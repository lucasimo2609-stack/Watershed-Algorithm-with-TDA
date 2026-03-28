if (!requireNamespace("imager", quietly = TRUE)) install.packages("imager")
library(imager)
set.seed(42)

# immagine sintetica integrata
IMG_PATH <- NULL

INVERT         <- FALSE   # TRUE  = inverte l'immagine prima del watershed

SIGMA          <- 2.0    # blur globale sull'immagine:
#   0   = nessun blur
#   1-2 = per foto normali
#   3-5 = immagini rumorose
SIGMA_GRADIENT <- 0      # blur prima del gradiente (usato solo se USE_GRADIENT=TRUE)
#   0   = nessun blur aggiuntivo (consigliato con USE_GRADIENT=FALSE)
#   1-3 = per foto naturali con USE_GRADIENT=TRUE

TAU          <- 0.2 # soglia persistenza esplicita; NULL = automatica
TAU_QUANTILE <-   NULL   # usato solo se TAU = NULL [0,1]

USE_GRADIENT <- FALSE    
# FALSE = watershed sull'intensita'
#         microscopia, blob, nuclei, cellule (con INVERT=TRUE)
# TRUE  = watershed sul gradiente
#         foto naturali, immagini con bordi netti

# immagini troppo pesanti, ridimentisono
RESIZE_TO <- 150L      # ridimensiona l'immagine a max RESIZE_TO pixel per lato
# prima dell'elaborazione (mantiene le proporzioni).
#   0 o NULL = nessun ridimensionamento


SYN_NR    <- 100L      # righe immagine sintetica
SYN_NC    <- 130L      # colonne immagine sintetica
SYN_NOISE <- 0.035

# =============================================================================
#  CARICAMENTO IMMAGINE
# =============================================================================

load_image_as_matrix <- function(path) {
  img <- load.image(path)
  if (spectrum(img) > 1L) img <- grayscale(img)
  # as.data.frame(cimg) restituisce colonne: x, y, cc, value
  # x = colonna (1..W), y = riga (1..H)
  df  <- as.data.frame(img)
  W   <- max(df$x)
  H   <- max(df$y)
  # Ricostruisce la matrice: riga = y, colonna = x
  mat <- matrix(df$value, nrow = W, ncol = H, byrow = FALSE)
  # as.data.frame di imager ordina per x poi y -> trasposta necessaria
  t(mat)
}

# =============================================================================
#  RIDIMENSIONAMENTO BILINEARE IN R PURO
# =============================================================================

# Ridimensiona mat a (nr_new x nc_new) con interpolazione bilineare
resize_bilinear <- function(mat, nr_new, nc_new) {
  nr_old <- nrow(mat);  nc_old <- ncol(mat)
  if (nr_new == nr_old && nc_new == nc_old) return(mat)
  
  # Coordinate sorgente per ogni pixel di destinazione (0-based)
  row_src <- (seq(0, nr_new - 1L) / (nr_new - 1L)) * (nr_old - 1L)
  col_src <- (seq(0, nc_new - 1L) / (nc_new - 1L)) * (nc_old - 1L)
  
  r0 <- pmin(floor(row_src), nr_old - 2L);  r1 <- r0 + 1L
  c0 <- pmin(floor(col_src), nc_old - 2L);  c1 <- c0 + 1L
  dr <- row_src - r0;  dc <- col_src - c0
  
  out <- matrix(0, nr_new, nc_new)
  for (i in seq_len(nr_new)) {
    ri0 <- r0[i] + 1L;  ri1 <- r1[i] + 1L;  w_r <- dr[i]
    top    <- mat[ri0, c0 + 1L] * (1 - dc) + mat[ri0, c1 + 1L] * dc
    bottom <- mat[ri1, c0 + 1L] * (1 - dc) + mat[ri1, c1 + 1L] * dc
    out[i, ] <- top * (1 - w_r) + bottom * w_r
  }
  out
}

# Calcola le nuove dimensioni mantenendo le proporzioni
compute_new_dims <- function(nr, nc, max_side) {
  if (is.null(max_side) || max_side <= 0L || max(nr, nc) <= max_side)
    return(c(nr, nc))
  if (nr >= nc) c(max_side, max(1L, round(nc * max_side / nr)))
  else          c(max(1L, round(nr * max_side / nc)), max_side)
}

# =============================================================================
#  PRE-PROCESSING IN R PURO
# =============================================================================

# Kernel gaussiano 1-D di lunghezza (2*r+1)
gauss_kernel1d <- function(sigma) {
  r   <- ceiling(3 * sigma)
  idx <- seq(-r, r)
  k   <- exp(-idx^2 / (2 * sigma^2))
  k / sum(k)
}

# Convoluzione 1-D lungo le righe (padding con valori bordo)
conv_rows <- function(mat, k) {
  nr  <- nrow(mat)
  nc  <- ncol(mat)
  r   <- (length(k) - 1L) %/% 2L
  pad <- cbind(mat[, rep(1L, r), drop = FALSE],
               mat,
               mat[, rep(nc, r), drop = FALSE])
  out <- matrix(0, nrow = nr, ncol = nc)
  for (i in seq_along(k))
    out <- out + k[i] * pad[, i:(i + nc - 1L), drop = FALSE]
  out
}

# Smoothing gaussiano 2-D separabile
gaussian_blur <- function(mat, sigma) {
  k   <- gauss_kernel1d(sigma)
  tmp <- conv_rows(mat, k)              # blur sulle righe
  t(conv_rows(t(tmp), k))              # blur sulle colonne
}

# Gradiente Sobel 3x3 migliore al rumore delle immagini reali
# rispetto alle differenze finite semplici.
# Sobel X:  [-1 0 +1]     Sobel Y:  [-1 -2 -1]
#           [-2 0 +2]               [ 0  0  0]
#           [-1 0 +1]               [+1 +2 +1]
sobel_magnitude <- function(mat) {
  nr <- nrow(mat);  nc <- ncol(mat)
  
  # Padding con replica del bordo (1 pixel su ogni lato)
  p <- matrix(0, nr + 2L, nc + 2L)
  p[2:(nr+1L), 2:(nc+1L)] <- mat
  p[1L,        ]  <- p[2L,       ]   # bordo superiore
  p[nr+2L,     ]  <- p[nr+1L,    ]   # bordo inferiore
  p[, 1L       ]  <- p[, 2L      ]   # bordo sinistro
  p[, nc+2L    ]  <- p[, nc+1L   ]   # bordo destro
  
  # Indici del centro dopo il padding
  r <- 2:(nr+1L);  c <- 2:(nc+1L)
  
  # Sobel X
  gx <- (-1*p[r-1L,c-1L] + 1*p[r-1L,c+1L]
         -2*p[r,   c-1L] + 2*p[r,   c+1L]
         -1*p[r+1L,c-1L] + 1*p[r+1L,c+1L])
  
  # Sobel Y
  gy <- (-1*p[r-1L,c-1L] - 2*p[r-1L,c] - 1*p[r-1L,c+1L]
         +1*p[r+1L,c-1L] + 2*p[r+1L,c] + 1*p[r+1L,c+1L])
  
  matrix(sqrt(gx^2 + gy^2), nrow = nr, ncol = nc)
}

# Normalizza una matrice in [0, 1]
normalize01 <- function(mat) {
  mn <- min(mat);  mx <- max(mat)
  matrix((mat - mn) / (mx - mn + 1e-12), nrow = nrow(mat), ncol = ncol(mat))
}

# =============================================================================
#  creo il VICINATO A 8
# =============================================================================

neighbors8 <- function(idx, nr, nc) {
  r  <- ((idx - 1L) %% nr) + 1L
  c  <- ((idx - 1L) %/% nr) + 1L
  dr <- c(-1L,-1L,-1L, 0L, 0L, 1L, 1L, 1L)
  dc <- c(-1L, 0L, 1L,-1L, 1L,-1L, 0L, 1L)
  r2 <- r + dr;  c2 <- c + dc
  ok <- r2 >= 1L & r2 <= nr & c2 >= 1L & c2 <= nc
  (c2[ok] - 1L) * nr + r2[ok]
}

# =============================================================================
#  UNION-FIND (path compression + union by rank + birth tracking)
# =============================================================================

#La Union-Find serve a tenere traccia di quali pixel appartengono allo stesso bacino 
#mentre l'algoritmo scansiona l'immagine in ordine crescente di intensità, 
#senza dover riesaminare tutti i pixel ogni volta

new_uf <- function() {
  e <- new.env(parent = emptyenv())
  e$parent <- integer(0);  e$rank <- integer(0)
  e$birth  <- numeric(0);  e$n    <- 0L
  e
}

uf_new_component <- function(uf, birth_val) {
  uf$n      <- uf$n + 1L
  uf$parent <- c(uf$parent, uf$n)
  uf$rank   <- c(uf$rank,   0L)
  uf$birth  <- c(uf$birth,  birth_val)
  uf$n
}

uf_find <- function(uf, x) {
  root <- x
  while (uf$parent[root] != root) root <- uf$parent[root]
  while (uf$parent[x] != root) { nxt <- uf$parent[x]; uf$parent[x] <- root; x <- nxt }
  root
}

uf_unite <- function(uf, x, y, saddle_val) {
  rx <- uf_find(uf, x);  ry <- uf_find(uf, y)
  if (rx == ry) return(NULL)
  if (uf$birth[rx] >= uf$birth[ry]) { dying <- rx; surviving <- ry
  } else                             { dying <- ry; surviving <- rx }
  pair <- list(dying_id = dying, surviving_id = surviving,
               birth = uf$birth[dying], death = saddle_val,
               persistence = saddle_val - uf$birth[dying])
  if (uf$rank[surviving] < uf$rank[dying]) { tmp <- surviving; surviving <- dying; dying <- tmp }
  uf$parent[dying] <- surviving
  if (uf$rank[surviving] == uf$rank[dying]) uf$rank[surviving] <- uf$rank[surviving] + 1L
  uf$birth[surviving] <- min(uf$birth[surviving], uf$birth[dying])
  pair
}

# =============================================================================
#  WATERSHED COMPLETO + COPPIE DI PERSISTENZA
# =============================================================================

compute_persistence_pairs <- function(fmat) {
  stopifnot(is.matrix(fmat), is.numeric(fmat))
  nr <- nrow(fmat);  nc <- ncol(fmat);  N <- nr * nc
  uf <- new_uf();  labels <- integer(N)
  seed_pixel <- integer(0);  plist <- vector("list", N);  n_pairs <- 0L
  
  for (idx in order(fmat)) {
    nbrs       <- neighbors8(idx, nr, nc)
    labeled_nb <- nbrs[labels[nbrs] > 0L]
    
    if (length(labeled_nb) == 0L) {
      cid <- uf_new_component(uf, fmat[idx])
      labels[idx] <- cid
      seed_pixel  <- c(seed_pixel, idx)
    } else {
      # Radici correnti dei bacini vicini gia' etichettati
      roots <- unique(vapply(labels[labeled_nb],
                             function(l) uf_find(uf, l), integer(1L)))
      
      # Valori di nascita: uf_find puo' restituire indici aggiornati
      # dopo la path-compression; ri-risolviamo le radici per sicurezza
      roots  <- unique(vapply(roots, function(r) uf_find(uf, r), integer(1L)))
      births <- vapply(roots, function(r) {
        b <- uf$birth[r]
        if (is.null(b) || length(b) == 0L || is.na(b)) Inf else b
      }, numeric(1L))
      
      # Scegli il bacino con il birth piu' antico come ricevente
      best_i    <- which.min(births)
      if (length(best_i) == 0L) best_i <- 1L   # fallback di sicurezza
      main_root <- roots[best_i]
      labels[idx] <- main_root
      
      # Fondi tutti gli altri bacini nel ricevente
      for (r in roots) {
        r1 <- uf_find(uf, main_root)
        r2 <- uf_find(uf, r)
        if (r1 != r2) {
          p <- uf_unite(uf, r1, r2, fmat[idx])
          if (!is.null(p)) { n_pairs <- n_pairs + 1L; plist[[n_pairs]] <- p }
        }
      }
    }
  }
  
  pairs_df <- if (n_pairs > 0L) {
    do.call(rbind, lapply(plist[seq_len(n_pairs)], as.data.frame))
  } else {
    data.frame(dying_id = integer(0), surviving_id = integer(0),
               birth = numeric(0), death = numeric(0), persistence = numeric(0))
  }
  list(pairs = pairs_df, seed_pixel = seed_pixel, nr = nr, nc = nc)
}

# =============================================================================
#  WATERSHED SEMINATO
# =============================================================================

seeded_watershed <- function(fmat, seeds_idx) {
  stopifnot(is.matrix(fmat))
  nr <- nrow(fmat);  nc <- ncol(fmat);  N <- nr * nc
  n_seeds <- length(seeds_idx)
  if (n_seeds == 0L) stop("Nessun seme sopravvissuto.")
  
  labels <- integer(N)
  for (k in seq_len(n_seeds)) labels[seeds_idx[k]] <- k
  
  pq <- data.frame(val = fmat[seeds_idx], idx = seeds_idx,
                   lbl = seq_len(n_seeds), stringsAsFactors = FALSE)
  
  while (nrow(pq) > 0L) {
    best <- which.min(pq$val)
    if (length(best) == 0L || is.na(best)) break
    cur  <- pq[best, ];  pq <- pq[-best, , drop = FALSE]
    cidx <- cur$idx;  clbl <- cur$lbl
    if (labels[cidx] != 0L && labels[cidx] != clbl) next
    labels[cidx] <- clbl
    nbrs      <- neighbors8(cidx, nr, nc)
    unlabeled <- nbrs[labels[nbrs] == 0L]
    if (length(unlabeled) == 0L) next
    labels[unlabeled] <- clbl
    pq <- rbind(pq, data.frame(val = fmat[unlabeled], idx = unlabeled,
                               lbl = clbl, stringsAsFactors = FALSE))
  }
  matrix(labels, nrow = nr, ncol = nc)
}

# =============================================================================
#  FUNZIONE PRINCIPALE
# =============================================================================

segment_persistence_watershed <- function(img_path       = NULL,
                                          img_mat        = NULL,
                                          invert         = FALSE,
                                          sigma          = 0,
                                          sigma_gradient = 0,
                                          resize_to      = NULL,
                                          tau            = NULL,
                                          tau_quantile   = 0.85,
                                          use_gradient   = FALSE) {
  # -- Input ------------------------------------------------------------------
  if (!is.null(img_path)) {
    if (!file.exists(img_path)) stop("File non trovato: ", img_path)
    mat <- load_image_as_matrix(img_path)
  } else if (!is.null(img_mat)) {
    mat <- img_mat
  } else {
    stop("Fornire img_path oppure img_mat.")
  }
  
  mat <- matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat))
  nr_orig <- nrow(mat);  nc_orig <- ncol(mat)
  cat(sprintf("\nImmagine originale: %d righe x %d colonne\n", nr_orig, nc_orig))
  
  # -- Ridimensionamento (se richiesto) --------------------------------------
  new_dims <- compute_new_dims(nr_orig, nc_orig, resize_to)
  if (new_dims[1] != nr_orig || new_dims[2] != nc_orig) {
    mat <- resize_bilinear(mat, new_dims[1], new_dims[2])
    cat(sprintf("  Ridimensionata a: %d x %d pixel\n", new_dims[1], new_dims[2]))
  }
  
  nr_img <- nrow(mat);  nc_img <- ncol(mat)
  
  # -- Normalizza [0,1] -------------------------------------------------------
  mat <- normalize01(mat)
  
  # -- Inversione opzionale ---------------------------------------------------
  # Il watershed trova i MINIMI. Se gli oggetti di interesse sono CHIARI
  # (microscopia a fluorescenza, blob bianchi) occorre invertire l'immagine
  # cosi' che ogni oggetto diventi una valle = un bacino separato.
  if (invert) {
    mat <- 1 - mat
    cat("  Immagine invertita (oggetti chiari -> minimi)\n")
  }
  
  # -- Smoothing globale opzionale (sigma = 0 -> salta) ---------------------
  if (sigma > 0) {
    mat_s <- matrix(as.numeric(gaussian_blur(mat, sigma)), nrow = nr_img, ncol = nc_img)
    cat(sprintf("  Smoothing globale: sigma = %.2f\n", sigma))
  } else {
    mat_s <- mat
  }
  
  # -- Funzione f da segmentare ----------------------------------------------
  if (use_gradient) {
    # Blur dedicato prima del gradiente (piu' importante per immagini reali)
    mat_for_grad <- if (sigma_gradient > 0) {
      cat(sprintf("  Blur pre-gradiente: sigma_gradient = %.2f\n", sigma_gradient))
      matrix(as.numeric(gaussian_blur(mat_s, sigma_gradient)), nrow = nr_img, ncol = nc_img)
    } else {
      cat("  Nessun blur pre-gradiente.\n")
      mat_s
    }
    # Sobel: robusto al rumore, bordi netti sulle immagini reali
    fmat <- sobel_magnitude(mat_for_grad)
    cat("  Gradiente: Sobel 3x3\n")
  } else {
    fmat <- mat_s
    cat("  Funzione f: intensita' (no gradiente)\n")
  }
  fmat <- normalize01(fmat)
  fmat <- matrix(as.numeric(fmat), nrow = nr_img, ncol = nc_img)
  
  # -- Watershed + persistenza -----------------------------------------------
  cat("Calcolo coppie di persistenza...\n")
  ws    <- compute_persistence_pairs(fmat)
  prs   <- ws$pairs
  n_min <- length(ws$seed_pixel)
  cat(sprintf("  Minimi locali:      %d\n", n_min))
  cat(sprintf("  Coppie persistenza: %d\n", nrow(prs)))
  
  # -- Soglia tau ------------------------------------------------------------
  if (is.null(tau)) {
    tau <- if (nrow(prs) > 0L) as.numeric(quantile(prs$persistence, tau_quantile)) else 0
    cat(sprintf("  Soglia auto (q=%.0f%%): tau = %.5f\n", tau_quantile * 100, tau))
  } else {
    cat(sprintf("  Soglia fornita: tau = %.5f\n", tau))
  }
  
  # -- Semplificazione -------------------------------------------------------
  cancelled     <- if (nrow(prs) > 0L) prs$dying_id[prs$persistence < tau] else integer(0)
  surviving_ids <- setdiff(seq_len(n_min), cancelled)
  seeds         <- ws$seed_pixel[surviving_ids]
  cat(sprintf("  Regioni finali: %d  (su %d minimi)\n", length(seeds), n_min))
  
  # -- Watershed seminato ----------------------------------------------------
  cat("Watershed seminato...\n")
  lbls <- seeded_watershed(fmat, seeds)
  
  list(original = mat, smoothed = mat_s, gradient = fmat,
       pairs = prs, tau = tau, n_segments = length(seeds), labels = lbls)
}

# =============================================================================
#  UTILITA' PER rasterImage
#  rasterImage vuole array[colonne, righe, RGB] con riga 1 = in basso.
#  Partendo da mat[nr x nc] (riga 1 = in alto):
#    flip verticale -> t() -> [nc x nr]  (colonne x righe) OK
# =============================================================================

gray_to_raster <- function(mat) {
  nr <- nrow(mat);  nc <- ncol(mat)
  t(mat[seq(nr, 1L), , drop = FALSE])
}

labels_to_raster <- function(lbl, pal) {
  nr  <- nrow(lbl);  nc  <- ncol(lbl)
  hex <- pal[as.vector(lbl[seq(nr, 1L), , drop = FALSE])]
  rv  <- col2rgb(hex) / 255
  a   <- array(0, dim = c(nc, nr, 3L))
  a[,,1L] <- t(matrix(rv[1L,], nrow = nr, ncol = nc))
  a[,,2L] <- t(matrix(rv[2L,], nrow = nr, ncol = nc))
  a[,,3L] <- t(matrix(rv[3L,], nrow = nr, ncol = nc))
  a
}

# =============================================================================
#  VISUALIZZAZIONE
# =============================================================================

plot_results <- function(res, title_prefix = "") {
  stopifnot(is.matrix(res$original), is.matrix(res$gradient), is.matrix(res$labels))
  op <- par(mfrow = c(2L, 3L), mar = c(2.5, 2.5, 3.2, 1.0), oma = c(0,0,2.5,0))
  on.exit(par(op))
  g256 <- gray.colors(256L, start = 0, end = 1)
  
  # Originale
  image(gray_to_raster(res$original), col = g256,
        axes = FALSE, asp = 1, main = paste(title_prefix, "Originale"))
  
  # Funzione f
  image(gray_to_raster(res$gradient), col = g256,
        axes = FALSE, asp = 1, main = "Funzione f  (gradiente / intensita')")
  
  # Persistence Diagram
  prs <- res$pairs
  if (nrow(prs) > 0L) {
    vivi <- prs$persistence >= res$tau
    plot(prs$birth, prs$death, pch = 16, cex = ifelse(vivi, 0.8, 0.45),
         col = ifelse(vivi, "#C0392B", "#7FB3D3"),
         xlab = "Birth", ylab = "Death", main = "Persistence Diagram",
         xlim = c(0,1), ylim = c(0,1))
    abline(0, 1, lty = 2, col = "grey55", lwd = 1.2)
    abline(a = res$tau, b = 1, lty = 2, col = "#27AE60", lwd = 2)
    legend("topleft", bty = "n", cex = 0.75,
           legend = c(sprintf("Cancellate   (p < %.3f)", res$tau),
                      sprintf("Sopravvissute (p >= %.3f)", res$tau)),
           col = c("#7FB3D3","#C0392B"), pch = 16)
  } else { plot.new(); title("Nessuna coppia") }
  
  # Barcode
  if (nrow(prs) > 0L) {
    ord  <- order(prs$persistence);  po <- prs[ord,];  n_p <- nrow(po)
    bcol <- ifelse(po$persistence >= res$tau, "#C0392B", "#7FB3D3")
    plot(NULL, xlim = range(c(po$birth, po$death)), ylim = c(0, n_p+1L),
         xlab = "Intensita'", ylab = "", main = "Barcode di Persistenza", yaxt = "n")
    for (i in seq_len(n_p))
      lines(c(po$birth[i], po$death[i]), c(i, i), col = bcol[i], lwd = 1.6)
    abline(v = res$tau, col = "#27AE60", lty = 2, lwd = 2)
  } else { plot.new(); title("Nessun barcode") }
  
  # Mappa segmentazione
  n_seg <- res$n_segments
  pal   <- colorRampPalette(c("#E74C3C","#3498DB","#2ECC71","#F39C12",
                              "#9B59B6","#1ABC9C","#E67E22","#34495E",
                              "#E91E63","#00BCD4"))(max(n_seg, 1L))
  plot(c(0,1), c(0,1), type="n", axes=FALSE, asp=1,
       main = sprintf("Segmentazione: %d regioni  (tau=%.3f)", n_seg, res$tau))
  rasterImage(labels_to_raster(res$labels, pal), 0, 0, 1, 1)
  
  # Overlay bordi
  lbl  <- res$labels;  nr_l <- nrow(lbl);  nc_l <- ncol(lbl)
  ri   <- seq_len(nr_l - 1L);  ci <- seq_len(nc_l - 1L)
  bnd  <- matrix(FALSE, nr_l, nc_l)
  bnd[ri,  ]      <- bnd[ri,  ]      | (lbl[ri,]      != lbl[ri+1L,])
  bnd[ri+1L,]     <- bnd[ri+1L,]     | (lbl[ri,]      != lbl[ri+1L,])
  bnd[,ci]        <- bnd[,ci]        | (lbl[,ci]       != lbl[,ci+1L])
  bnd[,ci+1L]     <- bnd[,ci+1L]     | (lbl[,ci]       != lbl[,ci+1L])
  
  ov <- array(rep(as.vector(res$original), 3L), dim = c(nr_l, nc_l, 3L))
  ov[,,1L][bnd] <- 1;  ov[,,2L][bnd] <- 0.08;  ov[,,3L][bnd] <- 0.08
  ov_r <- array(0, dim = c(nc_l, nr_l, 3L))
  for (ch in 1L:3L) ov_r[,,ch] <- t(ov[seq(nr_l,1L),,ch])
  
  plot(c(0,1), c(0,1), type="n", axes=FALSE, asp=1, main="Overlay bordi")
  rasterImage(ov_r, 0, 0, 1, 1)
  
  mtext(paste("Watershed + Persistenza Topologica", title_prefix),
        outer = TRUE, cex = 1.1, font = 2)
  invisible(res)
}

# =============================================================================
#  IMMAGINE SINTETICA  (3 valli evidenti + 1 no)
# =============================================================================

make_synthetic_image <- function(nr = 100L, nc = 130L, noise_sd = 0.035) {
  x <- matrix(seq(0,1, length.out=nc), nrow=nr, ncol=nc, byrow=TRUE)
  y <- matrix(seq(0,1, length.out=nr), nrow=nr, ncol=nc, byrow=FALSE)
  mat <- 1 -
    0.85*exp(-((x-0.20)^2+(y-0.28)^2)/0.011) -
    0.75*exp(-((x-0.72)^2+(y-0.32)^2)/0.010) -
    0.60*exp(-((x-0.45)^2+(y-0.75)^2)/0.009) -
    0.25*exp(-((x-0.78)^2+(y-0.80)^2)/0.005) +
    noise_sd * matrix(rnorm(nr*nc), nrow=nr, ncol=nc)
  matrix(pmax(0, pmin(1, mat)), nrow=nr, ncol=nc)
}

# =============================================================================
#  ESECUZIONE PRINCIPALE
# =============================================================================

cat("\n=====================================================\n")
cat("  Watershed con Persistenza Topologica\n")
cat("=====================================================\n")

if (!is.null(IMG_PATH)) {
  cat(sprintf("\nCaricamento: %s\n", IMG_PATH))
  res <- segment_persistence_watershed(
    img_path     = IMG_PATH,
    invert         = INVERT,
    sigma          = SIGMA,
    sigma_gradient = SIGMA_GRADIENT,
    resize_to      = RESIZE_TO,
    tau            = TAU,
    tau_quantile   = TAU_QUANTILE,
    use_gradient   = USE_GRADIENT
  )
  lbl_prefix <- basename(IMG_PATH)
  
} else {
  cat("\nIMG_PATH = NULL: uso immagine sintetica.\n")
  img_syn <- make_synthetic_image(nr=SYN_NR, nc=SYN_NC, noise_sd=SYN_NOISE)
  res <- segment_persistence_watershed(
    img_mat      = img_syn,
    invert         = INVERT,
    sigma          = SIGMA,
    sigma_gradient = SIGMA_GRADIENT,
    resize_to      = RESIZE_TO,
    tau            = TAU,
    tau_quantile   = TAU_QUANTILE,
    use_gradient   = USE_GRADIENT
  )
  lbl_prefix <- "Sintetica"
}

# -- pannelli principali ---------------------------------------------------
dev.new(width=15, height=9, noRStudioGD=TRUE)
plot_results(res, title_prefix = paste("--", lbl_prefix))

# =============================================================================
#  RIEPILOGO
# =============================================================================

cat("\n=====================================================\n")
cat("  RIEPILOGO\n")
cat("=====================================================\n")
prs <- res$pairs
cat(sprintf("  Sorgente              : %s\n", if (!is.null(IMG_PATH)) IMG_PATH else "sintetica"))
cat(sprintf("  Dimensione immagine   : %d x %d pixel\n", nrow(res$original), ncol(res$original)))
cat(sprintf("  Sigma globale         : %.2f\n",  SIGMA))
cat(sprintf("  Sigma pre-gradiente   : %.2f\n",  SIGMA_GRADIENT))
cat(sprintf("  Minimi locali totali  : %d\n",    nrow(prs) + res$n_segments))
cat(sprintf("  Coppie di persistenza : %d\n",    nrow(prs)))
cat(sprintf("  Soglia tau            : %.5f\n",  res$tau))
cat(sprintf("  Regioni finali        : %d\n",    res$n_segments))
if (nrow(prs) > 0L)
  cat(sprintf("  Persist. min/med/max  : %.4f / %.4f / %.4f\n",
              min(prs$persistence), median(prs$persistence), max(prs$persistence)))
cat("=====================================================\n")

