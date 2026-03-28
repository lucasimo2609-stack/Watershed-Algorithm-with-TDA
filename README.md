# Persistent Watershed Segmentation in R

This repository implements an advanced image segmentation framework that integrates the classical **Watershed** transformation with **Persistent Homology** (Topological Data Analysis). 

The project focuses on solving the common "over-segmentation" problem by filtering image basins based on their topological importance, following the principles described by **H. Edelsbrunner and J. Harer** in *Computational Topology: An Introduction*.

## Repository Structure

The project is organized into the following key files:

* **Project_Report.pdf**: A comprehensive technical report explaining the mathematical theory of persistence applied to image segmentation and an analysis of the results.
* **ESAME_PALLINE.R**: R script applied to a real-world case study (`PALLINE.jpg`). It demonstrates how persistence filters out noise and small reflections to correctly segment circular objects.
* **ESAME_SINTETICO.R**: R script that generates a synthetic digital image. It provides a "controlled environment" to visualize how the algorithm handles peaks and valleys without external noise.

## Core Concept

In standard Watershed, every local minimum (even a single noisy pixel) creates a region. This project utilizes **0-dimensional persistence** to:
1.  **Measure Stability**: Calculate the "life-span" (Persistence = Death - Birth) of each basin.
2.  **Topological Pruning**: Only basins that are deeper than a specific threshold $\tau$ are kept.
3.  **Guided Segmentation**: These stable basins are used as unique markers (seeds) for the Watershed flooding process.

## Requirements

To run the scripts, ensure you have R installed along with the following packages:

```R
install.packages(c("imager", "TDA", "magrittr"))
