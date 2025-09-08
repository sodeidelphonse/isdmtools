# Presentation of isdmtools
`isdmtools` is an R package designed to streamline the process of preparing, evaluating and visualizing spatial data for biodiversity species distribution modeling, with a specific focus on **integrated species distribution models (ISDMs)** with multi-source geospatial datasets within a Bayesian framework. It provides a robust and reproducible workflow for spatial cross-validation, data management, and model evaluation, leveraging the power of sf, terra, dplyr, purrr, and ggplot2 packages.

# Installation

### How to install the package?
You can install the development version of `isdmtools` directly from GitHub using `devtools`.

```R
install.packages("devtools") 
devtools::install_github("your-github-username/isdmtools")
```
### How can contributors manage the package dependencies?
This project uses `renv` to manage package dependencies and ensure reproducibility. To install all the necessary packages for this project, simply follow these steps:

- Make sure you have the `renv` package installed:

```R
install.packages("renv")
```

- With the project directory as your working directory, use `renv` to install all packages listed in the `renv.lock` file:

```R
renv::restore()
```

# Key Features
The package provides a set of core functions to handle common data preparation and evaluation tasks:

**Data Preparation**: Create DataFolds objects that bind multiple `sf` datasets and generate spatially-separated cross-validation folds using the constructor function `create_folds()` developed with the `blockCV` package. This ensures the resulting models are robust to spatial autocorrelation.

**Suitability Analysis**: Standardize model predictions for consistent mapping and compute a final habitat suitability index. The `suitability_index()` function transforms raw integrated model predictions into a suitability score using the generalized complementary log-log transform, `cloglog`.

**Model Evaluation**: Compute comprehensive evaluation metrics, including ROC-based and continuous metrics through the `compute_metrics()` function. The package also handles *weighted composite scores*, providing a holistic view of model performance.

**Mapping & Visualization**: Visualize model predictions and final habitat suitability maps. The plotting method `generate_maps()` is designed with `ggplot2` to be clear and informative and visualize multiple variables of model predictions (e.g. mean, median, standard deviation or quantiles), providing an easy way to interpret models' results.

**S3 Methods**: The package includes `print()` and `plot()` methods for the `DataFolds` class, providing a concise summary and a clear visualization of the cross-validation partitions.

**sample_bg_points()**: A constructor function for generating background points for presence-only data with `print()` and `plot()` methods for the `BackgroundPoints` class.

# Getting Started: A Complete Example
The core workflow of "isdmtools" involves creating a DataFolds object and then extracting specific folds for your modeling pipeline.

First, let's load the package and set up some dummy data for a hypothetical study region.

```R
library(isdmtools)
library(sf)
library(ggplot2)
library(dplyr)

### Create dummy presence-only and count data
set.seed(42)
presence_data <- data.frame(
  x = runif(100, 0, 4),
  y = runif(100, 6, 13),
  site = rbinom(100, 1, 0.6)
) %>%
  st_as_sf(coords = c("x", "y"), crs = 4326)

count_data <- data.frame(
  x = runif(50, 0, 4),
  y = runif(50, 6, 13),
  count = rpois(50, 5)
) %>%
  st_as_sf(coords = c("x", "y"), crs = 4326)

### Create a list of datasets
datasets_list <- list(Presence = presence_data, Count = count_data)

### Define a dummy study region
ben_utm_coords <- matrix(c(0, 6, 4, 6, 4, 13, 0, 13, 0, 6), ncol = 2, byrow = TRUE)
ben_utm <- st_sfc(st_polygon(list(ben_utm_coords)), crs = 4326)
ben_utm <- st_sf(data.frame(name = "Region"), ben_utm)

### Create the DataFolds object
my_folds <- create_folds(datasets_list, ben_utm, k = 5, seed = 23)

### Print a summary of the object
print(my_folds)

### Visualize the folds
plot_cv <- plot(my_folds)
print(plot_cv)

### Extract a specific fold (e.g., Fold 3) for modeling
splits_fold_3 <- extract_fold(my_folds, fold = 3)

You now have two clean lists of `sf` objects for training and testing
### Access the data for the 'Presence' dataset in the training set
train_data_presence <- splits_fold_3\$train\$Presence

### Access the data for the 'Count' dataset in the testing set
test_data_count <- splits_fold_3\$test\$count
```

# Usage with Prediction Models
The output of `isdmtools` is a set of clean sf objects, which makes it easy to integrate with various modeling tools. The extracted train and test data can be directly fed into your preferred modeling packages such as `inlabru`, `PointedSDMs`, `MCMC` software, or any `GLMs/GAMs` tools that can accommodate single or multiple responses. This ensures that your model predictions are validated using a robust spatial cross-validation approach.

# Contributing
We welcome contributions! If you encounter an issue or have a feature request, please open an issue on the GitHub repository.

# Citation
To cite this package in your research work, please run the command below which will give both the plain text and `LaTex` version of the citation in your R session: 

```R
citation("isdmtools")
```
