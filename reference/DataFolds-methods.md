# Methods for manipulating data from `DataFolds` object.

- `extract_fold`: Extract a specific fold from a data partition.

- `plot`: A method to visualize the spatial blocks and the corresponding
  train/test partitions of observations for block cross-validation (CV).

- `autoplot`: A method for the native visualisation of `spatialsample`
  objects.

- `print`: A method to print folds' information per dataset, including
  the species geometry and the number of points excluded via spatial
  buffering.

- `summary`: A method to print a concise summary of the `DataFolds`
  object.

## Usage

``` r
extract_fold(object, fold, ...)

# S3 method for class 'DataFolds'
extract_fold(object, fold, ...)

# S3 method for class 'DataFolds'
plot(x, nrow = 1, annotate = TRUE, ...)

# S3 method for class 'DataFolds'
print(x, ...)

# S3 method for class 'DataFolds'
summary(object, ...)

# S3 method for class 'DataFolds'
autoplot(object, ...)
```

## Arguments

- object:

  An object from which a spatial fold can be extracted (e.g. `DataFolds`
  object)

- fold:

  An integer specifying the fold ID to be extracted as the test set.

- ...:

  Additional arguments passed to specific methods.

- x:

  A `DataFolds` object.

- nrow:

  An integer specifying the number of rows needed for the panel plot.
  The default is 1.

- annotate:

  If `TRUE`, the north arrow and scale bar are added to the plot.

## Value

- `extract_fold`: A list containing two named elements (`train` and
  `test`). For the `DataFolds` class, each element is a named list of
  `sf` objects, with names corresponding to the original datasets.

- `plot`: Returns a `ggplot` object that can be modified.

- `autoplot`: Returns a `ggplot` object.

- `print`: Invisibly returns the original object.

- `summary`: Invisibly returns the original object with a `table` of
  observations count per fold and dataset type.

## Details

The `extract_fold` method is generic, meaning it provides a consistent
interface for different types of objects. The method dispatched depends
on the class of the `object` argument. The primary purpose is to
abstract the process of accessing training and testing data for a given
fold, making it easier to write generic cross-validation loops. This
method splits the data into training and testing sets for a given fold,
using original datasets to ensure each returned `sf` object contains
only original columns.

## See also

Other spatial blocking methods:
[`create_folds()`](https://sodeidelphonse.github.io/isdmtools/reference/create_folds.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a list of some dummy sf data with different columns
library(sf)
library(dplyr)
set.seed(42)

presence_data <- data.frame(
  x = runif(100, 0, 4),
  y = runif(100, 6, 13),
  site = rbinom(100, 1, 0.6)
) %>% st_as_sf(coords = c("x", "y"), crs = 4326)

count_data <- data.frame(
  x = runif(50, 0, 4),
  y = runif(50, 6, 13),
  count = rpois(50, 5)
) %>% st_as_sf(coords = c("x", "y"), crs = 4326)

datasets_list <- list(Presence = presence_data, Count = count_data)

# Create a dummy polygon for the study region
ben_coords <- matrix(c(0, 6, 4, 6, 4, 13, 0, 13, 0, 6), ncol = 2, byrow = TRUE)
ben_sf <- st_sfc(st_polygon(list(ben_coords)), crs = 4326)
ben_sf <- st_sf(data.frame(name = "Benin"), ben_sf)

# Create a DataFolds object with the default 'cluster' method
my_folds <- create_folds(datasets_list, ben_sf, k = 5)
print(my_folds)

#-- Extract the desired fold (e.g. third fold for k-fold CV)
splits_fold_3 <- extract_fold(my_folds, fold = 3)

# Access the clean datasets for modeling
train_data_presence <- splits_fold_3$train$Presence
test_data_count <- splits_fold_3$test$Count

# Check the columns to see that they are correct
names(train_data_presence)
names(test_data_count)

#-- Plot the previous folds
plot_cv <- plot(my_folds)
print(plot_cv)

# You can even customize the plot (e.g. adjusting the axes breaks)
plot_cv <- plot_cv +
  ggplot2::scale_x_continuous(breaks = seq(0, 4, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(6, 13, 2))
print(plot_cv)

#-- Summarise folds information
summary(my_folds)

#-- Run the native autoplot for a \code{spatialsample} blocking method
fold_ss <- create_folds(datasets_list, ben_sf, cv_method = "block", k = 5)
autoplot(fold_ss)
plot(fold_ss)
} # }
```
