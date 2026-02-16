# Create block cross-validation folds for multisource spatial datasets

A constructor function for the `DataFolds` S3 class. It binds multiple
`sf` datasets into a single object and generates spatially or
environmentally-separated cross-validation folds.

## Usage

``` r
create_folds(
  datasets,
  region_polygon = NULL,
  k = 5,
  seed = 23,
  cv_method = "cluster",
  ...
)

bind_datasets(datasets)
```

## Arguments

- datasets:

  A named list of `sf` objects. Each list element should be a spatial
  dataset with its name corresponding to the list element's name. The
  list elements must have the same geographical projection and
  coordinates reference system (CRS).

- region_polygon:

  An `sf` object representing the study area polygon.

- k:

  integer. It specifies the number of folds (k-fold cross-validation).

- seed:

  integer. It sets seed for reproducibility.

- cv_method:

  character. It specifies the spatial cross-validation method to use.
  Options are `"cluster"` (default) or `"spatial"`, see
  [cv_cluster](https://rdrr.io/pkg/blockCV/man/cv_cluster.html) or
  [cv_spatial](https://rdrr.io/pkg/blockCV/man/cv_spatial.html)
  functions. For `"block"`, `"buffer"`, `"location"`, or `"nndm"`, see
  the corresponding functions in the `spatialsample` package.

- ...:

  Additional arguments to be passed to the underlying blocking function
  in `blockCV` or `spatialsample`.

## Value

An S3 object of class `DataFolds` containing the combined data, fold
information, the region polygon, and the original datasets.

An `sf` object containing all datasets along with their corresponding
name.

## Details

This function first binds all datasets into a single `sf` object. It
then applies the chosen blocking method to create spatial folds. The
fold IDs are added to the combined data object, and the original
datasets and other relevant information are stored in the returned
`DataFolds` object.

The `"cluster"` method from the `blockCV` package, supports both spatial
and environmental clustering. Methods `"block"`, `"buffer"`, `"nndm"`,
and `"location"`are useful for distance-based exclusion (buffering) or
leaving out specific groups/locations (e.g., using
`group = "column_name"` with the
[spatial_leave_location_out_cv](https://spatialsample.tidymodels.org/reference/spatial_vfold.html)
method). Use `"spatial"` for `blockCV` grid-blocking, or `"block"` for
`spatialsample` grid-blocking.

The behavior of `create_folds` depends on the `cv_method` chosen.
Several methods require specific arguments passed via the ellipsis
(`...`):

- **`block` & `buffer`**: Accept `radius` and `buffer` arguments to
  define the size of the test/assessment areas and the width of the
  exclusion zones, respectively. Particularly, for grid-blocking with
  `"block"` scheme, additional arguments to
  [st_make_grid](https://r-spatial.github.io/sf/reference/st_make_grid.html)
  can be provided.

- **`nndm`**: Requires `prediction_sites`, an `sf` object representing
  the area where the model will be projected.

- **`location`**: Requires a `group` argument (character). This should
  be the name of a column in your datasets representing independent
  units like `"site_id"`, `"year"`, or `"observer"`. This tests model
  generalizability across these factors. It also accepts `radius` and
  `buffer` arguments for the size of assessment areas and the width of
  the exclusion zones, respectively.

## References

Mahoney MJ, Johnson LK, Silge J, Frick H, Kuhn M, Beier CM. Assessing
the performance of spatial cross-validation approaches for models of
spatially structured data. *arXiv* (2023)
[doi:10.48550/arXiv.2303.07334](https://doi.org/10.48550/arXiv.2303.07334)

Roberts DR, Bahn V, Ciuti S, Boyce MS, Elith J, Guillera-Arroita G,
Hauenstein S, Lahoz-Monfort JJ, Schröder B, Thuiller W, et al.
Cross-validation strategies for data with temporal, spatial,
hierarchical, or phylogenetic structure. *Ecography* (2017) 40:913–929.
[doi:10.1111/ecog.02881](https://doi.org/10.1111/ecog.02881)

Valavi R, Elith J, Lahoz-Monfort JJ, Guillera-Arroita G. blockCV: an R
package for generating spatially or environmentally separated folds for
k-fold cross-validation of species distribution models. *bioRxiv*
(2018). [doi:10.1101/357798](https://doi.org/10.1101/357798)

## See also

Other spatial blocking methods:
[`DataFolds-methods`](https://sodeidelphonse.github.io/isdmtools/reference/DataFolds-methods.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Create some dummy sf data with different columns
library(sf)
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

# Create a list of datasets
datasets_list <- list(Presence = presence_data, Count = count_data)

# Create a dummy polygon for the region (e.g. Benin's minimum bounding rectangle)
ben_coords <- matrix(c(0, 6, 4, 6, 4, 13, 0, 13, 0, 6), ncol = 2, byrow = TRUE)
ben_sf <- st_sfc(st_polygon(list(ben_coords)), crs = 4326)
ben_sf <- st_sf(data.frame(name = "Benin"), ben_sf)

# Create a DataFolds object using the default 'cluster' method
my_folds <- create_folds(datasets_list, ben_sf, k = 5)
print(my_folds)
} # }
```
