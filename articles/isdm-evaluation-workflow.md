# ISDM Evaluation Workflow

``` r
library(isdmtools)
library(sf)
library(terra)
library(fmesher)
library(ggplot2)
library(inlabru)
#library(INLA)  # required
```

## Introduction

As demonstrated in the [Get
started](https://sodeidelphonse.github.io/isdmtools/articles/isdmtools.md)
guide, the first output from the `isdmtools` package is a set of clean
`sf` objects, which makes it easy to integrate with various spatial
modeling tools using block cross-validation techniques. The extracted
training and testing data can be directly fed into your preferred
integrated modeling tools such as `inlabru`, `PointedSDMs`, or any
`GLMs/GAMs` tools that can accommodate multisource spatial data. This
ensures that your model predictions are validated using a robust spatial
cross-validation approach and comprehensive evaluation metrics.

This vignette shows how how the isdmtools can be used with other
predictive modelling tools such as `inlabru` for a complete workflow of
integrated species distribution modelling (ISDM) analysis.

### Step 1: Fitting a Bayesian integrated model with the `inlabru`

The `inlabru` package is a wrapper for the `R-INLA` package which is
designed for Bayesian Latent Gaussian Modelling using INLA (Integrated
Laplace Nested Approximations) and Extensions. Let’s develop a Bayesian
spatial model with the simulated data presented in [Get
started](https://sodeidelphonse.github.io/isdmtools/articles/isdmtools.md)
guide.

- **Model definition**

We assume the following basic joint model with a shared latent signal
$\xi(.)$ (i.e. a Gaussian random field):

``` math
 \begin{matrix} 
 Y_{\mathrm{count},i}|\xi(.) \sim \mathrm{Pois} \left(\mu_i \right), \quad i = 1,\ldots,n,\\
 \log(\mu_i) = \beta_{0,\mathrm{count}} + \xi(\mathbf{s}_i)\\[3mm]
 X_{\mathrm{presence}}|\xi(.) \sim \mathrm{IPP} \left(\lambda(\mathbf{s}) \right),\\
\log (\lambda(\mathbf{s})) = \beta_{0,\mathrm{presence}} + \xi(\mathbf{s})\\
\end{matrix}
```

where $IPP$ means a *Inhomogeneous Poisson Process* and $\mathbf{s}$ the
vector of a location coordinates.

You can now prepare the remaining data required to fit an integrated
model. First, we need to define the mesh to be used to approximate the
latent field as well as for the integration points in the IPP
likelihood. After this step, we set up the prior distributions for the
latent parameters such as the range and the marginal standard deviation.
Then we define the observation model for each data type to fuse them
using a joint likelihood estimation with INLA.

- **Model implementation**

``` r
# Create a "mesh" for the latent field 
mesh <- fmesher::fm_mesh_2d(
     boundary = ben_sf,
     max.edge = c(0.2, 0.5),
     offset = c(1e-3, 0.6),
     cutoff = 0.10,
     crs = "epsg:4326"
)
ggplot() + inlabru::gg(mesh)
```

``` r
# Set the PC-prior for the SPDE model. We estimate a longer range value as no spatial 
# autocorrelation was defined in the data generation process:
pcmatern <- INLA::inla.spde2.pcmatern(mesh,
                                      prior.range = c(1, 0.1), # P(spatial range < 1) = 0.1
                                      prior.sigma = c(1, 0.1)  # P(sigma > 1) = 0.1
                                      )
   
# The shared spatial latent component is denoted by 'spde'
jcmp <- ~ -1 + Presence_intercept(1) + Count_intercept(1) +
                  spde(geometry, model = pcmatern)
   
# Count observation model
obs_model_count <- inlabru::bru_obs(
     formula = count ~  + Count_intercept + spde,
     family = "poisson",
     data = train_data$Count
 )
   
# Presence-only observation model (LGCP)
obs_model_pres <- inlabru::bru_obs(
     formula = geometry ~ Presence_intercept + spde,
     family = "cp",
     data = train_data$Presence,
     domain = list(geometry = mesh),
     samplers = list(geometry = ben_sf)
)
   
# Model fit
jfit <- inlabru::bru(jcmp, obs_model_count, obs_model_pres,
                        options = list(control.inla = list(int.strategy = "eb"),
                                       bru_max_iter = 20)
                    )
```

We can collect model results after the fitting process.

``` r
 jfit$summary.fixed
 #>                     mean        sd      0.025quant  0.5quant   0.975quant  mode      kld
 #> Count_intercept    -0.2497590 0.3086958 -0.8547916 -0.2497590  0.3552737  -0.2497590  0
 #> Presence_intercept  0.9269141 0.2836352  0.3709992  0.9269141  1.4828289   0.9269141  0
 #>  
 jfit$summary.hyperpar
 #>                mean        sd      0.025quant  0.5quant   0.975quant  mode
 #> Range for spde 3.535334 2.5240208  0.9513318   2.8572241  10.2509603  1.9527898
 #> Stdev for spde 0.512346 0.1926203  0.2183487   0.4842595  0.9647017   0.4317979
```

As expected, the estimated *spatial range* is higher than 1. This is
because there is no strong spatial autocorrelation in the simulated
data.

### Step 2: Model prediction

``` r
# Define the prediction grids and projection system
grids      <- fmesher::fm_pixels(mesh, mask = ben_sf)
projection <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
```

``` r
# Model predictions
jpred <- predict(jfit, 
                 newdata = grids, 
                 formula = ~ spde + Presence_intercept,
                 n.samples = 500, 
                 seed = 24)
jpred <- prepare_predictions(jpred) 

jpred_count <- predict(jfit, 
                       newdata = grids, 
                       formula = ~ spde + Count_intercept ,
                       n.samples = 500, seed = 24)
jpred_count <- prepare_predictions(jpred_count)
```

### Step 3: Habitat suitability analysis

Once we have obtained the model predictions, we can then use the
`isdmtools` to perform habitat suitability analysis. This will allow us
to proceed with the evaluation of models and the visualisation of
results.

``` r
# Probability of presence
jt_prob <- suitability_index(jpred, 
                            post_stat = c("q0.025", "mean", "q0.975"), 
                            output_format = "prob",
                            response_type = "joint.po",
                            projection = projection,
                            scale_independent = TRUE
                            )
plot(jt_prob)
```

``` r
# Expected counts
jt_count <- suitability_index(jpred_count, 
                              post_stat = c("q0.025", "mean", "q0.975"), 
                              output_format = "response",
                              response_type = "count",
                              projection = projection
                              )
plot(jt_count)
```

### Step 4: Model performance evaluation using the test data

Various performance metrics can now be computed, including
dataset-specific and weighted composite scores.

``` r
 xy_observed <- rbind(st_coordinates(datasets_list$Presence)[, c("X","Y")], 
              st_coordinates(datasets_list$Count)[datasets_list$Count$count > 0, c("X","Y")])
   
 metrics <- c("auc", "tss", "accuracy", "rmse", "mae")
 eval_metrics <- compute_metrics(test_data, 
                                prob_raster = jt_prob$mean, 
                                expected_response = jt_count$mean,
                                xy_excluded = xy_observed, 
                                metrics = metrics,
                                overall_roc_metrics = c("auc", "tss", "accuracy"),
                                response_counts = "count"
                                )
print(eval_metrics)

#> ISDM Model Evaluation Results
#> ----------------------------------------------
#> Datasets Evaluated: Presence, Count 

#> Overall Performance:
#>  TOT ROC SCORE     : 0.8048
#>  TOT ERROR SCORE   : 1.9353
#> ----------------------------------------------
```

One can obtain detailed overview of the evaluation results via the
[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
method.

``` r
summary(eval_metrics)

#> ==============================================
#>        ISDM EVALUATION SUMMARY REPORT
#> ==============================================
#> Generated on: 2026-01-19 04:41:02 

#> --- Model Evaluation Settings ---
#> Random Seed         : 25
#> Background Points   : 1000
#> Spatial Context     : BackgroundPoints object attached
#> Threshold Logic     : best
#> Optimality Criterion: youden
#> Prediction Type     : Absolute Count (No Offset)

#> --- Detailed Metric Table ---
#>         Presence Count
#> AUC         0.917 0.750
#> TSS         0.791 0.750
#> ACCURACY    0.794 0.778
#> RMSE          N/A 2.119
#> MAE           N/A 1.752

#> --- Composite Scores (Weighted) ---
#>     AUC      TSS ACCURACY     RMSE      MAE 
#>    0.852    0.775    0.788    2.119    1.752 

#> --- Overall Performance ---
#>  TOT ROC SCORE     : 0.8048
#>  TOT ERROR SCORE   : 1.9353
#> ==============================================
```

As you will have noticed, continuous-outcome metrics such as MAE (mean
absolute error) and RMSE (root mean squared error) are not available for
presence-only data, which makes sense. Furthermore, the weighted
composite scores for continuous responses are identical to their
individual counterparts, since there is only one count response.
Moreover, you can check out the
[`get_background()`](https://sodeidelphonse.github.io/isdmtools/reference/get_background.md)
documentation for more details on background sample generated during the
model evaluation.

Next, you can iterate through all five spatial folds to obtain an
average model performance, then calculate the variation in metrics
between blocks. Finally, run a model on the full ‘datasets_list’ to make
the final prediction.

### Step 5: Prediction mapping

You can now generate a formal prediction map ready for publication.

``` r
map <- generate_maps(jt_prob, 
                   var_names = c("q0.025", "mean", "q0.975"), 
                   base_map = ben_sf,
                   legend_title = "suitability",  
                   panel_labels = c("(a) q2.5%", "(b) Mean", "(c) q97.5%"),
                   xaxis_breaks = seq(0, 4, 1),
                   yaxis_breaks = seq(6, 13, 2)
                   )
map
```

## Conclusion

You have successfully fused multi-source biodiversity data and generated
spatially independent partitions for robust model validation using
`isdmtools`. The toolkit has enabled the resampling of data, thereby
reducing spatial autocorrelation effects in the modelling process. This
has facilitated the analysis and a comprehensive evaluation of ISDM
using well-known metrics from the fields of statistics and machine
learning.
