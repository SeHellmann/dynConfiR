# dynConfiR 0.0.3 (April 2023)
## Improvements
- added the dynaViTE model (generalization of dynWEV with time-dependent 
confidence variable), which includes the additional lambda which can also be fitted
- improved fitting procedure for dynaViTE/dynWEV and 2DSD: confidence thresholds for 
starting parameters in the initial grid search are now estimated using simulations; 
this decreases the size of the initial grid drastically, while improving the initial
values
- adapted code to changes in the new dplyr release (1.1.1)


# dynConfiR 0.0.2 (December 2022)
## Improvements
- added the DDMConf model (density: `dDDMConf`, RNG: `rDDMConf`), also for fitting and prediction.
- `d2DSD` and `dWEV` added argument `lambda`. DynWEV and 2DSD are generalized with
the confidence measure explicitly depending on decision time. `lambda` controls the 
power of decision time by which the confidence measure is divided. 
- `d2DSD` and `dWEV` added argument `stop_on_zero`: For the calculation of likelihoods, it
is useful to stop as soon as one probability is 0, since then the likelihood is 0.
- included starting point and drift rate variation in IRMt (experimental!, see `dIRM2`)
- fitting 2DSD and dynWEV: improved the finding of stating values for confidence thresholds.
- many bug fixes and increased robustness


# dynConfiR 0.0.1 (May 2022)
First CRAN release Version 0.0.1

