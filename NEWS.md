# dynConfiR 1.0.0 (May 2025)
- change in the precision argument for all models. It is now an integer argument representing the approximate precision of the density in digits. Default is 6 for all models in the density functions.
- Added precision argument in the fitting functions. With default to 3
- renamed the DDMConf model to DDConf to align with published literature
- renamed the functions for dynaViTE model (ddynaViTE and rdynaViTE);
- Added functions for quantitative model comparison (see `?group_BMS?`)
- the functions dWEV and rWEV will be kept for now but maybe removed in future releases of the package, they will produce a deprecation warning now
- Minor changes: 
  - Bug fix in fitting Race Models with time-dependend confidence variable. Before, the weight 'wrt' was bound from above by the weight 'wx' because of a bug. 
  - Added a UserInterrupt() call in all longer-running C-functions. 

# dynConfiR 0.0.4 (January 2024)
Fixed a CRAN note and improved robustness in fitting functions, and simulations. 
Also improved initial parameter sets for model fitting.


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
- `d2DSD` and `ddynaViTE` added argument `lambda`. DynWEV and 2DSD are generalized with
the confidence measure explicitly depending on decision time. `lambda` controls the 
power of decision time by which the confidence measure is divided. 
- `d2DSD` and `ddynaViTE` added argument `stop_on_zero`: For the calculation of likelihoods, it
is useful to stop as soon as one probability is 0, since then the likelihood is 0.
- included starting point and drift rate variation in IRMt (experimental!, see `dIRM2`)
- fitting 2DSD and dynWEV: improved the finding of stating values for confidence thresholds.
- many bug fixes and increased robustness


# dynConfiR 0.0.1 (May 2022)
First CRAN release Version 0.0.1

