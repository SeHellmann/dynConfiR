
# dynConfiR 0.0.2 (December 2022)
## Improvements
- added the DDMConf model (density: `dDDMConf`, RNG: `rDDMConf`), also for fitting and prediction.
- `d2DSD` and `dWEV` added argument `omega`. DynWEV and 2DSD are generalized with
the confidence measure explicitly depending on decision time. `omega` controls the 
power of decision time by which the confidence measure is divided. 
- `d2DSD` and `dWEV` added argument `stop_on_zero`: For the calculation of likelihoods, it
is useful to stop as soon as one probability is 0, since then the likelihood is 0.
- included starting point and drift rate variation in IRMt (experimental!, see `dIRM2`)
- fitting 2DSD and dynWEV: improved the finding of stating values for confidence thresholds.
- many bug fixes and increased robustness


# dynConfiR 0.0.1 (May 2022)
First CRAN release Version 0.0.1

