# Analysis of 3rd Ontario Breeding Bird Atlas 

Spatio-temporal models of breeding bird atlas point-count data. Observations are counts of birds detected during 5-minute unlimited-distance point counts. Counts are small (typically 0–4) and are modeled with a Poisson GLM with log link.

The linear predictor includes:

- Intercept
- Land-cover covariates with log-linear effects
- Spatial random field (SPDE) representing baseline spatial variation in abundance
- 1D SPDE smooth for time of day (hours since sunrise) to account for detectability variation
- 1D SPDE smooth for day of year to capture seasonal detectability changes
- Random effects for 10×10 km atlas squares

To model change between two atlas periods (~20 years apart), the model includes:
- A fixed atlas-period effect (overall mean shift)
- A second spatial SPDE field interacted with atlas period to capture spatially structured change
- A second atlas-square random effect indexed by square-within-atlas to capture atlas-specific local deviations
    - The square-level random effects are identifiable because one indexes the same square ID across atlases, and the other indexes square ID within atlas.
- Habitat covariates that are time-referenced, meaning their pixel values can differ between atlases. Covariate coefficients are shared across atlases, but habitat change can still alter expected counts through changes in covariate values.

Thus atlas differences arise from:
- habitat change (via covariate values),
- a province-wide atlas effect,
- spatially structured smooth change fields (SPDE),
- atlas-specific square effects.

The estimand of interest is expected count during a standardized 5-minute survey. Predictions are generated for 1-km grid cells, fixing time-of-day and day-of-year values and omitting square-level random effects so predictions represent expected counts at an average square.
