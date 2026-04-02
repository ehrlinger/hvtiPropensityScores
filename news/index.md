# Changelog

## hvtiPropensityScores (development version)

### hvtiPropensityScores 0.0.0.9000

- Initial development scaffold.
- Added
  [`sample_ps_data()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/sample_ps_data.md)
  — reproducible synthetic cardiac-surgery dataset.
- Added
  [`ps_match()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/ps_match.md)
  — greedy nearest-neighbour 1:1 propensity score matching without
  replacement, with optional caliper.
- Added
  [`ps_weight()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/ps_weight.md)
  — IPTW weighting supporting ATE, ATT, and ATC estimands, with optional
  stabilisation and weight winsorisation.
- Added
  [`is_ps_data()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/is_ps_data.md)
  predicate and `print` / `summary` S3 methods for the common `ps_data`
  base class.
