# hvtiRpropensity (unreleased)

* Fixes an empty changelog on the pkgdown site. `NEWS.md` opened with a bare
  `# hvtiRpropensity` title, which made level one the title level and left the
  version headings at level two. pkgdown reads the top heading level as the
  version level, so it found no releases at all and warned "no version headings
  found" on every build. The title is gone and the version headings are at
  level one, matching the rest of the family. `utils::news()` was unaffected
  and still reports the same four versions.

# hvtiRpropensity 0.1.3

* Removed the explicit `Maintainer:` field from `DESCRIPTION` and moved the
  maintainer address to `john.ehrlinger@gmail.com` in `Authors@R`, matching the
  rest of the `hvtiR*` family. The two fields had disagreed — `Authors@R` named
  the `ehrlinj@ccf.org` address while `Maintainer:` named the gmail one — which
  `R CMD check --as-cran` reported as a note. `Maintainer:` is now derived from
  the `cre` role at build time, so the two cannot drift apart again.

# hvtiRpropensity 0.1.2

* Qualified the `rnorm()` calls in the `sample_ps_data*()` generators as
  `stats::rnorm()`, matching every other statistics call in the package, and
  declared `stats` under `Imports`. This clears the `R CMD check` note about
  an undefined global function. Generated data is unchanged — the seeds and
  the RNG draw order are identical.

# hvtiRpropensity 0.1.1

* Renamed from `hvtiPropensityScores` into the `hvtiR*` package family
  (`hvtiRtemplates`, `hvtiRutilities`, `hvtiRlifetables`, `hvtiRtables`).
  No user-facing function names or behaviour changed.

# hvtiRpropensity 0.1.0

* Initial development scaffold.
* Added `sample_ps_data()` — reproducible synthetic cardiac-surgery dataset.
* Added `ps_match()` — greedy nearest-neighbour 1:1 propensity score matching
  without replacement, with optional caliper.
* Added `ps_weight()` — IPTW weighting supporting ATE, ATT, and ATC estimands,
  with optional stabilisation and weight winsorisation.
* Added `is_ps_data()` predicate and `print` / `summary` S3 methods for the
  common `ps_data` base class.
