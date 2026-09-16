# gaofenglib

**Status: retired and archived.** This repository is read-only and receives no
further development, bug fixes, or dependency updates. It is retained so that
code which already installs it keeps working.

The package has been superseded by
**[gflab/clinstats](https://github.com/gflab/clinstats)**, which is actively
maintained and carries the same statistical routines with input validation,
tidy output, tests, documentation, and continuous integration.

## What this package was

A personal collection of R helpers for clinical and translational research,
written between 2016 and 2020: receiver operating characteristic summaries,
odds ratios, univariable and multivariable Cox and logistic regression tables,
time-dependent survival cut points, resampling-based Cox screening, TCGA
clinical extraction, expression-matrix cleaning, the Oncotype DX colon cancer
recurrence score, and two R Markdown report templates.

## Where to go instead

```r
install.packages("remotes")
remotes::install_github("gflab/clinstats")
```

| Task | gaofenglib (retired) | clinstats (maintained) |
| --- | --- | --- |
| ROC summary with optimal cut point and power | `calc_logit()` | `roc_summary()` |
| Odds ratio at the optimal threshold | `calc_or()` | `odds_ratio()` |
| Univariable and multivariable Cox table | `factor_analysis_cox()` | `cox_table()` |
| Univariable and multivariable logistic table | `factor_analysis_logit()` | `logistic_table()` |
| Train and test model evaluation | `eval_logit()` | `evaluate_model()` |
| Time-dependent survival cut point | `calc_cutoff_survivalroc()` | `survival_cutoff()` |
| Resampling-based Cox screening | `calc_resamp_cox()` | `resample_cox()` |
| Expression-matrix cleaning | `clean_dat()` | `clean_expression()` |
| Oncotype DX colon cancer score | `calc_oncotypedx_crc()` | `oncotype_crc()` |
| TCGA clinical extraction | `extract_tcga_clinical()` | `tcga_clinical()` |
| Survival response object | — | `survival_response()` |
| Curated CRC clinical data | from `curatedClinicalData` | `clin_crc`, `clin_crc_gf`, `clin_crc_cell` |

## How to migrate

The two packages are not drop-in compatible, so each call site needs a small
edit. Changes to expect:

* **Return values are tibbles, not base data frames.** Results can be
  filtered, joined, and exported directly; positional indexing and `rownames()`
  no longer apply.
* **Regression tables include every requested variable in the multivariable
  model by default.** Use `multivariable = "significant"` to restore the
  previous rule, in which only variables with a univariable likelihood-ratio
  p-value below 0.05 entered the model.
* **Expression matrices keep samples in rows, as before.** `clean_expression()`
  and `oncotype_crc()` chain directly; pass `gene_axis = "rows"` to
  `oncotype_crc()` for a transposed matrix.
* **`resample_cox()` reports selection frequencies** per marker instead of a
  matrix of p-values, and follows the `future` parallel plan.
* **Examples for every function are on its help page**, together with the
  parameters that changed.

Before:

```r
library(gaofenglib)

rfs <- get_survival(clin, type = "rfs")
factor_analysis_cox(clin_factors, rfs, limit = 60, string = TRUE)
```

After:

```r
library(clinstats)

cox_table(
  clin,
  time = "rfs.delay",
  event = "rfs.event",
  factors = names(clin_factors),
  max_time = 60,
  multivariable = "significant"
)
```

## Installing the archived package

The package still installs. `remotes` resolves the Bioconductor dependency
`survcomp` automatically:

```r
install.packages("remotes")
remotes::install_github("gaofeng21cn/gaofenglib")
```

Because the repository is archived, the code is frozen at its last commit and
no support is offered. Work that must reproduce a previous result should pin
that commit rather than track the branch.

## Why it was retired

The statistical routines were useful but had accumulated real defects: several
exported functions depended on packages that were never declared, so a clean
installation could fail at run time; `init_R()` and `update_all()` reinstalled
the user's environment, which a library should not do; and the code carried no
tests. Those functions and the bundled datasets were consolidated into
`clinstats` rather than patched here, so that the lab maintains one package
instead of several overlapping ones.

## License

Apache License 2.0. See [LICENSE](LICENSE).
