# gaofenglib

> **Retired.** This package is kept for reproducibility of code that already
> installs it. It is superseded by
> [gflab/clinstats](https://github.com/gflab/clinstats), which carries the
> same statistical routines with input validation, tidy output, tests, and
> documentation. See the migration table in the `clinstats` README.
>
> The install command below continues to work and no new development is
> planned here.

Useful R functions for personal use.

## Install

```
install.packages("devtools")
devtools::install_github("gaofeng21cn/gaofenglib")
```

## Successor

[gflab/clinstats](https://github.com/gflab/clinstats) replaces this package.
The function mapping is:

| gaofenglib | clinstats |
| --- | --- |
| `calc_logit()` | `roc_summary()` |
| `calc_or()` | `odds_ratio()` |
| `factor_analysis_cox()` | `cox_table()` |
| `factor_analysis_logit()` | `logistic_table()` |
| `eval_logit()` | `evaluate_model()` |
| `calc_cutoff_survivalroc()` | `survival_cutoff()` |
| `calc_resamp_cox()` | `resample_cox()` |
| `clean_dat()` | `clean_expression()` |
| `calc_oncotypedx_crc()` | `oncotype_crc()` |
| `extract_tcga_clinical()` | `tcga_clinical()` |

`init_R()` and `update_all()` were dropped: a package should not reinstall the
user's environment. The `clin_crc` datasets previously installed from
`curatedClinicalData` ship with `clinstats`.

#### Functions

- `init_R`


#### Templates
- `Analysis report(html)`
- `Analysis report(pdf)`
