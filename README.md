# SAMutils

`SAMutils` provides helper functions for Stock Synthesis Assessment Model (SAM) workflows built around the [`stockassessment`](https://github.com/fishfollower/SAM) R package.

The package focuses on:

- preparing SAM-ready input matrices from tidy data
- fitting and comparing model configurations
- extracting model outputs in tidy formats
- producing standard diagnostic plots and summary tables

## Installation

### Install from local source

From the project root:

```r
R CMD INSTALL .
```

Or from R:

```r
install.packages(".", repos = NULL, type = "source")
```

### Install development dependencies (optional)

```r
install.packages(c("devtools", "testthat"))
```

## Usage

Load the package:

```r
library(SAMutils)
```

### 1. Build SAM input matrices from tidy data

```r
model_dat <- tibble::tibble(
  year = c(2000, 2000, 2001, 2001),
  age = c(1, 2, 1, 2),
  M = c(0.20, 0.18, 0.21, 0.19)
)

M_mat <- sam.input(
  model_dat = model_dat,
  variable = "M",
  age_range = 1:2,
  year_range = 2000:2001,
  tail_f = mean
)
```

### 2. Run a SAM fit and extract tidy outputs

```r
# fit <- stockassessment::sam.fit(dat, conf)

# Yearly summaries (SSB, Fbar, recruitment, catch, etc.)
# by_year <- rby.sam(fit)

# Year-age summaries (N, F, optional stock attributes)
# by_year_age <- rbya.sam(fit)
```

### 3. Plot diagnostics

```r
# Observation residual bubble plot
# p_res <- model_resid_plot(list(res = residuals(fit)))

# Process residual bubble plot
# p_pres <- model_pres_resid_plot(list(process_res = safe_procres(fit)))
```

## Main function groups

- Input preparation: `sam.input()`, `format_SAM()`, `sam_build_data()`
- Fitting/comparison: `full_sam_fit()`, `compare_configs()`
- Output extraction: `rby.sam()`, `rbya.sam()`, `annotated_par_table()`
- Diagnostics/plots: `model_data_plot()`, `model_lik_plot()`, `model_retro_plot()`, `model_selectivity_plot()`, `model_resid_plot()`

## Testing

Run unit tests from the package root:

```r
testthat::test_local(".")
```
