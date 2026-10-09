# Central Park soil moisture data

This is the Central Park Soil data used in [Tree-Aggregated Predictive
Modeling of Microbiome
Data](https://www.biorxiv.org/content/10.1101/2020.09.01.277632v1.full),
with soil moisture as the response, aggregated to the family level. The
soil samples were collected in New York City's Central Park by [Ramirez
et al. (2014)](https://doi.org/10.1098/rspb.2014.1988). The counts were
summed to the family level with
[`aggregate_to_level()`](https://viettran.de/trac/reference/aggregate_to_level.md),
as in the trac paper.

## Usage

``` r
cps_moisture
```

## Format

A named list:

- y:

  Vector of n = 580 soil moisture values

- x:

  Matrix of soil 16S rRNA amplicon data, with n = 580 rows corresponding
  to soil samples and p = 1492 columns corresponding to families

- tree:

  Taxonomic tree of class `phylo`

- tax:

  A data frame containing the taxonomic information for each family

- A:

  A binary matrix encoding the tree structure with p = 1492 rows,
  corresponding to leaves, and 1634 columns, corresponding to all
  non-root nodes in the tree.

- covariates:

  A data frame with the n = 580 soil samples as rows and five numeric
  soil properties from the original sample data: `pH`, total carbon
  (`C`), total nitrogen (`N`), their ratio (`CN`) and `CO2_C`. `C`, `N`
  and `CN` are missing for 29 samples.
