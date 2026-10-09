# Aggregate data to a fixed level

This function aggregates a data set to a fixed level of a taxonomy. For
example, if aggregating to the Phylum level, it would sum up the columns
of x corresponding to OTUs/ASVs in that Phylum. It also adjusts the tree
information (in particular A, tree, and tax) for this new smaller
problem.

## Usage

``` r
aggregate_to_level(x, y, A, tax, level = 7, collapse = FALSE)
```

## Arguments

- x:

  n by p data matrix of counts, where n is the number of samples, p is
  the number of OTUs/ASVs.

- y:

  response vector of length n

- A:

  p by (t_size-1) binary matrix giving tree structure (t_size is the
  total number of nodes and the -1 is because we do not include the
  root). This is created by the function
  [`phylo_to_A`](https://viettran.de/trac/reference/phylo_to_A.md).

- tax:

  a tax table, which is a p by number-of-levels matrix, where cell ij
  gives the name of the level j ancestor of OTU/ASV i.

- level:

  which level are we aggregating to?

- collapse:

  see
  [`tax_table_to_phylo`](https://viettran.de/trac/reference/tax_table_to_phylo.md)
