# Convert a tax table to a phylo object

This is mostly the function
[`as.phylo.formula`](https://rdrr.io/pkg/ape/man/as.phylo.formula.html)
from the ape package but modified so that interior nodes get labeled.
Note: assumes that there are NOT single quotes around the labels.

## Usage

``` r
tax_table_to_phylo(x, data = parent.frame(), collapse = TRUE, ...)
```

## Arguments

- x:

  a right-side formula describing the taxonomic relationship:
  `~C1/C2/.../Cn`.

- data:

  the data.frame where to look for the variables (default to user's
  workspace).

- collapse:

  a logical value specifying whether to collapse single nodes in the
  returned tree (see details).

- ...:

  further arguments to be passed from other methods.

## Author

Julien Dutheil <dutheil@evolbio.mpg.de> and Eric Marcon, modified by
Jacob Bien
