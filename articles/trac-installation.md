# Installation

## Installation

trac internally uses the c-lasso solver (Simpson et al. 2021)
implemented in Python and reticulate to call c-lasso in R.

``` r

# install.packages("reticulate")
library(reticulate)
```

## Python

### Install c-lasso within R with reticulate

#### Install Miniconda

If you are not sure if you have python 3.6 (or later) installed use the
following code. If there is no installation of Python you will be asked
to install Miniconda and guided through the console.

``` r

py_config()
```

If know you do not have python installed install Miniconda and Python
3.6.3 with:

``` r

install_miniconda()
```

Alternatively consider following the instruction of the [Miniconda
installation](https://docs.conda.io/en/latest/miniconda.html) page.

#### Install c-lasso

First we need to install the dependencies of c-lasso, namely numpy
(Harris et al. 2020), scipy (Virtanen et al. 2020), matplotlib (Hunter
2007) and pandas (McKinney 2010). Afterwards we can install c-lasso.

``` r

conda_install(packages = c("numpy", "scipy", "matplotlib", "pandas"))
conda_install(packages = "c-lasso", pip = TRUE)
```

### Install c-lasso within the terminal

Install numpy, scipy, matplotlib, pandas and c-lasso within your virtual
/ conda environment. (Note: Depending on the operating system one has to
use commas or not)

``` bash
pip install numpy, scipy, matplotlib, pandas, c-lasso
```

Also consider pip3 if it does not work. Before loading trac specify the
the python path with

``` r

use_python(path)
```

## Install trac

``` r

# if devtools is not installed
# install.packages("devtools")
devtools::install_github("jacobbien/trac")
```

## References

Harris, Charles R, K Jarrod Millman, Stéfan J van der Walt, et al. 2020.
“Array Programming with NumPy.” *Nature* 585 (7825): 357–62.

Hunter, John D. 2007. “Matplotlib: A 2D Graphics Environment.” *IEEE
Annals of the History of Computing* 9 (03): 90–95.

McKinney, Wes. 2010. “Data Structures for Statistical Computing in
Python.” In *Proceedings of the 9th Python in Science Conference*,
edited by Stéfan van der Walt and Jarrod Millman.

Simpson, Léo, Patrick L. Combettes, and Christian L. Müller. 2021.
“C-Lasso - a Python Package for Constrained Sparse and Robust Regression
and Classification.” *Journal of Open Source Software* 6 (57): 2844.
<https://doi.org/10.21105/joss.02844>.

Virtanen, Pauli, Ralf Gommers, Travis E Oliphant, et al. 2020. “SciPy
1.0: Fundamental Algorithms for Scientific Computing in Python.” *Nature
Methods* 17 (3): 261–72.
