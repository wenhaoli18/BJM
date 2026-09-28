# Combine per-marker transposed design matrices into one block-diagonal matrix

Shared helper: replicates the `cbind`/`rbind` block-diagonal assembly
duplicated inline in `conditionalYT.R`/ `conditionalYDT.R`, so the
copula prediction functions can reuse it rather than re-copy it a third
time.

## Usage

``` r
build_block_diagonal_matrix(matrix_list)
```
