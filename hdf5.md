hdf5
================
2024-02-23

``` r
library(HDF5Array)
library(hdf5r)
library(scRNAseq)
library(tidyverse)
library(tictoc)
library(withr)
library(bench)
library(gt)
library(glue)
```

``` r
sce <- ZilionisLungData()
```

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

``` r
counts <- sce@assays@data$counts[,]
```

# manual hdf5

There are couple of R implementations of hdf5 for single cell gene
expression data:

But these are structured for single cell experiment / seurat / anndata
objects, which at minimum contain counts matrix, but can also include
multiple assays, cell metadata, projections, and other information.

Here I make an attempt to create my own hdf5 structure for taking a dgc
matrix as input.

# saving dgc matrix manually with hdf5r

Here I litter the function with tictoc::tic and tictoc::toc to
understand how long overall and each step takes.

The strategy here is to write actually a dense representation of the
sparase matrix.

The alternative would be to write a sparse representation which would be
faster to write (less data), and then convert to dense on the read in
endpoints. But here the goal is to be able to read fast - the writing
can be slow, so we but the processing (sparse -\> dense) on the write
end.

``` r
file_h5 <- here::here('data','counts.h5')
if (file.exists(file_h5)) file.remove(file_h5)
```

    ## [1] TRUE

``` r
write_dgc_to_h5 <- function(dgc, file, chunk_size = 500) {
  
  # open h5 connection and groups
  h5 <- H5File$new(file, mode = "w")
  on.exit(h5$close_all())
  h5_grp <- h5$create_group("grp")
  h5_grp_data <- h5_grp$create_dataset(
    "data",  
    dtype = h5types$H5T_NATIVE_FLOAT,
    space = H5S$new("simple", dims = dim(dgc), maxdims = dim(dgc)),
    chunk_dims = c(1, ncol(dgc))
  )
  
  tic('total')
  for (i in 1:floor((nrow(dgc) - 8)/chunk_size)) {
    
    index_start <- ((i - 1) * chunk_size) + 1
    index_end <- i * chunk_size
    tic(glue::glue('loop {i}, rows {index_start}:{index_end}'))
    
    h5_grp_data[index_start:index_end, ] <- dgc[index_start:index_end,] |> 
      as.matrix()
    toc(log = TRUE)
  }
  
  # final group
  index_start <- i*chunk_size + 1
  index_end <- nrow(dgc)
  
  tic(glue::glue('final loop, rows {index_start}:{index_end}'))
  h5_grp_data[index_start:index_end, ] <- as.matrix(dgc[index_start:index_end, ])
  toc(log = TRUE)
  
  # add rownames and colnames
  h5_grp[['rownames']] <- rownames(dgc)
  h5_grp[['colnames']] <- colnames(dgc)
  toc()
}
write_dgc_to_h5(counts, file_h5, chunk_size = 1000)
```

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 1, rows 1:1000: 7.5 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 2, rows 1001:2000: 7.06 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 3, rows 2001:3000: 7.08 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 4, rows 3001:4000: 7.14 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 5, rows 4001:5000: 7 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 6, rows 5001:6000: 6.36 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 7, rows 6001:7000: 7 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 8, rows 7001:8000: 6.61 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 9, rows 8001:9000: 6.62 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 10, rows 9001:10000: 6.33 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 11, rows 10001:11000: 6.24 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 12, rows 11001:12000: 7.2 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 13, rows 12001:13000: 6.48 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 14, rows 13001:14000: 7.32 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 15, rows 14001:15000: 6.73 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 16, rows 15001:16000: 6.92 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 17, rows 16001:17000: 6.68 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 18, rows 17001:18000: 7.34 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 19, rows 18001:19000: 6.75 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 20, rows 19001:20000: 6.83 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 21, rows 20001:21000: 6.79 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 22, rows 21001:22000: 7.33 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 23, rows 22001:23000: 6.46 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 24, rows 23001:24000: 6.43 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 25, rows 24001:25000: 6.64 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 26, rows 25001:26000: 6.82 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 27, rows 26001:27000: 7.17 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 28, rows 27001:28000: 11.36 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 29, rows 28001:29000: 7.21 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 30, rows 29001:30000: 7.84 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 31, rows 30001:31000: 6.2 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 32, rows 31001:32000: 7.52 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 33, rows 32001:33000: 11.05 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 34, rows 33001:34000: 8.48 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 35, rows 34001:35000: 6.8 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 36, rows 35001:36000: 7.79 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 37, rows 36001:37000: 6.43 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 38, rows 37001:38000: 8.31 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 39, rows 38001:39000: 6.53 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 40, rows 39001:40000: 6.62 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.3 GiB

    ## loop 41, rows 40001:41000: 13.19 sec elapsed

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 1.1 GiB

    ## final loop, rows 41001:41861: 7.74 sec elapsed
    ## total: 308.17 sec elapsed

# read 1 gene at a time

``` r
read_gene <- function(gene, file_h5) {
  #tic('whole thing')
  stopifnot(file.exists(file_h5))
  # open connections
  #tic('open')
  h5 <- H5File$new(
    file_h5, 
    mode = "r")
  
  #on.exit(h5$close_all())
  
  h5_data <- h5[['grp']][['data']]
  h5_rownames <- h5[['grp']][['rownames']]
  h5_colnames <- h5[['grp']][['colnames']]
  
  #toc()
  #tic('read gene')
  ind <- str_which(h5_rownames[], gene)
  #ind <- 18871
  #print(ind)
  gene <- h5_data[ind,]
  #toc()
  
  #tic('set name')
  gene <- setNames(gene, nm = h5_colnames[])
  #toc()
  
  #tic('close')
  h5$close_all()
  #toc()
  #toc()
  return(gene) 
}


gene <- read_gene('AC006486.2', file_h5)
```

# HDF5array

Bioconductor has the HDF5array R package that supports writing / loading
dense and sparse matrices from .h5 files.

Let’s see how this compares to my manual implementation.

First write to disk using `HDF5Array::writeHDF5Array`

``` r
file_h5array <- here::here('data', 'HDF5array.h5')
if (file.exists(file_h5array)) file.remove(file_h5array)
```

    ## [1] TRUE

``` r
tic();HDF5Array::writeHDF5Array(
  DelayedArray(counts), 
  file_h5array, as.sparse = FALSE, name = 'full', with.dimnames = TRUE);toc()
```

    ## <41861 x 173954> HDF5Matrix object of type "double":
    ##           bcIIOD bcHTNA bcDLAV ... bcELDH bcFGGM
    ##   5S_rRNA      0      0      0   .      0      0
    ## 5_8S_rRNA      0      0      0   .      0      0
    ##       7SK      0      0      0   .      0      0
    ##      A1BG      0      0      0   .      0      0
    ##  A1BG-AS1      0      0      0   .      0      0
    ##       ...      .      .      .   .      .      .
    ##   snoZ278      0      0      0   .      0      0
    ##    snoZ40      0      0      0   .      0      0
    ##     snoZ6      0      0      0   .      0      0
    ##  snosnR66      0      0      0   .      0      0
    ##    uc_338      0      0      0   .      0      0

    ## 443.03 sec elapsed

``` r
h5ls(file_h5array)
```

    ##             group           name       otype dclass            dim
    ## 0               / .full_dimnames   H5I_GROUP                      
    ## 1 /.full_dimnames              1 H5I_DATASET STRING          41861
    ## 2 /.full_dimnames              2 H5I_DATASET STRING         173954
    ## 3               /           full H5I_DATASET  FLOAT 41861 x 173954

``` r
hf5 <- HDF5Array(file_h5array,  name = 'full', as.sparse = TRUE)
hf5['AC006486.2',] |>  head()
```

    ## bcIIOD bcHTNA bcDLAV bcHNVA bcALZN bcEIYJ 
    ##      0      0      0      0      0      0

``` r
showtree(hf5)
```

    ## 41861x173954 double, sparse: HDF5Matrix object
    ## └─ 41861x173954 double, sparse: [seed] HDF5ArraySeed object

``` r
showtree(hf5[2,])
```

    ##  double: [seed] numeric object

``` r
class(hf5)
```

    ## [1] "HDF5Matrix"
    ## attr(,"package")
    ## [1] "HDF5Array"

``` r
is(hf5, 'DelayedMatrix')
```

    ## [1] TRUE

# compare HDF5array vs diy implementation

Compare performance

``` r
bench_read <- bench::mark(
  hf5['AC006486.2',],
  read_gene('AC006486.2', file_h5),
  counts['AC006486.2',]
)  

summary(bench_read) |> select(-memory, -result, -time, -gc) |>  gt()
```

<div id="rndavfqsrh" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#rndavfqsrh table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#rndavfqsrh thead, #rndavfqsrh tbody, #rndavfqsrh tfoot, #rndavfqsrh tr, #rndavfqsrh td, #rndavfqsrh th {
  border-style: none;
}
&#10;#rndavfqsrh p {
  margin: 0;
  padding: 0;
}
&#10;#rndavfqsrh .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#rndavfqsrh .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#rndavfqsrh .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#rndavfqsrh .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#rndavfqsrh .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#rndavfqsrh .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#rndavfqsrh .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#rndavfqsrh .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#rndavfqsrh .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#rndavfqsrh .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#rndavfqsrh .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#rndavfqsrh .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#rndavfqsrh .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#rndavfqsrh .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#rndavfqsrh .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#rndavfqsrh .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#rndavfqsrh .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#rndavfqsrh .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#rndavfqsrh .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#rndavfqsrh .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#rndavfqsrh .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#rndavfqsrh .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#rndavfqsrh .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#rndavfqsrh .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#rndavfqsrh .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#rndavfqsrh .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#rndavfqsrh .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#rndavfqsrh .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#rndavfqsrh .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#rndavfqsrh .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#rndavfqsrh .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#rndavfqsrh .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#rndavfqsrh .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#rndavfqsrh .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#rndavfqsrh .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#rndavfqsrh .gt_left {
  text-align: left;
}
&#10;#rndavfqsrh .gt_center {
  text-align: center;
}
&#10;#rndavfqsrh .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#rndavfqsrh .gt_font_normal {
  font-weight: normal;
}
&#10;#rndavfqsrh .gt_font_bold {
  font-weight: bold;
}
&#10;#rndavfqsrh .gt_font_italic {
  font-style: italic;
}
&#10;#rndavfqsrh .gt_super {
  font-size: 65%;
}
&#10;#rndavfqsrh .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#rndavfqsrh .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#rndavfqsrh .gt_indent_1 {
  text-indent: 5px;
}
&#10;#rndavfqsrh .gt_indent_2 {
  text-indent: 10px;
}
&#10;#rndavfqsrh .gt_indent_3 {
  text-indent: 15px;
}
&#10;#rndavfqsrh .gt_indent_4 {
  text-indent: 20px;
}
&#10;#rndavfqsrh .gt_indent_5 {
  text-indent: 25px;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="expression">expression</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="min">min</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="median">median</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="itr/sec">itr/sec</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="mem_alloc">mem_alloc</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="gc/sec">gc/sec</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="n_itr">n_itr</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="n_gc">n_gc</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="total_time">total_time</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="expression" class="gt_row gt_center">[, hf5, AC006486.2, </td>
<td headers="min" class="gt_row gt_center">1.84s</td>
<td headers="median" class="gt_row gt_center">1.84s</td>
<td headers="itr/sec" class="gt_row gt_right">0.5441028</td>
<td headers="mem_alloc" class="gt_row gt_center">10.11MB</td>
<td headers="gc/sec" class="gt_row gt_right">0</td>
<td headers="n_itr" class="gt_row gt_right">1</td>
<td headers="n_gc" class="gt_row gt_right">0</td>
<td headers="total_time" class="gt_row gt_center">1.84s</td></tr>
    <tr><td headers="expression" class="gt_row gt_center">read_gene, AC006486.2, file_h5</td>
<td headers="min" class="gt_row gt_center">769.14ms</td>
<td headers="median" class="gt_row gt_center">769.14ms</td>
<td headers="itr/sec" class="gt_row gt_right">1.3001609</td>
<td headers="mem_alloc" class="gt_row gt_center">4.97MB</td>
<td headers="gc/sec" class="gt_row gt_right">0</td>
<td headers="n_itr" class="gt_row gt_right">1</td>
<td headers="n_gc" class="gt_row gt_right">0</td>
<td headers="total_time" class="gt_row gt_center">769.14ms</td></tr>
    <tr><td headers="expression" class="gt_row gt_center">[, counts, AC006486.2, </td>
<td headers="min" class="gt_row gt_center">145.37ms</td>
<td headers="median" class="gt_row gt_center">152.82ms</td>
<td headers="itr/sec" class="gt_row gt_right">6.5416653</td>
<td headers="mem_alloc" class="gt_row gt_center">3.93MB</td>
<td headers="gc/sec" class="gt_row gt_right">0</td>
<td headers="n_itr" class="gt_row gt_right">4</td>
<td headers="n_gc" class="gt_row gt_right">0</td>
<td headers="total_time" class="gt_row gt_center">611.47ms</td></tr>
  </tbody>
  &#10;  
</table>
</div>

``` r
bench_read |>  autoplot(type = 'jitter')
```

![](hdf5_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

HDF5array slower than my manual method

Fastest is holding in memory, but not by that much.

# singlecellexperiment

``` r
sce_h5 <- SingleCellExperiment(assays = list(counts = hf5))
object.size(sce_h5) |>  print(units = 'auto')
```

    ## 9.1 Mb

``` r
object.size(counts) |>  print(units = 'auto')
```

    ## 624.4 Mb
