---
title: "Create Giotto object SI villi data"
author: "Sandy Kroh"
date: "July 16, 2025"
output:
  html_document:
    toc: yes
    number_sections: yes
    fig_caption: yes
    keep_md: yes
  pdf_document:
    toc: yes
params:
  Sample: Default
editor_options:
  chunk_output_type: inline
---



## Libraries


``` r
if (!requireNamespace("Giotto", quietly = TRUE))
  devtools::install_github("drieslab/Giotto@suite")
library(Giotto)
```

```
## Loading required package: GiottoClass
```

```
## Newer devel version of GiottoClass on GitHub: 0.4.9
```

```
## Giotto Suite 4.2.2
```

``` r
library(ggplot2)
```

```
## Warning: package 'ggplot2' was built under R version 4.4.3
```

``` r
# library(ggpubr)
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:GiottoClass':
## 
##     intersect, union
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

## Directories


``` r
output_dir <- here::here("1_spatial_setup", "import_Giotto_Villi_files")
dir.create(output_dir)
```

```
## Warning in dir.create(output_dir):
## 'D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\1_spatial_setup\import_Giotto_Villi_files'
## already exists
```

# Configuring Giotto package and raw data

Checking if miniconda and python is in place


``` r
# configure Giotto
genv_exists = checkGiottoEnvironment()
```

```
## Giotto can access environment found at:
##  'C:\Users\NieHau\AppData\Local\r-miniconda\envs\giotto_env/python.exe'
```

```
##  If this is the wrong environment, try specifying `envname` param
##  or set option "giotto.py_path" with the desired envname or path
```

``` r
print(genv_exists)
```

```
## [1] TRUE
```

``` r
if(!genv_exists){
  # The following command need only be run once to install the Giotto environment.
  installGiottoEnvironment()
}
```

## Import raw data


``` r
datax <- read.csv(paste0(here::here(), "/data/SO_arcsinh_si_imputed_Villi.csv"))

# inspect cell types
unique(datax$CellType)
```

```
##  [1] "Epithelia I"               "Epithelia II"             
##  [3] "Fibroblasts"               "Blood vessels"            
##  [5] "Lymphatics"                "Myeloid cells"            
##  [7] "B cells"                   "Plasma cells/Plasmablasts"
##  [9] "ILC2s"                     "CD8+ CD3- IEL"            
## [11] "Unresolved"                "NK cells/ILC1s/ILC3s"     
## [13] "T helper cells"            "T cytotox. cells"
```


``` r
samples <- unique(datax$Dataset)
gio_list <- NULL
for(samp in samples){

  # cur_datax
  cur_datax <- datax[datax$Dataset == samp,]
  
  # sample and image folder
  image_folder <- unique(cur_datax$Dataset)
  print(image_folder)

  # feature data
  cur_data <- cur_datax[,15:length(colnames(cur_datax))]
  rownames(cur_data) <- cur_datax$CellID
  cur_data <- t(cur_data)

  # metadata
  cur_metadata <- cur_datax[,c(2:7, 10:12)]
  rownames(cur_metadata) <- cur_datax$CellID

  # coordinates
  cur_coords <- cur_datax[,c("Loc_X", "Loc_Y")]
  rownames(cur_coords) <- cur_datax$CellID

  # main giotto object
  gio_list[[image_folder]] <- createGiottoObject(expression = cur_data, spatial_locs = cur_coords)

  # metadata
  gio_list[[image_folder]] <- addCellMetadata(gobject =  gio_list[[image_folder]], new_metadata = cur_metadata)

  # make delaunay graph
  gio_list[[image_folder]] <- createSpatialDelaunayNetwork(gio_list[[image_folder]])
}
```

```
## [1] "CTRL_FOV1_20210706"
```

```
## a giotto python environment was found
```

```
## Using python path:
##  "C:/Users/NieHau/AppData/Local/r-miniconda/envs/giotto_env/python.exe"
```

```
## Warning: Some of Giotto's expected python module(s) were not found:
## leidenalg
## (This is fine if python-based functions are not needed)
## 
## ** Python path used:
## "C:/Users/NieHau/AppData/Local/r-miniconda/envs/giotto_env/python.exe"
```

```
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "CTRL_FOV1_20210730"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
```

```
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "CTRL_FOV2_20210706"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "CTRL_FOV2_20210730"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "CTRL_FOV2_20210810"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "CTRL_FOV3_20210706"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "CTRL_FOV3_20210709"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "CTRL_FOV3_20210730"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "D1_FOV1_20211025"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "D1_FOV2_20211025"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "D1_FOV2_20220505"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "D1_FOV3_20211025"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "D1_FOV3_20220505"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "D3_FOV1_20210701"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "D3_FOV1_20210806"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "D3_FOV2_20210625"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "D3_FOV2_20210701"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "D3_FOV2_20210806"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "D3_FOV3_20210625"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "D3_FOV3_20210701"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "D3_FOV3_20210806"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

## Visualization

Visualization of Giotto objects


``` r
for (variable in 1:length(samples)) {
  print(spatPlot2D(gio_list[[variable]], cell_color = "CellType", point_size = 1)+
          ggtitle(names(gio_list[variable])))
}
```

![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-1.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-2.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-3.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-4.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-5.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-6.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-7.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-8.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-9.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-10.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-11.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-12.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-13.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-14.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-15.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-16.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-17.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-18.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-19.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-20.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-21.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-22.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-23.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-24.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-25.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-26.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-27.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-28.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-29.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-30.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-31.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-32.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-33.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-34.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-35.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-36.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-37.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-38.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-39.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-40.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-41.png)<!-- -->![](import_Giotto_Villi_files/figure-html/unnamed-chunk-6-42.png)<!-- -->

# Save data

## Session info


``` r
saveRDS(gio_list, file = paste0(here::here(), "/data/Giotto_data_Villi.rds"))
save.image(paste0(output_dir, "/environment.RData"))
sessionInfo()
```

```
## R version 4.4.2 (2024-10-31 ucrt)
## Platform: x86_64-w64-mingw32/x64
## Running under: Windows 10 x64 (build 19045)
## 
## Matrix products: default
## 
## 
## locale:
## [1] LC_COLLATE=English_Germany.utf8  LC_CTYPE=English_Germany.utf8   
## [3] LC_MONETARY=English_Germany.utf8 LC_NUMERIC=C                    
## [5] LC_TIME=English_Germany.utf8    
## 
## time zone: Europe/Berlin
## tzcode source: internal
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] dplyr_1.1.4       ggplot2_3.5.2     Giotto_4.2.2      GiottoClass_0.4.8
## 
## loaded via a namespace (and not attached):
##   [1] colorRamp2_0.1.0            deldir_2.0-4               
##   [3] gridExtra_2.3               rlang_1.1.5                
##   [5] magrittr_2.0.3              GiottoUtils_0.2.5          
##   [7] matrixStats_1.5.0           compiler_4.4.2             
##   [9] png_0.1-8                   vctrs_0.6.5                
##  [11] pkgconfig_2.0.3             SpatialExperiment_1.16.0   
##  [13] crayon_1.5.3                fastmap_1.2.0              
##  [15] backports_1.5.0             magick_2.8.7               
##  [17] XVector_0.46.0              labeling_0.4.3             
##  [19] ggraph_2.2.1                rmarkdown_2.29             
##  [21] UCSC.utils_1.2.0            purrr_1.0.4                
##  [23] xfun_0.51                   bluster_1.16.0             
##  [25] zlibbioc_1.52.0             cachem_1.1.0               
##  [27] GenomeInfoDb_1.42.3         jsonlite_1.9.1             
##  [29] DelayedArray_0.32.0         BiocParallel_1.40.2        
##  [31] tweenr_2.0.3                terra_1.8-54               
##  [33] parallel_4.4.2              cluster_2.1.6              
##  [35] R6_2.6.1                    bslib_0.9.0                
##  [37] RColorBrewer_1.1-3          reticulate_1.42.0          
##  [39] GenomicRanges_1.58.0        jquerylib_0.1.4            
##  [41] scattermore_1.2             Rcpp_1.0.14                
##  [43] SummarizedExperiment_1.36.0 knitr_1.50                 
##  [45] IRanges_2.40.1              Matrix_1.7-1               
##  [47] igraph_2.1.4                tidyselect_1.2.1           
##  [49] rstudioapi_0.17.1           abind_1.4-8                
##  [51] yaml_2.3.10                 viridis_0.6.5              
##  [53] codetools_0.2-20            lattice_0.22-6             
##  [55] tibble_3.2.1                Biobase_2.66.0             
##  [57] withr_3.0.2                 evaluate_1.0.4             
##  [59] polyclip_1.10-7             pillar_1.10.2              
##  [61] MatrixGenerics_1.18.1       checkmate_2.3.2            
##  [63] stats4_4.4.2                plotly_4.11.0              
##  [65] generics_0.1.4              rprojroot_2.0.4            
##  [67] S4Vectors_0.44.0            scales_1.4.0               
##  [69] gtools_3.9.5                glue_1.8.0                 
##  [71] lazyeval_0.2.2              tools_4.4.2                
##  [73] GiottoVisuals_0.2.12        BiocNeighbors_2.0.1        
##  [75] data.table_1.17.0           graphlayouts_1.2.2         
##  [77] tidygraph_1.3.1             cowplot_1.1.3              
##  [79] grid_4.4.2                  tidyr_1.3.1                
##  [81] colorspace_2.1-1            SingleCellExperiment_1.28.1
##  [83] GenomeInfoDbData_1.2.13     ggforce_0.5.0              
##  [85] cli_3.6.3                   rappdirs_0.3.3             
##  [87] S4Arrays_1.6.0              viridisLite_0.4.2          
##  [89] gtable_0.3.6                sass_0.4.10                
##  [91] digest_0.6.37               BiocGenerics_0.52.0        
##  [93] SparseArray_1.6.2           ggrepel_0.9.6              
##  [95] rjson_0.2.23                htmlwidgets_1.6.4          
##  [97] farver_2.1.2                memoise_2.0.1              
##  [99] htmltools_0.5.8.1           lifecycle_1.0.4            
## [101] httr_1.4.7                  here_1.0.1                 
## [103] MASS_7.3-61
```
