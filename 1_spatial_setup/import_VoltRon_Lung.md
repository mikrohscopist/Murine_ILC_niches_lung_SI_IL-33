---
title: "import VoltRon Lung"
date: "July 16, 2025"
output:
  html_document:
    toc: yes
    toc_float: yes
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
if (!requireNamespace("Seurat", quietly = TRUE))
  install.packages("Seurat")
if (!requireNamespace("VoltRon", quietly = TRUE))
  devtools::install_github("BIMSBbioinfo/VoltRon@dev")#  devtools::install_github("Artur-man/VoltRon")
if (!requireNamespace("magick", quietly = TRUE))
  install.packages("magick")
library(Seurat)
```

```
## Loading required package: SeuratObject
```

```
## Warning: package 'SeuratObject' was built under R version 4.4.3
```

```
## Loading required package: sp
```

```
## 
## Attaching package: 'SeuratObject'
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, t
```

``` r
library(VoltRon)
```

```
## 
## Attaching package: 'VoltRon'
```

```
## The following object is masked from 'package:Seurat':
## 
##     as.Seurat
```

```
## The following object is masked from 'package:SeuratObject':
## 
##     as.Seurat
```

``` r
library(magick)
```

```
## Warning: package 'magick' was built under R version 4.4.3
```

```
## Linking to ImageMagick 6.9.12.98
## Enabled features: cairo, freetype, fftw, ghostscript, heic, lcms, pango, raw, rsvg, webp
## Disabled features: fontconfig, x11
```

``` r
library(readr)
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
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

``` r
library(ggplot2)
```

```
## Warning: package 'ggplot2' was built under R version 4.4.3
```

``` r
library(ggpubr)
```

```
## Warning: package 'ggpubr' was built under R version 4.4.3
```

## Directories


``` r
output_dir <- here::here("1_spatial_setup", "import_VoltRon_Lung_files")
dir.create(output_dir)
```

```
## Warning in dir.create(output_dir):
## 'D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\1_spatial_setup\import_VoltRon_Lung_files'
## already exists
```

# Import Data


``` r
datax <- read.csv(paste0(here::here(), "/data/MELC_data_murine_lung_CTRL_D1_D2_D3_withfolders.csv"))

# inspect cell types
unique(datax$CellType)
```

```
##  [1] "Epithelia"               "EMCN CD31 Blood vessels"
##  [3] "LYVE1 CD31 vessels"      "LYVE1 CD90 Lymphatics"  
##  [5] "Myeloid cells"           "B cells & Plasma cells" 
##  [7] "NK cells/ILC1s"          "ILC3s"                  
##  [9] "T cytotox cells"         "T helper cells"         
## [11] "ILC2s"
```

The main loop to go over each FOV and create VoltRon assays


``` r
image_folders_path <- paste0(here::here(), "/data/MELC_data/Lung/")

samples <- unique(datax$Sample)

# do this for all samples
vr_list <- NULL
for(samp in samples[-35]){
  # current data fram cur_datax
  cur_datax <- datax[datax$Sample == samp,]
  rownames(cur_datax) <- cur_datax$CellID
  
  # current feature data frame
  cur_data <- cur_datax[,2:34]
  rownames(cur_data) <- cur_datax$CellID

  # sample and image folder
  image_folder <- unique(cur_datax$Sample)
  print(image_folder)

  # metadata
  cur_metadata <- cur_datax[,38:46]
  rownames(cur_metadata) <- cur_datax$CellID

  # coordinates
  cur_coords <- as.matrix(cur_datax[,c("Loc_X", "Loc_Y")])
  rownames(cur_coords) <- cur_datax$CellID

  # make voltron objects
  file_names <- list.files(paste0(image_folders_path, image_folder))
  file_names <- gsub(".png", "", file_names)
  image_names <- file_names
  image_list <- list()
  for(img in image_names){
    try({
      image_list[[img]] <- magick::image_read(paste0(image_folders_path,  image_folder, "/", img, ".png"))
    })
  }
  vr_list[[image_folder]] <- formVoltRon(t(cur_data), metadata = cur_metadata, image = image_list, coords = cur_coords, main.assay = "MELC", sample_name = image_folder)
  
  # flip and resize images
  vr_list[[image_folder]] <- flipCoordinates(vr_list[[image_folder]])
  vr_list[[image_folder]] <- resizeImage(vr_list[[image_folder]], size = 1000)
}
```

```
## [1] "20210910_1_lu_ctrl"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210910_1_lu_ctrl\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2964
## [1] "20210914_1_lu_ctrl"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210914_1_lu_ctrl\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2964
## [1] "20210922_1_lu_ctrl"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210922_1_lu_ctrl\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2964
## [1] "20210910_2_lu_ctrl"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210910_2_lu_ctrl\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2964
## [1] "20210914_2_lu_ctrl"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210914_2_lu_ctrl\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2964
## [1] "20210922_2_lu_ctrl"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210922_2_lu_ctrl\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2964
## [1] "20210910_3_lu_ctrl"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210910_3_lu_ctrl\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2964
## [1] "20210914_3_lu_ctrl"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210914_3_lu_ctrl\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2964
## [1] "20210922_3_lu_ctrl"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210922_3_lu_ctrl\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2964
## [1] "20220311_1"
## [1] "20220316_1"
## [1] "20220321_1"
## [1] "20220311_2"
## [1] "20220316_2"
## [1] "20220321_2"
## [1] "20220311_3"
## [1] "20220316_3"
## [1] "20220321_3"
## [1] "20220325_1"
## [1] "20220421_1"
## [1] "20220502_1"
## [1] "20220325_2"
## [1] "20220421_2"
## [1] "20220502_2"
## [1] "20220325_3"
## [1] "20220421_3"
## [1] "20220502_3"
## [1] "20210902_1_lu_d3"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210902_1_lu_d3\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2964
## [1] "20210906_1_lu_d3"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210906_1_lu_d3\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2964
## [1] "20210928_1_lu_d3"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210928_1_lu_d3\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2964
## [1] "20210902_2_lu_d3"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210902_2_lu_d3\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2964
## [1] "20210906_2_lu_d3"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210906_2_lu_d3\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2964
## [1] "20210928_2_lu_d3"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210928_2_lu_d3\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2964
## [1] "20210902_3_lu_d3"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210902_3_lu_d3\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2964
## [1] "20210928_3_lu_d3"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210928_3_lu_d3\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2964
```

## Manipulate VoltRon list

Here, we merge all VoltRon objects into one object


``` r
vr_merged <- merge(vr_list[[1]], vr_list[-1])
```

```
## Merging metadata ...
```

```
## Merging blocks and layers ...
```

Getting a full list of all images available in VoltRon


``` r
vrImageNames(vr_merged)
```

```
## [1] "image_1"
```

Getting all FOVs


``` r
SampleMetadata(vr_merged)
```

```
##         Assay    Layer             Sample
## Assay1   MELC Section1 20210910_1_lu_ctrl
## Assay2   MELC Section1 20210914_1_lu_ctrl
## Assay3   MELC Section1 20210922_1_lu_ctrl
## Assay4   MELC Section1 20210910_2_lu_ctrl
## Assay5   MELC Section1 20210914_2_lu_ctrl
## Assay6   MELC Section1 20210922_2_lu_ctrl
## Assay7   MELC Section1 20210910_3_lu_ctrl
## Assay8   MELC Section1 20210914_3_lu_ctrl
## Assay9   MELC Section1 20210922_3_lu_ctrl
## Assay10  MELC Section1         20220311_1
## Assay11  MELC Section1         20220316_1
## Assay12  MELC Section1         20220321_1
## Assay13  MELC Section1         20220311_2
## Assay14  MELC Section1         20220316_2
## Assay15  MELC Section1         20220321_2
## Assay16  MELC Section1         20220311_3
## Assay17  MELC Section1         20220316_3
## Assay18  MELC Section1         20220321_3
## Assay19  MELC Section1         20220325_1
## Assay20  MELC Section1         20220421_1
## Assay21  MELC Section1         20220502_1
## Assay22  MELC Section1         20220325_2
## Assay23  MELC Section1         20220421_2
## Assay24  MELC Section1         20220502_2
## Assay25  MELC Section1         20220325_3
## Assay26  MELC Section1         20220421_3
## Assay27  MELC Section1         20220502_3
## Assay28  MELC Section1   20210902_1_lu_d3
## Assay29  MELC Section1   20210906_1_lu_d3
## Assay30  MELC Section1   20210928_1_lu_d3
## Assay31  MELC Section1   20210902_2_lu_d3
## Assay32  MELC Section1   20210906_2_lu_d3
## Assay33  MELC Section1   20210928_2_lu_d3
## Assay34  MELC Section1   20210902_3_lu_d3
## Assay35  MELC Section1   20210928_3_lu_d3
```

# Visualization

Different visualization scripts.


``` r
# visualize all FOVs, default image is DAPI
vrSpatialPlot(vr_merged, group.by = "CellType", alpha = 0.2, ncol = 3, 
                  legend.loc = "bottom", legend.text.size = 11,legend.pt.size = 3, pt.size = 1) +
    theme(plot.margin = margin(0.1,0.1,0.1,0.1, "cm"), 
          title =element_text(size=4), 
          axis.title = element_blank(), 
          axis.ticks = element_blank(), 
          axis.text = element_blank(),
          legend.text = element_text(size = 8), 
          legend.title = element_blank())
```

```
## Warning: Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
## Duplicated `override.aes` is ignored.
```

<img src="import_VoltRon_Lung_files/figure-html/unnamed-chunk-8-1.png" width="100%" style="display: block; margin: auto;" />

# Save data

## Session info


``` r
saveRDS(vr_list, file = paste0(here::here(), "/data/VoltRon_data_lung.rds"))
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
## [1] ggpubr_0.6.1       ggplot2_3.5.2      dplyr_1.1.4        readr_2.1.5       
## [5] magick_2.8.7       VoltRon_0.2.0      Seurat_5.2.1       SeuratObject_5.1.0
## [9] sp_2.2-0          
## 
## loaded via a namespace (and not attached):
##   [1] RcppAnnoy_0.0.22        splines_4.4.2           later_1.4.1            
##   [4] bitops_1.0-9            tibble_3.2.1            BPCells_0.3.0          
##   [7] polyclip_1.10-7         fastDummies_1.7.5       lifecycle_1.0.4        
##  [10] rstatix_0.7.2           rprojroot_2.0.4         globals_0.18.0         
##  [13] lattice_0.22-6          MASS_7.3-61             backports_1.5.0        
##  [16] magrittr_2.0.3          plotly_4.11.0           sass_0.4.10            
##  [19] rmarkdown_2.29          jquerylib_0.1.4         yaml_2.3.10            
##  [22] httpuv_1.6.15           sctransform_0.4.1       spam_2.11-1            
##  [25] spatstat.sparse_3.1-0   reticulate_1.42.0       cowplot_1.2.0          
##  [28] pbapply_1.7-2           RColorBrewer_1.1-3      abind_1.4-8            
##  [31] zlibbioc_1.52.0         GenomicRanges_1.58.0    Rtsne_0.17             
##  [34] purrr_1.0.4             BiocGenerics_0.52.0     RCurl_1.98-1.17        
##  [37] rgl_1.3.18              GenomeInfoDbData_1.2.13 IRanges_2.40.1         
##  [40] S4Vectors_0.44.0        ggrepel_0.9.6           irlba_2.3.5.1          
##  [43] listenv_0.9.1           spatstat.utils_3.1-3    goftest_1.2-3          
##  [46] RSpectra_0.16-2         spatstat.random_3.3-3   fitdistrplus_1.2-2     
##  [49] parallelly_1.45.0       Rvcg_0.25               codetools_0.2-20       
##  [52] DelayedArray_0.32.0     tidyselect_1.2.1        UCSC.utils_1.2.0       
##  [55] farver_2.1.2            ScaledMatrix_1.14.0     matrixStats_1.5.0      
##  [58] stats4_4.4.2            base64enc_0.1-3         spatstat.explore_3.4-2 
##  [61] jsonlite_1.9.1          progressr_0.15.1        Formula_1.2-5          
##  [64] ggridges_0.5.6          survival_3.7-0          tools_4.4.2            
##  [67] ica_1.0-3               Rcpp_1.0.14             glue_1.8.0             
##  [70] gridExtra_2.3           SparseArray_1.6.2       here_1.0.1             
##  [73] xfun_0.51               MatrixGenerics_1.18.1   GenomeInfoDb_1.42.3    
##  [76] EBImage_4.48.0          withr_3.0.2             fastmap_1.2.0          
##  [79] shinyjs_2.1.0           caTools_1.18.3          digest_0.6.37          
##  [82] rsvd_1.0.5              R6_2.6.1                mime_0.13              
##  [85] colorspace_2.1-1        scattermore_1.2         gtools_3.9.5           
##  [88] tensor_1.5.1            jpeg_0.1-11             spatstat.data_3.1-6    
##  [91] tidyr_1.3.1             generics_0.1.4          data.table_1.17.0      
##  [94] httr_1.4.7              htmlwidgets_1.6.4       S4Arrays_1.6.0         
##  [97] scatterplot3d_0.3-44    uwot_0.2.3              pkgconfig_2.0.3        
## [100] gtable_0.3.6            lmtest_0.9-40           XVector_0.46.0         
## [103] ids_1.0.1               htmltools_0.5.8.1       carData_3.0-5          
## [106] dotCall64_1.2           fftwtools_0.9-11        scales_1.4.0           
## [109] png_0.1-8               spatstat.univar_3.1-2   knitr_1.50             
## [112] rstudioapi_0.17.1       tzdb_0.4.0              rjson_0.2.23           
## [115] reshape2_1.4.4          uuid_1.2-1              nlme_3.1-166           
## [118] cachem_1.1.0            zoo_1.8-13              Polychrome_1.5.4       
## [121] stringr_1.5.1           KernSmooth_2.23-24      parallel_4.4.2         
## [124] miniUI_0.1.2            pillar_1.11.0           grid_4.4.2             
## [127] vctrs_0.6.5             colorsGen_1.0.0         RANN_2.6.2             
## [130] gplots_3.2.0            promises_1.3.2          BiocSingular_1.22.0    
## [133] car_3.1-3               beachmat_2.22.0         xtable_1.8-4           
## [136] cluster_2.1.6           evaluate_1.0.4          cli_3.6.3              
## [139] locfit_1.5-9.12         compiler_4.4.2          rlang_1.1.5            
## [142] crayon_1.5.3            future.apply_1.20.0     ggsignif_0.6.4         
## [145] labeling_0.4.3          plyr_1.8.9              stringi_1.8.4          
## [148] viridisLite_0.4.2       deldir_2.0-4            BiocParallel_1.40.2    
## [151] lazyeval_0.2.2          tiff_0.1-12             spatstat.geom_3.3-6    
## [154] Matrix_1.7-1            RcppHNSW_0.6.0          hms_1.1.3              
## [157] patchwork_1.3.1         future_1.58.0           shiny_1.11.1           
## [160] ROCR_1.0-11             igraph_2.1.4            broom_1.0.8            
## [163] bslib_0.9.0             RCDT_1.3.0
```
