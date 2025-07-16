---
title: "Import VoltRon SI Villi"
author: "Sandy Kroh"
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
  devtools::install_github("Artur-man/VoltRon")
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

## Directories


``` r
output_dir <- here::here("1_spatial_setup", "import_VoltRon_Villi_files")
dir.create(output_dir)
```

```
## Warning in dir.create(output_dir):
## 'D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\1_spatial_setup\import_VoltRon_Villi_files'
## already exists
```

# Import Data

The main loop to go over each FOV and create VoltRon assays


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
image_folders_path <- paste0(here::here(), "/data/MELC_data/SI_Villi/")
samples <- unique(datax$Dataset)

# do this for all samples
vr_list <- NULL
for(samp in samples){
  # current data fram cur_datax
  cur_datax <- datax[datax$Dataset == samp,]
  rownames(cur_datax) <- cur_datax$CellID
  
  # current feature data frame
  cur_data <- cur_datax[,15:length(colnames(cur_datax))]
  rownames(cur_data) <- cur_datax$CellID

  # sample and image folder
  image_folder <- unique(cur_datax$Dataset)
  print(image_folder)

  # metadata
  cur_metadata <- cur_datax[,c(2:7, 10:12)]
  rownames(cur_metadata) <- cur_datax$CellID

  # coordinates
  cur_coords <- as.matrix(cur_datax[,c("Loc_X", "Loc_Y")])
  rownames(cur_coords) <- cur_datax$CellID

  # make voltron objects
  file_names <- list.files(paste0(image_folders_path, image_folder))
  file_names <- gsub(".png", "", file_names)
  image_names <- c("DAPI", file_names)
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
## [1] "CTRL_FOV1_20210706"
## [1] "CTRL_FOV1_20210730"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\SI_Villi\CTRL_FOV1_20210730\DAPI.png': No such file or directory @ error/blob.c/OpenBlob/2964
## [1] "CTRL_FOV2_20210706"
## [1] "CTRL_FOV2_20210730"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\SI_Villi\CTRL_FOV2_20210730\DAPI.png': No such file or directory @ error/blob.c/OpenBlob/2964
## [1] "CTRL_FOV2_20210810"
## [1] "CTRL_FOV3_20210706"
## [1] "CTRL_FOV3_20210709"
## [1] "CTRL_FOV3_20210730"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\SI_Villi\CTRL_FOV3_20210730\DAPI.png': No such file or directory @ error/blob.c/OpenBlob/2964
## [1] "D1_FOV1_20211025"
## [1] "D1_FOV2_20211025"
## [1] "D1_FOV2_20220505"
## [1] "D1_FOV3_20211025"
## [1] "D1_FOV3_20220505"
## [1] "D3_FOV1_20210701"
## [1] "D3_FOV1_20210806"
## [1] "D3_FOV2_20210625"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\SI_Villi\D3_FOV2_20210625\DAPI.png': No such file or directory @ error/blob.c/OpenBlob/2964
## [1] "D3_FOV2_20210701"
## [1] "D3_FOV2_20210806"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\SI_Villi\D3_FOV2_20210806\DAPI.png': No such file or directory @ error/blob.c/OpenBlob/2964
## [1] "D3_FOV3_20210625"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\SI_Villi\D3_FOV3_20210625\DAPI.png': No such file or directory @ error/blob.c/OpenBlob/2964
## [1] "D3_FOV3_20210701"
## [1] "D3_FOV3_20210806"
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
## Assay1   MELC Section1 CTRL_FOV1_20210706
## Assay2   MELC Section1 CTRL_FOV1_20210730
## Assay3   MELC Section1 CTRL_FOV2_20210706
## Assay4   MELC Section1 CTRL_FOV2_20210730
## Assay5   MELC Section1 CTRL_FOV2_20210810
## Assay6   MELC Section1 CTRL_FOV3_20210706
## Assay7   MELC Section1 CTRL_FOV3_20210709
## Assay8   MELC Section1 CTRL_FOV3_20210730
## Assay9   MELC Section1   D1_FOV1_20211025
## Assay10  MELC Section1   D1_FOV2_20211025
## Assay11  MELC Section1   D1_FOV2_20220505
## Assay12  MELC Section1   D1_FOV3_20211025
## Assay13  MELC Section1   D1_FOV3_20220505
## Assay14  MELC Section1   D3_FOV1_20210701
## Assay15  MELC Section1   D3_FOV1_20210806
## Assay16  MELC Section1   D3_FOV2_20210625
## Assay17  MELC Section1   D3_FOV2_20210701
## Assay18  MELC Section1   D3_FOV2_20210806
## Assay19  MELC Section1   D3_FOV3_20210625
## Assay20  MELC Section1   D3_FOV3_20210701
## Assay21  MELC Section1   D3_FOV3_20210806
```

# Visualization


``` r
# visualize all FOVs, default image is DAPI
vrSpatialPlot(vr_merged, group.by = "CellType", alpha = 0.2, ncol = 3, 
                  legend.loc = "bottom", legend.pt.size = 1, pt.size = 1) +
    theme(plot.margin = margin(0.1,0.1,0.1,0.1, "cm"), 
          title =element_text(size=4), 
          axis.title = element_blank(), 
          axis.ticks = element_blank(), 
          axis.text = element_blank(),
          legend.text = element_text(size = 8), 
          legend.title = element_text(size = 8))+
    ggplot2::guides(color=guide_legend(override.aes = list(size=20), ncol=4))
```

<img src="import_VoltRon_Villi_files/figure-html/unnamed-chunk-8-1.png" width="100%" style="display: block; margin: auto;" />

# Save data

## Session info


``` r
saveRDS(vr_list, file = paste0(here::here(), "/data/VoltRon_data_Villi.rds"))
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
## [1] ggpubr_0.6.0       ggplot2_3.5.2      dplyr_1.1.4        readr_2.1.5       
## [5] magick_2.8.7       VoltRon_0.2.0      Seurat_5.2.1       SeuratObject_5.1.0
## [9] sp_2.2-0          
## 
## loaded via a namespace (and not attached):
##   [1] RcppAnnoy_0.0.22       splines_4.4.2          later_1.4.1           
##   [4] bitops_1.0-9           tibble_3.2.1           polyclip_1.10-7       
##   [7] fastDummies_1.7.5      lifecycle_1.0.4        rstatix_0.7.2         
##  [10] rprojroot_2.0.4        globals_0.18.0         lattice_0.22-6        
##  [13] MASS_7.3-61            backports_1.5.0        magrittr_2.0.3        
##  [16] plotly_4.11.0          sass_0.4.10            rmarkdown_2.29        
##  [19] jquerylib_0.1.4        yaml_2.3.10            httpuv_1.6.15         
##  [22] sctransform_0.4.1      spam_2.11-1            spatstat.sparse_3.1-0 
##  [25] reticulate_1.42.0      cowplot_1.1.3          pbapply_1.7-2         
##  [28] RColorBrewer_1.1-3     abind_1.4-8            Rtsne_0.17            
##  [31] purrr_1.0.4            BiocGenerics_0.52.0    RCurl_1.98-1.17       
##  [34] rgl_1.3.18             IRanges_2.40.1         S4Vectors_0.44.0      
##  [37] ggrepel_0.9.6          irlba_2.3.5.1          listenv_0.9.1         
##  [40] spatstat.utils_3.1-3   goftest_1.2-3          RSpectra_0.16-2       
##  [43] spatstat.random_3.3-3  fitdistrplus_1.2-2     parallelly_1.45.0     
##  [46] Rvcg_0.25              codetools_0.2-20       tidyselect_1.2.1      
##  [49] farver_2.1.2           matrixStats_1.5.0      stats4_4.4.2          
##  [52] base64enc_0.1-3        spatstat.explore_3.4-2 jsonlite_1.9.1        
##  [55] progressr_0.15.1       Formula_1.2-5          ggridges_0.5.6        
##  [58] survival_3.7-0         tools_4.4.2            ica_1.0-3             
##  [61] Rcpp_1.0.14            glue_1.8.0             gridExtra_2.3         
##  [64] here_1.0.1             xfun_0.51              EBImage_4.48.0        
##  [67] withr_3.0.2            fastmap_1.2.0          shinyjs_2.1.0         
##  [70] caTools_1.18.3         digest_0.6.37          R6_2.6.1              
##  [73] mime_0.13              colorspace_2.1-1       scattermore_1.2       
##  [76] gtools_3.9.5           tensor_1.5.1           jpeg_0.1-11           
##  [79] spatstat.data_3.1-6    tidyr_1.3.1            generics_0.1.4        
##  [82] data.table_1.17.0      httr_1.4.7             htmlwidgets_1.6.4     
##  [85] S4Arrays_1.6.0         scatterplot3d_0.3-44   uwot_0.2.3            
##  [88] pkgconfig_2.0.3        gtable_0.3.6           lmtest_0.9-40         
##  [91] ids_1.0.1              htmltools_0.5.8.1      carData_3.0-5         
##  [94] dotCall64_1.2          fftwtools_0.9-11       scales_1.4.0          
##  [97] png_0.1-8              spatstat.univar_3.1-2  knitr_1.50            
## [100] rstudioapi_0.17.1      tzdb_0.4.0             reshape2_1.4.4        
## [103] rjson_0.2.23           uuid_1.2-1             nlme_3.1-166          
## [106] cachem_1.1.0           zoo_1.8-13             Polychrome_1.5.4      
## [109] stringr_1.5.1          KernSmooth_2.23-24     parallel_4.4.2        
## [112] miniUI_0.1.2           pillar_1.10.2          grid_4.4.2            
## [115] vctrs_0.6.5            colorsGen_1.0.0        RANN_2.6.2            
## [118] gplots_3.2.0           promises_1.3.2         car_3.1-3             
## [121] xtable_1.8-4           cluster_2.1.6          evaluate_1.0.4        
## [124] cli_3.6.3              locfit_1.5-9.12        compiler_4.4.2        
## [127] rlang_1.1.5            crayon_1.5.3           future.apply_1.20.0   
## [130] ggsignif_0.6.4         labeling_0.4.3         plyr_1.8.9            
## [133] stringi_1.8.4          viridisLite_0.4.2      deldir_2.0-4          
## [136] lazyeval_0.2.2         tiff_0.1-12            spatstat.geom_3.3-6   
## [139] Matrix_1.7-1           RcppHNSW_0.6.0         hms_1.1.3             
## [142] patchwork_1.3.1        future_1.58.0          shiny_1.10.0          
## [145] ROCR_1.0-11            igraph_2.1.4           broom_1.0.8           
## [148] bslib_0.9.0            RCDT_1.3.0
```
