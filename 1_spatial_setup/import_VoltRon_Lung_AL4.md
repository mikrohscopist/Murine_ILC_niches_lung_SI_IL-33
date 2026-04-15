---
title: "import VoltRon Lung - AL4 GATA3+ Th cells"
date: "April 15, 2026"
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
## Lade nötiges Paket: SeuratObject
```

```
## Lade nötiges Paket: sp
```

```
## 
## Attache Paket: 'SeuratObject'
```

```
## Die folgenden Objekte sind maskiert von 'package:base':
## 
##     intersect, t
```

``` r
library(VoltRon)
```

```
## 
## Attache Paket: 'VoltRon'
```

```
## Das folgende Objekt ist maskiert 'package:Seurat':
## 
##     as.Seurat
```

```
## Das folgende Objekt ist maskiert 'package:SeuratObject':
## 
##     as.Seurat
```

``` r
library(magick)
```

```
## Linking to ImageMagick 6.9.13.29
## Enabled features: cairo, freetype, fftw, ghostscript, heic, lcms, pango, raw, rsvg, webp
## Disabled features: fontconfig, x11
```

``` r
library(readr)
library(dplyr)
```

```
## 
## Attache Paket: 'dplyr'
```

```
## Die folgenden Objekte sind maskiert von 'package:stats':
## 
##     filter, lag
```

```
## Die folgenden Objekte sind maskiert von 'package:base':
## 
##     intersect, setdiff, setequal, union
```

``` r
library(ggplot2)
library(ggpubr)
```

## Directories


``` r
output_dir <- here::here("1_spatial_setup", "import_VoltRon_Lung_AL4_files")
dir.create(output_dir)
```

# Import Data


``` r
datax <- read.csv(paste0(here::here(), "/data/MELC_data_murine_lung_CTRL_D1_D2_D3_withfolders.csv"))
```

Data prep


``` r
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)

# 1. Define the target cell types and markers
target_types <- c("T helper cells", "T cytotox cells", "ILC2s", "ILC3s", "NK cells/ILC1s")
target_markers <- c("EOMES", "TBET", "GATA3eGFP", "RORgt", "CD4", "CD8a", "CD3")

# 2. Process the data
plot_data <- datax %>%
  # Filter for the specific cell types in AL3
  filter(AL3 %in% target_types) %>%
  # Select only the metadata and the markers of interest
  select(AL3, all_of(target_markers)) %>%
  # Pivot to long format: this puts all marker names in one column 
  # and all expression values in another, which is required for ridgeplots
  pivot_longer(cols = all_of(target_markers), 
               names_to = "Marker", 
               values_to = "Expression")

# 3. Create the Ridgeplot
ggplot(plot_data, aes(x = Expression, y = Marker, fill = AL3)) +
  geom_density_ridges(alpha = 0.7, scale = 1.5, rel_min_height = 0.01) +
  # You can facet by Marker or by Cell Type depending on how you want to compare
  facet_wrap(~ AL3, ncol = 5) + 
  theme_ridges() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10),
    axis.title.x = element_text(hjust = 0.5),
    axis.title.y = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Marker Expression Profiles by Cell Type (AL3)",
    x = "Expression Intensity",
    y = "Protein Marker"
  ) +
  # Optional: use a nice color palette
  scale_fill_viridis_d(option = "D")
```

```
## Picking joint bandwidth of 0.553
```

```
## Picking joint bandwidth of 0.868
```

```
## Picking joint bandwidth of 0.428
```

```
## Picking joint bandwidth of 0.429
```

```
## Picking joint bandwidth of 0.417
```

<img src="import_VoltRon_Lung_AL4_files/figure-html/unnamed-chunk-4-1.png" alt="" width="100%" style="display: block; margin: auto;" />


``` r
library(dplyr)

datax <- datax %>%
  mutate(AL4 = case_when(
    # --- T Helper Cell Logic ---
    AL3 == "T helper cells" & GATA3eGFP > 1.5 ~ "Th2",
    AL3 == "T helper cells" & RORgt > 1.5     ~ "Th17",
    AL3 == "T helper cells" & TBET > 1      ~ "Th1",
    AL3 == "T helper cells"                  ~ "Other Th",
    
    # --- T Cytotox Cell Logic ---
    AL3 == "T cytotox cells" & GATA3eGFP > 1.5 ~ "Tc2",
    AL3 == "T cytotox cells"                   ~ "Other Tc",
    
    # --- All other cells (Retain original AL3 label) ---
    TRUE ~ as.character(AL3)
  ))

# define AL2 as Celltypes and inspect cell types
datax$CellType <- datax$AL4

datax$CellType <- gsub("EMCN CD31 |LYVE1 CD90 |lasma", "", datax$CellType)
unique(datax$CellType)
```

```
##  [1] "Epithelia"          "Blood vessels"      "LYVE1 CD31 vessels"
##  [4] "Lymphatics"         "Myeloid cells"      "B cells & P cells" 
##  [7] "NK cells/ILC1s"     "ILC3s"              "Other Tc"          
## [10] "Tc2"                "Th17"               "Other Th"          
## [13] "Th2"                "Th1"                "ILC2s"
```

# VoltRon objects

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
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210910_1_lu_ctrl\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2967
## [1] "20210914_1_lu_ctrl"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210914_1_lu_ctrl\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2967
## [1] "20210922_1_lu_ctrl"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210922_1_lu_ctrl\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2967
## [1] "20210910_2_lu_ctrl"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210910_2_lu_ctrl\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2967
## [1] "20210914_2_lu_ctrl"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210914_2_lu_ctrl\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2967
## [1] "20210922_2_lu_ctrl"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210922_2_lu_ctrl\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2967
## [1] "20210910_3_lu_ctrl"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210910_3_lu_ctrl\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2967
## [1] "20210914_3_lu_ctrl"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210914_3_lu_ctrl\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2967
## [1] "20210922_3_lu_ctrl"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210922_3_lu_ctrl\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2967
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
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210902_1_lu_d3\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2967
## [1] "20210906_1_lu_d3"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210906_1_lu_d3\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2967
## [1] "20210928_1_lu_d3"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210928_1_lu_d3\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2967
## [1] "20210902_2_lu_d3"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210902_2_lu_d3\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2967
## [1] "20210906_2_lu_d3"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210906_2_lu_d3\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2967
## [1] "20210928_2_lu_d3"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210928_2_lu_d3\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2967
## [1] "20210902_3_lu_d3"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210902_3_lu_d3\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2967
## [1] "20210928_3_lu_d3"
## Error in eval(expr, envir) : 
##   Rterm.exe: UnableToOpenBlob `D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\data\MELC_data\Lung\20210928_3_lu_d3\Thumbs.db.png': No such file or directory @ error/blob.c/OpenBlob/2967
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

<img src="import_VoltRon_Lung_AL4_files/figure-html/unnamed-chunk-10-1.png" alt="" width="100%" style="display: block; margin: auto;" />

# Save data

## Session info


``` r
saveRDS(vr_list, file = paste0(output_dir, "/VoltRon_data_lung_AL4.rds"))
save.image(paste0(output_dir, "/environment.RData"))
sessionInfo()
```

```
## R version 4.5.2 (2025-10-31 ucrt)
## Platform: x86_64-w64-mingw32/x64
## Running under: Windows 11 x64 (build 26200)
## 
## Matrix products: default
##   LAPACK version 3.12.1
## 
## locale:
## [1] LC_COLLATE=German_Germany.utf8  LC_CTYPE=German_Germany.utf8   
## [3] LC_MONETARY=German_Germany.utf8 LC_NUMERIC=C                   
## [5] LC_TIME=German_Germany.utf8    
## 
## time zone: Europe/Berlin
## tzcode source: internal
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] ggridges_0.5.7     tidyr_1.3.1        ggpubr_0.6.2       ggplot2_4.0.1     
##  [5] dplyr_1.1.4        readr_2.1.6        magick_2.9.0       VoltRon_0.2.3     
##  [9] Seurat_5.3.1       SeuratObject_5.2.0 sp_2.2-0          
## 
## loaded via a namespace (and not attached):
##   [1] RcppAnnoy_0.0.22       splines_4.5.2          later_1.4.4           
##   [4] bitops_1.0-9           tibble_3.3.0           polyclip_1.10-7       
##   [7] fastDummies_1.7.5      lifecycle_1.0.5        rstatix_0.7.3         
##  [10] rprojroot_2.1.1        globals_0.19.0         lattice_0.22-7        
##  [13] MASS_7.3-65            backports_1.5.0        magrittr_2.0.4        
##  [16] plotly_4.12.0          sass_0.4.10            rmarkdown_2.30        
##  [19] jquerylib_0.1.4        yaml_2.3.10            httpuv_1.6.16         
##  [22] otel_0.2.0             sctransform_0.4.2      spam_2.11-1           
##  [25] spatstat.sparse_3.1-0  reticulate_1.44.0      cowplot_1.2.0         
##  [28] pbapply_1.7-4          RColorBrewer_1.1-3     abind_1.4-8           
##  [31] Rtsne_0.17             purrr_1.2.0            BiocGenerics_0.56.0   
##  [34] RCurl_1.98-1.17        rgl_1.3.34             IRanges_2.44.0        
##  [37] S4Vectors_0.48.0       ggrepel_0.9.6          irlba_2.3.5.1         
##  [40] listenv_0.10.0         spatstat.utils_3.2-0   goftest_1.2-3         
##  [43] RSpectra_0.16-2        spatstat.random_3.4-2  fitdistrplus_1.2-6    
##  [46] parallelly_1.45.1      Rvcg_0.25              codetools_0.2-20      
##  [49] DelayedArray_0.36.0    tidyselect_1.2.1       farver_2.1.2          
##  [52] ScaledMatrix_1.18.0    matrixStats_1.5.0      stats4_4.5.2          
##  [55] base64enc_0.1-6        spatstat.explore_3.5-3 jsonlite_2.0.0        
##  [58] progressr_0.18.0       Formula_1.2-5          survival_3.8-3        
##  [61] tools_4.5.2            ica_1.0-3              Rcpp_1.1.0            
##  [64] glue_1.8.0             gridExtra_2.3          SparseArray_1.10.1    
##  [67] here_1.0.2             xfun_0.56              MatrixGenerics_1.22.0 
##  [70] EBImage_4.52.0         withr_3.0.2            fastmap_1.2.0         
##  [73] shinyjs_2.1.1          caTools_1.18.3         digest_0.6.38         
##  [76] rsvd_1.0.5             R6_2.6.1               mime_0.13             
##  [79] colorspace_2.1-2       scattermore_1.2        gtools_3.9.5          
##  [82] tensor_1.5.1           jpeg_0.1-11            dichromat_2.0-0.1     
##  [85] spatstat.data_3.1-9    generics_0.1.4         data.table_1.17.8     
##  [88] httr_1.4.8             htmlwidgets_1.6.4      S4Arrays_1.10.0       
##  [91] scatterplot3d_0.3-45   uwot_0.2.4             pkgconfig_2.0.3       
##  [94] gtable_0.3.6           lmtest_0.9-40          S7_0.2.0              
##  [97] XVector_0.50.0         ids_1.0.1              htmltools_0.5.8.1     
## [100] carData_3.0-6          dotCall64_1.2          fftwtools_0.9-11      
## [103] scales_1.4.0           png_0.1-8              spatstat.univar_3.1-4 
## [106] knitr_1.51             rstudioapi_0.18.0      tzdb_0.5.0            
## [109] rjson_0.2.23           reshape2_1.4.5         uuid_1.2-2            
## [112] nlme_3.1-168           cachem_1.1.0           zoo_1.8-14            
## [115] Polychrome_1.5.4       stringr_1.6.0          KernSmooth_2.23-26    
## [118] parallel_4.5.2         miniUI_0.1.2           pillar_1.11.1         
## [121] grid_4.5.2             vctrs_0.6.5            colorsGen_1.0.0       
## [124] RANN_2.6.2             gplots_3.3.0           promises_1.5.0        
## [127] BiocSingular_1.26.0    car_3.1-5              beachmat_2.26.0       
## [130] xtable_1.8-8           cluster_2.1.8.1        evaluate_1.0.5        
## [133] cli_3.6.5              locfit_1.5-9.12        compiler_4.5.2        
## [136] rlang_1.1.6            future.apply_1.20.2    ggsignif_0.6.4        
## [139] labeling_0.4.3         plyr_1.8.9             stringi_1.8.7         
## [142] viridisLite_0.4.2      deldir_2.0-4           BiocParallel_1.44.0   
## [145] lazyeval_0.2.2         tiff_0.1-12            spatstat.geom_3.6-0   
## [148] Matrix_1.7-4           RcppHNSW_0.6.0         hms_1.1.4             
## [151] patchwork_1.3.2        future_1.69.0          shiny_1.13.0          
## [154] ROCR_1.0-12            igraph_2.2.1           broom_1.0.12          
## [157] bslib_0.10.0           RCDT_1.3.0
```
