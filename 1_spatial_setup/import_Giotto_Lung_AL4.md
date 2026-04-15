---
title: "Create Giotto object lung data - AL4 GATA3 Th cells"
author: "Sandy Kroh"
date: "April 15, 2026"
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
## Lade nötiges Paket: GiottoClass
```

```
## Giotto Suite 4.2.2
```

``` r
library(ggplot2)
library(ggpubr)
```

## Directories


``` r
output_dir <- here::here("1_spatial_setup", "import_Giotto_Lung_AL4_files")
dir.create(output_dir)
```

```
## Warning in dir.create(output_dir):
## 'D:\Repositories\2025_Kroh_et_al\Murine_ILC_niches_lung_SI_IL-33\1_spatial_setup\import_Giotto_Lung_AL4_files'
## existiert bereits
```

# Configuring Giotto package and raw data

Checking if miniconda and python is in place:


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

Load data frame with annotated cell types:


``` r
datax <- read.csv(paste0(here::here(), "/data/MELC_data_murine_lung_CTRL_D1_D2_D3_withfolders.csv"))
```

## Data prep

We will separate different Th subpopulations by thresholding:


``` r
library(dplyr)
```

```
## 
## Attache Paket: 'dplyr'
```

```
## Die folgenden Objekte sind maskiert von 'package:GiottoClass':
## 
##     intersect, union
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

![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-5-1.png)<!-- -->


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

## Giotto


``` r
samples <- unique(datax$FullInfo)
samples <- samples[!samples %in% "20210906_FOV3_D3"]
gio_list <- NULL
for(samp in samples){

  # cur_datax
  cur_datax <- datax[datax$FullInfo == samp,]

  # sample and image folder
  image_folder <- unique(cur_datax$Sample)
  print(image_folder)

  # data
  cur_data <- cur_datax[,2:34]
  rownames(cur_data) <- cur_datax$CellID
  cur_data <- t(cur_data)

  # metadata
  cur_metadata <- cur_datax[,38:46]
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
## [1] "20210910_1_lu_ctrl"
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
## [1] "20210914_1_lu_ctrl"
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
## [1] "20210922_1_lu_ctrl"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20210910_2_lu_ctrl"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20210914_2_lu_ctrl"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20210922_2_lu_ctrl"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20210910_3_lu_ctrl"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20210914_3_lu_ctrl"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20210922_3_lu_ctrl"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20220311_1"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20220316_1"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20220321_1"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20220311_2"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20220316_2"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20220321_2"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20220311_3"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20220316_3"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20220321_3"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20220325_1"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20220421_1"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20220502_1"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20220325_2"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20220421_2"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20220502_2"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20220325_3"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20220421_3"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20220502_3"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20210902_1_lu_d3"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20210906_1_lu_d3"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20210928_1_lu_d3"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20210902_2_lu_d3"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20210906_2_lu_d3"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20210928_2_lu_d3"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20210902_3_lu_d3"
```

```
## python already initialized in this session
##  active environment : 'giotto_env'
##  python version : 3.10
## Setting spatial network [cell] Delaunay_network
```

```
## [1] "20210928_3_lu_d3"
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

![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-1.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-2.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-3.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-4.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-5.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-6.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-7.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-8.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-9.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-10.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-11.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-12.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-13.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-14.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-15.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-16.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-17.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-18.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-19.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-20.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-21.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-22.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-23.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-24.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-25.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-26.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-27.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-28.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-29.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-30.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-31.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-32.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-33.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-34.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-35.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-36.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-37.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-38.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-39.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-40.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-41.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-42.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-43.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-44.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-45.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-46.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-47.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-48.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-49.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-50.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-51.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-52.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-53.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-54.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-55.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-56.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-57.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-58.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-59.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-60.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-61.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-62.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-63.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-64.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-65.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-66.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-67.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-68.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-69.png)<!-- -->![](import_Giotto_Lung_AL4_files/figure-html/unnamed-chunk-8-70.png)<!-- -->

# Save data

## Session information


``` r
saveRDS(gio_list, file = paste0(output_dir, "/Giotto_data_lung_AL4.rds"))
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
## [1] ggridges_0.5.7     tidyr_1.3.1        dplyr_1.1.4        ggpubr_0.6.2      
## [5] ggplot2_4.0.1      Giotto_4.2.2       GiottoClass_0.4.10
## 
## loaded via a namespace (and not attached):
##   [1] colorRamp2_0.1.0            deldir_2.0-4               
##   [3] gridExtra_2.3               rlang_1.1.6                
##   [5] magrittr_2.0.4              otel_0.2.0                 
##   [7] GiottoUtils_0.2.5           matrixStats_1.5.0          
##   [9] compiler_4.5.2              png_0.1-8                  
##  [11] vctrs_0.6.5                 pkgconfig_2.0.3            
##  [13] SpatialExperiment_1.20.0    fastmap_1.2.0              
##  [15] backports_1.5.0             magick_2.9.0               
##  [17] XVector_0.50.0              labeling_0.4.3             
##  [19] ggraph_2.2.2                rmarkdown_2.30             
##  [21] purrr_1.2.0                 xfun_0.56                  
##  [23] bluster_1.20.0              cachem_1.1.0               
##  [25] jsonlite_2.0.0              DelayedArray_0.36.0        
##  [27] BiocParallel_1.44.0         tweenr_2.0.3               
##  [29] terra_1.8-93                broom_1.0.12               
##  [31] parallel_4.5.2              cluster_2.1.8.1            
##  [33] R6_2.6.1                    bslib_0.10.0               
##  [35] RColorBrewer_1.1-3          reticulate_1.44.0          
##  [37] car_3.1-5                   GenomicRanges_1.62.0       
##  [39] jquerylib_0.1.4             scattermore_1.2            
##  [41] Rcpp_1.1.0                  Seqinfo_1.0.0              
##  [43] SummarizedExperiment_1.40.0 knitr_1.51                 
##  [45] IRanges_2.44.0              Matrix_1.7-4               
##  [47] igraph_2.2.1                tidyselect_1.2.1           
##  [49] rstudioapi_0.18.0           dichromat_2.0-0.1          
##  [51] abind_1.4-8                 yaml_2.3.10                
##  [53] viridis_0.6.5               codetools_0.2-20           
##  [55] lattice_0.22-7              tibble_3.3.0               
##  [57] Biobase_2.70.0              withr_3.0.2                
##  [59] S7_0.2.0                    evaluate_1.0.5             
##  [61] polyclip_1.10-7             pillar_1.11.1              
##  [63] carData_3.0-6               MatrixGenerics_1.22.0      
##  [65] checkmate_2.3.4             stats4_4.5.2               
##  [67] plotly_4.12.0               generics_0.1.4             
##  [69] rprojroot_2.1.1             S4Vectors_0.48.0           
##  [71] scales_1.4.0                gtools_3.9.5               
##  [73] glue_1.8.0                  lazyeval_0.2.2             
##  [75] tools_4.5.2                 GiottoVisuals_0.2.14       
##  [77] BiocNeighbors_2.4.0         data.table_1.17.8          
##  [79] ggsignif_0.6.4              graphlayouts_1.2.2         
##  [81] tidygraph_1.3.1             cowplot_1.2.0              
##  [83] grid_4.5.2                  colorspace_2.1-2           
##  [85] SingleCellExperiment_1.32.0 ggforce_0.5.0              
##  [87] Formula_1.2-5               cli_3.6.5                  
##  [89] rappdirs_0.3.4              S4Arrays_1.10.0            
##  [91] viridisLite_0.4.2           gtable_0.3.6               
##  [93] rstatix_0.7.3               sass_0.4.10                
##  [95] digest_0.6.38               BiocGenerics_0.56.0        
##  [97] SparseArray_1.10.1          ggrepel_0.9.6              
##  [99] rjson_0.2.23                htmlwidgets_1.6.4          
## [101] farver_2.1.2                memoise_2.0.1              
## [103] htmltools_0.5.8.1           lifecycle_1.0.5            
## [105] here_1.0.2                  httr_1.4.8                 
## [107] MASS_7.3-65
```
