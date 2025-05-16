---
title: "Lung and SI data"
author: "Sandy Kroh"
date: "May 15, 2025"
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
library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(here)
# library(stringr)
# library(glue)
# library(here)
library(readr)
# library(lubridate)
# library(data.table)
# library(clustree)
# library(magrittr)
# library(ggpubr)
# library(ggrepel)
# library(readxl)
# library(openxlsx)
```

## Parameters


``` r
set.seed(123)

input_dir <- here::here("data")

output_dir <- here::here("1_data_tidying", "Lung_SI_all_cells_all_ALs_files")
dir.create(output_dir)




main_markers <- c(
  "EpCAM", "EMCN", "LYVE1", "PDPN", "PDGFRa", "CD8a", "CD4",
  "CD45", "CD3", "IRF4", "Kappa", "CD11c", "CD127", "GATA3eGFP", "RORgt"
)

immune_markers <- c(
 "CD3", "CD4", "CD8a", "Kappa", "IRF4", "CD11c",
  "CD127", "CD90", "EOMES", "GATA3eGFP", "RORgt", "Ki67",  "KLRG1", "NKp46", "CD117", "Areg", "CCR6", "CD44", "MHCII", "Sca1"
)

ilc_markers <- c(
  "CD3", "CD4", "CD8a",
  "CD127", "CD90", "EOMES", "GATA3eGFP", "RORgt", "KLRG1", "NKp46", "CD117", "CCR6", "MHCII", "Ki67", "Areg", "IRF4", "Sca1", "CD44"
)


cols_nat <- c("magenta", "cyan", "blue", "purple", "green", 
                       "red", "yellow", "olivedrab1", "slateblue1", 
                       "darkcyan", "gold","indianred1", "seagreen", "deeppink", 
                       "orange", "brown", "violet",
                       "deeppink4", "pink", 
                       "grey", "black", "lightgreen", 
                       "#FF0066",  
                       "lightblue", "#FFCC99", "#CC00FF", 
                       "blueviolet",  "goldenrod4", 
                       "navy", "olivedrab", "lightcyan", "seagreen2", "darkviolet", "lightpink", "slateblue4", "olivedrab2")

colfunc <- colorRampPalette(c("darkcyan", "green", "yellow", "magenta", "purple"))
```

# Load data

## SI


``` r
SO.si <- readRDS(paste0(input_dir, "/SO_arcsinh_si_imputed__allALs_allCells.rds"))

SO.si@meta.data$Tissue.area <- gsub("Crypts", "Villi", SO.si@meta.data$Tissue.area)
SO.si <- subset(SO.si, subset = AL2 != "Unresolved")
# SO.villi <- subset(SO.si, subset = Tissue.area == "Villi")
```

Combine Mesenchymal cells I with Fibroblasts and rename mesenchymal cells II


``` r
SO.si$AL2 <- gsub("Epithelia", "Epithelia I" , SO.si$AL2)
SO.si$AL3 <- gsub("Epithelia", "Epithelia I" , SO.si$AL3)

SO.si$AL2 <- gsub("Mesenchymal stromal cells II", "Epithelia II" , SO.si$AL2)
SO.si$AL3 <- gsub("Mesenchymal stromal cells II", "Epithelia II" , SO.si$AL3)

SO.si$AL2 <- gsub("Mesenchymal stromal cells I", "Fibroblasts" , SO.si$AL2)
SO.si$AL3 <- gsub("Mesenchymal stromal cells I", "Fibroblasts" , SO.si$AL3)
```

## Lung


``` r
df <- read_csv(paste0(input_dir, "/20230808_SO_33M_arcsinh_lung_0.04_imputed_quantification_all_cells.csv"), 
    col_types = cols(...1 = col_skip()))

colnames(df) <- gsub("annotation_lvl", "AL", colnames(df))
df$AL3 <- gsub("EMCN CD31 ", "", df$AL3)
df$AL3 <- gsub("LYVE1 CD90 ", "", df$AL3)

CellIDs <- rownames(df)
# 2 TRANSPOSE DF -----
markers <- colnames(df)[12:43]
df_markers <- df %>%
  select(markers)
df_markers <- sapply(df_markers, as.numeric)
df_t <- t(df_markers)
colnames(df_t) <- CellIDs 
df_meta <- df %>%
  select(-markers) 
rownames(df_meta) <- CellIDs 
df_meta$Experiment <- as.factor(df_meta$Experiment)
df_meta$Dataset <- as.factor(df_meta$Dataset)
df_meta$CellID <- as.numeric(df_meta$CellID)
df_meta$Location_Center_X <- as.numeric(df_meta$Location_Center_X)
df_meta$Location_Center_Y <- as.numeric(df_meta$Location_Center_Y)
df_meta$UMAP_1 <- as.numeric(df_meta$UMAP_1)
df_meta$UMAP_2 <- as.numeric(df_meta$UMAP_2)
df_meta$AL1 <- as.factor(df_meta$AL1)
df_meta$AL2 <- as.factor(df_meta$AL2)
df_meta$AL3 <- as.factor(df_meta$AL3)
df_meta$AL4 <- as.factor(df_meta$AL4)

df_meta$max_features <- NULL

SO.lung <- CreateSeuratObject(
  counts = df_t,
  assay = "MELC",
  names.field = 0,
  names.delim = "_",
  meta.data = df_meta,
  project = "SeuratProject"
)


SO.lung@meta.data$AL3 <- factor(SO.lung@meta.data$AL3, levels = c(
  "NK cells/ILC1s", 
  "ILC2s", 
  "ILC3s", 
  "T helper cells", 
  "T cytotox cells", 
  "B cells & Plasma cells", 
  "Myeloid cells", 
  "Blood vessels", 
  "Lymphatics", 
  "LYVE1 CD31 vessels", 
  "Epithelia"
))

new.cluster.ids <- SO.lung@meta.data$AL3
```

# Save data

## Session Information


``` r
saveRDS(SO.si, paste0(output_dir, "/si_all_cells_all_ALs.rds"))
saveRDS(SO.lung, paste0(output_dir, "/lung_all_cells_all_ALs.rds"))
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
## [1] LC_COLLATE=English_Germany.utf8  LC_CTYPE=English_Germany.utf8    LC_MONETARY=English_Germany.utf8 LC_NUMERIC=C                     LC_TIME=English_Germany.utf8    
## 
## time zone: Europe/Berlin
## tzcode source: internal
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] readr_2.1.5        here_1.0.1         ggplot2_3.5.1      dplyr_1.1.4        Seurat_5.2.1       SeuratObject_5.0.2 sp_2.2-0          
## 
## loaded via a namespace (and not attached):
##   [1] deldir_2.0-4           pbapply_1.7-2          gridExtra_2.3          rlang_1.1.5            magrittr_2.0.3         RcppAnnoy_0.0.22       spatstat.geom_3.3-6    matrixStats_1.5.0      ggridges_0.5.6         compiler_4.4.2         png_0.1-8              vctrs_0.6.5            reshape2_1.4.4         stringr_1.5.1          crayon_1.5.3           pkgconfig_2.0.3        fastmap_1.2.0          promises_1.3.2         rmarkdown_2.29         tzdb_0.4.0             bit_4.6.0              purrr_1.0.4            xfun_0.51              cachem_1.1.0           jsonlite_1.9.1         goftest_1.2-3          later_1.4.1            spatstat.utils_3.1-3   irlba_2.3.5.1          parallel_4.4.2         cluster_2.1.6          R6_2.6.1               ica_1.0-3              bslib_0.9.0            stringi_1.8.4          RColorBrewer_1.1-3     spatstat.data_3.1-6    reticulate_1.42.0      parallelly_1.43.0      spatstat.univar_3.1-2  lmtest_0.9-40          jquerylib_0.1.4        scattermore_1.2        Rcpp_1.0.14            knitr_1.50             tensor_1.5             future.apply_1.11.3    zoo_1.8-13             sctransform_0.4.1      httpuv_1.6.15          Matrix_1.7-1          
##  [52] splines_4.4.2          igraph_2.1.4           tidyselect_1.2.1       abind_1.4-8            rstudioapi_0.17.1      yaml_2.3.10            spatstat.random_3.3-3  codetools_0.2-20       miniUI_0.1.2           spatstat.explore_3.4-2 listenv_0.9.1          lattice_0.22-6         tibble_3.2.1           plyr_1.8.9             withr_3.0.2            shiny_1.10.0           ROCR_1.0-11            evaluate_1.0.3         Rtsne_0.17             future_1.40.0          fastDummies_1.7.5      survival_3.7-0         polyclip_1.10-7        fitdistrplus_1.2-2     pillar_1.10.2          KernSmooth_2.23-24     plotly_4.10.4          generics_0.1.3         vroom_1.6.5            rprojroot_2.0.4        RcppHNSW_0.6.0         hms_1.1.3              munsell_0.5.1          scales_1.3.0           globals_0.17.0         xtable_1.8-4           glue_1.8.0             lazyeval_0.2.2         tools_4.4.2            data.table_1.17.0      RSpectra_0.16-2        RANN_2.6.2             dotCall64_1.2          cowplot_1.1.3          grid_4.4.2             tidyr_1.3.1            colorspace_2.1-1       nlme_3.1-166           patchwork_1.3.0        cli_3.6.3              spatstat.sparse_3.1-0 
## [103] spam_2.11-1            viridisLite_0.4.2      uwot_0.2.3             gtable_0.3.6           sass_0.4.10            digest_0.6.37          progressr_0.15.1       ggrepel_0.9.6          htmlwidgets_1.6.4      farver_2.1.2           htmltools_0.5.8.1      lifecycle_1.0.4        httr_1.4.7             mime_0.13              bit64_4.6.0-1          MASS_7.3-61
```
