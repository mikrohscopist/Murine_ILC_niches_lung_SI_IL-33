---
title: "Lung and SI data"
author: "Sandy Kroh"
date: "July 11, 2025"
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

SO.si$AL1 <- gsub("Epithelia & stroma", "Epithelia" , SO.si$AL1)

df_si <- FetchData(SO.si, c(
  "nCount_MELC",
  "nFeature_MELC",
  "Location_Center_X",
  "Location_Center_Y",
  "Treatment",
  "FOV",
  "Experiment",
  "Organ",
  "MELC.machine",
  "Tissue.area",
  "Dataset",
  "CellID",
  "AL1",
  "AL2",
  "AL3",
  nrow(SO.si)
))

# adapt order
df_si$AL1 <- factor(df_si$AL1, c(
  "Immune cells", "Vessels & stroma", "Epithelia"
))

df_si$AL2 <- factor(df_si$AL2, c(
  "ILCs",
  "CD8+ CD3- IEL",
  "T helper cells",
  "T cytotox. cells",
  "Myeloid cells",
  "B cells",
  "Plasma cells",
  "Fibroblasts",
  "Blood vessels",
  "Lymphatics",
  "Epithelia I",
  "Epithelia II"
))

df_si$AL3 <- factor(df_si$AL3, c(
  "ILC2s",
  "NK cells/ILC1s/ILC3s",
  "CD8+ CD3- IEL",
  "T helper cells",
  "T cytotox. cells",
  "Myeloid cells",
  "B cells",
  "Plasma cells",
  "Fibroblasts",
  "Blood vessels",
  "Lymphatics",
  "Epithelia I",
  "Epithelia II"
))
```

### SI villi

Quantification proportions


``` r
# info of AL is in long format
df_al1 <- df_si %>%
  reshape2::dcast(., Treatment + Dataset + Tissue.area ~ AL1) %>%
  mutate(TotalCellCountFOV = `Immune cells`+ `Vessels & stroma` + `Epithelia`) %>%
  select(Dataset, Treatment, Tissue.area, TotalCellCountFOV, 
         `Immune cells`, `Vessels & stroma`, `Epithelia`) %>%
  as.data.frame()


df_al2 <- df_si %>%
  reshape2::dcast(., Dataset ~ AL2) %>%
  as.data.frame()


df_al3 <- df_si %>%
  reshape2::dcast(., Dataset ~ AL3)%>%
  select(Dataset, `ILC2s`, `NK cells/ILC1s/ILC3s`) %>%
  as.data.frame()

# combine all ALs together
df_al <- merge(df_al1, df_al2, by = "Dataset")
df_al <- merge(df_al, df_al3, by = "Dataset")

# Treatment
df_al$Treatment <- factor(df_al$Treatment, 
                              levels = c("CTRL", "1", "3"))

# Convert to numeric
df_al[, c(4:length(colnames(df_al)))] <- lapply(df_al[, c(4:length(colnames(df_al)))], as.numeric)

# Check right formats
head(df_al)
```

```
##              Dataset Treatment Tissue.area TotalCellCountFOV Immune cells Vessels & stroma Epithelia ILCs CD8+ CD3- IEL T helper cells T cytotox. cells Myeloid cells B cells Plasma cells Fibroblasts Blood vessels Lymphatics Epithelia I Epithelia II ILC2s NK cells/ILC1s/ILC3s
## 1 CTRL_FOV1_20210706      CTRL       Villi              1341          870              106       365  100            11             18              168           365     158           50          38            24         44         203          162    79                   21
## 2 CTRL_FOV1_20210709      CTRL         ILF              3811         3094              230       487  396            22            190              492          1550     293          151         103            52         75         293          194   137                  259
## 3 CTRL_FOV1_20210730      CTRL       Villi              2348         1409              234       705  140            50             13              195           624     271          116          95           116         23         262          443    46                   94
## 4 CTRL_FOV1_20210810      CTRL       Villi              6218         5289              240       689  414            51            436              815          1802    1329          442          84            98         58         448          241   150                  264
## 5 CTRL_FOV2_20210706      CTRL       Villi              1739         1072               82       585  129            14             35              187           340     257          110          52            15         15         410          175    97                   32
## 6 CTRL_FOV2_20210709      CTRL         ILF              4116         3257              148       711  381            54            155              519          1565     363          220          54            29         65         315          396   248                  133
```


``` r
# Proportion of total cell count per FOV
for (element in c(colnames(df_al)[5:length(colnames(df_al))])) {
  number <- paste0("Prop_", element, "_perTotalCountFOV")
  df_al[number] <- round(df_al[element]/df_al["TotalCellCountFOV"]*100, 1)
}

# Proportion of total immune count per FOV
for (element in c(colnames(df_al)[8:14])) {
  number <- paste0("Prop_", element, "_perTotalImmuneFOV")
  df_al[number] <- round(df_al[element]/df_al["Immune cells"]*100, 1)
}
round(rowSums(df_al[39:length(colnames(df_al))], na.rm = TRUE), 0)
```

```
##  [1] 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100
```

``` r
# Proportion of total immune count per FOV
for (element in c(colnames(df_al)[20:21])) {
  number <- paste0("Prop_", element, "_perTotalImmuneFOV")
  df_al[number] <- round(df_al[element]/df_al["Immune cells"]*100, 1)
}

# Proportion of total ILC count per FOV
for (element in c(colnames(df_al)[20:21])) {
  number <- paste0("Prop_", element, "_perTotalILCsFOV")
  df_al[number] <- round(df_al[element]/df_al["ILCs"]*100, 1)
}
round(rowSums(df_al[48:49], na.rm = TRUE), 0)
```

```
##  [1] 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100
```

``` r
colnames(df_al)
```

```
##  [1] "Dataset"                                     "Treatment"                                   "Tissue.area"                                 "TotalCellCountFOV"                           "Immune cells"                                "Vessels & stroma"                            "Epithelia"                                   "ILCs"                                        "CD8+ CD3- IEL"                               "T helper cells"                              "T cytotox. cells"                            "Myeloid cells"                               "B cells"                                     "Plasma cells"                                "Fibroblasts"                                 "Blood vessels"                               "Lymphatics"                                  "Epithelia I"                                 "Epithelia II"                                "ILC2s"                                       "NK cells/ILC1s/ILC3s"                        "Prop_Immune cells_perTotalCountFOV"          "Prop_Vessels & stroma_perTotalCountFOV"      "Prop_Epithelia_perTotalCountFOV"             "Prop_ILCs_perTotalCountFOV"                  "Prop_CD8+ CD3- IEL_perTotalCountFOV"        
## [27] "Prop_T helper cells_perTotalCountFOV"        "Prop_T cytotox. cells_perTotalCountFOV"      "Prop_Myeloid cells_perTotalCountFOV"         "Prop_B cells_perTotalCountFOV"               "Prop_Plasma cells_perTotalCountFOV"          "Prop_Fibroblasts_perTotalCountFOV"           "Prop_Blood vessels_perTotalCountFOV"         "Prop_Lymphatics_perTotalCountFOV"            "Prop_Epithelia I_perTotalCountFOV"           "Prop_Epithelia II_perTotalCountFOV"          "Prop_ILC2s_perTotalCountFOV"                 "Prop_NK cells/ILC1s/ILC3s_perTotalCountFOV"  "Prop_ILCs_perTotalImmuneFOV"                 "Prop_CD8+ CD3- IEL_perTotalImmuneFOV"        "Prop_T helper cells_perTotalImmuneFOV"       "Prop_T cytotox. cells_perTotalImmuneFOV"     "Prop_Myeloid cells_perTotalImmuneFOV"        "Prop_B cells_perTotalImmuneFOV"              "Prop_Plasma cells_perTotalImmuneFOV"         "Prop_ILC2s_perTotalImmuneFOV"                "Prop_NK cells/ILC1s/ILC3s_perTotalImmuneFOV" "Prop_ILC2s_perTotalILCsFOV"                  "Prop_NK cells/ILC1s/ILC3s_perTotalILCsFOV"
```

Separate by tissue area:


``` r
# separate by tissue area
df_villi <- df_al %>%
  filter(Tissue.area == "Villi")
df_ilf <- df_al %>%
  filter(Tissue.area == "ILF")
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

Quantification


``` r
df_lung <- read_csv(here::here("data", "20230808_SO_33M_arcsinh_lung_0.04_imputed_quantification_all_ALs_counts.csv"), 
    col_types = cols(...1 = col_skip()))

# Tidying
df_lung <- df_lung %>% 
  rename(`Endothelial blood cells` = `EMCN CD31 Blood vessels`, 
         `Endothelia & stroma` = `Vessels`, 
         `Lymphatics` = `LYVE1 CD90 Lymphatics`) %>%
  select(-c(`ILC2s A`, `ILC2s B`))

# AL1
for (element in c(colnames(df_lung)[6:8])) {
  number <- paste0("Prop_", element, "_perTotalCountFOV")
  df_lung[number] <- round(df_lung[element]/df_lung["TotalCellCountFOV"]*100, 1)
}
round(rowSums(df_lung[20:22], na.rm = TRUE), 0)
```

```
##  [1] 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100
```

``` r
# AL2
for (element in c(colnames(df_lung)[12:19])) {
  number <- paste0("Prop_", element, "_perTotalImmuneFOV")
  df_lung[number] <- round(df_lung[element]/df_lung["Immune cells"]*100, 1)
}
round(rowSums(df_lung[23:27], na.rm = TRUE), 0)
```

```
##  [1] 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100
```

``` r
# AL3
for (element in c(colnames(df_lung)[17:19])) {
  number <- paste0("Prop_", element, "_perTotalILCsFOV")
  df_lung[number] <- round(df_lung[element]/df_lung["ILCs"]*100, 1)
}
round(rowSums(df_lung[31:33], na.rm = TRUE), 0)
```

```
##  [1] 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100
```

``` r
df_lung$`Tissue area` <- "Lung"
df_lung$Organ <- "Lung"

head(df_lung)
```

```
## # A tibble: 6 × 35
##   Dataset    Experiment   FOV Treatment TotalCellCountFOV `Endothelia & stroma` `Immune cells` Epithelia `Endothelial blood cells` `LYVE1 CD31 vessels` Lymphatics `Myeloid cells` `B cells & Plasma cells`  ILCs `T cytotox cells` `T helper cells` `NK cells/ILC1s` ILC2s ILC3s `Prop_Endothelia & stroma_perTotalCountFOV` `Prop_Immune cells_perTotalCountFOV` Prop_Epithelia_perTotalCountFOV `Prop_Myeloid cells_perTotalImmuneFOV` `Prop_B cells & Plasma cells_perTotalImmuneFOV` Prop_ILCs_perTotalImmuneFOV `Prop_T cytotox cells_perTotalImmuneFOV` `Prop_T helper cells_perTotalImmuneFOV` `Prop_NK cells/ILC1s_perTotalImmuneFOV` Prop_ILC2s_perTotalImmuneFOV Prop_ILC3s_perTotalImmuneFOV `Prop_NK cells/ILC1s_perTotalILCsFOV` Prop_ILC2s_perTotalILCsFOV Prop_ILC3s_perTotalILCsFOV `Tissue area` Organ
##   <chr>           <dbl> <dbl> <chr>                 <dbl>                 <dbl>          <dbl>     <dbl>                     <dbl>                <dbl>      <dbl>           <dbl>                    <dbl> <dbl>             <dbl>            <dbl>            <dbl> <dbl> <dbl>                                       <dbl>                                <dbl>                           <dbl>                                  <dbl>                                           <dbl>                       <dbl>                                    <dbl>                                   <dbl>                                   <dbl>                        <dbl>                        <dbl>                                 <dbl>                      <dbl>                      <dbl> <chr>         <chr>
## 1 20210910_1   20210910     1 CTRL                   1107                   816            244        47                       681                   81         54              46                      103    21                27               47               13     8     0                                        73.7                                 22                               4.2                                   18.9                                            42.2                         8.6                                     11.1                                    19.3                                     5.3                          3.3                          0                                    61.9                       38.1                          0 Lung          Lung 
## 2 20210914_1   20210914     1 CTRL                   1201                   847            243       111                       730                   64         53              50                       94    26                19               54               15    11     0                                        70.5                                 20.2                             9.2                                   20.6                                            38.7                        10.7                                      7.8                                    22.2                                     6.2                          4.5                          0                                    57.7                       42.3                          0 Lung          Lung 
## 3 20210922_1   20210922     1 CTRL                    625                   362            190        73                       290                   30         42              77                       53    25                 3               32                4    18     3                                        57.9                                 30.4                            11.7                                   40.5                                            27.9                        13.2                                      1.6                                    16.8                                     2.1                          9.5                          1.6                                  16                         72                           12 Lung          Lung 
## 4 20210910_2   20210910     2 CTRL                   1149                   732            216       201                       586                  104         42              55                       86    23                10               42                8    15     0                                        63.7                                 18.8                            17.5                                   25.5                                            39.8                        10.6                                      4.6                                    19.4                                     3.7                          6.9                          0                                    34.8                       65.2                          0 Lung          Lung 
## 5 20210914_2   20210914     2 CTRL                   1350                   959            286       105                       848                   64         47              50                      132    23                25               56               12    11     0                                        71                                   21.2                             7.8                                   17.5                                            46.2                         8                                        8.7                                    19.6                                     4.2                          3.8                          0                                    52.2                       47.8                          0 Lung          Lung 
## 6 20210922_2   20210922     2 CTRL                    640                   456            153        31                       358                   66         32              57                       37    21                 5               33                5    16     0                                        71.2                                 23.9                             4.8                                   37.3                                            24.2                        13.7                                      3.3                                    21.6                                     3.3                         10.5                          0                                    23.8                       76.2                          0 Lung          Lung
```

``` r
unique(df_lung$`Tissue area`)
```

```
## [1] "Lung"
```

``` r
df_lung[, c(1:4)] <- lapply(df_lung[, c(1:4)], as.character)

nrow(df_lung)
```

```
## [1] 36
```

``` r
colnames(df_lung)
```

```
##  [1] "Dataset"                                       "Experiment"                                    "FOV"                                           "Treatment"                                     "TotalCellCountFOV"                             "Endothelia & stroma"                           "Immune cells"                                  "Epithelia"                                     "Endothelial blood cells"                       "LYVE1 CD31 vessels"                            "Lymphatics"                                    "Myeloid cells"                                 "B cells & Plasma cells"                        "ILCs"                                          "T cytotox cells"                               "T helper cells"                                "NK cells/ILC1s"                                "ILC2s"                                         "ILC3s"                                         "Prop_Endothelia & stroma_perTotalCountFOV"     "Prop_Immune cells_perTotalCountFOV"            "Prop_Epithelia_perTotalCountFOV"               "Prop_Myeloid cells_perTotalImmuneFOV"          "Prop_B cells & Plasma cells_perTotalImmuneFOV"
## [25] "Prop_ILCs_perTotalImmuneFOV"                   "Prop_T cytotox cells_perTotalImmuneFOV"        "Prop_T helper cells_perTotalImmuneFOV"         "Prop_NK cells/ILC1s_perTotalImmuneFOV"         "Prop_ILC2s_perTotalImmuneFOV"                  "Prop_ILC3s_perTotalImmuneFOV"                  "Prop_NK cells/ILC1s_perTotalILCsFOV"           "Prop_ILC2s_perTotalILCsFOV"                    "Prop_ILC3s_perTotalILCsFOV"                    "Tissue area"                                   "Organ"
```

# Save data

## Session Information


``` r
write.csv(df_lung, paste0(output_dir, "/lung_proportions.csv"))
write.csv(df_si, paste0(output_dir, "/si_all_cells_all_ALs.csv"))
write.csv(df_villi, paste0(output_dir, "/si_villi_proportions.csv"))
write.csv(df_ilf, paste0(output_dir, "/si_ilf_proportions.csv"))
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
## [1] readr_2.1.5        here_1.0.1         ggplot2_3.5.2      dplyr_1.1.4        Seurat_5.2.1       SeuratObject_5.1.0 sp_2.2-0          
## 
## loaded via a namespace (and not attached):
##   [1] deldir_2.0-4           pbapply_1.7-2          gridExtra_2.3          rlang_1.1.5            magrittr_2.0.3         RcppAnnoy_0.0.22       spatstat.geom_3.3-6    matrixStats_1.5.0      ggridges_0.5.6         compiler_4.4.2         png_0.1-8              vctrs_0.6.5            reshape2_1.4.4         stringr_1.5.1          crayon_1.5.3           pkgconfig_2.0.3        fastmap_1.2.0          utf8_1.2.6             promises_1.3.2         rmarkdown_2.29         tzdb_0.4.0             bit_4.6.0              purrr_1.0.4            xfun_0.51              cachem_1.1.0           jsonlite_1.9.1         goftest_1.2-3          later_1.4.1            spatstat.utils_3.1-3   irlba_2.3.5.1          parallel_4.4.2         cluster_2.1.6          R6_2.6.1               ica_1.0-3              bslib_0.9.0            stringi_1.8.4          RColorBrewer_1.1-3     spatstat.data_3.1-6    reticulate_1.42.0      parallelly_1.45.0      spatstat.univar_3.1-2  lmtest_0.9-40          jquerylib_0.1.4        scattermore_1.2        Rcpp_1.0.14            knitr_1.50             tensor_1.5.1           future.apply_1.20.0    zoo_1.8-13             sctransform_0.4.1      httpuv_1.6.15         
##  [52] Matrix_1.7-1           splines_4.4.2          igraph_2.1.4           tidyselect_1.2.1       abind_1.4-8            rstudioapi_0.17.1      yaml_2.3.10            spatstat.random_3.3-3  codetools_0.2-20       miniUI_0.1.2           spatstat.explore_3.4-2 listenv_0.9.1          lattice_0.22-6         tibble_3.2.1           plyr_1.8.9             withr_3.0.2            shiny_1.10.0           ROCR_1.0-11            evaluate_1.0.4         Rtsne_0.17             future_1.58.0          fastDummies_1.7.5      survival_3.7-0         polyclip_1.10-7        fitdistrplus_1.2-2     pillar_1.10.2          KernSmooth_2.23-24     plotly_4.11.0          generics_0.1.4         vroom_1.6.5            rprojroot_2.0.4        RcppHNSW_0.6.0         hms_1.1.3              scales_1.4.0           globals_0.18.0         xtable_1.8-4           glue_1.8.0             lazyeval_0.2.2         tools_4.4.2            data.table_1.17.0      RSpectra_0.16-2        RANN_2.6.2             dotCall64_1.2          cowplot_1.1.3          grid_4.4.2             tidyr_1.3.1            colorspace_2.1-1       nlme_3.1-166           patchwork_1.3.1        cli_3.6.3              spatstat.sparse_3.1-0 
## [103] spam_2.11-1            viridisLite_0.4.2      uwot_0.2.3             gtable_0.3.6           sass_0.4.10            digest_0.6.37          progressr_0.15.1       ggrepel_0.9.6          htmlwidgets_1.6.4      farver_2.1.2           htmltools_0.5.8.1      lifecycle_1.0.4        httr_1.4.7             mime_0.13              bit64_4.6.0-1          MASS_7.3-61
```
