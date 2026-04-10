---
title: "Figure 3: Immune cells and ILCs in lung"
author: "Sandy Kroh"
date: "April 10, 2026"
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
# library(stringr)
# library(glue)
# library(here)
library(readr)
# library(lubridate)
# library(data.table)
# library(clustree)
# library(magrittr)
library(ggpubr)
# library(ggrepel)
# library(readxl)
# library(openxlsx)
```

## Parameters


``` r
set.seed(123)

input_dir <- here::here("1_data_tidying", "Lung_SI_all_cells_all_ALs_files")

output_dir <- here::here("2_visualizations_for_figures", "Fig_3_ILCs_in_lung_files")
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


``` r
SO.lung <- readRDS(paste0(input_dir, "/lung_all_cells_all_ALs.rds"))
dim(SO.lung)
```

```
## [1]    32 67537
```

``` r
SO.lung$AL1 <- gsub("Vessels", "Stromal cells", SO.lung$AL1)
```

Somehow TBET ended up in the metadata. We will put it back as a feature:


``` r
# 1. Extract the current expression matrix from the active assay
# Note: Using slot = "counts". If your MELC data is only stored in "data", change this.
current_matrix <- GetAssayData(SO.lung, assay = "MELC", slot = "counts")

# 2. Extract TBET from metadata and format it as a 1-row matrix
tbet_values <- matrix(SO.lung$TBET, nrow = 1)
rownames(tbet_values) <- "TBET"
colnames(tbet_values) <- colnames(SO.lung) # Ensure cell names match perfectly

# 3. Bind the TBET row to the bottom of the expression matrix
new_matrix <- rbind(current_matrix, tbet_values)

# 4. Overwrite the existing MELC assay with the updated matrix
SO.lung[["MELC"]] <- CreateAssayObject(counts = new_matrix)

# 5. Clean up: Remove TBET from the metadata so it doesn't cause confusion later
SO.lung$TBET <- NULL

# Verify the fix! TBET should now appear at the bottom of this list:
rownames(SO.lung)
```

```
##  [1] "Areg"      "B220"      "CCR6"      "CD117"     "CD11c"     "CD127"     "CD138"     "CD3"       "CD31"      "CD4"       "CD44"      "CD45"      "CD68"      "CD8a"      "CD90"      "EMCN"      "EpCAM"     "ICOS"      "KLRG1"     "Kappa"     "LYVE1"     "MHCII"     "NKp46"     "PDGFRa"    "PDPN"      "Sca1"      "EOMES"     "GATA3"     "GATA3eGFP" "IRF4"      "Ki67"      "RORgt"     "TBET"
```

# Visualization

## Dotplot AL2


``` r
SO.lung$AL2 <- factor(SO.lung$AL2, levels = c(
  "B cells & Plasma cells",
  "Myeloid cells", 
  "T cytotox. cells", 
  "T helper cells", 
  "ILCs"
))


dot_plot <- Seurat::DotPlot(subset(SO.lung, subset = AL1 == "Immune cells"), 
                group.by = "AL2",
                  features = c(
                    "CD90", 
                    "CD127", 
                    "GATA3eGFP", 
                    "KLRG1", 
                    "RORgt", 
                    "CD4", 
                    "CD3", 
                    "CD8a", 
                    "CD68",
                    "CD11c",
                    "MHCII",
                    "B220", 
                    "Kappa", 
                    "CD117"
                  ), 
                cols ="RdBu", assay = "MELC")+   
    RotatedAxis()+
    coord_flip()+
  ggtitle("AL2 - Immune cells\n")+
    theme(axis.text.x=element_text(size=10, angle = 45),
          axis.text.y=element_text(size=10), 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size=11), 
          plot.margin = margin(0.2, 2, 1, 0.2, "cm"), 
          plot.title = element_text(size=12, face = "bold", hjust = 0.5), 
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10, face = "bold"))+ 
  scale_color_gradient2(midpoint = 0, low = "gold", 
                            high = "blue", space = "Lab" )

dot_plot
```

<img src="Fig_3_ILCs_in_lung_files/figure-html/unnamed-chunk-5-1.png" alt="" width="100%" style="display: block; margin: auto;" />

## Ridge plots of ILC subytpes


``` r
ridge_plot <- RidgePlot(subset(SO.lung, subset = AL3 == c("NK cells/ILC1s", "ILC2s", "ILC3s")), 
          features = c("CD127", "CD90","GATA3eGFP", "KLRG1", "EOMES", "TBET",  "RORgt"),
          group.by = "AL3", ncol = 2, same.y.lims = TRUE, assay = "MELC", layer = "counts",
          cols = alpha(c(
            "navy",
            "seagreen2",
            "darkmagenta"
          ), 0.5))+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))

ridge_plot
```

<img src="Fig_3_ILCs_in_lung_files/figure-html/unnamed-chunk-6-1.png" alt="" width="100%" style="display: block; margin: auto;" />

``` r
RidgePlot(subset(SO.lung, subset = AL3 == c("NK cells/ILC1s", "ILC2s", "ILC3s")), 
          features = c("CD127", "CD90","GATA3eGFP",  "KLRG1", "EOMES", "TBET", "NKp46", "RORgt"),
          group.by = "AL3", ncol = 3, same.y.lims = TRUE, assay = "MELC", layer = "counts", y.max = 10, log = TRUE, 
          cols = alpha(c(
            "navy",
            "seagreen2",
            "darkmagenta"
          ), 0.5))+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))
```

<img src="Fig_3_ILCs_in_lung_files/figure-html/unnamed-chunk-6-2.png" alt="" width="100%" style="display: block; margin: auto;" />


``` r
fetched_data <- FetchData(SO.lung, vars = c("CD127", "CD90", "KLRG1", "ICOS",
                                       "TBET", "EOMES", "CCR6", "NKp46", "RORgt",
                                       "GATA3eGFP", "AL3", "Treatment", "Dataset"))

fetched_data <- fetched_data %>%
  filter(AL3 %in% c("NK cells/ILC1s", "ILC2s", "ILC3s", "T helper cells", "T cytotox cells")) %>%
  mutate(Treatment = factor(Treatment, levels = c(
    "CTRL", "1", "2", "3"
  ))) %>%
  group_by(AL3) %>%
  summarise(CD127_median = median(CD127), 
            GATA3eGFP_median = median(GATA3eGFP), 
            CD90_median = median(CD90), 
            ICOS_median = median(ICOS), 
            KLRG1_median = median(KLRG1), 
            RORgt_median = median(RORgt), 
            TBET_median = median(TBET), 
            EOMES_median = median(EOMES), 
            CCR6_median = median(CCR6), 
            NKp46_median = median(NKp46), 
            CD127_div = sd(CD127),
            CD127_skew = moments::skewness(CD127),
            CD127_kurt = moments::kurtosis(CD127),
            GATA3eGFP_div = sd(GATA3eGFP),
            GATA3eGFP_skew = moments::skewness(GATA3eGFP), 
            GATA3eGFP_kurt = moments::kurtosis(GATA3eGFP),
            CD90_div = sd(CD90), 
            CD90_skew = moments::skewness(CD90), 
            CD90_kurt = moments::kurtosis(CD90)
)

fetched_data
```

```
## # A tibble: 5 × 20
##   AL3             CD127_median GATA3eGFP_median CD90_median ICOS_median KLRG1_median RORgt_median TBET_median EOMES_median CCR6_median NKp46_median CD127_div CD127_skew CD127_kurt GATA3eGFP_div GATA3eGFP_skew GATA3eGFP_kurt CD90_div CD90_skew CD90_kurt
##   <fct>                  <dbl>            <dbl>       <dbl>       <dbl>        <dbl>        <dbl>       <dbl>        <dbl>       <dbl>        <dbl>     <dbl>      <dbl>      <dbl>         <dbl>          <dbl>          <dbl>    <dbl>     <dbl>     <dbl>
## 1 NK cells/ILC1s          0                5.59        0           4.31         4.22         0           6.54         6.78        5.45         5.50      2.21      0.562       1.67          2.49         -1.05            2.68     2.06     1.11       2.63
## 2 ILC2s                   5.61             5.96        5.43        5.43         4.81         0           0            0           5.09         0         2.48     -1.06        2.77          3.05         -0.706           1.82     2.81    -0.596      1.74
## 3 ILC3s                   5.65             5.43        5.79        6.15         0            5.27        4.00         0           6.26         4.00      1.52     -1.73        7.04          2.48         -1.07            2.55     2.30    -1.30       3.59
## 4 T helper cells          5.26             5.15        4.43        5.38         0            5.29        3.45         0           5.31         4.66      1.90     -1.43        4.46          2.79         -0.507           1.58     2.54    -0.394      1.68
## 5 T cytotox cells         5.14             0           5.49        4.90         0            0           2.43         3.87        5.40         0         2.03     -1.24        3.59          2.68          0.300           1.37     1.78    -1.40       4.73
```

Healthy


``` r
SO.ilc <- subset(SO.lung, subset = AL3 %in% c("NK cells/ILC1s", "ILC2s", "ILC3s"))
SO.ilc.ctrl <- subset(SO.ilc, subset = Treatment %in% c("CTRL"))

ridge_plot <- RidgePlot(SO.ilc.ctrl, 
          features = c("CD127", "CD90","GATA3eGFP", "KLRG1", "EOMES", "TBET",  "RORgt"),
          group.by = "AL3", ncol = 2, same.y.lims = TRUE, assay = "MELC", layer = "counts",
          cols = alpha(c(
            "navy",
            "seagreen2",
            "darkmagenta"
          ), 0.5))+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))

ridge_plot
```

<img src="Fig_3_ILCs_in_lung_files/figure-html/unnamed-chunk-8-1.png" alt="" width="100%" style="display: block; margin: auto;" />

``` r
RidgePlot(SO.ilc.ctrl, 
          features = c("CD127", "CD90","GATA3eGFP",  "KLRG1", "EOMES", "TBET", "NKp46", "RORgt"),
          group.by = "AL3", ncol = 3, same.y.lims = TRUE, assay = "MELC", layer = "counts", y.max = 10, log = TRUE, 
          cols = alpha(c(
            "navy",
            "seagreen2",
            "darkmagenta"
          ), 0.5))+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))
```

<img src="Fig_3_ILCs_in_lung_files/figure-html/unnamed-chunk-8-2.png" alt="" width="100%" style="display: block; margin: auto;" />

Differential expression


``` r
library(scater)
library(SpatialExperiment)
set.seed(8)

SO.sub <- subset(SO.lung, subset = AL3 %in% c("NK cells/ILC1s", "ILC2s", "ILC3s", 
                                              "T helper cells", "T cytotox cells"))

counts = GetAssayData(SO.sub, layer = "counts")

metadata = SO.sub@meta.data

spatial_coords <- as.matrix(metadata[, c("Location_Center_X", "Location_Center_Y")])

spe <- SpatialExperiment(
  assays = list(counts = counts),
  colData = metadata,
  spatialCoords = spatial_coords
)

# Validate the SpatialExperiment object
print(spe)
```

```
## class: SpatialExperiment 
## dim: 33 5577 
## metadata(0):
## assays(1): counts
## rownames(33): Areg B220 ... RORgt TBET
## rowData names(0):
## colnames(5577): 61961 61962 ... 67536 67537
## colData names(18): orig.ident nCount_MELC ... AL4 sample_id
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## spatialCoords names(2) : Location_Center_X Location_Center_Y
## imgData names(0):
```

``` r
# --- 3. Generate Aggregated Heatmap ---
# Step A: Pseudo-bulk aggregation by cell type
library(scuttle) # required for aggregateAcrossCells

ilc_markers <- c("CD127", "CD90","GATA3eGFP", "KLRG1", "EOMES", "TBET",  
                 "CD3", "CD8a", "CD4", "RORgt", "ICOS", 
                 "MHCII", "CD44", "Ki67", "NKp46")

pbs_annotated <- aggregateAcrossCells(spe, 
                                      ids = spe$AL3, 
                                      subset.row = ilc_markers, 
                                      use.assay.type = "counts", 
                                      statistics = "mean")

# Step B: Extract matrix and Z-score (scale) the data
library(pheatmap)
# Transpose so rows = Cell Types, columns = Markers
mat_annotated <- t(assay(pbs_annotated, "counts"))


# Step C: Plot Heatmap
pheatmap(mat = mat_annotated, 
         scale = "column", 
         clustering_method = "ward.D2",
         color = colorRampPalette(c("#648FFF", "white", "#FFB000"))(101), 
         display_numbers = TRUE, # Optional: shows Z-score values in the cells
         main = "Aggregated Lineage Signatures (Z-scored)")
```

<img src="Fig_3_ILCs_in_lung_files/figure-html/unnamed-chunk-9-1.png" alt="" width="100%" style="display: block; margin: auto;" />


``` r
library(Seurat)
library(ComplexHeatmap)
library(ggplotify) # Required for as.ggplot()

set.seed(8)

# --- 1. Subsetting ---
SO.sub <- subset(SO.lung, subset = AL3 %in% c("NK cells/ILC1s", "ILC2s", "ILC3s"))
SO.sub$AL3 <- droplevels(as.factor(SO.sub$AL3))

# --- 2. Direct Pseudo-bulk Aggregation ---
ilc_markers <- c("CD127", "CD90","GATA3eGFP", "KLRG1", "EOMES", "TBET",  
                 "CD3", "CD4", "RORgt", "ICOS", 
                 "MHCII", "CD44", "Ki67", "NKp46")

# Calculate the mean directly from Seurat
avg_exp <- AverageExpression(SO.sub, features = ilc_markers, slot = "counts", group.by = "AL3")
mat_annotated <- avg_exp[[DefaultAssay(SO.sub)]]

# --- 3. Manually Scale (Z-score) the Matrix ---
mat_scaled <- t(scale(t(mat_annotated)))

# Clean up math errors and cap outliers
mat_scaled[is.na(mat_scaled)] <- 0
mat_scaled[mat_scaled > 2] <- 2
mat_scaled[mat_scaled < -2] <- -2

# --- 4. Plot Heatmap using ComplexHeatmap ---
plot_heat <- ComplexHeatmap::pheatmap(
  mat = mat_scaled, 
  scale = "none", 
  clustering_method = "ward.D2",
  color = colorRampPalette(c("#648FFF", "white", "#FFB000"))(101), 
  breaks = seq(-2, 2, length.out = 102), 
  display_numbers = round(mat_scaled, 2), 
  number_color = "black",
  main = "ILC subtypes",
  treeheight_col = 10,
  treeheight_row = 20,
  name = "Z-Score" 
)

# --- 5. Convert to ggplot object for arrangement ---
# Using grid.grabExpr(draw()) ensures that all the internal sizing, 
# legends, and dendrograms from ComplexHeatmap are captured perfectly.
gg_heat <- as.ggplot(grid::grid.grabExpr(ComplexHeatmap::draw(plot_heat)))

# 'gg_heat' is now a standard ggplot object! 
# You can now combine it using patchwork, e.g.:
# combined_plot <- spatial_scatter_plot + gg_heat
print(gg_heat)
```

<img src="Fig_3_ILCs_in_lung_files/figure-html/unnamed-chunk-10-1.png" alt="" width="100%" style="display: block; margin: auto;" />

## IF overlays



Export TIFF images



### ILC1s/NK cells



### ILC2s


``` r
# load overview overlay image
img <- png::readPNG(
    "D:/Sandy/Promotion/Dissertation/Figures/20210902_1_li_d3/ILC2s/GATA3eGFP-c_CD45-m_CD3-y_ov_ROI.png"
  )

plot_ov <- ggplot() + background_image(img)+
  theme( panel.border = element_rect(colour = "white", fill=NA, size=2)) +
    theme(legend.position = "bottom", 
          plot.margin=margin(0.18,0,0.13,0,"cm"),
          text = element_blank(),
           axis.ticks = element_blank(),
         panel.grid =  element_blank(),
          legend.ticks = element_blank(),
         legend.title=element_blank(),
         legend.key = element_blank(),
        legend.text = element_text(size=14), 
         panel.background = element_rect(fill = 'black', 
                                         color = 'black', size = 1))+    
annotate("text", x=c(500, 1330, 1800, 600), y=c(2000, 2000, 2000, 1800), label= c("GATA3eGFP", "CD45", "CD3", "Identified ILC2s"), ymin = 0, ymax = 2048,
           xmin = 0, xmax = 2048, 
           col=c("cyan", "magenta", "yellow", "white"), size=6, parse=FALSE) +    
annotate("text", x=c(540, 850), y=c(1250, 950), label= c("ROI 1", "ROI 2"), ymin = 0, ymax = 2048,
           xmin = 0, xmax = 2048, 
           col=c("white", "white"), size=4, parse=FALSE)  


# load overview overlay image
img <- png::readPNG(
    "D:/Sandy/Promotion/Dissertation/Figures/20210902_1_li_d3/ILC2s/GATA3eGFP-c_CD45-m_CD3-y_ov_ROI-2.png"
  )

plot_1 <- ggplot() + background_image(img) +
    theme(legend.position = "bottom", 
          plot.margin=margin(0.18,0.1,0.13,0,"cm"),
          text = element_blank(),
           axis.ticks = element_blank(),
         panel.grid =  element_blank(),
          legend.ticks = element_blank(),
         legend.title=element_blank(),
         legend.key = element_blank(),
        legend.text = element_text(size=14), 
         panel.background = element_rect(fill = 'white', 
                                         color = 'white', size = 1)) + 
  annotate("text", x=2, y=9, 
           label= c("ROI 1"), ymin = 0, ymax = 10,
           xmin = 0, xmax = 10, 
           col=c("white"), size=4, parse=FALSE) 


# load overview overlay image
img <- png::readPNG(
    "D:/Sandy/Promotion/Dissertation/Figures/20210902_1_li_d3/ILC2s/GATA3eGFP-c_CD45-m_CD3-y_ov_ROI-1.png"
  )

plot_2 <- ggplot() + background_image(img)  +
    theme(legend.position = "bottom", 
          plot.margin=margin(0.18,0.1,0.13,0,"cm"),
          text = element_blank(),
           axis.ticks = element_blank(),
         panel.grid =  element_blank(),
          legend.ticks = element_blank(),
         legend.title=element_blank(),
         legend.key = element_blank(),
        legend.text = element_text(size=14), 
         panel.background = element_rect(fill = 'white', 
                                         color = 'white', size = 1)) + 
  annotate("text", x=2, y=9, 
           label= c("ROI 2"), ymin = 0, ymax = 10,
           xmin = 0, xmax = 10, 
           col=c("white"), size=4, parse=FALSE) 


zoomed <- ggarrange(plot_1, plot_2, nrow = 2, ncol = 1)

ilc2 <- ggarrange(plot_ov, zoomed, nrow = 1, ncol = 2, 
          heights = 5, 
          widths = c(4, 2))+
    theme(plot.margin=margin(0.5,0,0,0,"cm"))

ilc2
```

<img src="Fig_3_ILCs_in_lung_files/figure-html/unnamed-chunk-14-1.png" alt="" width="100%" style="display: block; margin: auto;" />

### ILC3s





## Combine plots for figure


``` r
figure <- ggarrange(dot_plot, gg_heat, 
          ncol = 2, nrow = 1, labels = c("A", "B"), widths = c(3, 2))+
  theme(plot.margin = margin(0, 0.1, 0, 0, "cm"))

ggarrange(figure, "NONE", ncol = 1, nrow = 2, heights = c(5, 0.3), labels = c("", "C"))
```

<img src="Fig_3_ILCs_in_lung_files/figure-html/unnamed-chunk-17-1.png" alt="" width="100%" style="display: block; margin: auto;" />

## Session Information


``` r
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
## [1] LC_COLLATE=German_Germany.utf8  LC_CTYPE=German_Germany.utf8    LC_MONETARY=German_Germany.utf8 LC_NUMERIC=C                    LC_TIME=German_Germany.utf8    
## 
## time zone: Europe/Berlin
## tzcode source: internal
## 
## attached base packages:
## [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] ggplotify_0.1.3             ComplexHeatmap_2.26.1       pheatmap_1.0.13             SpatialExperiment_1.20.0    scater_1.38.0               scuttle_1.20.0              SingleCellExperiment_1.32.0 SummarizedExperiment_1.40.0 Biobase_2.70.0              GenomicRanges_1.62.0        Seqinfo_1.0.0               IRanges_2.44.0              S4Vectors_0.48.0            BiocGenerics_0.56.0         generics_0.1.4              MatrixGenerics_1.22.0       matrixStats_1.5.0           ggpubr_0.6.2                readr_2.1.6                 ggplot2_4.0.1               dplyr_1.1.4                 Seurat_5.3.1                SeuratObject_5.2.0          sp_2.2-0                   
## 
## loaded via a namespace (and not attached):
##   [1] RcppAnnoy_0.0.22       splines_4.5.2          later_1.4.4            tibble_3.3.0           polyclip_1.10-7        fastDummies_1.7.5      lifecycle_1.0.5        rstatix_0.7.3          doParallel_1.0.17      rprojroot_2.1.1        globals_0.19.0         lattice_0.22-7         MASS_7.3-65            backports_1.5.0        magrittr_2.0.4         plotly_4.12.0          sass_0.4.10            rmarkdown_2.30         jquerylib_0.1.4        yaml_2.3.10            httpuv_1.6.16          otel_0.2.0             sctransform_0.4.2      spam_2.11-1            spatstat.sparse_3.1-0  reticulate_1.44.0      cowplot_1.2.0          pbapply_1.7-4          RColorBrewer_1.1-3     abind_1.4-8            Rtsne_0.17             purrr_1.2.0            yulab.utils_0.2.4      rappdirs_0.3.4         circlize_0.4.17        ggrepel_0.9.6          irlba_2.3.5.1          listenv_0.10.0         spatstat.utils_3.2-0   moments_0.14.1         goftest_1.2-3          RSpectra_0.16-2        spatstat.random_3.4-2  fitdistrplus_1.2-6     parallelly_1.45.1      codetools_0.2-20       DelayedArray_0.36.0    shape_1.4.6.1          tidyselect_1.2.1       farver_2.1.2           viridis_0.6.5         
##  [52] ScaledMatrix_1.18.0    spatstat.explore_3.5-3 jsonlite_2.0.0         GetoptLong_1.1.0       BiocNeighbors_2.4.0    progressr_0.18.0       Formula_1.2-5          iterators_1.0.14       ggridges_0.5.7         survival_3.8-3         foreach_1.5.2          tools_4.5.2            ica_1.0-3              Rcpp_1.1.0             glue_1.8.0             gridExtra_2.3          SparseArray_1.10.1     xfun_0.56              here_1.0.2             withr_3.0.2            fastmap_1.2.0          rsvd_1.0.5             digest_0.6.38          gridGraphics_0.5-1     R6_2.6.1               mime_0.13              colorspace_2.1-2       Cairo_1.7-0            scattermore_1.2        tensor_1.5.1           dichromat_2.0-0.1      spatstat.data_3.1-9    utf8_1.2.6             tidyr_1.3.1            data.table_1.17.8      httr_1.4.8             htmlwidgets_1.6.4      S4Arrays_1.10.0        uwot_0.2.4             pkgconfig_2.0.3        gtable_0.3.6           lmtest_0.9-40          S7_0.2.0               XVector_0.50.0         htmltools_0.5.8.1      carData_3.0-6          dotCall64_1.2          clue_0.3-67            scales_1.4.0           png_0.1-8              spatstat.univar_3.1-4 
## [103] knitr_1.51             rstudioapi_0.18.0      rjson_0.2.23           tzdb_0.5.0             reshape2_1.4.5         nlme_3.1-168           GlobalOptions_0.1.3    cachem_1.1.0           zoo_1.8-14             stringr_1.6.0          KernSmooth_2.23-26     parallel_4.5.2         miniUI_0.1.2           vipor_0.4.7            ggrastr_1.0.2          pillar_1.11.1          vctrs_0.6.5            RANN_2.6.2             promises_1.5.0         BiocSingular_1.26.0    car_3.1-5              beachmat_2.26.0        xtable_1.8-8           cluster_2.1.8.1        beeswarm_0.4.0         evaluate_1.0.5         magick_2.9.0           cli_3.6.5              compiler_4.5.2         crayon_1.5.3           rlang_1.1.6            future.apply_1.20.2    ggsignif_0.6.4         labeling_0.4.3         fs_1.6.6               plyr_1.8.9             ggbeeswarm_0.7.3       stringi_1.8.7          viridisLite_0.4.2      deldir_2.0-4           BiocParallel_1.44.0    lazyeval_0.2.2         spatstat.geom_3.6-0    Matrix_1.7-4           RcppHNSW_0.6.0         hms_1.1.4              patchwork_1.3.2        future_1.69.0          shiny_1.13.0           ROCR_1.0-12            igraph_2.2.1          
## [154] broom_1.0.12           bslib_0.10.0
```
