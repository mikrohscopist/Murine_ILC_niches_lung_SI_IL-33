---
title: "Figure 3: Visual validation AL1"
author: "Sandy Kroh"
date: "Mai 21, 2026"
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
library(rlang)
library(dplyr)
library(png)
library(grid)
library(ggplot2)
library(ggpubr)
```

## Parameters


``` r
set.seed(123)

input_dir <- here::here("data", "images")

output_dir <- here::here("2_visualizations_for_figures", "Fig_3_Visual_validation_AL1_files")
dir.create(output_dir)
```

# IF overlays


``` r
# 1. Define the path
img_path <- paste0(input_dir, "/Visual_validation_AL1_Immune.png")

# 2. Read the image and get dimensions
img <- readPNG(img_path)
h <- nrow(img)
w <- ncol(img)

# 3. Create the ggplot object
# We use coord_fixed to ensure the image doesn't stretch
img_plot_immune <- ggplot() +
  annotation_raster(img, xmin = 0, xmax = w, ymin = 0, ymax = h) +
  coord_fixed() +
  theme_void() +
  # Optional: adds a tiny margin to ensure the edges aren't clipped
  scale_x_continuous(expand = c(0, 0), limits = c(0, w)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, h)) +
  theme(plot.margin = margin(0, 0, 0, 0.7, "cm"))

# Display the plot
# print(img_plot_immune)
```


``` r
# 1. Define the path
img_path <- paste0(input_dir, "/Visual_validation_AL1_Vessels.png")

# 2. Read the image and get dimensions
img <- readPNG(img_path)
h <- nrow(img)
w <- ncol(img)

# 3. Create the ggplot object
# We use coord_fixed to ensure the image doesn't stretch
img_plot_vessels <- ggplot() +
  annotation_raster(img, xmin = 0, xmax = w, ymin = 0, ymax = h) +
  coord_fixed() +
  theme_void() +
  # Optional: adds a tiny margin to ensure the edges aren't clipped
  scale_x_continuous(expand = c(0, 0), limits = c(0, w)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, h)) +
  theme(plot.margin = margin(0, 0, 0, 0.7, "cm"))

# Display the plot
# print(img_plot_vessels)
```


``` r
# 1. Define the path
img_path <- paste0(input_dir, "/Visual_validation_AL1_Epithelia.png")

# 2. Read the image and get dimensions
img <- readPNG(img_path)
h <- nrow(img)
w <- ncol(img)

# 3. Create the ggplot object
# We use coord_fixed to ensure the image doesn't stretch
img_plot_epithelia <- ggplot() +
  annotation_raster(img, xmin = 0, xmax = w, ymin = 0, ymax = h) +
  coord_fixed() +
  theme_void() +
  # Optional: adds a tiny margin to ensure the edges aren't clipped
  scale_x_continuous(expand = c(0, 0), limits = c(0, w)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, h)) +
  theme(plot.margin = margin(0, 0, 0, 0.7, "cm"))

# Display the plot
# print(img_plot_epithelia)
```

# Visualization for figures

## Figure 3


``` r
ggarrange(img_plot_immune, img_plot_vessels, img_plot_epithelia, 
          nrow = 3, ncol = 1, 
          heights = c(4.1, 8.1, 4.1), 
          labels = "AUTO")
```

<embed src="D:/Repositories/2025_Kroh_et_al/Manuscript/EJI_submissions_for_publication/Figures_pdfs/Figure_3_Main-1.pdf" width="100%" style="display: block; margin: auto;" type="application/pdf" />

## Session Information


``` r
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
## [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] ggpubr_0.6.2  ggplot2_4.0.1 png_0.1-8     dplyr_1.1.4   rlang_1.1.6  
## 
## loaded via a namespace (and not attached):
##  [1] sass_0.4.10        generics_0.1.4     tidyr_1.3.1        rstatix_0.7.3      digest_0.6.38      magrittr_2.0.4     evaluate_1.0.5     RColorBrewer_1.1-3 fastmap_1.2.0      rprojroot_2.1.1    jsonlite_2.0.0     backports_1.5.0    Formula_1.2-5      purrr_1.2.0        scales_1.4.0       jquerylib_0.1.4    abind_1.4-8        cli_3.6.5          cowplot_1.2.0      withr_3.0.2        cachem_1.1.0       yaml_2.3.10        otel_0.2.0         tools_4.5.2        ggsignif_0.6.4     here_1.0.2         broom_1.0.12       vctrs_0.6.5        R6_2.6.1           lifecycle_1.0.5    car_3.1-5          pkgconfig_2.0.3    pillar_1.11.1      bslib_0.10.0       gtable_0.3.6       glue_1.8.0         xfun_0.56          tibble_3.3.0       tidyselect_1.2.1   rstudioapi_0.18.0  knitr_1.51         dichromat_2.0-0.1  farver_2.1.2       htmltools_0.5.8.1  rmarkdown_2.30     carData_3.0-6      labeling_0.4.3     compiler_4.5.2     S7_0.2.0
```
