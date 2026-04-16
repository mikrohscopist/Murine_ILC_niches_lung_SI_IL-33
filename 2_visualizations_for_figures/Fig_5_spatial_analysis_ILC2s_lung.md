---
title: "Figure 5: ILC2s localize"
author: "Sandy Kroh"
date: "April 16, 2026"
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
library(SeuratObject)
library(dplyr)
library(rstatix)

# remove.packages("rlang")
# remove.packages("dplyr")
# install.packages("rlang")
# install.packages("dplyr")
library(rlang)
library(dplyr)

if (!requireNamespace("Giotto", quietly = TRUE))
  devtools::install_github("drieslab/Giotto@suite")
if (!requireNamespace("VoltRon", quietly = TRUE))
#  devtools::install_github("Artur-man/VoltRon")
  devtools::install_github("BIMSBbioinfo/VoltRon@dev")

if (!requireNamespace("Seurat", quietly = TRUE))
  install.packages("Seurat")
library(Giotto)
library(Seurat)
library(VoltRon)
library(ggplot2)
library(ggpubr)
library(readr)
library(ggbeeswarm)
library(stringr)
```

## Parameters


``` r
set.seed(123)

input_dir <- here::here("1_data_tidying", "Lung_SI_all_cells_all_ALs_files")

output_dir <- here::here("2_visualizations_for_figures", "Fig_5_spatial_analysis_ILC2s_lung_files")
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

cols_ilcs_lung <- c("darkcyan", "seagreen2", "deeppink4")

cols_treat <- c("CTRL"  = "darkcyan", 
              "1"  = "gold",
              "2"  = "deeppink", 
              "3"  = "slateblue", 
              "D1"  = "gold",
              "D2"  = "deeppink", 
              "D3"  = "slateblue")


# colorblind-friendly palette for conditions
cols_con <- cols_treat
ColorsCellType_conditions <- cols_con

cols_heat <- c(
               "#648FFF", "white", "#FFB000")


order_al3 <- c(
      "NK cells/ILC1s",
    "ILC2s",
    "ILC3s",
    "T cytotox cells",
    "T helper cells",
    "B cells & Plasma cells",
    "Myeloid cells",
    "LYVE1 CD31 vessels",
    "LYVE1 CD90 Lymphatics",
    "EMCN CD31 Blood vessels",
    "Epithelia"

)
```

# Load data


``` r
# from import_Giotto.Rmd
gio_list <- readRDS(here::here("data", "Giotto_data_lung.rds"))

# from import_VoltRon.Rmd
vr_list <- readRDS(here::here("data", "VoltRon_data_lung.rds"))

# original data
metadatax <- read_csv(here::here("data", "MELC_data_murine_lung_CTRL_D1_D2_D3_withfolders.csv"))
metadatax <- metadatax %>%
  mutate(CellType = AL3)  %>%
  filter(Sample != "20210906_3_lu_d3")

unique(metadatax$CellType)
```

```
##  [1] "Epithelia"               "EMCN CD31 Blood vessels" "LYVE1 CD31 vessels"      "LYVE1 CD90 Lymphatics"   "Myeloid cells"           "B cells & Plasma cells"  "NK cells/ILC1s"          "ILC3s"                   "T cytotox cells"         "T helper cells"          "ILC2s"
```

``` r
vr_list_names <- unique(metadatax$Sample)


cell_proximities_list <- list()
for(samp in vr_list_names){
  print(samp)
  cell_proximities_list[[samp]] <-cellProximityEnrichment(
    gobject = gio_list[[samp]],
    cluster_column = 'CellType',
    spatial_network_name = 'Delaunay_network',
    adjust_method = 'fdr',
    number_of_simulations = 1000)
  cell_proximities_list[[samp]] <- cell_proximities_list[[samp]]$enrichm_res
}
```

```
## [1] "20210910_1_lu_ctrl"
## [1] "20210914_1_lu_ctrl"
## [1] "20210922_1_lu_ctrl"
## [1] "20210910_2_lu_ctrl"
## [1] "20210914_2_lu_ctrl"
## [1] "20210922_2_lu_ctrl"
## [1] "20210910_3_lu_ctrl"
## [1] "20210914_3_lu_ctrl"
## [1] "20210922_3_lu_ctrl"
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
## [1] "20210906_1_lu_d3"
## [1] "20210928_1_lu_d3"
## [1] "20210902_2_lu_d3"
## [1] "20210906_2_lu_d3"
## [1] "20210928_2_lu_d3"
## [1] "20210902_3_lu_d3"
## [1] "20210928_3_lu_d3"
```

``` r
vr_merged <- merge(vr_list[[1]], vr_list[-1])
vrImageNames(vr_merged)
```

```
## [1] "image_1"
```

``` r
unique(vr_merged$CellType)
```

```
##  [1] "Epithelia"               "EMCN CD31 Blood vessels" "LYVE1 CD31 vessels"      "LYVE1 CD90 Lymphatics"   "Myeloid cells"           "B cells & Plasma cells"  "NK cells/ILC1s"          "T cytotox cells"         "T helper cells"          "ILC2s"                   "ILC3s"
```

Calculate co-enrichment scores and plot them:


``` r
# FOVs for representative overview images
fovs <- vr_list_names

set_alpha <- 0.35
set_ptsize <- 2
set_nrows <- 1
set_ncols <- 3
cell_shape <- 20


cols_con <- c("darkcyan", "gold", "deeppink", 
                "slateblue")

cols_fov <- c("darkcyan", "gold", "deeppink", 
                "slateblue")

ColorsCellType <-  list(
  #`NK cells/ILC1s/ILC3s` = "cyan", 
  `ILC2s` = "magenta",
  #`ILC3s` = "magenta", 
  `EMCN CD31 Blood vessels` = "green")

uni_celltypes <- unique(vr_merged$CellType)
backgroundlist <- list("EpCAM","CD31","LYVE1","LYVE1","CD11c","B220","EOMES", "CD8a", "CD4", "GATA3eGFP", "RORgt")
names(backgroundlist) <- uni_celltypes
uni_celltypes <- uni_celltypes[!uni_celltypes %in% "ILC2s"]

g_master_list <- list()
for(unic in uni_celltypes){
  
    ### selected cell groups ####
  selected_celltypes <- c("ILC2s", unic)
  # interactions <- c("EMCN CD31 Blood vessels--ILC2s")
  interactions <- unique(cell_proximities_list[[vr_list_names[30]]]$unified_int)
  interactions <- interactions[grepl("ILC2s", interactions) & grepl(unic, interactions)]
  background_image <- backgroundlist[[unic]]
  if(length(interactions) > 0){
     ### get interaction results ####
    interaction_celltypes <- NULL
    for(samp in vr_list_names){
      cur_cell_proximities <- cell_proximities_list[[samp]]
      cur_cell_proximities <- cur_cell_proximities[cur_cell_proximities$unified_int %in% interactions,]
      sample <- unique(metadatax$FullInfo[metadatax$Sample==samp])
      if(nrow(cur_cell_proximities) > 0 & sample != "20210906_FOV3_D3"){
        interaction_celltypes <- rbind(interaction_celltypes,
                                       data.frame(cur_cell_proximities[cur_cell_proximities$unified_int %in% interactions,], 
                                                  experiment = strsplit(sample, split = "_")[[1]][1], fov = strsplit(sample, split = "_")[[1]][2], condition = strsplit(sample, split = "_")[[1]][3]))
      }
    }
    interaction_celltypes$p.adj <- ifelse(interaction_celltypes$enrichm > 0, interaction_celltypes$p.adj_higher, interaction_celltypes$p.adj_lower)

    
    sig_label <- as.character(ifelse(interaction_celltypes$p.adj < 0.1, paste0("*"), ""))

    g_test <- ggplot(interaction_celltypes, aes(x = condition, y = enrichm, fill = condition)) +
      geom_bar(stat = "identity", position = position_dodge2(width=0.9, preserve = "single")) +
      facet_grid(.~condition, scales = "free_x") +
      geom_text(aes(label=sig_label), position=position_dodge2(width=0.9, preserve = "single"), angle = 90, hjust = -0.02, size = 4) +
      ylim(-2,3.5)+
      NoLegend()+
      theme_classic2()+
      scale_fill_manual(values = cols_fov, name = "") +
      theme(axis.text.x = element_text(#angle = 50,
                                       vjust = 1, size = 12, hjust = 0.5, face = "bold"
                                       ),
            axis.text.y = element_text(hjust = 0.5, size = 12),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 12),
            plot.title = element_text(size =14, hjust = 0.5),
            plot.margin = margin(0, 0.5, 0.5, 0.5, "cm"),
            legend.position = "none",
            strip.background=element_blank(),
            strip.background.x= element_blank(),
            strip.text.x = element_text(size = 1, color = "white"),
            panel.grid.major.y = element_line())+
      NoLegend()+  
      ggtitle(gsub("LYVE1 CD90 ", "", interactions)) +
      ylab("Enrichment")
    
    
    g_master_list[[unic]] <- g_test+
        theme(plot.margin = margin(0, 0.5, 0, 0.5, "cm"))
  }
}

# ILC2s around ILC2s -----------------------------------------
interactions <- c("ILC2s--ILC2s")
unic <- "ILC2s"
background_image <- backgroundlist[["ILC2s"]]
interaction_celltypes <- NULL
for(samp in vr_list_names){
    cur_cell_proximities <- cell_proximities_list[[samp]]
    cur_cell_proximities <- cur_cell_proximities[cur_cell_proximities$unified_int %in% interactions,]
    sample <- unique(metadatax$FullInfo[metadatax$Sample==samp])
    if(nrow(cur_cell_proximities) > 0 & sample != "20210906_FOV3_D3"){
      interaction_celltypes <- rbind(interaction_celltypes,
                                     data.frame(cur_cell_proximities[cur_cell_proximities$unified_int %in% interactions,], 
                                                experiment = strsplit(sample, split = "_")[[1]][1], fov = strsplit(sample, split = "_")[[1]][2], condition = strsplit(sample, split = "_")[[1]][3]))
    }
  }
interaction_celltypes$p.adj <- ifelse(interaction_celltypes$enrichm > 0, interaction_celltypes$p.adj_higher, interaction_celltypes$p.adj_lower)
  # plot test results
  # sig_label <- as.character(ifelse(interaction_celltypes$p.adj < 0.1, paste0("p=",round(interaction_celltypes$p.adj,3)), ""))
  sig_label <- as.character(ifelse(interaction_celltypes$p.adj < 0.1, paste0("*"), ""))
    # print(sig_label)
g_test <- ggplot(interaction_celltypes, aes(x = condition, y = enrichm, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge2(width=0.9, preserve = "single")) +
  facet_grid(.~condition, scales = "free_x") +
  geom_text(aes(label=sig_label), position=position_dodge2(width=0.9, preserve = "single"), angle = 90, hjust = -0.02, size = 4) +
  ylim(-2,3.5)+
  NoLegend()+
  theme_classic2()+
  scale_fill_manual(values = cols_fov, name = "") +
  theme(axis.text.x = element_text(#angle = 50,
                                   vjust = 1, size = 12, hjust = 0.5, face = "bold"
                                   ),
        axis.text.y = element_text(hjust = 0.5, size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size =14, hjust = 0.5),
        plot.margin = margin(0, 0.5, 0.5, 0.5, "cm"),
        legend.position = "none",
        strip.background=element_blank(),
        strip.background.x= element_blank(),
        strip.text.x = element_text(size = 1, color = "white"),
        panel.grid.major.y = element_line())+
  NoLegend()+  
  ggtitle(interactions) +
  ylab("Enrichment")


g_master_list[[unic]] <- g_test+
    theme(plot.margin = margin(0, 0.5, 0, 0.5, "cm"))
```

# Visualization

## IF overlays


``` r
img <- png::readPNG(
    "D:/Repositories/2025_Kroh_et_al/Murine_ILC_niches_lung_SI_IL-33/data/images/Fig_5_ILC2_niche/20210906_1_ILC2s-w_CD90-c_CD31-m_LYVE1-y.png"
  )


my_colors <- c("cyan", "magenta", "yellow")

g <- grid::rasterGrob(img, interpolate=TRUE)



df_fov <- metadatax %>%
  mutate(AL1 = recode(
    AL1,
    "Immune cells" = "CD90.2",
    "Endothelia & stroma" = "CD31",
    "Epithelia" = "LYVE1"
  ))

plot_if_1 <- ggplot()+ 
  geom_point(data = df_fov, 
                       aes(x= Loc_X, y= Loc_Y, color = AL1), 
                       size = 1)+
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme(legend.title=element_blank())+ # exclude legend title
  scale_colour_discrete(name  ="MELC IF stainings",
                          breaks=c("CD45", "CD31",
                                   "EpCAM"),
                          labels=c("CD45", "CD31",
                                   "EpCAM")) +
  scale_color_manual(values = my_colors)+ 
    theme(legend.position = "bottom", 
          plot.margin=margin(1,0,0,0,"cm"),
          text = element_blank(),
           axis.ticks = element_blank(),
         panel.grid =  element_blank(),
          legend.ticks = element_blank(),
         legend.title=element_blank(),
         legend.key = element_blank(),
        legend.text = element_text(size=14), 
         panel.background = element_rect(fill = 'black', 
                                         color = 'black', size = 1))+    
  ggplot2::theme(legend.position = "bottom")+
  ggplot2::guides(color=guide_legend(override.aes = list(size=5), ncol=3),
                  fill=guide_legend(ncol = 1,byrow=TRUE))


plot_if_1 <- ggarrange(plot_if_1, 
            nrow = 1, ncol = 1, 
          #widths = c(4.5, 4.5), 
          align = "v",
          font.label=list(size=12),hjust=-0.5
          )+    
  ggplot2::theme(legend.position = "left")
  
plot_if_1
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-5-1.png" alt="" width="100%" style="display: block; margin: auto;" />

Plot ILC2s on overlay of endothelial markers


``` r
set_ptsize <- 2
cell_shape <- 18
set_alpha <- 1
```


``` r
# define markers for the Overlay

# CYAN
marker1 <- "CD90"
# MAGENTA
marker2 <- "FN"
# YELLOW
marker3 <- "LYVE1"
# GREEN
marker4 <- "EpCAM"
# RED
marker5 <- "CD31"

# define cell type of interested that should be plotted on the overlay
celltype_of_interest <- "ILC2s"

ColorsCellTypeSingle <-  list(
  #`NK cells/ILC1s/ILC3s` = "cyan", 
  `ILC2s` = "white",
  #`ILC3s` = "magenta", 
  `EMCN CD31 Blood vessels` = "green")


name_channel_key <- paste0(marker1, "-c_", marker2, "-m_", marker3, "-y_" 
                           # marker4, "-g_", 
                           # marker5, "-r" 
                           )
vr_merged <- combineChannels(vr_merged,
                             channels = c(marker1, marker2, marker3 
                                          # marker4, 
                                          # marker5
                                          ),
                             colors = c("cyan", "magenta", "yellow" 
                                        # "green",
                                        # "red"
                                        ),
                             channel_key = name_channel_key)



cell_type_of_interest <- "ILC2s"

# CTRL ------------------------------------------------------------
plot <- vrSpatialPlot(vr_merged, assay = paste0("Assay", 2), 
                        group.by = "CellType", 
                        group.ids = celltype_of_interest,
                        alpha = set_alpha, 
                        background = c("image_1", name_channel_key), 
                        pt.size = set_ptsize, cell.shape = cell_shape)+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(values = ColorsCellTypeSingle)+
  scale_fill_manual(values = ColorsCellTypeSingle)+
  theme_void()+ NoLegend()+ ggtitle(NULL)+
  ggtitle("CTRL")+
  theme(plot.title = element_text(hjust = 0.5, size = 14))



plot_if_ctrl <- plot +
  annotate("text", x=70, y=970, label= marker1,
           col="cyan", size=3.5, parse=TRUE) +
  annotate("text", x=180, y=970, label= marker2,
           col="magenta", size=3.5, parse=TRUE) +
  annotate("text", x=305, y=970, label= marker3,
           col="yellow", size=3.5, parse=TRUE)+
  # annotate("text", x=305, y=920, label= marker4,
  #          col="green", size=3.5, parse=TRUE)+
  # annotate("text", x=70, y=920, label= marker5,
  #          col="red", size=3.5, parse=TRUE)+
  annotate("segment", x = 680, xend = 985, y = 45, yend = 45, size = 1.6, 
  colour = "white")

plot_if_ctrl
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-7-1.png" alt="" width="100%" style="display: block; margin: auto;" />

``` r
# IL-33 D3 ------------------------------------------------------------
plot <- vrSpatialPlot(vr_merged, assay = paste0("Assay", 35), #2
                        group.by = "CellType", 
                        group.ids = celltype_of_interest,
                        alpha = set_alpha, 
                        background = c("image_1", name_channel_key), 
                        pt.size = set_ptsize, cell.shape = cell_shape)+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(values = ColorsCellTypeSingle)+
  scale_fill_manual(values = ColorsCellTypeSingle)+
  theme_void()+ NoLegend()+ ggtitle(NULL)+
  ggtitle("IL-33 day 3")+
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        text = element_text(size = 12))



plot_if_d3 <- plot +
  annotate("text", x=70, y=970, label= marker1,
           col="cyan", size=5, parse=TRUE) +
  annotate("text", x=180, y=970, label= marker2,
           col="magenta", size=5, parse=TRUE) +
  annotate("text", x=305, y=970, label= marker3,
           col="yellow", size=5, parse=TRUE)+
  # annotate("text", x=305, y=920, label= marker4,
  #          col="green", size=3.5, parse=TRUE)+
  # annotate("text", x=70, y=920, label= marker5,
  #          col="red", size=3.5, parse=TRUE)+
  annotate("segment", x = 680, xend = 985, y = 45, yend = 45, size = 1.6, 
  colour = "white")

plot_if_d3
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-7-2.png" alt="" width="100%" style="display: block; margin: auto;" />

## Coenrichment plots


``` r
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(ggplotify)
library(ggplot2)

# --- 1. Data Gathering (Remains the same) ---
all_interactions <- data.frame()
# We'll also capture the unique cell types here to use for our loop later
all_found_cell_types <- c() 

for (vr_name in vr_list_names) {
  cur_cell_proximities <- cell_proximities_list[[vr_name]]
  if(is.null(cur_cell_proximities) || nrow(cur_cell_proximities) == 0) next
  
  vr_samp <- subset(vr_merged, samples = vr_name)
  vr_meta <- Metadata(vr_samp)
  condition <- paste(unique(vr_meta$Condition), collapse = ", ")
  
  # Capture the cell types available in this dataset
  all_found_cell_types <- unique(c(all_found_cell_types, unique(vr_meta$CellType)))
  
  temp_df <- data.frame(
    Interaction = cur_cell_proximities$unified_int,
    Enrichment = cur_cell_proximities$enrichm,
    FOV = vr_name,
    Condition = condition
  )
  all_interactions <- rbind(all_interactions, temp_df)
}

# --- 2. Initialize Storage ---
all_heatmaps <- list()
target_cell_types <- sort(all_found_cell_types) # Loop through these

# --- 3. START THE LOOP ---
for (ct in target_cell_types) {
  
  # Filter for the current cell type in the loop
  wide_data <- all_interactions %>%
    filter(grepl(ct, Interaction)) %>%
    select(Interaction, FOV, Enrichment) %>%
    pivot_wider(names_from = FOV, values_from = Enrichment, values_fill = 0, values_fn = mean) 
  
  # Skip if this cell type has no interactions recorded
  if(nrow(wide_data) == 0) {
    message(paste("Skipping", ct, "- No interactions found."))
    next
  }

  # Convert to matrix
  enrichmat_global <- as.matrix(wide_data[, -1])
  rownames(enrichmat_global) <- wide_data$Interaction

  # Sync Annotations: Only keep FOVs that are present for THIS cell type
  anno_df_sub <- unique(all_interactions[, c("FOV", "Condition")])
  rownames(anno_df_sub) <- anno_df_sub$FOV
  anno_df_sub <- anno_df_sub[colnames(enrichmat_global), , drop = FALSE]

  # Create Annotation
  ha_top <- HeatmapAnnotation(
    Condition = anno_df_sub$Condition,
    col = list(Condition = cols_treat),
    annotation_name_side = "left"
  )

  # Create Heatmap Object
  curr_heatmap <- Heatmap(
    enrichmat_global,
    col = circlize::colorRamp2(c(-1, 0, 1), cols_heat),
    name = "Enrichment\nScore",
    heatmap_legend_param = list(direction = "vertical"),
    
    column_split = anno_df_sub$Condition,
    cluster_column_slices = FALSE, 
    cluster_columns = FALSE,        
    cluster_rows = TRUE,           
    
    row_title = "Interactions",
    row_title_gp = gpar(fontsize = 11),
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 9),
    show_row_dend = TRUE,
    
    top_annotation = ha_top,
    column_title = paste("Neighborhood Enrichment:", ct),
    column_title_gp = gpar(fontsize = 12, fontface = "bold")
  )

  # Convert to ggplot
  heat_grob <- grid::grid.grabExpr(
    draw(curr_heatmap, 
         heatmap_legend_side = "right",     
         annotation_legend_side = "right",  
         merge_legend = TRUE,
         padding = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
  )

  final_plot <- as.ggplot(heat_grob) + 
    labs(x = "FOV") + 
    theme(
      axis.title.x = element_text(size = 11, margin = margin(t = 10), hjust = 0.37)
    )

  # Store and Print
  all_heatmaps[[ct]] <- final_plot
  print(final_plot)
}
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-8-1.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-8-2.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-8-3.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-8-4.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-8-5.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-8-6.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-8-7.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-8-8.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-8-9.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-8-10.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-8-11.png" alt="" width="100%" style="display: block; margin: auto;" />


``` r
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)

# --- 1. Gather all data into one master dataframe ---
all_interactions <- data.frame()

for (vr_name in vr_list_names) {
  
  cur_cell_proximities <- cell_proximities_list[[vr_name]]
  
  # Skip if this FOV has no interaction data
  if(is.null(cur_cell_proximities) || nrow(cur_cell_proximities) == 0) next
  
  # Get metadata for condition
  vr_samp <- subset(vr_merged, samples = vr_name)
  vr_meta <- Metadata(vr_samp)
  condition <- paste(unique(vr_meta$Condition), collapse = ", ")
  
  # Create a temporary dataframe for this FOV
  temp_df <- data.frame(
    Interaction = cur_cell_proximities$unified_int,
    Enrichment = cur_cell_proximities$enrichm,
    FOV = vr_name,
    Condition = condition
  )
  
  # Append to master dataframe
  all_interactions <- rbind(all_interactions, temp_df)
}

# --- 2. Reshape data into a Wide Matrix ---
wide_data <- all_interactions %>%
  filter(grepl("ILC2s", Interaction)) %>%
  select(Interaction, FOV, Enrichment) %>%
  pivot_wider(names_from = FOV, values_from = Enrichment, values_fill = 0, values_fn = mean) 

# Convert to matrix and assign rownames
enrichmat_global <- as.matrix(wide_data[, -1])
rownames(enrichmat_global) <- wide_data$Interaction

# --- 3. Create Column Annotations (Condition mapping) ---
anno_df <- unique(all_interactions[, c("FOV", "Condition")])
rownames(anno_df) <- anno_df$FOV

# Reorder the rows
anno_df <- anno_df[colnames(enrichmat_global), , drop = FALSE]

condition_colors <- cols_treat

ha_top <- HeatmapAnnotation(
  Condition = anno_df$Condition,
  col = list(Condition = condition_colors),
  annotation_name_side = "left"
)

library(ggplotify)
library(ggplot2)

# --- 4. Draw the Unified Master Heatmap ---
master_heatmap <- Heatmap(
  enrichmat_global,
  
  col = circlize::colorRamp2(c(-1, 0, 1), cols_heat),
  name = "Enrichment\nScore",
  heatmap_legend_param = list(direction = "vertical"),
  
  column_split = anno_df$Condition,
  cluster_column_slices = FALSE, 
  cluster_columns = FALSE,        
  cluster_rows = TRUE,           
  
  # Y-Axis Label
  row_title = "Interactions",
  row_title_gp = gpar(fontsize = 11#, fontface = "bold"
                      ),

  # Aesthetics
  show_column_names = FALSE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 9),
  show_row_dend = TRUE,
  
  # Annotations and Titles
  top_annotation = ha_top,
  column_title = "Neighborhood Enrichment ILC2s",
  column_title_gp = gpar(fontsize = 12, fontface = "bold")
)

# --- 5. Convert to ggplot2 and explicitly position the X-axis ---

# First, capture the raw heatmap (with normal, minimal padding)
heat_grob <- grid::grid.grabExpr(
  draw(master_heatmap, 
       heatmap_legend_side = "right",     
       annotation_legend_side = "right",  
       merge_legend = TRUE,
       padding = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
)

# Convert to a standard ggplot object
gg_heat <- as.ggplot(heat_grob)

# Add the X-axis label and manually nudge it to center under the heatmap body
final_plot_ILC2s <- gg_heat + 
  labs(x = "FOV") + 
  theme(
    axis.title.x = element_text(
      size = 11, 
      # face = "bold", 
      margin = margin(t = 10), # Adds a nice gap between the plot and the text
      hjust = 0.37 
    )
  )

# Print your beautiful, perfectly aligned plot!
print(final_plot_ILC2s)
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-9-1.png" alt="" width="100%" style="display: block; margin: auto;" />


``` r
# fine tune the co-enrichment plot
# ILC2s around lymphatics -----------------------------------------
interactions <- c("ILC2s--LYVE1 CD90 Lymphatics")
unic <- "LYVE1 CD90 Lymphatics"
background_image <- backgroundlist[["LYVE1 CD90 Lymphatics"]]
interaction_celltypes <- NULL
for(samp in vr_list_names){
    cur_cell_proximities <- cell_proximities_list[[samp]]
    cur_cell_proximities <- cur_cell_proximities[cur_cell_proximities$unified_int %in% interactions,]
    sample <- unique(metadatax$FullInfo[metadatax$Sample==samp])
    if(nrow(cur_cell_proximities) > 0 & sample != "20210906_FOV3_D3"){
      interaction_celltypes <- rbind(interaction_celltypes,
                                     data.frame(cur_cell_proximities[cur_cell_proximities$unified_int %in% interactions,], 
                                                experiment = strsplit(sample, split = "_")[[1]][1], fov = strsplit(sample, split = "_")[[1]][2], condition = strsplit(sample, split = "_")[[1]][3]))
    }
  }
interaction_celltypes$p.adj <- ifelse(interaction_celltypes$enrichm > 0, interaction_celltypes$p.adj_higher, interaction_celltypes$p.adj_lower)
  # plot test results
  # sig_label <- as.character(ifelse(interaction_celltypes$p.adj < 0.1, paste0("p=",round(interaction_celltypes$p.adj,3)), ""))
  sig_label <- as.character(ifelse(interaction_celltypes$p.adj < 0.1, paste0("*"), ""))
    # print(sig_label)
plot_coenrichment_lymp <- ggplot(interaction_celltypes, aes(x = condition, y = enrichm, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge2(width=0.9, preserve = "single")) +
  facet_grid(.~condition, scales = "free_x") +
  geom_text(aes(label=sig_label), position=position_dodge2(width=0.9, preserve = "single"), angle = 90, hjust = -0.02, size = 4) +
  ylim(-1.6,3.6)+
  NoLegend()+
  theme_classic2()+
  scale_fill_manual(values = cols_fov, name = "") +
  ggtitle(gsub("LYVE1 CD90 ", "", interactions)) +
  theme(axis.text.x = element_text(#angle = 50,
                                   vjust = 1, size = 9, hjust = 0.5
                                   ),
        axis.text.y = element_text(hjust = 0.5, size = 9),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size =12, hjust = 0.5, face = "bold"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.position = "none",
        strip.background=element_blank(),
        strip.background.x= element_blank(),
        strip.text.x = element_text(size = 1, color = "white"),
        panel.grid.major.y = element_line())+
  NoLegend()+  
  ylab("Enrichment")


# fine tune the co-enrichment plot
# ILC2s around ILC2s -----------------------------------------
interactions <- c("ILC2s--ILC2s")
unic <- "ILC2s"
background_image <- backgroundlist[["ILC2s"]]
interaction_celltypes <- NULL
for(samp in vr_list_names){
    cur_cell_proximities <- cell_proximities_list[[samp]]
    cur_cell_proximities <- cur_cell_proximities[cur_cell_proximities$unified_int %in% interactions,]
    sample <- unique(metadatax$FullInfo[metadatax$Sample==samp])
    if(nrow(cur_cell_proximities) > 0 & sample != "20210906_FOV3_D3"){
      interaction_celltypes <- rbind(interaction_celltypes,
                                     data.frame(cur_cell_proximities[cur_cell_proximities$unified_int %in% interactions,], 
                                                experiment = strsplit(sample, split = "_")[[1]][1], fov = strsplit(sample, split = "_")[[1]][2], condition = strsplit(sample, split = "_")[[1]][3]))
    }
  }
interaction_celltypes$p.adj <- ifelse(interaction_celltypes$enrichm > 0, interaction_celltypes$p.adj_higher, interaction_celltypes$p.adj_lower)
  # plot test results
  # sig_label <- as.character(ifelse(interaction_celltypes$p.adj < 0.1, paste0("p=",round(interaction_celltypes$p.adj,3)), ""))
  sig_label <- as.character(ifelse(interaction_celltypes$p.adj < 0.1, paste0("*"), ""))
    # print(sig_label)
plot_coenrichment_ilc2 <- ggplot(interaction_celltypes, aes(x = condition, y = enrichm, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge2(width=0.9, preserve = "single")) +
  facet_grid(.~condition, scales = "free_x") +
  geom_text(aes(label=sig_label), position=position_dodge2(width=0.9, preserve = "single"), angle = 90, hjust = -0.02, size = 4) +
  ylim(-1.6,3.6)+
  NoLegend()+
  theme_classic2()+
  scale_fill_manual(values = cols_fov, name = "") +
  ggtitle(gsub("LYVE1 CD90 ", "", interactions)) +
  theme(axis.text.x = element_text(#angle = 50,
                                   vjust = 1, size = 9, hjust = 0.5
                                   ),
        axis.text.y = element_text(hjust = 0.5, size = 9),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size =12, hjust = 0.5, face = "bold"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.position = "none",
        strip.background=element_blank(),
        strip.background.x= element_blank(),
        strip.text.x = element_text(size = 1, color = "white"),
        panel.grid.major.y = element_line())+
  NoLegend()+  
  ylab("Enrichment")


# fine tune the co-enrichment plot
# ILC2s around Myeloid cells -----------------------------------------
interactions <- c("ILC2s--Myeloid cells")
unic <- "Myeloid cells"
background_image <- backgroundlist[["Myeloid cells"]]
interaction_celltypes <- NULL
for(samp in vr_list_names){
    cur_cell_proximities <- cell_proximities_list[[samp]]
    cur_cell_proximities <- cur_cell_proximities[cur_cell_proximities$unified_int %in% interactions,]
    sample <- unique(metadatax$FullInfo[metadatax$Sample==samp])
    if(nrow(cur_cell_proximities) > 0 & sample != "20210906_FOV3_D3"){
      interaction_celltypes <- rbind(interaction_celltypes,
                                     data.frame(cur_cell_proximities[cur_cell_proximities$unified_int %in% interactions,], 
                                                experiment = strsplit(sample, split = "_")[[1]][1], fov = strsplit(sample, split = "_")[[1]][2], condition = strsplit(sample, split = "_")[[1]][3]))
    }
  }
interaction_celltypes$p.adj <- ifelse(interaction_celltypes$enrichm > 0, interaction_celltypes$p.adj_higher, interaction_celltypes$p.adj_lower)
  # plot test results
  # sig_label <- as.character(ifelse(interaction_celltypes$p.adj < 0.1, paste0("p=",round(interaction_celltypes$p.adj,3)), ""))
  sig_label <- as.character(ifelse(interaction_celltypes$p.adj < 0.1, paste0("*"), ""))
    # print(sig_label)
plot_coenrichment_my <- ggplot(interaction_celltypes, aes(x = condition, y = enrichm, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge2(width=0.9, preserve = "single")) +
  facet_grid(.~condition, scales = "free_x") +
  geom_text(aes(label=sig_label), position=position_dodge2(width=0.9, preserve = "single"), angle = 90, hjust = -0.02, size = 4) +
  ylim(-1.6,3.6)+
  NoLegend()+
  theme_classic2()+
  scale_fill_manual(values = cols_fov, name = "") +
  ggtitle(gsub("LYVE1 CD90 ", "", interactions)) +
  theme(axis.text.x = element_text(#angle = 50,
                                   vjust = 1, size = 9, hjust = 0.5
                                   ),
        axis.text.y = element_text(hjust = 0.5, size = 9),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size =12, hjust = 0.5, face = "bold"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.position = "none",
        strip.background=element_blank(),
        strip.background.x= element_blank(),
        strip.text.x = element_text(size = 1, color = "white"),
        panel.grid.major.y = element_line())+
  NoLegend()+  
  ylab("Enrichment")

coenrichment_fig <- ggarrange(plot_coenrichment_ilc2, plot_coenrichment_my, plot_coenrichment_lymp, ncol = 3, nrow = 1)

coenrichment_fig_ann <- annotate_figure(coenrichment_fig, 
                left = text_grob("Coenrichment score", color = "black", rot = 90, 
                                 size = 11))

coenrichment_fig_ann
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-10-1.png" alt="" width="100%" style="display: block; margin: auto;" />

## CIN analysis


``` r
library(stringr)
tissuearea <- "Lung"
input_dir <- here::here("data", "CIN_analysis", tissuearea, "/")


# List all radii calculated 
list_files_rad <- list.files(path=input_dir, 
                         pattern=NULL, all.files=FALSE,
                         full.names=FALSE)

list_files_rad
```

```
## [1] "10_micm" "15_micm" "20_micm" "25_micm"
```

``` r
df_cin_lung = data.frame()

# define a dataframe to collect all data from one tissue area
df_radius <- data.frame()
for (defined_radius in list_files_rad) {
  print(defined_radius)
  list_files <- list.files(path=paste0(input_dir,
                                   defined_radius, 
                                   "/"), 
                       pattern=NULL, all.files=FALSE,
                       full.names=FALSE)
  # define a dataframe to collect the data for all datasets of one 
  df_cin = data.frame()
  for (element in list_files) {
    dir <- paste0(input_dir,
                  defined_radius, 
                  "/", element)
    df_sub <- read_csv(dir, 
      col_types = cols(...1 = col_skip()))
    df_sub$Filename <- paste0(defined_radius, 
                  "/", element)
    df_cin <- rbind(df_cin, df_sub)
  }

  df_cin$Radius <- defined_radius
  df_radius <- rbind(df_radius, df_cin)
}
```

```
## [1] "10_micm"
```

```
## [1] "15_micm"
```

```
## [1] "20_micm"
```

```
## [1] "25_micm"
```

``` r
# collect data from current tissue area to datafram df_all_cin
df_cin_lung <- df_radius
df_cin_lung$Radius <- gsub("_micm", " \u03BCm", df_cin_lung$Radius)
df_cin_lung$Experiment <- str_sub(df_cin_lung$Dataset,-8,-1)
df_cin_lung$FOV <- str_sub(df_cin_lung$Dataset,-10,-10)  
df_cin_lung$Treatment <- str_extract(df_cin_lung$Dataset, "[^_]+")
df_cin_lung$ExpID <- paste0(df_cin_lung$Experiment, "_", df_cin_lung$FOV)
colnames(df_cin_lung) <- gsub("Reference cell", "Reference", 
                                colnames(df_cin_lung))

colnames(df_cin_lung) <- gsub("NK cells & ILC1s", "NK cells/ILC1s", 
                              colnames(df_cin_lung))
df_cin_lung$Reference <- gsub("NK cells & ILC1s", "NK cells/ILC1s",
                              df_cin_lung$Reference)
head(df_cin_lung)
```

```
## # A tibble: 6 × 21
##   Dataset          Epithelia `EMCN CD31 Blood vessels` `LYVE1 CD31 vessels` `LYVE1 CD90 Lymphatics` `Myeloid cells` `B cells & Plasma cells` `NK cells/ILC1s` ILC3s `T cytotox cells` `T helper cells` ILC2s Reference               Organ Tissue.area Filename                                          Radius Experiment FOV   Treatment ExpID     
##   <chr>                <dbl>                     <dbl>                <dbl>                   <dbl>           <dbl>                    <dbl>            <dbl> <dbl>             <dbl>            <dbl> <dbl> <chr>                   <chr> <chr>       <chr>                                             <chr>  <chr>      <chr> <chr>     <chr>     
## 1 D3_FOV1_20210902        88                         4                    3                       1               3                        1                0     0                 0                0     0 Epithelia               Lung  Lung        10_micm/SPIAT_20210902_1_rad10_micm_Lung_freq.csv 10 μm  20210902   1     D3        20210902_1
## 2 D3_FOV1_20210902         2                        70                    6                       1               5                       11                1     0                 1                3     1 EMCN CD31 Blood vessels Lung  Lung        10_micm/SPIAT_20210902_1_rad10_micm_Lung_freq.csv 10 μm  20210902   1     D3        20210902_1
## 3 D3_FOV1_20210902         5                        18                   59                       2               7                        6                1     0                 1                2     1 LYVE1 CD31 vessels      Lung  Lung        10_micm/SPIAT_20210902_1_rad10_micm_Lung_freq.csv 10 μm  20210902   1     D3        20210902_1
## 4 D3_FOV1_20210902         5                        14                    8                      47              15                        8                0     0                 0                0     2 LYVE1 CD90 Lymphatics   Lung  Lung        10_micm/SPIAT_20210902_1_rad10_micm_Lung_freq.csv 10 μm  20210902   1     D3        20210902_1
## 5 D3_FOV1_20210902         5                        14                    6                       3              59                        4                1     0                 0                2     5 Myeloid cells           Lung  Lung        10_micm/SPIAT_20210902_1_rad10_micm_Lung_freq.csv 10 μm  20210902   1     D3        20210902_1
## 6 D3_FOV1_20210902         1                        27                    5                       1               3                       54                1     0                 2                4     1 B cells & Plasma cells  Lung  Lung        10_micm/SPIAT_20210902_1_rad10_micm_Lung_freq.csv 10 μm  20210902   1     D3        20210902_1
```

``` r
colnames(df_cin_lung)
```

```
##  [1] "Dataset"                 "Epithelia"               "EMCN CD31 Blood vessels" "LYVE1 CD31 vessels"      "LYVE1 CD90 Lymphatics"   "Myeloid cells"           "B cells & Plasma cells"  "NK cells/ILC1s"          "ILC3s"                   "T cytotox cells"         "T helper cells"          "ILC2s"                   "Reference"               "Organ"                   "Tissue.area"             "Filename"                "Radius"                  "Experiment"              "FOV"                     "Treatment"               "ExpID"
```

``` r
df_cin_lung$Treatment <- gsub("D", "", df_cin_lung$Treatment)
```


``` r
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggbeeswarm)
library(ggpubr)
library(patchwork)
library(rstatix) 

# 1. Target cells, References, and Radii
target_cells <- c("Epithelia", "EMCN CD31 Blood vessels", "LYVE1 CD31 vessels", 
                  "LYVE1 CD90 Lymphatics", "Myeloid cells", "B cells & Plasma cells", 
                  "NK cells/ILC1s", "ILC3s", "T cytotox cells", "T helper cells", "ILC2s")

all_refs <- unique(df_cin_lung$Reference)
all_radii <- unique(df_cin_lung$Radius)

# Initialize the global master list 
master_plot_list <- list()

# 2. Loop through Radii and References
for (set_radius in all_radii) {
  for (ref_type in all_refs) {
    
    base_df <- df_cin_lung %>%
      filter(Reference == ref_type) %>%
      filter(Radius == set_radius) %>%
      mutate(Treatment = factor(Treatment, levels = c("CTRL", "1", "2", "3")))
    
    if (nrow(base_df) == 0) next
    
    # CRITICAL: This temporary list MUST be outside the target_cells loop
    freq_plot_list <- list()
    
    # 3. Loop through the 11 target cell types
    for (tgt_cell in target_cells) {
      
      stats_df <- base_df %>%
        select(Treatment, FOV, Frequency = all_of(tgt_cell))
      
      global_max <- max(stats_df$Frequency, na.rm = TRUE)
      
      n_df <- stats_df %>%
        group_by(Treatment) %>%
        summarise(n = n(), .groups = "drop") %>%
        mutate(
          Label = paste0("n=", n),
          Y_pos = ifelse(is.finite(global_max), global_max * 1.05, 1) 
        )
      
      # --- STATISTICAL TESTING ---
      stat.test <- tryCatch({
        stats_df %>%
          wilcox_test(Frequency ~ Treatment, ref.group = "CTRL") %>%
          add_significance() %>% 
          add_xy_position(x = "Treatment")
      }, error = function(e) NULL)
      
      plot_title <- tgt_cell
      
      # --- BUILD THE PLOT ---
      p <- ggplot(stats_df, aes(x = Treatment, y = Frequency)) +
        geom_boxplot(outlier.colour = NA, fill = "white", alpha = 0.5) +
        geom_beeswarm(aes(color = Treatment), size = 1.5, cex = 2) + 
        # geom_text(data = n_df, aes(x = Treatment, y = Y_pos, label = Label), 
        #           size = 3.2, fontface = "italic", inherit.aes = FALSE) +
        ggtitle(gsub("LYVE1 CD90 ", "", plot_title)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.45))) +
        theme_classic() + 
        theme(axis.text.x = element_text(size = 9#, face = "bold"
                                         ),
              axis.text.y = element_text(size = 9),
              axis.title.x = element_blank(), 
              axis.title.y = element_text(size = 9), 
              plot.title = element_text(size = 11, hjust = 0.5),
              legend.position = "none") +
        ylab(paste0("Frequency within r=", set_radius))+
        scale_color_manual(values = cols_treat) 
      
      # --- THE BULLETPROOF STATS FIX ---
      if (!is.null(stat.test) && nrow(stat.test) > 0) {
        # Dynamically detect if rstatix used multiple testing (p.adj.signif) 
        # or a single comparison (p.signif) for this specific cell type!
        sig_label <- ifelse("p.adj.signif" %in% colnames(stat.test), "p.adj.signif", "p.signif")
        
        p <- p + stat_pvalue_manual(stat.test, 
                                    label = sig_label, # Uses the dynamically detected column
                                    hide.ns = TRUE, 
                                    y.position = global_max * 1.1, 
                                    step.increase = 0.1, 
                                    size = 4)
      }
      
      # Save to the temporary grid list
      freq_plot_list[[plot_title]] <- p
      
      # Save to the master list
      plot_key <- paste(set_radius, ref_type, tgt_cell, sep = "_")
      master_plot_list[[plot_key]] <- p
    }
    
    # --- PRINT THE PATCHWORK GRID ---
    if (length(freq_plot_list) > 0) {
      final_cin_grid <- wrap_plots(freq_plot_list, ncol = 4) +
        plot_annotation(
          title = paste0("Spatial Neighborhood (", set_radius, ") - Reference: ", ref_type),
          subtitle = 'Points represent individual FOVs; Stats: Wilcoxon vs CTRL',
          theme = theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 12, hjust = 0.5)
          )
        )
      
      message("Printing grid for: ", ref_type, " at ", set_radius)
      print(final_cin_grid)
    }
  }
}
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-1.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-2.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-3.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-4.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-5.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-6.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-7.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-8.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-9.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-10.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-11.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-12.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-13.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-14.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-15.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-16.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-17.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-18.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-19.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-20.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-21.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-22.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-23.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-24.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-25.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-26.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-27.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-28.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-29.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-30.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-31.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-32.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-33.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-34.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-35.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-36.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-37.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-38.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-39.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-40.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-41.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-42.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-43.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-44.png" alt="" width="100%" style="display: block; margin: auto;" />


``` r
cin_ilc <- 
  ggarrange(
  master_plot_list[["10 μm_ILC2s_ILC2s"]], 
  master_plot_list[["10 μm_ILC2s_Myeloid cells"]], 
  master_plot_list[["10 μm_ILC2s_LYVE1 CD90 Lymphatics"]], 
  ncol = 3, nrow = 1)+
  plot_annotation(
          title = paste0("ILC2s (CIN = 10 µm)"),
          theme = theme(
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
          )
        )


cin_my <- ggarrange(
  master_plot_list[["10 μm_Myeloid cells_ILC2s"]], 
  master_plot_list[["10 μm_Myeloid cells_Myeloid cells"]], 
  master_plot_list[["10 μm_Myeloid cells_LYVE1 CD90 Lymphatics"]], 
  ncol = 3, nrow = 1)+
  plot_annotation(
          title = paste0("Myeloid cells (CIN = 10 µm)"),
          theme = theme(
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
          )
        )

cin_ly <- ggarrange(
  master_plot_list[["10 μm_LYVE1 CD90 Lymphatics_ILC2s"]], 
  master_plot_list[["10 μm_LYVE1 CD90 Lymphatics_Myeloid cells"]], 
  master_plot_list[["10 μm_LYVE1 CD90 Lymphatics_LYVE1 CD90 Lymphatics"]], 
  ncol = 3, nrow = 1)+
  plot_annotation(
          title = paste0("Lymphatics (CIN = 10 µm)"),
          theme = theme(
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
          )
        )

plot_cin <- ggarrange(cin_ilc, cin_my, cin_ly, ncol = 1, nrow = 3, labels = c("C", "D", "E"))

plot_cin
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-13-1.png" alt="" width="100%" style="display: block; margin: auto;" />

## Minimum distance

### Load data


``` r
set.seed(8)
radius <- 10

cols_nat <- c("magenta", "cyan", "blue", "purple", "green", 
                       "red", "yellow", "olivedrab1", "slateblue1", 
                       "darkcyan", "seagreen", "deeppink", 
                       "orange", "brown", "violet",
                       "deeppink4", "pink", 
                       "grey", "black", "lightgreen", 
                       "#FF0066", "gold", 
                       "lightblue", "#FFCC99", "#CC00FF", 
                       "blueviolet",  "goldenrod4", 
                       "indianred1", "navy", "olivedrab", "lightcyan", "seagreen2", "darkviolet", "lightpink", "slateblue4", "olivedrab2")
                


ColorsCellType <-  list(`NK cells/ILC1s` = "darkcyan", 
                        `ILC2s` = "seagreen2", 
                        `ILC3s` = "darkmagenta", 
                        `T helper cells` = "deeppink",
                        `T cytotox cells` = "slateblue", 
                        `Myeloid cells` = "burlywood", 
                        `B cells & Plasma cells` = "indianred1",
                        `LYVE1 CD31 vessels` = "darkgreen", 
                        `LYVE1 CD90 Lymphatics` = "yellow", 
                        `EMCN CD31 Blood vessels` = "red", 
                        `Epithelia` = "green")
ColorsCellType 
```

```
## $`NK cells/ILC1s`
## [1] "darkcyan"
## 
## $ILC2s
## [1] "seagreen2"
## 
## $ILC3s
## [1] "darkmagenta"
## 
## $`T helper cells`
## [1] "deeppink"
## 
## $`T cytotox cells`
## [1] "slateblue"
## 
## $`Myeloid cells`
## [1] "burlywood"
## 
## $`B cells & Plasma cells`
## [1] "indianred1"
## 
## $`LYVE1 CD31 vessels`
## [1] "darkgreen"
## 
## $`LYVE1 CD90 Lymphatics`
## [1] "yellow"
## 
## $`EMCN CD31 Blood vessels`
## [1] "red"
## 
## $Epithelia
## [1] "green"
```

``` r
cols_treat <- c("darkcyan", "gold", "deeppink", "slateblue")

celltypes <- c(
  "NK cells/ILC1s", 
  "ILC2s", 
  "ILC3s", 
  "T helper cells", 
  "T cytotox cells", 
  "Myeloid cells", 
  "B cells & Plasma cells", 
  "LYVE1 CD31 vessels", 
  "LYVE1 CD90 Lymphatics", 
  "EMCN CD31 Blood vessels", 
  "Epithelia"
)

df_dist_all <- read_csv(here::here("data", "SPIAT_min_dist_all_lung.csv"), 
    col_types = cols(...1 = col_skip()))

df_dist_all$RefType <- factor(df_dist_all$RefType, levels = celltypes)
df_dist_all$NearestType <- factor(df_dist_all$NearestType, levels = celltypes)

df_dist_all <- df_dist_all %>%
  mutate(Treatment = paste0("D", Treatment)) %>%
  mutate(Treatment = gsub("DCTRL", "CTRL", Treatment))

head(df_dist_all)
```

```
## # A tibble: 6 × 11
##   ExpID      RefCell     RefType   NearestCell NearestType             Distance Pair                                 Name                                     Experiment   FOV Treatment
##   <chr>      <chr>       <fct>     <chr>       <fct>                      <dbl> <chr>                                <chr>                                         <dbl> <dbl> <chr>    
## 1 20210902_1 Cell_175084 Epithelia Cell_175108 EMCN CD31 Blood vessels     4.52 Epithelia -- EMCN CD31 Blood vessels SPIAT_20210902_1_min_distances_pairs.csv   20210902     1 D3       
## 2 20210902_1 Cell_175115 Epithelia Cell_175126 EMCN CD31 Blood vessels     4.53 Epithelia -- EMCN CD31 Blood vessels SPIAT_20210902_1_min_distances_pairs.csv   20210902     1 D3       
## 3 20210902_1 Cell_175154 Epithelia Cell_175123 EMCN CD31 Blood vessels     4.52 Epithelia -- EMCN CD31 Blood vessels SPIAT_20210902_1_min_distances_pairs.csv   20210902     1 D3       
## 4 20210902_1 Cell_175252 Epithelia Cell_175201 EMCN CD31 Blood vessels     6.98 Epithelia -- EMCN CD31 Blood vessels SPIAT_20210902_1_min_distances_pairs.csv   20210902     1 D3       
## 5 20210902_1 Cell_175423 Epithelia Cell_175453 EMCN CD31 Blood vessels     5.29 Epithelia -- EMCN CD31 Blood vessels SPIAT_20210902_1_min_distances_pairs.csv   20210902     1 D3       
## 6 20210902_1 Cell_175570 Epithelia Cell_175592 EMCN CD31 Blood vessels     7.93 Epithelia -- EMCN CD31 Blood vessels SPIAT_20210902_1_min_distances_pairs.csv   20210902     1 D3
```


``` r
# filter condition and reference cell types to compare ILC subtypes and 
# T cells
df_dist_ref <- df_dist_all %>%
  # filter(Treatment == condition) %>%
  filter(RefType == "ILC2s"
         |
           RefType == "NK cells/ILC1s"| RefType == "ILC3s"|
           RefType == "T helper cells"| RefType == "T cytotox cells"| RefType == "Myeloid cells"| RefType == "B cells & Plasma cells"
         )

unique(df_dist_ref$Treatment)
```

```
## [1] "D3"   "CTRL" "D1"   "D2"
```

``` r
df_dist_ref$Treatment <- factor(df_dist_ref$Treatment, levels = c(
  "CTRL", "D1", "D2", "D3"
))
unique(df_dist_ref$RefType)
```

```
## [1] Myeloid cells          B cells & Plasma cells NK cells/ILC1s         ILC3s                  T cytotox cells        T helper cells         ILC2s                 
## Levels: NK cells/ILC1s ILC2s ILC3s T helper cells T cytotox cells Myeloid cells B cells & Plasma cells LYVE1 CD31 vessels LYVE1 CD90 Lymphatics EMCN CD31 Blood vessels Epithelia
```

``` r
# LYMPHATICS ------------------------------------------------------------
celltype_of_interest <- "LYVE1 CD90 Lymphatics"
ypos <- 100


my_plot_list <- list()

celltype_of_interest <- "LYVE1 CD90 Lymphatics"

for (condition in c("CTRL", "D1", "D2", "D3")) {
  
  plot_data <- df_dist_ref %>%
    filter(Treatment == condition) %>%
    filter(NearestType == celltype_of_interest)
  
  # Test for statistical significance of ILC2s to the other cell types
  stat.test <- plot_data %>%
    dunn_test(Distance ~ RefType) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance()
  stat.test <- stat.test %>% add_xy_position(x = "RefType")
  
  stat.test <- stat.test %>%
    filter(group1 == "ILC2s"|group2 == "ILC2s"
           )
  
  # Create lav´bels that depict mean value 
  Labs = plot_data %>% 
    group_by(RefType) %>%
    summarise(lab_text = paste0(round(median(Distance), 0), " µm"), lab_pos = quantile(Distance)[3]
                )
  
  # create plot
  plot <- ggplot(plot_data, aes(x=RefType , 
            y = Distance, 
            fill = "RefType")) +
    geom_boxplot(fill="white", outliers = FALSE)+
    geom_beeswarm(aes(color = RefType), 
                  alpha = 0.2, size = 0.1, cex = 0.2)+
    stat_pvalue_manual(stat.test,
                       size = 6,
                       hide.ns = TRUE, 
                       y.position = 250, 
                       step.increase = 0.1
                       )+
    scale_y_continuous(expand = c(0, 0), limits = c(0,400))+
    rotate_x_text(angle = 30)+
  ggtitle(condition)+
  scale_color_manual(values = ColorsCellType)+
  theme_classic2()+
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, 
                                   # vjust = 1, hjust = 0.5, 
                                   size = 9),
    axis.title.x = element_blank(),
    axis.text.y = element_text(hjust = 0.5, size = 9),
    axis.title.y = element_text(size = 9),
          plot.title = element_text(size = 10, hjust = 0.5, 
                                    face = "bold"),
    plot.margin = margin(0.5, 1, 0.5, 1, "cm"),
    legend.title = element_text(size =11),
    legend.text = element_text(size =12))+
    NoLegend()+
    xlab("Cell types")+
    ylab("Minimal distance [µm]")+
    geom_text(aes(y = 150 , #lab_pos, 
                  label = lab_text, vjust = -0.5), angle = 30,
              data = Labs, size=3)
  plot
  
  assigned_name <- gsub(" ", "", paste0(celltype_of_interest, "_", condition))
  assign(assigned_name, plot )
  my_plot_list[[assigned_name]] <- plot
  

}

ggarrange(plotlist = my_plot_list[c(1:4)], nrow = 1, ncol = 4)
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-15-1.png" alt="" width="100%" style="display: block; margin: auto;" />

``` r
dist_lymph <- ggarrange(plotlist = my_plot_list[c(1, 4)], nrow = 1, ncol = 2, labels = c("C", "D"))

dist_lymph
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-15-2.png" alt="" width="100%" style="display: block; margin: auto;" />


``` r
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggbeeswarm)
library(ggpubr)
library(patchwork)
library(rstatix)

# --- 1. SETTINGS ---
celltype_of_interest <- "LYVE1 CD90 Lymphatics"
ref_cell_types <- unique(df_dist_ref$RefType)
treatment_levels <- c("CTRL", "D1", "D2", "D3")

master_dist_list <- list() 
grid_plot_list <- list()   

# --- 2. LOOP THROUGH REFERENCE CELL TYPES ---
for (ref_type in ref_cell_types) {
  
  # Step A: Filter data
  stats_df <- df_dist_ref %>%
    filter(RefType == ref_type) %>%
    filter(NearestType == celltype_of_interest) %>%
    mutate(Treatment = factor(Treatment, levels = treatment_levels))
  
  # Skip if there's no data for this specific pairing
  if (nrow(stats_df) < 5) next
  
  # Step B: Statistical Testing (vs CTRL only)
  stat.test <- tryCatch({
    stats_df %>%
      wilcox_test(Distance ~ Treatment, ref.group = "CTRL") %>%
      add_significance() %>%
      add_xy_position(x = "Treatment")
  }, error = function(e) NULL)
  
  # Step C: Calculate Median Distance for Labels
  median_df <- stats_df %>%
    group_by(Treatment) %>%
    summarise(med_dist = median(Distance, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      Label = paste0(round(med_dist, 0), " \u00b5m"), 
      Y_pos = 250 # Median labels kept at 400
    )
  
  # Step D: Build the Plot
  plot_title <- paste0(ref_type) 
  
  p <- ggplot(stats_df, aes(x = Treatment, y = Distance)) +
    geom_boxplot(outlier.colour = NA, fill = "white", alpha = 0.5) +
    geom_beeswarm(aes(color = Treatment), size = 0.1, cex = 0.2, alpha = 0.2) + 
    
    # Median Labels at y=400
    geom_text(data = median_df, aes(x = Treatment, y = Y_pos,
                                    label = Label), 
              size = 3, angle = 30, color = "black", 
              inherit.aes = FALSE) +
    ggtitle(plot_title) +
    # Y-axis fixed from 0 to 1000
    scale_y_continuous(limits = c(0, 500), expand = c(0, 0)) + 
    theme_classic() + 
    theme(axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 9),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 9),
          plot.title = element_text(size = 10, hjust = 0.5, 
                                    face = "bold"),
          plot.margin = margin(0.5, 1, 0.5, 1, "cm"),
          legend.position = "none") +
    ylab("Distance [\u00b5m]") +
    scale_color_manual(values = cols_treat)

  # Step E: Add Stat Brackets
  if (!is.null(stat.test) && nrow(stat.test) > 0) {
    sig_col <- if ("p.adj.signif" %in% colnames(stat.test)) "p.adj.signif" else "p.signif"
    
    p <- p + stat_pvalue_manual(
      stat.test, 
      label = sig_col, 
      hide.ns = TRUE, 
      tip.length = 0.01,
      # --- FIXED: Significance brackets start at y = 300 ---
      y.position = 400, 
      step.increase = 0.1, 
      size = 3.5
    )
  }
  
  # Save to lists
  grid_plot_list[[ref_type]] <- p
  plot_key <- paste0(ref_type, "_to_Lymphatics")
  master_dist_list[[plot_key]] <- p
}

# --- 3. PRINT THE FINAL GRID ---
if (length(grid_plot_list) > 0) {
  final_grid <- wrap_plots(grid_plot_list, ncol = 4) +
    plot_annotation(
      title = paste0("Distances to: ", celltype_of_interest),
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )
  
  print(final_grid)
}
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-16-1.png" alt="" width="100%" style="display: block; margin: auto;" />


``` r
# 1. Arrange the plots
p_arranged <- ggarrange(
  master_dist_list[["ILC2s_to_Lymphatics"]], 
  master_dist_list[["Myeloid cells_to_Lymphatics"]], 
  ncol = 2, nrow = 1  # Set to 2 rows for a vertical stack
)

# 2. Add the title using annotate_figure
final_plot <- annotate_figure(p_arranged,
                top = text_grob(label = "Minimal distance to Lymphatics", 
                                 size = 12, face = "bold"))

print(final_plot)
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-17-1.png" alt="" width="100%" style="display: block; margin: auto;" />


``` r
plot_dist <- ggarrange(dist_lymph, final_plot, nrow = 2, ncol = 1, heights = c(4, 3))

plot_dist
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-18-1.png" alt="" width="100%" style="display: block; margin: auto;" />

## Cellular Microenvironments (Tissue Hubs/Domains)

### Data preparation


``` r
df_all <- metadatax

df_all <- df_all %>%
  select(-c(X, ...1)) %>%
  mutate(Dataset = FullInfo, 
         Filename = Dataset, 
         Condition = str_split_i(Dataset, "_", -1))

# max values x and y before conversion
max(df_all$Loc_X)
```

```
## [1] 2033.556
```

``` r
max(df_all$Loc_Y)
```

```
## [1] 2034.286
```

``` r
df_all <- df_all %>%
  mutate(Loc_X = Loc_X * 0.325,
         Loc_Y = Loc_Y * 0.325)

# max values x and y after conversion
max(df_all$Loc_X)
```

```
## [1] 660.9056
```

``` r
max(df_all$Loc_Y)
```

```
## [1] 661.1429
```

``` r
head(df_all)
```

```
## # A tibble: 6 × 47
##    Areg  B220  CCR6 CD117 CD11c CD127 CD138   CD3  CD31   CD4  CD44  CD45  CD68  CD8a  CD90  EMCN EpCAM  ICOS KLRG1 Kappa LYVE1 MHCII NKp46 PDGFRa  PDPN  Sca1 EOMES GATA3 GATA3eGFP  IRF4  Ki67 RORgt  TBET Loc_X Loc_Y CellID FullInfo           Experiment FOV   Condition CellType  AL1       AL2       AL3       Sample             Dataset            Filename          
##   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>     <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl> <chr>                   <dbl> <chr> <chr>     <chr>     <chr>     <chr>     <chr>     <chr>              <chr>              <chr>             
## 1  5.49  0     0     0        0  0     7.81  0     0        0  0        0     0  0        0     0  0     0     0     0     0     7.31     0   0     8.34  7.24     0     0         0  5.31  0     0     0    239.   27.5  12029 20210910_FOV1_CTRL   20210910 FOV1  CTRL      Epithelia Epithelia Epithelia Epithelia 20210910_1_lu_ctrl 20210910_FOV1_CTRL 20210910_FOV1_CTRL
## 2  4.73  0     5.16  0        0  3.02  6.36  0     0        0  0        0     0  7.92     0     0  6.56  6.86  4.97  4.86  4.25  4.75     0   5.22  7.91  6.65     0     0         0  6.38  0     0     0     51.2  42.7  12050 20210910_FOV1_CTRL   20210910 FOV1  CTRL      Epithelia Epithelia Epithelia Epithelia 20210910_1_lu_ctrl 20210910_FOV1_CTRL 20210910_FOV1_CTRL
## 3  5.82  2.54  6.58  3.34     0  0     6.29  3.26  6.16     0  5.52     0     0  3.41     0     0  5.92  7.27  0     4.69  5.68  5.07     0   5.87  7.75  7.31     0     0         0  3.91  3.49  0     6.89 363.   46.6  12057 20210910_FOV1_CTRL   20210910 FOV1  CTRL      Epithelia Epithelia Epithelia Epithelia 20210910_1_lu_ctrl 20210910_FOV1_CTRL 20210910_FOV1_CTRL
## 4  6.40  0     6.77  0        0  5.39  0     0     0        0  0        0     0  0        0     0  0     7.80  0     0     0     0        0   0     8.65  0        0     0         0  0     5.21  0     0     36.8  74.4  12127 20210910_FOV1_CTRL   20210910 FOV1  CTRL      Epithelia Epithelia Epithelia Epithelia 20210910_1_lu_ctrl 20210910_FOV1_CTRL 20210910_FOV1_CTRL
## 5  6.52  0     7.35  5.86     0  3.60  5.30  4.29  6.75     0  5.06     0     0  0        0     0  6.29  7.26  0     1.82  4.86  5.61     0   4.84  7.38  6.96     0     0         0  6.26  4.70  5.78  0    447.  120.   12203 20210910_FOV1_CTRL   20210910 FOV1  CTRL      Epithelia Epithelia Epithelia Epithelia 20210910_1_lu_ctrl 20210910_FOV1_CTRL 20210910_FOV1_CTRL
## 6  4.99  0     6.13  6.10     0  0     4.89  0     7.59     0  0        0     0  0        0     0  6.06  6.44  3.61  4.47  0     0        0   0     8.13  7.52     0     0         0  5.91  4.13  0     0    477.  130.   12216 20210910_FOV1_CTRL   20210910 FOV1  CTRL      Epithelia Epithelia Epithelia Epithelia 20210910_1_lu_ctrl 20210910_FOV1_CTRL 20210910_FOV1_CTRL
```

### Niche analysis 20 µm

This analysis clusters the tissue based on local composition to find **functional zones**. For example, it might identify an "Inflamed Immune Hub" (dense T cells + B cells + ILCs) or a "Vascular Niche" (Vessels + Fibroblasts).

-   **The Biological Question:** Do periodontitis samples form tertiary lymphoid structures or specific immune clusters that healthy tissues lack?


``` r
set.seed(8)

library(dbscan)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(rstatix)
library(patchwork)
library(ggpubr)

# ========================================================================
# PART 1: COMPUTE GLOBAL CELLULAR NEIGHBORHOODS
# ========================================================================
# Parameters
search_radius <- 20  # Look at cells within 20 µm
num_niches <- 4        # Number of distinct spatial environments to identify (K)
all_cell_types <- unique(df_all$CellType)

# Create a globally unique ID for every cell so we can map them back flawlessly
df_all$Global_ID <- paste0(df_all$Filename, "_", df_all$CellID)
cn_matrix_list <- list()

print("Extracting local neighborhoods for all cells...")
```

```
## [1] "Extracting local neighborhoods for all cells..."
```

``` r
for (fov_name in unique(df_all$Filename)) {
  
  df_fov <- df_all %>% filter(Filename == fov_name)
  coords <- as.matrix(df_fov[, c("Loc_X", "Loc_Y")])
  
  # Find all neighbors within the radius (Extremely fast using kd-trees)
  neighbors <- dbscan::frNN(coords, eps = search_radius)$id
  
  # Empty matrix to store the counts (Rows = cells, Columns = cell types)
  fov_counts <- matrix(0, nrow = nrow(df_fov), ncol = length(all_cell_types))
  colnames(fov_counts) <- all_cell_types
  rownames(fov_counts) <- df_fov$Global_ID
  
  # Tabulate neighbors for each cell
  for (i in seq_along(neighbors)) {
    # Include the central cell itself in its own neighborhood
    neighbor_idx <- c(i, neighbors[[i]])
    neighbor_types <- df_fov$CellType[neighbor_idx]
    
    # Factor guarantees missing cell types get recorded as 0
    counts <- table(factor(neighbor_types, levels = all_cell_types))
    fov_counts[i, ] <- as.numeric(counts)
  }
  
  # Convert counts to fractions (proportions from 0 to 1)
  fov_freqs <- fov_counts / rowSums(fov_counts)
  cn_matrix_list[[fov_name]] <- fov_freqs
}

# Combine all images into one master matrix
master_cn_matrix <- do.call(rbind, cn_matrix_list)

print("Clustering microenvironments using K-means...")
```

```
## [1] "Clustering microenvironments using K-means..."
```

``` r
set.seed(8) # Lock seed for reproducibility
kmeans_res <- kmeans(master_cn_matrix, centers = num_niches, nstart = 25)

# Add the raw cluster assignments temporarily
df_all$Raw_Niche <- paste0("Niche ", kmeans_res$cluster[df_all$Global_ID])

# Map raw clusters to your specific biological names
# NOTE: The order of K-means clusters (1-5) can shift if the data changes slightly.
# Assuming Niche 1 = Epithelial, Niche 2 = B cell, etc. based on your prompt.
niche_names <- c(
  "Niche 1" = "Mixed Myeloid/LEC niche",
  "Niche 2" = "Mixed BPC/BEC",
  "Niche 3" = "Epithelial niche",
  "Niche 4" = "Blood endothelial niche"
)

df_all$Tissue_Niche <- unname(niche_names[df_all$Raw_Niche])
df_all$Tissue_Niche <- factor(df_all$Tissue_Niche, levels = c(
  "Mixed BPC/BEC",        
  "Mixed Myeloid/LEC niche",
  "Epithelial niche", 
  "Blood endothelial niche"       
))

niche_colors <- c(
  "Epithelial niche" = "#117733", 
  "Mixed Myeloid/LEC niche" = "#DDCC77", 
  "Blood endothelial niche" = "#882255", 
  "Mixed BPC/BEC" = "#332288"
)


# ========================================================================
# PART 2: WHAT IS IN EACH NICHE? (Composition Profiling)
# ========================================================================
niche_composition <- df_all %>%
  mutate(CellType = factor(CellType, levels = c(
    "NK cells/ILC1s",
    "ILC2s",
    "ILC3s",
    "T cytotox cells",
    "T helper cells",
    "B cells & Plasma cells",
    "Myeloid cells",
    "LYVE1 CD31 vessels",
    "LYVE1 CD90 Lymphatics",
    "EMCN CD31 Blood vessels",
    "Epithelia"
  ))) %>%
  group_by(Tissue_Niche, CellType) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Tissue_Niche) %>%
  mutate(Fraction = Count / sum(Count))

table_plot_al3 <- ggtexttable(niche_composition, 
                          rows = NULL
                          )

plot_comp_al3 <- ggplot(niche_composition, aes(x = Tissue_Niche, y = Fraction, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", color = "black", 
           linewidth = 0.05) +
  scale_fill_manual(values = ColorsCellType) +
  theme_classic() +
  labs(title = "Cellular Microenvironments",
       x = "Identified Tissue Niches", y = "Average Proportion of Cell Types", fill = "Cell Type") +
  theme(
    axis.text.x = element_text(angle = 45, face = "bold", size = 9 , 
                                   hjust = 1),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(face = "bold", hjust = -1, size = 12), 
    plot.margin = margin(0.2, 0, 0.5, 0, "cm"))

print(plot_comp_al3)
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-20-1.png" alt="" width="100%" style="display: block; margin: auto;" />

``` r
niche_composition <- df_all %>%
  group_by(Tissue_Niche, AL1) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Tissue_Niche) %>%
  mutate(Fraction = Count / sum(Count))


table_plot_al1 <- ggtexttable(niche_composition, 
                          rows = NULL
                          )

plot_comp_al1 <- ggplot(niche_composition, aes(x = Tissue_Niche, y = Fraction, fill = AL1)) +
  geom_bar(stat = "identity", position = "stack", color = "black", 
           linewidth = 0.05) +
  scale_fill_manual(values = c("yellow", "aquamarine3", "deeppink")) +
  theme_classic() +
  labs(title = "AL1 - Cellular Microenvironments",
       x = "Identified Tissue Niches", 
       y = "Average Proportion of Cell Types", 
       fill = "AL1 cell types") +
  theme(
    axis.text.x = element_text(angle = 45, size = 9, face = "bold" , 
                                   hjust = 1),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_blank())

print(plot_comp_al1)
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-20-2.png" alt="" width="100%" style="display: block; margin: auto;" />


``` r
# ========================================================================
# PART 2: WHAT IS IN EACH NICHE? (Composition Profiling)
# ========================================================================
niche_composition <- df_all %>%
  filter(AL2 == "ILCs") %>%
  mutate(CellType = factor(CellType, levels = c(
    "NK cells/ILC1s",
    "ILC2s",
    "ILC3s",
    "T cytotox cells",
    "T helper cells",
    "B cells & Plasma cells",
    "Myeloid cells",
    "LYVE1 CD31 vessels",
    "LYVE1 CD90 Lymphatics",
    "EMCN CD31 Blood vessels",
    "Epithelia"
  ))) %>%
  group_by(Tissue_Niche, CellType) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Tissue_Niche) %>%
  mutate(Fraction = Count / sum(Count))

table_plot_ilc <- ggtexttable(niche_composition, 
                          rows = NULL
                          )

plot_comp_ilc <- ggplot(niche_composition, 
                        aes(x = Tissue_Niche, y = Fraction, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", color = "black", 
           linewidth = 0.05) +
  scale_fill_manual(values = ColorsCellType) +
  theme_classic() +
  labs(title = "AL3 ILC-subtypes - Cellular Microenvironments",
       x = "Identified Tissue Niches", 
       y = "Average Proportion of Cell Types", 
       # CHANGED legend title here
       fill = "AL3 ILC subtypes") + 
  theme(
    axis.text.x = element_text(angle = 45, face = "bold", size = 9, 
                                   hjust = 1),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_blank())

print(plot_comp_ilc)
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-21-1.png" alt="" width="100%" style="display: block; margin: auto;" />


``` r
upper_supp_plot <- ggarrange(plot_comp_al1, "NONE", plot_comp_ilc, ncol = 3, nrow = 1, 
          widths = c(10, 1, 10), labels = c("A", "B", ""))

upper_supp_plot
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-22-1.png" alt="" width="100%" style="display: block; margin: auto;" />


``` r
print(table_plot_al1)
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-23-1.png" alt="" width="100%" style="display: block; margin: auto;" />

``` r
print(table_plot_al3)
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-23-2.png" alt="" width="100%" style="display: block; margin: auto;" />

``` r
print(table_plot_ilc)
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-23-3.png" alt="" width="100%" style="display: block; margin: auto;" />


``` r
# --- 1. Calculate Niche Composition per FOV ---
fov_composition <- df_all %>%
  group_by(FullInfo, Condition, Tissue_Niche) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(FullInfo) %>%
  mutate(Fraction = Count / sum(Count)) %>%
  ungroup()

# Ensure the Condition order is preserved for the X-axis sorting
# This keeps CTRL, D1, D2, D3 grouped together in the bar plot
fov_composition <- fov_composition %>%
  arrange(factor(Condition, levels = c("CTRL", "D1", "D2", "D3")), FullInfo) %>%
  mutate(FullInfo = factor(FullInfo, levels = unique(FullInfo)))

# --- Optional: Facet by Condition ---
# If you have 35+ FOVs, it is often better to facet them so they are readable
lower_supp_plot <- ggplot(fov_composition, aes(x = FullInfo, y = Fraction, fill = Tissue_Niche)) +
  geom_bar(stat = "identity", position = "stack", color = "black", 
           linewidth = 0.05) +
  facet_grid(~Condition, scales = "free_x", space = "free_x") +
  ggtitle("Cellular niches across acquired tissue regions") +
  scale_fill_manual(values = niche_colors) +
  theme_classic() +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  theme_classic() +
  labs(title = "Niche Composition across Conditions",
       y = "Proportion of Niche Area", 
       fill = "Tissue Niche") +
  theme(
    # --- REMOVE THE BOX ON TOP ---
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    
    # Cleaning up the X-axis
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    
    # General Aesthetics
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 10, face = "bold"),
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

lower_supp_plot
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-24-1.png" alt="" width="100%" style="display: block; margin: auto;" />


``` r
ggarrange(upper_supp_plot, "NONE", lower_supp_plot, ncol = 1, nrow = 3, labels = c("", "", "C"), heights = c(10, 1, 10))
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-25-1.png" alt="" width="100%" style="display: block; margin: auto;" />


``` r
# --- 1. Calculate CellType Proportions per Niche per FOV ---
cell_comp_per_niche <- df_all %>%
  mutate(CellType = factor(CellType, levels = c(
        "NK cells/ILC1s",
    "ILC2s",
    "ILC3s",
    "T cytotox cells",
    "T helper cells",
    "B cells & Plasma cells",
    "Myeloid cells",
    "LYVE1 CD31 vessels",
    "LYVE1 CD90 Lymphatics",
    "EMCN CD31 Blood vessels",
    "Epithelia"
  ))) %>%
  # Group by FOV, its Condition, the Niche, and the CellType
  group_by(FullInfo, Condition, Tissue_Niche, CellType) %>%
  summarise(Count = n(), .groups = "drop") %>%
  # Now group by FOV and Niche to calculate the internal percentage
  group_by(FullInfo, Tissue_Niche) %>%
  mutate(Fraction = Count / sum(Count)) %>%
  ungroup()

# 2. Ensure ordering for the X-axis (Grouping by Condition)
cell_comp_per_niche <- cell_comp_per_niche %>%
  arrange(factor(Condition, levels = c("CTRL", "D1", "D2", "D3")), FullInfo) %>%
  mutate(FullInfo = factor(FullInfo, levels = unique(FullInfo)))

# --- 3. Create the Multi-Faceted Stacked Bar Plot ---
ggplot(cell_comp_per_niche, aes(x = FullInfo, y = Fraction, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", color = "black", 
           linewidth = 0.05) +
  # Rows = Niches, Columns = Conditions
  facet_grid(Tissue_Niche ~ Condition, scales = "free_x", space = "fixed", 
             switch = "y") +
  scale_fill_manual(values = ColorsCellType) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  theme_classic() +
  labs(title = "Cell Type Composition per Tissue Niche",
       subtitle = "Each bar represents one FOV partitioned by cell type percentages",
       y = "Cell Type Proportion (%)", 
       x = "FOVs",
       fill = "AL3 cell types") +
  theme(
    # --- Clean up Facet Labels (No Boxes) ---
    strip.background = element_blank(),
    strip.text.x = element_text(face = "bold", size = 12),
    strip.text.y = element_text(face = "bold", size = 10, angle = 0), # Horizontal niche labels
    
    # Cleaning up the X-axis (FOVs)
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),

    # General Aesthetics
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 11, face = "bold"),
    axis.title.x = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
    panel.spacing.y = unit(1, "lines") # Add a bit of space between niche rows
  )
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-26-1.png" alt="" width="100%" style="display: block; margin: auto;" />

### Niche abundance


``` r
# ========================================================================
# PART 3: COMPARING CONDITIONS
# ========================================================================

# 1. Calculate the frequency of each niche per FOV
niche_abundance <- df_all %>%
  group_by(Filename, Condition, Tissue_Niche) %>%
  summarise(Cells_in_Niche = n(), .groups = "drop") %>%
  group_by(Filename) %>%
  mutate(Total_Cells = sum(Cells_in_Niche),
         Frequency = (Cells_in_Niche / Total_Cells) * 100) %>%
  ungroup()

# --- 1. SETTINGS & PREPARATION ---
# Ensure your Conditions are in the correct order for the x-axis
# Adjust levels to match your actual data (e.g., "CTRL", "D1", "D3")
target_conditions <- c("CTRL", "D1", "D2", "D3") 
niche_abundance$Condition <- factor(niche_abundance$Condition, levels = target_conditions)

all_niches <- unique(niche_abundance$Tissue_Niche)
niche_plot_list <- list()

# --- 2. LOOP THROUGH EACH NICHE ---
for (target_niche in all_niches) {
  
  # Step A: Filter data for this specific niche
  stats_df <- niche_abundance %>%
    filter(Tissue_Niche == target_niche)
  
  if (nrow(stats_df) == 0) next
  
  # Step B: Statistical Testing (Pairwise Wilcoxon)
  # This compares all conditions against each other
  stat.test <- tryCatch({
    stats_df %>%
      wilcox_test(Frequency ~ Condition) %>%
      add_significance() %>%
      add_xy_position(x = "Condition", step.increase = 0.1)
  }, error = function(e) NULL)
  
  # Step C: Calculate n labels (Number of FOVs per Condition)
  global_max <- max(stats_df$Frequency, na.rm = TRUE)
  
  n_df <- stats_df %>%
    group_by(Condition) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(
      Label = paste0("n=", n),
      # Position labels at the very top of the plot
      Y_pos = global_max * 1.1 
    )
  
  # Step D: Build the Plot
  p <- ggplot(stats_df, aes(x = Condition, y = Frequency)) +
    geom_boxplot(outlier.colour = NA, alpha = 0.5, fill = "white") +
    geom_beeswarm(aes(color = Condition), size = 1.5, cex = 2) + 
    # Add n= labels
    # geom_text(data = n_df, aes(x = Condition, y = 95, label = Label), 
    #           size = 3, fontface = "italic", inherit.aes = FALSE) +
    ggtitle(target_niche) +
    # Expand Y-axis for significance brackets
    scale_y_continuous(expand = c(0, 0.5), limits = c(0,90))+
    theme_classic() + 
    theme(
      axis.text.x = element_text(angle = 45, size = 9, face = "bold" , 
                                     hjust = 1),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 9),
      axis.title.y = element_text(size = 9),
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
      plot.margin = margin(0, 0.2, 0.3, 0.2, "cm"),
      legend.position = "none") +
    ylab("% of Total Cells per FOV") +
    scale_color_manual(values = cols_treat)

  # Step E: Add Stat Brackets
  if (!is.null(stat.test) && nrow(stat.test) > 0) {
    # Dynamically pick significance column
    sig_label <- ifelse("p.adj.signif" %in% colnames(stat.test), "p.adj.signif", "p.signif")
    
    p <- p + stat_pvalue_manual(stat.test, 
                                label = sig_label, 
                                hide.ns = TRUE, 
                                tip.length = 0.01,
                                step.increase = 0.05)
  }
  
  # Save to list
  niche_plot_list[[as.character(target_niche)]] <- p
}

# --- 3. COMBINE INTO GRID ---
if (length(niche_plot_list) > 0) {
  final_niche_grid <- wrap_plots(niche_plot_list, ncol = 4) +
    plot_annotation(
      title = "Cellular Niche Abundance across Conditions",
      subtitle = "Each point represents one FOV; Stats: Pairwise Wilcoxon Test",
      theme = theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
      )
    )
  
  print(final_niche_grid)
}
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-27-1.png" alt="" width="100%" style="display: block; margin: auto;" />


``` r
# --- 3. COMBINE INTO GRID ---
if (length(niche_plot_list) > 0) {
  final_niche_grid <- wrap_plots(niche_plot_list, ncol = 2) +
    plot_annotation(
      title = "Cellular Niche Abundance\nacross Conditions\n",
      theme = theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
        )
    )
  
  print(final_niche_grid)
}
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-28-1.png" alt="" width="100%" style="display: block; margin: auto;" />


``` r
upper_fig <- ggarrange(plot_comp_al3, final_niche_grid, ncol = 2, nrow = 1, 
          widths = c(3.5, 4) , labels = c("A", "B"))

upper_fig
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-29-1.png" alt="" width="100%" style="display: block; margin: auto;" />


``` r
# Calculate the frequency of each CellType within each Tissue_Niche per FOV
niche_composition <- df_all %>%
  group_by(Filename, Condition, Tissue_Niche, CellType) %>%
  summarise(Cell_Count = n(), .groups = "drop") %>%
  # Ensure all cell types are represented in every niche/FOV, even if count is 0
  complete(nesting(Filename, Condition, Tissue_Niche), CellType, fill = list(Cell_Count = 0)) %>%
  group_by(Filename, Tissue_Niche) %>%
  mutate(
    Total_Niche_Cells = sum(Cell_Count),
    Frequency = (Cell_Count / Total_Niche_Cells) * 100
  ) %>%
  ungroup() %>%
  # Remove NaNs if a niche had 0 cells total in an FOV
  filter(!is.na(Frequency))

niche_composition$Condition <- factor(niche_composition$Condition, levels = c("CTRL", "D1", "D2", "D3"))


all_niches <- unique(niche_composition$Tissue_Niche)
all_cell_types <- unique(niche_composition$CellType)
niche_plots_all <- list()

for (current_niche in all_niches) {
  
  message("Processing Niche: ", current_niche)
  niche_plots <- list()
  
  # Filter data for the specific niche
  niche_data <- niche_composition %>% filter(Tissue_Niche == current_niche)
  
  for (current_cell in all_cell_types) {
    
    stats_df <- niche_data %>% filter(CellType == current_cell)
    
    # Skip if there's no data at all for this cell type in this niche
    if(nrow(stats_df) == 0 || sum(stats_df$Frequency) == 0) next
    
    # Statistical Testing (Pairwise)
    stat.test <- tryCatch({
      stats_df %>%
        wilcox_test(Frequency ~ Condition) %>%
        add_significance() %>%
        add_xy_position(x = "Condition")
    }, error = function(e) NULL)
    
    # Calculate Max for Labeling
    y_max <- max(stats_df$Frequency, na.rm = TRUE)
    
    # Build the individual plot
    p <- ggplot(stats_df, aes(x = Condition, y = Frequency)) +
      geom_boxplot(outlier.colour = NA, fill = "white", alpha = 0.5) +
      geom_beeswarm(aes(color = Condition), size = 1, cex = 1.5, alpha = 0.6) +
      ggtitle(gsub("LYVE1 CD90 ", "", current_cell)) +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.4))) +
      scale_color_manual(values = cols_treat) +
      theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, size = 9, face = "bold" , 
                                     hjust = 1),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 9),
      axis.title.y = element_text(size = 9),
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
      plot.margin = margin(0, 0.2, 0.3, 0.2, "cm"),
      legend.position = "none") +
      NoLegend()
    
    # Add Stats
    if (!is.null(stat.test) && nrow(stat.test) > 0) {
      sig_col <- if("p.adj.signif" %in% colnames(stat.test)) "p.adj.signif" else "p.signif"
      p <- p + stat_pvalue_manual(stat.test, label = sig_col, hide.ns = TRUE, 
                                  tip.length = 0.01, step.increase = 0.1)
    }
    
    niche_plots[[current_cell]] <- p
    
    plot_name <- paste0(current_niche, "_", current_cell)
    niche_plots_all[[plot_name]] <- p
  }
  
  # Assemble the Grid for this Niche
  if (length(niche_plots) > 0) {
    combined_grid <- wrap_plots(niche_plots, ncol = 4) +
      plot_annotation(
        title = paste("Niche Composition:", gsub("LYVE1 CD90 ", "", current_niche)),
        subtitle = "Y-axis: % of total cells within this specific niche; Stats: Pairwise Wilcoxon",
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
      )
    
    print(combined_grid)
  }
}
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-30-1.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-30-2.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-30-3.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-30-4.png" alt="" width="100%" style="display: block; margin: auto;" />


``` r
middle_fig <- ggarrange(niche_plots_all[["Mixed Myeloid/LEC niche_ILC2s"]],
          niche_plots_all[["Mixed Myeloid/LEC niche_Myeloid cells"]],
          niche_plots_all[["Mixed Myeloid/LEC niche_LYVE1 CD90 Lymphatics"]], 
          ncol = 3, nrow = 1)

middle_fig <- annotate_figure(middle_fig, 
                              top = text_grob("Cellular abundance in Mixed Myeloid/LEC niche\n", 
                                              color = "black", face = "bold", 
                                              size = 12)
                              )

middle_fig
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-31-1.png" alt="" width="100%" style="display: block; margin: auto;" />


``` r
top_figure <- ggarrange(upper_fig, middle_fig, nrow = 2, ncol = 1, 
          heights = c(4, 2), labels = c("", "C"))

top_figure
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-32-1.png" alt="" width="100%" style="display: block; margin: auto;" />

### Spatial distribution of identified niches


``` r
library(ggplot2)
library(dplyr)
library(patchwork)

# ========================================================================
# PART 4: SPATIAL MAPPING (36 Plot Grid Assembly)
# ========================================================================

# 1. Define the Niche colors (Matches your specific clusters)
niche_colors <- c(
  "Epithelial niche" = "#117733", 
  "Mixed Myeloid/LEC niche" = "#DDCC77", 
  "Blood endothelial niche" = "#882255", 
  "Mixed BPC/BEC" = "#332288"
)

# 2. Grab and Sort exactly 35 FOVs by Condition
example_fovs_raw <- unique(df_all$Filename)[1:35] 

fov_order_df <- df_all %>%
  dplyr::filter(Filename %in% example_fovs_raw) %>%
  dplyr::select(Filename, Condition) %>%
  dplyr::distinct() %>%
  dplyr::arrange(factor(Condition, levels = c("CTRL", "D1", "D2", "D3")), Filename)

example_fovs <- fov_order_df$Filename
spatial_niche_plots <- list()

# 3. LOOP: Generate individual FOV maps
for(fov in example_fovs) {
  plot_data <- df_all %>% filter(Filename == fov)
  cond_label <- unique(plot_data$Condition)
  plot_title <- fov
  
  p <- ggplot(plot_data, aes(x = Loc_X, y = Loc_Y, color = Tissue_Niche)) +
    geom_point(size = 1, alpha = 0.5, shape = 19) + # Smaller points for high-density grid
    scale_color_manual(values = niche_colors) +
    scale_y_reverse() + # Standard spatial orientation
    coord_fixed() +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "black", color = NA),
      plot.background = element_rect(fill = "black", color = NA),
      plot.title = element_text(color = "white", hjust = 0.5, face = "bold", size = 7),
      legend.position = "none"
    ) +
    labs(title = plot_title)
  
  spatial_niche_plots[[fov]] <- p
}

# 4. THE LEGEND PLOT (The 36th Plot)
n_niches <- length(niche_colors)
text_legend_df <- data.frame(
  Niche = names(niche_colors),
  Y_pos = seq(n_niches, 1),
  X_pos = rep(0, n_niches)
)

legend_plot <- ggplot(text_legend_df, aes(x = X_pos, y = Y_pos)) +
  geom_point(aes(color = Niche), size = 5, shape = 15) + 
  geom_text(aes(label = Niche), color = "white", hjust = 0, nudge_x = 0.1, size = 4, fontface = "bold") +
  scale_color_manual(values = niche_colors) +
  xlim(-0.1, 1.5) + 
  ylim(0, n_niches + 1) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    legend.position = "none"
  )

# 5. ASSEMBLE THE 6x6 GRID
# Combine the 35 FOV plots + 1 legend plot
final_plot_list <- c(spatial_niche_plots, list(legend = legend_plot))

niche_map_grid <- wrap_plots(final_plot_list, ncol = 9) +
  plot_annotation(
    title = "Spatial Niche Atlas: Tissue Organization across Timepoints",
    theme = theme(
      plot.background = element_rect(fill = "black", color = NA),
      plot.title = element_text(color = "white", size = 20, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(color = "gray70", size = 14, hjust = 0.5),
      plot.margin = margin(1, 1, 1, 1)
    )
  )

# Render the atlas
print(niche_map_grid)
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-33-1.png" alt="" width="100%" style="display: block; margin: auto;" />


``` r
library(ggplot2)
library(dplyr)
library(patchwork)

# ========================================================================
# PART 4: SPATIAL MAPPING (3x3 CTRL vs 3x3 D3)
# ========================================================================

# 2. Function to generate a list of 9 plots for a specific condition
get_niche_plots <- function(data, target_condition, n_fovs = 9) {
  fovs <- data %>%
    filter(Condition == target_condition) %>%
    pull(Filename) %>%
    unique() %>%
    head(n_fovs)
  
  plot_list <- list()
  
  for(fov in fovs) {
    plot_data <- data %>% filter(Filename == fov)
    
    p <- ggplot(plot_data, aes(x = Loc_X, y = Loc_Y, color = Tissue_Niche)) +
      geom_point(size = 0.5, alpha = 0.6, shape = 19) + 
      scale_color_manual(values = niche_colors) +
      scale_y_reverse() + 
      coord_fixed() +
      theme_void() +
      theme(
        panel.background = element_rect(fill = "black", color = NA),
        plot.background = element_rect(fill = "black", color = NA),
        plot.title = element_text(color = "white", hjust = 0.5, size = 7),
        legend.position = "none"
      ) +
      labs(title = fov)
    
    plot_list[[fov]] <- p
  }
  return(plot_list)
}

# 3. Generate the plot lists
ctrl_plots <- get_niche_plots(df_all, "CTRL")
d3_plots   <- get_niche_plots(df_all, "D3")

# 4. Create the Legend Plot (HORIZONTAL 1 ROW, 4 COLUMNS)
# We spread the X coordinates out (1, 4, 7, 10) to give the text labels room
legend_df <- data.frame(
  Niche = names(niche_colors),
  Y = 1,
  X = c(1, 4, 7, 10) 
)

legend_plot <- ggplot(legend_df, aes(X, Y)) +
  geom_point(aes(color = Niche), size = 6, shape = 15) + 
  geom_text(aes(label = Niche), 
            color = "white", 
            hjust = 0, 
            nudge_x = 0.3, 
            size = 9/.pt, # Matches theme size 12pt
            fontface = "bold") +
  scale_color_manual(values = niche_colors) +
  # Expand limits so the last text label isn't cut off
  xlim(0.5, 13) + 
  theme_void() +
  theme(
    plot.background = element_rect(fill = "black", color = NA), 
    legend.position = "none",
    plot.margin = margin(t = 1, b = 1)
  )


# ========================================================================
# PART 4: SPATIAL MAPPING (Fixed Nested Titles)
# ========================================================================

# [Niche colors and get_niche_plots function remain the same as your code]

# 5. Assemble the Grids using wrap_elements to protect nested titles
grid_ctrl <- wrap_elements(full = (
  wrap_plots(ctrl_plots, ncol = 3) + 
    plot_annotation(
      title = "________________________  CTRL  ________________________",
      theme = theme(
        plot.background = element_rect(fill = "black", color = NA),
        plot.title = element_text(color = "white", size = 10, hjust = 0.5, margin = margin(b = 1))
      )
    ) & theme(plot.background = element_rect(fill = "black", color = NA))
))

grid_d3 <- wrap_elements(full = (
  wrap_plots(d3_plots, ncol = 3) + 
    plot_annotation(
      title = "______________________  IL-33 day 3  ______________________",
      theme = theme(
        plot.background = element_rect(fill = "black", color = NA),
        plot.title = element_text(color = "white", size = 10, hjust = 0.5, margin = margin(b = 1))
      )
    ) & theme(plot.background = element_rect(fill = "black", color = NA))
))

# 6. Assemble the Final Atlas
# The heights c(10, 1) usually works better to keep the legend readable
final_atlas <- (grid_ctrl | grid_d3) / legend_plot + 
  plot_layout(heights = c(10, 0.2)) +
  plot_annotation(
    title = "Spatial distribution of tissue niches",
    theme = theme(
      plot.background = element_rect(fill = "black", color = NA),
      plot.title = element_text(color = "white", size = 12, face = "bold", hjust = 0.5, margin = margin(t = 1, b = 1))
    )
  ) & theme(plot.background = element_rect(fill = "black", color = NA))

print(final_atlas)
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-34-1.png" alt="" width="100%" style="display: block; margin: auto;" />

Proliferation


``` r
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(rstatix)
library(patchwork)
library(ggpubr)

# --- 1. DATA PREPARATION ---
# Filter for target cells and calculate Mean Ki67 per FOV
ki67_df <- df_all %>%
  filter(CellType %in% c("NK cells/ILC1s", "ILC2s", "ILC3s")) %>%
  group_by(Filename, Condition, CellType) %>%
  summarise(Mean_Ki67 = mean(Ki67, na.rm = TRUE), .groups = "drop")

# Ensure Condition order
ki67_df$Condition <- factor(ki67_df$Condition, levels = c("CTRL", "D1", "D2", "D3"))

# --- 2. STATISTICAL TESTING ---
# Pairwise comparisons between conditions for each cell type
stat.test <- ki67_df %>%
  group_by(CellType) %>%
  wilcox_test(Mean_Ki67 ~ Condition) %>%
  add_significance() %>%
  add_xy_position(x = "Condition", dodge = 0.8)

# --- 3. PLOTTING LOOP ---
ki67_plots <- list()
target_cells <- c("NK cells/ILC1s", "ILC2s", "ILC3s")

for (ct in target_cells) {
  
  # Filter data for specific cell type
  stats_df <- ki67_df %>% filter(CellType == ct)
  
  # Filter stats for specific cell type
  curr_stat <- stat.test %>% filter(CellType == ct)
  
  # Detect correct p-value label
  sig_label <- if("p.adj.signif" %in% colnames(curr_stat)) "p.adj.signif" else "p.signif"
  
  # Build Plot
  p <- ggplot(stats_df, aes(x = Condition, y = Mean_Ki67)) +
    geom_boxplot(outlier.colour = NA, fill = "white", alpha = 0.5) +
    geom_beeswarm(aes(color = Condition), size = 1.5, cex = 2, alpha = 0.7) +
    
    # Add stats
    stat_pvalue_manual(curr_stat, 
                       label = sig_label, 
                       hide.ns = TRUE, 
                       step.increase = 0.1, 
                       size = 4) +
    
    # Aesthetic settings
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.4))) +
    scale_color_manual(values = cols_treat) +
    ggtitle(ct) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, size = 9, face = "bold", hjust = 1),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 9),
      axis.title.y = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
      plot.margin = margin(0.2, 0.2, 0.3, 0.2, "cm"),
      legend.position = "none"
    ) +
    ylab("Mean Ki67 Intensity per FOV")

  ki67_plots[[ct]] <- p
}

# --- 4. ASSEMBLE ---
final_ki67_plot <- wrap_plots(ki67_plots, ncol = 3) +
  plot_annotation(
    title = "Ki67 Expression Proliferation Analysis",
    subtitle = "Points represent mean intensity per FOV; Wilcoxon pairwise test",
    theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(hjust = 0.5))
  )

print(final_ki67_plot)
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-35-1.png" alt="" width="100%" style="display: block; margin: auto;" />


``` r
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)
library(rstatix)
library(patchwork)
library(ggpubr)

# --- 1. SETTINGS & PARAMETERS ---
target_markers <- c("Ki67", "GATA3eGFP", "EOMES", "TBET", "MHCII", "KLRG1")
target_cells   <- c("NK cells/ILC1s", "ILC2s", "ILC3s", "T helper cells", "T cytotox cells")

# --- 2. DATA AGGREGATION ---
# Calculate mean intensity for each marker per cell type per niche per FOV
niche_marker_data <- df_all %>%
  filter(CellType %in% target_cells) %>%
  group_by(Filename, Condition, Tissue_Niche, CellType) %>%
  summarise(across(all_of(target_markers), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  # Pivot markers to long format for easier looping
  pivot_longer(cols = all_of(target_markers), names_to = "Marker", values_to = "Intensity")

niche_marker_data$Condition <- factor(niche_marker_data$Condition, levels = c("CTRL", "D1", "D2", "D3"))
all_niches <- unique(niche_marker_data$Tissue_Niche)

# --- 3. PLOTTING LOOP ---
for (current_niche in all_niches) {
  
  message("Processing Niche: ", current_niche)
  
  # Filter for the current niche
  niche_sub <- niche_marker_data %>% filter(Tissue_Niche == current_niche)
  
  # This list will hold plots in the order: Row 1 (all cells), Row 2 (all cells), etc.
  grid_plots <- list()
  
  for (marker in target_markers) {
    for (cell in target_cells) {
      
      # Subset data for this specific panel
      stats_df <- niche_sub %>% filter(Marker == marker, CellType == cell)
      
      # Check if we have enough data to plot
      if(nrow(stats_df) < 2 || all(is.na(stats_df$Intensity))) {
        # Create an empty placeholder plot to maintain grid alignment
        p <- ggplot() + theme_void() 
      } else {
        
        # Statistical Testing (Pairwise Wilcoxon)
        stat.test <- tryCatch({
          stats_df %>%
            wilcox_test(Intensity ~ Condition) %>%
            add_significance() %>%
            add_xy_position(x = "Condition")
        }, error = function(e) NULL)
        
        # Build Panel
        p <- ggplot(stats_df, aes(x = Condition, y = Intensity)) +
          geom_boxplot(outlier.colour = NA, fill = "white", alpha = 0.5) +
          geom_beeswarm(aes(color = Condition), size = 0.8, cex = 1.2, alpha = 0.6) +
          scale_y_continuous(expand = expansion(mult = c(0.1, 0.4))) +
          scale_color_manual(values = cols_treat) +
          # Title logic: Only show CellType on top row and Marker on left column (optional)
          # Here we keep both for clarity in large grids
          labs(title = marker, subtitle = cell, 
               y = "Mean Intensity") +
          theme_classic() +
          theme(
            axis.text.x = element_text(angle = 45, size = 7, face = "bold", hjust = 1),
            axis.title.x = element_blank(),
            axis.text.y = element_text(size = 7),
            axis.title.y = element_text(size = 7),
            plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), 
            plot.subtitle = element_text(hjust = 0.5, size = 8),
            plot.margin = margin(2, 2, 2, 2, "pt"),
            legend.position = "none"
          )
        
        # Add Stats
        if (!is.null(stat.test) && nrow(stat.test) > 0) {
          sig_col <- if("p.adj.signif" %in% colnames(stat.test)) "p.adj.signif" else "p.signif"
          p <- p + stat_pvalue_manual(stat.test, label = sig_col, hide.ns = TRUE, 
                                      tip.length = 0.01, step.increase = 0.1, size = 2.5)
        }
      }
      grid_plots[[paste(marker, cell)]] <- p
    }
  }
  
  # --- 4. ASSEMBLE NICHE GRID ---
  # ncol = number of cell types (columns)
  if (length(grid_plots) > 0) {
    combined_grid <- wrap_plots(grid_plots, ncol = length(target_cells)) +
      plot_annotation(
        title = paste("Niche Protein Profile:", current_niche),
        subtitle = "Rows: Markers | Columns: Cell Types | Stats: Wilcoxon vs CTRL",
        theme = theme(
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 11, hjust = 0.5)
        )
      )
    
    # Recommend saving as a large PDF or printing
    # {r fig.height=12, fig.width=15}
    print(combined_grid)
  }
}
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-36-1.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-36-2.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-36-3.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-36-4.png" alt="" width="100%" style="display: block; margin: auto;" />

# Combine plots for figure


``` r
ggarrange(top_figure, "NONE", final_atlas, ncol = 1, nrow = 3, heights = c(6.5, 0.1, 5.2), 
          labels = c("", "", "D"), label.y = 1.06)
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-37-1.png" alt="" width="100%" style="display: block; margin: auto;" />


``` r
coenrichment <- ggarrange(final_plot_ILC2s, "NONE", coenrichment_fig_ann,  
          ncol = 1, nrow = 3, heights = c(1, 0.1, 0.9),
          labels = c("A", "", "B"))

coenrichment
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-38-1.png" alt="" width="100%" style="display: block; margin: auto;" />

``` r
plot_cin
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-38-2.png" alt="" width="100%" style="display: block; margin: auto;" />

``` r
# spiat <- ggarrange(dist_lymph, plot_cin,
#           ncol = 1, nrow = 2, heights = c(3, 1),
#           labels = c("", "E", "F", "G", "H", "I"))+
#   theme(plot.margin = margin(0, 0.1, 0, 0, "cm"))
# 
# ggarrange(coenrichment, spiat, 
#           ncol = 2, nrow = 1, widths = c(5.5, 3.5)
#           # labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I")
#           )+
#   theme(plot.margin = margin(0, 0.1, 0, 0, "cm"))
```


``` r
ggarrange(coenrichment, "NONE", plot_cin,  
          ncol = 1, nrow = 3, heights = c(5, 0.5, 6),
          labels = c("", "", ""))
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-39-1.png" alt="" width="100%" style="display: block; margin: auto;" />

Additional plots


``` r
# fine tune the co-enrichment plot
#  -----------------------------------------
interactions <- c("LYVE1 CD90 Lymphatics--T cytotox cells")
unic <- "LYVE1 CD90 Lymphatics"
background_image <- backgroundlist[["LYVE1 CD90 Lymphatics"]]
interaction_celltypes <- NULL
for(samp in vr_list_names){
    cur_cell_proximities <- cell_proximities_list[[samp]]
    cur_cell_proximities <- cur_cell_proximities[cur_cell_proximities$unified_int %in% interactions,]
    sample <- unique(metadatax$FullInfo[metadatax$Sample==samp])
    if(nrow(cur_cell_proximities) > 0 & sample != "20210906_FOV3_D3"){
      interaction_celltypes <- rbind(interaction_celltypes,
                                     data.frame(cur_cell_proximities[cur_cell_proximities$unified_int %in% interactions,], 
                                                experiment = strsplit(sample, split = "_")[[1]][1], fov = strsplit(sample, split = "_")[[1]][2], condition = strsplit(sample, split = "_")[[1]][3]))
    }
  }
interaction_celltypes$p.adj <- ifelse(interaction_celltypes$enrichm > 0, interaction_celltypes$p.adj_higher, interaction_celltypes$p.adj_lower)
  # plot test results
  # sig_label <- as.character(ifelse(interaction_celltypes$p.adj < 0.1, paste0("p=",round(interaction_celltypes$p.adj,3)), ""))
  sig_label <- as.character(ifelse(interaction_celltypes$p.adj < 0.1, paste0("*"), ""))
    # print(sig_label)
plot_coenrichment <- ggplot(interaction_celltypes, aes(x = condition, y = enrichm, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge2(width=0.9, preserve = "single")) +
  facet_grid(.~condition, scales = "free_x") +
  geom_text(aes(label=sig_label), position=position_dodge2(width=0.9, preserve = "single"), angle = 90, hjust = -0.02, size = 4) +
  ylim(-2,2.1)+
  NoLegend()+
  theme_classic2()+
  scale_fill_manual(values = cols_fov, name = "") +
  ggtitle(gsub("LYVE1 CD90 ", "", interactions)) +
  theme(axis.text.x = element_text(#angle = 50,
                                   vjust = 1, size = 12, hjust = 0.5, face = "bold"
                                   ),
        axis.text.y = element_text(hjust = 0.5, size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size =14, hjust = 0.5),
        plot.margin = margin(0.5, 1, 0.5, 1, "cm"),
        legend.position = "none",
        strip.background=element_blank(),
        strip.background.x= element_blank(),
        strip.text.x = element_text(size = 1, color = "white"),
        panel.grid.major.y = element_line())+
  NoLegend()+  
  ylab("Enrichment")

plot_coenrichment
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-40-1.png" alt="" width="100%" style="display: block; margin: auto;" />

``` r
# fine tune the co-enrichment plot
#  -----------------------------------------
interactions <- c("LYVE1 CD90 Lymphatics--Myeloid cells")
unic <- "LYVE1 CD90 Lymphatics"
background_image <- backgroundlist[["LYVE1 CD90 Lymphatics"]]
interaction_celltypes <- NULL
for(samp in vr_list_names){
    cur_cell_proximities <- cell_proximities_list[[samp]]
    cur_cell_proximities <- cur_cell_proximities[cur_cell_proximities$unified_int %in% interactions,]
    sample <- unique(metadatax$FullInfo[metadatax$Sample==samp])
    if(nrow(cur_cell_proximities) > 0 & sample != "20210906_FOV3_D3"){
      interaction_celltypes <- rbind(interaction_celltypes,
                                     data.frame(cur_cell_proximities[cur_cell_proximities$unified_int %in% interactions,], 
                                                experiment = strsplit(sample, split = "_")[[1]][1], fov = strsplit(sample, split = "_")[[1]][2], condition = strsplit(sample, split = "_")[[1]][3]))
    }
  }
interaction_celltypes$p.adj <- ifelse(interaction_celltypes$enrichm > 0, interaction_celltypes$p.adj_higher, interaction_celltypes$p.adj_lower)
  # plot test results
  # sig_label <- as.character(ifelse(interaction_celltypes$p.adj < 0.1, paste0("p=",round(interaction_celltypes$p.adj,3)), ""))
  sig_label <- as.character(ifelse(interaction_celltypes$p.adj < 0.1, paste0("*"), ""))
    # print(sig_label)
plot_coenrichment <- ggplot(interaction_celltypes, aes(x = condition, y = enrichm, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge2(width=0.9, preserve = "single")) +
  facet_grid(.~condition, scales = "free_x") +
  geom_text(aes(label=sig_label), position=position_dodge2(width=0.9, preserve = "single"), angle = 90, hjust = -0.02, size = 4) +
  ylim(-2,2.1)+
  NoLegend()+
  theme_classic2()+
  scale_fill_manual(values = cols_fov, name = "") +
  ggtitle(gsub("LYVE1 CD90 ", "", interactions)) +
  theme(axis.text.x = element_text(#angle = 50,
                                   vjust = 1, size = 12, hjust = 0.5, face = "bold"
                                   ),
        axis.text.y = element_text(hjust = 0.5, size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size =14, hjust = 0.5),
        plot.margin = margin(0.5, 1, 0.5, 1, "cm"),
        legend.position = "none",
        strip.background=element_blank(),
        strip.background.x= element_blank(),
        strip.text.x = element_text(size = 1, color = "white"),
        panel.grid.major.y = element_line())+
  NoLegend()+  
  ylab("Enrichment")

plot_coenrichment
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-40-2.png" alt="" width="100%" style="display: block; margin: auto;" />

Frequency of cell types per FOV/condition:


``` r
df_freq_cells <- metadatax %>%
  select(Condition, CellType, FullInfo) %>%
  janitor::tabyl(Condition, CellType) %>%
  janitor::adorn_percentages() %>%
  janitor::adorn_totals(c('row', 'col')) %>%
  janitor::adorn_pct_formatting(2)
```

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
## [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] dbscan_1.2.3          patchwork_1.3.2       ggplotify_0.1.3       tidyr_1.3.1           circlize_0.4.17       ComplexHeatmap_2.26.1 stringr_1.6.0         ggbeeswarm_0.7.3      readr_2.1.6           ggpubr_0.6.2          ggplot2_4.0.1         VoltRon_0.2.3         Seurat_5.3.1          Giotto_4.2.2          GiottoClass_0.4.10    rlang_1.1.6           rstatix_0.7.3         dplyr_1.1.4           SeuratObject_5.2.0    sp_2.2-0             
## 
## loaded via a namespace (and not attached):
##   [1] fs_1.6.6                    matrixStats_1.5.0           spatstat.sparse_3.1-0       bitops_1.0-9                lubridate_1.9.4             EBImage_4.52.0              doParallel_1.0.17           httr_1.4.8                  RColorBrewer_1.1-3          tools_4.5.2                 sctransform_0.4.2           backports_1.5.0             utf8_1.2.6                  R6_2.6.1                    lazyeval_0.2.2              uwot_0.2.4                  GetoptLong_1.1.0            withr_3.0.2                 gridExtra_2.3               GiottoUtils_0.2.5           progressr_0.18.0            cli_3.6.5                   Biobase_2.70.0              Cairo_1.7-0                 spatstat.explore_3.5-3      fastDummies_1.7.5           shinyjs_2.1.1               labeling_0.4.3              sass_0.4.10                 S7_0.2.0                    spatstat.data_3.1-9         ggridges_0.5.7              pbapply_1.7-4               yulab.utils_0.2.4           dichromat_2.0-0.1           parallelly_1.45.1           rstudioapi_0.18.0           gridGraphics_0.5-1          shape_1.4.6.1               generics_0.1.4              vroom_1.7.0                 gtools_3.9.5               
##  [43] ica_1.0-3                   spatstat.random_3.4-2       car_3.1-5                   Matrix_1.7-4                S4Vectors_0.48.0            abind_1.4-8                 terra_1.8-93                lifecycle_1.0.5             scatterplot3d_0.3-45        yaml_2.3.10                 snakecase_0.11.1            carData_3.0-6               SummarizedExperiment_1.40.0 gplots_3.3.0                SparseArray_1.10.1          Rtsne_0.17                  promises_1.5.0              crayon_1.5.3                miniUI_0.1.2                lattice_0.22-7              beachmat_2.26.0             cowplot_1.2.0               magick_2.9.0                pillar_1.11.1               knitr_1.51                  GenomicRanges_1.62.0        rjson_0.2.23                future.apply_1.20.2         codetools_0.2-20            glue_1.8.0                  spatstat.univar_3.1-4       data.table_1.17.8           vctrs_0.6.5                 png_0.1-8                   ids_1.0.1                   spam_2.11-1                 gtable_0.3.6                cachem_1.1.0                xfun_0.56                   S4Arrays_1.10.0             mime_0.13                   tidygraph_1.3.1            
##  [85] Seqinfo_1.0.0               survival_3.8-3              SingleCellExperiment_1.32.0 iterators_1.0.14            bluster_1.20.0              rgl_1.3.34                  fitdistrplus_1.2-6          ROCR_1.0-12                 colorsGen_1.0.0             nlme_3.1-168                bit64_4.6.0-1               RcppAnnoy_0.0.22            rprojroot_2.1.1             bslib_0.10.0                irlba_2.3.5.1               vipor_0.4.7                 KernSmooth_2.23-26          otel_0.2.0                  colorspace_2.1-2            BiocGenerics_0.56.0         tidyselect_1.2.1            bit_4.6.0                   compiler_4.5.2              BiocNeighbors_2.4.0         DelayedArray_0.36.0         plotly_4.12.0               checkmate_2.3.4             scales_1.4.0                caTools_1.18.3              lmtest_0.9-40               rappdirs_0.3.4              tiff_0.1-12                 SpatialExperiment_1.20.0    digest_0.6.38               goftest_1.2-3               fftwtools_0.9-11            spatstat.utils_3.2-0        rmarkdown_2.30              XVector_0.50.0              htmltools_0.5.8.1           GiottoVisuals_0.2.14        pkgconfig_2.0.3            
## [127] jpeg_0.1-11                 base64enc_0.1-6             MatrixGenerics_1.22.0       fastmap_1.2.0               GlobalOptions_0.1.3         htmlwidgets_1.6.4           shiny_1.13.0                Rvcg_0.25                   farver_2.1.2                jquerylib_0.1.4             zoo_1.8-14                  jsonlite_2.0.0              BiocParallel_1.44.0         BiocSingular_1.26.0         RCurl_1.98-1.17             magrittr_2.0.4              Formula_1.2-5               dotCall64_1.2               RCDT_1.3.0                  Rcpp_1.1.0                  viridis_0.6.5               reticulate_1.44.0           stringi_1.8.7               ggraph_2.2.2                MASS_7.3-65                 plyr_1.8.9                  parallel_4.5.2              listenv_0.10.0              ggrepel_0.9.6               deldir_2.0-4                graphlayouts_1.2.2          splines_4.5.2               tensor_1.5.1                hms_1.1.4                   locfit_1.5-9.12             colorRamp2_0.1.0            igraph_2.2.1                uuid_1.2-2                  spatstat.geom_3.6-0         ggsignif_0.6.4              RcppHNSW_0.6.0              reshape2_1.4.5             
## [169] stats4_4.5.2                ScaledMatrix_1.18.0         evaluate_1.0.5              foreach_1.5.2               tzdb_0.5.0                  tweenr_2.0.3                httpuv_1.6.16               RANN_2.6.2                  purrr_1.2.0                 polyclip_1.10-7             clue_0.3-67                 future_1.69.0               scattermore_1.2             ggforce_0.5.0               janitor_2.2.1               rsvd_1.0.5                  broom_1.0.12                xtable_1.8-8                RSpectra_0.16-2             later_1.4.4                 viridisLite_0.4.2           Polychrome_1.5.4            tibble_3.3.0                beeswarm_0.4.0              memoise_2.0.1               IRanges_2.44.0              cluster_2.1.8.1             timechange_0.3.0            globals_0.19.0              here_1.0.2
```
