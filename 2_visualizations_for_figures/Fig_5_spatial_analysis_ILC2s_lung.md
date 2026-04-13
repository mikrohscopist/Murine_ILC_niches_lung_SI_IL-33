---
title: "Figure 5: ILC2s localize"
author: "Sandy Kroh"
date: "April 13, 2026"
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

## ILC2s localize close to vessels

### B - IF overlays


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

### A - Coenrichment plots


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

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-8-1.png" alt="" width="100%" style="display: block; margin: auto;" />


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

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-9-1.png" alt="" width="100%" style="display: block; margin: auto;" />

### CIN analysis


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

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-1.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-2.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-3.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-4.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-5.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-6.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-7.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-8.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-9.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-10.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-11.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-12.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-13.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-14.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-15.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-16.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-17.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-18.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-19.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-20.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-21.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-22.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-23.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-24.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-25.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-26.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-27.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-28.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-29.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-30.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-31.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-32.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-33.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-34.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-35.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-36.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-37.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-38.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-39.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-40.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-41.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-42.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-43.png" alt="" width="100%" style="display: block; margin: auto;" /><img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-11-44.png" alt="" width="100%" style="display: block; margin: auto;" />


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

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-12-1.png" alt="" width="100%" style="display: block; margin: auto;" />

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
                        `Myeloid cells` = "gold", 
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
## [1] "gold"
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

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-14-1.png" alt="" width="100%" style="display: block; margin: auto;" />

``` r
dist_lymph <- ggarrange(plotlist = my_plot_list[c(1, 4)], nrow = 1, ncol = 2, labels = c("C", "D"))

dist_lymph
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-14-2.png" alt="" width="100%" style="display: block; margin: auto;" />


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

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-15-1.png" alt="" width="100%" style="display: block; margin: auto;" />


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

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-16-1.png" alt="" width="100%" style="display: block; margin: auto;" />


``` r
plot_dist <- ggarrange(dist_lymph, final_plot, nrow = 2, ncol = 1, heights = c(4, 3))

plot_dist
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-17-1.png" alt="" width="100%" style="display: block; margin: auto;" />

Check what is the frequency of ILC2s with a minimum distance to lymphatics less than 5 µm:


``` r
plot_data <- df_dist_ref %>%
  filter(Treatment == "CTRL") %>%
  filter(RefType == "ILC2s") %>%
  filter(NearestType == celltype_of_interest)

percentage_under_5 <- plot_data %>%
  summarise(pct = mean(Distance < 5) * 100) %>%
  pull(pct)

print(paste("Percentage of cells with Distance < 5 under healthy conditions:", round(percentage_under_5, 2), "%"))
```

```
## [1] "Percentage of cells with Distance < 5 under healthy conditions: 1.92 %"
```

``` r
plot_data <- df_dist_ref %>%
  filter(RefType == "ILC2s") %>%
  filter(NearestType == celltype_of_interest)

percentage_under_5 <- plot_data %>%
  summarise(pct = mean(Distance < 5) * 100) %>%
  pull(pct)

print(paste("Percentage of cells with Distance < 5:", round(percentage_under_5, 2), "%"))
```

```
## [1] "Percentage of cells with Distance < 5: 2.66 %"
```

### E - CIN plot



For comparison, check the CIN of :



## Combine plots for figure


``` r
coenrichment <- ggarrange(final_plot_ILC2s, "NONE", coenrichment_fig_ann,  
          ncol = 1, nrow = 3, heights = c(1, 0.1, 0.9),
          labels = c("A", "", "B"))

coenrichment
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-21-1.png" alt="" width="100%" style="display: block; margin: auto;" />

``` r
plot_cin
```

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-21-2.png" alt="" width="100%" style="display: block; margin: auto;" />

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

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-22-1.png" alt="" width="100%" style="display: block; margin: auto;" />

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

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-23-1.png" alt="" width="100%" style="display: block; margin: auto;" />

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

<img src="Fig_5_spatial_analysis_ILC2s_lung_files/figure-html/unnamed-chunk-23-2.png" alt="" width="100%" style="display: block; margin: auto;" />

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
##  [1] patchwork_1.3.2       ggplotify_0.1.3       tidyr_1.3.1           circlize_0.4.17       ComplexHeatmap_2.26.1 stringr_1.6.0         ggbeeswarm_0.7.3      readr_2.1.6           ggpubr_0.6.2          ggplot2_4.0.1         VoltRon_0.2.3         Seurat_5.3.1          Giotto_4.2.2          GiottoClass_0.4.10    rlang_1.1.6           rstatix_0.7.3         dplyr_1.1.4           SeuratObject_5.2.0    sp_2.2-0             
## 
## loaded via a namespace (and not attached):
##   [1] fs_1.6.6                    matrixStats_1.5.0           spatstat.sparse_3.1-0       bitops_1.0-9                lubridate_1.9.4             EBImage_4.52.0              doParallel_1.0.17           httr_1.4.8                  RColorBrewer_1.1-3          tools_4.5.2                 sctransform_0.4.2           backports_1.5.0             utf8_1.2.6                  R6_2.6.1                    lazyeval_0.2.2              uwot_0.2.4                  GetoptLong_1.1.0            withr_3.0.2                 gridExtra_2.3               GiottoUtils_0.2.5           progressr_0.18.0            cli_3.6.5                   Biobase_2.70.0              Cairo_1.7-0                 spatstat.explore_3.5-3      fastDummies_1.7.5           shinyjs_2.1.1               labeling_0.4.3              sass_0.4.10                 S7_0.2.0                    spatstat.data_3.1-9         ggridges_0.5.7              pbapply_1.7-4               yulab.utils_0.2.4           dichromat_2.0-0.1           parallelly_1.45.1           rstudioapi_0.18.0           gridGraphics_0.5-1          shape_1.4.6.1               generics_0.1.4              vroom_1.7.0                 gtools_3.9.5               
##  [43] ica_1.0-3                   spatstat.random_3.4-2       car_3.1-5                   Matrix_1.7-4                S4Vectors_0.48.0            abind_1.4-8                 terra_1.8-93                lifecycle_1.0.5             scatterplot3d_0.3-45        yaml_2.3.10                 snakecase_0.11.1            carData_3.0-6               SummarizedExperiment_1.40.0 gplots_3.3.0                SparseArray_1.10.1          Rtsne_0.17                  promises_1.5.0              crayon_1.5.3                miniUI_0.1.2                lattice_0.22-7              beachmat_2.26.0             cowplot_1.2.0               magick_2.9.0                pillar_1.11.1               knitr_1.51                  GenomicRanges_1.62.0        rjson_0.2.23                future.apply_1.20.2         codetools_0.2-20            glue_1.8.0                  spatstat.univar_3.1-4       data.table_1.17.8           vctrs_0.6.5                 png_0.1-8                   ids_1.0.1                   spam_2.11-1                 gtable_0.3.6                cachem_1.1.0                xfun_0.56                   S4Arrays_1.10.0             mime_0.13                   tidygraph_1.3.1            
##  [85] Seqinfo_1.0.0               survival_3.8-3              SingleCellExperiment_1.32.0 iterators_1.0.14            bluster_1.20.0              rgl_1.3.34                  fitdistrplus_1.2-6          ROCR_1.0-12                 colorsGen_1.0.0             nlme_3.1-168                bit64_4.6.0-1               RcppAnnoy_0.0.22            rprojroot_2.1.1             bslib_0.10.0                irlba_2.3.5.1               vipor_0.4.7                 KernSmooth_2.23-26          otel_0.2.0                  colorspace_2.1-2            BiocGenerics_0.56.0         tidyselect_1.2.1            bit_4.6.0                   compiler_4.5.2              BiocNeighbors_2.4.0         DelayedArray_0.36.0         plotly_4.12.0               checkmate_2.3.4             scales_1.4.0                caTools_1.18.3              lmtest_0.9-40               rappdirs_0.3.4              tiff_0.1-12                 SpatialExperiment_1.20.0    digest_0.6.38               goftest_1.2-3               fftwtools_0.9-11            spatstat.utils_3.2-0        rmarkdown_2.30              XVector_0.50.0              htmltools_0.5.8.1           GiottoVisuals_0.2.14        pkgconfig_2.0.3            
## [127] jpeg_0.1-11                 base64enc_0.1-6             MatrixGenerics_1.22.0       fastmap_1.2.0               GlobalOptions_0.1.3         htmlwidgets_1.6.4           shiny_1.13.0                Rvcg_0.25                   farver_2.1.2                jquerylib_0.1.4             zoo_1.8-14                  jsonlite_2.0.0              BiocParallel_1.44.0         BiocSingular_1.26.0         RCurl_1.98-1.17             magrittr_2.0.4              Formula_1.2-5               dotCall64_1.2               RCDT_1.3.0                  Rcpp_1.1.0                  viridis_0.6.5               reticulate_1.44.0           stringi_1.8.7               ggraph_2.2.2                MASS_7.3-65                 plyr_1.8.9                  parallel_4.5.2              listenv_0.10.0              ggrepel_0.9.6               deldir_2.0-4                graphlayouts_1.2.2          splines_4.5.2               tensor_1.5.1                hms_1.1.4                   locfit_1.5-9.12             colorRamp2_0.1.0            igraph_2.2.1                uuid_1.2-2                  spatstat.geom_3.6-0         ggsignif_0.6.4              RcppHNSW_0.6.0              reshape2_1.4.5             
## [169] stats4_4.5.2                ScaledMatrix_1.18.0         evaluate_1.0.5              foreach_1.5.2               tzdb_0.5.0                  tweenr_2.0.3                httpuv_1.6.16               RANN_2.6.2                  purrr_1.2.0                 polyclip_1.10-7             clue_0.3-67                 future_1.69.0               scattermore_1.2             ggforce_0.5.0               janitor_2.2.1               rsvd_1.0.5                  broom_1.0.12                xtable_1.8-8                RSpectra_0.16-2             later_1.4.4                 viridisLite_0.4.2           Polychrome_1.5.4            tibble_3.3.0                beeswarm_0.4.0              memoise_2.0.1               IRanges_2.44.0              cluster_2.1.8.1             timechange_0.3.0            globals_0.19.0              here_1.0.2
```
