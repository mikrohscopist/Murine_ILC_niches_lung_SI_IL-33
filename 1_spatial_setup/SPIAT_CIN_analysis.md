---
title: "SPIAT Cells in neighborhood"
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
library(SPIAT)
library(dplyr)
library(ggplot2)
library(stringr)
library(glue)
library(here)
library(readr)
library(lubridate)
library(data.table)
library(clustree)
library(magrittr)
library(ggpubr)
library(ggrepel)
library(moments)
library(plotrix)
library(rstatix)
library(forcats)
library(scales)
library(patchwork)
```

## Directories & parameters


``` r
output_env_dir <- here::here("1_spatial_setup", "SPIAT_CIN_analysis_files")
dir.create(output_env_dir)

output_dir <- here::here("data", "CIN_analysis")
dir.create(output_dir)

save_name <- "SPIAT"
```

# Lung

## Load data


``` r
set.seed(8)

# files names of the seurat object of the different annotation levels
filename <- here::here("data", "20230712_SO_33M_arcsinh_lung_0.04_imputed_SPIAT_data.csv")
  
df_all <- read_csv(filename, 
    col_types = cols(...1 = col_skip()))


# df_all <- df_all %>%
#   mutate(Loc_X = Location_Center_X, 
#          Loc_Y = Location_Center_Y, 
#          CellType = annotation_lvl3,
#          AL1 = annotation_lvl1, 
#          AL2 = annotation_lvl2, 
#          AL3 = annotation_lvl3
#          ) %>%
#   select(-c(annotation_lvl1, annotation_lvl2, annotation_lvl3, annotation_lvl4, 
#             Location_Center_X, Location_Center_Y, 
#             UMAP_1, UMAP_2, harmony_1, harmony_2))
# 

unique(df_all$Dataset)
```

```
##  [1] "CTRL_FOV1_20210910" "CTRL_FOV1_20210914" "CTRL_FOV1_20210922" "CTRL_FOV2_20210910" "CTRL_FOV2_20210914" "CTRL_FOV2_20210922" "CTRL_FOV3_20210910" "CTRL_FOV3_20210914" "CTRL_FOV3_20210922" "D1_FOV1_20220311"   "D1_FOV1_20220316"   "D1_FOV1_20220321"   "D1_FOV2_20220311"   "D1_FOV2_20220316"   "D1_FOV2_20220321"   "D1_FOV3_20220311"   "D1_FOV3_20220316"   "D1_FOV3_20220321"   "D2_FOV1_20220325"   "D2_FOV1_20220421"   "D2_FOV1_20220502"   "D2_FOV2_20220325"   "D2_FOV2_20220421"   "D2_FOV2_20220502"   "D2_FOV3_20220325"   "D2_FOV3_20220421"   "D2_FOV3_20220502"   "D3_FOV1_20210902"   "D3_FOV1_20210906"   "D3_FOV1_20210928"   "D3_FOV2_20210902"   "D3_FOV2_20210906"   "D3_FOV2_20210928"   "D3_FOV3_20210902"   "D3_FOV3_20210906"   "D3_FOV3_20210928"
```

## Data preparation

X and y coordinates are converted from pixel to µm by the factor 0.325.


``` r
# max values x and y before conversion
max(df_all$Loc_X)
```

```
## [1] 2034
```

``` r
max(df_all$Loc_Y)
```

```
## [1] 2034
```

``` r
df_all <- df_all %>%
  mutate(Loc_X = Loc_X * 0.325,
         Loc_Y = Loc_Y * 0.325)

# max values x and y after conversion
max(df_all$Loc_X)
```

```
## [1] 661.05
```

``` r
max(df_all$Loc_Y)
```

```
## [1] 661.05
```

Combine ILCs A and ILC2s B into ILC2s


``` r
df_all$ref_markers <- gsub("CD127,CD45,CD90,KLRG1,Kappa,MHCII,GATA3eGFP", 
                           "CD127,CD45,CD90,KLRG1,GATA3eGFP", df_all$ref_markers)
df_all$ref_markers <- gsub("CD127,CD45,CD90,KLRG1,GATA3eGFP,Ki67", 
                           "CD127,CD45,CD90,KLRG1,GATA3eGFP", df_all$ref_markers)

df_all$CellType <- gsub("NK cells/ILC1s", "NK cells & ILC1s", df_all$CellType)
```

Define markers, cell_types and reference:


``` r
# markers defined here have to be in the same order as in the df 
# you can also just simply use colnames of the respective columns
markers <-  c("Areg",
              "B220",
              "CCR6", 
              "CD117",
              "CD11c",
              "CD127",
              "CD138",
              "CD3",
              "CD31",
              "CD4",
              "CD44",
              "CD45",
              "CD68",
              "CD8a",
              "CD90",
              "EMCN",
              "EpCAM",
              "ICOS",
              "KLRG1",
              "Kappa",
              "LYVE1",
              "MHCII", 
              "NKp46",
              "PDGFRa",
              "PDPN",
              "Sca1",
              "EOMES",
              "GATA3",
              "GATA3eGFP",
              "IRF4",
              "Ki67", 
              "RORgt",
              "TBET")

# define cell types (also have to be in same order as in df)
cell_types <- unique(df_all$CellType)

# Specify the markers used for identification of the different cell types here
# format needs to be as shown here
# order needs to be the same as in the df
reference <- c(
      "EpCAM", # Epithelia
      "CD31,EMCN", # EMCN CD31 Blood vessels
      "CD31,LYVE1", # LYVE1 CD31 vessels
      "CD90,LYVE1", # LYVE1 CD90 Lymphatics
      "CD11c,CD45,CD68", # Myeloid cells
      "B220,CD138,CD45,IRF4,MHCII", # B cells & Plasma cells
      "CD117,CD45,NKp46,EOMES,TBET", # NK cells/ILC1s
      "CCR6,CD127,CD45,CD90,ICOS,RORgt", # ILC3s
      "CD3,CD45,CD8a", # T cytotox cells
      "CD3,CD4,CD45", # T helper cells
      "CD127,CD45,CD90,KLRG1,GATA3eGFP" # ILC2s
    )

my_colors <- c(
      "green", # Epithelia
      "red", # EMCN CD31 Blood vessels
      "yellow", # LYVE1 CD31 vessels
      "gold", # LYVE1 CD90 Lymphatics
      "darkviolet", # Myeloid cells
      "deeppink4", # B cells & Plasma cells
      "navy", # NK cells/ILC1s
      "darkcyan", # ILC3s
      "slateblue4", # T cytotox cells
      "deeppink", # T helper cells
      "seagreen2" # ILC2s 
    )
```

## Cells in neighborhood analysis

The next chunk will take the radius and unit provided in the header of the document and loop through all the available datasets in the Dataset column of the df_all dataframe. It will create csv files containing all the frequencies of neighboring cells within the radius provided and save the files in the folder provided.


``` r
for (radius in c(10, 15, 20, 25)) {
  unit <- "micm"
  datasets <- unique(df_all$Dataset)
  organ <- "Lung"
  tissuearea <- "Lung"
  cur_output_dir <- paste0(output_dir, "/", tissuearea, "/", radius, "_", unit)
  dir.create(cur_output_dir)
  
  for (variable in datasets) {
    # extract meta data from dataset name
    exp <- str_sub(variable,-8,-1)   
    fov <- str_sub(variable,-10,-10)   
    treat <- gsub("D", "", str_extract(variable, "[^_]+"))
    dataset <- paste0(exp, "_", fov)
    
  
    # filter dataframe for dataset
    df <- df_all %>%
      filter(Dataset == variable)
    
    
    # create gfi object from filtered dataframe
    ## rows are markers, columns are cells
    df_features <- t(df[, 12:44])
    intensity_matrix <- matrix(df_features, nrow = nrow(df_features), 
                               ncol = ncol(df_features))
    # define marker names as rownames
    rownames(intensity_matrix) <- rownames(df_features)
    
    # define cell IDs as colnames
    colnames(intensity_matrix) <- paste0("Cell_" , df$CellID) 
    
    # metadata (phenotypes, x/y coordinates, cell_types)
    # the order of the elements in these vectors correspond to the cell order 
    # in `intensity matrix`
    #phenotypes <- df$annotation_lvl2
    phenotypes <- df$ref_markers
    coord_x <- df$Loc_X
    coord_y <- df$Loc_Y
    sample_id <- df$Dataset
    condition <- df$Treatment 
    
    # gfi = general format image
    gfi <- format_image_to_spe(format = "general", 
                                                intensity_matrix = intensity_matrix,
                                                phenotypes = phenotypes, 
                                                coord_x = coord_x,
                                                coord_y = coord_y)
    
    # specify cell types
    gfi_form <- define_celltypes(
      gfi, 
      categories = reference, 
      category_colname = "Phenotype", 
      names = cell_types,
      new_colname = "Cell.Type")
    
    
    # Cells in neighborhood
    # Preparing a data frmae to collect all cells in neighborhood 
    #for the respective dataset
    # Reference cell type is assigned in Reference cell column
    # create data frame with 0 rows and 5 columns
    df_neigh <- data.frame(matrix(ncol = length(cell_types)+2, 
                                    nrow = length(cell_types)))
      
    #provide column names
    colnames(df_neigh) <- c("Dataset", cell_types, "Reference cell")
    df_neigh$Dataset <- as.character(variable)
    df_neigh$Organ <- organ
    df_neigh$Tissue.area <- tissuearea
  
    # loop through the cell types and create one data frame with 
    # all the neighboring frequencies for the resepctive reference cell type
    for (index in seq(length(cell_types))) {
        ref <- cell_types[index]
        ref_name <- gsub(" ", "", ref)
        df_neigh[index, "Reference cell"] <- ref
        for (element in cell_types) {
          apocwr <- average_percentage_of_cells_within_radius(spe_object = gfi_form, 
                                                  reference_celltype = ref, 
                                                  target_celltype = element, 
                                                  radius=radius, 
                                                  feature_colname="Cell.Type")
          df_neigh[index, element] <- as.numeric(round(apocwr))
          # print(paste0(ref, 
          #             " in radius of ", radius, " ", 
          #             unit, " around ", element, ": ",
          #             round(apocwr, 1)))
        }
        
        rowSums(df_neigh[2:(length(cell_types)+1)], na.rm = TRUE)
        head(df_neigh)
        }
    
    # save csv file of ILC2 A neighbor hood for respective dataset
    write.csv(df_neigh, 
              file = paste0(cur_output_dir, 
                                    "/", 
                                    save_name, 
                                    "_", 
                                    dataset, 
                                    "_rad", 
                                    radius,"_", unit, 
                            "_", 
                            tissuearea,
                                    "_freq.csv"))
    
    
    
  }

}
```

```
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
```

# Villi

## Load data


``` r
filename <- here::here("data", "SO_arcsinh_si_imputed_all_cells.csv")
  
df_all <- read_csv(filename, 
    col_types = cols(...1 = col_skip()))
head(df_all)
```

```
## # A tibble: 6 × 39
##   Dataset            Treatment CellID Experiment   FOV `Tissue area` Loc_X Loc_Y AL1       AL2       CellType    UMAP_1 UMAP_2  Areg  CCR6 CD117 CD11c CD127   CD3   CD4  CD44  CD45  CD8a  CD90  EMCN EpCAM KLRG1 Kappa LYVE1 MHCII NKp46 PDGFRa  PDPN  Sca1 EOMES GATA3eGFP  IRF4  Ki67 RORgt
##   <chr>              <chr>      <dbl>      <dbl> <dbl> <chr>         <dbl> <dbl> <chr>     <chr>     <chr>        <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl>     <dbl> <dbl> <dbl> <dbl>
## 1 CTRL_FOV1_20210706 CTRL           3   20210706     1 Villi          732.  19.2 Epithelia Epithelia Epithelia I   8.64  -2.33  5.82  0     0     0     0     0     0        0  0     0        0  0     9.08  0     0     0     0     0      0        0     0  0         5.00  0     0     6.55
## 2 CTRL_FOV1_20210706 CTRL           7   20210706     1 Villi          704.  26.3 Epithelia Epithelia Epithelia I  10.7   -1.82  6.22  0     0     0     0     0     0        0  0     0        0  4.34  8.99  4.95  0     0     0     0      0        0     0  0         0     0     0     7.12
## 3 CTRL_FOV1_20210706 CTRL           9   20210706     1 Villi          405.  30.4 Epithelia Epithelia Epithelia I   5.83  13.9   7.03  5.40  7.17  6.16  7.44  7.32  5.86     0  0     0        0  2.01  0     4.27  0     6.01  5.83  7.56   6.27     0     0  0         0     0     0     3.88
## 4 CTRL_FOV1_20210706 CTRL          17   20210706     1 Villi          354.  40.0 Epithelia Epithelia Epithelia I   5.72  13.2   6.84  4.24  7.23  6.61  7.00  7.33  5.57     0  5.40  0        0  4.31  0     5.31  0     6.05  6.10  7.24   6.87     0     0  4.57      0     0     0     5.21
## 5 CTRL_FOV1_20210706 CTRL          19   20210706     1 Villi         1711.  46.5 Epithelia Epithelia Epithelia I   5.54  -2.69  0     0     0     0     0     0     0        0  0     0        0  4.11  9.19  0     3.57  0     0     0      0        0     0  0         4.92  0     0     0   
## 6 CTRL_FOV1_20210706 CTRL          20   20210706     1 Villi          377.  50.0 Epithelia Epithelia Epithelia I   5.35  14.1   6.49  6.49  6.80  6.16  6.89  6.87  5.49     0  5.29  4.30     0  4.02  0     4.24  0     5.93  5.79  6.96   6.52     0     0  6.78      0     4.72  6.32  6.64
```

## Data preparation


``` r
# max values x and y before conversion
max(df_all$Loc_X)
```

```
## [1] 2043.708
```

``` r
max(df_all$Loc_Y)
```

```
## [1] 2042.244
```

``` r
df_all <- df_all %>%
  mutate(Loc_X = Loc_X * 0.325,
         Loc_Y = Loc_Y * 0.325)

# max values x and y after conversion
max(df_all$Loc_X)
```

```
## [1] 664.2052
```

``` r
max(df_all$Loc_Y)
```

```
## [1] 663.7294
```

``` r
unique(df_all$`Tissue area`)
```

```
## [1] "Villi" "ILF"
```


``` r
# markers defined here have to be in the same order as in the df 
# you can also just simply use colnames of the respective columns
markers <-  c("Areg",
              "CCR6", 
              "CD117",
              "CD11c",
              "CD127",
              "CD3",
              "CD4",
              "CD44",
              "CD45",
              "CD8a",
              "CD90",
              "EMCN",
              "EpCAM",
              "KLRG1",
              "Kappa",
              "LYVE1",
              "MHCII", 
              "NKp46",
              "PDGFRa",
              "PDPN",
              "Sca1",
              "EOMES",
              "GATA3eGFP",
              "IRF4",
              "Ki67", 
              "RORgt")

# define cell types (also have to be in same order as in df)
cell_types <- unique(df_all$CellType)

# Specify the markers used for identification of the different cell types here
# format needs to be as shown here
# order needs to be the same as in the df
reference <- c(
      "EpCAM", # Epithelia I
      "CD44,EpCAM,Sca1,Ki67", # Epithelia II
      "CD44,CD90,LYVE1,PDPN,Sca1", # Fibroblasts
      "CD90,EMCN,LYVE1,Sca1", # Blood vessels
      "CD90,LYVE1,PDPN", # Lymphatics
      "CD11c,CD45", # Myeloid cells
      "CD45,Kappa", # B cells 
      "CD45,Kappa,IRF4", #Plasma cells/Plasmablasts
      "CD45,GATA3eGFP", # ILC2s
      "CD45,CD8a,EpCAM", # IEL
      "CCR6,CD117,CD11c,CD127,CD3,CD4,CD90,EMCN,Kappa,LYVE1,MHCII,NKp46,PDGFRa,Sca1,EOMES,IRF4,RORgt", # Unresolved
      "CD127,CD45,EOMES,RORgt", # NK cells7ILC1s/ILC3s
      "CD3,CD45,CD8a", # T cytotox cells
      "CD3,CD4,CD45" # T helper cells
    )

df_all <- df_all %>%
  mutate(`ref_markers` = CellType) %>%
  mutate(`ref_markers` = recode(
    CellType, 
      "Epithelia I" = "EpCAM", # 
      "Epithelia II" = "CD44,EpCAM,Sca1,Ki67", # 
      "Fibroblasts" = "CD44,CD90,LYVE1,PDPN,Sca1", # 
      "Blood vessels" = "CD90,EMCN,LYVE1,Sca1", # 
      "Lymphatics" = "CD90,LYVE1,PDPN", # 
      "Myeloid cells" = "CD11c,CD45", # 
      "B cells" = "CD45,Kappa", #  
      "Plasma cells/Plasmablasts" = "CD45,Kappa,IRF4", #
      "ILC2s" = "CD45,GATA3eGFP", # 
      "CD8+ CD3- IEL" = "CD45,CD8a,EpCAM", # 
      "Unresolved" = "CCR6,CD117,CD11c,CD127,CD3,CD4,CD90,EMCN,Kappa,LYVE1,MHCII,NKp46,PDGFRa,Sca1,EOMES,IRF4,RORgt", # 
      "NK cells/ILC1s/ILC3s" = "CD127,CD45,EOMES,RORgt", # 
      "T cytotox. cells" = "CD3,CD45,CD8a", # 
      "T helper cells" = "CD3,CD4,CD45" # 
    
  ))
```

## Cells in neighborhood


``` r
for (radius in c(10, 15, 20, 25)) {
  unit <- "micm"
  organ <- "SI"
  tissuearea <- "Villi"
  
  # get only the datasets of the selected tissue area
  df_ta <- df_all %>%
    filter(`Tissue area` == tissuearea)
  
  datasets <- unique(df_ta$Dataset)
  cur_output_dir <- paste0(output_dir, "/", tissuearea, "/", radius, "_", unit)
  dir.create(cur_output_dir)
  
  for (variable in datasets) {
    # extract meta data from dataset name
    exp <- str_sub(variable,-8,-1)   
    fov <- str_sub(variable,-10,-10)   
    treat <- gsub("D", "", str_extract(variable, "[^_]+"))
    dataset <- paste0(exp, "_", fov)
    
  
    # filter dataframe for dataset
    df <- df_all %>%
      filter(Dataset == variable)
    
    
    ## rows are markers, columns are cells
    df_features <- t(df[, 14:39])
    intensity_matrix <- matrix(df_features, nrow = nrow(df_features), 
                               ncol = ncol(df_features))
    # define marker names as rownames
    rownames(intensity_matrix) <- rownames(df_features)
    
    # define cell IDs as colnames
    colnames(intensity_matrix) <- paste0("Cell_" , df$CellID) 
    
    # metadata (phenotypes, x/y coordinates, cell_types)
    # the order of the elements in these vectors correspond to the cell order 
    # in `intensity matrix`
    #phenotypes <- df$annotation_lvl2
    phenotypes <- df$ref_markers
    coord_x <- df$Loc_X
    coord_y <- df$Loc_Y
    sample_id <- df$Dataset
    condition <- df$Treatment 
    
    # gfi = general format image
    gfi <- format_image_to_spe(format = "general", 
                                                intensity_matrix = intensity_matrix,
                                                phenotypes = phenotypes, 
                                                coord_x = coord_x,
                                                coord_y = coord_y)  
      
  
    gfi_form <- define_celltypes(
        gfi, 
        categories = reference, 
        category_colname = "Phenotype", 
        names = cell_types,
        new_colname = "Cell.Type")
  
    # Cells in neighborhood
    # Preparing a data frmae to collect all cells in neighborhood 
    #for the respective dataset
    # Reference cell type is assigned in Reference cell column
    # create data frame with 0 rows and 5 columns
    df_neigh <- data.frame(matrix(ncol = length(cell_types)+2, 
                                    nrow = length(cell_types)))
      
    #provide column names
    colnames(df_neigh) <- c("Dataset", cell_types, "Reference cell")
    df_neigh$Dataset <- as.character(variable)
    df_neigh$Organ <- organ
    df_neigh$Tissue.area <- tissuearea
  
    # loop through the cell types and create one data frame with 
    # all the neighboring frequencies for the resepctive reference cell type
    for (index in seq(length(cell_types))) {
        ref <- cell_types[index]
        ref_name <- gsub(" ", "", ref)
        df_neigh[index, "Reference cell"] <- ref
        for (element in cell_types) {
          apocwr <- average_percentage_of_cells_within_radius(spe_object = gfi_form, 
                                                  reference_celltype = ref, 
                                                  target_celltype = element, 
                                                  radius=radius, 
                                                  feature_colname="Cell.Type")
          df_neigh[index, element] <- as.numeric(round(apocwr))
          # print(paste0(ref, 
          #             " in radius of ", radius, " ", unit, " around ", element, ": ",
          #             round(apocwr, 1)))
        }
        
        rowSums(df_neigh[2:(length(cell_types)+1)], na.rm = TRUE)
        head(df_neigh)
        }
    
    # save csv file of ILC2 A neighbor hood for respective dataset
    write.csv(df_neigh, 
              file = paste0(cur_output_dir, 
                                    "/", 
                                    save_name, 
                                    "_", 
                                    dataset, 
                                    "_rad", 
                                    radius,"_", unit, 
                            "_", 
                            tissuearea,
                                    "_freq.csv"))
    
    
    
  }
  

}
```

```
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
## [1] "There are no reference cells or no target cells"
```

# ILF

## Cells in neighborhood


``` r
for (radius in c(10, 15, 20, 25)) {
  unit <- "micm"
  organ <- "SI"
  tissuearea <- "ILF"
  
  # get only the datasets of the selected tissue area
  df_ta <- df_all %>%
    filter(`Tissue area` == tissuearea)
  
  datasets <- unique(df_ta$Dataset)
  cur_output_dir <- paste0(output_dir, "/", tissuearea, "/", radius, "_", unit)
  dir.create(cur_output_dir)
  
  for (variable in datasets) {
    # extract meta data from dataset name
    exp <- str_sub(variable,-8,-1)   
    fov <- str_sub(variable,-10,-10)   
    treat <- gsub("D", "", str_extract(variable, "[^_]+"))
    dataset <- paste0(exp, "_", fov)
    
  
    # filter dataframe for dataset
    df <- df_all %>%
      filter(Dataset == variable)
    
    
    ## rows are markers, columns are cells
    df_features <- t(df[, 14:39])
    intensity_matrix <- matrix(df_features, nrow = nrow(df_features), 
                               ncol = ncol(df_features))
    # define marker names as rownames
    rownames(intensity_matrix) <- rownames(df_features)
    
    # define cell IDs as colnames
    colnames(intensity_matrix) <- paste0("Cell_" , df$CellID) 
    
    # metadata (phenotypes, x/y coordinates, cell_types)
    # the order of the elements in these vectors correspond to the cell order 
    # in `intensity matrix`
    #phenotypes <- df$annotation_lvl2
    phenotypes <- df$ref_markers
    coord_x <- df$Loc_X
    coord_y <- df$Loc_Y
    sample_id <- df$Dataset
    condition <- df$Treatment 
    
    # gfi = general format image
    gfi <- format_image_to_spe(format = "general", 
                                                intensity_matrix = intensity_matrix,
                                                phenotypes = phenotypes, 
                                                coord_x = coord_x,
                                                coord_y = coord_y)  
      
  
    gfi_form <- define_celltypes(
        gfi, 
        categories = reference, 
        category_colname = "Phenotype", 
        names = cell_types,
        new_colname = "Cell.Type")
  
    # Cells in neighborhood
    # Preparing a data frmae to collect all cells in neighborhood 
    #for the respective dataset
    # Reference cell type is assigned in Reference cell column
    # create data frame with 0 rows and 5 columns
    df_neigh <- data.frame(matrix(ncol = length(cell_types)+2, 
                                    nrow = length(cell_types)))
      
    #provide column names
    colnames(df_neigh) <- c("Dataset", cell_types, "Reference cell")
    df_neigh$Dataset <- as.character(variable)
    df_neigh$Organ <- organ
    df_neigh$Tissue.area <- tissuearea
  
    # loop through the cell types and create one data frame with 
    # all the neighboring frequencies for the resepctive reference cell type
    for (index in seq(length(cell_types))) {
        ref <- cell_types[index]
        ref_name <- gsub(" ", "", ref)
        df_neigh[index, "Reference cell"] <- ref
        for (element in cell_types) {
          apocwr <- average_percentage_of_cells_within_radius(spe_object = gfi_form, 
                                                  reference_celltype = ref, 
                                                  target_celltype = element, 
                                                  radius=radius, 
                                                  feature_colname="Cell.Type")
          df_neigh[index, element] <- as.numeric(round(apocwr))
          # print(paste0(ref, 
          #             " in radius of ", radius, " ", unit, " around ", element, ": ",
          #             round(apocwr, 1)))
        }
        
        rowSums(df_neigh[2:(length(cell_types)+1)], na.rm = TRUE)
        head(df_neigh)
        }
    
    # save csv file of ILC2 A neighbor hood for respective dataset
    write.csv(df_neigh, 
              file = paste0(cur_output_dir, 
                                    "/", 
                                    save_name, 
                                    "_", 
                                    dataset, 
                                    "_rad", 
                                    radius,"_", unit, 
                            "_", 
                            tissuearea,
                                    "_freq.csv"))
    
    
    
  }
  

}
```

# Session Information


``` r
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
## [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] patchwork_1.3.1             scales_1.4.0                forcats_1.0.0               rstatix_0.7.2               plotrix_3.8-4               moments_0.14.1              ggrepel_0.9.6               ggpubr_0.6.0                magrittr_2.0.3              clustree_0.5.1              ggraph_2.2.1                data.table_1.17.0           lubridate_1.9.4             readr_2.1.5                 here_1.0.1                  glue_1.8.0                  stringr_1.5.1               ggplot2_3.5.2               dplyr_1.1.4                 SPIAT_1.8.1                 SpatialExperiment_1.16.0    SingleCellExperiment_1.28.1 SummarizedExperiment_1.36.0 Biobase_2.66.0              GenomicRanges_1.58.0        GenomeInfoDb_1.42.3         IRanges_2.40.1              S4Vectors_0.44.0            BiocGenerics_0.52.0         MatrixGenerics_1.18.1       matrixStats_1.5.0          
## 
## loaded via a namespace (and not attached):
##  [1] deldir_2.0-4            gridExtra_2.3           rlang_1.1.5             compiler_4.4.2          spatstat.geom_3.3-6     vctrs_0.6.5             pkgconfig_2.0.3         crayon_1.5.3            fastmap_1.2.0           backports_1.5.0         magick_2.8.7            XVector_0.46.0          utf8_1.2.6              rmarkdown_2.29          tzdb_0.4.0              UCSC.utils_1.2.0        bit_4.6.0               purrr_1.0.4             xfun_0.51               zlibbioc_1.52.0         cachem_1.1.0            jsonlite_1.9.1          goftest_1.2-3           DelayedArray_0.32.0     spatstat.utils_3.1-3    tweenr_2.0.3            parallel_4.4.2          broom_1.0.8             R6_2.6.1                bslib_0.9.0             stringi_1.8.4           RColorBrewer_1.1-3      spatstat.data_3.1-6     car_3.1-3               spatstat.univar_3.1-2   jquerylib_0.1.4         Rcpp_1.0.14             knitr_1.50              tensor_1.5.1            Matrix_1.7-1            igraph_2.1.4            timechange_0.3.0        tidyselect_1.2.1        rstudioapi_0.17.1       abind_1.4-8             yaml_2.3.10             viridis_0.6.5           spatstat.random_3.3-3   spatstat.explore_3.4-2 
## [50] lattice_0.22-6          tibble_3.2.1            withr_3.0.2             evaluate_1.0.4          polyclip_1.10-7         pillar_1.10.2           carData_3.0-5           dbscan_1.2.2            generics_0.1.4          vroom_1.6.5             rprojroot_2.0.4         hms_1.1.3               tools_4.4.2             ggsignif_0.6.4          graphlayouts_1.2.2      tidygraph_1.3.1         grid_4.4.2              tidyr_1.3.1             nlme_3.1-166            GenomeInfoDbData_1.2.13 ggforce_0.5.0           Formula_1.2-5           cli_3.6.3               spatstat.sparse_3.1-0   S4Arrays_1.6.0          viridisLite_0.4.2       gtable_0.3.6            sass_0.4.10             digest_0.6.37           SparseArray_1.6.2       rjson_0.2.23            farver_2.1.2            memoise_2.0.1           htmltools_0.5.8.1       lifecycle_1.0.4         httr_1.4.7              bit64_4.6.0-1           MASS_7.3-61
```
