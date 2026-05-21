---
title: "Figure 5: Quantification lung"
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
library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(rstatix)
library(ggbeeswarm)
library(readr)
library(ggpubr)
```

## Parameters


``` r
set.seed(123)

input_dir <- here::here("1_data_tidying", "Lung_SI_all_cells_all_ALs_files")

output_dir <- here::here("2_visualizations_for_figures", "Fig_5_quantification_lung_files")
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

cols_treat <- c("darkcyan", "gold", "deeppink", "slateblue")
```

# Load data


``` r
df_lung <- read_csv(paste0(input_dir, "/lung_proportions.csv"), 
    col_types = cols(...1 = col_skip()))
```

# Visualization

## 4A - Count of ILCs \@CTRL conditions


``` r
# filter for CTRL and convert to longer format
df <- df_lung %>%
  filter(Treatment == "CTRL") %>%
  select(Treatment, Dataset, `NK cells/ILC1s`, ILC2s, ILC3s) %>%
  tidyr::pivot_longer(cols = c(`NK cells/ILC1s`, ILC2s, ILC3s), names_to = "ILCtype") %>%
  mutate(ILCtype = factor(ILCtype, level = c(
    "NK cells/ILC1s", "ILC2s", "ILC3s"
  )))

# Testing for normal distribution
shapiro.test(df$value)
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  df$value
## W = 0.88672, p-value = 0.006732
```

``` r
moments::kurtosis(df$value)
```

```
## [1] 1.622249
```

``` r
moments::jarque.test(df$value)
```

```
## 
## 	Jarque-Bera Normality Test
## 
## data:  df$value
## JB = 2.4178, p-value = 0.2985
## alternative hypothesis: greater
```

``` r
# Kruskal-Wallis-test to check for significance between tested groups and effect size
res.kruskal <- df %>% kruskal_test(`value` ~ ILCtype)
res.kruskal
```

```
## # A tibble: 1 × 6
##   .y.       n statistic    df        p method        
## * <chr> <int>     <dbl> <int>    <dbl> <chr>         
## 1 value    27      18.1     2 0.000116 Kruskal-Wallis
```

``` r
df %>% kruskal_effsize(`value` ~ ILCtype)
```

```
## # A tibble: 1 × 5
##   .y.       n effsize method  magnitude
## * <chr> <int>   <dbl> <chr>   <ord>    
## 1 value    27   0.672 eta2[H] large
```

``` r
# Pairwise comparisons using Dunn's test
pwc <- df %>% 
  dunn_test(`value` ~ ILCtype, p.adjust.method = "bonferroni") 
pwc
```

```
## # A tibble: 3 × 9
##   .y.   group1         group2    n1    n2 statistic         p    p.adj p.adj.signif
## * <chr> <chr>          <chr>  <int> <int>     <dbl>     <dbl>    <dbl> <chr>       
## 1 value NK cells/ILC1s ILC2s      9     9     0.678 0.498     1        ns          
## 2 value NK cells/ILC1s ILC3s      9     9    -3.30  0.000962  0.00289  **          
## 3 value ILC2s          ILC3s      9     9    -3.98  0.0000690 0.000207 ***
```

``` r
# add N to plot
tab <- data.frame(xtabs(~ ILCtype, data = df))
head(tab)
```

```
##          ILCtype Freq
## 1 NK cells/ILC1s    9
## 2          ILC2s    9
## 3          ILC3s    9
```

``` r
# Add cell number per cluster to cluster labels
Labels = paste0("n = ", tab$Freq, "")



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "ILCtype")

plot_prop <- ggplot(df, aes(x = ILCtype, y = value, fill = "ILCtype"))+
  geom_boxplot(fill="white")+
  geom_beeswarm(aes(color = ILCtype), size = 2, cex = 3)+
  scale_color_manual(values = cols_ilcs_lung)+
  theme_classic() +
  theme(
      axis.text.x = element_text(angle = 45, size = 10, face = "bold", hjust = 1),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
      plot.margin = margin(0.5, 1, 0.5, 1, "cm"),
      legend.position = "none"
    ) +
  stat_pvalue_manual(pwc,
                       hide.ns = TRUE, size = 6,
                       step.increase = 0.1, y.position = 20) +
  xlab(NULL)+
  ylab("Count per FOV [#]")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,24))+
    # annotate(geom = 'text',
    #        x="ILC2s",
    #        y=28,
    #        label=Labels[1], 
    #        #angle = 90, 
    #        size = 10/.pt)+
  NoLegend()

# plot_prop
```


``` r
# filter for CTRL and convert to longer format
df <- df_lung %>%
  filter(Treatment == "3") %>%
  select(Treatment, Dataset, `NK cells/ILC1s`, ILC2s, ILC3s) %>%
  tidyr::pivot_longer(cols = c(`NK cells/ILC1s`, ILC2s, ILC3s), names_to = "ILCtype") %>%
  mutate(ILCtype = factor(ILCtype, level = c(
    "NK cells/ILC1s", "ILC2s", "ILC3s"
  )))

# Testing for normal distribution
shapiro.test(df$value)
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  df$value
## W = 0.81232, p-value = 0.0002255
```

``` r
moments::kurtosis(df$value)
```

```
## [1] 3.499666
```

``` r
moments::jarque.test(df$value)
```

```
## 
## 	Jarque-Bera Normality Test
## 
## data:  df$value
## JB = 7.4486, p-value = 0.02413
## alternative hypothesis: greater
```

``` r
# Kruskal-Wallis-test to check for significance between tested groups and effect size
res.kruskal <- df %>% kruskal_test(`value` ~ ILCtype)
res.kruskal
```

```
## # A tibble: 1 × 6
##   .y.       n statistic    df         p method        
## * <chr> <int>     <dbl> <int>     <dbl> <chr>         
## 1 value    27      20.9     2 0.0000291 Kruskal-Wallis
```

``` r
df %>% kruskal_effsize(`value` ~ ILCtype)
```

```
## # A tibble: 1 × 5
##   .y.       n effsize method  magnitude
## * <chr> <int>   <dbl> <chr>   <ord>    
## 1 value    27   0.787 eta2[H] large
```

``` r
# Pairwise comparisons using Dunn's test
pwc <- df %>% 
  dunn_test(`value` ~ ILCtype, p.adjust.method = "bonferroni") 
pwc
```

```
## # A tibble: 3 × 9
##   .y.   group1         group2    n1    n2 statistic          p     p.adj p.adj.signif
## * <chr> <chr>          <chr>  <int> <int>     <dbl>      <dbl>     <dbl> <chr>       
## 1 value NK cells/ILC1s ILC2s      9     9      1.87 0.0612     0.184     ns          
## 2 value NK cells/ILC1s ILC3s      9     9     -2.67 0.00748    0.0224    *           
## 3 value ILC2s          ILC3s      9     9     -4.55 0.00000544 0.0000163 ****
```

``` r
# add N to plot
tab <- data.frame(xtabs(~ ILCtype, data = df))
head(tab)
```

```
##          ILCtype Freq
## 1 NK cells/ILC1s    9
## 2          ILC2s    9
## 3          ILC3s    9
```

``` r
# Add cell number per cluster to cluster labels
Labels = paste0("n = ", tab$Freq, "")



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "ILCtype")

plot_d3 <- ggplot(df, aes(x = ILCtype, y = value, fill = "ILCtype"))+
  geom_boxplot(fill="white")+
  geom_beeswarm(aes(color = ILCtype), size = 2, cex = 3)+
  scale_color_manual(values = cols_ilcs_lung)+
  theme_classic() +
  theme(
      axis.text.x = element_text(angle = 45, size = 10, face = "bold", hjust = 1),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
      plot.margin = margin(0.5, 1, 0.5, 1, "cm"),
      legend.position = "none"
    ) +
  stat_pvalue_manual(pwc,
                       hide.ns = TRUE, size = 6,
                       step.increase = 0.1, y.position = 180) +
  xlab(NULL)+
  ylab("Count per FOV [#]")+
  ggtitle("ILC subtypes D3")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,250))+
    # annotate(geom = 'text',
    #        x="ILC2s",
    #        y=28,
    #        label=Labels[1], 
    #        #angle = 90, 
    #        size = 10/.pt)+
  NoLegend()

# plot_d3
```

## 4C - Frequencies of ILCs within ILC compartment \@ CTRL conditions


``` r
# filter for CTRL and convert to longer format
df <- df_lung %>%
  filter(Treatment == "CTRL") %>%
  select(Treatment, Dataset, `Prop_NK cells/ILC1s_perTotalILCsFOV`, Prop_ILC2s_perTotalILCsFOV, Prop_ILC3s_perTotalILCsFOV) %>%
  tidyr::pivot_longer(cols = c(`Prop_NK cells/ILC1s_perTotalILCsFOV`, Prop_ILC2s_perTotalILCsFOV, Prop_ILC3s_perTotalILCsFOV), names_to = "ILCtype") %>%
  mutate(ILCtype = gsub("Prop_|_perTotalILCsFOV", "", ILCtype), 
         ILCtype = factor(ILCtype, level = c(
    "NK cells/ILC1s", "ILC2s", "ILC3s"
  )))

# Testing for normal distribution
shapiro.test(df$value)
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  df$value
## W = 0.89454, p-value = 0.01006
```

``` r
moments::kurtosis(df$value)
```

```
## [1] 1.622854
```

``` r
moments::jarque.test(df$value)
```

```
## 
## 	Jarque-Bera Normality Test
## 
## data:  df$value
## JB = 2.2919, p-value = 0.3179
## alternative hypothesis: greater
```

``` r
# Kruskal-Wallis-test to check for significance between tested groups and effect size
res.kruskal <- df %>% kruskal_test(`value` ~ ILCtype)
res.kruskal
```

```
## # A tibble: 1 × 6
##   .y.       n statistic    df         p method        
## * <chr> <int>     <dbl> <int>     <dbl> <chr>         
## 1 value    27      18.5     2 0.0000974 Kruskal-Wallis
```

``` r
df %>% kruskal_effsize(`value` ~ ILCtype)
```

```
## # A tibble: 1 × 5
##   .y.       n effsize method  magnitude
## * <chr> <int>   <dbl> <chr>   <ord>    
## 1 value    27   0.686 eta2[H] large
```

``` r
# Pairwise comparisons using Dunn's test
pwc <- df %>% 
  dunn_test(`value` ~ ILCtype, p.adjust.method = "bonferroni") 
pwc
```

```
## # A tibble: 3 × 9
##   .y.   group1         group2    n1    n2 statistic         p    p.adj p.adj.signif
## * <chr> <chr>          <chr>  <int> <int>     <dbl>     <dbl>    <dbl> <chr>       
## 1 value NK cells/ILC1s ILC2s      9     9     0.812 0.417     1        ns          
## 2 value NK cells/ILC1s ILC3s      9     9    -3.25  0.00116   0.00347  **          
## 3 value ILC2s          ILC3s      9     9    -4.06  0.0000488 0.000146 ***
```

``` r
# add N to plot
tab <- data.frame(xtabs(~ ILCtype, data = df))
head(tab)
```

```
##          ILCtype Freq
## 1 NK cells/ILC1s    9
## 2          ILC2s    9
## 3          ILC3s    9
```

``` r
# Add cell number per cluster to cluster labels
Labels = paste0("n = ", tab$Freq, "")



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "ILCtype")

plot_freq <- ggplot(df, aes(x = ILCtype, y = value, fill = "ILCtype"))+
  geom_boxplot(fill="white")+
  geom_beeswarm(aes(color = ILCtype), size = 2, cex = 3)+
  scale_color_manual(values = cols_ilcs_lung)+
  theme_classic() +
  theme(
      axis.text.x = element_text(angle = 45, size = 10, face = "bold", hjust = 1),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
      plot.margin = margin(0.5, 1, 0.5, 1, "cm"),
      legend.position = "none"
    ) +
  stat_pvalue_manual(pwc,
                       hide.ns = TRUE, size = 6,
                       step.increase = 0.1, y.position = 90) +
  xlab(NULL)+
  ylab("Frequency per FOV/ILCs [%]")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,105))+
  NoLegend()


# plot_freq
```


``` r
# filter for CTRL and convert to longer format
df <- df_lung %>%
  filter(Treatment == "3") %>%
  select(Treatment, Dataset, `Prop_NK cells/ILC1s_perTotalILCsFOV`, Prop_ILC2s_perTotalILCsFOV, Prop_ILC3s_perTotalILCsFOV) %>%
  tidyr::pivot_longer(cols = c(`Prop_NK cells/ILC1s_perTotalILCsFOV`, Prop_ILC2s_perTotalILCsFOV, Prop_ILC3s_perTotalILCsFOV), names_to = "ILCtype") %>%
  mutate(ILCtype = gsub("Prop_|_perTotalILCsFOV", "", ILCtype), 
         ILCtype = factor(ILCtype, level = c(
    "NK cells/ILC1s", "ILC2s", "ILC3s"
  )))

# Testing for normal distribution
shapiro.test(df$value)
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  df$value
## W = 0.87901, p-value = 0.004573
```

``` r
moments::kurtosis(df$value)
```

```
## [1] 1.740977
```

``` r
moments::jarque.test(df$value)
```

```
## 
## 	Jarque-Bera Normality Test
## 
## data:  df$value
## JB = 2.4851, p-value = 0.2886
## alternative hypothesis: greater
```

``` r
# Kruskal-Wallis-test to check for significance between tested groups and effect size
res.kruskal <- df %>% kruskal_test(`value` ~ ILCtype)
res.kruskal
```

```
## # A tibble: 1 × 6
##   .y.       n statistic    df          p method        
## * <chr> <int>     <dbl> <int>      <dbl> <chr>         
## 1 value    27      23.1     2 0.00000943 Kruskal-Wallis
```

``` r
df %>% kruskal_effsize(`value` ~ ILCtype)
```

```
## # A tibble: 1 × 5
##   .y.       n effsize method  magnitude
## * <chr> <int>   <dbl> <chr>   <ord>    
## 1 value    27   0.881 eta2[H] large
```

``` r
# Pairwise comparisons using Dunn's test
pwc <- df %>% 
  dunn_test(`value` ~ ILCtype, p.adjust.method = "bonferroni") 
pwc
```

```
## # A tibble: 3 × 9
##   .y.   group1         group2    n1    n2 statistic          p      p.adj p.adj.signif
## * <chr> <chr>          <chr>  <int> <int>     <dbl>      <dbl>      <dbl> <chr>       
## 1 value NK cells/ILC1s ILC2s      9     9      2.41 0.0162     0.0485     *           
## 2 value NK cells/ILC1s ILC3s      9     9     -2.41 0.0162     0.0485     *           
## 3 value ILC2s          ILC3s      9     9     -4.81 0.00000150 0.00000451 ****
```

``` r
# add N to plot
tab <- data.frame(xtabs(~ ILCtype, data = df))
head(tab)
```

```
##          ILCtype Freq
## 1 NK cells/ILC1s    9
## 2          ILC2s    9
## 3          ILC3s    9
```

``` r
# Add cell number per cluster to cluster labels
Labels = paste0("n = ", tab$Freq, "")



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "ILCtype")

plot_freq_d3 <- ggplot(df, aes(x = ILCtype, y = value, fill = "ILCtype"))+
  geom_boxplot(fill="white")+
  geom_beeswarm(aes(color = ILCtype), size = 2, cex = 3)+
  scale_color_manual(values = cols_ilcs_lung)+
  theme_classic() +
  theme(
      axis.text.x = element_text(angle = 45, size = 10, face = "bold", hjust = 1),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
      plot.subtitle = element_text(hjust = 0.5, size = 11), 
      plot.margin = margin(0.5, 1, 0.5, 1, "cm"),
      legend.position = "none"
    ) +
  stat_pvalue_manual(pwc,
                       hide.ns = TRUE, size = 6,
                       step.increase = 0.1, y.position = 90) +
  xlab(NULL)+
  ylab("Frequency per FOV/ILCs [%]")+
  ggtitle("D3", subtitle = "ILC compartment")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,119))+
  NoLegend()


# plot_freq_d3
```

## 4B - Frequency of ILCs within immune compartment \@ CTRL


``` r
# filter for CTRL and convert to longer format
df <- df_lung %>%
  filter(Treatment == "CTRL") %>%
  select(Treatment, Dataset, `Prop_NK cells/ILC1s_perTotalImmuneFOV`, Prop_ILC2s_perTotalImmuneFOV, Prop_ILC3s_perTotalImmuneFOV) %>%
  tidyr::pivot_longer(cols = c(`Prop_NK cells/ILC1s_perTotalImmuneFOV`, Prop_ILC2s_perTotalImmuneFOV, Prop_ILC3s_perTotalImmuneFOV), names_to = "ILCtype") %>%
  mutate(ILCtype = gsub("Prop_|_perTotalImmuneFOV", "", ILCtype), 
         ILCtype = factor(ILCtype, level = c(
    "NK cells/ILC1s", "ILC2s", "ILC3s"
  )))

# Testing for normal distribution
shapiro.test(df$value)
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  df$value
## W = 0.90165, p-value = 0.0146
```

``` r
moments::kurtosis(df$value)
```

```
## [1] 2.567516
```

``` r
moments::jarque.test(df$value)
```

```
## 
## 	Jarque-Bera Normality Test
## 
## data:  df$value
## JB = 2.3173, p-value = 0.3139
## alternative hypothesis: greater
```

``` r
# Kruskal-Wallis-test to check for significance between tested groups and effect size
res.kruskal <- df %>% kruskal_test(`value` ~ ILCtype)
res.kruskal
```

```
## # A tibble: 1 × 6
##   .y.       n statistic    df        p method        
## * <chr> <int>     <dbl> <int>    <dbl> <chr>         
## 1 value    27      17.9     2 0.000132 Kruskal-Wallis
```

``` r
df %>% kruskal_effsize(`value` ~ ILCtype)
```

```
## # A tibble: 1 × 5
##   .y.       n effsize method  magnitude
## * <chr> <int>   <dbl> <chr>   <ord>    
## 1 value    27   0.661 eta2[H] large
```

``` r
# Pairwise comparisons using Dunn's test
pwc <- df %>% 
  dunn_test(`value` ~ ILCtype, p.adjust.method = "bonferroni") 
pwc
```

```
## # A tibble: 3 × 9
##   .y.   group1         group2    n1    n2 statistic         p    p.adj p.adj.signif
## * <chr> <chr>          <chr>  <int> <int>     <dbl>     <dbl>    <dbl> <chr>       
## 1 value NK cells/ILC1s ILC2s      9     9     0.692 0.489     1        ns          
## 2 value NK cells/ILC1s ILC3s      9     9    -3.26  0.00110   0.00329  **          
## 3 value ILC2s          ILC3s      9     9    -3.96  0.0000760 0.000228 ***
```

``` r
# add N to plot
tab <- data.frame(xtabs(~ ILCtype, data = df))
head(tab)
```

```
##          ILCtype Freq
## 1 NK cells/ILC1s    9
## 2          ILC2s    9
## 3          ILC3s    9
```

``` r
# Add cell number per cluster to cluster labels
Labels = paste0("n = ", tab$Freq, "")



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "ILCtype")

plot_freq_immune <- ggplot(df, aes(x = ILCtype, y = value, fill = "ILCtype"))+
  geom_boxplot(fill="white")+
  geom_beeswarm(aes(color = ILCtype), size = 2, cex = 3)+
  scale_color_manual(values = cols_ilcs_lung)+
  theme_classic() +
  theme(
      axis.text.x = element_text(angle = 45, size = 10, face = "bold", hjust = 1),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
      plot.margin = margin(0.5, 1, 0.5, 1, "cm"),
      legend.position = "none"
    ) +
  stat_pvalue_manual(pwc,
                       hide.ns = TRUE, size = 6,
                       step.increase = 0.1, y.position = 11) +
  xlab(NULL)+
  ylab("Frequency per FOV/Immune [%]")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,14))+
  NoLegend()

# plot_freq_immune
```


``` r
# filter for CTRL and convert to longer format
df <- df_lung %>%
  filter(Treatment == "3") %>%
  select(Treatment, Dataset, `Prop_NK cells/ILC1s_perTotalImmuneFOV`, Prop_ILC2s_perTotalImmuneFOV, Prop_ILC3s_perTotalImmuneFOV) %>%
  tidyr::pivot_longer(cols = c(`Prop_NK cells/ILC1s_perTotalImmuneFOV`, Prop_ILC2s_perTotalImmuneFOV, Prop_ILC3s_perTotalImmuneFOV), names_to = "ILCtype") %>%
  mutate(ILCtype = gsub("Prop_|_perTotalImmuneFOV", "", ILCtype), 
         ILCtype = factor(ILCtype, level = c(
    "NK cells/ILC1s", "ILC2s", "ILC3s"
  )))

# Testing for normal distribution
shapiro.test(df$value)
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  df$value
## W = 0.85373, p-value = 0.001367
```

``` r
moments::kurtosis(df$value)
```

```
## [1] 3.812105
```

``` r
moments::jarque.test(df$value)
```

```
## 
## 	Jarque-Bera Normality Test
## 
## data:  df$value
## JB = 7.3038, p-value = 0.02594
## alternative hypothesis: greater
```

``` r
# Kruskal-Wallis-test to check for significance between tested groups and effect size
res.kruskal <- df %>% kruskal_test(`value` ~ ILCtype)
res.kruskal
```

```
## # A tibble: 1 × 6
##   .y.       n statistic    df          p method        
## * <chr> <int>     <dbl> <int>      <dbl> <chr>         
## 1 value    27      23.2     2 0.00000917 Kruskal-Wallis
```

``` r
df %>% kruskal_effsize(`value` ~ ILCtype)
```

```
## # A tibble: 1 × 5
##   .y.       n effsize method  magnitude
## * <chr> <int>   <dbl> <chr>   <ord>    
## 1 value    27   0.883 eta2[H] large
```

``` r
# Pairwise comparisons using Dunn's test
pwc <- df %>% 
  dunn_test(`value` ~ ILCtype, p.adjust.method = "bonferroni") 
pwc
```

```
## # A tibble: 3 × 9
##   .y.   group1         group2    n1    n2 statistic          p      p.adj p.adj.signif
## * <chr> <chr>          <chr>  <int> <int>     <dbl>      <dbl>      <dbl> <chr>       
## 1 value NK cells/ILC1s ILC2s      9     9      2.41 0.0160     0.0481     *           
## 2 value NK cells/ILC1s ILC3s      9     9     -2.41 0.0160     0.0481     *           
## 3 value ILC2s          ILC3s      9     9     -4.82 0.00000146 0.00000438 ****
```

``` r
# add N to plot
tab <- data.frame(xtabs(~ ILCtype, data = df))
head(tab)
```

```
##          ILCtype Freq
## 1 NK cells/ILC1s    9
## 2          ILC2s    9
## 3          ILC3s    9
```

``` r
# Add cell number per cluster to cluster labels
Labels = paste0("n = ", tab$Freq, "")



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "ILCtype")

plot_freq_immune_d3 <- ggplot(df, aes(x = ILCtype, y = value, fill = "ILCtype"))+
  geom_boxplot(fill="white")+
  geom_beeswarm(aes(color = ILCtype), size = 2, cex = 3)+
  scale_color_manual(values = cols_ilcs_lung)+
  theme_classic() +
  theme(
      axis.text.x = element_text(angle = 45, size = 10, face = "bold", hjust = 1),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
      plot.subtitle = element_text(hjust = 0.5, size = 11), 
      plot.margin = margin(0.5, 1, 0.5, 1, "cm"),
      legend.position = "none"
    ) +
  stat_pvalue_manual(pwc,
                       hide.ns = TRUE, size = 6,
                       step.increase = 0.1, y.position = 15) +
  xlab(NULL)+
  ylab("Frequency per FOV/Immune [%]")+
  ggtitle("D3", subtitle = "Immune compartment")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,20))+
  NoLegend()

# plot_freq_immune_d3
```

## Total cell count IL-33


``` r
# filter for CTRL and convert to longer format
df <- df_lung %>%
  select(Treatment, Dataset, TotalCellCountFOV) %>%
  mutate(Treatment = paste0("D", Treatment)) %>%
  mutate(Treatment = gsub("DCTRL", "CTRL", Treatment)) %>%
  mutate(Treatment = factor(Treatment, level =c(
    "CTRL", "D1", "D2", "D3"
  )))

# Testing for normal distribution
shapiro.test(df$TotalCellCountFOV)
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  df$TotalCellCountFOV
## W = 0.82636, p-value = 5.688e-05
```

``` r
# Kruskal-Wallis-test to check for significance between tested groups and effect size
res.kruskal <- df %>% kruskal_test(`TotalCellCountFOV` ~ Treatment)
res.kruskal
```

```
## # A tibble: 1 × 6
##   .y.                   n statistic    df          p method        
## * <chr>             <int>     <dbl> <int>      <dbl> <chr>         
## 1 TotalCellCountFOV    36      27.0     3 0.00000597 Kruskal-Wallis
```

``` r
df %>% kruskal_effsize(`TotalCellCountFOV` ~ Treatment)
```

```
## # A tibble: 1 × 5
##   .y.                   n effsize method  magnitude
## * <chr>             <int>   <dbl> <chr>   <ord>    
## 1 TotalCellCountFOV    36   0.749 eta2[H] large
```

``` r
# Pairwise comparisons using Dunn's test
pwc <- df %>% 
  dunn_test(`TotalCellCountFOV` ~ Treatment, p.adjust.method = "bonferroni") 
pwc
```

```
## # A tibble: 6 × 9
##   .y.               group1 group2    n1    n2 statistic          p     p.adj p.adj.signif
## * <chr>             <chr>  <chr>  <int> <int>     <dbl>      <dbl>     <dbl> <chr>       
## 1 TotalCellCountFOV CTRL   D1         9     9     2.84  0.00449    0.0270    *           
## 2 TotalCellCountFOV CTRL   D2         9     9     0.515 0.607      1         ns          
## 3 TotalCellCountFOV CTRL   D3         9     9     4.56  0.00000502 0.0000301 ****        
## 4 TotalCellCountFOV D1     D2         9     9    -2.33  0.0200     0.120     ns          
## 5 TotalCellCountFOV D1     D3         9     9     1.72  0.0850     0.510     ns          
## 6 TotalCellCountFOV D2     D3         9     9     4.05  0.0000514  0.000308  ***
```

``` r
# add N to plot
tab <- data.frame(xtabs(~ Treatment, data = df))
head(tab)
```

```
##   Treatment Freq
## 1      CTRL    9
## 2        D1    9
## 3        D2    9
## 4        D3    9
```

``` r
# Add cell number per cluster to cluster labels
Labels = paste0("n = ", tab$Freq, "")



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "Treatment")

plot_count_all <- ggplot(df, aes(x = Treatment, y = TotalCellCountFOV, fill = "Treatment"))+
  geom_boxplot(fill="white")+
  geom_beeswarm(aes(color = Treatment), size = 2, cex = 3)+
  scale_color_manual(values = cols_treat)+
  theme_classic() +
  theme(
      axis.text.x = element_text(angle = 45, size = 10, face = "bold", hjust = 1),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
      plot.margin = margin(0.5, 1, 0.5, 1, "cm"),
      legend.position = "none"
    ) +
  stat_pvalue_manual(pwc,
                       hide.ns = TRUE, size = 6,
                       step.increase = 0.1, y.position = 5200) +
  xlab(NULL)+
  ylab("Total count per FOV [#]")+
  ggtitle("Total cells")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,6800))+
  NoLegend()
# +
#     annotate(geom = 'text',
#            x="ILC2s",
#            y=28,
#            label=Labels[1], 
#            #angle = 90, 
#            size = 10/.pt)

# plot_count_all
```

## Total immune count IL-33


``` r
# filter for CTRL and convert to longer format
df <- df_lung %>%
  select(Treatment, Dataset, `Immune cells`) %>%
  mutate(Treatment = paste0("D", Treatment)) %>%
  mutate(Treatment = gsub("DCTRL", "CTRL", Treatment)) %>%
  mutate(Treatment = factor(Treatment, level =c(
    "CTRL", "D1", "D2", "D3"
  )))

# Testing for normal distribution
shapiro.test(df$`Immune cells`)
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  df$`Immune cells`
## W = 0.76496, p-value = 3.524e-06
```

``` r
# Kruskal-Wallis-test to check for significance between tested groups and effect size
res.kruskal <- df %>% kruskal_test(`Immune cells` ~ Treatment)
res.kruskal
```

```
## # A tibble: 1 × 6
##   .y.              n statistic    df         p method        
## * <chr>        <int>     <dbl> <int>     <dbl> <chr>         
## 1 Immune cells    36      28.9     3 0.0000023 Kruskal-Wallis
```

``` r
df %>% kruskal_effsize(`Immune cells` ~ Treatment)
```

```
## # A tibble: 1 × 5
##   .y.              n effsize method  magnitude
## * <chr>        <int>   <dbl> <chr>   <ord>    
## 1 Immune cells    36   0.811 eta2[H] large
```

``` r
# Pairwise comparisons using Dunn's test
pwc <- df %>% 
  dunn_test(`Immune cells` ~ Treatment, p.adjust.method = "bonferroni") 
pwc
```

```
## # A tibble: 6 × 9
##   .y.          group1 group2    n1    n2 statistic           p      p.adj p.adj.signif
## * <chr>        <chr>  <chr>  <int> <int>     <dbl>       <dbl>      <dbl> <chr>       
## 1 Immune cells CTRL   D1         9     9      3.19 0.00143     0.00858    **          
## 2 Immune cells CTRL   D2         9     9      1.44 0.149       0.894      ns          
## 3 Immune cells CTRL   D3         9     9      5.08 0.000000379 0.00000228 ****        
## 4 Immune cells D1     D2         9     9     -1.75 0.0809      0.486      ns          
## 5 Immune cells D1     D3         9     9      1.89 0.0587      0.352      ns          
## 6 Immune cells D2     D3         9     9      3.64 0.000277    0.00166    **
```

``` r
# add N to plot
tab <- data.frame(xtabs(~ Treatment, data = df))
head(tab)
```

```
##   Treatment Freq
## 1      CTRL    9
## 2        D1    9
## 3        D2    9
## 4        D3    9
```

``` r
# Add cell number per cluster to cluster labels
Labels = paste0("n = ", tab$Freq, "")



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "Treatment")

plot_count_immune <- ggplot(df, aes(x = Treatment, y = `Immune cells`, fill = "Treatment"))+
  geom_boxplot(fill="white")+
  geom_beeswarm(aes(color = Treatment), size = 2, cex = 3)+
  scale_color_manual(values = cols_treat)+
  theme_classic() +
  theme(
      axis.text.x = element_text(angle = 45, size = 10, face = "bold", hjust = 1),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
      plot.margin = margin(0.5, 1, 0.5, 1, "cm"),
      legend.position = "none"
    ) +
  stat_pvalue_manual(pwc,
                       hide.ns = TRUE, size = 6,
                       step.increase = 0.1, y.position = 1750) +
  xlab(NULL)+
  ylab("Total count per FOV [#]")+
  ggtitle("Immune cells")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,2300))+
  NoLegend()
# +
#     annotate(geom = 'text',
#            x="ILC2s",
#            y=28,
#            label=Labels[1], 
#            #angle = 90, 
#            size = 10/.pt)

# plot_count_immune
```

## Total ILC count IL-33


``` r
# filter for CTRL and convert to longer format
df <- df_lung %>%
  select(Treatment, Dataset, ILCs) %>%
  mutate(Treatment = paste0("D", Treatment)) %>%
  mutate(Treatment = gsub("DCTRL", "CTRL", Treatment)) %>%
  mutate(Treatment = factor(Treatment, level =c(
    "CTRL", "D1", "D2", "D3"
  )))

# Testing for normal distribution
shapiro.test(df$ILCs)
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  df$ILCs
## W = 0.74981, p-value = 1.892e-06
```

``` r
# Kruskal-Wallis-test to check for significance between tested groups and effect size
res.kruskal <- df %>% kruskal_test(`ILCs` ~ Treatment)
res.kruskal
```

```
## # A tibble: 1 × 6
##   .y.       n statistic    df         p method        
## * <chr> <int>     <dbl> <int>     <dbl> <chr>         
## 1 ILCs     36      25.9     3 0.0000102 Kruskal-Wallis
```

``` r
df %>% kruskal_effsize(`ILCs` ~ Treatment)
```

```
## # A tibble: 1 × 5
##   .y.       n effsize method  magnitude
## * <chr> <int>   <dbl> <chr>   <ord>    
## 1 ILCs     36   0.714 eta2[H] large
```

``` r
# Pairwise comparisons using Dunn's test
pwc <- df %>% 
  dunn_test(`ILCs` ~ Treatment, p.adjust.method = "bonferroni") 
pwc
```

```
## # A tibble: 6 × 9
##   .y.   group1 group2    n1    n2 statistic           p      p.adj p.adj.signif
## * <chr> <chr>  <chr>  <int> <int>     <dbl>       <dbl>      <dbl> <chr>       
## 1 ILCs  CTRL   D1         9     9     3.14  0.00166     0.00998    **          
## 2 ILCs  CTRL   D2         9     9     2.20  0.0275      0.165      ns          
## 3 ILCs  CTRL   D3         9     9     4.99  0.000000600 0.00000360 ****        
## 4 ILCs  D1     D2         9     9    -0.940 0.347       1          ns          
## 5 ILCs  D1     D3         9     9     1.85  0.0648      0.389      ns          
## 6 ILCs  D2     D3         9     9     2.79  0.00533     0.0320     *
```

``` r
# add N to plot
tab <- data.frame(xtabs(~ Treatment, data = df))
head(tab)
```

```
##   Treatment Freq
## 1      CTRL    9
## 2        D1    9
## 3        D2    9
## 4        D3    9
```

``` r
# Add cell number per cluster to cluster labels
Labels = paste0("n = ", tab$Freq, "")



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "Treatment")

plot_count_ilc <- ggplot(df, aes(x = Treatment, y = ILCs, fill = "Treatment"))+
  geom_boxplot(fill="white")+
  geom_beeswarm(aes(color = Treatment), size = 2, cex = 3)+
  scale_color_manual(values = cols_treat)+
  theme_classic() +
  theme(
      axis.text.x = element_text(angle = 45, size = 10, face = "bold", hjust = 1),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
      plot.margin = margin(0.5, 1, 0.5, 1, "cm"),
      legend.position = "none"
    ) +
  stat_pvalue_manual(pwc,
                       hide.ns = TRUE, size = 6,
                       step.increase = 0.1, y.position = 180) +
  xlab(NULL)+
  ylab("Total count per FOV [#]")+
  ggtitle("ILCs")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,250))+
  NoLegend()
# +
#     annotate(geom = 'text',
#            x="ILC2s",
#            y=28,
#            label=Labels[1], 
#            #angle = 90, 
#            size = 10/.pt)

# plot_count_ilc
```

## Total count ILC1s/NK cells


``` r
# filter for CTRL and convert to longer format
df <- df_lung %>%
  select(Treatment, Dataset, `NK cells/ILC1s`) %>%
  mutate(Treatment = paste0("D", Treatment)) %>%
  mutate(Treatment = gsub("DCTRL", "CTRL", Treatment)) %>%
  mutate(Treatment = factor(Treatment, level =c(
    "CTRL", "D1", "D2", "D3"
  )))

# Testing for normal distribution
shapiro.test(df$`NK cells/ILC1s`)
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  df$`NK cells/ILC1s`
## W = 0.83781, p-value = 0.000101
```

``` r
# Kruskal-Wallis-test to check for significance between tested groups and effect size
res.kruskal <- df %>% kruskal_test(`NK cells/ILC1s` ~ Treatment)
res.kruskal
```

```
## # A tibble: 1 × 6
##   .y.                n statistic    df        p method        
## * <chr>          <int>     <dbl> <int>    <dbl> <chr>         
## 1 NK cells/ILC1s    36      20.8     3 0.000118 Kruskal-Wallis
```

``` r
df %>% kruskal_effsize(`NK cells/ILC1s` ~ Treatment)
```

```
## # A tibble: 1 × 5
##   .y.                n effsize method  magnitude
## * <chr>          <int>   <dbl> <chr>   <ord>    
## 1 NK cells/ILC1s    36   0.555 eta2[H] large
```

``` r
# Pairwise comparisons using Dunn's test
pwc <- df %>% 
  dunn_test(`NK cells/ILC1s` ~ Treatment, p.adjust.method = "bonferroni") 
pwc
```

```
## # A tibble: 6 × 9
##   .y.            group1 group2    n1    n2 statistic         p    p.adj p.adj.signif
## * <chr>          <chr>  <chr>  <int> <int>     <dbl>     <dbl>    <dbl> <chr>       
## 1 NK cells/ILC1s CTRL   D1         9     9     3.69  0.000227  0.00136  **          
## 2 NK cells/ILC1s CTRL   D2         9     9     1.82  0.0687    0.412    ns          
## 3 NK cells/ILC1s CTRL   D3         9     9     4.02  0.0000573 0.000344 ***         
## 4 NK cells/ILC1s D1     D2         9     9    -1.87  0.0621    0.373    ns          
## 5 NK cells/ILC1s D1     D3         9     9     0.337 0.736     1        ns          
## 6 NK cells/ILC1s D2     D3         9     9     2.20  0.0276    0.166    ns
```

``` r
# add N to plot
tab <- data.frame(xtabs(~ Treatment, data = df))
head(tab)
```

```
##   Treatment Freq
## 1      CTRL    9
## 2        D1    9
## 3        D2    9
## 4        D3    9
```

``` r
# Add cell number per cluster to cluster labels
Labels = paste0("n = ", tab$Freq, "")



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "Treatment")

plot_count_ilc1 <- ggplot(df, aes(x = Treatment, y = `NK cells/ILC1s`, fill = "Treatment"))+
  geom_boxplot(fill="white")+
  geom_beeswarm(aes(color = Treatment), size = 2, cex = 3)+
  scale_color_manual(values = cols_treat)+
  theme_classic() +
  theme(
      axis.text.x = element_text(angle = 45, size = 10, face = "bold", hjust = 1),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
      plot.margin = margin(0.5, 1, 0.5, 1, "cm"),
      legend.position = "none"
    ) +
  stat_pvalue_manual(pwc,
                       hide.ns = TRUE, size = 6,
                       step.increase = 0.1, y.position = 65) +
  xlab(NULL)+
  ylab("Total count per FOV [#]")+
  ggtitle("NK cells/ILC1s")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,78))+
  NoLegend()
# +
#     annotate(geom = 'text',
#            x="ILC2s",
#            y=28,
#            label=Labels[1], 
#            #angle = 90, 
#            size = 10/.pt)

# plot_count_ilc1
```

## Total count ILC2s


``` r
# filter for CTRL and convert to longer format
df <- df_lung %>%
  select(Treatment, Dataset, `ILC2s`) %>%
  mutate(Treatment = paste0("D", Treatment)) %>%
  mutate(Treatment = gsub("DCTRL", "CTRL", Treatment)) %>%
  mutate(Treatment = factor(Treatment, level =c(
    "CTRL", "D1", "D2", "D3"
  )))

# Testing for normal distribution
shapiro.test(df$`ILC2s`)
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  df$ILC2s
## W = 0.72409, p-value = 6.914e-07
```

``` r
# Kruskal-Wallis-test to check for significance between tested groups and effect size
res.kruskal <- df %>% kruskal_test(`ILC2s` ~ Treatment)
res.kruskal
```

```
## # A tibble: 1 × 6
##   .y.       n statistic    df        p method        
## * <chr> <int>     <dbl> <int>    <dbl> <chr>         
## 1 ILC2s    36      17.5     3 0.000546 Kruskal-Wallis
```

``` r
df %>% kruskal_effsize(`ILC2s` ~ Treatment)
```

```
## # A tibble: 1 × 5
##   .y.       n effsize method  magnitude
## * <chr> <int>   <dbl> <chr>   <ord>    
## 1 ILC2s    36   0.455 eta2[H] large
```

``` r
# Pairwise comparisons using Dunn's test
pwc <- df %>% 
  dunn_test(`ILC2s` ~ Treatment, p.adjust.method = "bonferroni") 
pwc
```

```
## # A tibble: 6 × 9
##   .y.   group1 group2    n1    n2 statistic         p    p.adj p.adj.signif
## * <chr> <chr>  <chr>  <int> <int>     <dbl>     <dbl>    <dbl> <chr>       
## 1 ILC2s CTRL   D1         9     9     1.33  0.183     1        ns          
## 2 ILC2s CTRL   D2         9     9     0.839 0.401     1        ns          
## 3 ILC2s CTRL   D3         9     9     3.96  0.0000743 0.000446 ***         
## 4 ILC2s D1     D2         9     9    -0.492 0.622     1        ns          
## 5 ILC2s D1     D3         9     9     2.63  0.00853   0.0512   ns          
## 6 ILC2s D2     D3         9     9     3.12  0.00179   0.0108   *
```

``` r
# add N to plot
tab <- data.frame(xtabs(~ Treatment, data = df))
head(tab)
```

```
##   Treatment Freq
## 1      CTRL    9
## 2        D1    9
## 3        D2    9
## 4        D3    9
```

``` r
# Add cell number per cluster to cluster labels
Labels = paste0("n = ", tab$Freq, "")



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "Treatment")

plot_count_ilc2 <- ggplot(df, aes(x = Treatment, y = `ILC2s`, fill = "Treatment"))+
  geom_boxplot(fill="white")+
  geom_beeswarm(aes(color = Treatment), size = 2, cex = 3)+
  scale_color_manual(values = cols_treat)+
  theme_classic() +
  theme(
      axis.text.x = element_text(angle = 45, size = 10, face = "bold", hjust = 1),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
      plot.margin = margin(0.5, 1, 0.5, 1, "cm"),
      legend.position = "none"
    ) +
  stat_pvalue_manual(pwc,
                       hide.ns = TRUE, size = 6,
                       step.increase = 0.1, y.position = 180) +
  xlab(NULL)+
  ylab("Total count per FOV [#]")+
  ggtitle("ILC2s")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,250))+
  NoLegend()
# +
#     annotate(geom = 'text',
#            x="ILC2s",
#            y=28,
#            label=Labels[1], 
#            #angle = 90, 
#            size = 10/.pt)

# plot_count_ilc2
```

## Total count ILC3s


``` r
# filter for CTRL and convert to longer format
df <- df_lung %>%
  select(Treatment, Dataset, `ILC3s`) %>%
  mutate(Treatment = paste0("D", Treatment)) %>%
  mutate(Treatment = gsub("DCTRL", "CTRL", Treatment)) %>%
  mutate(Treatment = factor(Treatment, level =c(
    "CTRL", "D1", "D2", "D3"
  )))

# Testing for normal distribution
shapiro.test(df$`ILC3s`)
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  df$ILC3s
## W = 0.83426, p-value = 8.433e-05
```

``` r
# Kruskal-Wallis-test to check for significance between tested groups and effect size
res.kruskal <- df %>% kruskal_test(`ILC3s` ~ Treatment)
res.kruskal
```

```
## # A tibble: 1 × 6
##   .y.       n statistic    df     p method        
## * <chr> <int>     <dbl> <int> <dbl> <chr>         
## 1 ILC3s    36      12.4     3 0.006 Kruskal-Wallis
```

``` r
df %>% kruskal_effsize(`ILC3s` ~ Treatment)
```

```
## # A tibble: 1 × 5
##   .y.       n effsize method  magnitude
## * <chr> <int>   <dbl> <chr>   <ord>    
## 1 ILC3s    36   0.295 eta2[H] large
```

``` r
# Pairwise comparisons using Dunn's test
pwc <- df %>% 
  dunn_test(`ILC3s` ~ Treatment, p.adjust.method = "bonferroni") 
pwc
```

```
## # A tibble: 6 × 9
##   .y.   group1 group2    n1    n2 statistic        p   p.adj p.adj.signif
## * <chr> <chr>  <chr>  <int> <int>     <dbl>    <dbl>   <dbl> <chr>       
## 1 ILC3s CTRL   D1         9     9     2.01  0.0444   0.266   ns          
## 2 ILC3s CTRL   D2         9     9     2.73  0.00625  0.0375  *           
## 3 ILC3s CTRL   D3         9     9     3.30  0.000976 0.00585 **          
## 4 ILC3s D1     D2         9     9     0.724 0.469    1       ns          
## 5 ILC3s D1     D3         9     9     1.29  0.198    1       ns          
## 6 ILC3s D2     D3         9     9     0.563 0.573    1       ns
```

``` r
# add N to plot
tab <- data.frame(xtabs(~ Treatment, data = df))
head(tab)
```

```
##   Treatment Freq
## 1      CTRL    9
## 2        D1    9
## 3        D2    9
## 4        D3    9
```

``` r
# Add cell number per cluster to cluster labels
Labels = paste0("n = ", tab$Freq, "")



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "Treatment")

plot_count_ilc3 <- ggplot(df, aes(x = Treatment, y = `ILC3s`, fill = "Treatment"))+
  geom_boxplot(fill="white")+
  geom_beeswarm(aes(color = Treatment), size = 2, cex = 3)+
  scale_color_manual(values = cols_treat)+
  theme_classic() +
  theme(
      axis.text.x = element_text(angle = 45, size = 10, face = "bold", hjust = 1),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
      plot.margin = margin(0.5, 1, 0.5, 1, "cm"),
      legend.position = "none"
    ) +
  stat_pvalue_manual(pwc,
                       hide.ns = TRUE, size = 6,
                       step.increase = 0.1, y.position = 11) +
  xlab(NULL)+
  ylab("Total count per FOV [#]")+
  ggtitle("ILC3s")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,14))+
  NoLegend()
# +
#     annotate(geom = 'text',
#            x="ILC2s",
#            y=28,
#            label=Labels[1], 
#            #angle = 90, 
#            size = 10/.pt)

# plot_count_ilc3
```

## Freq ILC1s of ILC compartment


``` r
# filter for CTRL and convert to longer format
df <- df_lung %>%
  select(Treatment, Dataset, `Prop_NK cells/ILC1s_perTotalILCsFOV`) %>%
  mutate(Treatment = paste0("D", Treatment)) %>%
  mutate(Treatment = gsub("DCTRL", "CTRL", Treatment)) %>%
  mutate(
    Treatment = factor(Treatment, level =c(
    "CTRL", "D1", "D2", "D3"
  )),
  `NK cells/ILC1s` = `Prop_NK cells/ILC1s_perTotalILCsFOV`)

# Testing for normal distribution
shapiro.test(df$`NK cells/ILC1s`)
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  df$`NK cells/ILC1s`
## W = 0.96616, p-value = 0.33
```

``` r
# Kruskal-Wallis-test to check for significance between tested groups and effect size
res.kruskal <- df %>% kruskal_test(`NK cells/ILC1s` ~ Treatment)
res.kruskal
```

```
## # A tibble: 1 × 6
##   .y.                n statistic    df      p method        
## * <chr>          <int>     <dbl> <int>  <dbl> <chr>         
## 1 NK cells/ILC1s    36      7.56     3 0.0559 Kruskal-Wallis
```

``` r
df %>% kruskal_effsize(`NK cells/ILC1s` ~ Treatment)
```

```
## # A tibble: 1 × 5
##   .y.                n effsize method  magnitude
## * <chr>          <int>   <dbl> <chr>   <ord>    
## 1 NK cells/ILC1s    36   0.143 eta2[H] large
```

``` r
# Pairwise comparisons using Dunn's test
pwc <- df %>% 
  dunn_test(`NK cells/ILC1s` ~ Treatment, p.adjust.method = "bonferroni") 
pwc
```

```
## # A tibble: 6 × 9
##   .y.            group1 group2    n1    n2 statistic       p  p.adj p.adj.signif
## * <chr>          <chr>  <chr>  <int> <int>     <dbl>   <dbl>  <dbl> <chr>       
## 1 NK cells/ILC1s CTRL   D1         9     9     0.850 0.395   1      ns          
## 2 NK cells/ILC1s CTRL   D2         9     9     0.313 0.754   1      ns          
## 3 NK cells/ILC1s CTRL   D3         9     9    -1.75  0.0810  0.486  ns          
## 4 NK cells/ILC1s D1     D2         9     9    -0.537 0.591   1      ns          
## 5 NK cells/ILC1s D1     D3         9     9    -2.60  0.00945 0.0567 ns          
## 6 NK cells/ILC1s D2     D3         9     9    -2.06  0.0396  0.237  ns
```

``` r
# add N to plot
tab <- data.frame(xtabs(~ Treatment, data = df))
head(tab)
```

```
##   Treatment Freq
## 1      CTRL    9
## 2        D1    9
## 3        D2    9
## 4        D3    9
```

``` r
# Add cell number per cluster to cluster labels
Labels = paste0("n = ", tab$Freq, "")



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "Treatment")

plot_freq_ilc1 <- ggplot(df, aes(x = Treatment, y = `NK cells/ILC1s`, fill = "Treatment"))+
  geom_boxplot(fill="white")+
  geom_beeswarm(aes(color = Treatment), size = 2, cex = 3)+
  scale_color_manual(values = cols_treat)+
  theme_classic() +
  theme(
      axis.text.x = element_text(angle = 45, size = 10, face = "bold", hjust = 1),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
      plot.margin = margin(0.5, 1, 0.5, 1, "cm"),
      legend.position = "none"
    ) +
  stat_pvalue_manual(pwc,
                       hide.ns = TRUE, size = 6,
                       step.increase = 0.1, y.position = 14) +
  xlab(NULL)+
  ylab("Frequency per FOV [%]")+
  ggtitle("NK cells/ILC1s")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,85))+
  NoLegend()
# +
#     annotate(geom = 'text',
#            x="ILC2s",
#            y=28,
#            label=Labels[1], 
#            #angle = 90, 
#            size = 10/.pt)

# plot_freq_ilc1
```

## Freq ILC2s of ILC compartment


``` r
# filter for CTRL and convert to longer format
df <- df_lung %>%
  select(Treatment, Dataset, Prop_ILC2s_perTotalILCsFOV) %>%
  mutate(Treatment = paste0("D", Treatment)) %>%
  mutate(Treatment = gsub("DCTRL", "CTRL", Treatment)) %>%
  mutate(Treatment = factor(Treatment, level =c(
    "CTRL", "D1", "D2", "D3"
  )),
  ILC2s = Prop_ILC2s_perTotalILCsFOV)

# Testing for normal distribution
shapiro.test(df$`ILC2s`)
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  df$ILC2s
## W = 0.96896, p-value = 0.3978
```

``` r
# Kruskal-Wallis-test to check for significance between tested groups and effect size
res.kruskal <- df %>% kruskal_test(`ILC2s` ~ Treatment)
res.kruskal
```

```
## # A tibble: 1 × 6
##   .y.       n statistic    df      p method        
## * <chr> <int>     <dbl> <int>  <dbl> <chr>         
## 1 ILC2s    36      8.72     3 0.0332 Kruskal-Wallis
```

``` r
df %>% kruskal_effsize(`ILC2s` ~ Treatment)
```

```
## # A tibble: 1 × 5
##   .y.       n effsize method  magnitude
## * <chr> <int>   <dbl> <chr>   <ord>    
## 1 ILC2s    36   0.179 eta2[H] large
```

``` r
# Pairwise comparisons using Dunn's test
pwc <- df %>% 
  dunn_test(`ILC2s` ~ Treatment, p.adjust.method = "bonferroni") 
pwc
```

```
## # A tibble: 6 × 9
##   .y.   group1 group2    n1    n2 statistic       p  p.adj p.adj.signif
## * <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl>  <dbl> <chr>       
## 1 ILC2s CTRL   D1         9     9    -1.13  0.259   1      ns          
## 2 ILC2s CTRL   D2         9     9    -0.984 0.325   1      ns          
## 3 ILC2s CTRL   D3         9     9     1.49  0.137   0.821  ns          
## 4 ILC2s D1     D2         9     9     0.145 0.884   1      ns          
## 5 ILC2s D1     D3         9     9     2.62  0.00885 0.0531 ns          
## 6 ILC2s D2     D3         9     9     2.47  0.0134  0.0806 ns
```

``` r
# add N to plot
tab <- data.frame(xtabs(~ Treatment, data = df))
head(tab)
```

```
##   Treatment Freq
## 1      CTRL    9
## 2        D1    9
## 3        D2    9
## 4        D3    9
```

``` r
# Add cell number per cluster to cluster labels
Labels = paste0("n = ", tab$Freq, "")



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "Treatment")

plot_freq_ilc2 <- ggplot(df, aes(x = Treatment, y = `ILC2s`, fill = "Treatment"))+
  geom_boxplot(fill="white")+
  geom_beeswarm(aes(color = Treatment), size = 2, cex = 3)+
  scale_color_manual(values = cols_treat)+
  theme_classic() +
  theme(
      axis.text.x = element_text(angle = 45, size = 10, face = "bold", hjust = 1),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
      plot.margin = margin(0.5, 1, 0.5, 1, "cm"),
      legend.position = "none"
    ) +
  stat_pvalue_manual(pwc,
                       hide.ns = TRUE, size = 6,
                       step.increase = 0.1, y.position = 14) +
  xlab(NULL)+
  ylab("Frequency per FOV [%]")+
  ggtitle("ILC2s")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,90))+
  NoLegend()
# +
#     annotate(geom = 'text',
#            x="ILC2s",
#            y=28,
#            label=Labels[1], 
#            #angle = 90, 
#            size = 10/.pt)

# plot_freq_ilc2
```

## Freq ILC3s of ILC compartment


``` r
# filter for CTRL and convert to longer format
df <- df_lung %>%
  select(Treatment, Dataset, Prop_ILC3s_perTotalILCsFOV) %>%
  mutate(Treatment = paste0("D", Treatment)) %>%
  mutate(Treatment = gsub("DCTRL", "CTRL", Treatment)) %>%
  mutate(Treatment = factor(Treatment, level =c(
    "CTRL", "D1", "D2", "D3"
  )),
  ILC3s = Prop_ILC3s_perTotalILCsFOV)

# Testing for normal distribution
shapiro.test(df$`ILC3s`)
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  df$ILC3s
## W = 0.84165, p-value = 0.000123
```

``` r
# Kruskal-Wallis-test to check for significance between tested groups and effect size
res.kruskal <- df %>% kruskal_test(`ILC3s` ~ Treatment)
res.kruskal
```

```
## # A tibble: 1 × 6
##   .y.       n statistic    df       p method        
## * <chr> <int>     <dbl> <int>   <dbl> <chr>         
## 1 ILC3s    36      13.3     3 0.00411 Kruskal-Wallis
```

``` r
df %>% kruskal_effsize(`ILC3s` ~ Treatment)
```

```
## # A tibble: 1 × 5
##   .y.       n effsize method  magnitude
## * <chr> <int>   <dbl> <chr>   <ord>    
## 1 ILC3s    36   0.321 eta2[H] large
```

``` r
# Pairwise comparisons using Dunn's test
pwc <- df %>% 
  dunn_test(`ILC3s` ~ Treatment, p.adjust.method = "bonferroni") 
pwc
```

```
## # A tibble: 6 × 9
##   .y.   group1 group2    n1    n2 statistic        p   p.adj p.adj.signif
## * <chr> <chr>  <chr>  <int> <int>     <dbl>    <dbl>   <dbl> <chr>       
## 1 ILC3s CTRL   D1         9     9     2.13  0.0330   0.198   ns          
## 2 ILC3s CTRL   D2         9     9     3.60  0.000316 0.00189 **          
## 3 ILC3s CTRL   D3         9     9     1.61  0.108    0.648   ns          
## 4 ILC3s D1     D2         9     9     1.47  0.141    0.849   ns          
## 5 ILC3s D1     D3         9     9    -0.524 0.600    1       ns          
## 6 ILC3s D2     D3         9     9    -1.99  0.0461   0.276   ns
```

``` r
# add N to plot
tab <- data.frame(xtabs(~ Treatment, data = df))
head(tab)
```

```
##   Treatment Freq
## 1      CTRL    9
## 2        D1    9
## 3        D2    9
## 4        D3    9
```

``` r
# Add cell number per cluster to cluster labels
Labels = paste0("n = ", tab$Freq, "")



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "Treatment")

plot_freq_ilc3 <- ggplot(df, aes(x = Treatment, y = `ILC3s`, fill = "Treatment"))+
  geom_boxplot(fill="white")+
  geom_beeswarm(aes(color = Treatment), size = 2, cex = 3)+
  scale_color_manual(values = cols_treat)+
  theme_classic() +
  theme(
      axis.text.x = element_text(angle = 45, size = 10, face = "bold", hjust = 1),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
      plot.margin = margin(0.5, 1, 0.5, 1, "cm"),
      legend.position = "none"
    ) +
  stat_pvalue_manual(pwc,
                       hide.ns = TRUE, size = 6,
                       step.increase = 0.1, y.position = 14) +
  xlab(NULL)+
  ylab("Frequency per FOV [%]")+
  ggtitle("ILC3s")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,16))+
  NoLegend()
# +
#     annotate(geom = 'text',
#            x="ILC2s",
#            y=28,
#            label=Labels[1], 
#            #angle = 90, 
#            size = 10/.pt)

# plot_freq_ilc3
```

## Freq of immune cells, Endothelia & stroma and epithelia within total cells


``` r
# IMMUNE CELLS ---------------------------------------------------------------
df <- df_lung %>%
  select(Treatment, Dataset, `Prop_Immune cells_perTotalCountFOV`) %>%
  mutate(Treatment = paste0("D", Treatment)) %>%
  mutate(Treatment = gsub("DCTRL", "CTRL", Treatment)) %>%
  mutate(Treatment = factor(Treatment, level =c(
    "CTRL", "D1", "D2", "D3"
  )),
  `Immune cells` = `Prop_Immune cells_perTotalCountFOV`)

# Testing for normal distribution
shapiro.test(df$`Immune cells`)
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  df$`Immune cells`
## W = 0.92858, p-value = 0.02267
```

``` r
# Kruskal-Wallis-test to check for significance between tested groups and effect size
res.kruskal <- df %>% kruskal_test(`Immune cells` ~ Treatment)
res.kruskal
```

```
## # A tibble: 1 × 6
##   .y.              n statistic    df        p method        
## * <chr>        <int>     <dbl> <int>    <dbl> <chr>         
## 1 Immune cells    36      21.8     3 0.000071 Kruskal-Wallis
```

``` r
df %>% kruskal_effsize(`Immune cells` ~ Treatment)
```

```
## # A tibble: 1 × 5
##   .y.              n effsize method  magnitude
## * <chr>        <int>   <dbl> <chr>   <ord>    
## 1 Immune cells    36   0.588 eta2[H] large
```

``` r
# Pairwise comparisons using Dunn's test
pwc <- df %>% 
  dunn_test(`Immune cells` ~ Treatment, p.adjust.method = "bonferroni") 
pwc
```

```
## # A tibble: 6 × 9
##   .y.          group1 group2    n1    n2 statistic         p     p.adj p.adj.signif
## * <chr>        <chr>  <chr>  <int> <int>     <dbl>     <dbl>     <dbl> <chr>       
## 1 Immune cells CTRL   D1         9     9     0.839 0.401     1         ns          
## 2 Immune cells CTRL   D2         9     9     1.44  0.149     0.893     ns          
## 3 Immune cells CTRL   D3         9     9     4.39  0.0000115 0.0000690 ****        
## 4 Immune cells D1     D2         9     9     0.604 0.546     1         ns          
## 5 Immune cells D1     D3         9     9     3.55  0.000389  0.00233   **          
## 6 Immune cells D2     D3         9     9     2.94  0.00325   0.0195    *
```

``` r
# add N to plot
tab <- data.frame(xtabs(~ Treatment, data = df))
head(tab)
```

```
##   Treatment Freq
## 1      CTRL    9
## 2        D1    9
## 3        D2    9
## 4        D3    9
```

``` r
# Add cell number per cluster to cluster labels
Labels = paste0("n = ", tab$Freq, "")



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "Treatment")

plot_immune <- ggplot(df, aes(x = Treatment, y = `Immune cells`, fill = "Treatment"))+
  geom_boxplot(fill="white")+
  geom_beeswarm(aes(color = Treatment), size = 2, cex = 3)+
  scale_color_manual(values = cols_treat)+
  theme_classic() +
  theme(
      axis.text.x = element_text(angle = 45, size = 10, face = "bold", hjust = 1),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
      plot.margin = margin(0.5, 1, 0.5, 1, "cm"),
      legend.position = "none"
    ) +
  stat_pvalue_manual(pwc,
                       hide.ns = TRUE, size = 6,
                       step.increase = 0.5, y.position = 40) +
  xlab(NULL)+
  ylab("Frequency/total cells per FOV [#]")+
  ggtitle("Immune cells")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,100))+
  NoLegend()


# Endothelia & stroma ---------------------------------------------------------------
df <- df_lung %>%
  select(Treatment, Dataset, `Prop_Endothelia & stroma_perTotalCountFOV`) %>%
  mutate(Treatment = paste0("D", Treatment)) %>%
  mutate(Treatment = gsub("DCTRL", "CTRL", Treatment)) %>%
  mutate(Treatment = factor(Treatment, level =c(
    "CTRL", "D1", "D2", "D3"
  )),
  `cells_of_interest` = `Prop_Endothelia & stroma_perTotalCountFOV`)

# Testing for normal distribution
shapiro.test(df$`cells_of_interest`)
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  df$cells_of_interest
## W = 0.97166, p-value = 0.4725
```

``` r
# Kruskal-Wallis-test to check for significance between tested groups and effect size
res.kruskal <- df %>% kruskal_test(`cells_of_interest` ~ Treatment)
res.kruskal
```

```
## # A tibble: 1 × 6
##   .y.                   n statistic    df        p method        
## * <chr>             <int>     <dbl> <int>    <dbl> <chr>         
## 1 cells_of_interest    36      18.3     3 0.000374 Kruskal-Wallis
```

``` r
df %>% kruskal_effsize(`cells_of_interest` ~ Treatment)
```

```
## # A tibble: 1 × 5
##   .y.                   n effsize method  magnitude
## * <chr>             <int>   <dbl> <chr>   <ord>    
## 1 cells_of_interest    36   0.479 eta2[H] large
```

``` r
# Pairwise comparisons using Dunn's test
pwc <- df %>% 
  dunn_test(`cells_of_interest` ~ Treatment, p.adjust.method = "bonferroni") 
pwc
```

```
## # A tibble: 6 × 9
##   .y.               group1 group2    n1    n2 statistic         p    p.adj p.adj.signif
## * <chr>             <chr>  <chr>  <int> <int>     <dbl>     <dbl>    <dbl> <chr>       
## 1 cells_of_interest CTRL   D1         9     9    -2.37  0.0177    0.106    ns          
## 2 cells_of_interest CTRL   D2         9     9    -2.26  0.0238    0.143    ns          
## 3 cells_of_interest CTRL   D3         9     9    -4.27  0.0000192 0.000115 ***         
## 4 cells_of_interest D1     D2         9     9     0.112 0.911     1        ns          
## 5 cells_of_interest D1     D3         9     9    -1.90  0.0572    0.343    ns          
## 6 cells_of_interest D2     D3         9     9    -2.01  0.0440    0.264    ns
```

``` r
# add N to plot
tab <- data.frame(xtabs(~ Treatment, data = df))
head(tab)
```

```
##   Treatment Freq
## 1      CTRL    9
## 2        D1    9
## 3        D2    9
## 4        D3    9
```

``` r
# Add cell number per cluster to cluster labels
Labels = paste0("n = ", tab$Freq, "")



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "Treatment")

plot_stroma <- ggplot(df, aes(x = Treatment, y = `cells_of_interest`, fill = "Treatment"))+
  geom_boxplot(fill="white")+
  geom_beeswarm(aes(color = Treatment), size = 2, cex = 3)+
  scale_color_manual(values = cols_treat)+
  theme_classic() +
  theme(
      axis.text.x = element_text(angle = 45, size = 10, face = "bold", hjust = 1),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
      plot.margin = margin(0.5, 1, 0.5, 1, "cm"),
      legend.position = "none"
    ) +
  stat_pvalue_manual(pwc,
                       hide.ns = TRUE, size = 6,
                       step.increase = 0.2, y.position = 80) +
  xlab(NULL)+
  ylab("Frequency/total cells per FOV [#]")+
  ggtitle("Endothelia & stroma")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,100))+
  NoLegend()


# Epithelia ---------------------------------------------------------------
df <- df_lung %>%
  select(Treatment, Dataset, `Prop_Epithelia_perTotalCountFOV`) %>%
  mutate(Treatment = paste0("D", Treatment)) %>%
  mutate(Treatment = gsub("DCTRL", "CTRL", Treatment)) %>%
  mutate(Treatment = factor(Treatment, level =c(
    "CTRL", "D1", "D2", "D3"
  )),
  `cells_of_interest` = `Prop_Epithelia_perTotalCountFOV`)

# Testing for normal distribution
shapiro.test(df$`cells_of_interest`)
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  df$cells_of_interest
## W = 0.91559, p-value = 0.009365
```

``` r
# Kruskal-Wallis-test to check for significance between tested groups and effect size
res.kruskal <- df %>% kruskal_test(`cells_of_interest` ~ Treatment)
res.kruskal
```

```
## # A tibble: 1 × 6
##   .y.                   n statistic    df      p method        
## * <chr>             <int>     <dbl> <int>  <dbl> <chr>         
## 1 cells_of_interest    36      9.28     3 0.0257 Kruskal-Wallis
```

``` r
df %>% kruskal_effsize(`cells_of_interest` ~ Treatment)
```

```
## # A tibble: 1 × 5
##   .y.                   n effsize method  magnitude
## * <chr>             <int>   <dbl> <chr>   <ord>    
## 1 cells_of_interest    36   0.196 eta2[H] large
```

``` r
# Pairwise comparisons using Dunn's test
pwc <- df %>% 
  dunn_test(`cells_of_interest` ~ Treatment, p.adjust.method = "bonferroni") 
pwc
```

```
## # A tibble: 6 × 9
##   .y.               group1 group2    n1    n2 statistic       p  p.adj p.adj.signif
## * <chr>             <chr>  <chr>  <int> <int>     <dbl>   <dbl>  <dbl> <chr>       
## 1 cells_of_interest CTRL   D1         9     9     2.66  0.00774 0.0465 *           
## 2 cells_of_interest CTRL   D2         9     9     2.50  0.0126  0.0755 ns          
## 3 cells_of_interest CTRL   D3         9     9     2.23  0.0260  0.156  ns          
## 4 cells_of_interest D1     D2         9     9    -0.168 0.867   1      ns          
## 5 cells_of_interest D1     D3         9     9    -0.436 0.663   1      ns          
## 6 cells_of_interest D2     D3         9     9    -0.269 0.788   1      ns
```

``` r
# add N to plot
tab <- data.frame(xtabs(~ Treatment, data = df))
head(tab)
```

```
##   Treatment Freq
## 1      CTRL    9
## 2        D1    9
## 3        D2    9
## 4        D3    9
```

``` r
# Add cell number per cluster to cluster labels
Labels = paste0("n = ", tab$Freq, "")



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "Treatment")

plot_epithelia <- ggplot(df, aes(x = Treatment, y = `cells_of_interest`, fill = "Treatment"))+
  geom_boxplot(fill="white")+
  geom_beeswarm(aes(color = Treatment), size = 2, cex = 3)+
  scale_color_manual(values = cols_treat)+
  theme_classic() +
  theme(
      axis.text.x = element_text(angle = 45, size = 10, face = "bold", hjust = 1),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
      plot.margin = margin(0.5, 1, 0.5, 1, "cm"),
      legend.position = "none"
    ) +
  stat_pvalue_manual(pwc,
                       hide.ns = TRUE, size = 6,
                       step.increase = 0.15, y.position = 40) +
  xlab(NULL)+
  ylab("Frequency/total cells per FOV [#]")+
  ggtitle("Epithelia")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,100))+
  NoLegend()


plot_major <- ggarrange(plot_immune, plot_stroma, plot_epithelia, 
          ncol = 3, nrow = 1, 
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"))+
  theme(plot.margin = margin(0, 0.1, 0, 0, "cm"))
```

## Differential marker expression of ILCs across conditions


``` r
SO.lung <- readRDS(paste0(here::here("1_data_tidying", "Lung_SI_all_cells_all_ALs_files"), "/lung_all_cells_all_ALs.rds"))
dim(SO.lung)
```

```
## [1]    32 67537
```

``` r
SO.lung$AL1 <- gsub("Vessels", "Stromal cells", SO.lung$AL1)

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


``` r
library(Seurat)
library(ComplexHeatmap)
library(ggplotify) # Required for as.ggplot()

set.seed(8)

# --- 1. Subsetting ---
SO.sub <- subset(SO.lung, subset = AL3 %in% c("ILC2s"))
SO.sub$AL3 <- droplevels(as.factor(SO.sub$AL3))

# --- 2. Direct Pseudo-bulk Aggregation ---
ilc_markers <- c("KLRG1", #"CD127", "CD90",
                 "ICOS", 
                 "MHCII", "CD44", "Ki67")

# Calculate the mean directly from Seurat
avg_exp <- AverageExpression(SO.sub, features = ilc_markers, slot = "counts", group.by = "Treatment")
mat_annotated <- avg_exp[[DefaultAssay(SO.sub)]]
# --- 3. Manually Scale (Z-score) the Matrix ---
mat_scaled <- t(scale(t(mat_annotated)))
colnames(mat_scaled) <- gsub("g", "D", colnames(mat_scaled))

# Clean up math errors and cap outliers
mat_scaled[is.na(mat_scaled)] <- 0
mat_scaled[mat_scaled > 2] <- 2
mat_scaled[mat_scaled < -2] <- -2

# --- 4. Plot Heatmap using ComplexHeatmap ---
plot_heat_ilc2s <- ComplexHeatmap::pheatmap(
  mat = mat_scaled, 
  scale = "none", 
  clustering_method = "ward.D2",
  color = colorRampPalette(c("#648FFF", "white", "#FFB000"))(101), 
  breaks = seq(-2, 2, length.out = 102), 
  display_numbers = round(mat_scaled, 2), 
  number_color = "black",
  treeheight_col = 10,
  treeheight_row = 20, 
  name = "Z-Score" 
)

# --- 5. Convert to ggplot object for arrangement ---
# Using grid.grabExpr(draw()) ensures that all the internal sizing, 
# legends, and dendrograms from ComplexHeatmap are captured perfectly.
gg_heat_ilc2s <- as.ggplot(grid::grid.grabExpr(ComplexHeatmap::draw(plot_heat_ilc2s, heatmap_legend_side = "bottom")))+
  ggtitle("ILC2s\nacross conditions\n")+
  theme(
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 9),
        plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
        plot.margin = margin(0.1, 0.25, 0, 0.25, "cm"),
        legend.position = "bottom"
      )

# 'gg_heat' is now a standard ggplot object! 
# You can now combine it using patchwork, e.g.:
# combined_plot <- spatial_scatter_plot + gg_heat
# print(gg_heat_ilc2s)
```


``` r
set.seed(8)

# --- 1. Subsetting ---
SO.sub <- subset(SO.lung, subset = AL3 %in% c("NK cells/ILC1s"))
SO.sub$AL3 <- droplevels(as.factor(SO.sub$AL3))

# --- 2. Direct Pseudo-bulk Aggregation ---
ilc_markers <- c("KLRG1", "CD127", "CD90",
                 "ICOS", "NKp46", "NK11",
                 "MHCII", "CD44", "Ki67")

# Calculate the mean directly from Seurat
avg_exp <- AverageExpression(SO.sub, features = ilc_markers, slot = "counts", group.by = "Treatment")
mat_annotated <- avg_exp[[DefaultAssay(SO.sub)]]
# --- 3. Manually Scale (Z-score) the Matrix ---
mat_scaled <- t(scale(t(mat_annotated)))
colnames(mat_scaled) <- gsub("g", "D", colnames(mat_scaled))

# Clean up math errors and cap outliers
mat_scaled[is.na(mat_scaled)] <- 0
mat_scaled[mat_scaled > 2] <- 2
mat_scaled[mat_scaled < -2] <- -2

# --- 4. Plot Heatmap using ComplexHeatmap ---
plot_heat_ilc1s <- ComplexHeatmap::pheatmap(
  mat = mat_scaled, 
  scale = "none", 
  clustering_method = "ward.D2",
  color = colorRampPalette(c("#648FFF", "white", "#FFB000"))(101), 
  breaks = seq(-2, 2, length.out = 102), 
  display_numbers = round(mat_scaled, 2), 
  number_color = "black",
  treeheight_col = 10,
  treeheight_row = 20, 
  name = "Z-Score" 
)

# --- 5. Convert to ggplot object for arrangement ---
# Using grid.grabExpr(draw()) ensures that all the internal sizing, 
# legends, and dendrograms from ComplexHeatmap are captured perfectly.
gg_heat_ilc1s <- as.ggplot(grid::grid.grabExpr(ComplexHeatmap::draw(plot_heat_ilc1s, heatmap_legend_side = "bottom")))+
  ggtitle("NK cells/ILC1s\nacross conditions\n")+
  theme(
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 9),
        plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
        plot.margin = margin(0.1, 0.25, 0, 0.25, "cm"),
        legend.position = "bottom"
      )

# 'gg_heat' is now a standard ggplot object! 
# You can now combine it using patchwork, e.g.:
# combined_plot <- spatial_scatter_plot + gg_heat
# print(gg_heat_ilc1s)


# --- 1. Subsetting ---
SO.sub <- subset(SO.lung, subset = AL3 %in% c("ILC3s"))
SO.sub$AL3 <- droplevels(as.factor(SO.sub$AL3))

# --- 2. Direct Pseudo-bulk Aggregation ---
ilc_markers <- c("KLRG1", "CD127", "CD90", "CD3", "CD4",
                 "ICOS", "NKp46", "NK11",
                 "MHCII", "CD44", "Ki67")

# Calculate the mean directly from Seurat
avg_exp <- AverageExpression(SO.sub, features = ilc_markers, slot = "counts", group.by = "Treatment")
mat_annotated <- avg_exp[[DefaultAssay(SO.sub)]]
# --- 3. Manually Scale (Z-score) the Matrix ---
mat_scaled <- t(scale(t(mat_annotated)))
colnames(mat_scaled) <- gsub("g", "D", colnames(mat_scaled))

# Clean up math errors and cap outliers
mat_scaled[is.na(mat_scaled)] <- 0
mat_scaled[mat_scaled > 2] <- 2
mat_scaled[mat_scaled < -2] <- -2

# --- 4. Plot Heatmap using ComplexHeatmap ---
plot_heat_ilc3s <- ComplexHeatmap::pheatmap(
  mat = mat_scaled, 
  scale = "none", 
  clustering_method = "ward.D2",
  color = colorRampPalette(c("#648FFF", "white", "#FFB000"))(101), 
  breaks = seq(-2, 2, length.out = 102), 
  display_numbers = round(mat_scaled, 2), 
  number_color = "black",
  treeheight_col = 10,
  treeheight_row = 20, 
  name = "Z-Score" 
)

# --- 5. Convert to ggplot object for arrangement ---
# Using grid.grabExpr(draw()) ensures that all the internal sizing, 
# legends, and dendrograms from ComplexHeatmap are captured perfectly.
gg_heat_ilc3s <- as.ggplot(grid::grid.grabExpr(ComplexHeatmap::draw(plot_heat_ilc3s, heatmap_legend_side = "bottom")))+
  ggtitle("ILC3s\nacross conditions\n")+
  theme(
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 9),
        plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
        plot.margin = margin(0.1, 0.25, 0, 0.25, "cm"),
        legend.position = "bottom"
      )

# 'gg_heat' is now a standard ggplot object! 
# You can now combine it using patchwork, e.g.:
# combined_plot <- spatial_scatter_plot + gg_heat
# print(gg_heat_ilc3s)


# --- 1. Subsetting ---
SO.sub <- subset(SO.lung, subset = AL3 %in% c("ILC3s"))
SO.sub$AL3 <- droplevels(as.factor(SO.sub$AL3))

# --- 2. Direct Pseudo-bulk Aggregation ---
ilc_markers <- c("ICOS", "CD11c", "CD68", "SiglecF", "CD4",
                 "MHCII", "CD44", "Ki67")

# Calculate the mean directly from Seurat
avg_exp <- AverageExpression(SO.sub, features = ilc_markers, slot = "counts", group.by = "Treatment")
mat_annotated <- avg_exp[[DefaultAssay(SO.sub)]]
# --- 3. Manually Scale (Z-score) the Matrix ---
mat_scaled <- t(scale(t(mat_annotated)))
colnames(mat_scaled) <- gsub("g", "D", colnames(mat_scaled))

# Clean up math errors and cap outliers
mat_scaled[is.na(mat_scaled)] <- 0
mat_scaled[mat_scaled > 2] <- 2
mat_scaled[mat_scaled < -2] <- -2

# --- 4. Plot Heatmap using ComplexHeatmap ---
plot_heat_my <- ComplexHeatmap::pheatmap(
  mat = mat_scaled, 
  scale = "none", 
  clustering_method = "ward.D2",
  color = colorRampPalette(c("#648FFF", "white", "#FFB000"))(101), 
  breaks = seq(-2, 2, length.out = 102), 
  display_numbers = round(mat_scaled, 2), 
  number_color = "black",
  treeheight_col = 10,
  treeheight_row = 20, 
  name = "Z-Score" 
)

# --- 5. Convert to ggplot object for arrangement ---
# Using grid.grabExpr(draw()) ensures that all the internal sizing, 
# legends, and dendrograms from ComplexHeatmap are captured perfectly.
gg_heat_my <- as.ggplot(grid::grid.grabExpr(ComplexHeatmap::draw(plot_heat_my, heatmap_legend_side = "bottom")))+
  ggtitle("Myeloid cells\nacross conditions\n")+
  theme(
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 9),
        plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
        plot.margin = margin(0.1, 0.25, 0, 0.25, "cm"),
        legend.position = "bottom"
      )

# 'gg_heat' is now a standard ggplot object! 
# You can now combine it using patchwork, e.g.:
# combined_plot <- spatial_scatter_plot + gg_heat
# print(gg_heat_my)
```


``` r
set.seed(8)

# --- 1. Subsetting ---
SO.sub <- subset(SO.lung, subset = AL3 %in% c("NK cells/ILC1s", "ILC2s", "ILC3s"))
SO.sub$AL3 <- droplevels(as.factor(SO.sub$AL3))

SO.sub <- subset(SO.sub, subset = Treatment %in% c("CTRL"))
SO.sub$Treatment <- droplevels(as.factor(SO.sub$Treatment))

# --- 2. Direct Pseudo-bulk Aggregation ---
ilc_markers <- c("KLRG1",   
                 "ICOS", 
                 "MHCII", "CD44", "Ki67", "NKp46")

# Calculate the mean directly from Seurat
avg_exp <- AverageExpression(SO.sub, features = ilc_markers, slot = "counts", group.by = "AL3")
mat_annotated <- avg_exp[[DefaultAssay(SO.sub)]]
# --- 3. Manually Scale (Z-score) the Matrix ---
mat_scaled <- t(scale(t(mat_annotated)))
colnames(mat_scaled) <- gsub("g", "D", colnames(mat_scaled))

# Clean up math errors and cap outliers
mat_scaled[is.na(mat_scaled)] <- 0
mat_scaled[mat_scaled > 2] <- 2
mat_scaled[mat_scaled < -2] <- -2

# --- 4. Plot Heatmap using ComplexHeatmap ---
plot_heat_ctrl <- ComplexHeatmap::pheatmap(
  mat = mat_scaled, 
  scale = "none", 
  clustering_method = "ward.D2",
  color = colorRampPalette(c("#648FFF", "white", "#FFB000"))(101), 
  breaks = seq(-2, 2, length.out = 102), 
  display_numbers = round(mat_scaled, 2), 
  number_color = "black",
  treeheight_col = 10,
  treeheight_row = 20, 
  name = "Z-Score" 
)

# --- 5. Convert to ggplot object for arrangement ---
# Using grid.grabExpr(draw()) ensures that all the internal sizing, 
# legends, and dendrograms from ComplexHeatmap are captured perfectly.
gg_heat_ctrl <- as.ggplot(grid::grid.grabExpr(ComplexHeatmap::draw(plot_heat_ctrl, heatmap_legend_side = "right")))+
  ggtitle("ILC subtypes\n@ CTRL\n")+
  theme(
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 9),
        plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
        plot.margin = margin(0.1, 0.5, 0, 0.5, "cm"),
        legend.position = "bottom"
      )

# 'gg_heat' is now a standard ggplot object! 
# You can now combine it using patchwork, e.g.:
# combined_plot <- spatial_scatter_plot + gg_heat


# --- 1. Subsetting ---
SO.sub <- subset(SO.lung, subset = AL3 %in% c("NK cells/ILC1s", "ILC2s", "ILC3s"))
SO.sub$AL3 <- droplevels(as.factor(SO.sub$AL3))

SO.sub <- subset(SO.sub, subset = Treatment %in% c("3"))
SO.sub$Treatment <- droplevels(as.factor(SO.sub$Treatment))

# Calculate the mean directly from Seurat
avg_exp <- AverageExpression(SO.sub, features = ilc_markers, slot = "counts", group.by = "AL3")
mat_annotated <- avg_exp[[DefaultAssay(SO.sub)]]
# --- 3. Manually Scale (Z-score) the Matrix ---
mat_scaled <- t(scale(t(mat_annotated)))
colnames(mat_scaled) <- gsub("g", "D", colnames(mat_scaled))

# Clean up math errors and cap outliers
mat_scaled[is.na(mat_scaled)] <- 0
mat_scaled[mat_scaled > 2] <- 2
mat_scaled[mat_scaled < -2] <- -2

# --- 4. Plot Heatmap using ComplexHeatmap ---
plot_heat_d3 <- ComplexHeatmap::pheatmap(
  mat = mat_scaled, 
  scale = "none", 
  clustering_method = "ward.D2",
  color = colorRampPalette(c("#648FFF", "white", "#FFB000"))(101), 
  breaks = seq(-2, 2, length.out = 102), 
  display_numbers = round(mat_scaled, 2), 
  number_color = "black",
  treeheight_col = 10,
  treeheight_row = 20, 
  name = "Z-Score" 
)

# --- 5. Convert to ggplot object for arrangement ---
# Using grid.grabExpr(draw()) ensures that all the internal sizing, 
# legends, and dendrograms from ComplexHeatmap are captured perfectly.
gg_heat_d3 <- as.ggplot(grid::grid.grabExpr(ComplexHeatmap::draw(plot_heat_d3, heatmap_legend_side = "right")))+
  ggtitle("ILC subtypes\n@ D3\n")+
  theme(
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 9),
        plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
        plot.margin = margin(0.1, 0.5, 0, 0.5, "cm"),
        legend.position = "bottom"
      )

# 'gg_heat' is now a standard ggplot object! 
# You can now combine it using patchwork, e.g.:
# combined_plot <- spatial_scatter_plot + gg_heat

# ggarrange(gg_heat_ctrl, gg_heat_d3, ncol = 2)
```

## Combine plots for figure


``` r
top_figure <- ggarrange(plot_prop, plot_freq_immune, plot_freq, 
          ncol = 3, nrow = 1, 
          labels = c("AUTO"), 
          label.x = 0.1)+
  theme(plot.margin = margin(0, 0, 0.2, 0, "cm"))


middle_figure_right <- ggarrange(
          plot_count_all, plot_count_immune, plot_count_ilc,
          plot_count_ilc2, 
          ncol = 2, nrow = 2, 
          labels = c("D", "E", "F",
                     "G"),
          label.x = 0.1)+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))

middle_figure <- ggarrange(
          middle_figure_right,
          gg_heat_ilc2s, 
          ncol = 2, nrow = 1, 
          widths = c(2, 1),
          labels = c("", "H"),
          label.x = 0.1)+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))


final_figure <- ggarrange(top_figure, middle_figure, 
          ncol = 1, nrow = 2, heights = c(3, 4.5))+
  theme(plot.margin = margin(0, 0.1, 0, 0, "cm"))

final_figure
```

<embed src="D:/Repositories/2025_Kroh_et_al/Manuscript/EJI_submissions_for_publication/Figures_pdfs/Figure_5_Main-1.pdf" width="100%" style="display: block; margin: auto;" type="application/pdf" />

``` r
# annotate_figure(final_figure,
#                bottom = text_grob("N: CTRL/D1/D2 = 9 FOVs, D3 = 8 FOVs; Kruskal-Wallis-test & Dunn’s test; Significance levels of adjusted p-values: **** = 1e-04, *** = 0.001, ** = 0.01, * = 0.05", 
#                                   color = "black", 
#                                   # hjust = -1,
#                                   face = "italic", size = 8))
```


``` r
suppl_figure <- ggarrange(plot_count_ilc1, plot_count_ilc3, plot_freq_immune_d3, plot_freq_d3,
          ncol = 2, nrow = 2, 
          heights = c(3, 4),
          labels = c("AUTO"),
          label.x = 0.1)+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))

suppl_figure
```

<embed src="D:/Repositories/2025_Kroh_et_al/Manuscript/EJI_submissions_for_publication/Figures_pdfs/Suppl_Figure_9-1.pdf" width="100%" style="display: block; margin: auto;" type="application/pdf" />

``` r
# annotate_figure(suppl_figure,
#                bottom = text_grob("N: CTRL/D1/D2 = 9 FOVs, D3 = 8 FOVs; Kruskal-Wallis-test & Dunn’s test; Significance levels of adjusted p-values: **** = 1e-04, *** = 0.001, ** = 0.01, * = 0.05",
#                                   color = "black",
#                                   # hjust = -1,
#                                   face = "italic", size = 8))
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
##  [1] ggplotify_0.1.3       ComplexHeatmap_2.26.1 ggpubr_0.6.2          readr_2.1.6           ggbeeswarm_0.7.3      rstatix_0.7.3         ggplot2_4.0.1         dplyr_1.1.4           Seurat_5.3.1          SeuratObject_5.2.0    sp_2.2-0             
## 
## loaded via a namespace (and not attached):
##   [1] RcppAnnoy_0.0.22       splines_4.5.2          later_1.4.4            tibble_3.3.0           polyclip_1.10-7        fastDummies_1.7.5      lifecycle_1.0.5        doParallel_1.0.17      rprojroot_2.1.1        globals_0.19.0         lattice_0.22-7         vroom_1.7.0            MASS_7.3-65            backports_1.5.0        magrittr_2.0.4         plotly_4.12.0          sass_0.4.10            rmarkdown_2.30         jquerylib_0.1.4        yaml_2.3.10            httpuv_1.6.16          otel_0.2.0             sctransform_0.4.2      spam_2.11-1            spatstat.sparse_3.1-0  reticulate_1.44.0      cowplot_1.2.0          pbapply_1.7-4          RColorBrewer_1.1-3     abind_1.4-8            Rtsne_0.17             purrr_1.2.0            BiocGenerics_0.56.0    yulab.utils_0.2.4      rappdirs_0.3.4         circlize_0.4.17        IRanges_2.44.0         S4Vectors_0.48.0       ggrepel_0.9.6          irlba_2.3.5.1          listenv_0.10.0         spatstat.utils_3.2-0   moments_0.14.1         goftest_1.2-3          RSpectra_0.16-2        spatstat.random_3.4-2  fitdistrplus_1.2-6     parallelly_1.45.1      codetools_0.2-20       tidyselect_1.2.1       shape_1.4.6.1         
##  [52] farver_2.1.2           matrixStats_1.5.0      stats4_4.5.2           spatstat.explore_3.5-3 jsonlite_2.0.0         GetoptLong_1.1.0       progressr_0.18.0       Formula_1.2-5          ggridges_0.5.7         survival_3.8-3         iterators_1.0.14       foreach_1.5.2          tools_4.5.2            ica_1.0-3              Rcpp_1.1.0             glue_1.8.0             gridExtra_2.3          xfun_0.56              here_1.0.2             withr_3.0.2            fastmap_1.2.0          digest_0.6.38          gridGraphics_0.5-1     R6_2.6.1               mime_0.13              colorspace_2.1-2       Cairo_1.7-0            scattermore_1.2        tensor_1.5.1           dichromat_2.0-0.1      spatstat.data_3.1-9    utf8_1.2.6             tidyr_1.3.1            generics_0.1.4         data.table_1.17.8      httr_1.4.8             htmlwidgets_1.6.4      uwot_0.2.4             pkgconfig_2.0.3        gtable_0.3.6           lmtest_0.9-40          S7_0.2.0               htmltools_0.5.8.1      carData_3.0-6          dotCall64_1.2          clue_0.3-67            scales_1.4.0           png_0.1-8              spatstat.univar_3.1-4  knitr_1.51             rstudioapi_0.18.0     
## [103] tzdb_0.5.0             reshape2_1.4.5         rjson_0.2.23           nlme_3.1-168           cachem_1.1.0           zoo_1.8-14             GlobalOptions_0.1.3    stringr_1.6.0          KernSmooth_2.23-26     parallel_4.5.2         miniUI_0.1.2           vipor_0.4.7            pillar_1.11.1          vctrs_0.6.5            RANN_2.6.2             promises_1.5.0         car_3.1-5              xtable_1.8-8           cluster_2.1.8.1        beeswarm_0.4.0         evaluate_1.0.5         magick_2.9.0           cli_3.6.5              compiler_4.5.2         rlang_1.1.6            crayon_1.5.3           future.apply_1.20.2    ggsignif_0.6.4         labeling_0.4.3         fs_1.6.6               plyr_1.8.9             stringi_1.8.7          viridisLite_0.4.2      deldir_2.0-4           lazyeval_0.2.2         spatstat.geom_3.6-0    Matrix_1.7-4           RcppHNSW_0.6.0         hms_1.1.4              patchwork_1.3.2        bit64_4.6.0-1          future_1.69.0          shiny_1.13.0           ROCR_1.0-12            igraph_2.2.1           broom_1.0.12           bslib_0.10.0           bit_4.6.0
```
