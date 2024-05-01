Microviz_PCOA
================
2024-04-30

``` r
library(tidyverse) ; packageVersion("tidyverse") # 1.3.2
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.0     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

    ## [1] '2.0.0'

``` r
library(phyloseq) ; packageVersion("phyloseq") # 1.42.0
```

    ## [1] '1.46.0'

``` r
library("pairwiseAdonis"); packageVersion("pairwiseAdonis") # 0.4
```

    ## Loading required package: vegan
    ## Loading required package: permute
    ## Loading required package: lattice
    ## This is vegan 2.6-4
    ## Loading required package: cluster

    ## [1] '0.4.1'

``` r
library(ggpubr)
library(microbiome)
```

    ## 
    ## microbiome R package (microbiome.github.com)
    ##     
    ## 
    ## 
    ##  Copyright (C) 2011-2022 Leo Lahti, 
    ##     Sudarshan Shetty et al. <microbiome.github.io>
    ## 
    ## 
    ## Attaching package: 'microbiome'
    ## 
    ## The following object is masked from 'package:vegan':
    ## 
    ##     diversity
    ## 
    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     alpha
    ## 
    ## The following object is masked from 'package:base':
    ## 
    ##     transform

``` r
library(ggplot2)
library(qiime2R)
```

    ## 
    ## Attaching package: 'qiime2R'
    ## 
    ## The following object is masked from 'package:ggpubr':
    ## 
    ##     mean_sd

``` r
library(microViz)
```

    ## microViz version 0.12.1 - Copyright (C) 2021-2024 David Barnett
    ## Attaching package: 'microViz'
    ## The following object is masked from 'package:ggpubr':
    ## 
    ## stat_chull
    ## ! Website: https://david-barnett.github.io/microViz
    ## ✔ Useful?  For citation details, run: `citation("microViz")`
    ## ✖ Silence? `suppressPackageStartupMessages(library(microViz))`

### 

# Sterling Butler

# NOAA/ UM

# RICA Nursery Project - qiime2 processed abundance PCOA with microViz

# 3/2024

### 

``` r
#set working directory
setwd("/Users/sterling.butler/Desktop/RICA_fastq/Merged/feature")
```

``` r
#import the feature table using the read_q2biom()
asv<- read_q2biom("table.biom")

#import the taxa file as a qiime2 artifact and then convert it into a usable tax table for phyloseq
setwd("/Users/sterling.butler/Desktop/RICA_fastq/Merged")
tax <- read_qza("taxonomy.qza")
taxa <- tax$data %>% as_tibble() %>% 
  separate(Taxon, sep=";", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
  select(-Confidence) %>% arrange(Feature.ID) %>% mutate(ASV = 1:n()) %>% 
  mutate(newcol = "ASV") %>%
  unite("ASVs", newcol:ASV)
```

    ## Warning: Expected 7 pieces. Missing pieces filled with `NA` in 1200 rows [2, 3, 5, 7, 8,
    ## 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, ...].

``` r
#make the taxa table into a data frame 
taxa<-as.data.frame(taxa)

#import the meta data
meta_df<-read.csv("qiime203112024_16S_RICA_META.csv")
```

``` r
#make the correct column names for the phyloseq object 
taxa1 = taxa %>% tibble::column_to_rownames("Feature.ID")
meta = meta_df %>% tibble::column_to_rownames("sample")

#convert the taxa table back into a matrix, so phyloseq can read
tax_mat = as.matrix(taxa1)

#make the individual phyloseq objects otu_table, tax_table, sample_data
otu<-otu_table(asv, taxa_are_rows = TRUE)
tax_table<-tax_table(tax_mat)
sample<-sample_data(meta)
```

``` r
#Check to see if the labels for samples and ASVs for the different 
sample_names(sample)
```

    ##  [1] "RICA2023-7"      "RICA2023-11"     "RICA2023-88"     "RICA2023-94"    
    ##  [5] "RICA2023-23"     "RICA2023-38"     "RICA2023-112"    "RICA2023-120"   
    ##  [9] "RICA2023-30"     "RICA2023-31"     "RICA2023-85"     "RICA2023-86"    
    ## [13] "RICA2023-50"     "RICA2023-55"     "RICA2023-141"    "RICA2023-160"   
    ## [17] "RICA2023-2"      "RICA2023-6"      "RICA2023-104"    "RICA2023-116"   
    ## [21] "RICA2023-18"     "RICA2023-19"     "RICA2023-84"     "RICA2023-89"    
    ## [25] "RICA2023-8"      "RICA2023-9"      "RICA2023-93"     "RICA2023-95"    
    ## [29] "RICA2023-57"     "RICA2023-60"     "RICA2023-154"    "RICA2023-155"   
    ## [33] "RICA2023-22"     "RICA2023-39"     "RICA2023-110"    "RICA2023-113"   
    ## [37] "RICA2023-24"     "RICA2023-32"     "RICA2023-114"    "RICA2023-118"   
    ## [41] "RICA2023-73"     "RICA2023-80"     "RICA2023-124"    "RICA2023-131"   
    ## [45] "RICA2023-21"     "RICA2023-37"     "RICA2023-108"    "RICA2023-109"   
    ## [49] "RICA2023-NTCP2"  "RICA2023-NTCP3"  "RICA2023-NTCPCR" "RICA2023-46"    
    ## [53] "RICA2023-47"     "RICA2023-145"    "RICA2023-157"    "RICA2023-5"     
    ## [57] "RICA2023-12"     "RICA2023-98"     "RICA2023-100"    "RICA2023-45"    
    ## [61] "RICA2023-56"     "RICA2023-156"    "RICA2023-29"     "RICA2023-40"    
    ## [65] "RICA2023-61"     "RICA2023-69"     "RICA2023-101"    "RICA2023-111"   
    ## [69] "RICA2023-128"    "RICA2023-129"    "RICA2023-63"     "RICA2023-70"    
    ## [73] "RICA2023-122"    "RICA2023-130"    "RICA2023-71"     "RICA2023-72"    
    ## [77] "RICA2023-132"    "RICA2023-133"    "RICA2023-42"     "RICA2023-51"    
    ## [81] "RICA2023-150"    "RICA2023-153"    "RICA2023-3"      "RICA2023-14"    
    ## [85] "RICA2023-82"     "RICA2023-83"     "RICA2023-65"     "RICA2023-68"    
    ## [89] "RICA2023-126"    "RICA2023-127"    "RICA2023-58"     "RICA2023-59"    
    ## [93] "RICA2023-146"    "RICA2023-149"

``` r
sample_names(otu)
```

    ##  [1] "RICA2023-100"    "RICA2023-101"    "RICA2023-104"    "RICA2023-108"   
    ##  [5] "RICA2023-109"    "RICA2023-11"     "RICA2023-110"    "RICA2023-111"   
    ##  [9] "RICA2023-112"    "RICA2023-113"    "RICA2023-114"    "RICA2023-116"   
    ## [13] "RICA2023-118"    "RICA2023-12"     "RICA2023-120"    "RICA2023-122"   
    ## [17] "RICA2023-124"    "RICA2023-126"    "RICA2023-127"    "RICA2023-128"   
    ## [21] "RICA2023-129"    "RICA2023-130"    "RICA2023-131"    "RICA2023-132"   
    ## [25] "RICA2023-133"    "RICA2023-14"     "RICA2023-141"    "RICA2023-145"   
    ## [29] "RICA2023-146"    "RICA2023-149"    "RICA2023-150"    "RICA2023-153"   
    ## [33] "RICA2023-154"    "RICA2023-155"    "RICA2023-156"    "RICA2023-157"   
    ## [37] "RICA2023-160"    "RICA2023-18"     "RICA2023-19"     "RICA2023-2"     
    ## [41] "RICA2023-21"     "RICA2023-22"     "RICA2023-23"     "RICA2023-24"    
    ## [45] "RICA2023-29"     "RICA2023-3"      "RICA2023-30"     "RICA2023-31"    
    ## [49] "RICA2023-32"     "RICA2023-37"     "RICA2023-38"     "RICA2023-39"    
    ## [53] "RICA2023-40"     "RICA2023-42"     "RICA2023-45"     "RICA2023-46"    
    ## [57] "RICA2023-47"     "RICA2023-5"      "RICA2023-50"     "RICA2023-51"    
    ## [61] "RICA2023-55"     "RICA2023-56"     "RICA2023-57"     "RICA2023-58"    
    ## [65] "RICA2023-59"     "RICA2023-6"      "RICA2023-60"     "RICA2023-61"    
    ## [69] "RICA2023-63"     "RICA2023-65"     "RICA2023-68"     "RICA2023-69"    
    ## [73] "RICA2023-7"      "RICA2023-70"     "RICA2023-71"     "RICA2023-72"    
    ## [77] "RICA2023-73"     "RICA2023-8"      "RICA2023-80"     "RICA2023-82"    
    ## [81] "RICA2023-83"     "RICA2023-84"     "RICA2023-85"     "RICA2023-86"    
    ## [85] "RICA2023-88"     "RICA2023-89"     "RICA2023-9"      "RICA2023-93"    
    ## [89] "RICA2023-94"     "RICA2023-95"     "RICA2023-98"     "RICA2023-NTCP2" 
    ## [93] "RICA2023-NTCP3"  "RICA2023-NTCPCR"

``` r
head(taxa_names(tax_table))
```

    ## [1] "0000f9a9a96f81bcfc9ab4393bf33e0e" "006b0e6060f45ba91f316426770dfa6c"
    ## [3] "006b922029778aad6b92f6d17e0ae056" "009cc48534879dcb5c792b9410603ba2"
    ## [5] "009feb0fdbe776e597506f303b0a8d33" "0118369d3233378cc756440e3dd46844"

``` r
head(taxa_names(otu))
```

    ## [1] "bd53e2483c985a090089d93fc039c4ac" "70d69f68d0230abfbb2f096379d6131e"
    ## [3] "d4a521b98cbc2f7183b70400d38c0216" "43c26c4de9813879819ff2893e668e67"
    ## [5] "3aa965194c47dee2ebb9ad04c942930a" "f7078305113a3daacae2e00137d0a262"

``` r
#lets make that phyloseq object now
phy= phyloseq(otu, tax_table, sample)

#remove Mitochondria/ chloroplast
phy <- phy %>% subset_taxa( Family!= "Mitochondria" | is.na(Family) & Class!="Chloroplast" | is.na(Class) ) 

#remove NTC and a sample with 20 total reads 
phy <-phy %>% subset_samples(Species!="NA")
phy <-phy %>% subset_samples(Bag_Number!="32")

#filter out the noise in the data, this filters out any taxa that have counts less than 5
phy_fill= filter_taxa(phy , function(x) sum(x > 5) > (0.05*length(x)), TRUE)

## Basic summary information for the phyloseq object that you just created
#number of taxa
ntaxa(phy_fill)
```

    ## [1] 138

``` r
#number of samples
nsamples(phy_fill)
```

    ## [1] 90

``` r
#number of variables
sample_variables(phy_fill)
```

    ##  [1] "qPCR_Plate_ID"                       "CAM_4.CT.mean"                      
    ##  [3] "TLC_8.CT.mean"                       "CAM_4.CT.sd"                        
    ##  [5] "TLC_8.CT.sd"                         "CAM_4.reps"                         
    ##  [7] "TLC_8.reps"                          "TLC_CAM_ratio"                      
    ##  [9] "ratio_mean"                          "delta_ratio"                        
    ## [11] "Genotype_Timepoint"                  "Extraction_Well"                    
    ## [13] "Extraction_Date"                     "Extraction_plate"                   
    ## [15] "Extraction_Kit"                      "Extraction_Kit_CatNumber"           
    ## [17] "Collection_Date"                     "Location"                           
    ## [19] "Species"                             "Bag_Number"                         
    ## [21] "Genotype"                            "Source_Lat"                         
    ## [23] "Source_Long"                         "Tree"                               
    ## [25] "Time_Point"                          "DHW"                                
    ## [27] "Order_Sampled"                       "Orgin_Nursery"                      
    ## [29] "Nursery_Lat"                         "Nursery_Long"                       
    ## [31] "conc_ngul"                           "X260.280"                           
    ## [33] "X260.230"                            "ddPCR_RICA_count"                   
    ## [35] "Bleach_index"                        "Bleach_index_rank"                  
    ## [37] "Growth_Specific"                     "Growth_index"                       
    ## [39] "Growth_6mon"                         "Growth_6mon_rank"                   
    ## [41] "symbiont_density_hemocytometer"      "symbiont_density_hemocytometer_rank"
    ## [43] "Heat_tolerance_CBASS"                "Heat_tolerance_CBASS_rank"          
    ## [45] "Notes"

``` r
#Use tax_fix() on your phyloseq data with default arguments to repair most tax_table problems (missing or uninformative values)
phy_fill <- tax_fix(phy_fill)
```

``` r
#####standard PCA for t1 
phy_fill %>% 
  subset_samples(Time_Point=="T1") %>%
  tax_transform("identity", rank = "Genus") %>% # don't transform!
  ord_calc("PCA") %>%
  ord_get() %>%
  phyloseq::plot_scree() + theme(axis.text.x = element_text(size = 10))
```

![](MicroViz_PCOA_MARKDOWN_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
#ordiations plot with ord_calc for standard PCA
phy_fill %>%
  subset_samples(Time_Point=="T1") %>%
  tax_transform("clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>%
  ord_plot(axes = c(1, 2),color = "Genotype", plot_taxa = 1:3, size = 2.5,
           tax_lab_style = tax_lab_style(size = 3, alpha = 0.5))+
  scale_colour_manual(values=c("#56B4E9","#CBD588","#5F7FC7", "darkorchid","#DA5724","#CD9BCD","#006700","yellow2",
                               "gray80", "#AD6F3B", "#673770","#8569D5", 
                               "#5E738F","#D1A33D", "orange","#ff6700","aquamarine4", "#652926",
                               "lightblue4", "lightpink","royalblue4","#D14285",
                               "palevioletred1", "#56B4E9","#CBD588", "#5F7FC7","#DA5724",
                               "#CD9BCD", "gray80", "darkorchid",
                               "#AD6F3B", "#673770","#D14285", "#652926"))
```

![](MicroViz_PCOA_MARKDOWN_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
#####standard PCA for t1 to t2 with arrows for taxa driving the seperation 

#look at distibution of the eigenvalues, for this data most was in PC1-PC2, so I am selecting those by using axes = c(1, 2) when plotting below
phy_fill %>%
  tax_transform("identity", rank = "Genus") %>% # don't transform!
  ord_calc("PCA") %>%
  ord_get() %>%
  phyloseq::plot_scree() + theme(axis.text.x = element_text(size = 10))
```

![](MicroViz_PCOA_MARKDOWN_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
#ordiations plot with ord_calc for standard PCA
phy_fill %>%
  tax_transform("clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>%
  ord_plot(axes = c(1, 2),color = "Time_Point", plot_taxa = 1:3, size = 2.5,
           tax_lab_style = tax_lab_style(size = 3.4, alpha = 0.5))+
  scale_colour_manual(values=c("blue", "red"))
```

![](MicroViz_PCOA_MARKDOWN_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

``` r
#test the if groups are different 
phy_clr = microbiome::transform(phy_fill, 'clr')
pairwise.adonis(t(otu_table(phy_clr)), sample_data(phy_clr)$Time_Point, sim.method = "euclidean",
                p.adjust.m = "bonferroni")
```

    ##      pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
    ## 1 T2 vs T1  1  2135.681 4.483397 0.04847786   0.001      0.001  **

``` r
#PCOA with time point as shape, genotype as color

phy_fill %>%
  tax_transform("clr", rank = "Genus") %>%
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc(method = "PCA") %>%
  ord_plot(axes = c(1, 2),color = "Genotype", shape="Time_Point", plot_taxa = 1:3, size = 2.9,
           tax_lab_style = tax_lab_style(size = 3.4, alpha = 0.5))+
  scale_colour_manual(values=c("#56B4E9","#CBD588","#5F7FC7", "darkorchid","#DA5724","#CD9BCD","#006700","yellow2",
                               "gray80", "#AD6F3B", "#673770","#8569D5", 
                               "#5E738F","#D1A33D", "orange","#ff6700","aquamarine4", "#652926",
                               "lightblue4", "lightpink","royalblue4","#D14285",
                               "palevioletred1", "#56B4E9","#CBD588", "#5F7FC7","#DA5724",
                               "#CD9BCD", "gray80", "darkorchid",
                               "#AD6F3B", "#673770","#D14285", "#652926"))
```

![](MicroViz_PCOA_MARKDOWN_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
######Aitchison Distance 
#look at distibution of the eigenvalues
phy_fill %>%
  tax_transform("identity", rank = "Genus") %>% # don't transform!
  dist_calc("aitchison") %>%
  ord_calc("PCoA") %>%
  ord_get() %>%
  phyloseq::plot_scree() + theme(axis.text.x = element_text(size = 6))
```

![](MicroViz_PCOA_MARKDOWN_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
#ordiations plot with ord_calc, Aitchison Distance 
phy_fill %>%
    tax_transform("identity", rank = "Genus") %>% # don't transform!
    dist_calc("aitchison") %>%
    ord_calc(method = "PCoA") %>%
    ord_plot(color  = "Time_Point", plot_taxa = 1:5, size = 2.9)
```

![](MicroViz_PCOA_MARKDOWN_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->
