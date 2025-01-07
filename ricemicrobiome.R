library(phyloseq)
library(microeco)
library(file2meco)
library(tidyr)
library(picante)
library(FSA)
library(rcompanion)
library(agricolae)
library(ggpubr)
library(WGCNA)
library(rgexf)
library(circlize)
library(ALDEx2)

##importing the qiime2 file and converting into the microeco package

ps <- qiime2meco(feature_table = "Desktop/Plant Microbiome/plant_denoise_table_for.qza",taxonomy_table = "Desktop/Plant Microbiome/SE-taxonomic_classification_gg.qza",sample_table = "Downloads/Sima-paper_1-16s rRNA-sequences/Sample_Metadata.tsv",phylo_tree = "Desktop/Plant Microbiome/rooted-tree.qza")
ps

## removing chloroplast and mitochondria reads 

ps$filter_pollution(taxa = c("mitochondria", "chloroplast"))

##make the samples consistent 
ps$tidy_dataset()

print(ps)

##rarefy the samples 

ps$sample_sums() %>% range

##as the minimum value is 64635, the minimum value is 60000

ps$rarefy_samples(sample.size = 60000)

##re-check the range

ps$sample_sums() %>% range

##caluclate the taxa_abundance

ps$cal_abund()

##Phylum-taxonomy-barplot

p1 <- trans_abund$new(dataset = ps, taxrank = "Phylum", ntaxa = 6)

p1$plot_bar(others_color = "grey70", facet = "Group_name", xtext_keep = FALSE, legend_text_italic = FALSE)

##Phylum-level-piechart

p1.1 <- trans_abund$new(dataset = ps, taxrank = "Phylum", ntaxa = 6, groupmean = "Group_name")

p1.1$plot_pie(facet_nrow = 1, add_label = TRUE)

##Class-level-boxplot

c1 <- trans_abund$new(dataset = ps, taxrank = "Class", ntaxa = 10)

c1$plot_box(group = "Group_name", xtext_angle = 30)

##Genus-level-heatmap

g1 <- trans_abund$new(dataset = ps, taxrank = "Genus", ntaxa = 20)

g1$plot_heatmap(facet = "Group_name", xtextkeep = FALSE, withmargin = FALSE)

##venn-diagram 

ps1 <- ps$merge_samples(use_group = "Group_name")

v1 <- trans_venn$new(ps1, ratio = NULL)

v1$plot_venn()

##alpha_diversity

ps$cal_alphadiv(PD = TRUE)

ad1 <- trans_alpha$new(dataset = ps, group = "Group_name")

head(ad1$data_stat)

ad1$cal_diff(method = "KW")

head(ad1$res_diff)

ad1$cal_diff(method = "KW_dunn")

head(ad1$res_diff)

ad1$cal_diff(method = "anova")

head(ad1$res_diff)

ad1$cal_diff(method = "wilcox")

ad1$plot_alpha(measure = "Chao1", shape = "Group_name", y_start = 0.1, y_increase = 0.1, add_sig_text_size = 6, boxplot_add = "jitter", order_x_mean = TRUE)

ad1$plot_alpha(measure = "Shannon", shape = "Group_name", y_start = 0.1, y_increase = 0.1, add_sig_text_size = 6, boxplot_add = "jitter", order_x_mean = TRUE)

##beta-diversity 

ps$cal_betadiv(unifrac = TRUE)

bd1 <- trans_beta$new(dataset = ps, group = "Group_name", measure = "bray")

bd1$cal_ordination(ordination = "PCoA")

bd1$plot_ordination(plot_color = "Group_name", plot_shape = "Group_name", plot_type = c("point", "ellipse"))

bd1$cal_manova(manova_all = TRUE)

bd1$res_manova

#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999

#adonis2(formula = use_formula, data = metadata)
#Df SumOfSqs      R2      F Pr(>F)    
#Group_name  3   1.1027 0.49618 2.6262  0.001 ***
#  Residual    8   1.1197 0.50382                  
#Total      11   2.2224 1.00000                  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

bd1$cal_manova(manova_all = FALSE)

bd1$res_manova

#Groups measure        F        R2 p.value p.adjusted Significance
#1   Karuppu Kuruvai vs Kaatu Yaanam    bray 3.257997 0.4488838     0.1        0.1             
#2 Karuppu Kuruvai vs Karuppu Kavuni    bray 2.359639 0.3710335     0.1        0.1             
#3     Karuppu Kuruvai vs Kala Namak    bray 2.407392 0.3757210     0.1        0.1             
#4    Kaatu Yaanam vs Karuppu Kavuni    bray 2.723105 0.4050368     0.1        0.1             
#5        Kaatu Yaanam vs Kala Namak    bray 3.355200 0.4561671     0.1        0.1             
#6      Karuppu Kavuni vs Kala Namak    bray 2.077648 0.3418506     0.1        0.1   

bd1$cal_ordination(ordination = "NMDS")

bd1$plot_ordination(plot_color = "Group_name", plot_shape = "Group_name", plot_type = c("point", "ellipse"))


bd2 <- trans_beta$new(dataset = ps, group = "Group_name", measure = "wei_unifrac")

bd2$cal_ordination(ordination = "PCoA")

bd2$plot_ordination(plot_color = "Group_name", plot_shape = "Group_name", plot_type = c("point", "ellipse"))

bd2$cal_ordination(ordination = "NMDS")

bd2$plot_ordination(plot_color = "Group_name", plot_shape = "Group_name", plot_type = c("point", "ellipse"))

bd2$cal_manova(manova_all = TRUE)

bd2$res_manova

#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999

#adonis2(formula = use_formula, data = metadata)
#Df SumOfSqs      R2      F Pr(>F)    
#Group_name  3 0.033107 0.53566 3.0763  0.001 ***
#  Residual    8 0.028698 0.46434                  
#Total      11 0.061805 1.00000                  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

bd2$cal_manova(manova_all = FALSE)

bd2$res_manova

#Groups     measure        F        R2 p.value p.adjusted Significance
#1   Karuppu Kuruvai vs Kaatu Yaanam wei_unifrac 4.668342 0.5385507     0.1       0.12             
#2 Karuppu Kuruvai vs Karuppu Kavuni wei_unifrac 2.582118 0.3922929     0.1       0.12             
#3     Karuppu Kuruvai vs Kala Namak wei_unifrac 2.086760 0.3428359     0.2       0.20             
#4    Kaatu Yaanam vs Karuppu Kavuni wei_unifrac 3.633543 0.4759970     0.1       0.12             
#5        Kaatu Yaanam vs Kala Namak wei_unifrac 3.544686 0.4698255     0.1       0.12             
#6      Karuppu Kavuni vs Kala Namak wei_unifrac 2.947998 0.4242946     0.1       0.12     

bd3 <- trans_beta$new(dataset = ps, group = "Group_name", measure = "unwei_unifrac")

bd3$cal_ordination(ordination = "PCoA")

bd3$plot_ordination(plot_color = "Group_name", plot_shape = "Group_name", plot_type = c("point", "ellipse"))

bd3$cal_ordination(ordination = "NMDS")

bd3$plot_ordination(plot_color = "Group_name", plot_shape = "Group_name", plot_type = c("point", "ellipse"))

bd3$cal_manova(manova_all = TRUE)

bd3$res_manova

#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999

#adonis2(formula = use_formula, data = metadata)
#Df SumOfSqs      R2      F Pr(>F)    
#Group_name  3  0.66314 0.40015 1.7789  0.001 ***
#  Residual    8  0.99410 0.59985                  
#Total      11  1.65725 1.00000                  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

bd3$cal_manova(manova_all = FALSE)

bd3$res_manova

#Groups       measure        F        R2 p.value p.adjusted Significance
#1   Karuppu Kuruvai vs Kaatu Yaanam unwei_unifrac 2.081707 0.3422899     0.1        0.1             
#2 Karuppu Kuruvai vs Karuppu Kavuni unwei_unifrac 1.893023 0.3212312     0.1        0.1             
#3     Karuppu Kuruvai vs Kala Namak unwei_unifrac 1.448751 0.2658868     0.1        0.1             
#4    Kaatu Yaanam vs Karuppu Kavuni unwei_unifrac 1.759201 0.3054592     0.1        0.1             
#5        Kaatu Yaanam vs Kala Namak unwei_unifrac 1.925843 0.3249906     0.1        0.1             
#6      Karuppu Kavuni vs Kala Namak unwei_unifrac 1.640058 0.2907874     0.1        0.1      

#network plot

t1 <- trans_network$new(dataset = ps, cor_method = "spearman", use_WGCNA_pearson_spearman = TRUE, filter_thres = 0.0001)

t1$cal_network(COR_p_thres = 0.01, COR_optimization = TRUE)

t1$cal_network(COR_p_thres = 0.01, COR_cut = 0.63)

t1$res_network

# invoke igraph cluster_fast_greedy function for this undirected network 
t1$cal_module(method = "cluster_fast_greedy")

t1$cal_network_attr()

t1$res_network_attr

#Vertex                 1.025000e+03
#Edge                   2.868000e+03
#Average_degree         5.596098e+00
#Average_path_length    2.181420e+00
#Network_diameter       1.300000e+01
#Clustering_coefficient 9.857816e-01
#Density                5.464939e-03
#Heterogeneity          1.371098e+00
#Centralization         2.773819e-02
#Modularity             9.276083e-01

t1$cal_sum_links(taxa_level = "Phylum")

t1$plot_sum_links(method = "circlize", transparency = 0.2, annotationTrackHeight = circlize::mm_h(c(2, 5)))

#for interactive_plot
t1$plot_sum_links(plot_pos = TRUE, plot_num = 10, color_values = RColorBrewer::brewer.pal(10, "Paired"))


#### Differential Abundance

df1 <- trans_diff$new(dataset = ps, method = "ALDEx2_kw", group = "Group_name", taxa_level = "Genus")

df1$plot_diff_abund(use_number = 1:10, group_order = c("Karuppu Kavuni", "Kaatu Yaanam", "Karuppu Kuruvai", "Kala Namak"), add_sig = TRUE)


