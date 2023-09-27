

library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
library(metagMisc)
library(ggpubr)
library(plyr)
library(miaTime)
library(knitr)
library(mia)

########## load and prepare ####

# load phyloseq object with adjsuted and emrged ASV names (fixing macrogen's inadequate bioinformatics)
load(file = "./Data/physeq_merged_renamed_pot_field_ASVs.RData")
physeq


#remove rare ASVs
otu_table(physeq) <- otu_table(physeq)[which(rowSums(otu_table(physeq)) > 7), ] 

#check library sizes
sample_sums(physeq)%>%hist
sample_sums(physeq)%>%sort(decreasing = FALSE)

#normalize data
physeq_filtered_df <- as.data.frame(otu_table(physeq))
rarefaction_curve<-
rarecurve(t(physeq_filtered_df),
          label = FALSE,
          step = 1000,
          main = "Rarefaction ", ylab = "Number of ASVs", xlab = "Number of DNA Sequences",
          abline(v = min(sample_sums(physeq)), lwd = 3, lty = 2))

ggsave(rarefaction_curve, filename = "./Results/Rarefaction_curve.pdf")


# now, let's rarefy the data
set.seed(100) 
physeq_rarefied <- rarefy_even_depth(
                                      physeq = physeq,
                                      sample.size = min(sample_sums(physeq)),
                                      rngseed = FALSE,
                                      replace = TRUE,
                                      trimOTUs = TRUE,
                                      verbose = TRUE)



#adjsut order of treatment factors, pallete
physeq_rarefied@sam_data$Treatment<-
factor(physeq_rarefied@sam_data$Treatment, levels = c("Control", "Manure", "Chitin", "BSF", "HC", "MW"))
pallete_els<-c("#000000", "#767171", "#AFABAB","#2F5597","#FFC000","#548235")

#separate phyloseq objects
physeq_rarefied_l<-phyloseq_sep_variable(physeq_rarefied, variable = c("Soil_type", "Time_point"))
physeq_rarefied_soiltype_l<-phyloseq_sep_variable(physeq_rarefied, variable = c("Soil_type"))

#save sepaated objects externaly
save(physeq_rarefied_l, file = "./Data/physeq_rarefied_l.RData")
save(physeq_rarefied_soiltype_l, file = "./Data/physeq_rarefied_soiltype_l.RData")








########## beta diversity plots ####

# Let's make a Non-Metric Multidimensional Scaling (NMDS) of all our samples 
nmds_rootAndSoil <- phyloseq::ordinate(physeq_rarefied,
                                       method = "NMDS", # this method has few assumptions and readly accepts different data structures
                                       distance = "bray", # bray-curtis distance is suitable for sparse data - such as the zero-inflated microbiome data we have
                                       trymax=200,
                                       try = 200,
                                       autotransform = TRUE) 

# plot the NMDS
NMDS_all_samples<-plot_ordination(
  physeq = physeq_rarefied,
  ordination = nmds_rootAndSoil,
  color = "Treatment",
  shape = "Soil_type"
) +
  theme_classic() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(aes(size = Time_point), alpha = 1) +
  theme(legend.position = "right")+
  scale_color_manual(values=pallete_els)

ggsave(NMDS_all_samples, filename = "./Results/NMDS_all_samples.pdf",
       width = 180,
       height = 180,
       units = "mm")

#plot an alternative NMDS visualization
NMDS_all_samples2<-plot_ordination(
  physeq = physeq_rarefied,
  ordination = nmds_rootAndSoil,
  color = "Soil_type",
  shape = "Treatment"
) +
  theme_classic() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(aes(size = Time_point), alpha = 1) +
  theme(legend.position = "right")+
  scale_color_manual(values=pallete_els)

ggsave(NMDS_all_samples2, filename = "./Results/NMDS_all_samples2.pdf",
       width = 180,
       height = 180,
       units = "mm")

#plot an alternative NMDS visualization
NMDS_all_samples3<-plot_ordination(
  physeq = physeq_rarefied,
  ordination = nmds_rootAndSoil,
  color = "Time_point",
  shape = "Treatment"
) +
  theme_classic() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  theme(legend.position = "right")+
  facet_wrap(~Soil_type)+
  scale_color_manual(values=pallete_els)

ggsave(NMDS_all_samples3, filename = "./Results/NMDS_all_samples3.pdf",
       width = 180,
       height = 180,
       units = "mm")



# define custom function to plot a list of NMDS results
NMDS_listing <- function(physeq_list) { # first the name of the new function you are saving, then you say you are going to define a function, and then you set it's arguments - here, a single list of phyloseq objects. the {} indicate where the contents of the custom function start and end.
  
  # here we perform the NMDS on all elements of the list using phyloseq::ordinate()
  NMDS_list <- lapply(physeq_list, function(x) {
    ordinate(
      physeq = x, # your phyloseq object, each element x in the list
      method = "NMDS", # ordination method to use
      distance = "bray", # distance metric to use
      try = 200, # number of permutations
      autotransform = TRUE
    )
  }) # automatic data transformation
  
  # mapply will perform one function on two lists. define the lists in the end of the function
  untitled_plot_list <- mapply(function(x, y) { # mapply will run one function on 2 lists
    plot_ordination(
      physeq = x, # your phyloseq object, each element x in the list
      ordination = y, # your phyloseq-generated ordination object, each element y in the list made above
      color = "Treatment"
    ) + # dot colors in the plot
      theme_classic() + # ggplot option to remove background
      labs(subtitle = paste("Stress:", round(y$stress, digits = 4))) + # this adds the NMDS stress to the plot as subtitle
      theme(plot.title = element_text(size = 10, face = "bold")) + # options for the title
      geom_point(aes(size = Time_point), alpha = 1) +
      theme(legend.position = "right")+
      scale_color_manual(values=pallete_els)
  }, # position for the legend
  x = physeq_list, # note that you only define x here, after you define your plot function
  y = NMDS_list, # note that you only definee y here, after you define your plot function
  SIMPLIFY = FALSE
  ) # if you simply, all resusts are saved in a vector
  
  # The plots above miss a title, so you don't really know who is root and who is soil. this code create a list of names to be used as titles in your plot
  tiles_list <- names(untitled_plot_list)
  
  Plot_list <- mapply(function(z, w) { # use mapply again to put those listed names in the list of plots
    z + ggtitle(w)
  },
  z = untitled_plot_list,
  w = tiles_list,
  SIMPLIFY = FALSE
  )
  
  return(list(NMDS_list, Plot_list)) # with this you specify the output of your custom function. Note we are saving both the NMDS calculations and the plot as a list, so we don't ahve to re-calculate the NMDS later
}


# run NMDS and put it in a single plot
NMDS_plot<-NMDS_listing(physeq_list = physeq_rarefied_l)

NMDS_plot_grid<-ggarrange(plotlist = NMDS_plot[[2]],ncol = 4, nrow = 2, common.legend = TRUE)  

ggsave(NMDS_plot_grid, filename = "./Results/NMDS_plot_grid.pdf",
       width = 400,
       height = 230,
       units = "mm")









########## beta dispersion ####

#run beta dispersion
#define funciton
beta_disp_plotAndTest <- function(phyloseq_list, group) {
  # phyloseq_list = a list of phyloseq objects
  # group = the variale you want to test the beta dispersion of, in quotes
  beta_disper_list <- lapply(phyloseq_list, function(x) {
    betadisper(phyloseq::distance(t(otu_table(x)), method = "bray"), sample_data(x)[[group]])
  }) # selects only column "group""
  
  # gets the names of the list
  tiles_list <- names(beta_disper_list)
  
  # runs anova on beta dispersions
  get_p <- lapply(beta_disper_list, function(x) {
    anova(x, permutations = 999)
  })
  
  
  p_dispersion <- map(get_p, 5) # gets the p value of the dispersion test
  p_dispersion <- p_dispersion[!is.na(p_dispersion)] # removes the generated NA
  
  
  # runs anova on beta dispersions
  bet_disp_PCOa_plot <- mapply(function(x, y, z) {
    plot(x,
         main = y,
         sub = z,
         xlab = "p value for homogeniety test:",
         ylab = "PCoA"
    )
  },
  x = beta_disper_list,
  y = tiles_list,
  z = p_dispersion,
  SIMPLIFY = FALSE
  )
  
  # runs anova on beta dispersions
  bet_disp_boxplot <- mapply(function(x, y) {
    boxplot(x, main = y)
  },
  x = beta_disper_list,
  y = tiles_list,
  SIMPLIFY = FALSE
  )
  
  
  return(bet_disp_PCOa_plot)
}

# now that we have a custom function, we can run it across all lists and variables
library(purrr)
set.seed(5235)
dip_result3 <- beta_disp_plotAndTest(physeq_rarefied_l, "Treatment")












########## beta diversity tests ####


# let's check now if sample type differ per species. note that this function uses Blocks as strata
permanova_with_blocks <- function(phyloseq_list, rhs_model) {
  # phyloseq_list = list of phyloseq objects
  # RHS_model = right hand side model, such as MeJA_treatment*Sample_type + Block
  lapply(phyloseq_list, function(x) {
    set.seed(1234)
    lhs_data <- phyloseq::distance(t(otu_table(x)), method = "bray")
    rhs_model_in <- paste(rhs_model)
    form <- as.formula(paste("lhs_data~", paste(rhs_model_in))) # getting the formulat properly evaluated as custom string is tricky
    output <- adonis2(form,
                      data = as(sample_data(x), "data.frame"), # changing with as.data.frame is insufficient
                      permutations = how(within = Within(type = "free"), nperm = 9999)
    ) # how defines the permutations, it is important to adjust it to the experimental design such as a time series
    return(output)
  })
}

permanova_tables<-permanova_with_blocks(phyloseq_list = physeq_rarefied_l,rhs_model = "Treatment")

# save and export as flextable
library("flextable")
library(tibble)
permanova_df_l<-lapply(permanova_tables, function(x){
  
  #save as dataframe, adjust rownames
  table_df<-as.data.frame(x)
  table_df<-rownames_to_column(table_df, var = "Factor")
  
  #round values, reduce numebr of digits
  table_df$SumOfSqs<-round(table_df$SumOfSqs, digits = 3)
  table_df$R2<-round(table_df$R2, digits = 3)
  table_df$F<-round(table_df$F, digits = 3)

  
  return(table_df)
  
  
})


flextable_l<-lapply(permanova_df_l, flextable)


names(flextable_l)<-names(permanova_df_l)

save_as_docx(flextable_l$Field.Week.4,
             flextable_l$Field.Week.8,
             flextable_l$Field.Week.12,
             flextable_l$Field.Week.16,
             flextable_l$Pot.Week.4,
             flextable_l$Pot.Week.8,
             flextable_l$Pot.Week.12,
             flextable_l$Pot.Week.16,
             values = flextable_l,
             path = "./Results/PERMANOVA_microbiome_tables_8plots.docx")















######### pairwise permanova
library("EcolUtils")
set.seed(303848)
pairwise_permanova <- lapply(physeq_rarefied_l, function (x)
    adonis.pair(dist.mat= phyloseq::distance(otu_table(x), method="bray"),
                Factor= as.factor(as(phyloseq::sample_data(x),"data.frame")$Treatment)))



# set time and site varaibles for mapply
site<-as.list(c(rep("Field", 4), rep("Pot", 4)))
# time<-as.list(rep(c("Week 4", "Week 8", "Week 12", "Week 16") , 2))
 time<-as.list(rep(c(4, 8, 12, 16) , 2))
# add site and time across mapply
pairwise_permanova_df<-
    mapply(function (x,y,z){
    x$Site<-z
    x$Time<-y
    return(x)
    },
  x = pairwise_permanova,
  z = site,
  y = time,
  SIMPLIFY = FALSE)

  # put all elements of the list into a single DF
pairwise_permanova_df<-do.call("rbind", pairwise_permanova_df)


library(stringr)
split_combiation<-str_split_fixed(string =pairwise_permanova_df$combination, pattern = " <-> ", n=2 )

pairwise_permanova_df$treatment_1<-split_combiation[,1]
pairwise_permanova_df$treatment_2<-split_combiation[,2]

pairwise_permanova_df$SumsOfSqs<-round(pairwise_permanova_df$SumsOfSqs, digits = 3)
pairwise_permanova_df$MeanSqs<-round(pairwise_permanova_df$MeanSqs, digits = 3)
pairwise_permanova_df$F.Model<-round(pairwise_permanova_df$F.Model, digits = 3)
pairwise_permanova_df$R2<-round(pairwise_permanova_df$R2, digits = 3)
pairwise_permanova_df$P.value<-round(pairwise_permanova_df$P.value, digits = 3)
pairwise_permanova_df$P.value.corrected<-round(pairwise_permanova_df$P.value.corrected, digits = 3)

pairwise_ft<-flextable(pairwise_permanova_df)

#flextable of pairwside comparsions
save_as_docx(pairwise_ft,
             path = "./Results/pairwise_PERMANOVA_tables.docx")


# filter df to make many plots
Control_df<-filter(pairwise_permanova_df, treatment_1 =="Control" | treatment_2 == "Control" )
Chitin_df<-filter(pairwise_permanova_df, treatment_1 =="Chitin" | treatment_2 == "Chitin" )
Manure_df<-filter(pairwise_permanova_df, treatment_1 =="Manure" | treatment_2 == "Manure" )
BSF_df<-filter(pairwise_permanova_df, treatment_1 =="BSF" | treatment_2 == "BSF" )
MW_df<-filter(pairwise_permanova_df, treatment_1 =="MW" | treatment_2 == "MW" )
HC_df<-filter(pairwise_permanova_df, treatment_1 =="HC" | treatment_2 == "HC" )


# line plot with one focal comparson at a time 
focal_line_plot<-
ggarrange(
  ggplot(Control_df, aes(x = Time ,y = R2))+
    geom_line(aes(color = combination))+
    facet_wrap(~Site)+
    theme_bw(),
  ggplot(Chitin_df, aes(x = Time ,y = R2))+
    geom_line(aes(color = combination))+
    facet_wrap(~Site)+
    theme_bw(),
  ggplot(Manure_df, aes(x = Time ,y = R2))+
    geom_line(aes(color = combination))+
    facet_wrap(~Site)+
    theme_bw(),
  ggplot(BSF_df, aes(x = Time ,y = R2))+
    geom_line(aes(color = combination))+
    facet_wrap(~Site)+
    theme_bw(),
  ggplot(MW_df, aes(x = Time ,y = R2))+
    geom_line(aes(color = combination))+
    facet_wrap(~Site)+
    theme_bw(),
  ggplot(HC_df, aes(x = Time ,y = R2))+
    geom_line(aes(color = combination))+
    facet_wrap(~Site)+
    theme_bw(), ncol = 1)


ggsave(focal_line_plot, filename = "./Results/comparison_lineplot_pairwise_permanova_microbiome.pdf",
       width = 120,
       height = 300,
       units = "mm")















##### permanova of mirobiome with treatment*time in pot and field samples in separate

permanova_tables2<-permanova_with_blocks(phyloseq_list = physeq_rarefied_soiltype_l ,rhs_model = "Treatment*Time_point")



permanova2_df_l<-lapply(permanova_tables2, function(x){
  
  #save as dataframe, adjust rownames
  table_df<-as.data.frame(x)
  table_df<-rownames_to_column(table_df, var = "Factor")
  
  #round values, reduce numebr of digits
  table_df$SumOfSqs<-round(table_df$SumOfSqs, digits = 3)
  table_df$R2<-round(table_df$R2, digits = 3)
  table_df$F<-round(table_df$F, digits = 3)
  
  
  return(table_df)
  
  
})


flextable_l2<-lapply(permanova2_df_l, flextable)


names(flextable_l2)<-names(permanova2_df_l)

save_as_docx(flextable_l2$Field,
             flextable_l2$Pot,
             values = flextable_l2,
             path = "./Results/PERMANOVA_microbiome_tables_2plots.docx")












########## alpha diversity plot ####

physeq_rarefied
library(microbiome)
# calculate diversity in topsoils

  taxon_diversity<-microbiome::diversity(physeq_rarefied,index = "Shannon")
  merg_to_ps<-sample_data(taxon_diversity) 
  taxon_diversity<-as(sample_data(merge_phyloseq(physeq_rarefied,merg_to_ps)),"data.frame") 


# plot diversity from topsoils
shannon_topsoil_plot<-
  ggplot(taxon_diversity, aes(x =Time_point , y = shannon, fill=Treatment))+
    geom_boxplot(alpha = 0.4, outlier.shape = NA)+
    geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(width = .75), dotsize = 0.3)+
    theme_bw()+
    scale_fill_manual(values=pallete_els)+
    labs(title = "Shannon diversity")+
    ylab("Shannon diversity index")+
    xlab("Time points")+
    facet_wrap(~Soil_type)
#save shannon plot
ggsave(shannon_topsoil_plot, filename = "./Results/shannon_plot.pdf",
       width = 180,
       height = 180,
       units = "mm")


########## alpha diversity test ####


# check homogeniety of variances in alpha diversity
library(car)
leveneTest((shannon) ~ Treatment*Soil_type*Time_point, data = taxon_diversity) 

#two-way anova, shanon diversity index
  tx <- with(taxon_diversity, interaction(Treatment,Soil_type,Time_point)) #needed for tukey to test interaction
  aovTukey<-aov(shannon ~ tx, data = taxon_diversity)#needed for tukey to test interaction
  anova_test<-Anova(lm((shannon) ~Treatment*Soil_type*Time_point, data = taxon_diversity, contrasts=list(Treatment=contr.sum, Soil_type=contr.sum, Time_point = contr.sum)), type = "2") 
  tukey_result<-TukeyHSD(aovTukey, "tx", group=TRUE, console=TRUE)#post-hoc
  

  
  anova_test$`Sum Sq`<-round(anova_test$`Sum Sq`, digits = 3)
  anova_test$`F value`<-round(anova_test$`F value`, digits = 3)
  anova_test$`Pr(>F)`<-round(anova_test$`Pr(>F)`, digits = 3)
  anova_test_df<-as.data.frame(anova_test)
  anova_test_df<- rownames_to_column(anova_test_df, var = "Factor")
  anova_test_ft<-flextable(anova_test_df)
  
  
  #flextable of pairwside comparsions
  save_as_docx(anova_test_ft,
               path = "./Results/shannon_anova.docx")
  
  
  
  
  
  
  
  #### matrix of heat trees ####
  
  library(metacoder)
  phyloseq_to_heat_tree_matrix2<-function(ps_object, sample_group){
    
    
    
    
    # this function output is a heat trees comparing metadata. the input is based on a phyloseq object
    
    # ps object =  a phyloseq object, containing sample metadata, otu table, and taxonomy table
    # sample_group = the name of a column in your emtadata that you want to compare in the heat tree. it has to be quoted.
    
    # this function will create a matrix of heat trees if your metadata has more than 2 groups. it should fail if it has only 1 group
    
    # turn NA into "" so it is not completety missing ofmr the tree  
    ps_object@tax_table@.Data[is.na(ps_object@tax_table)]<-""
    
    #remove unecessary taxonomic info (dada2id, "s__" and "ASV_id) by updating the tax table with a subset of the tax table
    tax_table(ps_object)<-tax_table(ps_object)[,1:6]
    
    
    # let's remove the "r__"ranks from the taxonomy, they can be useful but will pollute our plot
    tax_table(ps_object)[, colnames(tax_table(ps_object))] <- gsub(pattern = "[a-z]__", # regular expression pattern to search for
                                                                   x = tax_table(ps_object)[, colnames(tax_table(ps_object))], # "df"
                                                                   replacement = "") # replacement for pattern
    # transform from phyloseq to taxmap object
    taxmap_obj<-parse_phyloseq(ps_object)
    
    #get abundance per taxon
    taxmap_obj$data$tax_abund<-calc_taxon_abund(obj = taxmap_obj, 
                                                data = "otu_table",
                                                cols = taxmap_obj$data$sample_data$sample_id) 
    
    #get occurrence of ASVs per treatment
    # the sample groups needs some wrangling to fit within the soft code of the function
    taxmap_obj$data$tax_occ<- calc_n_samples(obj = taxmap_obj, 
                                             data = "tax_abund", 
                                             cols = taxmap_obj$data$sample_data$sample_id,
                                             groups = taxmap_obj$data$sample_data[colnames(taxmap_obj$data$sample_data)==sample_group][[1]]) 
    
    # calculate log2 median ratios and p values for a wilcoxcon test within taxas in this stress treatment groups
    # the sample groups needs some wrangling to fit within the soft code of the function
    taxmap_obj$data$diff_table <- compare_groups(obj = taxmap_obj,
                                                 data = "tax_abund",
                                                 cols = taxmap_obj$data$sample_data$sample_id, 
                                                 groups = taxmap_obj$data$sample_data[colnames(taxmap_obj$data$sample_data)==sample_group][[1]]) 
    
    # set differential log ratio to 0 based on adjusted p values
     taxmap_obj$data$diff_table$log2_median_ratio[taxmap_obj$data$diff_table$wilcox_p_value > 0.05] <- 0
    
    # define number of compared factors
    factors_compared<-taxmap_obj$data$sample_data[colnames(taxmap_obj$data$sample_data)==sample_group][[1]] 
    
    # draw the plot based on an if else statement: if there are 2 groups, plot a a heat tree comparing abundances between both groups, else plot a matrix of ehat trees. this function will fail if you only have 1 sample group! 
    
    if (length(unique(factors_compared)) == 2) {
      
      set.seed(1)
      output<- taxmap_obj %>%
        heat_tree(
          node_label = taxon_names,
                data = "diff_table", 
          node_size = n_obs,
          node_color = mean_diff,
          node_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
          node_color_range = c("cyan", "gray", "tan"), # The color palette used
          layout = "davidson-harel",
          initial_layout = "reingold-tilford")
      
    } else {
      
      set.seed(1)
      output<-heat_tree_matrix(taxmap_obj,
                               data = "diff_table", # this is the table with the data you want to plot
                               node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                               node_label = taxon_names,
                               node_color = log2_median_ratio, # A column from `taxmap_obj$data$diff_table_3treatments`
                               node_color_range = diverging_palette(), # The built-in palette for diverging data
                               node_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                               edge_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                               node_size_axis_label = "Number of OTUs",
                               node_color_axis_label = "Log2 ratio median proportions",
                               layout = "davidson-harel", # The primary layout algorithm
                               initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations
      
      
      
    }
    
    
    # clearly define the output object you will get from the function
    return(output)   
    
    
  }
  
  physeq_rarefied_filt_l<-lapply(physeq_rarefied_l, function(x) subset_taxa(physeq = x, Kingdom == "k__Bacteria"))


  # apply filter function on list of phyloseq objects - account for at least 0.5% of the reads in a sample & be present in at least 25% of the samples
  physeq_rarefied_filt_l<-lapply(physeq_rarefied_filt_l, function(x) filterPhyseq(x, 0.005, 25))
  
  
  heat_tree_microbe_l<-lapply(physeq_rarefied_filt_l[1:4], function(x)
    phyloseq_to_heat_tree_matrix2(ps_object = x, 
                                 sample_group = "Treatment"))
  

  heat_tree_array<-ggarrange(plotlist = heat_tree_microbe_l, labels = c("Week 4", "week 8", "week 12", "week 16"))

  #save shannon plot
  ggsave(plot = heat_tree_array, 
         filename = "./Results/microbiome_heat_tree_array_field.pdf",
         width = 600,
         height = 600,
         units = "mm")
  
  
  
  
  
  