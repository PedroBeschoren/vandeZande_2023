

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

# pacakge for microbiome time series: https://microbiome.github.io/miaTime/


#load inssect data
insect<-read.csv2(file = "./Data/Data_Sheets_exuviae_observations_complete1.csv")

#remove NA
insect<-na.omit(insect)
insect<-insect[complete.cases(insect), ]

#subset
insect_metadata<-(insect[,1:8])

#transform to pyoseq
insect_ps<-sample_data(insect_metadata)

#load insect species table and save as phyloseq
insect_sp<-insect[,9:17]
insect_ps_otu<-otu_table(t(insect_sp), taxa_are_rows = TRUE)

sample_names(insect_ps_otu)<-sample_names(insect_ps)

insect_ps<-merge_phyloseq(insect_ps_otu, insect_ps)

otu_table(insect_ps)[1:3,1:3]
insect_ps@otu_table

#NMDS of isnect species
set.seed(101)
nmds_insect <- phyloseq::ordinate(insect_ps,
                                       method = "PCoA", # this method has few assumptions and readly accepts different data structures
                                       distance = "bray", # bray-curtis distance is suitable for sparse data - such as the zero-inflated microbiome data we have
                                       autotransform = TRUE
) # automatically transforms your data, if needed. reduces weight of outliers
# weakties = FALSE prevests tress from colapsing into zero

NMDS_all_samples<-plot_ordination(
  physeq = physeq_rarefied,
  ordination = nmds_rootAndSoil,
  color = "Treatment",
  shape = "Soil_type"
) +
  theme_classic() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(aes(size = Time_point), alpha = 1) +
  theme(legend.position = "right")

#this is the inseq phyloseq object
insect_ps

#load data
physeq1 <- import_biom("./Data/featureTable_1.biom")
physeq2 <- import_biom("./Data/featureTable_2.biom")

#merge both field and potting soil
physeq<-merge_phyloseq(physeq1, physeq2)

#laod metadata
metadata <- read.csv2(file = "./Data/metadata.csv", row.names = 3)
metadata$Time_point<-factor(x = metadata$Time_point, levels = c("Week.4", "Week.8", "Week.12", "Week.16"))
metadata <-sample_data(metadata)



#merge metadata with  otu table
physeq<-merge_phyloseq(physeq, metadata)

#remove rare ASVs
otu_table(physeq) <- otu_table(physeq)[which(rowSums(otu_table(physeq)) > 7), ] # this drops ~3000 taxa.

#check library sizes
sample_sums(physeq)%>%hist
sample_sums(physeq)%>%sort(decreasing = FALSE)

#normalize data
physeq_filtered_df <- as.data.frame(otu_table(physeq))
rarecurve(t(physeq_filtered_df),
  #        col = sample_data(physeq)$Soil_type,
          label = FALSE,
          step = 600,
          main = "Rarefaction ", ylab = "Number of ASVs", xlab = "Number of DNA Sequences",
          abline(v = min(sample_sums(physeq)), lwd = 3, lty = 2))


# now, let's rarefy the data
set.seed(100) # set a random seed so that whenever you re-run this code you draw the same set of OTUs
# runt he rarefaction. with 1 argument per line it becomes easier to see what is going on. naming the arguments also help a lot!
physeq_rarefied <- rarefy_even_depth(
                                      physeq = physeq,
                                      sample.size = min(sample_sums(physeq)),
                                      rngseed = FALSE,
                                      replace = TRUE,
                                      trimOTUs = TRUE,
                                      verbose = TRUE
                                    )



#

# Let's make a Non-Metric Multidimensional Scaling (NMDS) of all our samples based on CSS normalization
set.seed(101)
nmds_rootAndSoil <- phyloseq::ordinate(physeq_rarefied,
                                       method = "NMDS", # this method has few assumptions and readly accepts different data structures
                                       distance = "bray", # bray-curtis distance is suitable for sparse data - such as the zero-inflated microbiome data we have
                                       trymax=200,
                                       try = 200,
                                       autotransform = TRUE
) # automatically transforms your data, if needed. reduces weight of outliers
# weakties = FALSE prevests tress from colapsing into zero

NMDS_all_samples<-plot_ordination(
  physeq = physeq_rarefied,
  ordination = nmds_rootAndSoil,
  color = "Treatment",
  shape = "Soil_type"
) +
  theme_classic() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(aes(size = Time_point), alpha = 1) +
  theme(legend.position = "right")

ggsave(NMDS_all_samples, filename = "./Results/NMDS_all_samples.pdf",
       width = 180,
       height = 180,
       units = "mm")

################# separate phyloseq objects

physeq_rarefied_l<-phyloseq_sep_variable(physeq_rarefied, variable = c("Soil_type", "Time_point"))

# define custom fucntion to plot a list of NMDS results
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
      theme(legend.position = "right")
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
  
  
  return(list(bet_disp_PCOa_plot, bet_disp_boxplot))
}

# now that we have a custom function, we can run it across all lists and variables
library(purrr)
set.seed(5235)
dip_result3 <- beta_disp_plotAndTest(physeq_rarefied_l, "Treatment")
dip_result3[[1]]

disppersion_result<-ggarrange(plotlist = dip_result3[[2]],ncol = 4, nrow = 2)  






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
             path = "./Results/PERMANOVA_tables_8plots.docx")



#pairwise permanova
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

# filter df
control_df<-filter(pairwise_permanova_df, treatment_1 =="Control" | treatment_2 == "Control" )
control_df$combination<-as.factor(control_df$combination)

#line 
ggplot(control_df, aes(x = Time ,y = R2))+
  geom_line(aes(color = combination))+
  facet_wrap(~Site)

########## alpha diversity

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
    labs(title = "Shannon diversity")+
    ylab("Shannon diversity index")+
    xlab("Time points")+
    facet_wrap(~Soil_type)

#save shannon plot
ggsave(shannon_topsoil_plot, filename = "./Results/shannon_plot.pdf",
       width = 180,
       height = 180,
       units = "mm")





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
  anova_test<-flextable(anova_test)
  
  #flextable of pairwside comparsions
  save_as_docx(shannon_anova_ft,
               path = "./Results/shannon_anova.docx")
  
  ###########
  # run miaTime to have a better view of time effects
  
  physeq_rarefied_soiltype_l<-phyloseq_sep_variable(physeq_rarefied, variable = c("Soil_type"))
  
  #convert phyloseq into summarized tree experiment
  tse_l<-lapply(physeq_rarefied_soiltype_l, function (x) makeTreeSEFromPhyloseq(x))
  
  mia::taxonomyRanks(tse_l$Field) 
  
  #adjust time formating in the table
  tse_l<- 
  lapply(tse_l, function (x){
    x@colData$Time_point<-   as.numeric(gsub(pattern = "Week.", replacement = "", x = x@colData$Time_point))
return(x)                                        
  })
  
  
  
  
  tse_field <- getStepwiseDivergence(tse_l$Field,
                                     group = "Treatment",
                                     time_field = "Time_point",
                                     time_interval = 4,
                                     name_divergence = "time_divergence",
                                     name_timedifference = "time_difference",
                                     assay.type="counts",
                                     FUN = vegan::vegdist,
                                     method="bray") 
  
  
  field_stepwise_plot<-as.data.frame(tse_field@colData) %>%
    ggplot(aes(x=Time_point, y=time_divergence))+
    geom_point(aes(color=Treatment), position = position_jitterdodge(jitter.width = 0.1,
                                                                     dodge.width = 0.7))+
    geom_smooth()+
    theme_bw()+
    ylim(0.3, 0.9)+
    labs(x="Time (weeks) in field", y="Stepwise time divergence") +
    facet_wrap(~Treatment)
  

  
  
  
  
  
  
 #pot 
  tse_pot <- getStepwiseDivergence(tse_l$Pot,
                                     group = "Treatment",
                                     time_field = "Time_point",
                                     time_interval = 4,
                                     name_divergence = "time_divergence",
                                     name_timedifference = "time_difference",
                                     assay.type="counts",
                                     FUN = vegan::vegdist,
                                     method="bray") 
  
  
  pot_stepwise_plot<-as.data.frame(tse_pot@colData) %>%
    ggplot(aes(x=Time_point, y=time_divergence))+
    geom_point(aes(color=Treatment), position = position_jitterdodge(jitter.width = 0.1,
                                                                     dodge.width = 0.7))+
    geom_smooth()+
    theme_bw()+
    labs(x="Time (weeks) in pots", y="Stepwise time divergence") +
    ylim(0.3, 0.9)+
    facet_wrap(~Treatment)
  
  #put both plots in one pannel adn export
  stepwise_time_plot<-
  ggarrange(field_stepwise_plot, pot_stepwise_plot, labels = "AUTO", common.legend = TRUE)

 ggsave(stepwise_time_plot, file = "./stepwise_time_plot.pdf", width = 180, height = 120, units = "mm")

 
 # check specific microbes of interest
 tse_pot$ASV1
 
 library(miaViz)
 plotSeries(tse_pot,
            x = "Time_point",
            y = c("ASV1", "ASV2"),
            assay.type = "counts")+
   theme_bw()
 
 
 
  