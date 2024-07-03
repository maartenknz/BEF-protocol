#### BEF Course Analysis data Leipzig Canopy Crane 2023 ####

##this is not detailed skript as you are supposed to programm your code with the help of other scripts alone for this project. Yet some suggestions

#### load data and do first checks ####
setwd("C:/Uni/Master/BEF/practical/Analysis")

library(ggplot2) #for making graphics
library(cowplot) #for making multiplots of your graphics
library(tidyr) #we are going to use this for going from long to wide datasets
library(viridis) ##for fancy color palettes
library(emmeans) #for conducting post hoc tests
library(tidyverse)
#source("SummaryFunction.R") #for calculating SE and SD
library(patchwork) # for facetted graphs
library(lmerTest)
library(performance)
library(Rcpp)
library(vegan)
library(ggfortify)# needed to plot PCA
library(plotly)# needed to plot PCA
library(lme4)
library(magrittr)

rm(list=ls()) # clear history

#### data transformation ####

weights<-read.csv("Data input - weights.csv",dec=",")
names(weights)<-c("treeID","spec","species","total_LDW","total_LA","mean_LA")
stomata<-read.csv("Data input - stomata.csv",dec=",")
names(stomata)<-c("treeID",	"species",	"Leaf #",	"stom_dens_1",	"stom_dens_2",	"stom_dens_3",	"stom_length_1",	"stom_length_2",	"stom_length_3",	"stom_length_4",	"stom_length_5")
all<-read.csv("Data input - data all.csv",dec=",")
all[,9]<-as.numeric(gsub(",", ".", gsub("\\.", "", all[,9])))

colnames(all)
names(all)<-c("treeID",	"spec",	"species",	"leaf_nr",	"press_bar",	"Comments",	"stom_conduct_gsw","Transp_E",	"Quant_eff_dark_Fv_Fm",	"ETR",	"stom_dens_1",	"stom_dens_2",	"stom_dens_3",	"stom_length_1",	"stom_length_2",	"stom_length_3",	"stom_length_4",	"stom_length_5","total_LDW","total_LA","mean_LA")



lcc_tree <- read.csv("LCC_data_inventory_BEF_2024.csv",sep=";",header=T, dec=",")[1:15,]
lcc_tree <- lcc_tree[,1:13]
names(lcc_tree)
names(lcc_tree)<- c("FID","ID","treeID","spec","species","DBH_2014","hcbasal_2014","hctotal_2014" ,"DBH_2020","Growth14_20","hcbasal_2020","hctotal_2020","drought_estimate")



#calculating SLA from leaf area and dry weight 
all$sla<-all$total_LA/(all$total_LDW/1000)
all$stomatal_dens<-round(apply(all[,11:13],MARGIN=1,FUN=mean),1)
all$stomatal_length<-round(apply(all[,14:18],MARGIN=1,FUN=mean),1)

colnames(all)

#lcc_leafmean <-alternative way to calculate with summarize the leaf a data
leaf_all_a_new <- all[,c(1,5,7,8,9,10,21:24)]

lcc_leafmeannew_a <- leaf_all_a_new %>%
  group_by(treeID)%>%
  summarise_all(funs(mean), na.rm=T) %>%
  ungroup()

#### no combine to full data set on tree individual level
lcc_all_tree <-merge(lcc_tree,lcc_leafmeannew_a,by=c("treeID"),all.x=T)
names(lcc_all_tree)

fov <- 18 /40
fov.area <- (fov/2)^2 * pi


lcc_all_tree$stomatal_density <- lcc_all_tree$stomatal_dens / fov.area
lcc_all_tree$stomatal_length <- lcc_all_tree$stomatal_length * 2.5
lcc_all_tree$species[lcc_all_tree$species == "Til_pla"] <- "Til_cor"
# write.table(lcc_all_tree,"data_tree_mean.csv",sep=",")##speichert jede Tabelle

###################################################################################
#### Data analysis and plotting ####
##################################################################################

# lcc_all_tree <-read.csv("data_tree_mean.csv",sep=",",header=T, dec=".")
#### one version of a box plot

boxplot_DBH14 <- ggplot(lcc_all_tree, aes(x=species, y=DBH_2014))+
    geom_boxplot()+
    theme_bw()+
    ylab("DBH 2014\n")+
    xlab("\nSpecies")+
    theme(axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title = element_text(size =14),
          panel.grid=element_blank(),
          plot.margin=unit(c(1,1,1,1),units=,"cm"),
          legend.position = "none")

### analysis species differences with AOV####
model1 <- aov(DBH_2014~species,data=lcc_all_tree)
summary(model1)
TukeyHSD(model1)
plot(TukeyHSD(model1))


#### linear model of traits predicting drought responses in 2028 ####
model_drought1 <- lm(drought_estimate ~ stomatal_length, data = lcc_all_tree)
summary(model_drought1)

### ggplot example graph ####
ggplot(aes(x=stomatal_length, y=drought_estimate), data = lcc_all_tree)+  # this is the base command to create a ggplot 
  geom_point()+ # add points for each datapoint
  geom_smooth(method=lm)+ #add a regression line, 'lm' specifies a linear regression
  xlab("Stomatal Length (?m)")+ #add x-axis title
  ylab("Individual Sapflow")+ #add y-axis title
  theme_bw()+ #this is preset for design of a ggplot
  theme(axis.text = element_text(size=14), #size of the labels on the axes
        axis.title = element_text(size=18), # size of the titles of the axes
        title = element_text(size= 22), #size of th emain title
        plot.title.position = "panel") # centre the title above the panel


#### PCA with new input table and means ####
tree_means <- lcc_all_tree[order(lcc_all_tree$species),] 
str(tree_means)
tree_means$species

#### Principal component analysis ####
#you need ggfortify to execute the autoplot command here
## check here for help to visualize PCA results: https://plotly.com/ggplot2/pca-visualization/
names(tree_means)
results_pca <- prcomp(tree_means[,c(8,9,15:16,18:19,21:22)], center=T, scale = TRUE)
summary(results_pca)
p <- autoplot(results_pca)
p <- autoplot(results_pca,data=tree_means,colour="species",size=3,
         loadings=T,loadings.colour="black",loadings.label=T, loadings.label.size=3,loadings.label.colour="black")

#two options for plotting
p
ggplotly(p)


#### how to extract the aces scores of the PCA ####
tree_means$PCA1 <- predict(results_pca,newdata=tree_means)[,1]
tree_means$PCA2 <- predict(results_pca,newdata=tree_means)[,2]

#### PCA can only work with a fully filled matrix. If you have NA in the data (like we have for water potential) you have to remove the entire data point (also for the other traits) ####
names(tree_means)
results_pca <- prcomp(na.omit(tree_means[,c(8,9,15:16,18:22)]), center=T, scale = TRUE)
summary(results_pca)
p <- autoplot(results_pca)
p <- autoplot(results_pca,data=tree_means,colour="species",size=3,
              loadings=T,loadings.colour="black",loadings.label=T, loadings.label.size=3,loadings.label.colour="black")

#two options for plotting
p
ggplotly(p)

#### aggregate data over species ####
tree_means <-read.csv("data_tree_mean.csv",sep=",",header=T, dec=".")
str(tree_means)
lcc_treemean <- aggregate(cbind(stom_conduct_gsw, Transp_E,Quant_eff_dark_Fv_Fm, press_bar, stomatal_dens, stomatal_length, sla, DBH_2014, Growth14_20,drought_estimate)~species, data=tree_means, mean) 

model_drought2 <- lm( drought_estimate ~ stomatal_length, data = lcc_treemean)
summary(model_drought2)

ggplot(aes(x=stomatal_length, y=drought_estimate), data = lcc_treemean)+  # this is the base command to create a ggplot 
  geom_point()+ # add points for each datapoint
  geom_smooth(method=lm)+ #add a regression line, 'lm' specifies a linear regression
  xlab("Stomatal Length (?m)")+ #add x-axis title
  ylab("Sapflow_reduced")+ #add y-axis title
  theme_bw()+ #this is preset for design of a ggplot
  theme(axis.text = element_text(size=14), #size of the labels on the axes
        axis.title = element_text(size=18), # size of the titles of the axes
        title = element_text(size= 22), #size of th emain title
        plot.title.position = "panel") # centre the title above the panel

##############################################################

#### Group data analysis ####

#### linear modelling ####

# model growth ~ height
model.height <- lm(Growth14_20 ~ hctotal_2014, data = lcc_all_tree)
summary(model.height)
plot(model.height) # cooks distance problem #2

# model growth ~ DBH
model.DBH <- lm(Growth14_20 ~ DBH_2014, data = lcc_all_tree)
summary(model.DBH)
plot(model.DBH) # cooks distance problem #2

# model growth ~ Drought
model.growth <- lm(drought_estimate ~ Growth14_20, data = lcc_all_tree)
summary(model.growth)
plot(model.growth) # cooks distance problem #4

# model Drought ~ SLA
model.sla <-  lm(drought_estimate ~ sla, data = lcc_all_tree)
summary(model.sla) # marginally significant
plot(model.sla)

# model Drought ~ Stomatal density
model.density <-  lm(drought_estimate ~ stomatal_density, data = lcc_all_tree)
summary(model.density)
plot(model.density)

# model Drought ~ Stomatal size
model.length <-  lm(drought_estimate ~ stomatal_length, data = lcc_all_tree)
summary(model.length)
plot(model.length) # cooks distance problem #1 and 8

# model Drought ~ Stomatal conductance
model.conductance <-  lm(drought_estimate ~ stom_conduct_gsw, data = lcc_all_tree)
summary(model.conductance)
plot(model.conductance)

# model Drought ~ Leaf water potential
model.lwp <-  lm(drought_estimate ~ press_bar, data = lcc_all_tree)
summary(model.lwp) # significant
plot(model.lwp) # cooks distance problem #8 and 12

stats <- summary(model.lwp)$coefficients

#### plots ####

# Drought ~ SLA

ggplot(aes(x=sla, y=drought_estimate), data = lcc_all_tree)+  # this is the base command to create a ggplot 
  geom_point(size=3, aes(color=species))+ # add points for each datapoint
  geom_smooth (method=lm)+ #add a regression line, 'lm' specifies a linear regression
  xlab("SLA")+ #add x-axis title
  ylab("Drought resistance") + #add y-axis title
  ggtitle("Influence of specific leaf area on drought resistance \nof forest tree species") +
  theme_bw()+ #this is preset for design of a ggplot
  theme(axis.text = element_text(size=14), #size of the labels on the axes
        axis.title = element_text(size=18), # size of the titles of the axes
        title = element_text(size= 24), #size of the main title
        legend.text = element_text(size=14))


# Drought ~ Leaf water potential

ggplot(aes(x=press_bar, y=drought_estimate), data = lcc_all_tree)+  # this is the base command to create a ggplot 
  geom_point(size=3, aes(color=species))+ # add points for each datapoint
  geom_smooth(method=lm)+ #add a regression line, 'lm' specifies a linear regression
  xlab("Leaf water potential (bar)")+ #add x-axis title
  ylab("Drought resistance")+ #add y-axis title
  ggtitle("Influence of leaf water potential on drought resistance \nof forest tree species") +
  theme_bw()+ #this is preset for design of a ggplot
  theme(axis.text = element_text(size=14), #size of the labels on the axes
        axis.title = element_text(size=18), # size of the titles of the axes
        title = element_text(size= 24), #size of the main title
        legend.text = element_text(size=14))


# Drought ~ leaf water potential with species effect

ggplot(aes(x=press_bar, y=drought_estimate ,colour = spec), data = lcc_all_tree)+  # this is the base command to create a ggplot 
  geom_point(size=3)+ # add points for each datapoint
  geom_smooth(method=lm, )+ #add a regression line, 'lm' specifies a linear regression
  xlab("Leaf water potential (bar)")+ #add x-axis title
  ylab("Drought resistance")+ #add y-axis title
  ggtitle("Influence of leaf water potential on drought resistance \nof forest tree species") +
  theme_bw()+ #this is preset for design of a ggplot
  theme(axis.text = element_text(size=14), #size of the labels on the axes
        axis.title = element_text(size=18), # size of the titles of the axes
        title = element_text(size= 24), #size of the main title
  )


#### PCA for all traits ####

names(lcc_all_tree)
results_pca <- prcomp(lcc_all_tree[,c(6,8,10,14:15,20,22:23)], center=T, scale = TRUE)
summary(results_pca)
p <- autoplot(results_pca)
p <- autoplot(results_pca,data=lcc_all_tree,colour="species",size=3,
              loadings=T,loadings.colour="black",loadings.label=T, loadings.label.size=3,loadings.label.colour="black")

lcc_no_frax <- subset(lcc_all_tree, species != "Fra_exc")

names(lcc_all_tree)
results_pca <- prcomp(lcc_all_tree[,c(6,8,14:15,20,22:23)], center=T, scale = TRUE)
summary(results_pca)
p <- autoplot(results_pca)
p <- autoplot(results_pca,data=lcc_all_tree,colour="species",size=3,
              loadings=T,loadings.colour="black",loadings.label=T, loadings.label.size=3,loadings.label.colour="black")

#two options for plotting
p
ggplotly(p)

test <- na.omit(tree_means[,c(5,9:18,20:22)])

lcc_all_tree$PCA1 <- predict(results_pca,newdata=lcc_all_tree)[,1]
lcc_all_tree$PCA2 <- predict(results_pca,newdata=lcc_all_tree)[,2]

names(lcc_no_frax)
results_pca <- prcomp(lcc_no_frax[,c(6,8,10,14:15,20,22:23)], center=T, scale = TRUE)
summary(results_pca)
p <- autoplot(results_pca)
p <- autoplot(results_pca,data=lcc_no_frax,colour="species",size=3,
              loadings=T,loadings.colour="black",loadings.label=T, loadings.label.size=3,loadings.label.colour="black")
p
