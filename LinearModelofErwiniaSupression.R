library(dplyr)
library(ggplot2)
library(ggpubr)
library(FactoMineR)
library(factoextra)
library(corrr)
library(ggcorrplot)
library(clusterSim)
library(tidyverse)
library(tidyr)
library(ggpmisc)
library(EcoSimR)
req <- substitute(require(x, character.only = TRUE))
libs<-c("sjPlot", "stargazer", "coefplot", "dotwhisker", "visreg", "jtools", "ggfortify", "olsrr", "DescTools", "interactions")
sapply(libs, function(x) eval(req) || {install.packages(x); eval(req)})

##################
## Read in data ##
##################

daway <- "~/Erwinia2022/ErwiniaR-2.csv"
CSV <- read.csv(daway)
data<-CSV
head(data)

#create the response variables 
data$log10raw<-log10(data$cfus+1)
data$log10av<-log10(data$cfu+1)

# growth metrics are significantly correlated, so we can collapse growth into one metric
Growthdata <- data %>% dplyr::select(r:t_mid)
GrowCorMat <- cor(Growthdata, method = "spearman")
ggcorrplot(GrowCorMat)

# resource use matrix
res.use <- data %>% dplyr::select(GluNetChange_mg.uL:SucNetChange_mg.uL) 
ResCorrMat <- cor(res.use, method= "spearman")
ggcorrplot(ResCorrMat)
ggscatter(res.use, x = "SucNetChange_mg.uL", y = "FruNetChange_mg.uL", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Glucose", ylab = "Sucrose")


# normalize the data 
Growthdata.norm <- data.Normalization(Growthdata, type = "n1", normalization = "column")
qqnorm(Growthdata.norm$t_mid)
hist(Growthdata.norm$r)
Resuse.norm <- data.Normalization(res.use, type = "n1", normalization = "column")


# perform the PCA 
Growthdata.PCA <- prcomp(Growthdata.norm, center = TRUE, scale = TRUE)
summary(Growthdata.PCA)


Resuse.PCA <- prcomp(Resuse.norm, center = TRUE, scale = TRUE)
summary(Resuse.PCA)


# Checking that there is no longer correlation among the variables 
res1 <- cor(Growthdata.PCA$x, method="pearson")
corrplot::corrplot(res1, method= "color", order = "hclust", tl.pos = 'n')

res2 <- cor(Resuse.PCA$x, method="pearson")
corrplot::corrplot(res1, method= "color", order = "hclust", tl.pos = 'n')

# Loading plots, this tells us how strong each variable is correlated to each PC 
fviz_pca_var(Growthdata.PCA,axes = c(1, 2))
fviz_pca_var(Resuse.PCA,axes = c(1, 2))

#Pull out the PCs into a data frame to later use in the linear model 
Growth.PCs <- as.data.frame(Growthdata.PCA$x)
Resuse.PCs <- as.data.frame(Resuse.PCA$x) %>% rename("R.PC1" = "PC1", "R.PC2" = "PC2", "R.PC3" = "PC3")
Growth.PCs <- Growth.PCs %>% rename("G.PC1"= "PC1" ,"G.PC2"= "PC2","G.PC3"= "PC3")




########################
## Read in amino data ##
########################

## Amino data 
setwd("C:/Users/A02342347/Documents/Erwinia2022/Amino data")
ConsumptionData <- read.csv("ConsumptionR.csv")

df1 <- ConsumptionData %>% mutate(background = as.factor(background), ABA = as.numeric(AmountABA.nmol.), 
                                   ASN = as.numeric(AmountASN), ASP = as.numeric(ASP),
                                   GLN.HIS = as.numeric(GLN.HIS),GLU = as.numeric(GLU),
                                   GLY = as.numeric(GLY),ILE = as.numeric(ILE),
                                   LEU = as.numeric(LEU),LYS = as.numeric(LYS),
                                   PHE = as.numeric(PHE),SER = as.numeric(SER),
                                   THR = as.numeric(THR),TRP = as.numeric(TRP),
                                   VAL.MET = as.numeric(VAL.MET)) 
df2 <- df1 %>%  
  mutate(TotalAminos = rowSums(across(c("ASP":"ASN")), na.rm=TRUE)) 

##############################################################################
#### Linear Models ####
##############################################################################

# Subset/merge dataframes
##### May have to change directory here ####
setwd("~/Erwinia2022")
phydist <- read.csv("Phydist.csv")

data2 <- cbind(Growth.PCs, Resuse.PCs, data)
data_merge <- merge.data.frame(data2, df2, by = c("name", "background", "id"))
data_merge$background <- as.factor(data_merge$background)

data_bact <- data_merge %>% filter(kingdom == "Bacteria")
data_bact2 <- merge.data.frame(data_bact, phydist, by = "name")
data_bact2$background <- as.factor(data_bact2$background)

data_fungi <- data_merge %>% filter(kingdom == "Fungi")
data_fungi$background <- as.factor(data_fungi$background)




##################
## Creating LMs ##
##################

# Overall model

SugarMoc <- lm(SucNetChange_mg.uL ~ name, data = data_merge)
summary(SugarMoc)
pHModel <- lm(pH ~ name, data = data_merge)
summary(pHModel)

   
ModelOverlap <- lm(Overlap ~ASP + GLN.HIS + GLU + ASN  + ILE + LEU  + SER + THR  + GLY  + PHE  + VAL.MET, data = data_merge)
summary(ModelOverlap)

Model1 <- lm(log10raw ~ Overlap + TotalAminos + pH + G.PC1 + R.PC1 + R.PC2 + kingdom + background, data = data_merge)
summary(Model1)

# Model subset for bacteria
ModelBact1 <- lm(log10raw ~  Overlap + Dist + TotalAminos + pH + G.PC1 + R.PC1 + R.PC2  + background, data = data_bact2)
summary(ModelBact1)

# Model subset for fungi
ModelFungi1 <- lm(log10raw ~ TotalAminos + Overlap + pH + G.PC2  + R.PC1 + R.PC2 + background , data = data_fungi)
summary(ModelFungi1)

overl <- data_bact2 %>% dplyr::select(Overlap, Dist)
cor(over)
ggscatter(overl, x = "Overlap", y = "Dist", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Overlap", ylab = "Dist")



### Model selection 
library(olsrr)
ggplot2::autoplot(Model1, which = 1:6, ncol = 2) ### Assumptions do not seem to be met!
olsrr::ols_test_breusch_pagan(Model1) # Test of Homo/Heteroskedasticity
a <- ols_step_best_subset(Model1)
a

ggplot2::autoplot(ModelBact1, which = 1:6, ncol = 2) ### Assumptions do not seem to be met!
olsrr::ols_test_breusch_pagan(ModelBact1) # Test of Homo/Heteroskedasticity
b <- ols_step_best_subset(ModelBact1)
b

ggplot2::autoplot(ModelFungi1, which = 1:6, ncol = 2) ### Assumptions do not seem to be met!
olsrr::ols_test_breusch_pagan(ModelFungi1) # Test of Homo/Heteroskedasticity
c <- ols_step_best_subset(ModelFungi1)
c

######### Best subset models 
# Overall model
Model2 <- lm(log10raw ~ TotalAminos + pH + R.PC2 + background, data = data_merge)
summary(Model2)


# Model subset for bacteria
ModelBact2 <- lm(log10raw ~ Overlap + Dist + background + R.PC1 + R.PC2, data = data_bact2)
summary(ModelBact2)



# Model subset for fungi
ModelFungi2 <- lm(log10raw ~  Overlap+ pH + G.PC2  + R.PC1 + background, data = data_fungi)
summary(ModelFungi2)

##### Model Comparison
sjPlot::tab_model(Model1, Model2, pred.labels = c("Intercept", "iv1 label", "iv2 label"), 
                  dv.labels = c("Model 1 overall", "Model 2 overall"), string.pred = "Estimates", string.ci = "95% CIs", 
                  string.p = "p-values")

sjPlot::tab_model(ModelBact1, ModelBact2, pred.labels = c("Intercept", "iv1 label", "iv2 label"), 
                  dv.labels = c("Model 1 Bacteria", "Model 2 Bacteria"), string.pred = "Estimates", string.ci = "95% CIs", 
                  string.p = "p-values")

sjPlot::tab_model(ModelFungi1, ModelFungi2, pred.labels = c("Intercept", "iv1 label", "iv2 label"), 
                  dv.labels = c("Model 1 Fungi", "Model 2 Fungi"), string.pred = "Estimates", string.ci = "95% CIs", 
                  string.p = "p-values")


####### Two way interaction models 
Modelint <- lm(log10raw ~ (Overlap + TotalAminos + pH  + R.PC1 + R.PC2 + background)^2 + G.PC1 + G.PC2 + kingdom , data = data_merge)
summary(Modelint)

ModelBactInt <- lm(log10raw ~ (Overlap + TotalAminos + pH  + R.PC1 + R.PC2 + background)^2 + Dist + G.PC1 + G.PC2, data = data_bact3)
summary(ModelBactInt)


ModelFungiInt <- lm(log10raw ~ (Overlap + TotalAminos + pH  + R.PC1 + R.PC2 + background)^2 + G.PC1 + G.PC2   , data = df_fungi)
summary(ModelFungiInt)


###########################################
## Figures for predictors for each model ##
###########################################
library(patchwork)
library(grid)
library(gridExtra)
library(lemon)
### Overall model 

g <- guides(fill = guide_legend(override.aes = list(color = scales::hue_pal()(3),
                                                    shape = c(16, 16, 16), 
                                                    size = c(1, 1, 1),
                                                    alpha = c(1, 1, 1)),))
SucConsumption <- ggplot(data = data_merge, aes(x = background, y = SucNetChange_mg.uL)) + geom_boxplot() 

pHlin <- ggplot(data = data_merge, aes(x = pH, y = log10raw, color = background))+ geom_point(size = 3) + geom_smooth(method="lm", se=FALSE, span=0.95) + scale_color_manual(name="Nectar background (w/v)", values=c("#DCE319FF", "#287D8EFF")) + labs(x="pH") + ylab(expression("Erwinia growth (Log"["10"] ~"CFUs + 1)")) + theme(axis.text.x=element_text(size=14, color="black"), axis.title.x=element_text(size=14, color="black"),axis.text.y=element_text(size=14, color="black"), axis.title.y=element_text(size=16)) + theme(legend.position="none")+g


RPC2lin <- ggplot(data = data_merge, aes(x = R.PC2, y = log10raw, color = background))+ geom_point(size = 3) +
  geom_smooth(method="lm", se=FALSE, span=0.95) + scale_color_manual(name="Nectar background (w/v)", values=c("#DCE319FF", "#287D8EFF"))+ labs(x="PC2 (Sucrose concentration)") + ylab(expression("Erwinia growth (Log"["10"] ~"CFUs + 1)")) + theme(axis.text.x=element_text(size=14, color="black"), axis.title.x=element_text(size=14, color="black"), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none")+g

background <- ggplot(data = data_merge, aes(x = background, y = log10raw, color = background))+ geom_boxplot() + scale_color_manual(name="Nectar background (w/v)", values=c("#DCE319FF", "#287D8EFF")) + 
  labs(x="Background") + ylab(expression("Erwinia growth (Log"["10"] ~"CFUs + 1)")) + theme(axis.text.x=element_text(size=14, color="black"), axis.title.x=element_text(size=14, color="black"), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none")+g

TotalAminoslin <- ggplot(data = data_merge, aes(x = TotalAminos, y = log10raw, color=background))+ geom_point(size = 3) + geom_smooth(method="lm", se= FALSE, span=0.95) + scale_color_manual(name="Nectar background (w/v)", values=c("#DCE319FF", "#287D8EFF"))+
  labs(x="Total Amino Consumption") + ylab(expression("Erwinia growth (Log"["10"] ~"CFUs + 1)"))+ theme(axis.text.x=element_text(size=12, color="black"), axis.title.x=element_text(size=16, color="black"), axis.text.y=element_text(size=14, color="black"), axis.title.y=element_text(size=16))+theme(legend.position="none")+g


combined <- background + Overlaplin + pHlin + RPC2lin 
combined +  plot_layout(nrow=1, guides="collect") + theme(legend.position = c(1,1))

grid_arrange_shared_legend(TotalAminoslin, background, ncol = 2, nrow = 1, position='bottom')
grid_arrange_shared_legend(pHlin, RPC2lin, ncol = 2, nrow = 1, position='bottom')

####### Bacterial Model

Overlaplin <- ggplot(data = data_bact2, aes(x = Overlap, y = log10raw, color = background, shape = background))+ geom_point(size = 3) +
  geom_smooth(method="lm", se = FALSE, span=0.95) + geom_smooth(method="lm", se= FALSE, span=0.95) + scale_color_manual(name="Nectar background", values=c("#DCE319FF", "#287D8EFF")) + scale_shape_manual(name="Nectar background", values=c(16,18))+
  stat_poly_eq()+ labs(x="Niche overlap (shared amino use)") + ylab(expression("Erwinia growth (Log"["10"] ~"CFUs + 1)"))+ theme(axis.text.x=element_text(size=12, color="black"), axis.title.x=element_text(size=16, color="black"), axis.text.y=element_text(size=14, color="black"), axis.title.y=element_text(size=16))+theme(legend.position="none")+g
Overlaplin


Distlin <- ggplot(data = data_bact2, aes(x = Dist, y = log10raw, color = background, shape = background))+ geom_point(size = 3) +
  geom_smooth(method="lm", se=FALSE, span=0.95) +  scale_color_manual(name="Nectar background", values=c("#DCE319FF", "#287D8EFF")) + scale_shape_manual(name="Nectar background", values=c(16,18))+
  stat_poly_eq()+ labs(x="Phylogenetic Distance") + ylab(expression("Erwinia growth (Log"["10"] ~"CFUs + 1)")) + theme(axis.text.x=element_text(size=14, color="black"), axis.title.x=element_text(size=14, color="black"), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none")+g
Distlin


RPC1Bactlin <- ggplot(data = data_bact2, aes(x = R.PC1, y = log10raw, color = background, shape = background))+ geom_point(size = 3) +geom_smooth(method="lm", se=FALSE, span=0.95) + scale_color_manual(name="Nectar background", values=c("#DCE319FF", "#287D8EFF")) + scale_shape_manual(name="Nectar background", values=c(16,18))+
  stat_poly_eq()+ labs(x="PC1 (Glucose/Fructose concentration)") + ylab(expression("Erwinia growth (Log"["10"] ~"CFUs + 1)")) + theme(axis.text.x=element_text(size=14, color="black"), axis.title.x=element_text(size=14, color="black"), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none")+g
RPC1Bactlin

RPC2Bactlin <- ggplot(data = data_bact2, aes(x = R.PC2, y = log10raw, color = background, shape = background))+ geom_point(size = 3) +geom_smooth(method="lm", se=FALSE, span=0.95) + scale_color_manual(name="Nectar background", values=c("#DCE319FF", "#287D8EFF")) + scale_shape_manual(name="Nectar background", values=c(16,18))+
  stat_poly_eq()+ labs(x="PC2 Sucrose concentration") + ylab(expression("Erwinia growth (Log"["10"] ~"CFUs + 1)")) + theme(axis.text.x=element_text(size=12, color="black"), axis.title.x=element_text(size=16, color="black"), axis.text.y=element_text(size=14, color="black"), axis.title.y=element_text(size=16))+theme(legend.position="none")+g

RPC2Bactlin

backgroundbact <- ggplot(data = data_bact2, aes(x = background, y = log10raw, color = background, shape = background))+ geom_boxplot()  + scale_color_manual(name="Nectar background", values=c("#DCE319FF", "#287D8EFF")) + scale_shape_manual(name="Nectar background", values=c(16,18))+
  labs(x="background") + ylab(expression("Erwinia growth (Log"["10"] ~"CFUs + 1)")) + theme(axis.text.x=element_text(size=14, color="black"), axis.title.x=element_text(size=14, color="black"), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none")+g
backgroundbact


grid_arrange_shared_legend(Overlaplin, Distlin,  ncol = 2, nrow = 1, position='bottom')
grid_arrange_shared_legend(RPC2Bactlin,backgroundbact, ncol = 2, nrow = 1, position='bottom')



###### 4 panel figure 



over <- data_bact2 %>% dplyr::select(Dist, Overlap)

o <- cor(over)
ggscatter(over, x = "Overlap", y = "Dist", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Niche Overlap", ylab = "Phylogenetic distance")




########## Fungi plots

pHFungilin <- ggplot(data = data_fungi, aes(x = pH, y = log10raw,  color = background, shape = background))+ geom_point(size = 3) +
  geom_smooth(method="lm", se=FALSE, span=0.95) +  scale_color_manual(name="Nectar background", values=c("#DCE319FF", "#287D8EFF")) + scale_shape_manual(name="Nectar background", values=c(16,18))+
  stat_poly_eq()+ labs(x="pH") + ylab(expression("Erwinia growth (Log"["10"] ~"CFUs + 1)")) 
pHFungilin

OverlapFungilin <- ggplot(data = data_fungi, aes(x = Overlap, y = log10raw, color = background, shape = background))+ geom_point(size = 3) +geom_smooth(method="lm", se=FALSE, span=0.95) +  scale_color_manual(name="Nectar background", values=c("#DCE319FF", "#287D8EFF")) + scale_shape_manual(name="Nectar background", values=c(16,18))+
  stat_poly_eq()+ labs(x="Niche Overlap") + ylab(expression("Erwinia growth (Log"["10"] ~"CFUs + 1)")) 
OverlapFungilin

RPC1Fungilin <- ggplot(data = data_fungi, aes(x = R.PC1, y = log10raw, color = background, shape = background))+ geom_point(size = 3) +geom_smooth(method="lm", se=FALSE, span=0.95) + scale_color_manual(name="Nectar background", values=c("#DCE319FF", "#287D8EFF")) + scale_shape_manual(name="Nectar background", values=c(16,18))+
  stat_poly_eq()+ labs(x="Sugar use PC 1") + ylab(expression("Erwinia growth (Log"["10"] ~"CFUs + 1)")) 
RPC1Fungilin ## zeros pulling this one

GPC2Fungilin <- ggplot(data = data_fungi, aes(x = G.PC2, y = log10raw, color = background, shape = background))+ geom_point(size = 3) +geom_smooth(method="lm", se=FALSE, span=0.95) + scale_color_manual(name="Nectar background", values=c("#DCE319FF", "#287D8EFF")) + scale_shape_manual(name="Nectar background", values=c(16,18))+
  stat_poly_eq()+ labs(x="Growth PC 2") + ylab(expression("Erwinia growth (Log"["10"] ~"CFUs + 1)")) 
GPC2Fungilin ## zeros pulling this one

backgroundfung <- ggplot(data = data_fungi, aes(x = background, y = log10raw, color = background, shape = background))+ geom_boxplot()  + scale_color_manual(name="Nectar background", values=c("#DCE319FF", "#287D8EFF")) + scale_shape_manual(name="Nectar background", values=c(16,18))+
  labs(x="background") + ylab(expression("Erwinia growth (Log"["10"] ~"CFUs + 1)")) 
backgroundfung

######### Erwinia growth in the competition assay
data3 <- data %>% filter(name != "Erwinia amylovora")

                                                  
#########################################################################################################


p <- ggplot(data3, aes(log10av, y=factor(name, levels=rev(levels(factor(name)))))) +
  geom_line(aes(group = name)) +
  geom_point(aes(color = factor(background)), size=4) +
  labs(x=expression("Erwinia growth (Log"["10"] ~"CFUs + 1)")) +
  ylab("Nectar microbe") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.title.y=element_text(size=18, color = "black")) +
  theme(axis.text.y=element_text(size = 10, face="italic", color = "black")) +
  xlim(0, 10) +
  scale_color_manual(name="Nectar background", values=c("#DCE319FF", "#287D8EFF")) +
  theme(legend.position="bottom")

p <- p + geom_point(data=data3, aes(log10raw, y=factor(name, levels=rev(levels(factor(name)))), color=factor(background)), size=2, shape=21) + scale_shape(solid=FALSE) 
 
p

#####################################
## Pianka Overlap value generation ##
#####################################


### Data wrangling 
df_average7.5 <- df2 %>% group_by(id) %>% filter(background == "7.5") %>%
  summarise(asp = mean(ASP),gln_his = mean(GLN.HIS),glu = mean(GLU),gly = mean(GLY),ile = mean(ILE),leu = mean(LEU),lys = mean(LYS),
            phe = mean(PHE),ser = mean(SER),thr = mean(THR),trp = mean(TRP),val_met = mean(VAL.MET)) %>% replace(is.na(.), 0)

df_average15 <- df2 %>% group_by(id) %>% filter(background == "15")  %>%  
  summarise(asp = mean(ASP),gln_his = mean(GLN.HIS),glu = mean(GLU), gly = mean(GLY),ile = mean(ILE), leu = mean(LEU),lys = mean(LYS),
  phe = mean(PHE),ser = mean(SER), thr = mean(THR),trp = mean(TRP),val_met = mean(VAL.MET))%>% replace(is.na(.), 0)

df_average15ABS <- abs(df_average15[,2:13])
df_average7.5ABS <- abs(df_average7.5[,2:13])

aa15<- cbind(df_average15$id, df_average15ABS)
aa7.5 <- cbind(df_average7.5$id, df_average7.5ABS)

aa7.5 <- rename(aa7.5, "id" = "df_average7.5$id")
aa15 <- rename(aa15, "id" = "df_average15$id")



#### Creating pairwise matricies 
AsaY2_7.5 <- as.matrix(subset(aa7.5, id == "Ea8R"|id == "AsaY2"))
Eas7.5 <- aa7.5 %>% filter(id == "Ea8R"|id == "Ea8") %>% dplyr::select(2:13)%>% as.matrix()
AsaY2_7.5 <- aa7.5 %>% filter(id == "AsaY2" | id == "Ea8R") %>% dplyr::select(2:13)%>% as.matrix()
Asay1_7.5 <- aa7.5 %>% filter(id == "AsaY1" | id == "Ea8R") %>% dplyr::select(2:13)%>% as.matrix()
B.mega7.5<- aa7.5 %>% filter(id == "Ea8R"|id == "B. mega") %>% dplyr::select(2:13)%>% as.matrix()
B.sub7.5<- aa7.5 %>% filter(id == "Ea8R"|id == "B. sub") %>% dplyr::select(2:13)%>% as.matrix()
Baugh7.5 <- aa7.5 %>% filter(id == "Ea8R"|id == "Baugh_Org_Int_N_Y1") %>% dplyr::select(2:13)%>% as.matrix()
Champs7.5 <- aa7.5 %>% filter(id == "Ea8R"|id == "Champs_P_B2") %>% dplyr::select(2:13)%>% as.matrix()
EC005_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id == "EC_005") %>% dplyr::select(2:13)%>% as.matrix()
EC006_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id =="EC_006") %>% dplyr::select(2:13)%>% as.matrix()
EC033_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id == "EC_033") %>% dplyr::select(2:13)%>% as.matrix()
EC052_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id == "EC_052") %>% dplyr::select(2:13)%>% as.matrix()
EC061_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id =="EC_061") %>% dplyr::select(2:13)%>% as.matrix()
EC064_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id =="EC_064") %>%dplyr::select(2:13)%>% as.matrix()
EC067_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id =="EC_067") %>% dplyr::select(2:13)%>% as.matrix()
EC072_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id == "EC_072") %>% dplyr::select(2:13)%>% as.matrix()
EC088_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id =="EC_088") %>% dplyr::select(2:13)%>% as.matrix()
EC089_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id =="EC_089") %>% dplyr::select(2:13)%>% as.matrix()
EC094_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id =="EC_094") %>% dplyr::select(2:13)%>% as.matrix()
EC112_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id == "EC_112") %>% dplyr::select(2:13)%>% as.matrix()
EC115_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id =="EC_115") %>% dplyr::select(2:13)%>% as.matrix()
EC116_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id == "EC_116") %>% dplyr::select(2:13)%>% as.matrix()
EC124_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id == "EC_124") %>% dplyr::select(2:13)%>% as.matrix()
EC126_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id == "EC_126") %>% dplyr::select(2:13)%>% as.matrix()
EC137_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id == "EC_137") %>% dplyr::select(2:13)%>% as.matrix()
Gall7.5<- aa7.5 %>% filter(id == "Ea8R"|id =="Gall_Forb_Edge_N_Y1") %>% dplyr::select(2:13)%>% as.matrix()
Mavb1_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id =="MAV_P_B1") %>% dplyr::select(2:13)%>% as.matrix()
Mavb2_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id =="MAV_P_B2") %>% dplyr::select(2:13)%>% as.matrix()
Mueller7.5 <- aa7.5 %>% filter(id == "Ea8R"|id =="Muller_Con_Int_N_Y1") %>% dplyr::select(2:13)%>% as.matrix()
P2_P_B2_7.5<- aa7.5 %>% filter(id == "Ea8R"|id =="P2_P_B2") %>% dplyr::select(2:13)%>% as.matrix()
P8_P_B1_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id =="P8_P_B1") %>% dplyr::select(2:13)%>% as.matrix()
P8_P_B2_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id =="P8_P_B2") %>% dplyr::select(2:13)%>% as.matrix()
P8_P_B3_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id =="P8_P_B3") %>% dplyr::select(2:13)%>% as.matrix()
Roll_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id == "Roll_Baugh_Org_Int_N_Y1") %>%dplyr::select(2:13)%>% as.matrix()
UCDFST02.252_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id =="UCDFST 02-252") %>% dplyr::select(2:13)%>% as.matrix()
UCDFST02.257_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id =="UCDFST 02-257") %>% dplyr::select(2:13)%>% as.matrix()
UCDFST02.286_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id =="UCDFST 02-286") %>% dplyr::select(2:13)%>% as.matrix()
UCDFST02.305_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id =="UCDFST 02-305") %>% dplyr::select(2:13)%>% as.matrix()
UCDFST02.329_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id =="UCDFST 02-329") %>% dplyr::select(2:13)%>% as.matrix()
UCDFST10.272_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id =="UCDFST 10-272") %>% dplyr::select(2:13)%>% as.matrix()
UCDFST67.77_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id =="UCDFST 67-77") %>% dplyr::select(2:13)%>% as.matrix()
UCDFST67.79_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id =="UCDFST 67-79") %>% dplyr::select(2:13)%>% as.matrix()
UCDFST68.105_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id =="UCDFST 68-105") %>% dplyr::select(2:13)%>% as.matrix()
UCDFST68.140_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id =="UCDFST 68-140") %>% dplyr::select(2:13)%>% as.matrix()
## Special cases 
B.mega7.5 <- df2 %>% filter(id == "B. mega"| id == "Ea8R", background == "7.5", No. != "297")%>% group_by(id)%>%summarise(asp = mean(ASP),gln_his = mean(GLN.HIS),glu = mean(GLU),gly = mean(GLY),ile = mean(ILE),leu = mean(LEU),lys = mean(LYS),
                                                                                                                          phe = mean(PHE),ser = mean(SER),thr = mean(THR),trp = mean(TRP),val_met = mean(VAL.MET)) %>% replace(is.na(.), 0)%>% dplyr::select(asp:val_met) %>% as.matrix()%>% abs()
EC137_7.5 <- df2 %>% filter(id == "EC_137"| id == "Ea8R", background == "7.5", No. != "210")%>% group_by(id)%>%summarise(asp = mean(ASP),gln_his = mean(GLN.HIS),glu = mean(GLU),gly = mean(GLY),ile = mean(ILE),leu = mean(LEU),lys = mean(LYS),
                                                                                                                         phe = mean(PHE),ser = mean(SER),thr = mean(THR),trp = mean(TRP),val_met = mean(VAL.MET)) %>% replace(is.na(.), 0)%>% dplyr::select(asp:val_met) %>% as.matrix()%>% abs()
AsaY1_7.5 <- df2 %>% filter(id == "Asa_Int_AS_Y1"| id == "Ea8R", background == "7.5")%>% group_by(id)%>%summarise(asp = mean(ASP),gln_his = mean(GLN.HIS),glu = mean(GLU),gly = mean(GLY),ile = mean(ILE),leu = mean(LEU),lys = mean(LYS),
                                                                                                                  phe = mean(PHE),ser = mean(SER),thr = mean(THR),trp = mean(TRP),val_met = mean(VAL.MET)) %>% replace(is.na(.), 0)%>% dplyr::select(asp:val_met) %>% as.matrix()%>% abs()
AsaY2_7.5 <- df2 %>% filter(id == "Asa_Edge_AS_Y2"| id == "Ea8R", background == "7.5")%>% group_by(id)%>%summarise(asp = mean(ASP),gln_his = mean(GLN.HIS),glu = mean(GLU),gly = mean(GLY),ile = mean(ILE),leu = mean(LEU),lys = mean(LYS),
                                                                                                                   phe = mean(PHE),ser = mean(SER),thr = mean(THR),trp = mean(TRP),val_met = mean(VAL.MET)) %>% replace(is.na(.), 0)%>% dplyr::select(asp:val_met) %>% as.matrix()%>% abs()


Eas15 <- aa15 %>% filter(id == "Ea8"|id =="Ea8R") %>% dplyr::select(2:13)%>% as.matrix()
AsaY2_15 <- aa15 %>% filter(id == "Ea8R"|id =="Asa_Edge_AS_Y2") %>% dplyr::select(2:13)%>% as.matrix()
Asay1_15 <- aa15 %>% filter(id == "Ea8R"|id =="Asa_Int_AS_Y1") %>% dplyr::select(2:13)%>% as.matrix()
B.mega15<- aa15 %>% filter(id == "Ea8R"|id =="B. mega") %>% dplyr::select(2:13)%>% as.matrix()
B.sub15<- aa15 %>% filter(id == "Ea8R"|id =="B. sub") %>% dplyr::select(2:13)%>% as.matrix()
Baugh15 <- aa15 %>% filter(id == "Ea8R"|id =="Baugh_Org_Int_N_Y1") %>% dplyr::select(2:13)%>% as.matrix()
Champs15 <- aa15 %>% filter(id == "Ea8R"|id == "Champs_P_B2") %>% dplyr::select(2:13)%>% as.matrix()
EC005_15 <- aa15 %>% filter(id == "Ea8R"|id =="EC_005") %>% dplyr::select(2:13)%>% as.matrix()
EC006_15 <- aa15 %>% filter(id == "Ea8R"|id == "EC_006") %>%dplyr::select(2:13)%>% as.matrix()
EC033_15 <- aa15 %>% filter(id == "Ea8R"|id =="EC_033") %>% dplyr::select(2:13)%>% as.matrix()
EC052_15 <- aa15 %>% filter(id == "Ea8R"|id == "EC_052") %>% dplyr::select(2:13)%>% as.matrix()
EC061_15 <- aa15 %>% filter(id == "Ea8R"|id =="EC_061") %>% dplyr::select(2:13)%>% as.matrix()
EC064_15 <- aa15 %>% filter(id == "Ea8R"|id == "EC_064") %>% dplyr::select(2:13)%>% as.matrix()
EC067_15 <- aa15 %>% filter(id == "Ea8R"|id =="EC_067") %>% dplyr::select(2:13)%>% as.matrix()
EC072_15 <- aa15 %>% filter(id == "Ea8R"|id =="EC_072") %>% dplyr::select(2:13)%>% as.matrix()
EC088_15 <- aa15 %>% filter(id == "Ea8R"|id == "EC_088") %>% dplyr::select(2:13)%>% as.matrix()
EC089_15 <- aa15 %>% filter(id == "Ea8R"|id =="EC_089") %>% dplyr::select(2:13)%>% as.matrix()
EC094_15 <- aa15 %>% filter(id == "Ea8R"|id =="EC_094") %>% dplyr::select(2:13)%>% as.matrix()
EC112_15 <- aa15 %>% filter(id == "Ea8R"|id =="EC_112") %>% dplyr::select(2:13)%>% as.matrix()
EC115_15 <- aa15 %>% filter(id == "Ea8R"|id =="EC_115") %>% dplyr::select(2:13)%>% as.matrix()
EC116_15 <- aa15 %>% filter(id == "Ea8R"|id =="EC_116") %>% dplyr::select(2:13)%>% as.matrix()
EC124_15 <- aa15 %>% filter(id == "Ea8R"|id == "EC_124") %>% dplyr::select(2:13)%>% as.matrix()
EC126_15 <- aa15 %>% filter(id == "Ea8R"|id == "EC_126") %>% dplyr::select(2:13)%>% as.matrix()
EC137_15 <- aa15 %>% filter(id == "Ea8R"|id == "EC_137") %>% dplyr::select(2:13)%>% as.matrix()
Gall15<- aa15 %>% filter(id == "Ea8R"|id =="Gall_Forb_Edge_N_Y1") %>% dplyr::select(2:13)%>% as.matrix()
Mavb1_15 <- aa15 %>% filter(id == "Ea8R"|id =="MAV_P_B1") %>% dplyr::select(2:13)%>% as.matrix()
Mavb2_15 <- aa15 %>% filter(id == "Ea8R"|id == "MAV_P_B2") %>% dplyr::select(2:13)%>% as.matrix()
Mueller15 <- aa15 %>% filter(id == "Ea8R"|id == "Muller_Con_Int_N_Y1") %>% dplyr::select(2:13)%>% as.matrix()
P2_P_B2_15<- aa15 %>% filter(id == "Ea8R"|id =="P2_P_B2") %>% dplyr::select(2:13)%>% as.matrix()
P8_P_B1_15 <- aa15 %>% filter(id == "Ea8R"|id == "P8_P_B1") %>% dplyr::select(2:13)%>% as.matrix()
P8_P_B2_15 <- aa15 %>% filter(id == "Ea8R"|id =="P8_P_B2") %>% dplyr::select(2:13)%>% as.matrix()
P8_P_B3_15 <- aa15 %>% filter(id == "Ea8R"|id == "P8_P_B3") %>% dplyr::select(2:13)%>% as.matrix()
Roll_15 <- aa15 %>% filter(id == "Ea8R"|id =="Roll_Baugh_Org_Int_N_Y1") %>% dplyr::select(2:13)%>% as.matrix()
UCDFST02.252_15 <- aa15 %>% filter(id == "Ea8R"|id =="UCDFST 02-252") %>% dplyr::select(2:13)%>% as.matrix()
UCDFST02.257_15 <- aa15 %>% filter(id == "Ea8R"|id =="UCDFST 02-257") %>% dplyr::select(2:13)%>% as.matrix()
UCDFST02.286_15 <- aa15 %>% filter(id == "Ea8R"|id =="UCDFST 02-286") %>% dplyr::select(2:13)%>% as.matrix()
UCDFST02.305_15 <- aa15 %>% filter(id == "Ea8R"|id =="UCDFST 02-305") %>% dplyr::select(2:13)%>% as.matrix()
UCDFST02.329_15 <- aa15 %>% filter(id == "Ea8R"|id == "UCDFST 02-329") %>% dplyr::select(2:13)%>% as.matrix()
UCDFST10.272_15 <- aa15 %>% filter(id == "Ea8R"|id == "UCDFST 10-272") %>% dplyr::select(2:13)%>% as.matrix()
UCDFST67.77_15 <- aa15 %>% filter(id == "Ea8R"|id =="UCDFST 67-77") %>% dplyr::select(2:13)%>% as.matrix()
UCDFST67.79_15 <- aa15 %>% filter(id == "Ea8R"|id =="UCDFST 67-79")%>% dplyr::select(2:13)%>% as.matrix()
UCDFST68.105_15 <- aa15 %>% filter(id == "Ea8R"|id == "UCDFST 68-105") %>% dplyr::select(2:13)%>% as.matrix()
UCDFST68.140_15 <- aa15 %>% filter(id == "Ea8R"|id =="UCDFST 68-140") %>% dplyr::select(2:13)%>% as.matrix()
EC84_7.5 <- aa7.5 %>% filter(id == "Ea8R"|id == "EC_84") %>% dplyr::select(2:13)%>% as.matrix()
EC84_15 <- aa15 %>% filter(id == "Ea8R"|id == "EC_84") %>% dplyr::select(2:13)%>% as.matrix()




PEas7.5 <- pianka(Eas7.5)
PEas15 <- pianka(Eas15)
PAsaY2_7.5 <- pianka(AsaY2_7.5)
PAsaY2_15 <- pianka(AsaY2_15)
PAsaY1_7.5 <- pianka(AsaY1_7.5)
PAsaY1_15 <- pianka(Asay1_15)
Pbmega7.5 <- pianka(B.mega7.5)
Pbmega15 <- pianka(B.mega15)
Pbsub7.5 <- pianka(B.sub7.5)
Pbsub15 <- pianka(B.sub15)
PBaugh7.5 <- pianka(Baugh7.5)
PBaugh15 <- pianka(Baugh15)
PChamps7.5 <- pianka(Champs7.5)
PChamps15 <- pianka(Champs15)
PEC005_7.5 <- pianka(EC005_7.5)
PEC005_15 <- pianka(EC005_15)
PEC006_7.5 <- pianka(EC006_7.5)
PEC006_15 <- pianka(EC006_15)
PEC033_7.5 <- pianka(EC033_7.5)
PEC033_15 <- pianka(EC033_15)
PEC052_7.5 <- pianka(EC052_7.5)
PEC052_15 <- pianka(EC052_15)
PEC061_7.5 <- pianka(EC061_7.5)
PEC061_15 <- pianka(EC061_15)
PEC064_7.5 <- pianka(EC064_7.5)
PEC064_15 <- pianka(EC064_15)
PEC067_7.5 <- pianka(EC067_7.5)
PEC067_15 <- pianka(EC067_15)
PEC072_7.5 <- pianka(EC072_7.5)
PEC072_15 <- pianka(EC072_15)
PEC088_7.5 <- pianka(EC088_7.5)
PEC088_15 <- pianka(EC088_15)
PEC089_7.5 <- pianka(EC089_7.5)
PEC089_15 <- pianka(EC089_15)
PEC094_7.5 <- pianka(EC094_7.5)
PEC094_15 <- pianka(EC094_15)
PEC112_7.5 <- pianka(EC112_7.5)
PEC112_15 <- pianka(EC112_15)
PEC115_7.5 <- pianka(EC115_7.5)
PEC115_15 <- pianka(EC115_15)
PEC116_7.5 <- pianka(EC116_7.5)
PEC116_15 <- pianka(EC116_15)
PEC124_7.5 <- pianka(EC124_7.5)
PEC124_15 <- pianka(EC124_15)
PEC126_7.5 <- pianka(EC126_7.5)
PEC126_15 <- pianka(EC126_15)
PEC137_7.5 <- pianka(EC137_7.5)
PEC137_15 <- pianka(EC137_15)
PEC84_7.5 <- pianka(EC84_7.5)
PEC84_15 <- pianka(EC84_15)
PGall7.5 <- pianka(Gall7.5)
PGall15 <- pianka(Gall15)
PMavB1_7.5 <- pianka(Mavb1_7.5)
PMavB1_15 <- pianka(Mavb1_15)
PMavB2_7.5 <- pianka(Mavb2_7.5)
PMavB2_15 <- pianka(Mavb2_15)
PMueller_7.5 <- pianka(Mueller7.5)
PMueller_15 <- pianka(Mueller15)
PP2_P_B2_7.5 <- pianka(P2_P_B2_7.5)
PP2_P_B2_15 <- pianka(P2_P_B2_15)
PP8_P_B1_7.5 <- pianka(P8_P_B1_7.5)
PP8_P_B1_15 <- pianka(P8_P_B1_15)
PP8_P_B2_7.5 <- pianka(P8_P_B2_7.5)
PP8_P_B2_15 <- pianka(P8_P_B2_15)
PP8_P_B3_7.5 <- pianka(P8_P_B3_7.5)
PP8_P_B3_15 <- pianka(P8_P_B3_15)
PRoll_7.5 <- pianka(Roll_7.5)
PRoll_15 <- pianka(Roll_15)
PUCDFST02.252_7.5 <- pianka(UCDFST02.252_7.5)
PUCDFST02.252_15 <- pianka(UCDFST02.252_15)
PUCDFST02.257_7.5 <- pianka(UCDFST02.257_7.5)
PUCDFST02.257_15 <- pianka(UCDFST02.257_15)
PUCDFST02.286_7.5 <- pianka(UCDFST02.286_7.5)
PUCDFST02.286_15 <- pianka(UCDFST02.286_15)
PUCDFST02.305_7.5 <- pianka(UCDFST02.305_7.5)
PUCDFST02.305_15 <- pianka(UCDFST02.305_15)
PUCDFST02.329_7.5 <- pianka(UCDFST02.329_7.5)
PUCDFST02.329_15 <- pianka(UCDFST02.329_15)
PUCDFST10.272_7.5 <- pianka(UCDFST10.272_7.5)
PUCDFST10.272_15 <- pianka(UCDFST10.272_15)
PUCDFST067.77_7.5 <- pianka(UCDFST67.77_7.5)
PUCDFST67.77_15 <- pianka(UCDFST67.77_15)
PUCDFST67.79_7.5 <- pianka(UCDFST67.79_7.5)
PUCDFST67.79_15 <- pianka(UCDFST67.79_15)
PUCDFST68.105_7.5 <- pianka(UCDFST68.105_7.5)
PUCDFST68.105_15 <- pianka(UCDFST68.105_15)
PUCDFST68.140_7.5 <- pianka(UCDFST68.140_7.5)
PUCDFST68.140_15 <- pianka(UCDFST68.140_15)




