
# R Script Graf_et_al_2017_R: R-Code to reproduce the results of
# Nadin Graf, Roman Bucher, Ralf B. Schäfer, Martin H. Entling

#----------------------------------------------------------------------------------------
# Contrasting short-term effects of aquatic subsidies on a recipient and simplified food web
#----------------------------------------------------------------------------------------
# submitted to Biology Letters

# The code has been written by:
# Nadin Graf
# University of Koblenz-Landau 
# Fortstrasse 7 
# 76829 Landau
# GERMANY
# graf-nadin@uni-landau.de

# Revised by Ralf B. Schäfer

# The code has been written to reproduce
# See 'graf_dataset.xls' for details

#############################
#      important notes      #
#############################

#Structure of the code:
# I # STATISTICS -------------------------------------------------------> line 45 to 156
# I.1 # indirect effects of aquatic subsidies on terrestrial plants-----> line 88 to 116
# I.2 # indirect effects of aquatic subsidies on terrestrial herbivores-> line 117 to 136
# I.3 # indirect effects of aquatic subsidies on terrestrial spiders----> line 137 to 158
# II # FIGURES ---------------------------------------------------------> line 159 to 288
# II.1 # Treatment vs. plants ------------------------------------------> line 168 to 238
# II.2 # Treatment vs. herbivores --------------------------------------> line 238 to 283
# II.3 # combine plots -------------------------------------------------> line 284 to 238
#----------------------------------------------------------------------------------------

# Set the path to your working directory here
pfad = "~/Literatur/Publications/2017/Nadin"
setwd(pfad)

#############################
#        STATISTICS         #
#############################

library(MASS)
library(doBy)
library(multcomp)
library(XLConnect)
# Function to assess standard error and confidence interval 
# conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE, conf.interval=.95) 
{  
  require(doBy)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) 
  {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # Collapse the data
  formula <- as.formula(paste(measurevar, paste(groupvars, collapse=" + "), sep=" ~ "))
  datac <- summaryBy(formula, data=data, FUN=c(length2,mean,sd), na.rm=na.rm)
  
  # Rename columns
  names(datac)[ names(datac) == paste(measurevar, ".mean",    sep="") ] <- measurevar
  names(datac)[ names(datac) == paste(measurevar, ".sd",      sep="") ] <- "sd"
  names(datac)[ names(datac) == paste(measurevar, ".length2", sep="") ] <- "N"
  # sqrt = ^1/2
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df = N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  # return value of interest
  return(datac)
}


#### indirect effects of aquatic subsdies on terrestrial plants ####
wb <- loadWorkbook("graf_dataset.xls")
dat <- readWorksheet(wb, sheet = 2, colTypes = c("character", "numeric", "character", "numeric", "numeric", "numeric", "numeric", "numeric"))

# herbivory: leaf area eaten by herbivores
# Leaf surface from start of experiment minus leaf area left at the end of the experiment

herb_1 <- dat$herb + 0.5/10^6  # a small value is for the glmm in the gamma family to run
# 0s cannot be calculated with a gamma-distribution
dat$Treatment <- factor(dat$Treatment)
pm1 <- glmmPQL(herb_1 ~ Treatment, random=~1|Replicate, family = Gamma(link ='log'), dat)
summary(glht(model = pm1, linfct = mcp(Treatment = "Tukey")))

# size of plants at the beginning of the experiment
pm2 <- glmmPQL(Start ~ Treatment, random=~1|Replicate, family = poisson, dat) 
summary(glht(model = pm2, linfct = mcp(Treatment = "Tukey")))

# size of plants at the end of the experiment
pm3 <- glmmPQL(end ~ Treatment, random=~1|Replicate, family = poisson, dat) 
summary(glht(model = pm3, linfct = mcp(Treatment = "Tukey")))

# growth of plants in size at the end of the experiment
pm4 <- glmmPQL(Growth ~ Treatment, random=~1|Replicate, family = poisson, dat)
summary(glht(model = pm4, linfct = mcp(Treatment = "Tukey")))

# dry mass of plants at the end of the experiment
pm5 <- glmmPQL(Mass~Treatment, random=~1|Replicate, family = poisson, dat)
summary(glht(model = pm5, linfct = mcp(Treatment = "Tukey")))



#### indirect effects of aquatic subsidies on terrestrial herbivores ####

dath<-read.table("Data Mesocosm_prey predator.csv",header=TRUE, sep=",")
datp<-dath[dath$prey == 'weevil', ]  # per mesocosm

# survival of Phyllobius sp.
hm1<-glm(Phyllobius~Treatment, family = poisson, datp) 
anova(hm1, test = "Chisq")

# survival of leafhoppers
hm2<-glm(Cicadina~Treatment, family = poisson, datp) 
anova(hm2, test = "Chisq")

# difference between the indirect effect of aquatic subsidies on herbivore response
dath$resp <- 100/dath$introduced*dath$survived
xy<-glmmPQL(resp~Treatment*prey,random=~1|Replicate,  family = negative.binomial(theta = 1), dath)
summary(xy)

#### direct effects of aquatic subsidies on terrestrial spiders ####

# survival of Pisaura mirablils
sm1<-glm(SurvPis~Treatment, family = poisson, datp)
anova(sm1, test = "Chisq")

# growth of mass of Pisaura mirablis             
sm2<-lm(PMG~Treatment, datp)
summary(sm2)
      
# width of prosoma of Pisaura mirablis      
sm3<-glm(PPWE~Treatment, family = poisson, datp)
anova(sm3, test = "Chisq")        

# Replacement of Tetragnatha sp.      
sm4<-glm(ReplTet~Treatment, family = poisson, datp)
anova(sm4, test = "Chisq")





#############################
#          FIGURES          #
#############################

library(ggplot2) 
library(cowplot) 

pd <- position_dodge(0) # manually specify position 

##### Treatment vs. plant ########

# Herbivory
herb.dat<- summarySE(dat, measurevar="herb", groupvars="Treatment")
herb<-ggplot(herb.dat, aes(Treatment, herb, cex=1)) +
  geom_errorbar(aes(ymin=herb-se, ymax=herb+se), width=0.1,
                size = 0.15, colour = "black", position = pd) +
  geom_line(linetype="blank", position = pd) +
  geom_point(size = 2)+
  scale_x_discrete(limits = c('control', 'aquatic', 'terrestrial'))+
  labs(y = expression('Herbivory ['~cm^2~'/day]'), x = '')+
  theme(axis.title.x = element_text(colour="black", size=13),
        axis.title.y = element_text(colour="black", size=13),
        axis.text.x  = element_blank(),
        axis.text.y  = element_text(colour="black", size=11),        
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill="white"), 
        panel.border = element_rect(colour="black", fill=NA, size=0.5),
        legend.position="none")+
  annotate("text",x=1,y=0.0185,label="a")+
  annotate("text",x=2,y=0.0185,label="b")+
  annotate("text",x=3,y=0.0185,label="b")
herb


Mass.dat<- summarySE(dat, measurevar="Mass", groupvars="Treatment")
Mass<-ggplot(Mass.dat, aes(Treatment, Mass, cex=1)) +
  geom_errorbar(aes(ymin=Mass-se, ymax=Mass+se), width=0.1,
                size = 0.15, colour = "black", position = pd) +
  geom_line(linetype="blank", position = pd) +
  geom_point(size = 2)+
  scale_x_discrete(limits = c("control", "aquatic", "terrestrial"))+
  labs(y = expression('Dry Mass [g]'), x = '')+
  theme(axis.title.x = element_text(colour="black", size=13),
        axis.title.y = element_text(colour="black", size=13),
        axis.text.x  = element_text(colour="black", size=11),
        axis.text.y  = element_text(colour="black", size=11),        
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill="white"), 
        panel.border = element_rect(colour="black", fill=NA, size=0.5),
        legend.position="none")+
  annotate("text",x=1,y=2.85,label="a b")+
  annotate("text",x=2,y=2.85,label="a")+
  annotate("text",x=3,y=2.85,label="b")
Mass


Growth.dat<- summarySE(dat, measurevar="Growth", groupvars="Treatment")
Growth<-ggplot(Growth.dat, aes(Treatment, Growth, cex=1)) +
  geom_errorbar(aes(ymin=Growth-se, ymax=Growth+se), width=0.1,
                size = 0.15, colour = "black", position = pd) +
  geom_line(linetype="blank", position = pd) +
  geom_point(size = 2)+
  scale_x_discrete(limits = c("control", "aquatic", "terrestrial"))+
  labs(y = expression('Growth [mm]'), x = '')+
  theme(axis.title.x = element_text(colour="black", size=13),
        axis.title.y = element_text(colour="black", size=13),
        axis.text.x  = element_blank(),
        axis.text.y  = element_text(colour="black", size=11),        
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill="white"), 
        panel.border = element_rect(colour="black", fill=NA, size=0.5),
        legend.position="none")+
  annotate("text",x=1,y=290,label="a")+
  annotate("text",x=2,y=290,label="a")+
  annotate("text",x=3,y=290,label="a")
Growth


##### Treatment vs. herbivore #####

Phyl.dat<- summarySE(datp, measurevar="Phyllobius", groupvars="Treatment")
Phyl<-ggplot(Phyl.dat, aes(Treatment, Phyllobius, cex=1)) +
  geom_errorbar(aes(ymin=Phyllobius-se, ymax=Phyllobius+se), width=0.1,
                size = 0.15, colour = "black", position = pd) +
  geom_line(linetype="blank", position = pd) +
  geom_point(size=2)+
  scale_x_discrete(limits = c("aquatic", "terrestrial"))+
  labs(y = expression(paste('Weevil')), x = '')+
  theme(axis.title.x = element_text(colour="black", size=13),
        axis.title.y = element_text(colour="black", size=13),
        axis.text.x  = element_text(colour="black", size=11),
        axis.text.y  = element_text(colour="black", size=11),        
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill="white"), 
        panel.border = element_rect(colour="black", fill=NA, size=0.5),
        legend.position="none")+
  scale_y_continuous(breaks=seq(0,4,by=1),labels=abs(seq(0,4,by=1)))+
  annotate("text",x=1,y=4,label="a")+
  annotate("text",x=2,y=4,label="a")
Phyl



Cic.dat<- summarySE(datp, measurevar="Cicadina", groupvars="Treatment")
Cic<-ggplot(Cic.dat, aes(Treatment, Cicadina, cex=1)) +
  geom_errorbar(aes(ymin=Cicadina-se, ymax=Cicadina+se), width=0.1,
                size = 0.15, colour = "black", position = pd) +
  geom_line(linetype="blank", position = pd) +
  geom_point(size=2)+
    scale_x_discrete(limits = c("aquatic", "terrestrial"))+
  labs(y = expression('Leafhopper'), x = '')+
  theme(axis.title.x = element_text(colour="black", size=13),
        axis.title.y = element_text(colour="black", size=13),
        axis.text.x  = element_blank(),
        axis.text.y  = element_text(colour="black", size=11),        
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill="white"), 
        panel.border = element_rect(colour="black", fill=NA, size=0.5),
        legend.position="none")+
  annotate("text",x=1,y=4.2,label="a")+
  annotate("text",x=2,y=4.2,label="b")
Cic


##### combine plots ########
res<-plot_grid(Cic, herb, Phyl, Growth, NULL, Mass, 
               labels = c("A", "C", "B", "D",' ' , "E"), 
               align='v', nrow=3, ncol = 2)
res
