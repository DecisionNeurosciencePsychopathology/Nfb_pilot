#Load in libraries
library(readxl)
library(nlme)
library(xtable)
library(lme4)
library(readr)
library(ggplot2)
library(ggpubr)
library(scales)
library(lsmeans)
library(car)
library(stargazer)

#Data import
nfb_db_ROI_merge <- read_excel("C:/kod/placebo_neurofeedback/R/nfb_db_ROI_merge.xlsx")
df<-nfb_db_ROI_merge
colnames(df)

#First plot box plot for frequency of responses -- see creating box_plots.R


#Run the model codes feom LMER_code_NFB.R to get the ls_nmX models to plot
if(!require(psych)){install.packages("psych")}
if(!require(ordinal)){install.packages("ordinal")}
if(!require(car)){install.packages("car")}
if(!require(RVAideMemoire)){install.packages("RVAideMemoire")}
if(!require(lsmeans)){install.packages("lsmeans")}
if(!require(multcompView)){install.packages("multcompView")}
library(psych)
library(ordinal)
library(car)
library(multcompView)




ls_mean_return_obj <- function(ls_in){
  CLD = cld(ls_in,
            alpha=0.05,
            Letters=letters,
            adjust="tukey")
  return(CLD)
}

pd = position_dodge(0.4) 


# CLD = cld(ls_nm1,
#           alpha=0.05,
#           Letters=letters,
#           adjust="tukey")

#####STIM RATINGS####
CLD_nm1a <- ls_mean_return_obj(ls_nm1)
CLD_nm1b <- ls_mean_return_obj(ls_nm1)



####FEEDBACK RATINGS####
CLD_nm2_a <- ls_mean_return_obj(ls_nm2)
CLD_nm2_b <- ls_mean_return_obj(ls_nm2)
CLD_nm2_c <- ls_mean_return_obj(ls_nm2)


####PCA####
CLD_nm4_a <- ls_mean_return_obj(ls_nm4) #Stim*PCA
CLD_nm4_b <- ls_mean_return_obj(ls_nm4) #Feedback*PCA


##STIM PLotting ##

#Plot 1 - with stim
ggplot(CLD_nm1a,
       aes(x = stim,
           y = lsmean)) +
  geom_point(shape = 15,
             size = 4,
             position = pd,
             color="blue") +
  geom_errorbar(aes(ymin = lower.CL,
                    ymax = upper.CL),
                width = 0.2,
                size = 1.0,
                position = pd,
                color="blue") +
  theme_bw() +
  
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  
  ylab("Expectancy ratings") + xlab("Infusion Cue") + scale_x_discrete(labels=c("No infusion", "Infusion")) +
  ggtitle("GGtitle") +
  labs(caption=paste0("\nYou can put a caption here ",
                      "\nand even here!"), hjust=0.5)



#Plot 2 - with feedback_ratinglag
p <-ggplot(CLD_nm1b,
       aes(x = feedback_ratinglag,
           y = lsmean)) +
  geom_point(shape = 15,
             size = 4,
             position = pd,
             color="blue") +
  geom_errorbar(aes(ymin = lower.CL,
                    ymax = upper.CL),
                width = 0.2,
                size = 1.0,
                position = pd,
                color="blue") +
  theme_bw() +
  
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  
  ylab("Expectancy ratings") + xlab("Feedback ratings") + 
  ggtitle("GGtitle") +
  labs(caption=paste0("\nYou can put a caption here ",
                      "\nand even here!"), hjust=0.5)


p + scale_x_continuous(breaks=c(-2:2),labels=c("Worse", "Slightly Worse", "No Change", "Slightly Better", "Better")) 

##Feedback PLotting ##

#Plot 1 - With stim
ggplot(CLD_nm2_a,
       aes(x = stim,
           y = lsmean)) +
  geom_point(shape = 15,
             size = 4,
             position = pd,
             color="blue") +
  geom_errorbar(aes(ymin = lower.CL,
                    ymax = upper.CL),
                width = 0.2,
                size = 1.0,
                position = pd,
                color="blue") +
  theme_bw() +
  
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  
  ylab("Feedback ratings") + xlab("Infusion Cue") + scale_x_discrete(labels=c("No infusion", "Infusion")) +
  ggtitle("GGtitle") +
  labs(caption=paste0("\nYou can put a caption here ",
                      "\nand even here!"), hjust=0.5)





#Plot 2 - With stim rating
ggplot(CLD_nm2_b,
       aes(x = stim_rating,
           y = lsmean)) +
  geom_point(shape = 15,
             size = 4,
             position = pd,
             color="blue") +
  geom_errorbar(aes(ymin = lower.CL,
                    ymax = upper.CL),
                width = 0.2,
                size = 1.0,
                position = pd,
                color="blue") +
  theme_bw() +
  
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  
  ylab("Feedback ratings") + xlab("Expectancy Ratings") + 
  scale_x_continuous(breaks=c(-2:2),labels=c("Worse", "Slightly Worse", "No Change", "Slightly Better", "Better")) +
  ggtitle("GGtitle") +
  labs(caption=paste0("\nYou can put a caption here ",
                      "\nand even here!"), hjust=0.5)





#Plot 3 - With feedback
CLD_nm2_c$.group=gsub(" ", "", CLD_nm2_c$.group)
ggplot(CLD_nm2_c,
       aes(x = feedback,
           y = lsmean)) +
  geom_point(shape = 15,
             size = 4,
             position = pd,
             color="blue") +
  geom_errorbar(aes(ymin = lower.CL,
                    ymax = upper.CL),
                width = 0.2,
                size = 1.0,
                position = pd,
                color="blue") +
  theme_bw() +
  
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  
  ylab("Feedback ratings") + xlab("Feedback") + 
  scale_x_continuous(breaks=c(-2:-1,1:2),labels=c("100% Negative", "50% Negative", "50% Positive", "100% Positive")) +
  ggtitle("GGtitle") +
  labs(caption=paste0("\nYou can put a caption here ",
                      "\nand even here!"), hjust=0.5)
  


#PCA plotting
CLD_nm4_a$PCA <- as.factor(CLD_nm4_a$PCA)
p<-ggplot(CLD_nm4_a,
       aes(x = stim,
           y = lsmean,
           colour = factor(PCA,labels=c("Decreased activation", "No change in activation","Increased activation")))) +
  geom_point(shape = 15,
             size = 4,
             position = pd) +
  geom_errorbar(aes(ymin = lower.CL,
                    ymax = upper.CL),
                width = 0.2,
                size = 1.0,
                position = pd) +
  theme_bw() +
  
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  
  ylab("Mood ratings") + xlab("Feedback") + 
  scale_x_discrete(labels=c("No infusion", "Infusion")) +
  ggtitle("GGtitle") +
  labs(caption=paste0("\nYou can put a caption here ",
                      "\nand even here!"), hjust=0.5)

p$labels$colour <- "Corticostriatothalamic\nnetwork"
p










CLD_nm4_b$PCA <- as.factor(CLD_nm4_b$PCA)
#Replace the values for plotting - hacky but works
CLD_nm4_b$feedback[CLD_nm4_b$feedback == -2] <-5
CLD_nm4_b$feedback[CLD_nm4_b$feedback == -1] <-6
CLD_nm4_b$feedback[CLD_nm4_b$feedback == 1] <-7
CLD_nm4_b$feedback[CLD_nm4_b$feedback == 2] <-8

p<-ggplot(CLD_nm4_b,
       aes(x = feedback,
           y = lsmean,
           colour = factor(PCA,labels=c("Decreased activation", "No change in activation","Increased activation")))) +
  geom_point(shape = 15,
             size = 4,
             position = pd) +
  geom_errorbar(aes(ymin = lower.CL,
                    ymax = upper.CL),
                width = 0.2,
                size = 1.0,
                position = pd) +
  theme_bw() +
  
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  
  ylab("Mood ratings") + xlab("Feedback") + 
  scale_x_continuous(labels=c("100% Negative", "50% Negative", "50% Positive", "100% Positive")) +
  ggtitle("GGtitle") +
  labs(caption=paste0("\nYou can put a caption here ",
                      "\nand even here!"), hjust=0.5)
p$labels$colour <- "Corticostriatothalamic\nnetwork"
p



