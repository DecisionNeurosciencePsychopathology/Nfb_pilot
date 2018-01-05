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
nfb_db_ROI_merge <- read_excel("C:/kod/placebo_neurofeedback/Paper/nfb_db_ROI_merge.xlsx")
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
if(!require(extrafont)){install.packages("extrafont")}
library(psych)
library(ordinal)
library(car)
library(multcompView)
library(extrafont)

#fonts
font_import()
loadfonts(device="win")
#fonts()


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



####Mood ratings####
CLD_nm2_a <- ls_mean_return_obj(ls_nm2b)
CLD_nm2_b <- ls_mean_return_obj(ls_nm2b)
CLD_nm2_c <- ls_mean_return_obj(ls_nm2b)
CLD_nm2_d <- ls_mean_return_obj(ls_nm2b)


####PCA####
CLD_nm4_a <- ls_mean_return_obj(ls_nm4) #Stim*PCA
CLD_nm4_b <- ls_mean_return_obj(ls_nm4) #Feedback*PCA

###Endorphins###
CLD_nm5_a <- ls_mean_return_obj(ls_nm5) #Stim Endoph
CLD_nm6_a <- ls_mean_return_obj(ls_nm6) #Feedback Endoph

###HAMILTON###
CLD_nm7_a <- ls_mean_return_obj(ls_nm7a) #Infuction Exp Ratings Hamilton
CLD_nm8_a <- ls_mean_return_obj(ls_nm8a) #Feedback Mood ratings hamilton
CLD_nm8_b <- ls_mean_return_obj(ls_nm8b) #Mood ratings vs Exp Ratings Hamilton

###Education###
CLD_nm9_a <- ls_mean_return_obj(ls_nm9) #Infusion Expectancy Ratings Edu
CLD_nm10_a <- ls_mean_return_obj(ls_nm10a) #Feedback Mood ratings Edu
CLD_nm10_b <- ls_mean_return_obj(ls_nm10) #Infusion Mood ratings Edu



#Bar vars
bar_top_width = 0.25
bar_col_width = 4
point_size = 7
text_size = 55
my_colors = c("dodgerblue1", "green2", "firebrick2")
hgt = 12
wdth= 25

##STIM PLotting ##

#Plot 1 - with stim
tiff('figure_1_stimVsexp.tiff', units="in", width=wdth, height=hgt, res=300)

ggplot(CLD_nm1a,
       aes(x = stim,
           y = lsmean)) +
  geom_point(shape = 15,
             size = point_size,
             position = pd,
             color="blue") +
  geom_errorbar(aes(ymin = lower.CL,
                    ymax = upper.CL),
                width = bar_top_width,
                size = bar_col_width,
                position = pd,
                color="blue") +
  theme_bw() +
  
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0),
        text = element_text(family = "Times New Roman", size=text_size),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour="black"), 
        axis.text.y = element_text(colour="black")) +
  
  ylab("Expectancy ratings\n") + xlab("\nPlacebo Infusion Cue") + scale_x_discrete(labels=c("No infusion", "Infusion"))

dev.off()



#Plot 2 - with feedback_ratinglag
tiff('figure_2_feedVsexp.tiff', units="in", width=30, height=hgt, res=300)
p <-ggplot(CLD_nm1b,
       aes(x = feedback_ratinglag,
           y = lsmean)) +
  geom_point(shape = 15,
             size = point_size,
             position = pd,
             color="blue") +
  geom_errorbar(aes(ymin = lower.CL,
                    ymax = upper.CL),
                width = bar_top_width,
                size = bar_col_width,
                position = pd,
                color="blue") +
  theme_bw() +
  
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0),
        text = element_text(family = "Times New Roman", size=text_size),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour="black"), 
        axis.text.y = element_text(colour="black")) +
  
  ylab("Expectancy ratings\n") + xlab("\nMood ratings")


p + scale_x_continuous(breaks=c(-3:3),labels=c("Much Worse", "Worse", "Slightly Worse", "No Change", "Slightly Better", "Better", "Much Better")) 

dev.off()

##Feedback PLotting ##

#Plot 1 - With stim
tiff('figure_3_stimVsFeedRat.tiff', units="in", width=wdth, height=hgt, res=300)

ggplot(CLD_nm2_a,
       aes(x = stim,
           y = lsmean)) +
  geom_point(shape = 15,
             size = point_size,
             position = pd,
             color="blue") +
  geom_errorbar(aes(ymin = lower.CL,
                    ymax = upper.CL),
                width = bar_top_width,
                size = bar_col_width,
                position = pd,
                color="blue") +
  theme_bw() +
  
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0),
        text = element_text(family = "Times New Roman", size=text_size),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour="black"), 
        axis.text.y = element_text(colour="black")) +
  
  ylab("Mood ratings\n") + xlab("\nPlacebo Infusion Cue") + scale_x_discrete(labels=c("No infusion", "Infusion"))

dev.off()



#Plot 2 - With stim rating
tiff('figure_4_stimRatVsFeedRat.tiff', units="in", width=30, height=hgt, res=300)
ggplot(CLD_nm2_b,
       aes(x = stim_rating,
           y = lsmean)) +
  geom_point(shape = 15,
             size = point_size,
             position = pd,
             color="blue") +
  geom_errorbar(aes(ymin = lower.CL,
                    ymax = upper.CL),
                width = bar_top_width,
                size = bar_col_width,
                position = pd,
                color="blue") +
  theme_bw() +
  
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0),
        text = element_text(family = "Times New Roman", size=text_size),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour="black"), 
        axis.text.y = element_text(colour="black")) +
  
  ylab("Mood ratings\n") + xlab("\nExpectancy Ratings") + 
  scale_x_continuous(breaks=c(-2:2),labels=c("Worse", "Slightly Worse", "No Change", "Slightly Better", "Better"))

dev.off()



#Plot 3 - With feedback
#Replace the values for plotting - hacky but works
#CLD_nm2_c$feedback[CLD_nm2_c$feedback == -2] <-5
#CLD_nm2_c$feedback[CLD_nm2_c$feedback == -1] <-6
#CLD_nm2_c$feedback[CLD_nm2_c$feedback == 1] <-7
#CLD_nm2_c$feedback[CLD_nm2_c$feedback == 2] <-8
tiff('figure_5_feedVsfeedRat.tiff', units="in", width=wdth, height=hgt, res=300)
CLD_nm2_c$.group=gsub(" ", "", CLD_nm2_c$.group)
ggplot(CLD_nm2_c,
       aes(x = feedback,
           y = lsmean)) +
  geom_point(shape = 15,
             size = point_size,
             position = pd,
             color="blue") +
  geom_errorbar(aes(ymin = lower.CL,
                    ymax = upper.CL),
                width = bar_top_width,
                size = bar_col_width,
                position = pd,
                color="blue") +
  theme_bw() +
  
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0),
        text = element_text(family = "Times New Roman", size=text_size),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour="black"), 
        axis.text.y = element_text(colour="black")) +
  
  ylab("Mood ratings\n") + xlab("\nFeedback") + 
  scale_x_continuous(labels=c("100% Negative", "50% Negative", "", "50% Positive", "100% Positive"))
  
dev.off()


#Plot 4 - X = feedback Y = feedback rating by stim rating (much like PCA plots below)
CLD_nm2_d$stim_rating <- as.factor(CLD_nm2_d$stim_rating)
#CLD_nm2_d$feedback[CLD_nm2_d$feedback == -2] <-5
#CLD_nm2_d$feedback[CLD_nm2_d$feedback == -1] <-6
#CLD_nm2_d$feedback[CLD_nm2_d$feedback == 1] <-7
#CLD_nm2_d$feedback[CLD_nm2_d$feedback == 2] <-8
tiff('figure_6_feedVsfeedRatByStimrating.tiff', units="in", width=wdth+5, height=12, res=300)
p<-ggplot(CLD_nm2_d,
          aes(x = feedback,
              y = lsmean,
              colour = factor(stim_rating,labels=c("Worse", "No change","Better")))) +
  geom_point(shape = 15,
             size = point_size,
             position = pd) +
  geom_errorbar(aes(ymin = lower.CL,
                    ymax = upper.CL),
                width = bar_top_width,
                size = bar_col_width,
                position = pd) +
  theme_bw() +
  
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0),
        text = element_text(family = "Times New Roman", size=text_size)) +
  
  ylab("Mood ratings\n") + xlab("\nFeedback") + 
  scale_x_continuous(labels=c("100% Negative", "50% Negative", "", "50% Positive", "100% Positive"))

p$labels$colour <- "Expectancy ratings"
p<-p + scale_color_manual(values=my_colors)
p + theme(axis.text=element_text(size=text_size),
          axis.title=element_text(size=text_size,face="bold"),
          legend.text = element_text(size=text_size),
          title = element_text(size = text_size),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"),
          legend.position = c(0.2, 0.75),
          legend.key.height = unit(0.5, 'inches')) 

dev.off()



#PCA plotting
CLD_nm4_a$PCA <- as.factor(CLD_nm4_a$PCA)
tiff('figure_7_PCA_stim.tiff', units="in", width=wdth+5, height=hgt+5, res=300)
p<-ggplot(CLD_nm4_a,
       aes(x = stim,
           y = lsmean,
           colour = factor(PCA,labels=c("Low Response", "Mean Response","High Response")))) +
  geom_point(shape = 15,
             size = point_size,
             position = pd) +
  geom_errorbar(aes(ymin = lower.CL,
                    ymax = upper.CL),
                width = bar_top_width,
                size = bar_col_width,
                position = pd) +
  theme_bw() +
  
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0),
        text = element_text(family = "Times New Roman", size=text_size)) +
  
  ylab("Mood ratings\n") + xlab("\nInfusion") + 
  scale_x_discrete(labels=c("No infusion", "Infusion"))

p$labels$colour <- "Corticostriatothalamic\nloop"
p<-p + scale_color_manual(values=my_colors)
p + theme(axis.text=element_text(size=text_size),
          axis.title=element_text(size=text_size,face="bold"),
          legend.text = element_text(size=text_size),
          title = element_text(size = text_size),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"),
          legend.position = "none") 

dev.off()



CLD_nm4_b$PCA <- as.factor(CLD_nm4_b$PCA)
#Replace the values for plotting - hacky but works
#CLD_nm4_b$feedback[CLD_nm4_b$feedback == -2] <-5
#CLD_nm4_b$feedback[CLD_nm4_b$feedback == -1] <-6
#CLD_nm4_b$feedback[CLD_nm4_b$feedback == 1] <-7
#CLD_nm4_b$feedback[CLD_nm4_b$feedback == 2] <-8
tiff('figure_7_PCA_feedback.tiff', units="in", width=wdth+5, height=hgt+5, res=300)
p<-ggplot(CLD_nm4_b,
       aes(x = feedback,
           y = lsmean,
           colour = factor(PCA,labels=c("Low Response", "Mean Response","High Response")))) +
  geom_point(shape = 15,
             size = point_size,
             position = pd) +
  geom_errorbar(aes(ymin = lower.CL,
                    ymax = upper.CL),
                width = bar_top_width,
                size = bar_col_width,
                position = pd) +
  theme_bw() +
  
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0),
        text = element_text(family = "Times New Roman", size=text_size)) +
  
  ylab("Mood ratings\n") + xlab("\nFeedback") + 
  scale_x_continuous(labels=c("100% Negative", "50% Negative", "", "50% Positive", "100% Positive"))

p$labels$colour <- "Corticostriatothalamic loop"
p<-p + scale_color_manual(values=my_colors)
p + theme(axis.text=element_text(size=text_size),
        axis.title=element_text(size=text_size,face="bold"),
        legend.text = element_text(size=text_size),
        title = element_text(size = text_size),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour="black"), 
        axis.text.y = element_text(colour="black"),
        #legend.key.height = unit(0.6, 'inches')
        legend.position = "none") 

dev.off()
    


ggplot_nf_data <- function(df,x_data,y_data,factor_data, factor_labels,text_size,y_title,x_title,legend_title){
  p <- ggplot(df,
            aes(x = x_data,
                y = y_data,
                colour = factor(factor_data,labels=factor_labels))) +
    geom_point(shape = 15,
               size = point_size,
               position = pd) +
    geom_errorbar(aes(ymin = lower.CL,
                      ymax = upper.CL),
                  width = bar_top_width,
                  size = bar_col_width,
                  position = pd) +
    theme_bw() +
    
    theme(axis.title   = element_text(face = "bold"),
          axis.text    = element_text(face = "bold"),
          plot.caption = element_text(hjust = 0),
          text = element_text(family = "Times New Roman", size=text_size)) +
    ylab(paste(y_title,"\n")) + xlab(paste("\n",x_title))
  
    p$labels$colour <- legend_title
  p + theme(axis.text=element_text(size=text_size),
            axis.title=element_text(size=text_size,face="bold"),
            legend.text = element_text(size=text_size),
            title = element_text(size = text_size),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(colour="black"), 
            axis.text.y = element_text(colour="black"),
            legend.key.height = unit(0.7, 'inches'))
}



#New plots ENDORPH
#Stim
tiff('figure_8_Endorf_Stim.tiff', units="in", width=wdth+5, height=hgt+8, res=300)
p<-ggplot_nf_data(CLD_nm5_a,
            CLD_nm5_a$stim,
            CLD_nm5_a$lsmean,
            CLD_nm5_a$endorf2_1,
            c("Decreases (-150 pg/ml)", "No change","Increases (150 pg/ml)"),
            text_size,
            "Expectancy Ratings",
            "Infusion Cue",
            "Change in circulating β-Endorphins (pg/ml)")
p<-p + scale_color_manual(values=my_colors) + theme(legend.position = "none")
p + scale_x_discrete(labels=c("No infusion", "Infusion")) 
dev.off()

#Feedback
tiff('figure_8_Endorf_Feed.tiff', units="in", width=wdth+5, height=hgt+8, res=300)
#CLD_nm6_a$feedback[CLD_nm6_a$feedback == -2] <-5
#CLD_nm6_a$feedback[CLD_nm6_a$feedback == -1] <-6
#CLD_nm6_a$feedback[CLD_nm6_a$feedback == 1] <-7
#CLD_nm6_a$feedback[CLD_nm6_a$feedback == 2] <-8
p <- ggplot_nf_data(CLD_nm6_a,
               CLD_nm6_a$feedback,
               CLD_nm6_a$lsmean,
               CLD_nm6_a$endorf2_1,
               c("Decreases (-150 pg/ml)", "No change","Increases (150 pg/ml)"),
               text_size,
               "Mood ratings",
               "Feedback",
               "Change in circulating β-Endorphins (pg/ml)")
p<-p + scale_color_manual(values=my_colors) + theme(legend.position = "none")
p + scale_x_continuous(labels=c("100% Negative", "50% Negative", "", "50% Positive", "100% Positive"))
dev.off()

###HAMILTON###
#New plots Hamiliton
#Stim
tiff('figure_9_HAM_stim.tiff', units="in", width=wdth, height=hgt+7, res=300)
p<-ggplot_nf_data(CLD_nm7_a,
                  CLD_nm7_a$stim,
                  CLD_nm7_a$lsmean,
                  CLD_nm7_a$SNAithHamilton,
                  c("No anhedonia (0)", "Low anhedonia (5)","High anhedonia (10)"),
                  text_size,
                  "Expectancy Ratings",
                  "Infusion Cue",
                  "Anhedonia Scores (SHAPS)")
p<-p + scale_color_manual(values=my_colors) + theme(legend.position = "none")
p + scale_x_discrete(labels=c("No infusion", "Infusion"))
dev.off()

#Feedback
#CLD_nm8_a$feedback[CLD_nm8_a$feedback == -2] <-5
#CLD_nm8_a$feedback[CLD_nm8_a$feedback == -1] <-6
#CLD_nm8_a$feedback[CLD_nm8_a$feedback == 1] <-7
#CLD_nm8_a$feedback[CLD_nm8_a$feedback == 2] <-8
tiff('figure_10_HAM_feed.tiff', units="in", width=wdth+5, height=hgt+8, res=300)
p <- ggplot_nf_data(CLD_nm8_a,
                    CLD_nm8_a$feedback,
                    CLD_nm8_a$lsmean,
                    CLD_nm8_a$SNAithHamilton,
                    c("No anhedonia (0)", "Low anhedonia (5)","High anhedonia (10)"),
                    text_size,
                    "Mood ratings",
                    "Feedback",
                    "Anhedonia Scores (SHAPS)")
p<-p + scale_color_manual(values=my_colors) + theme(legend.position = "none")
p + scale_x_continuous(labels=c("100% Negative", "50% Negative", "", "50% Positive", "100% Positive"))
dev.off()

#Stim rating vs feedback rating
tiff('figure_11_HAM_stimvsfeedRat.tiff', units="in", width=wdth+5, height=hgt+8, res=300)
p <- ggplot_nf_data(CLD_nm8_b,
                    CLD_nm8_b$stim_rating,
                    CLD_nm8_b$lsmean,
                    CLD_nm8_b$SNAithHamilton,
                    c("No anhedonia (0)", "Low anhedonia (5)","High anhedonia (10)"),
                    text_size,
                    "Mood ratings",
                    "Expectancy Rating",
                    "Anhedonia Scores (SHAPS)")
p<-p + scale_color_manual(values=my_colors) + theme(legend.position = "none")
p + scale_x_continuous(labels=c("Worse", "Slightly Worse", "No Change", "Slightly Better", "Better"))
dev.off()

###Education###
#New education graphs
tiff('figure_12_EDU_stim.tiff', units="in", width=wdth, height=hgt, res=300)
p<-ggplot_nf_data(CLD_nm9_a,
                  CLD_nm9_a$stim,
                  CLD_nm9_a$lsmean,
                  CLD_nm9_a$YearsOfEd,
                  c("9 years", "13 years","18 years"),
                  text_size,
                  "Expectancy Ratings",
                  "Infusion Cue",
                  "Education")
p<-p + scale_color_manual(values=my_colors) + theme(legend.position = "none")
p + scale_x_discrete(labels=c("No infusion", "Infusion"))
dev.off()

#Feedback
#CLD_nm10_a$feedback[CLD_nm10_a$feedback == -2] <-5
#CLD_nm10_a$feedback[CLD_nm10_a$feedback == -1] <-6
#CLD_nm10_a$feedback[CLD_nm10_a$feedback == 1] <-7
#CLD_nm10_a$feedback[CLD_nm10_a$feedback == 2] <-8
tiff('figure_13_EDU_feed.tiff', units="in", width=wdth+5, height=hgt+7, res=300)
p <- ggplot_nf_data(CLD_nm10_a,
                    CLD_nm10_a$feedback,
                    CLD_nm10_a$lsmean,
                    CLD_nm10_a$YearsOfEd,
                    c("9 years", "13 years","18 years"),
                    text_size,
                    "Mood ratings",
                    "Feedback",
                    "Education")
p<-p + scale_color_manual(values=my_colors) + theme(legend.position = "none")
p + scale_x_continuous(labels=c("100% Negative", "50% Negative", "", "50% Positive", "100% Positive"))
dev.off()


#Infusion Mood ratings
tiff('figure_14_EDU_stim_feedRat.tiff', units="in", width=wdth+5, height=hgt+7, res=300)
p<-ggplot_nf_data(CLD_nm10_b,
                  CLD_nm10_b$stim,
                  CLD_nm10_b$lsmean,
                  CLD_nm10_b$YearsOfEd,
                  c("9 years", "13 years","18 years"),
                  text_size,
                  "Mood ratings",
                  "Infusion Cue",
                  "Education")
p<-p + scale_color_manual(values=my_colors) + theme(legend.position = "none")
p + scale_x_discrete(labels=c("No infusion", "Infusion"))
dev.off()
