library(lme4)
library(lsmeans)
library(data.table)
library(nlme)
library(xtable)
library(ggplot2)
# load data
library(readr)


b <- read_csv("~/Google Drive/skinner/Pecina_k/new_task_design/Behavioral/pilot_behav.txt")
b$WillImpResp <- gsub("NaN","NA",b$WillImpResp)
b$WillImpRt <- gsub("NaN","NA",b$WillImpRt)
b$ImprovedResp<- gsub("NaN","NA",b$ImprovedResp)
b$ImprovedRt<- gsub("NaN","NA",b$ImprovedRt)

View(b)

b$drug <- b$InfusionNum==1 | b$InfusionNum==2
b$hi_reinf <- b$InfusionNum==1 | b$InfusionNum==3
as.factor(b$Participant)
as.factor(b$InfusionNum)

# convert responses to numeric
b$exp_rating <- as.numeric(as.character(b$WillImpResp))
b$exp_rt <- as.numeric(as.character(b$WillImpRt))

b$feed_rating <- as.numeric(as.character(b$ImprovedResp))
b$feed_rt <- as.numeric(as.character(b$ImprovedRt))


class(b$exp_rating)



exp_m1 <- lme(exp_rating ~ drug + hi_reinf, data = b, ~ 1|Participant, na.action = na.omit)
anova(exp_m1)

exp_m2 <- lme(exp_rating ~ drug*hi_reinf, data = b, ~ 1|Participant, na.action = na.omit)
anova(exp_m2)
anova(exp_m1,exp_m2)

exp_m3 <- lme(exp_rating ~ drug*hi_reinf + hi_reinf*TrialNum, data = b, ~ 1|Participant, na.action = na.omit)
anova(exp_m3)
anova(exp_m2,exp_m3)

f_m1 <- lme(feed_rating ~ drug*hi_reinf*exp_rating, data = b, ~ 1|Participant, na.action = na.omit)
anova(f_m1)




f_m1 <- lme(feedback_rating ~ stim_rating + feedback + feedback_ratinglag, data = placeboR, ~stim_rating + feedback|subject,na.action = na.omit)
anova(f_m1)

