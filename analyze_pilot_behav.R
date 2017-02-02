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
b$exp_rating <- 2-b$exp_rating
b$exp_rt <- as.numeric(as.character(b$WillImpRt))



b$feed_rating <- as.numeric(as.character(b$ImprovedResp))
b$feed_rating <- 2 - b$feed_rating
b$feed_rt <- as.numeric(as.character(b$ImprovedRt))
b$feed_rt_log <- log(b$feed_rt)

b <- b %>% group_by(Participant) %>% mutate(exp_rating_lag = lag(exp_rating, n=1, order_by=TrialNum),
                                            feed_rating_lag = lag(feed_rating, n=1, order_by=TrialNum),
                                            exp_rt_lag = lag(exp_rt, n=1, order_by=TrialNum),

                                            feed_rt_lag = lag(feed_rt, n=1, order_by=TrialNum),
                                             feed_rating_lag = lag(feed_rating, n=1, order_by=TrialNum)
                                            )


b[, "exp_rating_lag":=c(NA, exp_rating[-.N]), by=Participant]


class(b$exp_rating)

exp_m1 <- glmer(exp_rating ~ drug*hi_reinf + (1|Participant), binomial(link = "logit"), data = b,na.action = na.omit)
car::Anova(exp_m1)


exp_m2 <- glmer(exp_rating ~ drug*hi_reinf + exp_rating_lag + (1|Participant), binomial(link = "logit"), data = b,na.action = na.omit)
car::Anova(exp_m2)


exp_rt_m1 <- lme(exp_rt ~ drug*hi_reinf + TrialNum, data = b, ~ 1|Participant, na.action = na.omit)
summary(exp_rt_m1)
anova(exp_rt_m1)

exp_rt_m1r <- lmer(exp_rt ~ drug*hi_reinf + TrialNum + (1|Participant), data = b, na.action = na.omit)
summary(exp_rt_m1r)
car::Anova(exp_rt_m1r)

exp_rt_m2r <- lmer(exp_rt ~ drug*hi_reinf + TrialNum + exp_rt_lag + (1|Participant), data = b, na.action = na.omit)
summary(exp_rt_m2r)
car::Anova(exp_rt_m2r)



f_m1 <-  glmer(feed_rating ~ drug*hi_reinf + (1|Participant), binomial(link = "logit"), data = b,na.action = na.omit)
car::Anova(f_m1)


f_m2 <-  glmer(feed_rating ~ drug + hi_reinf + (1|Participant), binomial(link = "logit"), data = b,na.action = na.omit)
summary(f_m2)
car::Anova(f_m2)

f_m2 <-  glmer(feed_rating ~ drug + hi_reinf +  (1|Participant), binomial(link = "logit"), data = b,na.action = na.omit)
summary(f_m2)
car::Anova(f_m2)


f_rt_m1 <- lme(feed_rt ~ drug*hi_reinf + TrialNum, data = b, ~ 1|Participant, na.action = na.omit)
summary(f_rt_m1)
anova(f_rt_m1)

f_rt_m2 <- lme(feed_rt ~ drug*hi_reinf + TrialNum, data = b, ~ 1|Participant, na.action = na.omit)
summary(f_rt_m2)
anova(f_rt_m2)


f_rt_m3 <- lme(feed_rt ~ drug + hi_reinf + TrialNum + feed_rt_lag, data = b, ~ 1|Participant, na.action = na.omit)
summary(f_rt_m3)
anova(f_rt_m3)

f_rt_m1r <- lmer(feed_rt ~ drug*hi_reinf + TrialNum + feed_rt_log + (1|Participant), data = b, na.action = na.omit)
summary(f_rt_m1r)
car::Anova(f_rt_m1r)


