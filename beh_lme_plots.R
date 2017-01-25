## analysize and plot estimated marginal means in Marta's placebo/NF data


library(nlme)
library(xtable)
library(lme4)
# load data
placeboR <- read_delim("~/Dropbox/data_projects/placebo_marta/placeboR.txt",
                       " ", escape_double = FALSE, trim_ws = TRUE)
# check variables
View(placeboR)


exp_m1 <- lme(stim_rating ~ stim*feedback_ratinglag, data = placeboR, ~stim + feedback_ratinglag|subject,na.action = na.omit)
anova(exp_m1)

# added AR1 -- not much different
exp_m2 <- lme(stim_rating ~ stim*feedback_ratinglag, data = placeboR, ~stim + feedback_ratinglag|subject,na.action = na.omit, correlation = corAR1())
anova(exp_m2)
anova(exp_m1,exp_m2)

f_m1 <- lme(feedback_rating ~ stim_rating + feedback + feedback_ratinglag, data = placeboR, ~stim_rating + feedback|subject,na.action = na.omit)
anova(f_m1)

f_m2 <- lme(feedback_rating ~ stim_rating + stim + feedback + feedback_ratinglag, data = placeboR, ~ stim + feedback|subject,na.action = na.omit)
anova(f_m2)

f_m3 <- lme(feedback_rating ~ stim_rating*stim + feedback, data = placeboR, ~ stim + feedback|subject,na.action = na.omit)
anova(f_m3)

# added AR1 -- not much different
f_m4 <- lme(feedback_rating ~ stim_rating*stim + feedback, data = placeboR, ~ stim + feedback|subject,na.action = na.omit, correlation = corAR1())
anova(f_m4)
anova(f_m3, f_m4)


# plot estimated marginal means
#install.packages("lsmeans")

library("lsmeans")
f_rg3 = ref.grid(f_m3)



ls_f_m3 <- lsmeans(f_m3,"feedback", at = list(feedback = c(-2,-1,1,2)))
ls_f_m4 <- lsmeans(f_m3,"stim_rating",by = "stim", at = list(stim_rating = c(-2,-1,0,1,2)))

print.xtable(f_m3, type = "HTML")

plot(ls_f_m3, type ~ feedback, horiz=F,ylab = "Mood rating", xlab = "Neurofeedback")

plot(ls_f_m4, type ~ stim_rating, horiz=F, ylab = "Mood ratings", xlab = "Expectancy")
strip=strip.custom(factor.levels=c("No infusion", "Placebo infusion"))
lsmip(ls_f_m4, stim_rating~stim)

# interaction plot
object = f_m3;
formula =  ~ stim_rating*stim;
type = "response"
lsmip(object, formula, type,
      pch = c(1,2,6,7,9,10,15:20),
      lty = 1, plotit = TRUE)

# Read in additional ind. differences and brain variables
bdf <- read_csv("~/Google Drive/skinner/Pecina_k/alex_behav_full_merge_with_betas.csv")
colnames(bdf)
bf_m1 <- lme(feedback_rating ~ stim_rating + feedback + feedback_ratinglag, data = bdf, ~stim_rating + feedback|subject,na.action = na.omit)
anova(bf_m1)

# import individual differences and betas data
df <- read_csv("~/Google Drive/skinner/Pecina_k/alex_behav_full_merge_with_betas.csv")
df$Sex <- as.factor(df$Sex)
df$subject<- as.factor(df$subject)
# rename ROIs
df$visual1 <- df$CuneusVisual
df$motor2 <- df$PostcentralMotor
df$cing3 <- df$LRCingulateMedFrontalGy
df$crbl4 <- df$Cerebellum
df$str5 <- df$LLentiformNucPutamen
df$striatum <- df$str5
df$postc6 <-df$RPostcentral
df$crbl7 <- df$Cerebellum2
df$thal8 <- df$LThalamus

ggplot(data = df, mapping = aes(df$cing3,df$str5)) +geom_point(mapping = aes(df$cing3,df$str5))


nm1 <- (lmer(feedback_rating ~ stim_rating + feedback +feedback_ratinglag + (1 +feedback + stim_rating |subject), df))
summary(nm1)
car::Anova(nm1)

nm2 <- (lmer(feedback_rating ~ stim_rating*stim + feedback +feedback_ratinglag + (1 +feedback + stim_rating |subject), df))
summary(nm2)
car::Anova(nm2)
anova(nm1,nm2)

ls_nm2 <- lsmeans(nm2,"stim_rating",by = "stim", at = list(stim_rating = c(-2,-1,0,1,2)))
plot(ls_nm2, type ~ stim_rating, horiz=F, ylab = "Mood ratings", xlab = "Expectancy")

nm3 <- (lmer(feedback_rating ~ stim_rating*stim + feedback +feedback_ratinglag + visual1*stim + (1 +feedback + stim_rating |subject), df))
summary(nm3)
car::Anova(nm3)

nm4 <- (lmer(feedback_rating ~ stim_rating*stim + feedback +feedback_ratinglag + motor2*stim + (1 +feedback + stim_rating |subject), df))
summary(nm4)
car::Anova(nm4)
anova(nm3,nm4)




nm5 <- (lmer(feedback_rating ~ stim_rating*stim + feedback +feedback_ratinglag + striatum*stim + (1 +feedback + stim_rating |subject), df))
summary(nm5)
car::Anova(nm5)
anova(nm4,nm5)

ls_nm5 <- lsmeans(nm5,"stim",by = "striatum", at = list(striatum = c(1,5,10)))

setwd("/Users/localadmin/Google Drive/skinner/Pecina_k/results/")

pdf("str_modulation.pdf", width=5, height=3.5)
plot(ls_nm5, type ~ stim, horiz=F, ylab = "Mood ratings", xlab = "infusion", main = "Striatal activity modulates effect of infusion")
dev.off()

# nm6 <- (lmer(feedback_rating ~ stim_rating*stim + feedback +feedback_ratinglag + visual1*stim  + (1 +feedback + stim_rating |subject), df))
# summary(nm6)
# car::Anova(nm6)
# anova(nm4,nm6)

# expectancy rating models

em1<- (lmer(stim_rating ~ stim + feedback_ratinglag*striatum  + (1 +stim  |subject), df))
summary(em1)
car::Anova(em1)

em2<- (lmer(stim_rating ~ stim + feedback_ratinglag*striatum + stim_diff  + (1 +stim  |subject), df))
summary(em2)
car::Anova(em2)

em3<- (lme(stim_rating ~ stim + feedback_ratinglag*striatum + stim_diff  + (1 +stim  |subject), df))
summary(em2)
car::Anova(em2)


ls_em2 <- lsmeans(em2,"feedback_ratinglag",by = "striatum", at = list(striatum = c(1,5,10),feedback_ratinglag = c(-2,0,2)))
pdf("str_modulation_feedback_rating_lag.pdf", width=5, height=3.5)
plot(ls_em2, type ~ stim, horiz=F, ylab = "Expectancy ratings", xlab = "Previous mood ratings", main = "Striatal activity modulates effect of earlier mood")
dev.off()

# include indiv. differences in feedback rating models

# rescale
numcols <- grep("^c\\.",names(df))

# df$apathy <- scale(df$ApathyEvaluationScale)

im1 <- (lmer(feedback_rating ~ stim_rating*stim + feedback +feedback_ratinglag + (1 +feedback + stim_rating |subject), df))
summary(im1)
car::Anova(im1)

im2 <- (lmer(feedback_rating ~ stim_rating*stim + feedback +feedback_ratinglag + df$DAS7*feedback + (1 +feedback + stim_rating |subject), df))
summary(im2)
car::Anova(im2)
anova(im1,im2)


im3 <- (lmer(feedback_rating ~ stim_rating*stim + feedback +feedback_ratinglag + df$YearsOfEd*feedback*stim + (1 +feedback + stim_rating |subject), df))
summary(im3)
car::Anova(im3)
anova(im1,im3)
ls_im3 <- lsmeans(im3,"stim",by = "df$YearsOfEd", at = list(df$YearsOfEd = c(8,13,18)))

ls_nm2 <- lsmeans(nm2,"stim_rating",by = "stim", at = list(stim_rating = c(-2,-1,0,1,2)))
plot(ls_nm2, type ~ stim_rating, horiz=F, ylab = "Mood ratings", xlab = "Expectancy")


im4 <- (lmer(feedback_rating ~ stim_rating*stim + feedback +feedback_ratinglag + df$striatum*stim + df$YearsOfEd*stim + (1 +feedback + stim_rating |subject), df))
summary(im4)
car::Anova(im4)
anova(im3,im4)



ls_nm5 <- lsmeans(nm5,"stim",by = "striatum", at = list(striatum = c(1,5,10)))




## check Apathy ratings
ggplot(data = df, mapping = aes(df$ApathyEvaluationScale,df$subject)) +geom_point(mapping = aes(df$ApathyEvaluationScale,df$subject))

ggplot(data = df, mapping = aes(df$ApathyEvaluationScale,df$Age)) +geom_point(mapping = aes(df$ApathyEvaluationScale,df$Age))

# troubleshoot non-convergence

tt <- getME(im1,"theta")
ll <- getME(im1,"lower")
min(tt[ll==0])

hist(df$Age)

#,correlation = corAR1()

# try nlme instead of lme4






