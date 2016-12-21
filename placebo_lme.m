%% Multi-level models of reinforcement predicting subject ratings on Marta's placebo/neurofeedback data

%% Read in data
cd('/Users/localadmin/Dropbox/data_projects/placebo_marta/');
load('p');
% read_nf_data
% cs = [];
% us = [];
% for i = 1:length(feedback)
% cs(i) = strcmpi(char(stim(i)),'Infusion');
% if strcmpi(char(feedback(i)),'100% Pos Feedback');
% us(i) = 2;
% elseif strcmpi(char(feedback(i)),'50% Pos Feedback');
% us(i) = 1;
% elseif strcmpi(char(feedback(i)),'50% Neg Feedback');
% us(i) = -1;
% elseif strcmpi(char(feedback(i)),'100% Neg Feedback');
% us(i) = -2;
% end
% end
% cs = logical(cs);
% % nf = table(us',cs',stim,feedback,feedback_ratings,stim_ratings);
% 
% %% write to table for LME analysis
% clear p; p = table;
% trials = 72;
% start = 1;
% for sub = 1:size(feedback_ratings,2)
%     p.feedback_rating(start:start+trials-1,1) = feedback_ratings(:,sub);
%     p.subject(start:start+trials-1,1) = sub*ones(trials,1);
%     p.stim_rating(start:start+trials-1,1) = stim_ratings(:,sub);
%     p.trial(start:start+trials-1,1) = [1:72]';
%     p.stim(start:start+trials-1,1)  = cs';
%     p.feedback(start:start+trials-1,1) = us';
%     start = start + trials;
% end
% p.subject = nominal(p.subject);

% Sanity check
figure(1); boxplot(p.feedback_rating,p.feedback); xlabel('Feedback'); ylabel('Rating');
figure(2); boxplot(p.stim_rating, p.stim); xlabel('Infusion'); ylabel('Rating');

% 
% % separate feedback magnitude and valence
% p.feedback_mag = abs(p.feedback);
% p.feedback_valence = p.feedback > 0;
% 
% 
% % calculate lagged feedback representing the reward rate
% p.feedbacklag = [NaN; p.feedback(1:end-1)];
% p.feedbacklag(p.trial==1) = NaN; % make sure there is no carryover from previous subject
% 
% p.feedback_ratinglag = [NaN; p.feedback_rating(1:end-1)];
% p.feedback_ratinglag(p.trial==1) = NaN; % make sure there is no carryover from previous subject
% 
% p.feedback_maglag = [NaN; p.feedback_mag(1:end-1)];
% p.feedback_maglag(p.trial==1) = NaN; % make sure there is no carryover from previous subject
% 
% p.feedback_valencelag = [NaN; p.feedback_valence(1:end-1)];
% p.feedback_valencelag(p.trial==1) = NaN; % make sure there is no carryover from previous subject
% 
% % sigmoid transform
% p.feedback_rating_sigm = sigm(p.feedback_rating)';
% p.stim_rating_sigm = sigm(p.stim_rating)';

%% create lagged stimulus IV
% p.stimlag = [NaN; p.stim(1:end-1)];
% p.stim_diff = p.stimlag==p.stim;

%% make sure stim is nominal

p.stim = nominal(p.stim);
p.stimlag = nominal(p.stimlag);

% model1 = fitlme(p,'feedback_rating ~ 1 + feedback + trial + stim*trial + (feedback*trial + stim|subject)')
% anova(model1)
% model2 = fitlme(p,'feedback_rating ~ 1 + feedback + trial + stim*trial + (feedback + stim|subject)')
% anova(model2)

%% Feedback ratings 1: the best model to date shows that their feedback ratings are influenced by both neurofeedback and infusion

% stim_1 - infusion

feed_model = fitlme(p,'feedback_rating ~ 1 + feedback + trial + stim + (feedback + stim|subject)')
% anova(feed_model)

% % try with sigmoid transform -- same results
% feed_model_sigm = fitlme(p,'feedback_rating_sigm ~ 1 + feedback + trial + stim + (feedback + stim|subject)')
% anova(feed_model_sigm)


%% Feedback ratings 2: slight suggestion that people track their reward rate 
% indicated by a NS effect of lagged feedback rating
feed_model_back_main = fitlme(p,'feedback_rating ~ 1 + feedback + feedback_ratinglag + trial + stim + (feedback + feedback_ratinglag + stim|subject)')
% anova(feed_model_back)

% based on Marta's suggestion of adding stim_rating
feed_model_simple = fitlme(p,'feedback_rating ~ 1 + feedback*trial + stim_rating*stim + stim*trial + stim*feedback + (feedback + stim_rating + stim |subject)')

% currently best model
feed_model_lagged = fitlme(p,'feedback_rating ~ 1 + feedback*trial + stim_rating*stim + feedback_ratinglag*stim_diff*trial + stim*trial + (feedback + stim_rating + stim  |subject)')



feed_model_back = fitlme(p,'feedback_rating ~ 1 + feedback+trial + feedback_ratinglag*trial*stim_diff +  stim + (feedback + feedback_ratinglag + stim|subject)')

% feed_model_back_pure = fitlme(p,'feedback_rating ~ 1 + feedback+trial + feedbacklag*trial*stim_diff +  stim + (feedback + feedbacklag + stim|subject)')


compare(feed_model_back_pure, feed_model_back)

% 
% % try separating valence and magnitude
% feed_model_val_mag = fitlme(p,'feedback_rating ~ 1 + feedback_mag*feedback_valence + trial + stim + (feedback_mag*feedback_valence + stim|subject)')
% anova(feed_model_val_mag)

% compare(feed_model,feed_model_val_mag)


%% Expectancy ratings

% they prefer infusion, but there is no evidence of learning as indicated by NS trial*stim_1 interaction

exp_model_simple = fitlme(p,'stim_rating ~ 1 + stim*trial + feedback_ratinglag + (stim + feedback_ratinglag|subject)')


exp_model = fitlme(p,'stim_rating ~ 1 + stim*trial + stim*feedback_ratinglag + (stim + feedback_ratinglag + trial|subject)')

exp_model_congr = fitlme(p,'stim_rating ~ 1 + stim_diff*feedback_ratinglag*stim + stim*trial + stim*feedback_ratinglag + (stim + feedback_ratinglag + trial|subject)')

anova(exp_model)
compare(exp_model,exp_model_congr)

% check sigmoid transform -- same results
% exp_model = fitlme(p,'stim_rating_sigm ~ 1 + trial*stim + feedback_ratinglag + (stim + feedback_ratinglag|subject)')
% anova(exp_model)

%% get predicted responses for plotting interactions
%replicate data
t = p;
%fix predictors of no interest
t.subject = nominal(ones(size(t.subject)));

% will prediction run with undefined variables?  NO
% t.subject = nominal(NaN(size(t.subject)));


% t.trial = 30*ones(size(t.subject));
t.trial = p.trial;
t.feedback = ones(size(t.subject));



% get predicted response

[yhat, yhatCI, df] = predict(feed_model_simple,t);

% plot effects
diff = t.stim_diff==1;
same = t.stim_diff==0;
infus = eq(p.stim,'true');
noinf = eq(p.stim,'false');

% split by previous feedback
hi = t.feedback_ratinglag>0;
lo = t.feedback_ratinglag<0;

figure(99); clf;
subplot(1,2,1)
errorbar(t.trial(same & hi),yhat(same & hi),yhatCI(same & hi,1)-yhat(same & hi),yhatCI(same & hi,2)-yhat(same & hi)); hold on; 
errorbar(t.trial(diff & hi),yhat(diff & hi),yhatCI(diff & hi,1)-yhat(diff & hi),yhatCI(diff & hi,2)-yhat(diff & hi)); hold on; 

errorbar(t.trial(same & lo),yhat(same & lo),yhatCI(same & lo,1)-yhat(same & lo),yhatCI(same & lo,2)-yhat(same & lo)); hold on; 
errorbar(t.trial(diff & lo),yhat(diff & lo),yhatCI(diff & lo,1)-yhat(diff & lo),yhatCI(diff & lo,2)-yhat(diff & lo)); hold on; 


legend('fear', 'happy', 'scrambled', 'Location', 'southeast');
ylabel('predicted response time')
ax = gca;
ax.XTick = [0 1]; ax.XTickLabel = {'Omission'; 'Reward'};
title('Reward by emotion interaction on clock task response times (error bars = 95% CI)'


% [B,Bnames,stats] = randomEffects(model1);

%% Model diagnostics
figure(99); clf;
subplot(1,2,1)
F = fitted(feed_model);
R = response(feed_model);
plot(R,F,'rx')
xlabel('Response')
ylabel('Fitted feedback ratings')
subplot(1,2,2)
F = fitted(exp_model);
R = response(exp_model);
plot(R,F,'bx')
xlabel('Response')
ylabel('Fitted infusion ratings')


figure(100);
plotResiduals(exp_model,'fitted')
