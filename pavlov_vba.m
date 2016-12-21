% pavlov_vba
% fits Pavlovian RW-like model to Pecina neurofeedback placebo task
% first run read_nf_data.m
close all
% clear variables




f_fname = @f_pavlov; % evolution function (Q-learning)
% g_fname = @g_Id; % observation function (softmax mapping)
% g_fname = @g_cr; %scales associative strength linearly

fit_subjects = 1;
sanity_checks = 0;
censor_zeros = 0;
stop_after_each_subject = 1;
y_feed_ratings = 0; %% predict only feedback ratings
y_feed_and_exp_ratings = 1; %% this is a mixed-response model predicting both expectancy (infusion/no infusion) and feedback ratingss

y_sigmoid = 1;
decay = 0;
biases = 0;
infusion_expectancy = 1;
learning = 1;


%% read in data
% data_dir = '/Users/localadmin/Dropbox/data_projects/placebo_marta/';
% fname = '2subjectratings.xlsx';
% s = readtable(sprintf('%s%s', data_dir, fname));
% s = table2struct(s);
%% make CS, US
cs = [];
us = [];
for i = 1:length(feedback)
    cs(i) = strcmpi(char(stim(i)),'Infusion');
    if strcmpi(char(feedback(i)),'100% Pos Feedback');
    us(i) = 2;
    elseif strcmpi(char(feedback(i)),'50% Pos Feedback');
    us(i) = 1;
    elseif strcmpi(char(feedback(i)),'50% Neg Feedback');
    us(i) = -1;
    elseif strcmpi(char(feedback(i)),'100% Neg Feedback');
    us(i) = -2;
    end
end

cs = logical(cs);
congruent = (cs == 1 & us > 0) |  (cs == 0 & us < 0);  % code congruent trials to make sure credit is assigned to the right stimulus
check1 = FeedbackRatins;
%% sanity checks
if sanity_checks
figure(99); clf;
plot(smooth(mean(~isnan(stim_ratings(cs,:)),2)));
% plot(smooth(y1(cs))); legend Subject 1;hold on; plot(smooth(y2(cs))); legend('subject 1', 'subject 2');  
clear title xlabel ylabel;
xlabel('blue: Infusion, Red: no infusion')
ylabel('Ratings, moving average')
hold on;
plot(smooth(mean(~isnan(stim_ratings(~cs,:)),2)));
% plot(smooth(y1(~cs))); legend Subject 1;hold on; plot(smooth(y2(~cs))); legend('subject 1', 'subject 2');  
hold off;
end
% see if they learn the infusion value

if censor_zeros
feedback_ratings(feedback_ratings==0) = NaN;
stim_ratings(stim_ratings==0) = NaN;
end

%% diagnose feedback ratings distribution
diagnose = 0;
if diagnose
figure(101); clf;
for sub = 12:size(feedback_ratings,2)
    clear Fit;
    subplot(5,5,sub)
    scatter(us,feedback_ratings(:,sub)+rand./5); hold on;
    Fit = polyfit(us',feedback_ratings(:,sub),1);
    plot(polyval(us, Fit)); hold off;
%     %     scatter(us(cs),feedback_ratings(cs,sub)+rand./5)
%     hold on;
%     %     scatter(us(~cs),feedback_ratings(~cs,sub)+rand./5)
%     Fit = polyfit(us(~cs)',feedback_ratings(~cs,sub)+rand./5,1);
%     plot(polyval(Fit,us(~cs)));
% hold off; 

end
end

for ct=9:size(feedback_ratings,2) % skip the 1st subject
%% define response
if y_feed_ratings
y = feedback_ratings(:,ct);
g_fname = @g_feedback_ratings; % predicts feedback ratings as a function of expectancy and actual feedback;
elseif y_feed_and_exp_ratings
y = [feedback_ratings(:,ct) stim_ratings(:,ct)]'; % make a mixed response of feedback and infusion ratings

%% with sigmoid transform
if y_sigmoid
options.sources(1).out  = 1;
options.sources(1).type = 1;
options.sources(2).out  = 2;
options.sources(2).type = 1;
else
options.sources(1).out  = 1;
options.sources(1).type = 0;
options.sources(2).out  = 2;
options.sources(2).type = 0;
end
    
g_fname = @g_mixed;
else
y = stim_ratings(:,ct);
g_fname = @g_logistic;
end
% y2 = feedback_ratings(:,3);


options.inF.noCS = 0;
options.inG.noCS = 0;

if y_sigmoid
%     y = (sigmoid(y));
    
     y(1,:) =  y(1,:)./(max(y(1,:)) - min(y(1,:))) + .5;
     y(2,:) =  y(2,:)./(max(y(2,:)) - min(y(2,:))) + .5;

% options.binomial = 0; 
%% if fitting a linear response, rescale it to be between 
else
     y(1,:) =  y(1,:)./(max(y(1,:)) - min(y(1,:))) + .5;
     y(2,:) =  y(2,:)./(max(y(2,:)) - min(y(2,:))) + .5;
end
% y = y';

%% the u is not shifted for the feedback ratings model!!!
u = [cs 0; ...  % 1 infusion on current trial
     us 0; ...  % 2 feedback on current trial
     0 cs; ...% 3 infusion on previous trial (leading to the current, pre-update value)
     0 us; ... %4 feedback on previous trial
     congruent 0; %5 congruent on current trial
     0 congruent];     %6 congruent on previous trial

options.isYout      = zeros(size(y)) ;
options.isYout=isnan(y) ;
options.isYout(1) = 1;
options.inF.binomial_override = 1;
options.inF.decay = decay;
options.inG.decay = decay;
options.inG.biases = biases;
options.inG.infusion_expectancy = infusion_expectancy;

dim = struct('n',2,'n_theta',1,'n_phi',2);

if biases
    dim.n_phi = dim.n_phi + 2;
end
if options.inF.noCS
    dim.n = dim.n-1;
end

if decay
    dim.n = dim.n+1;
end
if infusion_expectancy
    dim.n_phi = dim.n_phi +1;
end

priors.muPhi = zeros(dim.n_phi,1);
priors.muTheta = zeros(dim.n_theta,1);

% priors.muTheta = 0.1*ones(dim.n_theta,1); %% fix LR
priors.muX0 = zeros(dim.n,1);
if decay
    priors.muTheta(2) = 3;          % initialize feedback decay parameter
    priors.muX0(3,1) = 1;           % initialize feedback salience at 1
end
priors.SigmaPhi = 1e1*eye(dim.n_phi);
priors.SigmaTheta = 1e1*eye(dim.n_theta); 
% priors.SigmaTheta = 0*eye(dim.n_theta);  %% fix LR
priors.SigmaX0 = zeros*eye(dim.n);
% priors.SigmaX0 = 1e1*eye(dim.n); %% try fitting prior expectancy
priors.a_alpha = Inf;
priors.b_alpha = 0;
options.priors = priors;
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
% displayResults(posterior,out,y,x,x0,theta,phi,Inf,Inf);
L(ct) = out.F;
if stop_after_each_subject
keyboard
end
end

% VBA_ReDisplay(posterior,out,1);
