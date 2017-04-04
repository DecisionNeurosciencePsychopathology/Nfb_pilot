% fits Pavlovian RW-like model or a no-learning model to Pecina neurofeedback placebo task
% first run read_nf_data.m
close all
% clear variables
% read_nf_data




fit_subjects = 1;
sanity_checks = 0;
diagnose = 0;
censor_zeros = 0;
stop_after_each_subject = 0;
y_feed_ratings = 0; %% predict only feedback ratings -- currently broken 
y_feed_and_exp_ratings = 1; %% this is a mixed-response model predicting both expectancy (infusion/no infusion) and feedback ratingss

%If we want to fix the learning rate and by what value
options.inF.fixed_learning_rate = 1;
options.inF.fixed_alpha = -1.3801; %This will force the exp to be .2009

%To save or not to save results
save_results=1;

y_sigmoid = 1;
decay = 0; 
biases = 0;
infusion_expectancy = 1;
learning = 1; %Execute the learning rule
graphics = 0; %Tuen graphics on and off
track_pe = 1; %If you want to track pe as a hidden state

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
    hold off;
end

if censor_zeros
    feedback_ratings(feedback_ratings==0) = NaN;
    stim_ratings(stim_ratings==0) = NaN;
end

%% diagnose feedback ratings distribution
if diagnose
    figure(101); clf;
    for sub = 1:size(feedback_ratings,2)
        clear Fit;
        subplot(5,5,sub)
        scatter(us,feedback_ratings(:,sub)+rand./5); hold on;
        Fit = polyfit(us',feedback_ratings(:,sub),1);
        plot(polyval(us, Fit)); hold off;
        
    end
end

for ct=1:size(feedback_ratings,2) % skip the 1st subject
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
    
    if y_sigmoid && y_feed_ratings
        %     y = (sigmoid(y));
        y =  y./(max(y) - min(y)) + .5;
    elseif y_sigmoid && y_feed_and_exp_ratings
        
        y(1,:) =  y(1,:)./(max(y(1,:)) - min(y(1,:))) + .5;
        y(2,:) =  y(2,:)./(max(y(2,:)) - min(y(2,:))) + .5;
        
        % options.binomial = 0;
    end
    % y = y';
    
    %% the u is not shifted for the feedback ratings model!!!
    u = [cs 0; ...  % 1 infusion on current trial
        us 0; ...  % 2 feedback on current trial
        0 cs; ...% 3 infusion on previous trial (leading to the current, pre-update value)
        0 us; ... %4 feedback on previous trial
        congruent 0; %5 congruent on current trial
        0 congruent];     %6 congruent on previous trial
    
    %options.isYout      = zeros(size(y)) ;
    options.isYout=isnan(y) ;
    options.isYout(1) = 1;
    options.inF.binomial_override = 1;
    options.inF.decay = decay;
    options.inG.decay = decay;
    options.inG.biases = biases;
    options.inG.infusion_expectancy = infusion_expectancy;
    options.inF.track_pe = track_pe;
    options.inG.track_pe = track_pe;
    
    
    if learning
        f_fname = @f_pavlov; % evolution function (Q-learning)
        dim = struct('n',2,'n_theta',1,'n_phi',2);
        priors.muTheta = zeros(dim.n_theta,1);
        priors.SigmaTheta = 1e1*eye(dim.n_theta);
        % priors.muTheta = 0.1*ones(dim.n_theta,1); %% fix LR
        priors.muTheta(1) = -1.3801; %% fix LR at .2
        % priors.SigmaTheta = 0*eye(dim.n_theta);  %% fix LR
        
        if decay
            dim.n = dim.n+1;
            priors.muTheta(2) = 3;          % initialize feedback decay parameter
            priors.muX0(3,1) = 1;           % initialize feedback salience at 1
        end
        
    else
        dim = struct('n',0,'n_theta',0,'n_phi',2);
        g_fname = @g_mixed_no_learn;
        
        
    end
    
    if biases
        dim.n_phi = dim.n_phi + 2;
    end
    if options.inF.noCS
        dim.n = dim.n-1;
    end
    
    if infusion_expectancy
        dim.n_phi = dim.n_phi +1;
    end
    
    if track_pe
        dim.n = dim.n+1;
    end
    
    priors.muPhi = zeros(dim.n_phi,1);
    priors.SigmaPhi = 1e1*eye(dim.n_phi);
    
    priors.a_alpha = Inf;
    priors.b_alpha = 0;
    options.priors = priors;
    
    %Turn graphics on or off
    
    options.DisplayWin = graphics;
    options.GnFigs = graphics;
    
    %Set prior covariance of x_0
    priors.SigmaX0 = zeros*eye(dim.n);
    % priors.SigmaX0 = 1e1*eye(dim.n); %% try fitting prior expectancy
    
    %Set prior mean for x_0
    priors.muX0 = zeros(dim.n,1);
    
    
    if learning
        [posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
    else
        [posterior,out] = VBA_NLStateSpaceModel(y,u,[],g_fname,dim,options);
    end
    % displayResults(posterior,out,y,x,x0,theta,phi,Inf,Inf);
    L(ct) = out.F;
    
    if save_results
        if ~exist('pavlov_lookup_table','var')
            load('pavlov_lookup_table')
        end
        %Just to make sure get the tmp id of the subject
        subj = subjs{ct};
        subj = subj(end-2:end);
        idx = find(~cellfun('isempty',strfind(pavlov_lookup_table.scan_id,subj)));
        tmp_id = pavlov_lookup_table.tmp_id(idx);
        
        %Finally save the data
        if ~isemtpy(idx)
            save(sprintf('vba_data/subj_%d_vba_data.mat',tmp_id),'out','posterior')
        end
        
    end
    
    if stop_after_each_subject
        keyboard
    end
end

% VBA_ReDisplay(posterior,out,1);
