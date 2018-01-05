function s=write_placebo_regressors
%This function will take the data generated from the
%generate_regressor_data.m script and put them in FSL styled .dat files to
%be used in fMRI regressor model analysis via AFNI


%This variable name is currently s -- only load is arg is not passed
%load('pavlovian_data_pe_and_feedback_times.mat');
load('pavlovian_proc_data.mat');

fprintf('\nCreating subject specific regressor files\n\n');

%Make regs folder if it doesn't exist
data_dump_str = 'regs/';
if ~exist(data_dump_str,'file')
    mkdir(data_dump_str)
end


%Grab the subjects in the data set
subjs = fieldnames(s.subjects);

%Make salience binary, and mean corrected
s.Salience_binary = s.Salience == 1;
s.Salience_binary_MC = s.Salience_binary-mean(s.Salience_binary);
s.Mag12(s.Magnitude==1) = 2;
s.Mag12(s.Magnitude==-1) = 1;




for i = 1:length(subjs)
    s.id = subjs{i};
    try
        
        %This is when the bar appears to increase
        %Infusion_cs2_onset = s.subjects.(s.id).feedback_onset_all_runs-4;
        %For now we don't even have to worry about it since we are using a stick
        %Infusion_cs2_offset = s.subjects.(s.id).feedback_onset-4;
        
        %Create the censor regressor TR = 2
        s=createPlaceboCensorRegressor(s,s.total_blocks);
        s.feeback_censor_value='pos';
        s_pos=createPlaceboCensorRegressor(s,s.total_blocks);
        s.feeback_censor_value='neg';
        s_neg=createPlaceboCensorRegressor(s,s.total_blocks);
        
        %% Trial
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).infusion_onset_all_runs,s.subjects.(s.id).feedbackresponse_offset_all_runs,sprintf('%s_trial',s.id),true(s.n_t,1),0,s);
        
        %% Infusion + feedback ratings
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).any_rating,s.subjects.(s.id).any_rating+1,sprintf('%s_all_ratings',s.id),s.subjects.(s.id).any_rating_censor,0,s); %Stick no duration
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).all_resp_event_onsets,s.subjects.(s.id).any_rating,sprintf('%s_all_ratings_rt_convolv',s.id),s.subjects.(s.id).any_rating_censor,0,s); %rt convolved
        
        %% Using the actual ratings taken from the response column of MDF all subjs
%         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).all_resp_event_onsets,s.subjects.(s.id).any_rating,sprintf('%s_all_responses_rt_convolv',s.id),s.subjects.(s.id).all_actual_ratings,0,s); %rt convolved
         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).infusionresponse_onset_all_runs,s.subjects.(s.id).infusionresponse_offset_all_runs,sprintf('%s_infusion_ratings_parametric',s.id),s.subjects.(s.id).infusion_ratings_all_runs,0,s); %rt convolved
         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedbackresponse_onset_all_runs,s.subjects.(s.id).feedbackresponse_offset_all_runs,sprintf('%s_feedback_ratings_parametric',s.id),s.subjects.(s.id).feedback_ratings_all_runs,0,s); %rt convolved
%         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedbackresponse_onset_all_runs,s.subjects.(s.id).feedback_resp_rt,sprintf('%s_Mag12_feedback_responses_aligned_rt_convolv',s.id),s.Mag12',0,s); %rt convolved
%         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedbackresponse_onset_all_runs,s.subjects.(s.id).feedback_resp_rt,sprintf('%s_Mag12_feedback_responses_aligned_rt_convolv_MC',s.id),s.Mag12'-mean(s.Mag12),0,s); %rt convolved
        
        %% Create a mag*valence regressor  fb rating convolved with RT
%         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedbackresponse_onset_all_runs,s.subjects.(s.id).feedback_resp_rt,sprintf('%s_MagVal_feedback_responses_aligned_rt_convolv',s.id),s.Mag12'.*s.Valence,0,s); %rt convolved
        
        
%         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).infusion_onset_all_runs,s.subjects.(s.id).infusionresponse_onset_all_runs,sprintf('%s_infusion_responses_inf_aligned_RTboxcar',s.id),s.subjects.(s.id).infusion_ratings_all_runs,0,s); %rt convolved
%         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedback_onset_all_runs,s.subjects.(s.id).feedbackresponse_onset_all_runs,sprintf('%s_feedback_responses_feed_aligned_RTboxcar',s.id),s.subjects.(s.id).feedback_ratings_all_runs,0,s); %rt convolved
        
        
        %% Feedback_flux- i.e. onset2 = feedback fluctuation onset, not the actual onset of when feedback is shown
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedback_onset_all_runs,s.subjects.(s.id).feedback_offset_all_runs,sprintf('%s_feedback',s.id),s.subjects.(s.id).feedback_censor_array',0,s);
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedback_onset_all_runs,s.subjects.(s.id).feedback_offset_all_runs,sprintf('%s_feedback',s.id),true(s.n_t,1),0,s);
        write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedback_onset_all_runs,s.subjects.(s.id).feedback_offset_all_runs,sprintf('%s_feedback_flux',s.id),true(s.n_t,1),0,s); %Stick no duration
%         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedback_onset_all_runs,s.subjects.(s.id).feedback_onset_all_runs+.5,sprintf('%s_feedback_flux',s.id),true(s.n_t,1),0,s); %Stick no duration
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedback_onset_all_runs,s.subjects.(s.id).feedback_offset_all_runs,sprintf('%sPE_feedback_flux_aligned',s.id),s.subjects.(s.id).pe',0,s); %Stick no duration
        %feedback_inflection (flux) * msg and valence
%         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedback_onset_all_runs,s.subjects.(s.id).feedback_onset_all_runs+.5,sprintf('%s_mag12_convolv_feedback_flux',s.id),s.Mag12',0,s); %Stick no duration
%         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedback_onset_all_runs,s.subjects.(s.id).feedback_onset_all_runs+.5,sprintf('%s_valence_convolv_feedback_flux',s.id),s.Valence,0,s); %Stick no duration
        
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedback_onset_all_runs,s.subjects.(s.id).feedback_offset_all_runs,sprintf('%s_Salience_feedback_flux',s.id),s.Salience,0,s); %Infusion=1 no Infusion=-1
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedback_onset_all_runs,s.subjects.(s.id).feedback_offset_all_runs,sprintf('%s_Valence_feedback_flux',s.id),s.Valence,0,s); %Positive=1 Negative=-1
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedback_onset_all_runs,s.subjects.(s.id).feedback_offset_all_runs,sprintf('%s_Magnitude_feedback_flux',s.id),s.Magnitude,0,s); %High Feedback=1 Low Feedback=-1
%         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedback_onset_all_runs,s.subjects.(s.id).feedback_offset_all_runs,sprintf('%s_Congruent_feedback_flux',s.id),s.Congruent,0,s); %Congruent=1 Not Congruent=-1
        
        %% Feedback parametric regressors
        %RT convolv aligned with feedback rating
        write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedbackresponse_onset_all_runs,s.subjects.(s.id).feedbackresponse_offset_all_runs,sprintf('%s_Valence_feedback_rating',s.id),s.Valence,0,s); %Positive=1 Negative=-1
        write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedbackresponse_onset_all_runs,s.subjects.(s.id).feedbackresponse_offset_all_runs,sprintf('%s_Magnitude_feedback_rating',s.id),s.Magnitude,0,s); %High Feedback=1 Low Feedback=-1
        
        %Aligned with feedback rating stick
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedbackresponse_onset_all_runs,s.subjects.(s.id).feedbackresponse_onset_all_runs+.5,sprintf('%s_Valence_feedback_rating_stick',s.id),s.Valence,0,s); %Positive=1 Negative=-1
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedbackresponse_onset_all_runs,s.subjects.(s.id).feedbackresponse_onset_all_runs+.5,sprintf('%s_Magnitude_feedback_rating_stick',s.id),s.Magnitude,0,s); %High Feedback=1 Low Feedback=-1
        
        %Mag and valence by feedback ratings
        write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedbackresponse_onset_all_runs,s.subjects.(s.id).feedbackresponse_offset_all_runs,sprintf('%s_Valence_X_ratings',s.id),s.Valence.*s.subjects.(s.id).feedback_ratings_all_runs,0,s); %Positive=1 Negative=-1
        write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedbackresponse_onset_all_runs,s.subjects.(s.id).feedbackresponse_offset_all_runs,sprintf('%s_Magnitude_X_ratings',s.id),s.Magnitude.*s.subjects.(s.id).feedback_ratings_all_runs,0,s); %High Feedback=1 Low Feedback=-1

        %% Feedback ratings
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedback_resp_rt,s.subjects.(s.id).feedback_offset_all_runs,sprintf('%s_feedback_ratings',s.id),s.subjects.(s.id).feedback_censor_array',0,s);
%         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedback_resp_rt,s.subjects.(s.id).feedback_resp_rt+1,sprintf('%s_feedback_rating_times',s.id),s.subjects.(s.id).feedback_censor_array,0,s); %Stick (Previously %s_feedback_ratings)
         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedbackresponse_onset_all_runs,s.subjects.(s.id).feedbackresponse_offset_all_runs,sprintf('%s_feedback_rating_times',s.id),true(s.n_t,1),0,s); %rt convolved
         %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedbackresponse_onset_all_runs,s.subjects.(s.id).feedbackresponse_offset_all_runs,sprintf('%s_feedback_ratings',s.id),s.subjects.(s.id).feedback_ratings_all_runs,0,s); %rt convolved
         %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedbackresponse_onset_all_runs,s.subjects.(s.id).feedbackresponse_onset_all_runs+.5,sprintf('%s_feedback_ratings_times_stick',s.id),true(s.n_t,1),0,s); %rt convolved
%         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedbackresponse_onset_all_runs,s.subjects.(s.id).feedback_resp_rt,sprintf('%s_Valence_feedback_rating_times_rt_convolv',s.id),s.Valence,0,s); %Positive=1 Negative=-1
%         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedbackresponse_onset_all_runs,s.subjects.(s.id).feedback_resp_rt,sprintf('%s_Magnitude_feedback_rating_times_rt_convolv',s.id),s.Magnitude,0,s); %High Feedback=1 Low Feedback=-1
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedbackresponse_onset_all_runs,s.subjects.(s.id).feedback_resp_rt,sprintf('%s_infusion_ratingsXinf_no_inf_feedback_ratings_aligned_rt_convolv',s.id),s.subjects.(s.id).feedback_censor_array.*s.subjects.(s.id).infusion_ratings_all_runs.*s.Salience_binary_MC,0,s); %Expectancy ratings x inf/noInf aligned with feedback ratings, RT convolved
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedbackresponse_onset_all_runs,s.subjects.(s.id).feedback_resp_rt,sprintf('%s_infusion_NoInfusion_feedback_ratings_aligned_rt_convolved_MC',s.id),s.Salience_binary_MC,0,s); %mean corrected inf no infusion, rt convolved
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedbackresponse_onset_all_runs,s.subjects.(s.id).feedback_resp_rt,sprintf('%s_infusion_ratings_feedback_ratings_aligned_rt_convolved',s.id),s.subjects.(s.id).feedback_censor_array.*s.subjects.(s.id).infusion_ratings_all_runs,0,s); %infusion ratings aligned with feedback ratings rt convolved, multiply by censor array since we censor 0's
        
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedbackresponse_onset_all_runs,s.subjects.(s.id).feedbackresponse_offset_all_runs,sprintf('%sPE_feedback_ratings_aligned',s.id),s.subjects.(s.id).pe',0,s); %Stick no duration
        %s.subjects.(s.id).pe_mc = s.subjects.(s.id).pe - mean(s.subjects.(s.id).pe);
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedbackresponse_onset_all_runs,s.subjects.(s.id).feedbackresponse_offset_all_runs,sprintf('%sPE_feedback_ratings_aligned_mc',s.id),s.subjects.(s.id).pe_mc',0,s); %Stick no duration
        
        %% Infusion
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedback_onset_all_runs,s.subjects.(s.id).feedback_offset_all_runs,sprintf('%s_infusion',s.id),s.subjects.(s.id).infusion_censor_array',0,s);
        write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).infusion_onset_all_runs,s.subjects.(s.id).infusion_offset_all_runs,sprintf('%s_infusion',s.id),true(s.n_t,1),0,s);
%         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).infusion_onset_all_runs,s.subjects.(s.id).infusion_onset_all_runs+.5,sprintf('%s_infusion_times',s.id),s.subjects.(s.id).infusion_censor_array,0,s);
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).infusion_onset_all_runs,s.subjects.(s.id).infusion_onset_all_runs+.5,sprintf('%s_infusion_NoInfusion',s.id),s.Salience_binary,0,s); %-1 is no infusion
        write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).infusionresponse_onset_all_runs,s.subjects.(s.id).infusionresponse_offset_all_runs,sprintf('%s_infusion_NoInfusion_inf_rating',s.id),s.Salience_binary,0,s); %-1 is no infusion
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).infusionresponse_onset_all_runs,s.subjects.(s.id).infusionresponse_onset_all_runs+.5,sprintf('%s_infusion_NoInfusion_inf_rating_stick',s.id),s.Salience_binary,0,s); %-1 is no infusion
        
        %% Infusion cs2 -- just use FS1 in the model
%         write3Ddeconv_startTimes(data_dump_str,Infusion_cs2_onset,Infusion_cs2_onset+.5,sprintf('%s_infusion_cs2_times',s.id),s.subjects.(s.id).infusion_censor_array,0,s);
        
        %% Infusion rating times
        write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).infusionresponse_onset_all_runs,s.subjects.(s.id).infusionresponse_offset_all_runs,sprintf('%s_infusion_rating_times',s.id),true(s.n_t,1),0,s); %RT convolved
        write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).infusionresponse_onset_all_runs,s.subjects.(s.id).infusionresponse_offset_all_runs,sprintf('%s_infusion_NoInfusion_X_ratings',s.id),s.subjects.(s.id).infusion_ratings_all_runs.*s.Salience_binary,0,s); %infNoinf x inf ratings
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).infusionresponse_onset_all_runs,s.subjects.(s.id).infusionresponse_onset_all_runs+.5,sprintf('%s_infusion_rating_times_stick',s.id),true(s.n_t,1),0,s); %Stick
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).infusionresponse_onset_all_runs,s.subjects.(s.id).infusionresponse_offset_all_runs,sprintf('%s_expectancy_ratings',s.id),s.subjects.(s.id).infusion_ratings_all_runs,0,s); %RT convolved with actual responses
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).infusion_resp_rt,s.subjects.(s.id).infusion_offset_all_runs,sprintf('%s_infusion_ratings',s.id),s.subjects.(s.id).infusion_censor_array',0,s);
%         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).infusion_resp_rt,s.subjects.(s.id).infusion_resp_rt+1,sprintf('%s_infusion_rating_times',s.id),s.subjects.(s.id).infusion_censor_array,0,s); %Stick (Previously %s_infusion_ratings) aka motor as well
%         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).infusionresponse_onset_all_runs,s.subjects.(s.id).infusion_resp_rt,sprintf('%s_infusion_rating_times_rt_convolv',s.id),s.subjects.(s.id).infusion_censor_array,0,s); %Stick
%         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).infusionresponse_onset_all_runs,s.subjects.(s.id).infusion_resp_rt,sprintf('%s_Salience_infusion_rating_times_rt_convolv',s.id),s.Salience,0,s); %Infusion=1 no Infusion=-1
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).infusion_resp_rt,s.subjects.(s.id).infusion_resp_rt+1,sprintf('%s_infusion_rating_times_Inf_NoInf',s.id),s.subjects.(s.id).infusion_censor_array.*s.Salience_binary,0,s); %-1 is no infusion
%         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).infusionresponse_onset_all_runs,s.subjects.(s.id).infusion_resp_rt,sprintf('%s_infusion_rating_times_Inf_NoInf_rt_convolv',s.id),s.subjects.(s.id).infusion_censor_array.*s.Salience_binary,0,s); %rt convolved
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).infusionresponse_onset_all_runs,s.subjects.(s.id).infusion_resp_rt,sprintf('%s_infusion_NoInfusion_infusion_ratings_aligned_rt_convolved_MC',s.id),s.Salience_binary_MC,0,s); %mean corrected inf no infusion, rt convolved
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).infusionresponse_onset_all_runs,s.subjects.(s.id).infusion_resp_rt,sprintf('%s_infusion_ratings_infusion_ratings_aligned_rt_convolved',s.id),s.subjects.(s.id).infusion_censor_array.*s.subjects.(s.id).infusion_ratings_all_runs,0,s); %infusion ratings aligned with inf ratings rt convolved, multiply by censor array since we censor 0's
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).infusionresponse_onset_all_runs,s.subjects.(s.id).infusion_resp_rt,sprintf('%s_infusion_ratingsXinf_no_inf_infusion_ratings_aligned_rt_convolved',s.id),s.subjects.(s.id).infusion_censor_array.*s.subjects.(s.id).infusion_ratings_all_runs.*s.Salience_binary_MC,0,s); %Expectancy ratings x inf/noInf aligned with infusion ratings, RT convolved
        
        %% PE
        %         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedback_onset_all_runs,s.subjects.(s.id).feedback_onset_all_runs+1,sprintf('%s_PE',s.id),s.subjects.(s.id).pe',0,s); %stick
        %         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedback_resp_rt,s.subjects.(s.id).feedback_resp_rt+1,sprintf('%s_PE_feedbackRating_aligned',s.id),s.subjects.(s.id).pe',0,s); %Aligned with fb ratings, stick
        %         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedbackresponse_onset_all_runs,s.subjects.(s.id).feedback_resp_rt,sprintf('%s_PE_feedbackRating_aligned_rt_convolv',s.id),s.subjects.(s.id).pe',0,s); %Aligned with fb ratings, stick
        
        %         %only have parametic values of 0,1 or 2
        %         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedback_onset_all_runs,s.subjects.(s.id).feedback_onset_all_runs+1,sprintf('%s_PE_plus12',s.id),s.subjects.(s.id).pe_pos_only_rounded',0,s); %stick
        %         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedback_resp_rt,s.subjects.(s.id).feedback_resp_rt+1,sprintf('%s_PE_plus12_feedbackRating_aligned',s.id),s.subjects.(s.id).pe_pos_only_rounded',0,s); %stick
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedbackresponse_onset_all_runs,s.subjects.(s.id).feedback_resp_rt,sprintf('%s_PE_plus12_feedbackRating_rt_convolv',s.id),s.subjects.(s.id).pe_pos_only_rounded',0,s); %rt convolved
        %         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedbackresponse_onset_all_runs,s.subjects.(s.id).feedback_resp_rt,sprintf('%s_PE_plus12_feedbackRating_rt_convolv_InfnoInf_contrast',s.id),s.subjects.(s.id).pe_pos_only_rounded'.*s.Salience_binary,0,s); %rt convolved
        
        %         %only have parametic values of 0,1 or -1.
        %         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedback_onset_all_runs,s.subjects.(s.id).feedback_onset_all_runs+1,sprintf('%s_PE_plus_minus1',s.id),s.subjects.(s.id).pe_neg_pos_ones_only,0,s); %stick
        %write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedback_resp_rt,s.subjects.(s.id).feedback_resp_rt+1,sprintf('%s_PE_plus_minus1_feedbackRating_aligned',s.id),s.subjects.(s.id).pe_neg_pos_ones_only,0,s); %stick
        %         write3Ddeconv_startTimes(data_dump_str,s.subjects.(s.id).feedbackresponse_onset_all_runs,s.subjects.(s.id).feedback_resp_rt,sprintf('%s_PE_plus_minus1_feedbackRating_aligned_rt_convolv',s.id),s.subjects.(s.id).pe_neg_pos_ones_only,0,s); %rt convolved
        
        
        %% Censor
        gdlmwrite([data_dump_str sprintf('%sCensorOnly.regs',s.id)],s.hrf_regs.to_censor');
%         gdlmwrite([data_dump_str sprintf('%sPosFeedbackCensor.regs',s.id)],s_pos.hrf_regs.to_censor');
%         gdlmwrite([data_dump_str sprintf('%sNegFeedbackCensor.regs',s.id)],s_neg.hrf_regs.to_censor');
        gdlmwrite([data_dump_str sprintf('%sPosFeedbackCensor_inclusive.regs',s.id)],s_pos.hrf_regs.to_censor'); %Inclusive meaing that we are including the missed response trials
        gdlmwrite([data_dump_str sprintf('%sNegFeedbackCensor_inclusive.regs',s.id)],s_neg.hrf_regs.to_censor');
    catch
        s.ids_skipped{i} =  s.id;
        continue
    end
end


function [x,y]=write3Ddeconv_startTimes(file_loc,event_beg,event_end,fname,modulator,noFSL,b)
% Function will write FSL styled regressors in dat files for fMRI analysis
% Inputs:
% file_loc: file location (str)
% event_beg: the time in miliseconds of the event beginning
% event_end: the time in milliseconds of the event ending
% fname: the file names
% censor: the censor vector or parametric value vector depending on the regressor
% noFSL either to write a FSL file or a different single line version (see 3dDeconvolve help for more info)
% trial_index: the position of when a new block starts (trialwise)

if nargin <6
    %censor = 1;
    noFSL=0;
end
format long

x(:,1) = event_beg';
x(:,2) = event_end'-event_beg';
%x=x./1000; %Convert to seconds (not for clock already in seconds)
x(:,3) = ones(length(x),1).*modulator; %originally was modulator'
%write the -stim_times_FSL

if ~noFSL
    %Save to regs folder
    %dlmwrite([file_loc fname '.dat'],x,'delimiter','\t','precision','%.6f')
    c = asterisk(x,b); %Add in asterisks and clean up data
    dlmcell([file_loc fname '.dat'],c,'delimiter','\t')
    %dlmcell([data_dump_str filename],c,'\t');
    y=0;
else
    %write the -stim_times file
    fname = [fname '_noFSL'];
    y = x(logical(x(:,3)),1)';
    %Quick fix hack for just first ten trials troubleshoot SPMG2
    %y = y(1:10);
    dlmwrite([file_loc fname '.dat'],y,'delimiter','\t','precision','%.6f')
end
return

function c = asterisk(x,b)
%adding asterisk to existing .dat files also removes any nans present

c=[];
ast = {'*', '*', '*'};
for i = 1:length(b.trial_index)
    %Set up trial ranges
    trial_index_1 = b.trial_index(i);
    trial_index_2 = trial_index_1 + b.trials_per_block-1;
    block_data = num2cell(x(trial_index_1: trial_index_2,:));
    if i<length(b.trial_index)
        c = [c; block_data; ast];
    else
        c = [c; block_data;];
    end
end

%clean up any nans
%fh = @(y) all(isnan(y(:)));
c = c(~any(cellfun(@isnan,c),2),:);
%c(cellfun(fh, c)) = [];
%Check on c!

return


function foo = createSimpleRegressor(event_begin,event_end,epoch_window,conditional_trials)
% this was not a problem earlier, but for some reason it is now: find indices that would
% result in a negative value and set them both to 0
qbz = ( event_begin == 0 ); qez = ( event_end == 0 );
event_begin( qbz | qez ) = 0; event_end( qbz | qez ) = 0;

% check if optional censoring variable was used
if(~exist('conditional_trials','var') || isempty(conditional_trials))
    conditional_trials = true(length(event_begin),1);
elseif(~islogical(conditional_trials))
    % needs to be logical format to index cells
    conditional_trials = logical(conditional_trials);
end

% this only happened recently, but it's weird
if(any((event_end(conditional_trials)-event_begin(conditional_trials)) < 0))
    error('MATLAB:bandit_fmri:time_travel','feedback is apparently received before RT');
end

% create epoch windows for each trial
epoch = arrayfun(@(a,b) a:b,event_begin,event_end,'UniformOutput',false);

% for each "epoch" (array of event_begin -> event_end), count events
% per_event_histcs = cellfun(@(h) histc(h,epoch_window),epoch(conditional_trials),'UniformOutput',false);
% foo = logical(sum(cell2mat(per_event_histcs),1));

foo = zeros(size(epoch_window));

for n = 1:numel(epoch)
    if(conditional_trials(n))
        foo = logical(foo + histc(epoch{n},epoch_window));
    end
end


% createAndCatRegs(event_begin,event_end,epoch_window);

return