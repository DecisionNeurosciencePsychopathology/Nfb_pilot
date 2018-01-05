function s=generate_regressor_data
%Grab the pes from subjects' vba output, currently we skipped the first
%subject and the pe hidden state is muX(3), aditionally grab the feedback
%durations from the csv file provided ('MDFBaseline.csv')

%Grab files
files = sort(glob('vba_data/subj_*_vba*.mat'));

%Load look up table
load('pavlov_lookup_table.mat');

%Load rt matrices
baseline_rts = readtable('baseline_rts.csv');
baseline_rts.RT=str2double(baseline_rts.RT); %Convert to double real quick (save later)
no_baseline_rts = readtable('no_baseline_rts.csv');

%Load in the feedback times
%subj_data = readtable('MDFBaseline.csv');
subj_data = readtable('non_interpolated_files/MdfResponseAll.csv');

%total number of runs
runs=6;

%Initalize struct
s = struct;

%Generate plots or not
make_plots=0;

%If we want to censor the baseline trials or not
censor_baseline=0;

s.n_t = 72; %Num trials for vba analysis
s.n_events = 360; %Total number of events in task. Events:( Inf. | Inf. Resp | Feed | Feed. Resp. | Baseline ) x 12 (trials) x 6 (runs) = 360
s.total_blocks = 6; %Total blocks
s.trials_per_block = s.n_t/s.total_blocks; %Trials per block
s.trial_index = 1:s.trials_per_block:s.total_blocks*s.trials_per_block; %Trial indexing

for i = 1:length(files)
    load(files{i})
    tmp_id=regexp(files{i},'\d+','match');
    tmp_id = tmp_id{:};
    
    [scan_id, block_len, cond] = look_up_id(tmp_id,pavlov_lookup_table);
    
    %Let user know current subject id
    fprintf('Constructing regressor data for Subject ID %s\n', scan_id{:})
    
    %Kick out if no scan id found
    if isempty(scan_id)
        warning('Subject not found in lookup table!')
        continue
    end
    
    %Keep cond in subject struct
    s.subjects.(scan_id{:}).cond = cond{:};
    
    %Generate any VBA plots here
    if make_plots
        gen_vba_plots(posterior,out,i,scan_id{:})
    end
    
    
    %Comment out PE stuff for now
    s.subjects.(scan_id{:}).pe = posterior.muX(3,:);
    %
    %     %Create PE parameterization using only 0,1,or -1
    %     neg_index = s.subjects.(scan_id{:}).pe<0;
    %     pos_index = s.subjects.(scan_id{:}).pe>0;
    %     s.subjects.(scan_id{:}).pe_neg_pos_ones_only=zeros(length(s.subjects.(scan_id{:}).pe),1);
    %     s.subjects.(scan_id{:}).pe_neg_pos_ones_only(neg_index)=-1;
    %     s.subjects.(scan_id{:}).pe_neg_pos_ones_only(pos_index)=1;
    %
    %     %Create PE parameterization using only 0,1,or 2
    %     s.subjects.(scan_id{:}).pe_pos_only_rounded = round(abs(s.subjects.(scan_id{:}).pe)); %Round the absolute values (or floor them?)
    %     s.subjects.(scan_id{:}).pe_pos_only_rounded(s.subjects.(scan_id{:}).pe_pos_only_rounded>=3)=2; %set max at 2
    
    s.subjects.(scan_id{:}).org_id = tmp_id;
    
    %Initialize feedback struct
    s.subjects.(scan_id{:}).feedback_onset_all_runs = [];
    s.subjects.(scan_id{:}).feedback_offset_all_runs = [];
    
    %Grab the block length as some do not have baseline condition
    s.subjects.(scan_id{:}).block_length = block_len;
    
    %If processing baseline subjects the index value is 5 otherwise it's 3
    if block_len>200
        t_idx = 5;
    else
        t_idx = 4;
    end
    
    %Grab the data needed per subject per run
    s=pull_event_times(s,subj_data,runs,scan_id{:},'Feedback');
    s=pull_event_times(s,subj_data,runs,scan_id{:},'Infusion');
    s=pull_event_times(s,subj_data,runs,scan_id{:},'Baseline'); %needed for censoring
    
    %Compile all onsets just to have data on hand
    s.subjects.(scan_id{:}).all_resp_event_onsets = zeros(s.n_t*2,1);
    s.subjects.(scan_id{:}).all_resp_event_onsets(1:2:end) = s.subjects.(scan_id{:}).infusionresponse_onset_all_runs;
    s.subjects.(scan_id{:}).all_resp_event_onsets(2:2:end) = s.subjects.(scan_id{:}).feedbackresponse_onset_all_runs;
    
        
    %Get the contrasts values (Salience,Valence, ect)
    if strcmp(scan_id{:},'NF122')
        s=pull_contrasts(s,scan_id{:},subj_data);
    end
    
    %Initialize censor struct
    s.subjects.(scan_id{:}).run_censor_full = [];
    
    %Initialize rating array
    s.subjects.(scan_id{:}).all_actual_ratings = zeros((length(s.subjects.(scan_id{:}).infusion_ratings_all_runs)...
        + length(s.subjects.(scan_id{:}).feedback_ratings_all_runs)),1);
    
    s=pull_event_times(s,subj_data,runs,scan_id{:},'CENSOR',censor_baseline);
    s.subjects.(scan_id{:}).run_censor_full = ones(s.n_events,1); %Use this line if you need to do NOT want to censor the entire runs, i.e. keep the volumes with baseline trials
    %s.subjects.(scan_id{:}).run_censor_vba = [];
    
    %Initialize the censor array that only takes NAN into account
    s.subjects.(scan_id{:}).response_nan_censor_full = ones(length(s.subjects.(scan_id{:}).run_censor_full),1);
    
    
    %Grab trial-wise censor -- censor the trials we dont want! i.e. make them 0
%     s.subjects.(scan_id{:}).feedback_censor_array =  s.subjects.(scan_id{:}).feedback_ratings_all_runs~=0; % ~isnan(s.subjects.(scan_id{:}).feedback_ratings_all_runs)
%     s.subjects.(scan_id{:}).infusion_censor_array = s.subjects.(scan_id{:}).infusion_ratings_all_runs~=0;  %~isnan(s.subjects.(scan_id{:}).infusion_ratings_all_runs)
    
%     s.subjects.(scan_id{:}).feedback_nan_censor_array = ~isnan(s.subjects.(scan_id{:}).feedback_ratings_all_runs);
%     s.subjects.(scan_id{:}).infusion_nan_recensor_array = ~isnan(s.subjects.(scan_id{:}).infusion_ratings_all_runs);
    
    %Updated
%     s.subjects.(scan_id{:}).feedback_nan_censor_array = ~isnan(s.subjects.(scan_id{:}).feedback_ratings_raw(:,1));
%     s.subjects.(scan_id{:}).infusion_nan_recensor_array = ~isnan(s.subjects.(scan_id{:}).infusion_ratings_raw(:,1));
    
    %To include 0's as well
    s.subjects.(scan_id{:}).feedback_nan_censor_array = ~isnan(s.subjects.(scan_id{:}).feedback_ratings_raw(:,1)) .* s.subjects.(scan_id{:}).feedback_ratings_all_runs(:,1)~=0;
    s.subjects.(scan_id{:}).infusion_nan_recensor_array = ~isnan(s.subjects.(scan_id{:}).infusion_ratings_raw(:,1)) .* s.subjects.(scan_id{:}).infusion_ratings_all_runs(:,1)~=0;
    
    try
        %Update full run for censor vector
%         s.subjects.(scan_id{:}).run_censor_full(3:t_idx:end) = s.subjects.(scan_id{:}).feedback_censor_array;
%         s.subjects.(scan_id{:}).run_censor_full(4:t_idx:end) = s.subjects.(scan_id{:}).feedback_censor_array;
%         s.subjects.(scan_id{:}).run_censor_full(1:t_idx:end) = s.subjects.(scan_id{:}).infusion_censor_array;
%         s.subjects.(scan_id{:}).run_censor_full(2:t_idx:end) = s.subjects.(scan_id{:}).infusion_censor_array;
        
        %Combine all ratings into a regressor just in case we need it for the future
%         s.subjects.(scan_id{:}).all_actual_ratings(1:2:end) = s.subjects.(scan_id{:}).infusion_ratings_all_runs;
%         s.subjects.(scan_id{:}).all_actual_ratings(2:2:end) = s.subjects.(scan_id{:}).feedback_ratings_all_runs;
        
        %Update full run for censor vector
        s.subjects.(scan_id{:}).response_nan_censor_full(3:t_idx:end) =  s.subjects.(scan_id{:}).feedback_nan_censor_array;
        s.subjects.(scan_id{:}).response_nan_censor_full(4:t_idx:end) =  s.subjects.(scan_id{:}).feedback_nan_censor_array;
        s.subjects.(scan_id{:}).response_nan_censor_full(1:t_idx:end) = s.subjects.(scan_id{:}).infusion_nan_recensor_array;
        s.subjects.(scan_id{:}).response_nan_censor_full(2:t_idx:end) = s.subjects.(scan_id{:}).infusion_nan_recensor_array;
        
    catch
        continue
    end
    
    
end

save pavlovian_proc_data s

%Remove the extra subjects! ie anything under <121



function [id_out, block_length, cond]=look_up_id(tmp_id,lookup_table)
%Grab the scan id
id_out=lookup_table.scan_id(lookup_table.tmp_id==str2double(tmp_id));
block_length = lookup_table.block_length(lookup_table.tmp_id==str2double(tmp_id));
cond = lookup_table.condition(lookup_table.tmp_id==str2double(tmp_id));

function s=pull_event_times(s,subj_data,runs,id,event,censor_baseline)
%This function will pull the times from MDFBaseline.csv which was
%previosuly converted to a table (subj_data)

try censor_baseline; catch,censor_baseline=0; end

%Initialize
s.subjects.(id).([lower(event) '_onset_all_runs'])=[];
s.subjects.(id).([lower(event) '_offset_all_runs'])=[];
if strcmpi(event,'INFUSION') || strcmpi(event,'FEEDBACK')
    s.subjects.(id).([lower(event) '_ratings_all_runs']) = [];
end

%Taken from non_interpolated_files/ResponseModelBasline.m, each of these
%corresponds to TrialModelCond in the main xlsx file (MdfResponseAll at the
%time of writing this). The inf/feedback and respective responses should
%always be 3 apart from one another ex: if TrialModCond is 1 for infusion
%it should be 4 for InfusionResponse and so forth. Coding for resgressor is
%as follows 1 = pos 0 = no change -1 = negative
% ConditionName =
% 1) 'infusion_positive';
% 2) 'infusion_zero';
% 3) 'infusion_negative';
% 4) 'expectation_positive';
% 5) 'expectation_zero';
% 6) 'expectation_negative';
% 7) 'feedback_positive';
% 8) 'feedback_zero';
% 9) 'feedback_negative';
% 10) 'outcome_positive';
% 11) 'outcome_zero';
% 12) 'outcome_negative';

%Pull the idx for the specific subject
subj_indices = strcmp(subj_data.x_Participant,id);

if ~strcmp(event,'CENSOR')
    %We will need to grab the infusion/feedback and infusion/feedback response pairs
    event_response = [event 'Response']; %Concat respone string
    
    %Grab the event times, in which the duration is RT by default
    s.subjects.(id).([lower(event) '_onset_all_runs']) = subj_data.Onset2(subj_indices & strcmpi(subj_data.ModelTrialPm,event));
    s.subjects.(id).([lower(event_response) '_onset_all_runs']) = subj_data.Onset2(subj_indices & strcmpi(subj_data.ModelTrialPm,event_response));
    s.subjects.(id).([lower(event) '_offset_all_runs']) = s.subjects.(id).([lower(event) '_onset_all_runs']) + subj_data.TrialModelDur(subj_indices & strcmpi(subj_data.ModelTrialPm,event));
    s.subjects.(id).([lower(event_response) '_offset_all_runs']) = s.subjects.(id).([lower(event_response) '_onset_all_runs']) + subj_data.TrialModelDur(subj_indices & strcmpi(subj_data.ModelTrialPm,event_response));
    
    %To get subject specific ratings
    
    %Pull the data for specific id
    trial_mod_cond=[subj_data.TrialModelCond((subj_indices & strcmpi(subj_data.ModelTrialPm,event)))...
        subj_data.TrialModelCond((subj_indices & strcmpi(subj_data.ModelTrialPm,event_response)))];
    
    %Store it for cesnsoring later
     s.subjects.(id).([lower(event) '_ratings_raw'])=trial_mod_cond;
    
    
    tmp = zeros(1,length(trial_mod_cond));
    for i = 1:length(trial_mod_cond)
        pairs = trial_mod_cond(i,:);
        if isnan(pairs(1))==1 || any(mod(pairs-11,3)==0)
            tmp(i)=0;
        elseif mod(pairs,3)==0
            tmp(i)=-1;
        elseif mod(pairs-10,3)==0
            tmp(i)=1;
        else
            %Just as a sanity check, throw error is the difference between the
            %event and the event resposne is something other than 3 units
            %apart.
            error('\nEvent paris for %s do not match up, further inspection required!\n',event);
        end
    end
    
    s.subjects.(id).([lower(event) '_ratings_all_runs']) = tmp';
else
    stop=0;
    %Censor only the trials in which there is no RT data for (i.e. Nan)
    %In this case 1 means we want to keep it
%     if censor_baseline && strcmp(s.subjects.(id).cond,'baseline')
%         trials_to_censor_indices = (subj_indices & isnan(subj_data.TrialModelDur));
%     else
%         trials_to_censor_indices = (subj_indices & isnan(subj_data.TrialModelDur) & ~strcmpi(subj_data.ModelTrialPm,'Baseline'));
%     end
%     (subj_data.RunNum(trials_to_censor_indices).*12) - (12-subj_data.TrialNum(trials_to_censor_indices));
    
    %s.subjects.(id).run_censor_full = subj_data
end




% %Can probably just del this for loop
% for j = 1:runs
%     if ~strcmp(event,'CENSOR')
%     else
%         %Censor out full runs if needed -- pretty sure this will break now
%         %additionally I believe it is incorrect
%         tmp_censor_data = subj_data.IncludeRun(strcmp(subj_data.x_Participant,id) &...
%             strcmp(subj_data.Run,['Run_0' num2str(j)]));
%         s.subjects.(id).run_censor_full = [s.subjects.(id).run_censor_full; strcmp(tmp_censor_data,'TRUE')];
%     end
%     
% end


function gen_vba_plots(post,out,fig_number,subj_id)
%This function will take in the outputs from the VBA toolbox and
%generate any plots needed to diagnose parametric regressor, and or
%subjects behavior.

%Figure 1 plot PE's
figure(fig_number)
clf;
plot(post.muX(3,:));
title(['PE for subject ' subj_id])

function s=pull_contrasts(s,id,subj_data)
%This function will grab the contrasts values, assuming they are all the
%same! We can just use one subject and remove the Nan's

contrasts = {'Salience' 'Valence' 'Magnitude' 'Congruent'};

for contrast=contrasts
    %Grab the data
    tmp_idx = strcmp(subj_data.x_Participant,{id});
    tmp_contrast_data = subj_data.(contrast{:});
    s.(contrast{:}) = tmp_contrast_data(tmp_idx);
    
    rm_nans = @(A) A(~isnan(A)); %Remove any Nans
    s.(contrast{:}) = rm_nans(s.(contrast{:}));
end



%
%     files = glob(['subjectRatings 2' filesep id filesep id '_run?.txt']);
%     formatSpec = '%f %d';
%     sizeData = [2 Inf];
%     event_times = [];
%     num_responses = 24;
%     for i = 1:length(files)
%         fileID = fopen(files{i},'r');
%         tmp_event_times = fscanf(fileID,formatSpec,sizeData);
%         %Some runs are skipped for various reasons
%         if isempty(tmp_event_times)
%             tmp_event_times = zeros(1,num_responses);
%         end
%         event_times = [event_times; tmp_event_times(1,:)']; %Grab only the times for now skip the ratings.
%         fclose(fileID);
%     end
%
%     if length(event_times)~=(24*6)
%         return
%         %error('Check the subjects .txt files for outliers...')
%     else
%         s.subjects.(id).any_rating = event_times;
%         s.subjects.(id).any_rating_censor = s.subjects.(id).any_rating==0;
%         s.subjects.(id).infusion_resp_rt = event_times(1:2:end);
%         s.subjects.(id).feedback_resp_rt = event_times(2:2:end);
%     end



%Old code
%Break up the censor logical array into censor events -- On second
%thought these might not be needed?
% 1 = infusion
% 2 = inf. response
% 3 = feedback
% 4 = feedback response
% 5 = baseline
%     events = {'infusion'; 'infusion_resp'; 'feedback'; 'feedback_resp'; 'baseline'};
%     n_events = length(events);
%     for j = 1:n_events
%         s.subjects.(scan_id{:}).([(events{j}) '_censor'])=s.subjects.(scan_id{:}).run_censor_full(j:n_events:end);
%     end


%         %Just grab if we include the run or not
%         tmp_censor_data = subj_data.IncludeRun(strcmp(subj_data.x_Participant,scan_id) &...
%             strcmp(subj_data.Run,['Run_0' num2str(j)]));
%         %tmp_censor_data_vba = subj_data.IncludeRun(strcmp(subj_data.x_Participant,scan_id) &...
%         %   strcmp(subj_data.Run,['Run_0' num2str(j)]) & strcmp(subj_data.ModelTrialPm,'Feedback'));
%
%         s.subjects.(scan_id{:}).feedback_onset_all_runs = [s.subjects.(scan_id{:}).feedback_onset_all_runs; s.subjects.(scan_id{:}).feedback_onset.(['run' num2str(j)])];
%         s.subjects.(scan_id{:}).feedback_offset_all_runs = [s.subjects.(scan_id{:}).feedback_offset_all_runs; s.subjects.(scan_id{:}).feedback_offset.(['run' num2str(j)])];
%         s.subjects.(scan_id{:}).run_censor_full = [s.subjects.(scan_id{:}).run_censor_full; strcmp(tmp_censor_data,'TRUE')];
%         %s.subjects.(scan_id{:}).run_censor_vba = [s.subjects.(scan_id{:}).run_censor_vba; strcmp(tmp_censor_data,'TRUE')];



%%%% -------ARCHIVE OF GETTING RT DURATION REGS------- %%%%

%%PULL IN DATA
    %Pull the detailed timings from the txt files
%     if strcmp(cond,'baseline')
%         rt_data = baseline_rts;
%     else
%         rt_data = no_baseline_rts;
%     end


%%CALL TO FUNCTION
%Not needed now as RT is the default duration
%     s=get_detailed_event_times(s,scan_id{:},rt_data,runs);

%%FUNCTION
% function s=get_detailed_event_times(s,id,rt_data,runs)
% %Initialize
% s.subjects.(id).event_rt_onset_all_runs=[];
% s.subjects.(id).event_rt_offset_all_runs=[];
%
% for j = 1:runs
%
%     %Pull event on/offset times
%     s.subjects.(id).event_rt_onset.(['run' num2str(j)])=...
%         rt_data.Onset2(strcmp(rt_data.Participant,{id}) &...
%         strcmp(rt_data.Run,['Run_0' num2str(j)])) + ...
%         rt_data.RT(strcmp(rt_data.Participant,{id}) &...
%         strcmp(rt_data.Run,['Run_0' num2str(j)])) ;
%
%     %Worry about offsets later
%     %          s.subjects.(id).event_offset.(['run' num2str(j)]) =...
%     %             s.subjects.(id).event_onset.(['run' num2str(j)])...
%     %             + rt_data.Duration(strcmp(rt_data.Participant,{id}) &...
%     %             strcmp(rt_data.Run,['Run_0' num2str(j)]));
%
%     %Compile into one vecotr
%     s.subjects.(id).event_rt_onset_all_runs = [ s.subjects.(id).event_rt_onset_all_runs;  s.subjects.(id).event_rt_onset.(['run' num2str(j)])];
%     %s.subjects.(id).event_offset_all_runs = [s.subjects.(id).event_offset_all_runs; s.subjects.(id).event_offset.(['run' num2str(j)])];
%
% end
% %Might have to make Nan's 0's
% s.subjects.(id).any_rating = s.subjects.(id).event_rt_onset_all_runs;
% s.subjects.(id).any_rating_censor = ~isnan(s.subjects.(id).any_rating);
% s.subjects.(id).infusion_resp_rt = s.subjects.(id).event_rt_onset_all_runs(1:2:end);
% s.subjects.(id).feedback_resp_rt = s.subjects.(id).event_rt_onset_all_runs(2:2:end);
