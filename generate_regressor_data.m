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
subj_data = readtable('MDFAllSubjects.csv');

%total number of runs
runs=6;

%Initalize struct
s = struct;

%Generate plots or not
make_plots=0;

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
    
    %Generate any VBA plots here
    if make_plots
        gen_vba_plots(posterior,out,i,scan_id{:})
    end
    
    s.subjects.(scan_id{:}).pe = posterior.muX(3,:);
    
    %Create PE parameterization using only 0,1,or -1
    neg_index = s.subjects.(scan_id{:}).pe<0;
    pos_index = s.subjects.(scan_id{:}).pe>0;
    s.subjects.(scan_id{:}).pe_neg_pos_ones_only=zeros(length(s.subjects.(scan_id{:}).pe),1);
    s.subjects.(scan_id{:}).pe_neg_pos_ones_only(neg_index)=-1;
    s.subjects.(scan_id{:}).pe_neg_pos_ones_only(pos_index)=1;
    
    %Create PE parameterization using only 0,1,or 2    
    s.subjects.(scan_id{:}).pe_pos_only_rounded = round(abs(s.subjects.(scan_id{:}).pe)); %Round the absolute values (or floor them?) 
    s.subjects.(scan_id{:}).pe_pos_only_rounded(s.subjects.(scan_id{:}).pe_pos_only_rounded>=3)=2; %set max at 2
    
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
    s=pull_event_times(s,subj_data,runs,scan_id{:},'FeedbackResponse');
    s=pull_event_times(s,subj_data,runs,scan_id{:},'Infusion');
    s=pull_event_times(s,subj_data,runs,scan_id{:},'InfusionResponse');
    s=pull_event_times(s,subj_data,runs,scan_id{:},'Baseline');
    
    s.subjects.(scan_id{:}).all_resp_event_onsets = zeros(s.n_t*2,1);
    s.subjects.(scan_id{:}).all_resp_event_onsets(1:2:end) = s.subjects.(scan_id{:}).infusionresponse_onset_all_runs;
    s.subjects.(scan_id{:}).all_resp_event_onsets(2:2:end) = s.subjects.(scan_id{:}).feedbackresponse_onset_all_runs;
    
    
    %Pull the detailed timings from the txt files
    if strcmp(cond,'baseline')
        rt_data = baseline_rts;
    else
        rt_data = no_baseline_rts;
    end
    
    s=get_detailed_event_times(s,scan_id{:},rt_data,runs);
    
    
    %Get the contrasts values (Salience,Valence, ect)
    if strcmp(scan_id{:},'NF122')
        s=pull_contrasts(s,scan_id{:},subj_data);
    end
            
    %Initialize censor struct
    s.subjects.(scan_id{:}).run_censor_full = [];
    s=pull_event_times(s,subj_data,runs,scan_id{:},'CENSOR');
    %s.subjects.(scan_id{:}).run_censor_full = ones(s.n_events,1); %Use this line if you need to do NOT want to censor the entire runs, i.e. keep the volumes with baseline trials
    %s.subjects.(scan_id{:}).run_censor_vba = [];
    
    %Grab trial-wise censor -- censor the trials we dont want! i.e. make them 0
    s.subjects.(scan_id{:}).feedback_censor_array = ~isnan(out.y(1,:));
    s.subjects.(scan_id{:}).infusion_censor_array = ~isnan(out.y(2,:));
    
    try
        %Update full run for censor vector
        s.subjects.(scan_id{:}).run_censor_full(3:t_idx:end) = s.subjects.(scan_id{:}).feedback_censor_array;
        s.subjects.(scan_id{:}).run_censor_full(4:t_idx:end) = s.subjects.(scan_id{:}).feedback_censor_array;
        s.subjects.(scan_id{:}).run_censor_full(1:t_idx:end) = s.subjects.(scan_id{:}).infusion_censor_array;
        s.subjects.(scan_id{:}).run_censor_full(2:t_idx:end) = s.subjects.(scan_id{:}).infusion_censor_array;
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

function s=pull_event_times(s,subj_data,runs,id,event)
%This function will pull the times from MDFBaseline.csv which was
%previosuly converted to a table (subj_data)

%Initialize
s.subjects.(id).([lower(event) '_onset_all_runs'])=[];
s.subjects.(id).([lower(event) '_offset_all_runs'])=[];
if strcmpi(event,'INFUSION') || strcmpi(event,'FEEDBACK')
    s.subjects.(id).([lower(event) '_ratings_all_runs']) = [];
end

for j = 1:runs
    if ~strcmp(event,'CENSOR')
        %Pull event on/offset times
        s.subjects.(id).([lower(event) '_onset']).(['run' num2str(j)])=...
            subj_data.Onset2(strcmp(subj_data.x_Participant,{id}) &...
            strcmp(subj_data.Run,['Run_0' num2str(j)]) & strcmp(subj_data.ModelTrialPm,event));
        
        s.subjects.(id).([lower(event) '_offset']).(['run' num2str(j)]) =...
            s.subjects.(id).([lower(event) '_onset']).(['run' num2str(j)])...
            + subj_data.Duration(strcmp(subj_data.x_Participant,{id}) &...
            strcmp(subj_data.Run,['Run_0' num2str(j)]) & strcmp(subj_data.ModelTrialPm,event));
        
        %Compile into one vecotr
        s.subjects.(id).([lower(event) '_onset_all_runs']) = [s.subjects.(id).([lower(event) '_onset_all_runs']);  s.subjects.(id).([lower(event) '_onset']).(['run' num2str(j)])];
        s.subjects.(id).([lower(event) '_offset_all_runs']) = [s.subjects.(id).([lower(event) '_offset_all_runs']); s.subjects.(id).([lower(event) '_offset']).(['run' num2str(j)])];
        
        %Pull the infusion or feedback ratings depedning on event supplied
        if strcmpi(event,'INFUSION') || strcmpi(event,'FEEDBACK')
            s.subjects.(id).([lower(event) '_ratings']).(['run' num2str(j)])=...
                subj_data.Response(strcmp(subj_data.x_Participant,{id}) &...
                strcmp(subj_data.Run,['Run_0' num2str(j)]) & strcmp(subj_data.ModelTrialPm,event));
            
            
            %NOTE we currently censor the 0's, therefore in the future if
            %we are to include them we cannot just make the Nan response
            %values 0, just a warning!
        if any(isnan(s.subjects.(id).([lower(event) '_ratings']).(['run' num2str(j)])))
            warning('Nans were found in eitherfeedback or infusion response variable! Currently this is bad, will now set all Nans to zeros!')
            nan_index=isnan(s.subjects.(id).([lower(event) '_ratings']).(['run' num2str(j)]));
            s.subjects.(id).([lower(event) '_ratings']).(['run' num2str(j)])(nan_index)=0;
        end
            
            
        s.subjects.(id).([lower(event) '_ratings_all_runs']) = [s.subjects.(id).([lower(event) '_ratings_all_runs']);  s.subjects.(id).([lower(event) '_ratings']).(['run' num2str(j)])];
        end
    else
        %Censor out full runs if needed
        tmp_censor_data = subj_data.IncludeRun(strcmp(subj_data.x_Participant,id) &...
            strcmp(subj_data.Run,['Run_0' num2str(j)]));
        s.subjects.(id).run_censor_full = [s.subjects.(id).run_censor_full; strcmp(tmp_censor_data,'TRUE')];
    end
end

function s=get_detailed_event_times(s,id,rt_data,runs)
%Initialize
 s.subjects.(id).event_rt_onset_all_runs=[];
 s.subjects.(id).event_rt_offset_all_runs=[];

for j = 1:runs

        %Pull event on/offset times
        s.subjects.(id).event_rt_onset.(['run' num2str(j)])=...
            rt_data.Onset2(strcmp(rt_data.Participant,{id}) &...
            strcmp(rt_data.Run,['Run_0' num2str(j)])) + ...
            rt_data.RT(strcmp(rt_data.Participant,{id}) &...
            strcmp(rt_data.Run,['Run_0' num2str(j)])) ;
        
        %Worry about offsets later
%          s.subjects.(id).event_offset.(['run' num2str(j)]) =...
%             s.subjects.(id).event_onset.(['run' num2str(j)])...
%             + rt_data.Duration(strcmp(rt_data.Participant,{id}) &...
%             strcmp(rt_data.Run,['Run_0' num2str(j)]));
        
        %Compile into one vecotr
        s.subjects.(id).event_rt_onset_all_runs = [ s.subjects.(id).event_rt_onset_all_runs;  s.subjects.(id).event_rt_onset.(['run' num2str(j)])];
        %s.subjects.(id).event_offset_all_runs = [s.subjects.(id).event_offset_all_runs; s.subjects.(id).event_offset.(['run' num2str(j)])];

end
        %Might have to make Nan's 0's
        s.subjects.(id).any_rating = s.subjects.(id).event_rt_onset_all_runs;
        s.subjects.(id).any_rating_censor = ~isnan(s.subjects.(id).any_rating);
        s.subjects.(id).infusion_resp_rt = s.subjects.(id).event_rt_onset_all_runs(1:2:end);
        s.subjects.(id).feedback_resp_rt = s.subjects.(id).event_rt_onset_all_runs(2:2:end);
        
        
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