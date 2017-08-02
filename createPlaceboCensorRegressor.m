function b=createPlaceboCensorRegressor(b,num_blocks)
%This is a file of 1s and 0s, indicating which     
%time points are to be included (1) and which are  
%to be excluded (0).


frequency_scale_hz = 10;
scan_tr = 2; %According to 3dinfo -verb!
block_length = b.subjects.(b.id).block_length; 
%b.stim_NextOnsetTime=[b.stim_OnsetTime(2:end); b.stim_RTTime(end)];
b.stim_OnsetTime = b.subjects.(b.id).infusion_onset_all_runs*1000; %Convert to ms
%b.stim_NextOnsetTime=[b.subjects.(b.id).infusion_onset_all_runs(2:11); b.subjects.(b.id).baseline_offset_all_runs(end)]*1000; %Convert to ms
%b.missed_responses = ( b.stim_RT == 0 ); %Only included no responses
%b.missed_responses = ( b.chosen_stim == 999 ); %Includes wrong button presses and rt == 0
b.trials_to_censor = ~b.subjects.(b.id).run_censor_full; %DID I HAVE THIS BACKWARDS!? there was no '~' before 11/4

%I was actually censoring the 0 trials :( let's try what we originally
%thought was going on
b.trials_to_censor = zeros(length(b.subjects.(b.id).run_censor_full),1);


%Censor only the NANs in MDF file
b.trials_to_censor  =~b.subjects.(b.id).response_nan_censor_full;

%b.trials_to_censor(b.rewardVec==0) = 1; %Adding in the trials in which there was no reward as well.
b.stim_NextOnsetTime = []; %Hacky work around to get matrices to play nice
b.no_baseline_end_time = 394; %seconds (I don't really like hard coding this BE CAREFUL!)

for block_n = 1:num_blocks
        %Set up trial ranges
    trial_index_1 = b.trial_index(block_n);
    trial_index_2 = trial_index_1 + b.trials_per_block-1;
    
    if block_length > 200
        b.stim_NextOnsetTime=[b.stim_NextOnsetTime; b.subjects.(b.id).infusion_onset_all_runs(trial_index_1+1:trial_index_2)*1000; b.subjects.(b.id).baseline_offset_all_runs(trial_index_2)*1000]; %Convert to ms
    else
        b.stim_NextOnsetTime=[b.stim_NextOnsetTime; b.subjects.(b.id).infusion_onset_all_runs(trial_index_1+1:trial_index_2)*1000; b.no_baseline_end_time*1000]; %Convert to ms
    end
    bin_size = 1/frequency_scale_hz*1000; % convert Hz to msec
    epoch_window = b.stim_OnsetTime(trial_index_1:trial_index_2):bin_size:b.stim_OnsetTime(trial_index_1:trial_index_2)+scan_tr*block_length*1000;    
    event_beg = b.stim_OnsetTime(trial_index_1:trial_index_2); event_end = b.stim_NextOnsetTime(trial_index_1:trial_index_2);
    
    tmp_reg.(['regressors' num2str(block_n)]).to_censor = ...
        createSimpleRegressor(event_beg, event_end, epoch_window, b.trials_to_censor(trial_index_1:trial_index_2));
    tmp_reg.(['regressors' num2str(block_n)]).to_censor = ones(size(tmp_reg.(['regressors' num2str(block_n)]).to_censor)) - tmp_reg.(['regressors' num2str(block_n)]).to_censor;
    
    
    % NB: the first 5s are censored because they capture HRF to events
    % preceding the first trial
    tmp_reg.(['hrfreg' num2str(block_n)]).to_censor = ...
        gsresample( ...
            [zeros(50,1)' tmp_reg.(['regressors' num2str(block_n)]).to_censor(1:end-51)], ...
        10,1./scan_tr);
        
end


fnm = fieldnames(tmp_reg.regressors1)';
ct=1:length(fnm);
try
    b.hrf_regs.(fnm{ct}) = [tmp_reg.hrfreg1.(fnm{ct}) tmp_reg.hrfreg2.(fnm{ct}) tmp_reg.hrfreg3.(fnm{ct}) tmp_reg.hrfreg4.(fnm{ct}) tmp_reg.hrfreg5.(fnm{ct}) tmp_reg.hrfreg6.(fnm{ct})];
catch
    error('Censor is bad new bears, check the number of runs')
end
b.hrf_regs.to_censor = 1-(ceil(b.hrf_regs.to_censor));
b.hrf_regs.to_censor = ~b.hrf_regs.to_censor;
%b.hrf_regs.to_censor = b.hrf_regs.to_censor(1:length(b.hrf_regs.RT));



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