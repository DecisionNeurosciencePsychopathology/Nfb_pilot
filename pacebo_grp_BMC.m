%  runs group model comparison on SCEPTIC log model evidence data from sceptic_fit_group_vba.m

% clear variables;
close all;

%Set path to data storage on bek
file_path = '/Users/localadmin/Dropbox/data_projects/placebo_marta/';
%data_file = 'group_model_comparison_L';
%data_file = 'grp_kalman_16_nbasis_40_nsteps_1_uasversion_with_q';
%data_file = 'grp_L_kalman_logistic16_nbasis40_nsteps1_uaversion';
% data_file = 'L9_fixed_ksoft_kprocnoise_kuvsum_ksigvol_kuvlog_fixeduv_q_fixeddecay.mat';
data_file = 'L2';

%Load in proper data file
group_data=load([file_path data_file]);

%Do you want to run models in families?
family_flag =0;

L = group_data.L2;
%sceptic.L = group_data.L9;
%sceptic.modelnames = models;

modelnames = group_data.modelnames;
%modelnames = {'fixed' 'kalmanSoftmax' 'kalmanProcessnoise' 'kalmanUVsum' 'kalmanSigmavolatility' 'kalmanLogistic' 'fixedUV' 'Qstep' 'fixedDecay'};
options.modelnames = modelnames;

if family_flag
    % %% Define 'fixed' as the first family, 'kalman', as the second, q third
    idx = strfind(modelnames, 'fixed');
    %options.families{1} = find(strcmpi(modelnames, 'fixed')); 
    options.families{1} = find(not(cellfun('isempty', idx))); 
    idx = strfind(modelnames, 'kalman');
    options.families{2} = find(not(cellfun('isempty', idx)));
    options.families{3} = find(strcmpi(modelnames, 'Qstep'));
    [posterior,out] = VBA_groupBMC(sceptic.L,options);
end



[posterior,out] = VBA_groupBMC(L);

%% Clean up model strings to make graphic look pretty
% modelnames = cellfun(@strrep, modelnames', repmat({'_'},length(modelnames),1), repmat({' '},length(modelnames),1), 'UniformOutput', false);
for i = 1:length(out.options.handles.ha)-2
    xlabel(out.options.handles.ha(i),'models')
    set(out.options.handles.ha(i),'xtick',1:length(modelnames), 'XTickLabel',char(modelnames))
end


%% save output figure
h = figure(1);
file_str=input('What do you want save the figure as? ', 's');
saveas(h,[file_path file_str])
%saveas(h,[file_path 'VBA_group_BMC_nbasis16_nsteps_kalman_40_1_uaversion_q'])

%savefig(h,sprintf('VBA_group_BMC_nbasis4_nsteps10_m1=aversion_supertight_tau_m2=aversion_tight_tau_m3_no_aversion'))
%savefig(h,sprintf('VBA_group_BMC_nbasis4_nsteps10_m1=aversion_supertight_tau_m2=aversion_tight_tau_m3_no_aversion'))
