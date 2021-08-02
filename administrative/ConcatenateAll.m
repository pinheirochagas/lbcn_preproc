function data_all = ConcatenateAll(sbj_name, project_name, block_names, dirs,elecs, datatype, freq_band, locktype, concat_params)
%% Define electrodes

% Load subjectVar
load([dirs.original_data filesep sbj_name filesep 'subjVar_' sbj_name '.mat'])

if strcmp(sbj_name, 'S20_151_HT') && strcmp(project_name, 'MMR') 
%     original = load('/Volumes/LBCN8T/Stanford/data/neuralData/originalData/S20_151_HT/global_Calculia_production_S20_151_HT_E20-540_0011.mat');
%     mmr = load([dirs.original_data, filesep, sbj_name,'/global_',project_name,'_',sbj_name,'_',block_names{1},'.mat']);
%     elecs = 1:mmr.globalVar.nchan;
%     [~,idx] = setdiff(mmr.globalVar.channame, original.globalVar.channame);
%     elecs(idx) = [];
else
    load([dirs.original_data, filesep, sbj_name,'/global_',project_name,'_',sbj_name,'_',block_names{1},'.mat'])
    
    if isempty(elecs)
        % load globalVar (just to get ref electrode, # electrodes)
        load([dirs.original_data, filesep, sbj_name,'/global_',project_name,'_',sbj_name,'_',block_names{1},'.mat'])
        elecs = setdiff(1:size(subjVar.elinfo,1),globalVar.refChan);
        elecs = 1:size(subjVar.elinfo,1);
    end
end



if isempty(concat_params)
    concat_params = genConcatParams(false); % default: no downsampling
end


if ~isfield(subjVar, 'elinfo')
    data_format = GetFSdataFormat(sbj_name, 'Stanford');
    subjVar = CreateSubjVar(sbj_name);
else
end


if strcmp(datatype,'Spec')
    tdim = 4; % time dimension after concatenating
    tag = [locktype,'lock_bl_corr']; % specifies type of data to load
    
elseif strcmp(datatype,'Band')
    tdim = 3;
    if strcmp(project_name, 'GradCPT')
        tag = [locktype,'lock']; % specifies type of data to load
    else
        tag = [locktype,'lock_bl_corr']; % specifies type of data to load
    end
elseif strcmp(datatype,'CAR')
    tdim = 3;
    tag = [locktype,'lock']; % specifies type of data to load
end
%% loop through electrodes
data_all.trialinfo = [];
concatfields = {'wave'}; % type of data to concatenate
for ei = 1:length(elecs)
    el = elecs(ei);
    
    data_bn = concatBlocks(sbj_name,project_name,block_names,dirs,el,freq_band,datatype,concatfields,tag);
    elecnans(ei) = sum(sum(isnan(data_bn.wave)));
    
    % Add extra spike detector
    % if task is active, set sample after RT to nan
    if (~strcmp(project_name, 'Rest') && ~strcmp(project_name, 'EglyDriver'))
        for iout = 1:size(data_bn.wave, 1)
            data_tmp(iout,:) = data_bn.wave(iout,:);
            data_tmp(iout,data_bn.time>data_bn.trialinfo.RT(iout)) = nan;
        end
        max_val = max(data_tmp, [], 2);
        if strcmp(subjVar.elinfo.DK_lobe{el}, 'Occipital')
            thrhold = 4;
        else
            thrhold = 3;
        end
        data_bn.trialinfo.spike_hfb = zscore(log(max_val))>thrhold ;
    end
    %% Add extra spike detector
%     for it = 1:size(data_bn.wave,1)
%         data_tmp_out = data_tmp;
%         max_val = max(data_tmp_out, [], 2);
% %         var_val = nanvar(data_tmp_out, [], 2)
% %         boxplot(max_val)
% %         plot(data_tmp_out')
% % %         out_max = isoutlier(max(data_tmp_out, [], 2), 'gesd')
% % %         f_out_max = find(max_val == max(max_val(out_max ==1)));
% %         mad_m = mad_median(max_val)
% % 
% %         f_out_max = find(max_val>median(max_val)*3);
%         if strcmp(subjVar.elinfo.DK_lobe{el}, 'Occipital')
%             thrhold = 4;
%         else
%             thrhold = 3;
%         end
%         
%         f_out_max = find(zscore(log(max_val))>thrhold);
% %         f_out_var = find(zscore(log(var_val))>thrhold);
%         
%         data_tmp_out(f_out_max,:) = [];
%         
%                 plot(data_tmp_out')
% 
%         
%         
%         s = data_bn.wave(it,:);
%         data_bn.trialinfo.spike_hfb(it) = sum(s>prctile(s,99)*5);        
%     end
    
    if strcmp(concat_params.noise_method,'timepts')
        data_bn = removeBadTimepts(data_bn,concat_params.noise_fields_timepts);
    elseif strcmp(concat_params.noise_method,'none')
        
    elseif strcmp(concat_params.noise_method,'trials')
        bad_trials = [];
        for i = 1:length(concat_params.noise_fields_trials)
            bad_trials = union(bad_trials,find(data_bn.trialinfo.(concat_params.noise_fields_trials{i})));
            bad_trials = reshape(bad_trials,length(bad_trials),1);
        end
        % dont do this... seriously... 
%         if strcmp(datatype,'Band') || strcmp(datatype,'CAR')
%             data_bn.wave(bad_trials,:) = NaN;
%         else
%             data_bn(bad_trials,:,:) = NaN;
%         end
        % Define bad channels as a function of number of good trials.
               
    end
    
    if concat_params.decimate % smooth and downsample (optional)
        if data_bn.fsample == 999
            ds_rate = 2; % FIX THIS, it assumes fs = 1000Hz.
            data_all.fsample =concat_params.fs_targ ;
            data_all.time = data_bn.time(1:ds_rate:end);
        else
            ds_rate = floor(data_bn.fsample/concat_params.fs_targ); % FIX THIS, it assumes fs = 1000Hz.
            data_all.fsample = data_bn.fsample/ds_rate;
            data_all.time = data_bn.time(1:ds_rate:end);
        end

        %         if concat_params.sm_win > 0 % if smoothing first
        %             winSize = floor(data_bn.fsample*concat_params.sm_win);
        %             gusWin= gausswin(winSize)/sum(gausswin(winSize));
        %             data_bn.wave = convn(data_bn.wave,shiftdim(gusWin,-tdim),'same'); % convolve data w/gaussian along time dimension
        %         else
        %         end
        % downsample
        if strcmp(datatype,'Band') || strcmp(datatype,'CAR')
            data_bn.wave = data_bn.wave(:,1:ds_rate:end);
        elseif strcmp(datatype,'Spec')
            data_bn.wave = data_bn.wave(:,:,1:ds_rate:end);
        end
        
    else
        data_all.time = data_bn.time;
        data_all.fsample = data_bn.fsample;
    end
    
    % Concatenate all subjects all trials
    if strcmp(datatype,'Band') || strcmp(datatype,'CAR')
        data_all.wave(:,ei,:) = data_bn.wave;
    elseif strcmp(datatype,'Spec')
        data_all.wave(:,:,ei,:) = data_bn.wave;
    end
    
    %     data_all.label = data_bn.label;
    
    data_all.trialinfo = [data_bn.trialinfo];
    data_all.trialinfo_all{el} = [data_bn.trialinfo];
    %     data_all.labels{ei} = data_bn.label;
    disp(['concatenating elec ',num2str(el)])
    data_all.label(ei) = subjVar.elinfo.FS_label(ei); % just modified that
    data_all.chan_num(ei) = subjVar.elinfo.chan_num(ei); % just modified that
    if strcmp(concat_params.noise_method,'trials')
            data_all.bad_trials{el} = bad_trials;
    else
    end
end


%% Correct for the actual recorded channels
% if size(data_all.wave, 2) > 1
%     if size(data_all.wave, 2) ~= size(subjVar.elinfo,1)
%         if ~isempty(str2num(globalVar.channame{1}))
%             disp('Correctig for the actual recorded channels')
%             nchan_fs = size(subjVar.elinfo,1);
%             in_chan_cmp = false(1,nchan_fs);
%             for i = 1:nchan_fs
%                 in_chan_cmp(i) = ismember(subjVar.elinfo.FS_label(i),globalVar.channame);
%             end
%             % If TDT, channels are all in freesurfer?
%         else
%             nchan_cmp = size(globalVar.channame,2);
%             in_fs = false(1,nchan_cmp);
%             for i = 1:nchan_cmp
%                 in_fs(i) = ismember(globalVar.channame(i),subjVar.elinfo.FS_label);
%             end
%             data_all.wave = data_all.wave(:, in_fs, :);
%             data_all.trialinfo_all = data_all.trialinfo_all(in_fs);
%             data_all.label = subjVar.elinfo.FS_label;
%         end
%     else
%         data_all.label = subjVar.elinfo.FS_label;
%     end
% else
% end




% Concatenate bad channels
badChan = [];
for bi = 1:length(block_names)
    % Load globalVar
    load([dirs.data_root,'/OriginalData/',sbj_name,'/global_',project_name,'_',sbj_name,'_',block_names{bi},'.mat'])
    badChan = [globalVar.badChan badChan];
end
% Finalize
% data_all.time = data.time;
% data_all.fsample = data.fsample;
data_all.badChan = unique(badChan);
data_all.project_name = project_name;


%% Exclude or not channels
if isfield(concat_params, 'exclude_nan_chan') && concat_params.exclude_nan_chan
    % check for channels with nan
    nan_channel = [];
    for i = 1:size(data_all.wave,2)
        for ii = 1:size(data_all.wave,1)
            n_nan = sum(isnan(data_all.wave(ii,i,:)));
            if n_nan > 0
                if n_nan > floor(size(data_all.wave,3)*0.1) % if number of nans exceed more than 10% of the trial time points, simply replace it with zeros.
                    data_all.wave(ii,i,:) = zeros(size(data_all.wave,3),1)';
                    n_trial_nan(i,ii) = 1;
                else % if less than 10% interpolate
                    sprintf('interpolating nan values of chan %s, trial%s.png', num2str(i), num2str(ii))
                    data_all.wave(ii,i,:) = fillmissing(data_all.wave(ii,i,:),'linear');
                    n_trial_nan(i,ii) = 0;
                end
            else
                
            end
            nan_channel(i) = sum(sum(isnan(data_all.wave(:,i,:))));
        end
    end
    % Exclude channels with more than 10% of nan trials.
    sum_n_trial_nan = sum(n_trial_nan,2);
    good_chans = logical(sum_n_trial_nan< size(data_all.wave,1)*0.05);
    data_all.wave = data_all.wave(:,good_chans,:);
    data_all.trialinfo_all = data_all.trialinfo_all{good_chans};
    data_all.label = data_all.label{good_chans};
end

if strcmp(concat_params.data_format, 'fieldtrip_raw')
    
    %     Reshape to trials and then channelsXtimes
    if size(data_all.wave,2) == 1
        for i = 1:size(data_all.wave,1)
            waveOrg{i} = squeeze(data_all.wave(i,:,:))';
        end
    else
        for i = 1:size(data_all.wave,1)
            waveOrg{i} = squeeze(data_all.wave(i,:,:));
        end
    end
    
    data_all.trial =  waveOrg;
    %     data_all.trial =  data_all.wave;
    
    data_all = rmfield(data_all, 'wave');
    data_all = rmfield(data_all, 'badChan');
    data_all = rmfield(data_all, 'project_name');
    data_all = rmfield(data_all, 'trialinfo_all');
    
    %% Filterout bad trials 
    bad_trials = find(data_all.trialinfo.bad_epochs == 99999999999); % temp
    data_all.trialinfo(bad_trials,:)=[];
    data_all.trial(bad_trials)=[];

    trialinfo =  data_all.trialinfo.(concat_params.trialinfo_var); % be carefull with that, simple solution for EglyDriver, only including one column
    %     trialinfo = data_all.trialinfo; % be carefull with that, simple solution for EglyDriver, only including one column
    %     trialinfo = [data_all.trialinfo.RT data_all.trialinfo.isCalc]; % be carefull with that, simple solution for EglyDriver, only including one column
    time = data_all.time;
    ntrials = size(data_all.trialinfo,1);
    data_all =  rmfield(data_all, 'trialinfo');
    data_all =  rmfield(data_all, 'time');
    for i = 1:ntrials
        %           data_all.trialinfo{i} = trialinfo;
        data_all.time{i} = time;
    end
    data_all.trialinfo = trialinfo;
    %      data_all.trialinfo = trialinfo;
    %     data_all.time = data_all.time;
    %     data_all.label = data_all.label';
    
elseif strcmp(concat_params.data_format, 'fieldtrip_fq')
    
    data_all.trial = data_all.wave;
    data_all.trialinfo = data_all.trialinfo.(concat_params.trialinfo_var); % be carefull with that, simple solution for EglyDriver, only including one column
    data_all = rmfield(data_all, 'wave');
    data_all = rmfield(data_all, 'badChan');
    data_all = rmfield(data_all, 'project_name');
    data_all = rmfield(data_all, 'trialinfo_all');
    
    data_all = rmfield(data_all, 'bad_trials');

    data_all.dimord = 'rpt_chan_time';
else
    
end

end
