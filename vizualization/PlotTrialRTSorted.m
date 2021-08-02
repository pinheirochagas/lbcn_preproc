
function h = PlotTrialRTSorted(data,column,conds,plot_params)

% plots average timecourse for each condition, separately for each electrode
% INPUTS:
%       data: can be either freq x trial x time  or trial x time
%       column: column of data.trialinfo by which to sort trials for plotting
%       conds:  cell containing specific conditions to plot within column (default: all of the conditions within column)
%               can group multiple conds together by having a cell of cells
%               (e.g. conds = {{'math'},{'autobio','self-internal'}})
%       col:    colors to use for plotting each condition (otherwise will
%               generate randomly)
%       plot_params:    controls plot features (see genPlotParams.m script)

% will only work for non-spectral (i.e. Band) daata
datatype = 'NonSpec';

if isempty(plot_params)
    plot_params = genPlotParams(project_name,'timecourse');
end

ncategs = length(conds);

cmap = cbrewer2('RdBu');
cmap = cmap(end:-1:1,:); % flip colormap
% add white to deal with the nan post RT. 
cmap(2:end+1, :) = cmap;
cmap(1, :) = [1 1 1];


%%
winSize = floor(data.fsample*plot_params.sm);
gusWin= gausswin(winSize)/sum(gausswin(winSize));

plot_data = cell(1,ncategs);

if strcmp(datatype,'Spec')
    freq_inds = data.freqs >= plot_params.freq_range(1) & data.freqs <= plot_params.freq_range(2);
    data.wave = squeeze(nanmean(data.wave(freq_inds,:,:)));  % avg across freq. domain
end

data.wave = convn(data.wave,gusWin','same');

% Dirty fix CHANGE THAT!
% Dirty fix CHANGE THAT!
% Dirty fix CHANGE THAT!
data.wave(data.trialinfo.RT>5,:) = [];
data.trialinfo(data.trialinfo.RT>5,:) = [];
% Dirty fix CHANGE THAT!
% Dirty fix CHANGE THAT!
% Dirty fix CHANGE THAT!

% group data by conditions
% if plot_params.multielec
%     groupall = true;
% else
%     groupall = false;
% end

% if plotting single trials, include noisy trials so can plot in different color
groupall = false;
[grouped_trials_all,~] = groupConds(conds,data.trialinfo,column,'none',[],groupall);

[grouped_trials,cond_names] = groupConds(conds,data.trialinfo,column,plot_params.noise_method,plot_params.noise_fields_trials,groupall);
% if eliminating noisy trials, keep track of how many clean trials remain for
% each condition (and include in figure legend)
if (strcmp(plot_params.noise_method,'trials'))
    for gi = 1:length(cond_names)
        cond_names{gi} = [cond_names{gi},' (',num2str(length(grouped_trials{gi})),' of ',num2str(length(grouped_trials_all{gi})), ' trials)'];
    end
end

plot_data = cell(1,ncategs); % with noisy epochs excluded
plot_data_all = cell(1,ncategs); %including noisy epochs
for ci = 1:ncategs
    plot_data{ci} = data.wave(grouped_trials{ci},:);
    trial_data{ci} = data.trialinfo(grouped_trials{ci},:);
    plot_data_all{ci} = data.wave(grouped_trials_all{ci},:);
    trial_data_all{ci} = data.trialinfo(grouped_trials_all{ci},:);
    
    nstim(ci) = size(trial_data{ci}.allonsets,2);
    
    RTLock{ci}=nan(1,length(grouped_trials{ci}));
    for i = 1:length(grouped_trials{ci})
        RTLock{ci}(i) = trial_data{ci}.allonsets(i,nstim(ci))-trial_data{ci}.allonsets(i,1)+trial_data{ci}.RT(i);        
        postRTinds = find(data.time>RTLock{ci}(i));
        plot_data{ci}(i,postRTinds)=nan;
    end
    trial_data{ci}.RTLock = RTLock{ci}';
    % define which colomn to sort
    
%    if strcmp(plot_params.sort_column, 'RT')
%         [~,sortInds] = sortrows(trial_data{ci},{'RTLock', ''});        
%    elseif strcmp(plot_params.sort_column, 'operand2')        
%         [~,sortInds] = sortrows(trial_data{ci},{'operand2', 'RTLock'});        
%    end
   
    [~,sortInds] = sortrows(trial_data{ci},plot_params.sort_columns);        


    plot_data{ci}=plot_data{ci}(sortInds,:);
    trial_data{ci} = trial_data{ci}(sortInds,:);
end

% smooth and plot data
if plot_params.set_figure == 1
    figureDim = [0 0 .2 .4];
    figure('units', 'normalized', 'outerposition', figureDim)
else
end

% hold on
for ci = 1:ncategs
    %     if plot_params.single_trial
    if ncategs > 1
        subplot(ncategs,1,ci)
    else
    end
    clims = [-prctile(plot_data{ci}(:),97.5) prctile(plot_data{ci}(:),97.5)];
    if size(plot_data{ci},1)>1
%         h = imagesc(data.time,1:size(plot_data{ci},1),plot_data{ci},clims);
        imagesc(data.time,1:size(plot_data{ci},1),plot_data{ci},clims);
        colormap(cmap)
        hold on
        plot(trial_data{ci}.RTLock,1:size(plot_data{ci},1),'.', 'MarkerSize', 6, 'Color', 'k')
        
        % Plot stimuli list
        if plot_params.plot_slist == 1
            % Here configure which field you wanna plot
            for i = 1:size(plot_data{ci},1)
                if strcmp(trial_data{ci}.keys(i),'1')
                    text(5,i,[trial_data{ci}.wlist{i},' (True)'])
                elseif strcmp(trial_data{ci}.keys(i),'2')
                    text(5,i,[trial_data{ci}.wlist{i},' (False)'])
                else
                    text(5,i,trial_data{ci}.wlist{i})
                end
            end
            
        else
        end
        
        
        plot([0 0],ylim,'k-','LineWidth',2)
%         title(cond_names{ci})
        xlabel(plot_params.xlabel)
        ylabel('RT-sorted trials')
        set(gca,'fontsize',plot_params.textsize)
        set(gca, 'xlim', plot_params.xlim)
        
        if size(data.trialinfo.allonsets,2) > 1
            time_events = cumsum(nanmean(diff(data.trialinfo.allonsets(:,1:nstim(ci)),1,2)));
            for i = 1:length(time_events)
                plot([time_events(i) time_events(i)],ylim,'k-','LineWidth',1)
            end
        end
        
        % Add horizontal lines to separate conditions
%         if ~iscell(trial_data{ci}.(plot_params.group_conds))
%             trial_data{ci}.(plot_params.group_conds) = cellstr(num2str(trial_data{ci}.(plot_params.group_conds)));
%         else
%         end
%         tab_tmp = tabulate(trial_data{ci}.(plot_params.group_conds));
%         % FIX this, now only working for conditions with 2 names! 
%         for i = 1:size(tab_tmp,1)
%             if i == 1
%                 line_loc = cell2mat(tab_tmp(i,2));
%                 if line_loc == 1
%                     line_loc_lab = line_loc;
%                 else
%                     line_loc_lab = line_loc/2;
%                 end
%             else
%                 line_loc = cell2mat(tab_tmp(i,2)) + cell2mat(tab_tmp(i-1,2));
%                     if line_loc == 1
%                         line_loc_lab = (cell2mat(tab_tmp(i,2))) + cell2mat(tab_tmp(i-1,2));
%                     else
%                         line_loc_lab = (cell2mat(tab_tmp(i,2)))/2 + cell2mat(tab_tmp(i-1,2));
%                     end
%             end
%             plot(xlim,[line_loc+0.5 line_loc+0.5],'k-','LineWidth',3)
%             lin_locs(i) = line_loc_lab;
%         end
%         yticks(lin_locs)
%         yticklabels(tab_tmp(:,1))
%         
        
        
        box off
        hold on
%         colorbar('south')
    else
        h = imagesc(0);
    end
end


set(gcf,'color','w')



