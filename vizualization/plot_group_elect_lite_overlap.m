function plot_group_elect_lite_overlap(data,task,cond_names, time, e)

% plot_data = [];
% for i = 1:length(data.wave)
%     for ic = 1:length(cond_names)
%         trials = find(strcmp(data.trialinfo_all{i}.(column), cond_names{ic}));
%         plot_data{ic}(i,:) = nanmean(data.wave{i}(trials,:));
%     end
% end

plot_params = genPlotParams(task,'timecourse');
lineprops = [];
lineprops.style= '-';
lineprops.width = plot_params.lw;
lineprops.edgestyle = '-';

load('cdcol_2018.mat')

plot_params.col = cbrewer2('Reds', 11);
plot_params.col(1:2,:) = [];
% plot_params.col = flip(plot_params.col);

try
    h = [];
    lineprops.col{1} = plot_params.col(e,:);
%     mseb(time,nanmean(data.plot_data),nanste(data.plot_data),lineprops,0);
    hold on
    h=plot(time,nanmean(data.plot_data),'LineWidth',plot_params.lw,'Color',plot_params.col(e,:));
    
    if strcmp(task, 'Memoria')
        time_events = [1 2.1 3.2 4.3];
        for ii = 1:length(time_events)
            plot([time_events(ii) time_events(ii)],ylim,'Color', [.5 .5 .5], 'LineWidth',1)
        end
    elseif strcmp(task, 'Calculia')
        time_events = [0 0.9 1.8 2.8 3.7 4.6];
        for ii = 1:length(time_events)
            plot([time_events(ii) time_events(ii)],ylim,'Color', [.5 .5 .5], 'LineWidth',1)
        end
    end
catch
end

try
    
    
    xlim(plot_params.xlim)
    xlabel(plot_params.xlabel)
    ylabel(plot_params.ylabel)
    set(gca,'fontsize',plot_params.textsize)
    box off
    xline(0, 'LineWidth',1.5);
    yline(0, 'LineWidth',1.5);
    set(gcf,'color','w');
    leg = legend(h,cond_names,'Location','Northeast', 'AutoUpdate','off', 'Interpreter', 'none');
    legend boxoff
    set(leg,'fontsize',14, 'Interpreter', 'none')
catch
end



end

