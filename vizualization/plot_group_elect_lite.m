function plot_group_elect_lite(data,task,cond_names, time)

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

figure('units', 'normalized', 'outerposition', [0 0 .4  .3])

try
    subplot(1,2,1)
    h = [];
    for i = 1:length(data.plot_data)
        lineprops.col{1} = plot_params.col(i,:);
        mseb(time,nanmean(data.plot_data{i}),nanste(data.plot_data{i}),lineprops,0);
        hold on
        h(i)=plot(time,nanmean(data.plot_data{i}),'LineWidth',plot_params.lw,'Color',plot_params.col(i,:));
        
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
    end
catch
end

try
    
    
     xlim(plot_params.xlim)
%     xlim([-3 .5])

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


subplot(1,2,2)

colors = hsv(size(data.plot_data{2},1));
for ie = 1:size(data.plot_data{2},1)
    hold on
    plot(time, data.plot_data{2}(ie,:), 'LineWidth', 2, 'Color', colors(ie, :))
    [v,idx] = max(data.plot_data{2}(ie,:));
    text(time(idx),v,[data.subj_name{ie} ' - ' data.label{ie}], 'Interpreter', 'none', 'Color', colors(ie, :))
end

if strcmp(task, 'Memoria')
    hold on
    time_events = [1 2.1 3.2 4.3];
    for ii = 1:length(time_events)
        plot([time_events(ii) time_events(ii)],ylim,'Color', [.5 .5 .5], 'LineWidth',1)
    end
elseif strcmp(task, 'Calculia')
    hold on
    time_events = [0 0.9 1.8 2.8 3.7 4.6];
    for ii = 1:length(time_events)
        plot([time_events(ii) time_events(ii)],ylim,'Color', [.5 .5 .5], 'LineWidth',1)
    end
end
xlim(plot_params.xlim)
xlabel(plot_params.xlabel)
ylabel(plot_params.ylabel)
set(gca,'fontsize',plot_params.textsize)
box off
xline(0, 'LineWidth',1.5);
yline(0, 'LineWidth',1.5);
set(gcf,'color','w');




end

