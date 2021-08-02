function plotPAC(PAC,conds,phase_elec,amp_elec)

if isempty(amp_elec) % assume same as phase_elec
    amp_elec = phase_elec;
end

p_ticks = 1:3:length(PAC.phase_freq);
a_ticks = 1:3:length(PAC.amp_freq);

p_label = cell(1,length(p_ticks));
a_label = cell(1,length(a_ticks));

for i = 1:length(p_ticks)
    p_label{i} = num2str(round(PAC.phase_freq(p_ticks(i))));
end

for i = 1:length(a_ticks)
    a_label{i} = num2str(round(PAC.amp_freq(a_ticks(i))));
end

figure('Position',[200 200 400*length(conds) 500])


for ci = 1:length(conds)
    data_plot = abs(PAC.(conds{ci}).(['p',phase_elec,'_a',amp_elec]))';
    max_caxis(ci) = prctile(data_plot(:), 97)
end

for ci = 1:length(conds)
    subplot(1,length(conds),ci)
    data_plot = abs(PAC.(conds{ci}).(['p',phase_elec,'_a',amp_elec]))';
    imagesc(data_plot)
    axis xy
    hold on
    set(gca,'XTick',p_ticks)
    set(gca,'YTick',a_ticks)
    set(gca,'XTickLabel',p_label)
    set(gca,'YTickLabel',a_label)
    set(gca,'FontSize',20)
    colormap(cbrewer2('Reds'))
    caxis([0 max(max_caxis)])

    % title('Math Trials')
    % colormap(jet)
    set(gcf,'color','w')
    xlabel('Phase freq. (Hz)')
    ylabel('Amp freq. (Hz)')
    colorbar
    title(conds{ci})
    
end
suptitle(['Phase: ',phase_elec,'; Amp: ',amp_elec])