
%% |PIPELINE FOR THE MATH PAPER|

%% Define paths to directories
[server_root, comp_root, code_root] = AddPaths('Pedro_iMAC');
dirs = InitializeDirs(' ', ' ', comp_root, server_root, code_root); % 'Pedro_NeuroSpin2T'


%% Tasks:
% * *Localizers:* VTC, Scrambled, AllCateg, Logo, 7Heaven
% * *Calculation Simultaneous:* MMR, UCLA, MFA
% * *Calculation Sequential:* Calculia, Memoria

%% Paper folder
result_dir = '/Users/pinheirochagas/Pedro/Stanford/papers/spatiotempoal_dynamics_math/results/';
figure_dir = '/Users/pinheirochagas/Pedro/Stanford/papers/spatiotempoal_dynamics_math/figures/';
dirs = InitializeDirs(' ', ' ', comp_root, server_root, code_root); % 'Pedro_NeuroSpin2T'
dirs.result_dir = result_dir;
%% Define final cohorts
% Read the google sheets
[DOCID,GID] = getGoogleSheetInfo('math_network','cohort');
sinfo = GetGoogleSpreadsheet(DOCID, GID);
subject_names = sinfo.sbj_name;
sinfo = sinfo(strcmp(sinfo.subjVar, '1'),:);
% Filer for usable data
sinfo = sinfo(strcmp(sinfo.behavior, '1'),:);
subjects = unique(sinfo.sbj_name);
subjects_unique = unique(subjects);

% if offline
sinfo = readtable('/Users/pinheirochagas/Pedro/Stanford/code/lbcn_personal/projects/math_network/math_network_cohort.csv')
subject_names = sinfo.sbj_name;
sinfo = sinfo(sinfo.subjVar == 1,:);
% Filer for usable data
sinfo = sinfo(sinfo.behavior == 1,:);
subjects = unique(sinfo.sbj_name);
subjects_unique = unique(subjects);



%% Plot coverage
outdir = '/Volumes/LBCN8T/Stanford/data/Results/coverage/';
% Filer for subjVar data


% Vizualize task per subject
vizualize_task_subject_coverage(sinfo, 'task_group')
savePNG(gcf, 600, [figure_dir, 'tasks_break_group_down_subjects.png'])

% Plot full coverage
vars = {'LvsR','MNI_coord', 'WMvsGM', 'sEEG_ECoG', 'DK_lobe', 'Yeo7', 'Yeo17', 'DK_long_josef'};
subjVar_all_all = ConcatSubjVars(subjects, dirs, vars);
subjVar_all = subjVar_all_all(strcmp(subjVar_all_all.WMvsGM, 'GM') | strcmp(subjVar_all_all.WMvsGM, 'WM'), :);
sort_tabulate(subjVar_all.WMvsGM)
sort_tabulate(subjVar_all.sEEG_ECoG)
sort_tabulate(subjVar_all.DK_lobe_generic)
sort_tabulate(subjVar_all.Yeo7)
sort_tabulate(subjVar_all.Yeo17)

data_ecog = subjVar_all(strcmp(subjVar_all.sEEG_ECoG, 'ECoG'), :);
data_seeg = subjVar_all(strcmp(subjVar_all.sEEG_ECoG, 'sEEG'), :);
sort_tabulate(data_ecog.WMvsGM)
sort_tabulate(data_ecog.LvsR)
sort_tabulate(data_seeg.WMvsGM)
sort_tabulate(data_seeg.LvsR)


subjVar_all(strcmp(subjVar_all.DK_long_josef, 'OUT OF BRAIN'), :) = [];
subjVar_all(strcmp(subjVar_all.DK_long_josef, 'EXCLUDE'), :) = [];

subjVars = struct
subjVars.elinfo = subjVar_all;

cfg = getPlotCoverageCFG('full');
cfg.views ={'lateral', 'lateral', 'ventral', 'ventral'}
cfg.hemis = {'left', 'right', 'left', 'right'}
cfg.subplots = [2, 2]
cfg.alpha = 0.5
PlotModulation(dirs, subjVars, cfg)



% Plot coverage by task group
col_group = 'task_group'; % or task_group
task_group = unique(sinfo.(col_group));
cols = hsv(length(task_group));
for i = 1:length(task_group)
    disp(task_group{i})
    sinfo_tmp = sinfo(strcmp(sinfo.(col_group), task_group{i}),:);
    subjects = unique(sinfo_tmp.sbj_name);
    subjVar_all = ConcatSubjVars(subjects, dirs);
    cfg = getPlotCoverageCFG('tasks_group');
    cfg.MarkerColor = cols(i,:);
    PlotModulation(dirs, subjVar_all, cfg)
end





%% Univariate Selectivity
tag = 'stim';
tasks = unique(sinfo.task);
tasks = {'MMR', 'UCLA', 'Memoria', 'Calculia'};
dirs = InitializeDirs(tasks{1}, sinfo.sbj_name{1}, comp_root, server_root, code_root); % 'Pedro_NeuroSpin2T'
dirs.result_dir = result_dir;
for it = 1:length(tasks)
    task = tasks{it};
    sinfo_tmp = sinfo(strcmp(sinfo.task, task),:);
    parfor i = 1:size(sinfo_tmp,1)
        ElecSelectivityAll(sinfo_tmp.sbj_name{i}, dirs, task, 'stim', 'Band', 'HFB')
    end
end

      
%% Proportion selectivity Calc simultaneous


vars = {'chan_num', 'FS_label', 'LvsR','MNI_coord', 'WMvsGM', 'sEEG_ECoG', 'DK_lobe', 'Yeo7', 'Yeo17', 'DK_long_josef', ...
    'elect_select', 'act_deact_cond1', 'act_deact_cond2', 'sc1c2_FDR', 'sc1b1_FDR' , 'sc2b2_FDR', ...
    'sc1c2_Pperm', 'sc1b1_Pperm', 'sc2b2_Pperm', 'sc1c2_tstat', 'sc1b1_tstat', 'sc2b2_tstat'};

task = 'MMR';
sinfo_MMR = sinfo(strcmp(sinfo.task, task),:);
el_selectivity_MMR = concat_elect_select(sinfo_MMR.sbj_name, task, dirs, vars);
task = 'UCLA';
sinfo_UCLA = sinfo(strcmp(sinfo.task, task),:);
el_selectivity_UCLA = concat_elect_select(sinfo_UCLA.sbj_name, task, dirs, vars);

el_selectivity_calc_sim = [el_selectivity_MMR;el_selectivity_UCLA]
el_selectivity_calc_sim = el_selectivity_calc_sim(strcmp(el_selectivity_calc_sim.WMvsGM, 'GM') | strcmp(el_selectivity_calc_sim.WMvsGM, 'WM'), :);
el_selectivity_calc_sim = el_selectivity_calc_sim(~strcmp(el_selectivity_calc_sim.Yeo7, 'FreeSurfer_Defined_Medial_Wall'),:)



selectivities = {{'math only'}, {'math selective'}, {'math deact'}, {'no selectivity'}, {'episodic only', 'autobio only'}, {'episodic selective', 'autobio selective'}, {'autobio deact',  'episodic deact'}, {'math and episodic', 'math and autobio'}};
columns = {'Yeo7', 'DK_lobe'};
figure('units', 'normalized', 'outerposition', [0 0 1 1])
for ic = 1:length(columns)
    for i = 1:length(selectivities)
        selectivity = selectivities{i};
        column = columns{ic};
        subplot(2,4,i)
        if strcmp(selectivity, 'math deact') == 1
            el_temp = el_selectivity_calc_sim(el_selectivity_calc_sim.act_deact_cond1 == -1,:);
        elseif contains(selectivity, {'autobio deact',  'episodic deact'}) == 2
            el_temp = el_selectivity_calc_sim(el_selectivity_calc_sim.act_deact_cond2 == -1,:);
        else
            el_temp = el_selectivity_calc_sim(contains(el_selectivity_calc_sim.elect_select, selectivity),:);
        end
        plot_frequency(el_temp, column, 'ascend', 'horizontal')
        title(selectivity)
    end
    savePNG(gcf, 300, [figure_dir, ['selectivity ', column, '.png']])
end


selectivities = {{'math only'},{'episodic only', 'autobio only'}};
columns = {'DK_long_josef'};
figure('units', 'normalized', 'outerposition', [0 0 1 1])
for ic = 1:length(columns)
    for i = 1:length(selectivities)
        selectivity = selectivities{i};
        column = columns{ic};
        subplot(1,2,i)
        if strcmp(selectivity, 'math deact') == 1
            el_temp = el_selectivity_calc_sim(el_selectivity_calc_sim.act_deact_cond1 == -1,:);
        elseif contains(selectivity, {'autobio deact',  'episodic deact'}) == 2
            el_temp = el_selectivity_calc_sim(el_selectivity_calc_sim.act_deact_cond2 == -1,:);
        else
            el_temp = el_selectivity_calc_sim(contains(el_selectivity_calc_sim.elect_select, selectivity),:);
        end
        plot_frequency(el_temp, column, 'ascend', 'horizontal')
        title(selectivity)
    end
    savePNG(gcf, 300, [figure_dir, ['selectivity ', column, '.png']])
end

%% Plot coverage selectivity

% Plot coverage by task group
el_selectivity_only = el_selectivity_calc_sim(contains(el_selectivity_calc_sim.elect_select, 'only'),:);


cfg = getPlotCoverageCFG('tasks_group');
cfg.MarkerSize = 10;
cfg.alpha = 0.5;
cfg.views = {'lateral', 'lateral', 'ventral', 'ventral'};
cfg.hemis = {'left', 'right', 'right', 'left'};
cfg.subplots = [2,2];
cfg.figureDim = [0 0 1 1];

load('cdcol_2018.mat')
for i = 1:size(el_selectivity_only,1)
    if contains(el_selectivity_only.elect_select{i}, 'math') == 1
        cfg.MarkerColor(i,:) = cdcol.light_cadmium_red;
    else
        cfg.MarkerColor(i,:) = cdcol.sapphire_blue;
    end
end
PlotModulation(dirs, el_selectivity_only, cfg)
savePNG(gcf, 600, [figure_dir, 'MMR_selectivity_brain_only.png'])



el_selectivity_selective = el_selectivity_calc_sim(contains(el_selectivity_calc_sim.elect_select, 'selective'),:);
load('cdcol_2018.mat')
for i = 1:size(el_selectivity_selective,1)
    if contains(el_selectivity_selective.elect_select{i}, 'math selective') == 1
        cfg.MarkerColor(i,:) = cdcol.light_cadmium_red;
    else
        cfg.MarkerColor(i,:) = cdcol.sapphire_blue;
    end
end
PlotModulation(dirs, el_selectivity_selective, cfg)
savePNG(gcf, 300, [figure_dir, 'MMR_selectivity_brain_selective.png'])



el_selectivity_math_deact = el_selectivity_calc_sim(el_selectivity_calc_sim.act_deact_cond1 == -1,:);
el_selectivity_math_deact(~strcmp(el_selectivity_math_deact.Yeo7, 'Somatomotor'),:)
load('cdcol_2018.mat')
for i = 1:size(el_selectivity_math_deact,1)
    if contains(el_selectivity_math_deact.elect_select{i}, 'math selective') == 1
        cfg.MarkerColor(i,:) = cdcol.light_cadmium_red;
    else
        cfg.MarkerColor(i,:) = cdcol.sapphire_blue;
    end
end
PlotModulation(dirs, el_selectivity_math_deact, cfg)
savePNG(gcf, 300, [figure_dir, 'MMR_el_selectivity_math_deact.png'])





%% MEMORIA

task = 'Memoria';
sinfo_Memoria = sinfo(strcmp(sinfo.task, task),:);
el_selectivity_Memoria = concat_elect_select(sinfo_Memoria.sbj_name, task, dirs, vars);
el_selectivity_Memoria = el_selectivity_Memoria(strcmp(el_selectivity_Memoria.WMvsGM, 'GM') | strcmp(el_selectivity_Memoria.WMvsGM, 'WM'), :);
el_selectivity_Memoria = el_selectivity_Memoria(~strcmp(el_selectivity_Memoria.Yeo7, 'FreeSurfer_Defined_Medial_Wall'),:)

selectivities = {{'math only'}, {'math selective'}, {'math deact'}, {'no selectivity'}, {'autobio only'}, {'autobio selective'}, {'autobio deact'}, {'math and autobio'}};
columns = {'Yeo7', 'DK_lobe'};
figure('units', 'normalized', 'outerposition', [0 0 1 1])
for ic = 1:length(columns)
    for i = 1:length(selectivities)
        selectivity = selectivities{i};
        column = columns{ic};
        subplot(2,4,i)
        if strcmp(selectivity, 'math deact') == 1
            el_temp = el_selectivity_Memoria(el_selectivity_Memoria.act_deact_cond1 == -1,:);
        elseif contains(selectivity, {'autobio deact',  'episodic deact'}) == 2
            el_temp = el_selectivity_Memoria(el_selectivity_Memoria.act_deact_cond2 == -1,:);
        else
            el_temp = el_selectivity_Memoria(contains(el_selectivity_Memoria.elect_select, selectivity),:);
        end
        plot_frequency(el_temp, column, 'ascend', 'horizontal')
        title(selectivity)
    end
    savePNG(gcf, 300, [figure_dir, ['Memoria_selectivity ', column, '.png']])
end


selectivities = {{'math only'},{'episodic only', 'autobio only'}};
columns = {'DK_long_josef'};
figure('units', 'normalized', 'outerposition', [0 0 1 1])
for ic = 1:length(columns)
    for i = 1:length(selectivities)
        selectivity = selectivities{i};
        column = columns{ic};
        subplot(1,2,i)
        if strcmp(selectivity, 'math deact') == 1
            el_temp = el_selectivity_Memoria(el_selectivity_Memoria.act_deact_cond1 == -1,:);
        elseif contains(selectivity, {'autobio deact',  'episodic deact'}) == 2
            el_temp = el_selectivity_Memoria(el_selectivity_Memoria.act_deact_cond2 == -1,:);
        else
            el_temp = el_selectivity_Memoria(contains(el_selectivity_Memoria.elect_select, selectivity),:);
        end
        plot_frequency(el_temp, column, 'ascend', 'horizontal')
        title(selectivity)
    end
    savePNG(gcf, 300, [figure_dir, ['Memoria_selectivity ', column, '.png']])
end



el_selectivity_only = el_selectivity_Memoria(contains(el_selectivity_Memoria.elect_select, 'only'),:);


cfg = getPlotCoverageCFG('tasks_group');
cfg.MarkerSize = 10;
cfg.alpha = 0.5;
cfg.views = {'lateral', 'lateral', 'ventral', 'ventral'};
cfg.hemis = {'left', 'right', 'right', 'left'};
cfg.subplots = [2,2];
cfg.figureDim = [0 0 1 1];

load('cdcol_2018.mat')
for i = 1:size(el_selectivity_only,1)
    if contains(el_selectivity_only.elect_select{i}, 'math') == 1
        cfg.MarkerColor(i,:) = cdcol.light_cadmium_red;
    else
        cfg.MarkerColor(i,:) = cdcol.sapphire_blue;
    end
end
PlotModulation(dirs, el_selectivity_only, cfg)
savePNG(gcf, 600, [figure_dir, 'Memoria_selectivity_brain_only.png'])


el_selectivity_selective = el_selectivity_Memoria(contains(el_selectivity_Memoria.elect_select, 'selective'),:);


load('cdcol_2018.mat')
for i = 1:size(el_selectivity_selective,1)
    if contains(el_selectivity_selective.elect_select{i}, 'math selective') == 1
        cfg.MarkerColor(i,:) = cdcol.light_cadmium_red;
    else
        cfg.MarkerColor(i,:) = cdcol.sapphire_blue;
    end
end
PlotModulation(dirs, el_selectivity_selective, cfg)
savePNG(gcf, 300, [figure_dir, 'Memoria_selectivity_brain_selective.png'])



%% ALL MATH
el_selectivity_all_calc = [el_selectivity_Memoria; el_selectivity_calc_sim];


selectivities = {{'math only'}, {'math selective'}, {'math deact'}, {'no selectivity'}, {'episodic only', 'autobio only'}, {'episodic selective', 'autobio selective'}, {'autobio deact',  'episodic deact'}, {'math and episodic', 'math and autobio'}};
columns = {'Yeo7', 'DK_lobe'};
figure('units', 'normalized', 'outerposition', [0 0 1 1])
for ic = 1:length(columns)
    for i = 1:length(selectivities)
        selectivity = selectivities{i};
        column = columns{ic};
        subplot(2,4,i)
        if strcmp(selectivity, 'math deact') == 1
            el_temp = el_selectivity_all_calc(el_selectivity_all_calc.act_deact_cond1 == -1,:);
        elseif contains(selectivity, {'autobio deact',  'episodic deact'}) == 2
            el_temp = el_selectivity_all_calc(el_selectivity_all_calc.act_deact_cond2 == -1,:);
        else
            el_temp = el_selectivity_all_calc(contains(el_selectivity_all_calc.elect_select, selectivity),:);
        end
        plot_frequency(el_temp, column, 'ascend', 'horizontal')
        title(selectivity)
    end
    savePNG(gcf, 300, [figure_dir, ['Calc_all_selectivity ', column, '.png']])
end
close all


selectivities = {{'math only'},{'episodic only', 'autobio only'}};
columns = {'DK_long_josef'};
figure('units', 'normalized', 'outerposition', [0 0 1 1])
for ic = 1:length(columns)
    for i = 1:length(selectivities)
        selectivity = selectivities{i};
        column = columns{ic};
        subplot(1,2,i)
        if strcmp(selectivity, 'math deact') == 1
            el_temp = el_selectivity_all_calc(el_selectivity_all_calc.act_deact_cond1 == -1,:);
        elseif contains(selectivity, {'autobio deact',  'episodic deact'}) == 2
            el_temp = el_selectivity_all_calc(el_selectivity_all_calc.act_deact_cond2 == -1,:);
        else
            el_temp = el_selectivity_all_calc(contains(el_selectivity_all_calc.elect_select, selectivity),:);
        end
        plot_frequency(el_temp, column, 'ascend', 'horizontal')
        title(selectivity)
    end
    savePNG(gcf, 300, [figure_dir, ['Calc_all_selectivity ', column, '.png']])
end
close all




el_selectivity_only = el_selectivity_all_calc(contains(el_selectivity_all_calc.elect_select, 'only'),:);
sort_tabulate(el_selectivity_only.elect_select, 'descend')

cfg = getPlotCoverageCFG('tasks_group');
cfg.MarkerSize = 10;
cfg.alpha = 0.8;
cfg.views = {'lateral', 'lateral', 'ventral', 'medial', 'medial', 'ventral',};
cfg.hemis = {'left', 'right', 'left', 'left', 'right', 'right', };
cfg.subplots = [2,3];
cfg.figureDim = [0 0 1 1];

load('cdcol_2018.mat')
for i = 1:size(el_selectivity_only,1)
    if contains(el_selectivity_only.elect_select{i}, 'math') == 1
        cfg.MarkerColor(i,:) = cdcol.light_cadmium_red;
    else
        cfg.MarkerColor(i,:) = cdcol.sapphire_blue;
    end
end
PlotModulation(dirs, el_selectivity_only, cfg)
savePNG(gcf, 600, [figure_dir, 'All_calc_selectivity_brain_only.png'])


el_selectivity_selective = el_selectivity_all_calc(contains(el_selectivity_all_calc.elect_select, 'selective'),:);


load('cdcol_2018.mat')
for i = 1:size(el_selectivity_selective,1)
    if contains(el_selectivity_selective.elect_select{i}, 'math selective') == 1
        cfg.MarkerColor(i,:) = cdcol.light_cadmium_red;
    else
        cfg.MarkerColor(i,:) = cdcol.sapphire_blue;
    end
end
PlotModulation(dirs, el_selectivity_selective, cfg)
savePNG(gcf, 300, [figure_dir, 'All_calc_selectivity_brain_selective.png'])


%% Calculia

task = 'Calculia';
sinfo_Calculia = sinfo(strcmp(sinfo.task, task),:);
sinfo_Calculia(contains(sinfo_Calculia.sbj_name, {'S14_74_OD', 'S15_87_RL', 'S16_95_JOB', 'S15_83_RR', 'S16_96_LF'}),:) = [];
el_selectivity_Calculia = concat_elect_select(sinfo_Calculia.sbj_name, task, dirs, vars);
el_selectivity_Calculia = el_selectivity_Calculia(strcmp(el_selectivity_Calculia.WMvsGM, 'GM') | strcmp(el_selectivity_Calculia.WMvsGM, 'WM'), :);
el_selectivity_Calculia = el_selectivity_Calculia(~strcmp(el_selectivity_Calculia.Yeo7, 'FreeSurfer_Defined_Medial_Wall'),:)

selectivities = {{'math only'}, {'math selective'}, {'math deact'}, {'no selectivity'}, {'autobio only'}, {'autobio selective'}, {'autobio deact'}, {'math and autobio'}};
columns = {'Yeo7', 'DK_lobe'};
figure('units', 'normalized', 'outerposition', [0 0 1 1])
for ic = 1:length(columns)
    for i = 1:length(selectivities)
        selectivity = selectivities{i};
        column = columns{ic};
        subplot(2,4,i)
        el_temp = el_selectivity_Memoria(contains(el_selectivity_Memoria.elect_select, selectivity),:);
        plot_frequency(el_temp, column, 'ascend', 'horizontal')
        title(selectivity)
    end
    savePNG(gcf, 300, [figure_dir, ['Memoria_selectivity ', column, '.png']])
end


el_selectivity_only = el_selectivity_Calculia(~contains(el_selectivity_Calculia.elect_select, {'no selectivity', 'digit_active and letter_active'}),:);

el_selectivity_only = el_selectivity_Calculia(contains(el_selectivity_Calculia.elect_select, {'only'}),:);
sort_tabulate(el_selectivity_only.elect_select, 'descend')

cfg = getPlotCoverageCFG('tasks_group');
cfg.MarkerSize = 10;
cfg.alpha = 0.5;
cfg.views = {'lateral', 'lateral', 'ventral', 'medial', 'medial', 'ventral',};
cfg.hemis = {'left', 'right', 'left', 'left', 'right', 'right', };
cfg.subplots = [2,3];
cfg.figureDim = [0 0 1 1];

load('cdcol_2018.mat')
for i = 1:size(el_selectivity_only,1)
    if contains(el_selectivity_only.elect_select{i}, {'digit_active selective', 'digit_active only'}) == 1
        cfg.MarkerColor(i,:) = cdcol.light_cadmium_red;
    elseif contains(el_selectivity_only.elect_select{i}, {'letter_active selective', 'letter_active only'}) == 1
        cfg.MarkerColor(i,:) = cdcol.sapphire_blue;
    else
        cfg.MarkerColor(i,:) = [0 0 0];
    end
end
PlotModulation(dirs, el_selectivity_only, cfg)
savePNG(gcf, 300, [figure_dir, 'Calculia_only_brain_selective.png'])


%% Viz proportions
el_selectivity = simplify_selectivity(el_selectivity_all_calc, 'MMR');
sort_tabulate(el_selectivity.elect_select, 'descend')

el_selectivity_only = el_selectivity(contains(el_selectivity.elect_select, 'only'), :)
el_selectivity_only = el_selectivity_only(~contains(el_selectivity_only.Yeo7, 'Depth'),:)


conditions = {'math', 'memory'};
Yeo7_networks = {'Frontoparietal', 'Dorsal Attention', 'Default', 'Limbic',  'Ventral Attention','Visual', 'Somatomotor'};

frequencies = [];
for i = 1:length(conditions)
    tmp_Yeo7 = el_selectivity_only(contains(el_selectivity_only.elect_select, conditions{i}),:);
    tmp_Yeo7 = sort_tabulate(tmp_Yeo7.Yeo7, 'descend');
    for in = 1:length(Yeo7_networks)
        frequencies(i,in) = tmp_Yeo7{strcmp(tmp_Yeo7.value, Yeo7_networks{in}), 2};
    end
end
frequencies = frequencies'

% frequencies = flip(frequencies');
[frequencies, idx] = sortrows(frequencies, 1, 'descend')
frequencies = flip(frequencies);
Yeo7_networks = Yeo7_networks(idx)
Yeo7_networks = flip(Yeo7_networks)

ba = barh(frequencies, 'stacked' ,'EdgeColor', 'k','LineWidth',2)
ba(1).FaceColor = cdcol.light_cadmium_red;
ba(2).FaceColor = cdcol.sapphire_blue

set(gca,'fontsize',16)
xlabel('Number of electrodes')
yticks(1:length(Yeo7_networks))
ylim([0, length(Yeo7_networks)+1])
yticklabels(Yeo7_networks)
set(gca,'TickLabelInterpreter','none')

for i = 1:length(conditions)
    for in = 1:length(Yeo7_networks)
        if i == 1
            txt = text(frequencies(in,i)-1,in, num2str(frequencies(in,i)), 'FontSize', 20, 'HorizontalAlignment', 'right', 'Color', 'w');
        elseif i > 1
            txt = text(frequencies(in,i)+sum(frequencies(in,1:i-1))-1,in, num2str(frequencies(in,i)), 'FontSize', 20, 'HorizontalAlignment', 'right', 'Color', 'w');
        end
    end
end
title('Frequency of math vs. memory only sites per intrinsic network')
savePNG(gcf, 300, [figure_dir, 'math_all_frequencies_Yeo7_stacked.png'])






conditions = {'math', 'memory'};
hemis = {'L', 'R'};

frequencies = [];
for i = 1:length(conditions)
    LvsR = el_selectivity_only(contains(el_selectivity_only.elect_select, conditions{i}),:);
    LvsR = sort_tabulate(LvsR.LvsR, 'descend');
    for in = 1:length(hemis)
        frequencies(i,in) = LvsR{strcmp(LvsR.value, hemis{in}), 2};
    end
end
frequencies = frequencies'

% frequencies = flip(frequencies');
[frequencies, idx] = sortrows(frequencies, 1, 'descend')
frequencies = frequencies;
hemis = hemis(idx)
hemis = flip(hemis)

ba = bar(frequencies, 'stacked' ,'EdgeColor', 'k','LineWidth',2)
ba(1).FaceColor = cdcol.light_cadmium_red;
ba(2).FaceColor = cdcol.sapphire_blue

set(gca,'fontsize',16)
ylabel('Number of electrodes')
xlim([0, length(hemis)+1])
xlabel('Hemispheres')
xticklabels({'Left', 'Right'})
set(gca,'TickLabelInterpreter','none')

for i = 1:length(conditions)
    for in = 1:length(Yeo7_networks)
        if i == 1
            txt = text(in, frequencies(in,i)-1, num2str(frequencies(in,i)), 'FontSize', 20, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', 'w');
        elseif i > 1
            txt = text(in, frequencies(in,i)+sum(frequencies(in,1:i-1))-1, num2str(frequencies(in,i)), 'FontSize', 20, 'HorizontalAlignment', 'center',  'VerticalAlignment', 'top', 'Color', 'w');
        end
    end
end
title('Frequency of math vs. memory only sites per hemi network')
savePNG(gcf, 300, [figure_dir, 'math_all_frequencies_hemi_stacked.png'])


%%
[DOCID,GID] = getGoogleSheetInfo('math_network','cohort');
sinfo = GetGoogleSpreadsheet(DOCID, GID);
subject_names = sinfo.sbj_name;
sinfo = sinfo(strcmp(sinfo.subjVar, '1'),:);
% Filer for usable dat
subjects = unique(sinfo.sbj_name);


for i = 1:length(subjects)
    try
        CreateSubjVar(subjects{i}, comp_root, server_root, code_root)
    catch
        fname = sprintf('%s/%s_subjVar_error.csv',dirs.comp_root, subjects{i});
        csvwrite(fname, 's')
    end
end


parfor i = 1:length(subjects)
    CreateSubjVar(subjects{i}, comp_root, server_root, code_root)
end



CreateSubjVar('S16_100_AF', comp_root, server_root, code_root)



%% VTCLoc
%% Univariate Selectivity
vars = {'chan_num', 'FS_label', 'LvsR','MNI_coord', 'WMvsGM', 'sEEG_ECoG', 'DK_lobe', 'Yeo7', 'Yeo17', 'DK_long_josef', ...
    'elect_select', 'act_deact_cond1', 'act_deact_cond2', 'sc1c2_FDR', 'sc1b1_FDR' , 'sc2b2_FDR', ...
    'sc1c2_Pperm', 'sc1b1_Pperm', 'sc2b2_Pperm', 'sc1c2_tstat', 'sc1b1_tstat', 'sc2b2_tstat'};


task = 'VTCLoc';
sinfo_VTCLoc = sinfo(strcmp(sinfo.task, task),:);
sinfo_VTCLoc(strcmp(sinfo_VTCLoc.sbj_name,'S17_118_TW'),:) = [];
sinfo_VTCLoc(strcmp(sinfo_VTCLoc.sbj_name,'S15_91_RP'),:) = [];




el_selectivity_VTCloc_faces = concat_elect_select(sinfo_VTCLoc.sbj_name, task, dirs, vars);
el_selectivity_VTCloc_numbers = concat_elect_select(sinfo_VTCLoc.sbj_name, task, dirs, vars);
el_selectivity_VTCloc_words = concat_elect_select(sinfo_VTCLoc.sbj_name, task, dirs, vars);

el_selectivity_VTC = el_selectivity_VTCloc_faces
el_selectivity_VTC.elect_select_faces = el_selectivity_VTC.elect_select
el_selectivity_VTC.elect_select_words = el_selectivity_VTCloc_words.elect_select
el_selectivity_VTC.elect_select_numbers = el_selectivity_VTCloc_numbers.elect_select


VTC_selective = el_selectivity_VTC(contains(el_selectivity_VTC.elect_select_faces, {'faces only', 'faces selective'}) | contains(el_selectivity_VTC.elect_select_numbers, {'numbers only', 'numbers selective'}) | contains(el_selectivity_VTC.elect_select_words, {'words only', 'words selective'}),:);
sort_tabulate(VTC_selective.elect_select, 'descend')

cfg = getPlotCoverageCFG('tasks_group');
cfg.MarkerSize = 10;
cfg.alpha = 0.4;
cfg.views = {'lateral', 'lateral', 'ventral', 'medial', 'medial', 'ventral',};
cfg.hemis = {'left', 'right', 'left', 'left', 'right', 'right', };
cfg.subplots = [2,3];
cfg.figureDim = [0 0 1 1];
cfg.CorrectFactor = 10;

load('cdcol_2018.mat')
for i = 1:size(VTC_selective,1)
    if contains(VTC_selective.elect_select_faces{i}, {'faces only', 'faces selective'}) == 1
        cfg.MarkerColor(i,:) = cdcol.grass_green;
    elseif contains(VTC_selective.elect_select_words{i}, {'words only', 'words selective'}) == 1
        cfg.MarkerColor(i,:) = cdcol.orange;
    elseif contains(VTC_selective.elect_select_numbers{i}, {'numbers only', 'numbers selective'}) == 1
        cfg.MarkerColor(i,:) = cdcol.light_cadmium_red;
    end
    
end
PlotModulation(dirs, VTC_selective, cfg)
savePNG(gcf, 300, [figure_dir, 'Calculia_only_brain_selective.png'])




%% ReadNumWord
%% Univariate Selectivity


vars = {'chan_num', 'FS_label', 'LvsR','MNI_coord', 'WMvsGM', 'sEEG_ECoG', 'DK_lobe', 'Yeo7', 'Yeo17', 'DK_long_josef', ...
    'elect_select', 'act_deact_cond1', 'act_deact_cond2', 'sc1c2_FDR', 'sc1b1_FDR' , 'sc2b2_FDR', ...
    'sc1c2_Pperm', 'sc1b1_Pperm', 'sc2b2_Pperm', 'sc1c2_tstat', 'sc1b1_tstat', 'sc2b2_tstat'};


task = 'ReadNumWord';
sinfo_ReadNumWord_numbers = sinfo(strcmp(sinfo.task, task),:);
sinfo_ReadNumWord_numbers(strcmp(sinfo_ReadNumWord_numbers.sbj_name,'S12_36_SrS'),:) = [];





el_selectivity_ReadNumWord_numbers = concat_elect_select(sinfo_VTCLoc.sbj_name, task, dirs, vars);


ReadNumWord_selective = el_selectivity_ReadNumWord_numbers(contains(el_selectivity_ReadNumWord_numbers.elect_select, {'numbers selective', 'numbers only'}),:);
sort_tabulate(ReadNumWord_selective.elect_select, 'descend')

cfg = getPlotCoverageCFG('tasks_group');
cfg.MarkerSize = 10;
cfg.alpha = 0.4;
cfg.views = {'lateral', 'lateral', 'ventral', 'medial', 'medial', 'ventral',};
cfg.hemis = {'left', 'right', 'left', 'left', 'right', 'right', };
cfg.subplots = [2,3];
cfg.figureDim = [0 0 1 1];
cfg.CorrectFactor = 10;

load('cdcol_2018.mat')
for i = 1:size(VTC_selective,1)
    if contains(VTC_selective.elect_select_faces{i}, {'faces only', 'faces selective'}) == 1
        cfg.MarkerColor(i,:) = cdcol.grass_green;
    elseif contains(VTC_selective.elect_select_words{i}, {'words only', 'words selective'}) == 1
        cfg.MarkerColor(i,:) = cdcol.orange;
    elseif contains(VTC_selective.elect_select_numbers{i}, {'numbers only', 'numbers selective'}) == 1
        cfg.MarkerColor(i,:) = cdcol.light_cadmium_red;
    end
    
end
PlotModulation(dirs, VTC_selective, cfg)
savePNG(gcf, 300, [figure_dir, 'Calculia_only_brain_selective.png'])



%%
vars = {'sbj_name', 'LvsR','MNI_coord', 'WMvsGM', 'sEEG_ECoG', 'DK_lobe', 'Yeo7', 'Yeo17'};
subjVars = ConcatSubjVars(subjects, dirs, vars);

sort_tabulate(subjVars.Yeo7, 'descend')





for i = 1:size(subjVars,1)
    str_tmp = strsplit(subjVars.sbj_name{i}, '_');
    
    subjVars.sbj_number(i) = str2num(str_tmp{2});
end

for i = 1:size(subjVars_old,1)
    subjVars_old.sbj_number(i) = str2num(subjVars_old.sbj_name{i});
end


subjVars = subjVars(ismember(subjVars.sbj_number, subjVars_old.sbj_number),:);




subjVars(contains(subjVars.WMvsGM, {'empty', 'FreeSurfer_Defined_Medial_Wall', 'EMPTY'} ),:) = []
subjVars_old(contains(subjVars_old.WMvsGM, {'empty', 'FreeSurfer_Defined_Medial_Wall', 'EMPTY'} ),:) = []



sort_tabulate(subjVars.WMvsGM, 'descend')
sort_tabulate(subjVars_old.WMvsGM, 'descend')


%% Group analyses

vars = {'chan_num', 'FS_label', 'LvsR','MNI_coord', 'WMvsGM', 'sEEG_ECoG', 'DK_lobe', 'Yeo7', 'Yeo17', 'DK_long_josef', ...
    'elect_select', 'act_deact_cond1', 'act_deact_cond2', 'sc1c2_FDR', 'sc1b1_FDR' , 'sc2b2_FDR', ...
    'sc1c2_Pperm', 'sc1b1_Pperm', 'sc2b2_Pperm', 'sc1c2_tstat', 'sc1b1_tstat', 'sc2b2_tstat'};


regions = {'ITG','IPS', 'SPL', 'MFG'};
task = 'MMR';
sinfo_task = sinfo(strcmp(sinfo.task, task),:);
el_selectivity = concat_elect_select(sinfo_task.sbj_name, task, dirs, vars);
el_selectivity = el_selectivity(strcmp(el_selectivity.WMvsGM, 'GM') | strcmp(el_selectivity.WMvsGM, 'WM'), :);
el_selectivity = el_selectivity(~strcmp(el_selectivity.Yeo7, 'FreeSurfer_Defined_Medial_Wall'),:);
el_selectivity = el_selectivity(contains(el_selectivity.elect_select, 'math only'),:);


for i = 1:length(regions)
    el_selec_tmp = el_selectivity(strcmp(el_selectivity.DK_long_josef, regions{i}),[1,end-1]);
    if isempty(el_selec_tmp)
    else
        data_all = concatenate_multiple_elect(el_selec_tmp, task, dirs, 'Band', 'HFB', 'stim');
        cond_names = {'autobio', 'math'};
        column = 'condNames';
        subplot(length(regions),1, i)
        plot_group_elect(data_all,task, cond_names, column);
        if strcmp(regions{i} , 'ITG')
            ylim([-0.4 2])
        else
            ylim([-0.2 1])
        end
        title([regions{i} ': ' num2str(length(data_all.wave)), ' electrodes'])
    end
end
savePNG(gcf, 600, [figure_dir, task, '_group_regions_selective.png'])


data_all = concatenate_multiple_elect(el_selectivity, task, dirs, 'Band', 'HFB', 'stim');


el_selec_tmp = el_selectivity(strcmp(el_selectivity.DK_long_josef, regions{2}),[1,end-1]);
data_all = concatenate_multiple_elect(el_selec_tmp, task, dirs, 'Band', 'HFB', 'stim');
cond_names = {'autobio', 'math'};
column = 'condNames';
stats_params = genStatsParams(task);
STATS = stats_group_elect(data_all,data_all.time, task,cond_names, column, stats_params);


title([regions{2} ': ' num2str(length(data_all.wave)), ' electrodes'])


% Group analyses

1. Select which electrodes to include
-Statistical
Compare all trails agains baseline
permutation test between avg baseline period whithin trial vs. avg 1s period within trial
-Anatomical
ROIS
After these steps you should have the electrode_list

2. Concatenate all electrodes of interest
data_all = concatenate_multiple_elect(electrode_list, task, dirs, 'Band', 'HFB', 'stim');

3. Compare two conditions across electrodes
-Define conditions and time window
cond_names = {'autobio', 'math'};
column = 'condNames';
stats_params = genStatsParams(task);
-Average each electrode per condition within the 1s period
(now you have a single value per electrode per condition)
Ready to compare electrodes with independent sample t-test
STATS = stats_group_elect(data_all,data_all.time, task,cond_names, column, stats_params);




%% Vizualize proportions


vars = {'chan_num', 'FS_label', 'LvsR','MNI_coord', 'WMvsGM', 'sEEG_ECoG', 'DK_lobe', 'Yeo7', 'Yeo17', 'DK_long_josef', 'elect_select'};
elec_select = concat_selectivity_tasks(sinfo, 'calc_simultaneous', vars, dirs);



cond_names = {'math only', 'memory only'};
brain_group_list = {'Frontoparietal', 'Dorsal Attention', 'Default', 'Limbic',  'Ventral Attention','Visual', 'Somatomotor'};
brain_group_list = {'Depth', 'Frontoparietal', 'Dorsal Attention', 'Default', 'Limbic',  'Ventral Attention','Visual', 'Somatomotor'};

brain_group = 'Yeo7';


plot_proportion_selectivity(elec_select, 'MMR', cond_names, brain_group, brain_group_list)

brain_group_list = {'L','R'};
brain_group = 'LvsR';
plot_proportion_selectivity(el_selectivity_calc_sim, 'MMR', cond_names, brain_group, brain_group_list)



%% MMR and memoria comparison electrode by electrode
sinfo_MMR = sinfo(strcmp(sinfo.task, 'MMR'),:);
sinfo_Memoria = sinfo(strcmp(sinfo.task, 'Memoria'),:);
subjects_MMR_Memoria = intersect(sinfo_MMR.sbj_name, sinfo_Memoria.sbj_name);

el_selectivity_MMR = concat_elect_select(subjects_MMR_Memoria, 'MMR', dirs, vars);
el_selectivity_Memoria = concat_elect_select(subjects_MMR_Memoria, 'Memoria', dirs, vars);
elect_select_MMR = el_selectivity_MMR.elect_select;
elect_select_Memoria = el_selectivity_Memoria.elect_select;

elec_select = el_selectivity_MMR(:,[1:12,end-1])
elec_select.elect_select_MMR = elect_select_MMR;
elec_select.elect_select_Memoria = elect_select_Memoria;
elec_select = elec_select(strcmp(elec_select.WMvsGM, 'GM') | strcmp(elec_select.WMvsGM, 'WM'), :);
elec_select = elec_select(~strcmp(elec_select.Yeo7, 'FreeSurfer_Defined_Medial_Wall'),:);

sort_tabulate(elec_select.elect_select_Memoria(strcmp(elec_select.elect_select_MMR, 'math only')), 'descend')
sort_tabulate(elec_select.elect_select_MMR(strcmp(elec_select.elect_select_Memoria, 'math only')), 'descend')



corrcoef([el_selectivity_MMR.sc1c2_tstat, el_selectivity_Memoria.sc1c2_tstat], 'rows','complete')


scatter(plotvals(:,1),plotvals(:,2),40,c,'filled'),colorbar;



scatter_kde(el_selectivity_MMR.sc1c2_tstat, el_selectivity_Memoria.sc1c2_tstat,  'filled', 'MarkerSize', 50)
colormap viridis

sum_t = nansum([el_selectivity_MMR.sc1c2_tstat, el_selectivity_Memoria.sc1c2_tstat], 2);
rgb = vals2colormap(sum_t*-1, 'cmRedBlue', [-10 10]);

scatter(el_selectivity_MMR.sc1c2_tstat, el_selectivity_Memoria.sc1c2_tstat, 100, rgb, 'filled')
xlabel('T-value Calc simultaneous')
ylabel('T-value Calc sequential')
set(gca,'fontsize',16)
box on
axis square

savePNG(gcf, 600, [figure_dir, task, 'MMR_Memoria_correspondence.png'])








%% ROL
%% ROL from Jessica and Omri NC

project_name = 'UCLA';

parfor i = 1:size(sinfo_UCLA,1)
    getROLALL_NC(sinfo_UCLA.sbj_name{i},project_name,[],dirs,[],'HFB',[],'condNames',{'autobio', 'math'}) ;% condNames
end

plot_ROL_scatter


ROL_var = {'onsets', 'peaks'};
elecs = [105, 16, 61];
col = cbrewer2('Blues',6)
col = col(3:end,:)
col = cool(3)
cond_names = {'math'};
plot_ROL_scatter(ROL, ROL_var, elecs, cond_names, 'each column separately', 2,col)
plot_ROL_scatter(ROL, ROL_var, elecs, cond_names, 'by one column', 3,col)



getROLALL_NC('S20_151_HT',project_name,[],dirs,[1:156],'HFB',[],'condNames',{'autobio', 'math'}) ;% condNames

getROLALL_NC('S13_57_TVD',project_name,[],dirs,[105,16,61],'HFB',[],'condNames',{'math'}) ;% condNames


%% Integrate selectivity and ROL
vars = {'chan_num', 'FS_label', 'LvsR','MNI_coord', 'fsaverageINF_coord', 'WMvsGM', 'sEEG_ECoG', 'DK_lobe', 'Yeo7', 'Yeo17', 'DK_long_josef', ...
    'elect_select', 'act_deact_cond1', 'act_deact_cond2', 'sc1c2_FDR', 'sc1b1_FDR' , 'sc2b2_FDR', ...
    'sc1c2_Pperm', 'sc1b1_Pperm', 'sc2b2_Pperm', 'sc1c2_tstat', 'sc1b1_tstat', 'sc2b2_tstat'};

task = 'MMR';
sinfo_MMR = sinfo(strcmp(sinfo.task, task),:);
el_selectivity_MMR = concat_elect_select(sinfo_MMR.sbj_name, task, dirs, vars);


el_selectivity_MMR_ROL = concat_elect_select_rol(sinfo_MMR.sbj_name, task, dirs, vars);


el_selectivity = el_selectivity_MMR_ROL(strcmp(el_selectivity_MMR_ROL.elect_select, 'math only'), :)
el_selectivity = el_selectivity(~strcmp(el_selectivity.Yeo7, 'FreeSurfer_Defined_Medial_Wall'),:);
el_selectivity = el_selectivity(~contains(el_selectivity.DK_long_josef, {'OUT OF BRAIN', 'WHITE MATTER', 'EXCLUDE', 'POSTCENTRAL GYRUS', 'PRECENTRAL GYRUS'}),:);

el_selectivity.ROL_math_avg = cellfun(@nanmean, el_selectivity.ROL_math_onsets)
el_selectivity = el_selectivity(~isnan(el_selectivity.ROL_math_avg),:);


dir_save = '/Volumes/LBCN8T/Stanford/data/electrode_localization/plots/';

col_group = 'task_group'; % or task_group
cols = hsv(length(labels_josef));
cfg = getPlotCoverageCFG('tasks_group');
cfg.figureDim = [0 0 1 0.7];
cfg.views = {'lateral', 'lateral', 'ventral', 'ventral', 'medial', 'medial', 'posterior', 'posterior'}; %{'lateral', 'lateral', 'ventral', 'ventral'};
cfg.hemis = {'left', 'right', 'left', 'right', 'left', 'right', 'left', 'right'}; %{'left', 'right', 'left', 'right'};
cfg.subplots = [2,4];  % 2,2
cfg.plot_label = 0;
cfg.colum_label = {'sbj_name', 'FS_label'};
cfg.alpha = 0.1;
cfg.MarkerSize = 10;
cfg.MarkerColor = cols(i,:);

PlotModulation(dirs, el_selectivity, cfg)
fname = sprintf('%selectrodes_%s.png', dir_save, labels_josef{i});
savePNG(gcf, 300, fname)
close all


el_selectivity.LvsR = repmat({'L'}, size(el_selectivity,1),1,1)
el_selectivity.MNI_coord(:,1) = abs(el_selectivity.MNI_coord(:,1))*-1

cfg.ind = cellfun(@nanmean, el_selectivity.ROL_math_onsets)
cfg.MarkerColor = [];
cfg.views = {'lateral', 'ventral'};
cfg.hemis = {'left', 'left'};
cfg.subplots = [1,2];  % 2,2
cfg.MarkerSize = 20;
cfg.Colormap = 'Reds'
cfg.Cortex = 'MNI';
PlotModulation(dirs, el_selectivity, cfg)



sort_tabulate(el_selectivity.DK_long_josef, 'descend')
labels_plot = {'MFG', 'IPS', 'ITG', 'SPL', 'FG', 'IFG'}
el_selectivity = el_selectivity(contains(el_selectivity.DK_long_josef, labels_plot) & ~strcmp(el_selectivity.DK_long_josef, {'mSFG'}),:);
mean_labels = varfun(@median,el_selectivity,'InputVariables','ROL_math_avg', 'GroupingVariables','DK_long_josef');
[mean_labels, idx] = sortrows(mean_labels, 3);
% rol per region
el_selectivity.DK_long_josef_sort = el_selectivity.DK_long_josef;

for i = 1:length(mean_labels)
    labels_tmp = el_selectivity.DK_long_josef_sort(strcmp(el_selectivity.DK_long_josef_sort, mean_labels.DK_long_josef{i})) =
    
end

ROL_vars = {'ROL_math_onsets', 'ROL_math_peaks'};
ROL_vars = {'sc1c2_tstat', 'sc1b1_tstat'};



for i = 1:length(ROL_vars)
    subplot(1,2,i)
    %     el_selectivity.ROL_tmp = cellfun(@nanmean, el_selectivity.(ROL_vars{i}))
    %     ROL_tmp = varfun(@median,el_selectivity,'InputVariables','ROL_tmp', 'GroupingVariables','DK_long_josef');
    ROL_tmp = varfun(@median,el_selectivity,'InputVariables',ROL_vars{i}, 'GroupingVariables','DK_long_josef');
    
    [~, idx] = sortrows(ROL_tmp, 3);
    vi = violinplot(el_selectivity.(ROL_vars{i}), el_selectivity.DK_long_josef, idx); % accuracy Min Max Result Abs deviant decade cross order;
    for iiii = 1:length(vi)
        cols_vi = cbrewer2('Greens', length(vi)+2);
        cols_vi(1:2,:) = [];
        %         cols_vi = flip(cols_vi);
        vi(iiii).ViolinPlot.FaceColor = 'none';
        vi(iiii).ViolinAlpha =  1;
        vi(iiii).ScatterPlot.MarkerFaceColor = cols_vi(iiii,:);
        vi(iiii).EdgeColor = cols_vi(iiii,:);
        vi(iiii).BoxColor = cols_vi(iiii,:);
        vi(iiii).ViolinPlot.LineWidth = 2;
    end
    set(gca,'FontSize', 16)
    xtickangle(45)
    if contains(ROL_vars{i}, 'onsets')
        ylabel('ROL (ms)');
    elseif contains(ROL_vars{i}, 'peaks')
        
        ylabel('Time to peak (ms)');
    end
    if contains(ROL_vars{i}, 'sc1c2_tstat')
        ylabel('T-stat math vs. memory');
    else
        ylabel('T-stat math vs. baseline');
    end
    axis square
end


suptitle('Math only electrodes timing', 'FontSize', 20)

%% Single subject simultaneous
task = 'MMR';
sinfo_MMR = sinfo(strcmp(sinfo.task, task),:);

vars = {'chan_num', 'FS_label', 'LvsR','MNI_coord', 'fsaverageINF_coord', 'WMvsGM', 'sEEG_ECoG', 'DK_lobe', 'Yeo7', 'Yeo17', 'DK_long_josef', ...
    'elect_select', 'act_deact_cond1', 'act_deact_cond2', 'sc1c2_FDR', 'sc1b1_FDR' , 'sc2b2_FDR', ...
    'sc1c2_Pperm', 'sc1b1_Pperm', 'sc2b2_Pperm', 'sc1c2_tstat', 'sc1b1_tstat', 'sc2b2_tstat'};

el_selectivity_MMR_ROL = concat_elect_select_rol(sinfo_MMR.sbj_name, 'MMR', dirs, vars);


labels_plot = {'FG', 'ITG', 'SPL', 'IPS', 'IFG', 'MFG', 'SFG'}
lobes_plot = {'Temporal', 'Parietal', 'Frontal'}
el_selectivity = el_selectivity_MMR_ROL(strcmp(el_selectivity_MMR_ROL.elect_select, {'math only'}), :)
el_selectivity = el_selectivity(contains(el_selectivity.DK_long_josef, labels_plot) & ~strcmp(el_selectivity.DK_long_josef, {'mSFG'}),:);

el_selectivity = el_selectivity(~strcmp(el_selectivity.Yeo7, 'FreeSurfer_Defined_Medial_Wall'),:);
el_selectivity = el_selectivity(~contains(el_selectivity.DK_long_josef, {'OUT OF BRAIN', 'WHITE MATTER', 'EXCLUDE', 'POSTCENTRAL GYRUS', 'PRECENTRAL GYRUS'}),:);
el_selectivity = el_selectivity(contains(el_selectivity.DK_long_josef, labels_plot) & ~strcmp(el_selectivity.DK_long_josef, {'mSFG'}),:);
el_selectivity.ROL_math_avg = cellfun(@nanmean, el_selectivity.ROL_math_onsets)


el_selectivity.labels_num = [];
el_selectivity.lobes_num = [];
for i = 1:size(el_selectivity,1)
    el_selectivity.labels_num(i) = find(strcmp(el_selectivity.DK_long_josef{i}, labels_plot)); 
end

el_tmp_all = [];
for i = 1:length(subjects)
    el_temp = el_selectivity(strcmp(el_selectivity.sbj_name, subjects{i}),:);
    if size(el_temp,1) > 1 & std(el_temp.labels_num) > 0
        el_tmp_all = [el_tmp_all; el_temp];
    else
    end
end

subjects = unique(el_tmp_all.sbj_name);
cols = hsv(length(subjects));


for i = 1:length(subjects)
    subplot(6,5,i)
    el_temp = el_tmp_all(strcmp(el_tmp_all.sbj_name, subjects{i}),:);
    plot(el_temp.labels_num, el_temp.ROL_math_avg, 'o', 'MarkerFaceColor', cols(i,:), 'MarkerEdgeColor', cols(i,:), 'Color', cols(i,:), 'MarkerSize', 10)
    boxplot([el_temp.labels_num, el_temp.ROL_math_avg],'PlotStyle','compact')
    
    hold on
    xlim([0.5 7.5])
    ylim([0 .4])
    set(gca,'xticklabels', labels_plot)
    set(gca,'FontSize', 12)
    ylabel('ROL (ms)')
end


cols = viridis(length(subjects));


sub_oder = []
for i = 1:length(subjects)
    el_temp = el_selectivity(strcmp(el_selectivity.sbj_name, subjects{i}),:);
    
    means = varfun(@nanmean,el_temp,'InputVariables','ROL_math_avg', 'GroupingVariables','labels_num');
    stds =  varfun(@nanstd,el_temp,'InputVariables','ROL_math_avg', 'GroupingVariables','labels_num');
    sub_oder(i,1) = round(sum(means.labels_num)+std(means.labels_num));
    
end
[~, idx] = sort(sub_oder, 'ascend');
subjects = subjects(idx);
subjects(16) = []


for i = 1:length(subjects)
    subplot(2,13,i)
    el_temp = el_selectivity(strcmp(el_selectivity.sbj_name, subjects{i}),:);
    
    %     plot(el_tmp.labels_num, el_tmp.ROL_math_avg, '-o', 'MarkerFaceColor', cols(i,:), 'MarkerEdgeColor', cols(i,:), 'Color', cols(i,:), 'MarkerSize', 10)
    hold on
    
    means = varfun(@nanmean,el_temp,'InputVariables','ROL_math_avg', 'GroupingVariables','labels_num');
    stds =  varfun(@nanstd,el_temp,'InputVariables','ROL_math_avg', 'GroupingVariables','labels_num');
    plot(means.labels_num, means.nanmean_ROL_math_avg, '-o', 'MarkerFaceColor', cols(i,:), 'MarkerEdgeColor', cols(i,:), 'Color', cols(i,:), 'MarkerSize', 7, 'LineWIdth', 2)
    for ii = 1:size(means,1)
        mt = means.nanmean_ROL_math_avg(ii);
        st = stds.nanstd_ROL_math_avg(ii);
        %         plot(means.labels_num(ii), mt, '-o', 'MarkerFaceColor', cols(i,:), 'MarkerEdgeColor', cols(i,:), 'Color', cols(i,:), 'MarkerSize', 10)
        if st > 0
            line([means.labels_num(ii) means.labels_num(ii)], [mt-st mt+st], 'Color', cols(i,:), 'LineWidth', 2)
        else
        end
    end
    xlim([0.5 7.5])
    ylim([0 0.45])
    
    set(gca,'xticklabels', labels_plot)
    
    set(gca,'FontSize', 8)
    if i == 1 || i == 14
        ylabel('ROL (ms)')
    else
        set(gca,'ytick',[])
        xtickangle(45)
        set(gca,'xtick',1:7)
    end
    
    if i < 14
        set(gca,'xtick',[])
        
    else
    end
    %     set(gca,'color', [.9 .9 .9]);
    %     set(gcf,'color', [.9 .9 .9]);
    set(gcf,'color', 'w')
    %     grid on
end


cfg.chan_highlight = 1;
cfg.highlight_col = [1 0 0];

cfg.views = {'ventral', 'ventral', 'lateral', 'lateral'};
cfg.hemis = {'right', 'left', 'left', 'right'};
cfg.figureDim = [0 0 .5 1];
cfg.subplots = [2, 2];
cfg.alpha = 0.6;
cfg.MarkerSize = 15;
cfg.MarkerSize_chan_highlight = 10;

for i = 1:length(subjects)
    s = subjects{i};
    el_temp = el_selectivity(strcmp(el_selectivity.sbj_name, s),:);
    load([dirs.original_data filesep  s filesep 'subjVar_'  s '.mat']);
    for ii = 1:size(el_temp,1)
        cfg.chan_highlight = el_temp.chan_num(ii);
        PlotCoverageElect(subjVar, cfg)
        fname = sprintf('%scoverage/%s_%s_%s_coverage_math_%s.png',dirs.result_dir, el_temp.DK_long_josef{ii}, num2str(el_temp.chan_num(ii)), s, task);
        savePNG(gcf, 300, fname)
    end
end



%% Temporal dynamics
el_selectivity_MMR_ROL.ROL_math_onsets




%% Cross correlation
labels_plot = {'FG', 'ITG', 'SPL', 'IPS', 'IFG', 'MFG', 'SFG'}

el_selectivity_MMR_ROL.DK_long_josef(strcmp(el_selectivity_MMR_ROL.DK_long_josef, 'mSFG')) = {'MFG'};

el_select_math = el_selectivity_MMR_ROL(contains(el_selectivity_MMR_ROL.DK_long_josef, labels_plot),:)
el_select_math = el_select_math(contains(el_select_math.elect_select, {'math only', 'math selective'}),:);
subjects = unique(el_select_math.sbj_name)


regions = {'ITG', 'IPS'};
el_pair = el_select_math(contains(el_select_math.DK_long_josef, regions),:);



region_pres = [];
s_sim = [];
for i = 1:length(subjects)
    s_tmp = el_pair(strcmp(el_pair.sbj_name, subjects{i}),:);
    for ii = 1:length(regions)
        region_pres{i}(ii) = sum(strcmp(s_tmp.DK_long_josef, regions{ii}));
    end
    
    if sum(region_pres{i}>0) == length(regions)
        s_sim(i) = 1;
    else
        s_sim(i) = 0;
    end
    
end
subjects = subjects(s_sim == 1);
subjects


dirs.figures_dir = '/Users/pinheirochagas/Pedro/Stanford/papers/spatiotempoal_dynamics_math/figures/';

% pairs
% single
for i = 1:length(subjects)
    s_tmp = el_select_math(strcmp(el_select_math.sbj_name, subjects{i}),:);
    electrodes = nchoosek(s_tmp.chan_num, 2);
    if size(electrodes,1) > 1
        
        for ii = 1:size(electrodes)
            
            electrodes_labels = [s_tmp.DK_long_josef(s_tmp.chan_num == electrodes(ii,1)) s_tmp.DK_long_josef(s_tmp.chan_num == electrodes(ii,2))];
            [~, elec_order] = sort(electrodes_labels);
            electrodes_labels = electrodes_labels(elec_order);
            elect_plot = electrodes(ii,elec_order);
            
            corr_gen_time(subjects{i}, 'MMR', dirs, elect_plot, electrodes_labels, [-0.1 3], 25, 'condNames', {'math'});
            fdir = sprintf('%s%s/', dirs.figures_dir, [electrodes_labels{1} '_' electrodes_labels{2}]);
            if ~exist(fdir, 'dir')
                mkdir(fdir)
            else
            end
            fname_gen_corr = sprintf('%s%s/%s_%s_%s_%s_%s_corr_gen.png', dirs.figures_dir, [electrodes_labels{1} '_' electrodes_labels{2}], subjects{i}, num2str(elect_plot(1)), electrodes_labels{1}, num2str(elect_plot(2)), electrodes_labels{2});
            savePNG(gcf, 300, fname_gen_corr)
            
            cfg.chan_highlight = electrodes(ii,:);
            cfg.highlight_col = repmat([1 0 0], length(cfg.chan_highlight),1);
            cfg.views = {'ventral', 'lateral'};
            if strcmp(subjVar.elinfo.LvsR(cfg.chan_highlight(1)), 'R')
                cfg.hemis = {'right', 'right'};
            else
                cfg.hemis = {'left', 'left'};
            end
            cfg.figureDim = [0 0 .7 .5];
            cfg.subplots = [1,2];
            cfg.alpha = 0.6;
            cfg.MarkerSize = 15;
            cfg.MarkerSize_chan_highlight = 10;
            cfg.correction_factor = 10;
            
            load([dirs.original_data filesep  subjects{i} filesep 'subjVar_'  subjects{i} '.mat']);
            PlotCoverageElect(subjVar, cfg)
            fname_coverage = sprintf('%s%s/%s_%s_%s_%s_%s_coverage_gen.png', dirs.figures_dir, [electrodes_labels{1} '_' electrodes_labels{2}], subjects{i},  num2str(elect_plot(1)), electrodes_labels{1},  num2str(elect_plot()), electrodes_labels{2});
            savePNG(gcf, 300, fname_coverage)
            
            
            close all
        end
    else
    end
end




for i = 1:length(subjects)
    s_tmp = el_select_math(strcmp(el_select_math.sbj_name, subjects{i}),:);
    electrodes = s_tmp.chan_num;
    if size(electrodes,1) > 1
        
        for ii = 1:size(electrodes)
            electrodes_labels = [s_tmp.DK_long_josef(s_tmp.chan_num == electrodes(ii)) s_tmp.DK_long_josef(s_tmp.chan_num == electrodes(ii))];
            corr_gen_time(subjects{i}, 'MMR', dirs, [electrodes(ii) electrodes(ii)], electrodes_labels, [-0.1 3], 25, 'condNames', {'math'});
            
            
            
            
            
            
            fdir = sprintf('%s%s/', dirs.figures_dir, [electrodes_labels{1} '_' electrodes_labels{2}]);
            if ~exist(fdir, 'dir')
                mkdir(fdir)
            else
            end
            fname = sprintf('%s%s/%s_%s_%s_%s_%s_corr_gen.png', dirs.figures_dir, [electrodes_labels{1} '_' electrodes_labels{2}], subjects{i}, num2str(electrodes(ii)), electrodes_labels{1}, num2str(electrodes(ii)), electrodes_labels{2});
            savePNG(gcf, 300, fname)
            close all
        end
        
    else
    end
    
end

% organize


%% state space
corr_els = corr(data_el{1}, data_el{2});

[coeff,score,latent,tsquared] = pca(corr_els,'NumComponents',3);
plot3(score(:,1),score(:,2),score(:,3), 'Color', 'b', 'LineWidth', 3)
axis square
xlabel('1st PC')
ylabel('2nd PC')
zlabel('3rd PC')
grid on
set(gcf,'color', 'w')
set(gca,'FontSize', 20)




%% Calculate x_corr multiple subjects

for i = 1:length(subjects)
    try
        xcorr_multiple_subj(subjects{i}, el_select_math, task, dirs)
        subj_error(i) = 0;
    catch
        subj_error(i) = 1;
    end
end


%% Plot x_corr multiple subjects
for i = 1:length(subjects)
    s = subjects{i};
    fname = sprintf('%s%s_%s_x_corr.mat', dirs.result_dir, s, task);
    load(fname,'xcorr_all')
    s_tmp = el_select_math(strcmp(el_select_math.sbj_name, s),:);
    block_names = BlockBySubj(s, task);
    electrodes = nchoosek(s_tmp.chan_num, 2);
    col = viridis(size(electrodes,1));
    for ii = 1:size(electrodes,1)
        sp = numSubplots(size(electrodes,1));
        subplot(sp(1),sp(2),ii)
        e1 = electrodes(ii,1);
        e2 = electrodes(ii,2);
        
        plot(xcorr_all.lags,xcorr_all.trace_mn.math{e2 ,e1,1},'Color',col(i,:),'lineWidth',3)
        hold on
        plot(xcorr_all.lags,xcorr_all.permtrace_mn.math{e2 ,e1,1},'color',col(i,:),'lineWidth',1)
        plot(xcorr_all.lags,xcorr_all.permtrace_mn.math{e2 ,e1,1}+xcorr_all.permtrace_sd.math{e1 ,e2,1},'color',col(i,:),'lineWidth',1)
        plot(xcorr_all.lags,xcorr_all.permtrace_mn.math{e2 ,e1,1}-xcorr_all.permtrace_sd.math{e1 ,e2,1},'color',col(i,:),'lineWidth',1)
        xline(0)
        yline(0)
        set(gcf,'color','w')
        set(gca,'fontsize',16)
        xlabel('Lag (s)')
        ylabel('z-scored cross correlation')
        ylim([0 max(xcorr_all.trace_mn.math{e2 ,e1,1}) + 1])
        title([num2str(e2) '_' s_tmp.DK_long_josef{(s_tmp.chan_num == e2)} ' - ' num2str(e1) '_' s_tmp.DK_long_josef{(s_tmp.chan_num == e1)}], 'Interpreter', 'none')
        %     ylim([-1 6])
        %         title(['IPS > MFG'])
        %         savePNG(gcf, 300, fname)
        %         close all
    end
    fname = sprintf('%s%s_%s_%s_cross_corr_gen.png', dirs.figures_dir, subjects{i}, num2str(electrodes(1)),num2str(electrodes(2)));
    
end








for i = 1:length(subjects)
    s_tmp = el_pair(strcmp(el_pair.sbj_name, subjects{i}),:);
    xcorr_all = laggedCorrPerm(subjects{i},project_name,block_names,dirs,e1,e2,'HFB',[],'condNames',{'math'});
    
    
    
    
    electrodes = [];
    for ii = 1:length(regions)
        els = s_tmp(strcmp(s_tmp.DK_long_josef, regions{ii}),:);
        [~, idx] = max(els.sc1c2_tstat);
        electrodes(ii) = els.chan_num(idx);
    end
    e1 = electrodes(2);
    e2 = electrodes(1);
    block_names = BlockBySubj(subjects{i}, 'MMR');
    xcorr_all = laggedCorrPerm(subjects{i},project_name,block_names,dirs,e1,e2,'HFB',[],'condNames',{'math'});
    plot(xcorr_all.lags,xcorr_all.trace_mn.math{e1 ,e2,1},'Color',col(i,:),'lineWidth',4)
    hold on
    plot(xcorr_all.lags,xcorr_all.permtrace_mn.math{e1 ,e2,1},'color',col(i,:),'lineWidth',1)
    plot(xcorr_all.lags,xcorr_all.permtrace_mn.math{e1 ,e2,1}+xcorr_all.permtrace_sd.math{e1 ,e2,1},'color',col(i,:),'lineWidth',1)
    plot(xcorr_all.lags,xcorr_all.permtrace_mn.math{e1 ,e2,1}-xcorr_all.permtrace_sd.math{e1 ,e2,1},'color',col(i,:),'lineWidth',1)
    plot([0 0],ylim,'k-','linewidth',1)
    plot(xlim,[0 0],'k-','linewidth',1)
    set(gcf,'color','w')
    set(gca,'fontsize',16)
    xlabel('Lag (s)')
    ylabel('z-scored cross correlation')
    %     ylim([-1 6])
    title(['IPS > MFG'])
    fname = sprintf('%s%s_%s_%s_cross_corr_gen.png', dirs.figures_dir, subjects{i}, num2str(electrodes(1)),num2str(electrodes(2)));
    savePNG(gcf, 300, fname)
    close all
end


% Single trials


elect1 = [105 15 61 81 70];
for i1 = 1:length(elect1)
    figure('units', 'normalized', 'outerposition', [0 0 1 .4])
    e1 = elect1(i1);
    if e1 == 105
        elect2 = [15, 61, 81 70];
    elseif e1 == 15
        elect2 = [61, 81 70];
    elseif e1 == 61
        elect2 = [81 70];
    elseif e1 == 81
        elect2 = 70;
    end
    
    for i2 = 1:length(elect2)
        e2 = elect2(i2);
        if e1 ~= e2
            subplot(1,4,i2)
            col_plot = viridis(size(xcorr_all.trace.math{e1, e2, 1},1));
            for i = 1:size(xcorr_all.trace.math{e1, e2, 1},1)
                
                [V, I] = max(xcorr_all.trace.math{e1, e2, 1}(i,:));
                if ~isnan(V)
                    lag_per_trial(i) = xcorr_all.lags(I);
                else
                    lag_per_trial(i) = nan;
                end
                hold on
                plot(xcorr_all.lags,xcorr_all.trace.math{e1, e2, 1}(i,:), 'Color', col_plot(i,:))
                plot(lag_per_trial(i),V,'o', 'MarkerSize', 10, 'MarkerFaceColor', col_plot(i,:), 'MarkerEdgeColor', 'k')
            end
            plot(xcorr_all.lags, xcorr_all.trace_mn.math{e2, e1, 1}, 'LineWidth', 5, 'color', 'k')
            xlabel('Lag (s)')
            ylabel('single trials crosscorr with peak')
            set(gca,'fontsize',16)
            xticks([-5:1:5])
            box on
            grid on
        else
        end
    end
    save2pdf(['/Users/pinheirochagas/Pedro/drive/Stanford/grants/R03/resubmission/lagged_MMR_single_trials_elect_', num2str(e1) '.pdf'], gcf, 300)
    savePNG(gcf, 300, ['/Users/pinheirochagas/Pedro/drive/Stanford/grants/R03/resubmission/lagged_MMR_single_trials_elect_', num2str(e1) '.png'])
end








%%

for i = 1:size(el,1)
    elt = el(i,:);
    sbj_name = elt.names{1}
    block_names = BlockBySubj(sbj_name, project_name);
    
    [pairs_tmp idx] =  sort([elt.DK_lobe_1 elt.DK_lobe_2]);
    el1 = elt.(['chan_num_' num2str(idx(1))]);
    el2 = elt.(['chan_num_' num2str(idx(2))]);
    xcorr = laggedCorrPerm(sbj_name,project_name,block_names,dirs,el1, el2,'HFB',[],'condNames',{'math'});
    [x,y] = max(xcorr.trace_mn.math{el2,el1,1});
    xcorr_lag(i) = xcorr.lags(y);
    figure('units', 'normalized', 'outerposition', [0 0 0.5 0.25])
    hold on
    plot(xcorr.lags,xcorr.trace_mn.math{el2, el1},'Color','b','lineWidth',3)
    xlim([min(xcorr.lags) max(xcorr.lags)])
    xlabel('Lags (sec.)')
    ylabel('zscored correlation')
    xline(0, 'Color', [.5 .5 .5], 'LineWidth', 2)
    plot(xcorr.lags(y), xcorr.trace_mn.math{el2, el1}(y), 'r.', 'MarkerSize', 30)
    set(gca,'FontSize', 12)
    set(gca,'xtick',-4:0.25:4)
    grid on
    title_plot = [pairs_tmp{2},'_',num2str(el2),'_', pairs_tmp{1}, '_', num2str(el1)];
    title(title_plot, 'Interpreter', 'none')
    fname = [outdir, sbj_name, '_' title_plot, '_peak_xcorr.png'];
    savePNG(gcf,300,fname)
end
el.xcorr_lag = xcorr_lag';




%%


el_selectivity_MMR_ROL.ROL_math_avg = cellfun(@nanmean, el_selectivity_MMR_ROL.ROL_math_onsets)
el_selectivity_MMR_ROL.peak_math_avg = cellfun(@nanmean, el_selectivity_MMR_ROL.ROL_math_peaks)



labels_plot = {'FG', 'ITG', 'IPS', 'SPL', 'IFG', 'MFG', 'SFG'};
colors = hsv(7);
for i = 1:length(labels_plot)
    ax = subplot(1,7,i)
    el_temp = el_selectivity_MMR_ROL(strcmp(el_selectivity_MMR_ROL.DK_long_josef, labels_plot{i}),:);
    %     el_temp = el_temp(contains(el_temp.elect_select, {'math'}),:);
    el_temp = el_temp(~isnan(el_temp.peak_math_avg),:);
    scatter(el_temp.peak_math_avg, el_temp.sc1c2_tstat, 50, 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:))
    h = lsline(ax)
    h.LineWidth = 3;
    h.Color = 'k';
    xlabel('ROL (ms)')
    ylabel('t-stat')
    title(labels_plot{i})
    %     xlim([0 max(el_temp.ROL_math_avg)])
    [r, p] = corrcoef(el_temp.peak_math_avg, el_temp.sc1c2_tstat)
    if p(2) < 0.05
        r = num2str(r(2));
        text(max(el_temp.peak_math_avg), max(el_temp.sc1c2_tstat), r(1:6), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 16)
    else
    end
end
set(gcf,'color','w')
fname = sprintf('%scorr_time_to_peak_t-stat.png', dirs.figures_dir);
savePNG(gcf, 300, fname)q





%% Memory activations
el = el_selectivity_all_calc;
el = el(contains(el.elect_select, {'autobio only', 'autobio selective', 'episodic selective', 'episodic only'}),:)
el = el(:,contains(el.Properties.VariableNames, vars))

el_hipp_vatc = el(contains(el.DK_long_josef, {'HIPPOCAMPUS ANTERIOR', 'HIPPOCAMPUS MID', 'HIPPOCAMPUS POSTERIOR', 'vaTC'}),:);
el_dmn = el(strcmp(el.Yeo7, 'Default'),:);
el_memory = el(contains(el.elect_select, {'autobio', 'episodic'}) & ~contains(el.elect_select, {'autobio deact', 'episodic deact'} ),:)

vars = {'sbj_name', 'FS_label', 'chan_num', 'elect_select', 'DK_long_josef', 'Yeo7', 'LvsR', 'sEEG_ECoG'}


vars = {'LvsR','MNI_coord', 'WMvsGM', 'sEEG_ECoG', 'DK_lobe', 'Yeo7', 'Yeo17', 'DK_long_josef'};


%%

vars = {'chan_num', 'FS_label', 'LvsR','MNI_coord', 'LEPTO_coord', 'WMvsGM', 'sEEG_ECoG', 'DK_lobe', 'Yeo7', 'Yeo17', 'DK_long_josef', ...
    'elect_select', 'act_deact_cond1', 'act_deact_cond2', 'sc1c2_FDR', 'sc1b1_FDR' , 'sc2b2_FDR', ...
    'sc1c2_Pperm', 'sc1b1_Pperm', 'sc2b2_Pperm', 'sc1c2_tstat', 'sc1b1_tstat', 'sc2b2_tstat'};

task = 'MMR';
sinfo_MMR = sinfo(strcmp(sinfo.task, task),:);
el_selectivity_MMR = concat_elect_select(sinfo_MMR.sbj_name, task, dirs, vars);
task = 'UCLA';
sinfo_UCLA = sinfo(strcmp(sinfo.task, task),:);
el_selectivity_UCLA = concat_elect_select(sinfo_UCLA.sbj_name, task, dirs, vars);

el_selectivity_calc_sim = [el_selectivity_MMR;el_selectivity_UCLA]
el_selectivity_calc_sim = el_selectivity_calc_sim(strcmp(el_selectivity_calc_sim.WMvsGM, 'GM') | strcmp(el_selectivity_calc_sim.WMvsGM, 'WM'), :);
el_selectivity_calc_sim = el_selectivity_calc_sim(~strcmp(el_selectivity_calc_sim.Yeo7, 'FreeSurfer_Defined_Medial_Wall'),:)

task = 'Memoria';
sinfo_Memoria = sinfo(strcmp(sinfo.task, task),:);
el_selectivity_Memoria = concat_elect_select(sinfo_Memoria.sbj_name, task, dirs, vars);
el_selectivity_Memoria = el_selectivity_Memoria(strcmp(el_selectivity_Memoria.WMvsGM, 'GM') | strcmp(el_selectivity_Memoria.WMvsGM, 'WM'), :);
el_selectivity_Memoria = el_selectivity_Memoria(~strcmp(el_selectivity_Memoria.Yeo7, 'FreeSurfer_Defined_Medial_Wall'),:)


el_selectivity_all_calc = [el_selectivity_Memoria; el_selectivity_calc_sim];
el_selectivity = simplify_selectivity(el_selectivity_all_calc, 'MMR');


labels_plot = {'ITG', 'SPL', 'IPS', 'IFG', 'MFG', 'SFG'}
el_selectivity_math = el_selectivity(contains(el_selectivity.elect_select, {'math only', 'math selective'}),:);
coord_labels = {'X', 'Y', 'Z'}

colors = hsv(3);
count = 1

for i =  1:length(labels_plot)
    el_temp = el_selectivity_math(strcmp(el_selectivity_math.DK_long_josef, labels_plot{i}),:);
    el_temp.LEPTO_coord(:,1) = abs(el_temp.LEPTO_coord(:,1))
    for ii = 1: size(el_temp.LEPTO_coord,2)
        ax = subplot(6,3,count)
        scatter(el_temp.LEPTO_coord(:,ii), el_temp.sc1c2_tstat, 50, 'MarkerFaceColor', colors(ii,:), 'MarkerEdgeColor', colors(ii,:))
        h = lsline(ax);
        h.LineWidth = 3;
        h.Color = 'k';
        xlabel(coord_labels{ii})
        ylabel('t-stat')
        title(labels_plot{i})
        %     xlim([0 max(el_temp.ROL_math_avg)])
        [r, p] = corrcoef(el_temp.LEPTO_coord(:,ii), el_temp.sc1c2_tstat);
        if p(2) < 0.05
            r = num2str(r(2));
            text(max(el_temp.LEPTO_coord(:,ii)), max(el_temp.sc1c2_tstat), r(1:6), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 16)
        else
        end
        count = count + 1
    end
    
    
    
    %% 
end



%% Salience error


%% Univariate Selectivity
tag = 'stim';
task = 'MMR';
dirs = InitializeDirs(task, sinfo_MMR.sbj_name{1}, comp_root, server_root, code_root); % 'Pedro_NeuroSpin2T'
dirs.result_dir = result_dir;

cfg = [];
[cfg.column,cfg.conds] = getCondNames(task);

cfg.column = 'correctness';
cfg.conds = {'correct', 'incorrect'};

cfg.stats_params = genStatsParams(task);
cfg.stats_params.task_win = [1 2];

parfor i = 1:size(sinfo_MMR,1)
    ElecSelectivityAll(sinfo_MMR.sbj_name{i}, dirs, task, 'stim', 'Band', 'HFB', cfg)
end


subjects = sinfo_MMR.sbj_name;
subjects(strcmp(subjects, 'S19_137_AF')) = [];
subjects(strcmp(subjects, 'S20_149_DR')) = [];
subjects(strcmp(subjects, 'S19_137_AF')) = [];
subjects(strcmp(subjects, 'S19_137_AF')) = [];
vars = {'chan_num', 'FS_label', 'LvsR','MNI_coord', 'LEPTO_coord', 'WMvsGM', 'sEEG_ECoG', 'DK_lobe', 'Yeo7', 'Yeo17', 'DK_long_josef', ...
    'elect_select', 'act_deact_cond1', 'act_deact_cond2', 'sc1c2_FDR', 'sc1b1_FDR' , 'sc2b2_FDR', ...
    'sc1c2_Pperm', 'sc1b1_Pperm', 'sc2b2_Pperm', 'sc1c2_tstat', 'sc1b1_tstat', 'sc2b2_tstat'};

elect_select_all = concat_elect_select(subjects, task, dirs, vars, 'correctness');


elect = elect_select_all(~strcmp(elect_select_all.elect_select, 'no selectivity'),:);

elect_inc = elect_select_all(~strcmp(elect_select_all.elect_select, 'no selectivity'),:);











