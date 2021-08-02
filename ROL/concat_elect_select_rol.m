function elect_select_all = concat_elect_select_rol(subjects, task, hypothesis, dirs, vars)

vars_selec = {'elect_select', 'act_deact_cond1', 'act_deact_cond2', 'sc1c2_FDR', 'sc1b1_FDR' , 'sc2b2_FDR', ...
    'sc1c2_Pperm', 'sc1b1_Pperm', 'sc2b2_Pperm', 'sc1c2_tstat', 'sc1b1_tstat', 'sc2b2_tstat'};


elect_select_all = table;
for is = 1:length(subjects)
    s = subjects{is};
    fname = sprintf('%scomputation/%s/%s_%s_%s.mat',dirs.paper_results, hypothesis, task, s, hypothesis);
    
    
    if exist(fname) == 2
        load(fname)
        load([dirs.original_data filesep  s filesep 'subjVar_'  s '.mat']);
        subjVar = subjVar.elinfo(:, contains(subjVar.elinfo.Properties.VariableNames, vars));
        if sum(strcmp(subjVar.Properties.VariableNames, 'DK_long_josef')) == 0
            subjVar.DK_long_josef = repmat({'not done'},size(subjVar,1),1,1);
        else
        end
        el_selectivity_tmp = el_selectivity(:, contains(el_selectivity.Properties.VariableNames, vars_selec));
        el_selectivity_tmp.sbj_name = repmat({s},size(el_selectivity_tmp,1),1);
        
        
        % Load subject variable and include variables
        el_selectivity_tmp = [subjVar, el_selectivity_tmp];
        
        % add ROL
        fname = sprintf('%sROL/%s_%s_ROL.mat',dirs.paper_results, s, task);
        
        load(fname)
        cond_names = fieldnames(ROL);
        for ir1 = 1:length(cond_names)
            cond = cond_names{ir1};
            ROL_idx = fieldnames(ROL.(cond));
            for ir2 = 1:length(ROL_idx)
                ROL_tmp_name = sprintf('ROL_%s_%s', cond, ROL_idx{ir2});
                el_selectivity_tmp.(ROL_tmp_name) = ROL.(cond).(ROL_idx{ir2});
            end
        end
        
        elect_select_all = [elect_select_all; el_selectivity_tmp];
        
    else
    end
    
end

if sum(contains(vars, 'DK_lobe')) > 0
    elect_select_all.DK_lobe_generic = DK_lobe_generic(elect_select_all);
else
end


end



