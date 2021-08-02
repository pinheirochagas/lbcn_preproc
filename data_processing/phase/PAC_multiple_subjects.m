function PAC_multiple_subjects(s,el_selectivity, pac_params, dirs)
    el_temp = el_selectivity(strcmp(el_selectivity.sbj_name, s),:);
    block_names = BlockBySubj(s, el_temp.task{1});
    electrodes = el_temp.chan_num;
    all_pairs = allcomb(electrodes, electrodes);
    phase_elecs = all_pairs(:,1)';
    amp_elecs =   all_pairs(:,2)';
    PAC = computePACAll(s,el_temp.task{1},block_names,dirs,phase_elecs,amp_elecs,[],'SpecDense','stim','condNames',{'autobio', 'math'},pac_params);
    fname = sprintf('%sPAC_MMR_%s.mat',dirs.pac, s);
    save(fname, 'PAC')
end

