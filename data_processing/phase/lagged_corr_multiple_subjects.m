function lagged_corr_multiple_subjects(s,el_selectivity, dirs)

    
    fname = sprintf('%slagged_corr_MMR_%s.mat',dirs.lagged_corr, s);
    if ~exist(fname, 'file')
        
        el_temp = el_selectivity(strcmp(el_selectivity.sbj_name, s),:);
        
        block_names = BlockBySubj(s, el_temp.task{1});
        e1 = el_temp.chan_num;
        e2 = el_temp.chan_num;
        
        xcorr = laggedCorrPerm(s,el_temp.task{1},block_names,dirs,e1,e2,'HFB',[],'condNames',{'math'});
        
        
   
    save(fname, 'xcorr')
    else
    end
end

