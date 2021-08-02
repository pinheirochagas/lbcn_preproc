 function PlotModulation(dirs, subjVar, cfg)

%% Load comon brain
% elinfo = subjVar.elinfo;
elinfo = subjVar.elinfo;


% Decide if project left
if cfg.project_left
    elinfo.MNI_coord(:,1) = -abs(elinfo.MNI_coord(:,1));
    elinfo.LvsR = repmat({'L'}, size(elinfo,1),1);
end

fsaverage_dir = '/Applications/freesurfer/subjects/fsaverage/surf'; % correct that:'/Applications/freesurfer/freesurfer/subjects/fsaverage/surf/rh.pial'
if strcmp(cfg.Cortex, 'MNI')
    [cmcortex.right.vert cmcortex.right.tri]=read_surf(fullfile('/Applications/freesurfer/subjects/fsaverage/surf',['rh.' 'pial']));
    [cmcortex.left.vert cmcortex.left.tri]=read_surf(fullfile('/Applications/freesurfer/subjects/fsaverage/surf',['lh.' 'pial']));    
    coords_plot = elinfo.MNI_coord;
    
elseif  strcmp(cfg.Cortex, 'MNI_inflated')
    [cmcortex.right.vert cmcortex.right.tri]=read_surf(fullfile('/Applications/freesurfer/subjects/fsaverage/surf',['rh.' 'inflated_avg']));
    [cmcortex.left.vert cmcortex.left.tri]=read_surf(fullfile('/Applications/freesurfer/subjects/fsaverage/surf',['lh.' 'inflated_avg']));
    coords_plot = elinfo.fsaverageINF_coord;
    
elseif  strcmp(cfg.Cortex, 'native')
    cmcortex = subjVar.cortex;
    coords_plot = elinfo.LEPTO_coord;
else
    error('you must specify the cortical space to plot, either MNI or native.')
end



% basic parameters:
decimate = true;
final_fs = 50;

% Get color indices
% [col_idx,colors_plot] = colorbarFromValues(ind, 'RedBlue', [], true);

if ~isempty(cfg.ind)
%     [col_idx,colors_plot] = colorbarFromValues(cfg.ind, cfg.Colormap, cfg.clim, cfg.color_center_zero);
    [col_idx, colors_plot] = vals2colormap(cfg.ind, cfg.Colormap, cfg.clim);
    col_idx = 1:size(colors_plot,1);

%      col_idx(col_idx==0)=1; % dirty fix
    MarkerEdgeColor = [.3 .3 .3];
%     colors_plot = flip(colors_plot);
%     MarkerEdgeColor = 'none';
elseif size(cfg.MarkerColor,1) == 1
    col_idx = ones(size(elinfo,1),1);
    colors_plot = repmat(cfg.MarkerColor, size(elinfo,1), 1);
    MarkerEdgeColor = cfg.MarkerEdgeColor;
elseif size(cfg.MarkerColor,1) > 1 
    col_idx = 1:size(elinfo,1);
    colors_plot = cfg.MarkerColor;
    MarkerEdgeColor = cfg.MarkerEdgeColor;
%     MarkerEdgeColor =  'none';
end




%% Plot electrodes as dots in native space
if cfg.MarkerSize_mod
%     marker_size = abs(cfg.ind).*cfg.MarkerSize; 
    marker_size = abs(cfg.ind).*round(cfg.MarkSizeEffect*abs(cfg.ind));
    marker_size(marker_size<=0) = 0.0001;
else
%     marker_size = repmat(cfg.MarkerSize,size(elinfo,1),1);
     marker_size = cfg.MarkerSize;

end


% figureDim = [0 0 1 .4];

% col_label = repmat(cfg.col_label, size(elinfo,1),1);
if cfg.bad_chans
    cfg.col_label(cfg.bad_chans,:) = repmat([1 0 0], length(cfg.bad_chans),1);
    colors_plot(cfg.bad_chans,:) = repmat([1 1 1], length(cfg.bad_chans),1);
else
    
end


% f1 = figure('units', 'normalized', 'outerposition', cfg.figureDim);
views =  cfg.views;
hemis = cfg.hemis;

for i = 1:length(views)
     subplot(cfg.subplots(1),cfg.subplots(2),i)
    coords_plot = CorrectElecLoc(coords_plot, views{i}, hemis{i}, elinfo.sEEG_ECoG, cfg.CorrectFactor); %
    ctmr_gauss_plot(cmcortex.(hemis{i}),[0 0 0], 0, hemis{i}, views{i})
    
    for ii = 1:size(coords_plot,1)
%         if ii ~= cfg.chan_highlight
            % Only plot on the relevant hemisphere
            if ~isnan(col_idx(ii))
                if (strcmp(hemis{i}, 'left') == 1 && strcmp(elinfo.LvsR(ii), 'R') == 1) || (strcmp(hemis{i}, 'right') == 1 && strcmp(elinfo.LvsR(ii), 'L') == 1)
                else
%                     plot3(coords_plot(ii,1),coords_plot(ii,2),coords_plot(ii,3), 'o', 'MarkerSize', marker_size, 'MarkerFaceColor', colors_plot(col_idx(ii),:), 'MarkerEdgeColor', MarkerEdgeColor);
                    plot3(coords_plot(ii,1),coords_plot(ii,2),coords_plot(ii,3), 'o', 'MarkerSize', marker_size(ii), 'MarkerFaceColor', colors_plot(col_idx(ii),:), 'MarkerEdgeColor', MarkerEdgeColor);
                end
            else
            end
%         else
%         end
    end
    
    if ~isempty(cfg.chan_highlight)
        for hi = 1:length(cfg.chan_highlight)
            for ii = 1:size(coords_plot,1)
                % Only plot on the relevant hemisphere
                if ~isnan(col_idx(ii))
                    if (strcmp(hemis{i}, 'left') == 1 && strcmp(elinfo.LvsR(ii), 'R') == 1) || (strcmp(hemis{i}, 'right') == 1 && strcmp(elinfo.LvsR(ii), 'L') == 1)
                    else
                        plot3(coords_plot(cfg.chan_highlight(hi),1),coords_plot(cfg.chan_highlight(hi),2),coords_plot(cfg.chan_highlight(hi),3), 'o', 'MarkerSize', cfg.marker_size_high, 'MarkerFaceColor', cfg.highlight_face_col(hi,:), 'MarkerEdgeColor', cfg.highlight_edge_col(hi,:));
                        %                 plot3(coords_plot(ii,1),coords_plot(ii,2),coords_plot(ii,3), 'o', 'MarkerSize', marker_size(ii), 'MarkerFaceColor', colors_plot(col_idx(ii),:), 'MarkerEdgeColor', MarkerEdgeColor);
                    end
                else
                end
            end
        end
    else
    end
%     if ~isempty(cfg.chan_highlight)
%         for hi = 1:length(cfg.chan_highlight)
%             if (strcmp(hemis{i}, 'left') && strcmpi(elinfo.LvsR{cfg.chan_highlight(hi)},'L'))
%                 plot3(coords_plot(cfg.chan_highlight(hi),1),coords_plot(cfg.chan_highlight(hi),2),coords_plot(cfg.chan_highlight(hi),3), 'o', 'MarkerSize', cfg.marker_size_high, 'MarkerFaceColor', cfg.highlight_face_col(hi,:), 'MarkerEdgeColor', cfg.highlight_edge_col(hi,:));
% %                 text(coords_plot(cfg.chan_highlight(hi),1),coords_plot(cfg.chan_highlight(hi),2),coords_plot(cfg.chan_highlight(hi),3), cfg.colum_label{hi}, 'FontSize', cfg.label_font_size, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Interpreter', 'none');
%             elseif (strcmp(hemis{i}, 'right') && strcmpi(elinfo.LvsR{cfg.chan_highlight(hi)},'R'))
%                 plot3(coords_plot(cfg.chan_highlight(hi),1),coords_plot(cfg.chan_highlight(hi),2),coords_plot(cfg.chan_highlight(hi),3), 'o', 'MarkerSize', cfg.marker_size_high, 'MarkerFaceColor', cfg.highlight_face_col(hi,:), 'MarkerEdgeColor', cfg.highlight_edge_col(hi,:));
% %                 text(coords_plot(cfg.chan_highlight(hi),1),coords_plot(cfg.chan_highlight(hi),2),coords_plot(cfg.chan_highlight(hi),3), cfg.colum_label{hi}, 'FontSize', cfg.label_font_size, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Interpreter', 'none');
%             else
%             end
%         end
%     end
    
    
        if cfg.plot_label 
            for ii = 1:size(coords_plot,1)
                hold on
                if (strcmp(hemis{i}, 'left') == 1 && strcmp(elinfo.LvsR(ii), 'R') == 1) || (strcmp(hemis{i}, 'right') == 1 && strcmp(elinfo.LvsR(ii), 'L') == 1)
                else
                    if strcmp(cfg.colum_label, 'chan_num')
                        label = num2str(elinfo.(cfg.colum_label)(ii));
                    else
                        if length(cfg.colum_label) > 1
                            label = [];
                            for il = 1:length(cfg.colum_label)
                                label{il}= elinfo.(cfg.colum_label{il}){ii};
                            end
                            label = strjoin(label, '_');
                        else
                        end
                    end
                    text(coords_plot(ii,1),coords_plot(ii,2),coords_plot(ii,3), label, 'FontSize', cfg.label_font_size, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Interpreter', 'none');
                end
            end
        else
        end
    
    alpha(cfg.alpha)
end
if cfg.save
    savePNG(gcf, 300, [dirs.result_root filesep 'selectivity' filesep 'group_selectivity4_' cfg.project_name '_' cortex_space '.png']); % ADD TASK AND CONDITION
    close all
else
end

end