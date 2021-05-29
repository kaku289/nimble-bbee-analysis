%% Loading data and creating plots for Introduction figure
% Plots are plotted for each track and highlights rref segment(s), entry
% segment(s), filtered data, and data estimated from system identification
clc; clear; close all;
if isunix
    dataFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A2_rrefEntryEstimation_everything_f1.mat';
    dataFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A2_rrefEntryEstimation_everything_f1pt5_all.mat';
    DirPlots = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/plots/BlindLandingTracks/articleA2_plots';
elseif ispc
    dataFile = 'D:/light_intensity_experiments/postprocessing/BlindLandingtracks_A2_rrefEntryEstimation_everything_f1.mat';
    dataFile = 'D:/light_intensity_experiments/postprocessing/BlindLandingtracks_A2_rrefEntryEstimation_everything_f1pt5_all.mat';
    DirPlots = 'D:/light_intensity_experiments/postprocessing/plots/BlindLandingTracks/articleA2_plots';
end
% load(dataFile);

%% For f=1 (based on track numbers)
delPreviousPlots = true; % BE CAREFUL - when set to true, all previously saved plots are deleted
savePlots = true;
savePDFs = true;

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
chosen_fac = 1;
% behaviour = {'rising','constant','sleeping'};
selected_tracks = [151 166 388];

assert(length(pattern)==size(data_all,1) && length(light)==size(data_all,2));
for ct_pattern=2%1:length(pattern)
    for ct_light=3%1:length(light)
        % Delete previous plots
        if delPreviousPlots && savePlots && exist(fullfile(DirPlots, [pattern{data_all(ct_pattern,ct_light).ct_pattern} '_' light{data_all(ct_pattern,ct_light).ct_light}]), 'dir')
            rmdir(fullfile(DirPlots, [pattern{data_all(ct_pattern,ct_light).ct_pattern} '_' light{data_all(ct_pattern,ct_light).ct_light}]),'s');
        end
        
        % Create sub-directory for each treatment if it doesn't exist
        if savePlots && ~exist(fullfile(DirPlots, [pattern{data_all(ct_pattern,ct_light).ct_pattern} '_' light{data_all(ct_pattern,ct_light).ct_light}]), 'dir')
            mkdir(fullfile(DirPlots, [pattern{data_all(ct_pattern,ct_light).ct_pattern} '_' light{data_all(ct_pattern,ct_light).ct_light}]));
        end
        DirPlots_treatment = fullfile(DirPlots, [pattern{data_all(ct_pattern,ct_light).ct_pattern} '_' light{data_all(ct_pattern,ct_light).ct_light}]);
        
        for ct2=1:length(selected_tracks)
            ct1 = selected_tracks(ct2);
%             % create tf vector in order of rrefEntrySegments
            tfs = data4est_lowpass{ct_pattern,ct_light}.tfest(:,:,data4est_lowpass{ct_pattern,ct_light}.track_indexes == ct1, 2);
            sysiddata = data4est_lowpass{ct_pattern,ct_light}.iddata(data4est_lowpass{ct_pattern,ct_light}.track_indexes == ct1);
            
            % plot data
%             plotHandles = data_all(ct_pattern,ct_light).tracks_fac(ct1).plot_rrefsEntry_with_acc(chosen_fac);
%             data_all(ct_pattern,ct_light).tracks_fac(ct1).plot_rrefsEntry_withSimulatedData(chosen_fac, tfs);
            plotHandles = data_all(ct_pattern,ct_light).tracks_fac(ct1).plot_rrefsEntry_withActualFilteredEstimatedData2(chosen_fac, sysiddata, tfs)
            if ~isempty(plotHandles)
                
                % Resizing the figures
                for i=1:length(plotHandles)
%                     plotHandles(i).Position(3) = 560;
%                     plotHandles(i).Position(4) = 560;
                    
                    % For sys-id simulation plots
                    plotHandles(i).Position(3) = 450;
                    plotHandles(i).Position(4) = 650;
                    
                    if i==1
                        figureName = ['fac_' num2str(chosen_fac,'%0.2f') ...
                            '_track' ...
                            num2str(ct1) ...
                            '.png'];
                    end
                    
                    saveas(plotHandles(i), fullfile(DirPlots_treatment, figureName) ,'png');
                    
                    if savePDFs
                        print(plotHandles(i), strrep(fullfile(DirPlots_treatment, figureName), ...
                            '.png','.pdf'), '-dpdf');
                    end
                end
                
                close(plotHandles);
            end
        end
    end
end

%% For f=1.5 (based on track identity - more general/easier than above)
delPreviousPlots = true; % BE CAREFUL - when set to true, all previously saved plots are deleted
savePlots = true;
savePDFs = true;

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
chosen_fac = 1.5;
% behaviour = {'rising','constant','sleeping'};
selected_tracks = {'20190704_153000_170000_obj7984_track234_excerpt1', '20190704_153000_170000_obj9821_track290_excerpt1', '20190708_140000_153000_obj664_track13_excerpt1'};

assert(length(pattern)==size(data_all,1) && length(light)==size(data_all,2));
for ct_pattern=2%1:length(pattern)
    for ct_light=3%1:length(light)
        % Delete previous plots
        if delPreviousPlots && savePlots && exist(fullfile(DirPlots, [pattern{data_all(ct_pattern,ct_light).ct_pattern} '_' light{data_all(ct_pattern,ct_light).ct_light}]), 'dir')
            rmdir(fullfile(DirPlots, [pattern{data_all(ct_pattern,ct_light).ct_pattern} '_' light{data_all(ct_pattern,ct_light).ct_light}]),'s');
        end
        
        % Create sub-directory for each treatment if it doesn't exist
        if savePlots && ~exist(fullfile(DirPlots, [pattern{data_all(ct_pattern,ct_light).ct_pattern} '_' light{data_all(ct_pattern,ct_light).ct_light}]), 'dir')
            mkdir(fullfile(DirPlots, [pattern{data_all(ct_pattern,ct_light).ct_pattern} '_' light{data_all(ct_pattern,ct_light).ct_light}]));
        end
        DirPlots_treatment = fullfile(DirPlots, [pattern{data_all(ct_pattern,ct_light).ct_pattern} '_' light{data_all(ct_pattern,ct_light).ct_light}]);
        
        % Gather details of all landing tracks
        datenum = [data_all(ct_pattern,ct_light).landingTrack.datenum];
        obj_id = [data_all(ct_pattern,ct_light).landingTrack.obj_id];
        
        for ct2=1:length(selected_tracks)
            selected_track_str = selected_tracks{ct2};
            selected_track_str_parts = strsplit(selected_track_str,'_');
            
            indx = find(datenum == str2num(selected_track_str_parts{1}) & obj_id == str2num(selected_track_str_parts{4}(4:end)));
            assert(length(indx)==1); % If not, tighten the indx selection with starting time of the treatment
            
            ct1 = indx;
%             % create tf vector in order of rrefEntrySegments
            tfs = data4est_lowpass{ct_pattern,ct_light}.tfest(:,:,data4est_lowpass{ct_pattern,ct_light}.track_indexes == ct1, 2);
            sysiddata = data4est_lowpass{ct_pattern,ct_light}.iddata(data4est_lowpass{ct_pattern,ct_light}.track_indexes == ct1);
            
            % plot data
%             plotHandles = data_all(ct_pattern,ct_light).tracks_fac(ct1).plot_rrefsEntry_with_acc(chosen_fac);
%             data_all(ct_pattern,ct_light).tracks_fac(ct1).plot_rrefsEntry_withSimulatedData(chosen_fac, tfs);
            plotHandles = data_all(ct_pattern,ct_light).tracks_fac(ct1).plot_rrefsEntry_withActualFilteredEstimatedData2(chosen_fac, sysiddata, tfs)
            if ~isempty(plotHandles)
                
                % Resizing the figures
                for i=1:length(plotHandles)
%                     plotHandles(i).Position(3) = 560;
%                     plotHandles(i).Position(4) = 560;
                    
                    % For sys-id simulation plots
                    plotHandles(i).Position(3) = 450;
                    plotHandles(i).Position(4) = 650;
                    
                    if i==1
                        figureName = ['fac_' num2str(chosen_fac,'%0.2f') ...
                            '_track' ...
                            num2str(ct1) ...
                            '.png'];
                    end
                    
                    saveas(plotHandles(i), fullfile(DirPlots_treatment, figureName) ,'png');
                    
                    if savePDFs
                        print(plotHandles(i), strrep(fullfile(DirPlots_treatment, figureName), ...
                            '.png','.pdf'), '-dpdf');
                    end
                end
                
                close(plotHandles);
            end
        end
    end
end