%%
clc; 
clear all; close all;

% to open/modify the Matlab text editor sessions
addpath('./lib/EditorSessionManager');

% to generate colormaps
addpath('./lib/DrosteEffect-BrewerMap-221b913');

% to add definition of classes being used for this analysis
addpath('./lib/li_analysis');

% to include class definitions used in videoscorer_data
addpath('./lib/video_scorer/');
% rmpath('./lib/video_scorer/Source');

% to include higher order accurate differentiation function
addpath('./lib/diffxy');

% to include fmf reader
addpath('./lib/flymovieformat');

% to include hline and vline function
addpath('./lib/hline_vline');

% to include boxplot2
addpath('./lib/boxplot2-pkg/boxplot2');
addpath('./lib/boxplot2-pkg/minmax');
%% Find r* entry intervals

% Inputs
close all; clc;
% clear;

if isunix
    inputFile = '/media/reken001/Disk_12/honeybee_experiments/postprocessing/BlindLandingtracks_A4_LDF_rref.mat';
    outputFile = '/media/reken001/Disk_12/honeybee_experiments/postprocessing/BlindLandingtracks_A4_rrefEntry_with_rdot.mat';
    DirPlots = '/media/reken001/Disk_12/honeybee_experiments/postprocessing/plots/BlindLandingTracks/rref_entry_with_time';
    % DirPlots = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/plots/BlindLandingTracks/rref_entry';
    % DirPlots = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/plots/BlindLandingTracks/rref_estimate_3dspeed';
elseif ispc
    inputFile = 'D:/light_intensity_experiments/postprocessing/BlindLandingtracks_A1_rref.mat';
    outputFile = 'D:/light_intensity_experiments/postprocessing/BlindLandingtracks_A2_rrefEntry_with_rdot.mat';
%     DirPlots = 'D:/light_intensity_experiments/postprocessing/plots/BlindLandingTracks/rref_entry_with_time';
    DirPlots = 'D:/light_intensity_experiments/postprocessing/plots/BlindLandingTracks/rref_entry_with_time2';
    DirPlots = 'D:/light_intensity_experiments/postprocessing/plots/BlindLandingTracks/rref_entry_rdotsim';
end

load(inputFile);


delPreviousPlots = true; % BE CAREFUL - when set to true, all previously saved plots are deleted
savePlots = true;
savePDFs = false;


factors = [0.25:0.25:2.5];
factors = [1];

interval_for_rdot_estimate = [0.2 0.8]; % in percentage of rref
for ct_folder = 1:length(landingTracks)
            if isempty(landingTracks{ct_folder})
                continue;
            end
            
            % Delete previous plots
            if delPreviousPlots && savePlots && exist(fullfile(DirPlots, landingTracks{ct_folder}(1).pattern), 'dir')
                rmdir(fullfile(DirPlots, landingTracks{ct_folder}(1).pattern),'s');
            end
        
            % Create sub-directory for each treatment if it doesn't exist
            if savePlots && ~exist(fullfile(DirPlots, landingTracks{ct_folder}(1).pattern), 'dir')
                mkdir(fullfile(DirPlots, landingTracks{ct_folder}(1).pattern));
            end
            DirPlots_treatment = fullfile(DirPlots, landingTracks{ct_folder}(1).pattern);
        
            for ct_track=1:length(landingTracks{ct_folder}) % for each landing track
                track = landingTracks{ct_folder}(ct_track);
                for ct_excerpt=1:length(track.state_LDF) % for each track excerpt
                    excerpt = track.state_LDF(ct_excerpt);
                    
                    excerpt.find_rrefEntry();
                    excerpt.find_rdot_estimate_in_rrefEntry(interval_for_rdot_estimate);
                    
                    %%%% Check - Comment for normal running %%%
                    %                               if any((excerpt.rrefEntrySegments(4).delta_r > -excerpt.rrefEntrySegments(4).rref) & excerpt.rrefEntrySegments(4).isRise)
                    %                                   keyboard;
                    %                               end
                    %%%% Check finished %%%%%%%%%%%%%%%%%%%%%%%
                    
                    if savePlots
                        for ct_factor=1:length(factors)
                            plotHandles = excerpt.plot_rdotSimulation_with_actualdata(factors(ct_factor));
                            %                                       plotHandles = excerpt.plot_rrefsEntry_with_rdotestimate(factors(ct_factor));
                            
                            %                                       plotHandles = excerpt.plot_rrefsEntry_with_time(factors(ct_factor));
                            %                                       plotHandles = excerpt.plot_rrefsEntry(factors(ct_factor));
                            %                                       plotHandles = excerpt.plot_rrefs_with3dspeed(factors(ct_factor));
                            
                            if ~isempty(plotHandles)
                                
                                % Resizing the figures
                                for i=1:length(plotHandles)
                                    % for plot_rdotSimulation_with_actualdata
                                    plotHandles(i).Position(3) = 885;
                                    plotHandles(i).Position(4) = 820;
                                    
                                    % for rest
                                    %                                               plotHandles(i).Position(3) = 680;
                                    %                                               plotHandles(i).Position(4) = 545;
                                    
                                    if i==1
                                        figureName = ['fac_' num2str(factors(ct_factor),'%0.2f') '_' ...
                                            '_track' ...
                                            num2str(ct_track) '_excerpt' num2str(ct_excerpt) ...
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
                    
                end % for excerpt
            end % for track
            
end
keyboard
save(outputFile, 'landingTracks');
keyboard;