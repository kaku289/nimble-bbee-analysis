%% Display info about captured video instances

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
%% Find and save information about all video recordings
DataDir = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/';
VideoDir = '/media/reken001/Disk_07/steady_wind_experiments/Videos';

inputFile = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/BlindLandingtracks_A3_LDF_rref_rrefEntry.mat';
% load(inputFile);

winds = unique([treatments.wind]);
behaviour = {'rising','constant','sleeping'};

for ct_wind = 6%1:length(winds)
    for ct_behaviour = 2%1:length(behaviour)
        disp(' ');
        disp([' wind: ' num2str(winds(ct_wind)) ...
            ', behaviour: ' behaviour{ct_behaviour}]);
        
        % Selecting relevant treatments
        if strcmpi(behaviour{ct_behaviour}, 'rising')
            relevantTreatments = treatments(rem(1:length(treatments), 8)==1);
        elseif strcmpi(behaviour{ct_behaviour}, 'constant')
            relevantTreatments = treatments( [treatments.wind] == winds(ct_wind) & ...
                rem(1:length(treatments), 8)>1 & ...
                rem(1:length(treatments), 8)<8);
        elseif strcmpi(behaviour{ct_behaviour}, 'sleeping')
            relevantTreatments = treatments(rem(1:length(treatments), 8)==0);
        else
            error('What other treatments did you perform dude?')
        end
        
        % Only analyse treatments that have uniform wind measurements
        % avaliable throughout their course of running
        hasUniformHwData = arrayfun(@(x) x.hwData.hasUniformHwData,relevantTreatments);
        relevantTreatments = relevantTreatments(hasUniformHwData);
        
        %             appData_filename = ['GUIDE_appData_Blindtracks_' pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour} '.mat'];
        %             if exist(fullfile(DataDir, appData_filename), 'file')
        %                 load(fullfile(DataDir, appData_filename));
        %             end
        
        for ct_treatment=1:length(relevantTreatments) % for each relevant treatment
            disp(['Into treatment: ' num2str(ct_treatment)]);
            treatment = relevantTreatments(ct_treatment);
            
            % Go to the folder where videos should be stored and
            % extract relevant information about them
            
            datenum_str = num2str(treatment.datenum);
            treatmentVideosDir = fullfile(VideoDir, datenum_str(1:4), datenum_str(5:6), datenum_str(7:8));
            if exist(treatmentVideosDir, 'dir')
                allVideoFiles_cam84 = dir(fullfile(treatmentVideosDir, '*84.fmf'));
                %                     allVideoFiles_cam85 = dir(fullfile(treatmentVideosDir, '*85.fmf'));
                isVideoUseful_cam84 = cellfun(@(x) str2double(x(1:6))>=treatment.startTime && ...
                    str2double(x(1:6))<=treatment.endTime,{allVideoFiles_cam84.name});
                %                     isVideoUseful_cam85 = cellfun(@(x) str2double(x(1:6))>=treatment.startTime && ...
                %                                             str2double(x(1:6))<=treatment.endTime,{allVideoFiles_cam85.name});
                treatmentVideos_cam84 = allVideoFiles_cam84(isVideoUseful_cam84);
                %                     treatmentVideos_cam85 = allVideoFiles_cam85(isVideoUseful_cam85);
                %                     treatmentVideos = [treatmentVideos_cam84; treatmentVideos_cam85];
                treatmentVideos = treatmentVideos_cam84;
                
                treatment.videosInfo = recordedVideosInformation.empty;
                
                if ~isempty(treatmentVideos)
                    treatment.videosInfo(length(treatmentVideos)) = recordedVideosInformation();
                    for ct_video=1:length(treatmentVideos)
                        treatment.videosInfo(ct_video).extractInformation(treatmentVideos(ct_video));
                    end
                end
            else
                disp(['Videos do NOT exist for treatment: ' datenum_str ', ' num2str(treatment.startTime) '-' num2str(treatment.endTime) ]);
            end
            %                 keyboard;
        end
    end
end
%% Save file
outputFile = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/BlindLandingtracks_A3_LDF_rref_rrefEntry_videos.mat';
save(outputFile, 'treatments', '-v7.3', '-nocompression');

% save(outputFile, 'treatments');
