%% Display info about captured video instances

%%
clc; 
clear; close all;

% to generate colormaps
addpath('/home/reken001/Pulkit/lib/DrosteEffect-BrewerMap-221b913');

% to add definition of classes being used for this analysis
addpath('/home/reken001/Pulkit/lib/li_analysis');

% to include class definitions used in videoscorer_data
addpath('/home/reken001/Pulkit/lib/video_scorer/');
% rmpath('/home/reken001/Pulkit/lib/video_scorer/Source');

% to include higher order accurate differentiation function
addpath('/home/reken001/Pulkit/lib/diffxy');

% to include fmf reader
addpath('/home/reken001/Pulkit/lib/flymovieformat');

% to include hline and vline function
addpath('/home/reken001/Pulkit/lib/hline_vline');

% to include all custom defined class definitions
addpath('/home/reken001/Pulkit/MATLAB');

%% Find and save information about all video recordings
DataDir = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/';
VideoDir = '/media/reken001/Disk_08_backup/light_intensity_experiments/Videos';

inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A1.mat';
load(inputFile);

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
behaviour = {'rising','constant','sleeping'};

for ct_pattern = 2%1:length(pattern)
    for ct_light = 3%1:length(light)
        for ct_behaviour = 2%1:length(behaviour)
        
           % Selecting relevant treatments
            if strcmpi(behaviour{ct_behaviour}, 'rising')
                relevantTreatments = treatments(strcmpi({treatments.pattern}, pattern{ct_pattern}) & ...
                                     strcmpi({treatments.light}, light{ct_light}) & ...
                                     rem(1:length(treatments), 8)==1);
            elseif strcmpi(behaviour{ct_behaviour}, 'constant')
                relevantTreatments = treatments(strcmpi({treatments.pattern}, pattern{ct_pattern}) & ...
                                     strcmpi({treatments.light}, light{ct_light}) & ...
                                     rem(1:length(treatments), 8)>1 & ...
                                     rem(1:length(treatments), 8)<8);
            elseif strcmpi(behaviour{ct_behaviour}, 'sleeping')
                relevantTreatments = treatments(strcmpi({treatments.pattern}, pattern{ct_pattern}) & ...
                                     strcmpi({treatments.light}, light{ct_light}) & ...
                                     rem(1:length(treatments), 8)==0);
            else
                error('What other treatments did you perform dude?')
            end
            
%             appData_filename = ['GUIDE_appData_Blindtracks_' pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour} '.mat'];
%             if exist(fullfile(DataDir, appData_filename), 'file')
%                 load(fullfile(DataDir, appData_filename));
%             end
            
            for ct_treatment=1:length(relevantTreatments) % for each relevant treatment
                treatment = relevantTreatments(ct_treatment);
                
                % Go to the folder where videos should be stored and
                % extract relevant information about them
                
                datenum_str = num2str(treatment.datenum);
                treatmentVideosDir = fullfile(VideoDir, datenum_str(1:4), datenum_str(5:6), datenum_str(7:8));
                if exist(treatmentVideosDir, 'dir')
                    allVideoFiles_cam84 = dir(fullfile(treatmentVideosDir, '*84.fmf'));
                    allVideoFiles_cam85 = dir(fullfile(treatmentVideosDir, '*85.fmf'));
                    isVideoUseful_cam84 = cellfun(@(x) str2double(x(1:6))>=treatment.startTime && ...
                                            str2double(x(1:6))<=treatment.endTime,{allVideoFiles_cam84.name});
                    isVideoUseful_cam85 = cellfun(@(x) str2double(x(1:6))>=treatment.startTime && ...
                                            str2double(x(1:6))<=treatment.endTime,{allVideoFiles_cam85.name});
                    treatmentVideos_cam84 = allVideoFiles_cam84(isVideoUseful_cam84);
                    treatmentVideos_cam85 = allVideoFiles_cam85(isVideoUseful_cam85);
                    treatmentVideos = [treatmentVideos_cam84; treatmentVideos_cam85];
                    
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
end
%% Save file
save(inputFile, 'treatments');
