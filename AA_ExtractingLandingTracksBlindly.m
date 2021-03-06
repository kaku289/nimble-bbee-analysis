%% To "blindly" select landing tracks from the whole dataset

%%
clc; clear; close all;

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

%% User inputs
rootDir = '/media/reken001/Disk_08_backup/light_intensity_experiments/';

%% Initializing variables
metadataDir = 'Metadata';
calibDir = 'Calibration/_xml';
tempRHDir = 'TempRH';
trackingDir = 'postprocessing';
videosDir = 'Videos';
videoScoringDir = 'postprocessing';


file2d = 'data2d_distorted.csv';
fileDataAssociation = 'data_association.csv';
file3d = 'h5_flying_kalman_estimates.csv';
filecalib = 'calibration.xml';

%% Creating treatments array
% % Read treatment data from xlsx file
treatmentSchedule = readcell(fullfile(rootDir, trackingDir, 'treatment_schedule.xlsx'),'FileType','spreadsheet','Sheet','light','Range','C2:R11');
startTimes = readmatrix(fullfile(rootDir, trackingDir, 'treatment_schedule.xlsx'),'FileType','spreadsheet','Sheet','light','Range','A4:A11');
endTimes = readmatrix(fullfile(rootDir, trackingDir, 'treatment_schedule.xlsx'),'FileType','spreadsheet','Sheet','light','Range','B4:B11');

% % Read manually clicked 2-D disc centers 
discCenters_2d.Hive = readcell(fullfile(rootDir, trackingDir, 'Disc_centers.xlsx'),'FileType','spreadsheet','Sheet','light','Range','A2:I20');
discCenters_2d.Feeder = readcell(fullfile(rootDir, trackingDir, 'Disc_centers.xlsx'),'FileType','spreadsheet','Sheet','light','Range','A24:I42');
cameras = {'Basler_22549584', 'Basler_22549585', 'Basler_22549587', 'Basler_22956425'};

% Creating treatments array
treatments = Lighttreatment.empty;
for ct_day=1:size(treatmentSchedule,2)
    for ct_treatment=1:8 % 8 treatments in a day including rising and sleeping
        
        treatments(end+1) = Lighttreatment(treatmentSchedule{1,ct_day}, ...
                            treatmentSchedule{2,ct_day}, ...
                            treatmentSchedule{ct_treatment+2,ct_day}, ...
                            startTimes(ct_treatment), endTimes(ct_treatment));     
                        
        % Adding landing discs to each treatment
        sides = fieldnames(discCenters_2d); % cell array containing names of the sides where landing platform are present (i.e. Hive or Feeder)
        nSides = length(sides);
        for i=1:nSides
            side = sides{i};
            cam_ids_and_points2D = {4,2};
            for ct_cam = 1:4
                indx_row = [false false false [discCenters_2d.(side){4:end,1}] == treatments(end).datenum];
                indx_col = [false strcmpi(discCenters_2d.(side)(2,2:end), cameras{ct_cam})];
                cam_ids_and_points2D{ct_cam, 1} = cameras{ct_cam};
                cam_ids_and_points2D{ct_cam, 2} = [discCenters_2d.(side){indx_row, indx_col}];
            end
            treatments(end).addLandingDisc(cam_ids_and_points2D, side);
%             keyboard;
        end
        
        
    end
    
end

%% Finding landing tracks in each treatment
y_margin = 0.04; % in m
time_margin = 1.5; % in s

% Getting list of tracking files
tracking_files = dir(fullfile(rootDir, trackingDir, ['*mainbrain.offline']));

filenums = num2cell(cellfun(@(x) str2double(x(10:15)), {tracking_files(:).name}));
[tracking_files.filenum] = filenums{:}; % length of folders is not same as treatments array!

datenums = num2cell(cellfun(@(x) str2double(x(1:8)), {tracking_files(:).name}));
[tracking_files.datenum] = datenums{:}; % length of folders is not same as treatments array!

for ct_treatment=1:length(treatments) % for each treatment
    treatment = treatments(ct_treatment); % Always remember working by reference here
    
    
    
    % Find tracking files (also in case multiple files exist)
    % corresponding to this treatment
    treatment.trackingFiles = tracking_files(...
        [tracking_files(:).filenum] >= treatment.startTime & ...
        [tracking_files(:).filenum] <= treatment.endTime & ...
        [tracking_files(:).datenum] == treatment.datenum);
    
    if isempty(treatment.trackingFiles)
        continue;
    end
    
    % Pick any treatment tracking file - as all of them will
    % have same calibration file
    treatment.calib = FlydraCalibration(fullfile(treatment.trackingFiles(1).folder, treatment.trackingFiles(1).name, filecalib));
    
    % compute center of all landing discs for this treatment
    for ct_landingDiscs = 1:length(treatment.landingDiscs)
        treatment.landingDiscs(ct_landingDiscs).computeCenter(treatment.calib);
        if strcmpi(treatment.landingDiscs(ct_landingDiscs).side,'Hive')
           y_hiveDisc = treatment.landingDiscs(ct_landingDiscs).center(2);
           center_hiveDisc = treatment.landingDiscs(ct_landingDiscs).center;
           radius_hiveDisc = treatment.landingDiscs(ct_landingDiscs).radius;
        elseif strcmpi(treatment.landingDiscs(ct_landingDiscs).side,'Feeder')
           y_feederDisc = treatment.landingDiscs(ct_landingDiscs).center(2);
           center_feederDisc = treatment.landingDiscs(ct_landingDiscs).center;
           radius_feederDisc = treatment.landingDiscs(ct_landingDiscs).radius;
        end
    end
    disc_centers = [treatment.landingDiscs.center];
    landingDiscs = treatment.landingDiscs;
    boundingBox = [min(disc_centers(1,:))-2*treatment.landingDiscs(1).radius ...
                   max(disc_centers(1,:))+2*treatment.landingDiscs(1).radius ... 
                   min(disc_centers(2,:)) ...
                   max(disc_centers(2,:)) ... 
                   0.0 ...
                   0.48]; 
    
                   
    for ct_file=1:length(treatment.trackingFiles) % for each tracking file
        trackingFile = treatment.trackingFiles(ct_file);
        
        % load h5_flying_kalman_estimates
        data3d = readmatrix(fullfile(rootDir, trackingDir, trackingFile.name, file3d), 'Range', 2); % loading from second row onwards
        
        unique_obj_ids = unique(data3d(:,1));
        nObjects = length(unique_obj_ids);
        
        obj_ids = data3d(:,1);
        
         % Initializing the landing tracks
        landingTracks = BlindLandingtrack.empty(nObjects, 0);
        
        % % % Find landing tracks
        % do parallel search for all object ids in this data3d
        pattern = treatment.pattern;
        light = treatment.light;
        foldername = trackingFile.name;
        datenum = treatment.datenum;
        parfor ct_id=1:nObjects
%         for ct_id=1:nObjects
            obj_id = unique_obj_ids(ct_id);
            track = sortrows(data3d(obj_ids==obj_id, :), 3);
            track = track(IsInsideBox(track(:,4:6), boundingBox), :);
            
            landingTracks(ct_id) = BlindLandingtrack(pattern, light, foldername, datenum, obj_id);
            landingTracks_per_object = landingTracks(ct_id);
            
            % % Extract track excerpt that corresponds to landing tracks
            % find indices of all points that intersect 0.24 m
            % This only selects landing tracks that go through y=0.24 m
            % plane and precludes the ones that have shorter length of
            % approach
%             intersect_indices = find((track(1:end-1,5) < 0.24 & ...
%                                       track(2:end,5) >= 0.24 & ...
%                                       track(1:end-1,6) >= 0.02) | ...
%                                       (track(1:end-1,5) > 0.24 & ...
%                                       track(2:end,5) <= 0.24 & ...
%                                       track(1:end-1,6) >= 0.02));
                                  
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% OR create union of intersect
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% indices closer to both sides (10
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% cms from the platforms)
            
            intersect_indices_nearHive = find((track(1:end-1,5) < y_hiveDisc-2.5*y_margin & ...
                                        track(2:end,5) >= y_hiveDisc-2.5*y_margin & ...
                                        track(1:end-1,6) >= 0.02) | ...
                                        (track(1:end-1,5) > y_hiveDisc-2.5*y_margin & ...
                                        track(2:end,5) <= y_hiveDisc-2.5*y_margin & ...
                                        track(1:end-1,6) >= 0.02));
                                    
            intersect_indices_nearFeeder = find((track(1:end-1,5) > y_feederDisc+2.5*y_margin & ...
                                        track(2:end,5) <= y_feederDisc+2.5*y_margin & ...
                                        track(1:end-1,6) >= 0.02) | ...
                                        (track(1:end-1,5) < y_feederDisc+2.5*y_margin & ...
                                        track(2:end,5) >= y_feederDisc+2.5*y_margin & ...
                                        track(1:end-1,6) >= 0.02));

            intersect_indices = sort(unique([intersect_indices_nearHive; intersect_indices_nearFeeder]));
                                    
                                  
            % Look into 1.5s ahead and behind for each index to see if they cross y_margin away from either platforms 
            for ct=1:length(intersect_indices)
                time_intersect = track(intersect_indices(ct),3);
                y_track_ahead = track(track(:,3) >= time_intersect & track(:,3) <= time_intersect+time_margin, 5);
                y_track_behind = track(track(:,3) >= time_intersect-time_margin & track(:,3) <= time_intersect, 5);
%                 z_track_ahead = track(track(:,3) >= time_intersect & track(:,3) <= time_intersect+time_margin, 6);
%                 z_track_behind = track(track(:,3) >= time_intersect-time_margin & track(:,3) <= time_intersect, 6);
                if any(y_track_ahead > y_hiveDisc-y_margin)
                    % track excerpt moving towards hive
                    
                    % % Now choosing track excerpt
                    
                    % min_indx based on min y
                    [~, min_indx] = min(y_track_behind);
                    min_indx = intersect_indices(ct) - (length(y_track_behind) - min_indx);
                    
%                     % min_indx based on first point that is off the ground
%                     min_indx1 = find(z_track_behind>0.015,1);
%                     if isempty(min_indx2) % flight starts in track_ahead
%                         min_indx1 = find(z_track_ahead>0.015,1);
%                         if isempty(min_indx1) % flight never started
%                             disp('Entire track is a walking track');
%                             continue;
%                         else 
%                             min_indx = intersect_indices(ct) + min_indx1;
%                         end
%                     else % flight starts in track behind
%                         min_indx1 = intersect_indices(ct) - (length(z_track_behind) - min_indx1);
%                         if min_indx1 < min
%                         end
%                     end
                        
                            
                        
                    
                    
                    % max_indx based on max y
                    [~, max_indx1] = max(y_track_ahead);
                    
                    % max_indx based on first point that is within 1 cm
                    % distance from the landing disc
                    max_indx2 = find(y_track_ahead>y_hiveDisc-0.01,1);
                    
                    if ~isempty(max_indx2) && max_indx2 < max_indx1
                        max_indx = intersect_indices(ct) + max_indx2 - 1;
                    else
                        max_indx = intersect_indices(ct) + max_indx1 - 1;
                    end              
                    
                    rawState =  track(min_indx:max_indx,:);
                    indx_offground = find(rawState(:,6)>0.015, 1);
                    if ~isempty(indx_offground)
                        rawState = rawState(indx_offground:end,:);
                    end
                    
                    side = 'Hive';
                    
                    if size(rawState,1) > 30 && ...
                       all(IsInsideCylinder(center_hiveDisc', center_hiveDisc'-[0 2.5*y_margin 0], radius_hiveDisc, rawState(end-10:end,4:6))) && ... % mean(rawState(end-10:end,6)) >= 0.05 && mean(rawState(end-10:end,6)) <= 0.35 ...
                       sum(rawState(:,6)>0.05)/size(rawState,1) > 0.3
                        landingTracks_per_object.rawTrack(end+1) = rawState_BlindLandingtrack(rawState, side);
                    end
                    

                    
                elseif any(y_track_ahead < y_feederDisc+y_margin)
                    % track excerpt moving towards feeder
                    
                    % min_indx based on max y
                    [~, min_indx] = max(y_track_behind);
                    min_indx = intersect_indices(ct) - (length(y_track_behind) - min_indx);
                    
                    % max_indx based on min y
                    [~, max_indx1] = min(y_track_ahead);
                    
                    % max_indx based on first point that is within 1 cm
                    % distance from the landing disc
                    max_indx2 = find(y_track_ahead<y_feederDisc+0.01,1);
                    
                    if ~isempty(max_indx2) && max_indx2 < max_indx1
                        max_indx = intersect_indices(ct) + max_indx2 - 1;
                    else
                        max_indx = intersect_indices(ct) + max_indx1 - 1;
                    end  
                    
                    rawState =  track(min_indx:max_indx,:);
                    indx_offground = find(rawState(:,6)>0.015, 1);
                    if ~isempty(indx_offground)
                        rawState = rawState(indx_offground:end,:);
                    end
                    
                    side = 'Feeder';
                    
                    if size(rawState,1) > 30 && ...
                       all(IsInsideCylinder(center_feederDisc', center_feederDisc'+[0 2.5*y_margin 0], radius_feederDisc, rawState(end-10:end,4:6))) && ... % mean(rawState(end-10:end,6)) >= 0.05 && mean(rawState(end-10:end,6)) <= 0.35 ...
                       sum(rawState(:,6)>0.05)/size(rawState,1) > 0.3
                        landingTracks_per_object.rawTrack(end+1) = rawState_BlindLandingtrack(rawState, side);
                    end
                    
                end
            end
            
            % Remove duplicate track excerpts
            landingTracks_per_object.removeDuplicateTracks();
            
            
%             if isempty(landingTracks_per_object.rawTrack)
%                 landingTracks(ct_id) = BlindLandingtrack.empty(1,0);
%             end
            
            
%             for ct_landingDiscs = 1:length(landingDiscs)
%                 y_disc = disc_centers(2, ct_landingDiscs);
%                 if strcmpi(landingDiscs(ct_landingDiscs).side, 'Feeder')
%                     intersect_indices_feeder = find((track(1:end-1,5) < y_disc - y_margin & ...
%                         track(2:end,5) >= y_disc - y_margin) | (track(1:end-1,5) > y_disc - y_margin & ...
%                         track(2:end,5) <= y_disc - y_margin));
%                 else
%                     intersect_indices_hive = find(track(1:end-1,5) < y_disc + y_margin & ...
%                         track(2:end,5) >= y_disc + y_margin);
%                 end
%                 
%             end
            
            
            
        end
        
        % Select only non-empty ones
        treatment.landingTracks = [treatment.landingTracks landingTracks(arrayfun(@(x) ~isempty(x.rawTrack), landingTracks))];
        
        
        
    end
%     keyboard
    if rem(ct_treatment,8) == 0
        % save data file at the end of each day
        save(fullfile(rootDir, trackingDir, 'BlindLandingtracks.mat'), 'treatments');
%         keyboard
    end
    
    
end
save(fullfile(rootDir, trackingDir, 'BlindLandingtracks.mat'), 'treatments');
keyboard;

%% Display information about the # of tracks
clc; close all; 
clear; 

% Inputs
inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks.mat';
load(inputFile);

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
behaviour = {'rising','constant','sleeping'};

for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        for ct_behaviour = 1:length(behaviour)
            
            disp(['Pattern: ' pattern{ct_pattern} ...
                  ', light: ' light{ct_light} ...
                  ', behaviour: ' behaviour{ct_behaviour}]);
              
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
            landingTracks = [relevantTreatments.landingTracks];
            
            % # of landing tracks
            disp(['# of distinct Flydra objects (landingTracks): ' num2str(length(landingTracks))]);
  
        end
    end
end
keyboard;

%% Plot raw vs filtered data plots and save them on hard disk

clc; clear;
close all; 

DirPlots = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/plots/BlindLandingTracks/raw_vs_filtered';
savePlots = true;
delPreviousPlots = true;

% Inputs
inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks.mat';
load(inputFile);

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
behaviour = {'rising','constant','sleeping'};
% 
% pattern = {'spokes'};
% light = {'high'};

for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        for ct_behaviour = 1:length(behaviour)
        
            % Delete previous plots
            if delPreviousPlots && exist(fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour}]), 'dir')
                rmdir(fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour}]),'s');
            end
        
            % Create sub-directory for each treatment if it doesn't exist
            if savePlots && ~exist(fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour}]), 'dir')
                mkdir(fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour}]));
            end
            DirPlots_treatment = fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour}]);
            
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
            
            % Plot data 
            for ct_treatment=1:length(relevantTreatments)
                treatment = relevantTreatments(ct_treatment);
                if ~isempty(treatment.landingTracks)
                    disp(['Plotting, day: ' num2str(treatment.datenum) ...
                          ', pattern: ' treatment.pattern ...
                          ', light: ' treatment.light ...
                          ', behaviour: ' behaviour{ct_behaviour}]);
                    
                      for ct_track=1:length(treatment.landingTracks)
                          plotHandles = treatment.landingTracks(ct_track).plotData();
                          if savePlots
                              % Resizing the figures
                              for i=1:length(plotHandles)
                                  plotHandles(i).Position(3) = 680;
                                  plotHandles(i).Position(4) = 545;
                                  
                                  if i==1
                                      figureName = ['traj_' num2str(treatment.datenum) '_' ...
                                          num2str(treatment.startTime) '_' num2str(treatment.endTime) '_' ...
                                          num2str(treatment.landingTracks(ct_track).rawTrack(1).rawState(1,3)) ...
                                          '.png'];
                                  elseif i==2
                                      figureName = ['vel_' num2str(treatment.datenum) '_' ...
                                          num2str(treatment.startTime) '_' num2str(treatment.endTime) '_' ...
                                          num2str(treatment.landingTracks(ct_track).rawTrack(1).rawState(1,3)) ...
                                          '.png'];
                                  elseif i==3
                                      figureName = ['acc_' num2str(treatment.datenum) '_' ...
                                          num2str(treatment.startTime) '_' num2str(treatment.endTime) '_' ...
                                          num2str(treatment.landingTracks(ct_track).rawTrack(1).rawState(1,3)) ...
                                          '.png'];
                                  end
                                  
                                  saveas(plotHandles(i), fullfile(DirPlots_treatment, figureName) ,'png');
                              end
                          end
%                           keyboard;
                          close(plotHandles);
                      end
                    
                end
            end
            
            
        end
    end
end
keyboard;

%% Compute landing tracks in reference frame of the corresponding landing discs
% NOTE - Landing disc reference frame at Hive-disc center is left-handed and frame at 
% Feeder-disc center is right-handed
% Advantage - This makes comparison of landing parameters (such as Vy at LE) across discs much
% easier

clc; close all; clear;

% % % Inputs
inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks.mat';
load(inputFile);

for ct_treatment=1:length(treatments)
    treatment = treatments(ct_treatment);
    landingDiscs = treatment.landingDiscs;
    
    for ct_track=1:length(treatment.landingTracks)
        landingTrack = treatment.landingTracks(ct_track);
        
        % Computing states in landing disc reference frame
        landingTrack.compute_states_in_LDF(landingDiscs);
    end
end
save(inputFile, 'treatments');

% pattern = {'checkerboard', 'spokes'};
% light = {'low', 'medium', 'high'};
% 
% for ct_pattern = 1:length(pattern)
%     for ct_light = 1:length(light)
%         disp(['%%%%%%%%%% Treatment: ' pattern{ct_pattern} ' and ' light{ct_light} ' light condition %%%%%%%%%%%%']);
%         
%         relevantTreatments = treatments(strcmpi({treatments.pattern}, pattern{ct_pattern}) & ...
%                                         strcmpi({treatments.light}, light{ct_light}));
%                          
%         for ct_treatment = 1:length(relevantTreatments)
%             treatment = relevantTreatments(ct_treatment);
%             landingDiscs = treatment.landingDiscs;
%             landingTracks = treatment.landingTracks;
%             for ct_track=1:length(landingTracks)
%                 landingTrack = landingTracks(ct_track);
%                 
%                 % Computing states in landing disc reference frame
%                 landingTrack.compute_states_in_LDF(landingDiscs);
% %                 keyboard;
%             end
% %             keyboard;
%         end
% 
%     end
% end
% save(inputFile, 'treatments');

%% Plot optical expansion figures (V vs y ones) and save them on hard disk

clc; close all; 
clear; 

DirPlots = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/plots/BlindLandingTracks/oe_info_Vvsy';
savePlots = true;
delPreviousPlots = true;

% Inputs
inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks.mat';
load(inputFile);

% pattern = {'checkerboard', 'spokes'};
% light = {'low', 'medium', 'high'};
% behaviour = {'rising','constant','sleeping'};
% 
pattern = {'spokes'};
% light = {'high'};
behaviour = {'constant'};
light = {'low', 'medium', 'high'};


for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        for ct_behaviour = 1:length(behaviour)
        
            % Delete previous plots
            if delPreviousPlots && exist(fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour}]), 'dir')
                rmdir(fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour}]),'s');
            end
        
            % Create sub-directory for each treatment if it doesn't exist
            if savePlots && ~exist(fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour}]), 'dir')
                mkdir(fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour}]));
            end
            DirPlots_treatment = fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour}]);
            
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
            
            if isempty(relevantTreatments)
                error('relevantTreatments length is zero. Not possible!'); 
            end
            
            % Initializing figure on which data needs to be plotted
            opticalExpansionPlot_Time = figure;
            figure(opticalExpansionPlot_Time);
            subplotHandles_Time(1) = subplot(5,1,1); hold on;
            ylabel('V_{gy} (m/s)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
            
            subplotHandles_Time(2) = subplot(5,1,2); hold on;
            ylabel('y (m)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
            
            subplotHandles_Time(3) = subplot(5,1,3); hold on;
            ylabel('r (1/s)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
            ylim([-8 2]);
            
            subplotHandles_Time(4) = subplot(5,1,4); hold on;
            ylabel('r dot (1/s^2)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
            ylim([-50 30]);
            
            subplotHandles_Time(5) = subplot(5,1,5); hold on;
            ylabel('a_y (m/s^2)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
            xlabel('Time (s)', 'FontSize', 15);
            
            opticalExpansionPlot_Time.Position(3) = 1050;
            opticalExpansionPlot_Time.Position(4) = 800;
            opticalExpansionPlot_Time.Visible = 'off';
            
            
            % Initializing figure on which data needs to be plotted
            opticalExpansionPlot_Distance = figure;
            figure(opticalExpansionPlot_Distance);
            subplotHandles_Distance(1) = subplot(4,1,1); hold on;
            ylabel('V_{gy} (m/s)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
            
            subplotHandles_Distance(2) = subplot(4,1,2); hold on;
            ylabel('r (1/s)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
            ylim([-8 2]);
            
            subplotHandles_Distance(3) = subplot(4,1,3); hold on;
            ylabel('r dot (1/s^2)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
            ylim([-50 30]);
            
            subplotHandles_Distance(4) = subplot(4,1,4); hold on;
            ylabel('a_y (m/s^2)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
            xlabel('Distance from the platform (m)', 'FontSize', 15);
            
            opticalExpansionPlot_Distance.Position(3) = 1050;
            opticalExpansionPlot_Distance.Position(4) = 800;
            opticalExpansionPlot_Distance.Visible = 'off';
            
            % Plot data 
            for ct_treatment=1:length(relevantTreatments)
                treatment = relevantTreatments(ct_treatment);
                if ~isempty(treatment.landingTracks)
                    disp(['Plotting, day: ' num2str(treatment.datenum) ...
                          ', pattern: ' treatment.pattern ...
                          ', light: ' treatment.light ...
                          ', behaviour: ' behaviour{ct_behaviour}]);
                    
                      for ct_track=1:length(treatment.landingTracks) % for each landing track
                          for ct_excerpt=1:length(treatment.landingTracks(ct_track).state_LDF) % for each track excerpt
                              BlindLandingtrack.plotDataLDF_Time(treatment.landingTracks(ct_track).state_LDF(ct_excerpt), subplotHandles_Time);
                              BlindLandingtrack.plotDataLDF_Distance(treatment.landingTracks(ct_track).state_LDF(ct_excerpt), subplotHandles_Distance);
                              if savePlots   
                                  figureName = [treatment.landingTracks(ct_track).foldername(1:15) '_' ...
                                      num2str(treatment.landingTracks(ct_track).obj_id) '_' ...
                                      num2str(treatment.landingTracks(ct_track).state_LDF(ct_excerpt).filteredState(1,1)) ...
                                      '.png'];
                                      
                                  saveas(opticalExpansionPlot_Time, fullfile(DirPlots_treatment, ['time_' figureName]) ,'png');
                                  saveas(opticalExpansionPlot_Distance, fullfile(DirPlots_treatment, ['dist_' figureName]) ,'png');
                                  
                              end
%                               keyboard;
                              % clear all axes
                              arrayfun(@(x) cla(x), [subplotHandles_Time, subplotHandles_Distance]);
                          end
                      end
                    
                end
            end
            close all;
            
            
        end
    end
end
keyboard;

%% Plot optical expansion figures (a vs V ones) and save them on hard disk

clc; close all; 
clear; 

DirPlots = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/plots/BlindLandingTracks/oe_info_avsV';
savePlots = true;
delPreviousPlots = true;

% Inputs
inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks.mat';
load(inputFile);

% pattern = {'checkerboard', 'spokes'};
% light = {'low', 'medium', 'high'};
% behaviour = {'rising','constant','sleeping'};
% 
pattern = {'checkerboard'};
% light = {'high'};
behaviour = {'constant'};
light = {'low', 'medium', 'high'};


for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        for ct_behaviour = 1:length(behaviour)
        
            % Delete previous plots
            if delPreviousPlots && exist(fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour}]), 'dir')
                rmdir(fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour}]),'s');
            end
        
            % Create sub-directory for each treatment if it doesn't exist
            if savePlots && ~exist(fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour}]), 'dir')
                mkdir(fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour}]));
            end
            DirPlots_treatment = fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour}]);
            
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
            
            if isempty(relevantTreatments)
                error('relevantTreatments length is zero. Not possible!'); 
            end
            
            % Initializing figures on which data needs to be plotted
            % Figure 1 - Plotting with time
            opticalExpansionPlot_Time = figure;
            figure(opticalExpansionPlot_Time);
            subplotHandles_Time(1) = subplot(3,1,1); hold on;
            ylabel('V_{gy} / y (1/s)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
            ylim([-8 2]);
            
            subplotHandles_Time(2) = subplot(3,1,2); hold on;
            ylabel('a_y / V_{gy} (1/s)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
            ylim([-60 60]);
            
            subplotHandles_Time(3) = subplot(3,1,3); hold on;
            ylabel('a_y (m/s^2)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
            xlabel('Time (s)', 'FontSize', 15);
            
            opticalExpansionPlot_Time.Position(3) = 1050;
            opticalExpansionPlot_Time.Position(4) = 800;
            pause(2);
            opticalExpansionPlot_Time.Visible = 'off';
            
            
            % Figure 2 - Plotting with distance from the platform
            opticalExpansionPlot_Distance = figure;
            figure(opticalExpansionPlot_Distance);
            subplotHandles_Distance(1) = subplot(4,1,1); hold on;
            ylabel('V_{gy} (m/s)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
            
            subplotHandles_Distance(2) = subplot(4,1,2); hold on;
            ylabel('V_{gy} / y (1/s)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
            ylim([-8 2]);
            
            subplotHandles_Distance(3) = subplot(4,1,3); hold on;
            ylabel('a_y / V_{gy} (1/s)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
            ylim([-60 60]);
            
            subplotHandles_Distance(4) = subplot(4,1,4); hold on;
            ylabel('a_y (m/s^2)', 'FontSize', 15);
            set(gca, 'FontSize', 15); grid on;
            xlabel('y (m)', 'FontSize', 15);
%             ylim([-10 10]);
            
            opticalExpansionPlot_Distance.Position(3) = 1050;
            opticalExpansionPlot_Distance.Position(4) = 800;
            pause(2);
            opticalExpansionPlot_Distance.Visible = 'off';
            
%             % Figure 3 - Plotting with velocity
%             opticalExpansionPlot_Velocity = figure;
%             figure(opticalExpansionPlot_Velocity); hold on;
%             ylabel('a_y (m/s^2)', 'FontSize', 15);
%             xlabel('V_{gy} (m/s)', 'FontSize', 15);
%             set(gca, 'FontSize', 15); grid on;
%             opticalExpansionPlot_Velocity.Visible = 'off';
            
            pause(5);
            
            % Plot data 
            for ct_treatment=1:length(relevantTreatments)
                treatment = relevantTreatments(ct_treatment);
                if ~isempty(treatment.landingTracks)
                    disp(['Plotting, day: ' num2str(treatment.datenum) ...
                          ', pattern: ' treatment.pattern ...
                          ', light: ' treatment.light ...
                          ', behaviour: ' behaviour{ct_behaviour}]);
                    
                      for ct_track=1:length(treatment.landingTracks) % for each landing track
                          for ct_excerpt=1:length(treatment.landingTracks(ct_track).state_LDF) % for each track excerpt
                              BlindLandingtrack.plotDataLDF_Time(treatment.landingTracks(ct_track).state_LDF(ct_excerpt), subplotHandles_Time);
                              BlindLandingtrack.plotDataLDF_Distance(treatment.landingTracks(ct_track).state_LDF(ct_excerpt), subplotHandles_Distance);
%                               BlindLandingtrack.plotDataLDF_abyV(treatment.landingTracks(ct_track).state_LDF(ct_excerpt), opticalExpansionPlot_Velocity);
                              if savePlots   
                                  figureName = [treatment.landingTracks(ct_track).foldername(1:15) '_' ...
                                      num2str(treatment.landingTracks(ct_track).obj_id) '_' ...
                                      num2str(treatment.landingTracks(ct_track).state_LDF(ct_excerpt).filteredState(1,1)) ...
                                      '.png'];
                                      
                                  saveas(opticalExpansionPlot_Time, fullfile(DirPlots_treatment, ['with_time_' figureName]) ,'png');
                                  saveas(opticalExpansionPlot_Distance, fullfile(DirPlots_treatment, ['with_dist_' figureName]) ,'png');
%                                   saveas(opticalExpansionPlot_Velocity, fullfile(DirPlots_treatment, ['with_vel_' figureName]) ,'png');
                                  
                              end
%                               keyboard;
                              % clear all axes
                              arrayfun(@(x) cla(x), [subplotHandles_Time, subplotHandles_Distance]);
%                               arrayfun(@(x) cla(x), [subplotHandles_Time, subplotHandles_Distance, opticalExpansionPlot_Velocity.Children]);
                          end
                      end
                    
                end
            end
            close all;
            
            
        end
    end
end
keyboard;

%% function used in the script
function drawline(center, slope_deg, figure_handle)
    figure(figure_handle); hold on;
    viscircles(center, 1,'Color','red');
    
    
    x = [center(1)+20*cos(slope_deg), center(1)-20*cos(slope_deg)];
    y = [center(2)+20*sin(slope_deg), center(2)-20*sin(slope_deg)];
    line(x,y,'Color','green','LineWidth', 2);
    
    x = [center(1)+20*cos(-slope_deg), center(1)-20*cos(-slope_deg)];
    y = [center(2)+20*sin(-slope_deg), center(2)-20*sin(-slope_deg)];
    line(x,y,'Color','red','LineWidth', 2);
%     
%     x = [center(1)+20*cos(-slope_deg+pi/2), center(1)-20*cos(-slope_deg+pi/2)];
%     y = [center(2)+20*sin(-slope_deg+pi/2), center(2)-20*sin(-slope_deg+pi/2)];
%     line(x,y,'Color','white','LineWidth', 2);
%     
%     x = [center(1)+20*cos(slope_deg+pi/2), center(1)-20*cos(slope_deg+pi/2)];
%     y = [center(2)+20*sin(slope_deg+pi/2), center(2)-20*sin(slope_deg+pi/2)];
%     line(x,y,'Color','magenta','LineWidth', 2);
    

end
function in = IsInsideBox(data,box)
    % nans are counted outside box
    % data = [X Y Z]
    % box = [xmin xmax ymin ymax zmin zmax]
    in = [(data(:,1) >= box(1) & data(:,1) <= box(2)) & ...
          (data(:,2) >= box(3) & data(:,2) <= box(4)) & ...
          (data(:,3) >= box(5) & data(:,3) <= box(6))];        
end
function in = NanIsInsideBox(data,box)
    % Nans are counted inside box
    % data = [X Y Z]
    % box = [xmin xmax ymin ymax zmin zmax]
    in = [(data(:,1) >= box(1) & data(:,1) <= box(2)) & ...
          (data(:,2) >= box(3) & data(:,2) <= box(4)) & ...
          (data(:,3) >= box(5) & data(:,3) <= box(6))]; 
    in(isnan(data(:,1)) | isnan(data(:,2)) | isnan(data(:,1))) = true;
end
function in = IsInsideCylinder(p1, p2, radius, q)
    % p1 - 1X3
    % p2 - 1X3
    % r - 1X1
    % q - query points (NX3)
    vec = p2-p1;
    const = radius*norm(vec);
    
    in = dot(q-p1, repmat(vec,size(q,1),1), 2) >= 0 & ...
         dot(q-p2, repmat(vec,size(q,1),1), 2) <= 0 & ...
         sum((cross(q-p1, repmat(vec, size(q,1), 1))).^2,2).^0.5 <= const;

end

% Functions used in this script file

function [folders] = GetFolders(direc,pattern)
    % Get a list of all folders in a directory containing patterns in their names.
    
    % direc - path where to look
    % pattern - cell array containing string patterns that must be in the
    % names of the folders

    all    = dir(direc);
    folders = all([all(:).isdir]);
    names    = {folders.name};

    % Get a logical vector that discards folders that aren't required
    flags = ~strcmp(names, '.') & ...
              ~strcmp(names, '..') & ~strcmp(cellfun(@(x){x(1)},names),'~');%~cellfun('isempty',strfind(A,B));
    for ct=1:length(pattern)
       flags = flags & contains(names,pattern{ct});
    end

    % Extract only those that contain pattern
    folders = names(flags);
end