%% Landing tracks analysis

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

%% Compute parameters representing landing performance
clc; close all; clear;

% Inputs
inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/videoscoring_and_tracks.mat';
load(inputFile);

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};

for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        disp(['%%%%%%%%%% Treatment: ' pattern{ct_pattern} ' and ' light{ct_light} ' light condition %%%%%%%%%%%%']);
        
        relevantTreatments = treatments([treatments.videosParsed] & ...
                             strcmpi({treatments.pattern}, pattern{ct_pattern}) & ...
                             strcmpi({treatments.light}, light{ct_light}));
        
        
        relevantVideos = [relevantTreatments.usefulVideoFiles];
        
        for ct_video=1:length(relevantVideos)
            videoFile = relevantVideos(ct_video);
            track = videoFile.landingTrack;
            
            if ~isempty(track.state) && track.confirmedChecks_le
                % Handling leg extension
                if ~isnan(videoFile.scoredData.landingParams.legextension_framenumber)
                    time_le = videoFile.scoredData.timestamps(...
                                    videoFile.scoredData.landingParams.legextension_framenumber);
                    track.time = track.state(:,1) - time_le;

                    if max(track.state(:,1)) >= time_le && min(track.state(:,1)) <= time_le
                        track.state_le = [time_le interp1(track.state(:,1)-time_le, track.state(:,2:10), 0, 'makima')];
    %                     track.state_le = track.state(abs(track.time)<1e-5, :);

%                         track.displacement = (sum((track.state(:,2:4)-track.state_le(:,2:4)).^2, 2)).^0.5;
                    else
                        disp('Not computing LE state');
                    end
                end
                
                % Handling touchdown 1
                if ~isnan(videoFile.scoredData.landingParams.touchdown_framenumber)
                    time_t1 = videoFile.scoredData.timestamps(...
                                    videoFile.scoredData.landingParams.touchdown_framenumber);

                    assert(time_t1 >= time_le, 'Tracks are not sorted with time rows!');
                    
%                     track.time = track.state(:,1) - time_t1;
                    
                    track.duration_le_t1 = time_t1 - time_le;

                    if max(track.state(:,1)) > time_t1 && min(track.state(:,1)) <= time_t1
                        track.state_t1 = [time_t1 interp1(track.state(:,1)-time_t1, track.state(:,2:10), 0, 'makima')];
    %                     track.state_t1 = track.state(abs(track.state(:,1)-time_t1)<1e-5, :);
                        
                        track.displacement = (sum((track.state(:,2:4)-track.state_t1(:,2:4)).^2, 2)).^0.5;
        
                        pos_with_endpts = [track.state_le(:, 2:4); ...
                                           track.state(track.state(:,1)-track.state_le(1,1) > 1e-5 & ...
                                                      track.state(:,1)-track.state_t1(1,1) < -1e-5, 2:4); ...
                                           track.state_t1(:,2:4)];

                        track.dist_le_t1 = sum((sum((diff(pos_with_endpts)).^2, 2)).^0.5);

                        track.disp_le_t1 = (sum((track.state_t1(:,2:4)-track.state_le(:, 2:4)).^2, 2)).^0.5;

                        track.turtuosity_le_t1 = track.dist_le_t1/track.disp_le_t1;
                    else
                        disp('Not computing T1 state');
                    end
                end
                
                % Handling touchdown 2
                if ~isnan(videoFile.scoredData.landingParams.touchdown2_framenumber)
                    time_t2 = videoFile.scoredData.timestamps(...
                                    videoFile.scoredData.landingParams.touchdown2_framenumber);

                    assert(time_t2 >= time_t1, 'Tracks are not sorted with time rows!');

                    track.duration_t1_t2 = time_t2 - time_t1;

                    track.duration_le_t2 = time_t2 - time_le;

                    if max(track.state(:,1)) >= time_t2 && min(track.state(:,1)) <= time_t2
                        track.state_t2 = [time_t2 interp1(track.state(:,1)-time_t2, track.state(:,2:10), 0, 'makima')];
    %                     track.state_t2 = track.state(abs(track.state(:,1)-time_t2)<1e-5, :);
                    else
                        disp('Not computing T2 state');
                    end

                end
                
                
            end
        end

    end
end
save(inputFile, 'treatments');

%% Compute landing tracks in reference frame of the corresponding landing discs
% NOTE - Landing disc reference frame at Hive-disc center is left-handed and frame at 
% Feeder-disc center is right-handed
% Advantage - This makes comparison of landing parameters (such as Vy at LE) across discs much
% easier

% This can only be done once state_le, t1, t2 are computed using above
% section

clc; close all; clear;

% Inputs
inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/videoscoring_and_tracks.mat';
load(inputFile);

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};

for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        disp(['%%%%%%%%%% Treatment: ' pattern{ct_pattern} ' and ' light{ct_light} ' light condition %%%%%%%%%%%%']);
        
        relevantTreatments = treatments([treatments.videosParsed] & ...
                             strcmpi({treatments.pattern}, pattern{ct_pattern}) & ...
                             strcmpi({treatments.light}, light{ct_light}));
                         
        for ct_treatment = 1:length(relevantTreatments)
            treatment = relevantTreatments(ct_treatment);
            landingDiscs = treatment.landingDiscs;
            relevantVideos = treatment.usefulVideoFiles;
            for ct_video=1:length(relevantVideos)
                videoFile = relevantVideos(ct_video);
                track = videoFile.landingTrack;
                
                if ~isempty(track.state) && track.confirmedChecks_le
                    % Computing states in landing disc reference frame
                    track.compute_states_in_LDF(landingDiscs);
                end
                
            end
        end

    end
end
save(inputFile, 'treatments');


%% Checks - Comparison of number of tracks computed with number of videos scored
clc; 
pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};

for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        disp(['Treatment: ' pattern{ct_pattern} ' and ' light{ct_light} ' light condition']);
        
        relevantTreatments = treatments([treatments.videosParsed] & ...
                             strcmpi({treatments.pattern}, pattern{ct_pattern}) & ...
                             strcmpi({treatments.light}, light{ct_light}));
        
        disp(['# of repititions: ' num2str(length(relevantTreatments))]);                 
        
        relevantVideos = {relevantTreatments.usefulVideoFiles};
        arraysize = cellfun(@(x) length(x), relevantVideos, 'UniformOutput', false);
        disp(['# of useful videos: ' num2str(length([relevantVideos{:}])) ' (' num2str([arraysize{:}]) ')']);
        

        relevantTracks = cellfun(@(x) [x.landingTrack], relevantVideos, 'UniformOutput', false);
        relevantTracks = cellfun(@(x) x(...
            arrayfun(@(y) ~isempty(y.state), x)), relevantTracks, 'UniformOutput', false);

        arraysize = cellfun(@(x) length(x), relevantTracks, 'UniformOutput', false);
        disp(['# of landing tracks for which state is computed: ' num2str(length([relevantTracks{:}])) ' (' num2str([arraysize{:}]) ')']);
        
        
        relevantTracks = cellfun(@(x) [x.landingTrack], relevantVideos, 'UniformOutput', false);
        relevantTracks = cellfun(@(x) x(...
            arrayfun(@(y) (~isempty(y.state) & y.confirmedChecks_le), x)), relevantTracks, 'UniformOutput', false);

        arraysize = cellfun(@(x) length(x), relevantTracks, 'UniformOutput', false);
        disp(['# of landing tracks for which state is computed and both checks are confirmed: ' num2str(length([relevantTracks{:}])) ' (' num2str([arraysize{:}]) ')']);
        
        
        relevantVideos = [relevantTreatments.usefulVideoFiles];
        relevantTracks = [relevantVideos.landingTrack];
        
        relevantTracks = relevantTracks(arrayfun(@(x) ~isempty(x.state), relevantTracks));
        relevantVideos = relevantVideos(arrayfun(@(x) ~isempty(x.state), relevantTracks));
        
        disp(['# of tracks containing non-empty LE state: ' num2str(sum(arrayfun(@(x) ~isempty(x.state_le), relevantTracks)))]);
        disp(['# of tracks containing non-empty T1 state: ' num2str(sum(arrayfun(@(x) ~isempty(x.state_t1), relevantTracks)))]);
        disp(['# of tracks containing non-empty T2 state: ' num2str(sum(arrayfun(@(x) ~isempty(x.state_t2), relevantTracks)))]);
%         sum(arrayfun(@(x) x.state_le(1) <= x.state(end,1), relevantTracks))

%         for ct_video=1:length(relevantVideos)
%             plot3Dand2Ddata(relevantVideos(ct_video));
%         end
        
    end
end


%% Plot landing performance parameters
% Plots are per pattern i.e., to see difference among different light
% conditions for each pattern

clc;
close all; clear;

inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/videoscoring_and_tracks.mat';
load(inputFile);


pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};

labels = cell(0, 1); % for x-axis of boxplots
pattern_label = {'+', 'x'}; % + is checkerboard, x is spokes
light_label = {'L', 'M', 'H'};

duration_le_t1 = cell(0, 1);
duration_t1_t2 = cell(0, 1);
duration_le_t2 = cell(0, 1);

dist_le_t1 = cell(0, 1);
disp_le_t1 = cell(0, 1);
turtuosity_le_t1 = cell(0, 1);

speed_le = cell(0,1);
speed_t1 = cell(0,1);
speed_t2 = cell(0,1);

y_platform_le = cell(0,1); % y distance from platform at LE
y_platform_t1 = cell(0,1); % y distance from platform at T1
y_platform_t2 = cell(0,1); % y distance from platform at T2

Vgy_le = cell(0,1); % Ground velocity in y-direction at LE in LDF (landing disc reference frame)
Vgy_t1 = cell(0,1);
Vgy_t2 = cell(0,1);

for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        relevantTreatments = treatments([treatments.videosParsed] & ...
                             strcmpi({treatments.pattern}, pattern{ct_pattern}) & ...
                             strcmpi({treatments.light}, light{ct_light}));
        
        
        relevantVideos = [relevantTreatments.usefulVideoFiles];
        
        relevantTracks = [relevantVideos.landingTrack];
        
        relevantTracks = relevantTracks(arrayfun(@(x) ~isempty(x.state) & x.confirmedChecks_le, relevantTracks));
        
        labels{end+1, 1} = [pattern_label{ct_pattern} ' ' light_label{ct_light}];
        
        duration_le_t1{end+1, 1} = [relevantTracks.duration_le_t1];
        duration_t1_t2{end+1, 1} = [relevantTracks.duration_t1_t2];
        duration_le_t2{end+1, 1} = [relevantTracks.duration_le_t2];
        
        dist_le_t1{end+1, 1} = [relevantTracks.dist_le_t1];
        disp_le_t1{end+1, 1} = [relevantTracks.disp_le_t1];
        turtuosity_le_t1{end+1, 1} = [relevantTracks.turtuosity_le_t1];
        
        state_le = vertcat(relevantTracks.state_le);
        speed_le{end+1, 1} = (sum(state_le(:,5:7).^2, 2)).^0.5;
        
        state_t1 = vertcat(relevantTracks.state_t1);
        speed_t1{end+1, 1} = (sum(state_t1(:,5:7).^2, 2)).^0.5; 
        
        state_t2 = vertcat(relevantTracks.state_t2);
        speed_t2{end+1, 1} = (sum(state_t2(:,5:7).^2, 2)).^0.5; 
        
        state_le_LDF = vertcat(relevantTracks.state_le_LDF);
        y_platform_le{end+1, 1} = state_le_LDF(:,3);
        Vgy_le{end+1, 1} = state_le_LDF(:,6);
        
        state_t1_LDF = vertcat(relevantTracks.state_t1_LDF);
        y_platform_t1{end+1, 1} = state_t1_LDF(:,3);
        Vgy_t1{end+1, 1} = state_t1_LDF(:,6);
        
        state_t2_LDF = vertcat(relevantTracks.state_t2_LDF);
        y_platform_t2{end+1, 1} = state_t2_LDF(:,3);
        Vgy_t2{end+1, 1} = state_t2_LDF(:,6);

    end
end

approach_duration_figHandle = createBoxPlot(duration_le_t1, labels, 'Approach duration (s)');
landing_duration_figHandle = createBoxPlot(duration_t1_t2, labels, 'Landing duration (s)');
total_landing_duration_figHandle = createBoxPlot(duration_le_t2, labels, 'Total landing duration (s)');

dist_le_t1_figHandle = createBoxPlot(dist_le_t1, labels, 'Dist LE and T1 (m)');
disp_le_t1_figHandle = createBoxPlot(disp_le_t1, labels, 'Disp LE and T1 (m)');
turtuosity_le_t1_figHandle = createBoxPlot(turtuosity_le_t1, labels, 'Turtuosity LE and T1 (m)');

speed_le_figHandle = createBoxPlot(speed_le, labels, 'Speed LE (m/s)');
speed_t1_figHandle = createBoxPlot(speed_t1, labels, 'Speed T1 (m/s)');
speed_t2_figHandle = createBoxPlot(speed_t2, labels, 'Speed T2 (m/s)');

y_platform_le_figHandle = createBoxPlot(y_platform_le, labels, 'y from platform at LE (m)');
y_platform_t1_figHandle = createBoxPlot(y_platform_t1, labels, 'y from platform at T1 (m)');
y_platform_t2_figHandle = createBoxPlot(y_platform_t2, labels, 'y from platform at T2 (m)');

Vgy_le_figHandle = createBoxPlot(Vgy_le, labels, 'V_{gy} at LE (m)');
Vgy_t1_figHandle = createBoxPlot(Vgy_t1, labels, 'V_{gy} at T1 (m)');
Vgy_t2_figHandle = createBoxPlot(Vgy_t2, labels, 'V_{gy} at T2 (m)');

%% Plot trajectory components with time and distance from platform for each track and save them on hard disk
% Plots are produced in the landing discs reference frame
DirPlots = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/plots/trajectories';
savePlots = true;

% clc; clear;
close all; 
% Inputs
% inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/videoscoring_and_tracks.mat';
% load(inputFile);

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
% 
% pattern = {'spokes'};
% light = {'high'};

for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        
        % Create sub-directory for each treatment if it doesn't exist
        if savePlots && ~exist(fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light}]), 'dir')
            mkdir(fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light}]));
        end
        DirPlots_treatment = fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light}]);
        
        relevantTreatments = treatments([treatments.videosParsed] & ...
                             strcmpi({treatments.pattern}, pattern{ct_pattern}) & ...
                             strcmpi({treatments.light}, light{ct_light}));
        
        
        relevantVideos = [relevantTreatments.usefulVideoFiles];
        
        for ct_video=1:length(relevantVideos)
            videoFile = relevantVideos(ct_video);
            track = videoFile.landingTrack;
            
            if ~isempty(track.state) && track.confirmedChecks_le
                disp(['Analysing, day: '  videoFile.scoredData.year ...
                    videoFile.scoredData.month  videoFile.scoredData.day ...
                    ', video: ' videoFile.scoredData.name]);
                
                % Plotting the data
                plotHandles = track.plotData_LDF();
                
                if savePlots
                    % Resizing the figures
                    for i=1:length(plotHandles)
                        plotHandles(i).Position(3) = 680;
                        plotHandles(i).Position(4) = 545;
                        
                        if i==1
                            figureName = ['traj_' videoFile.datenum '_' ...
                                videoFile.scoredData.name(1:6) '_' ...
                                videoFile.scoredData.landingParams.landing_side '_' ...
                                videoFile.scoredData.landingParams.landing_location '.png'];
                        elseif i==2
                            figureName = ['vel_' videoFile.datenum '_' ...
                                videoFile.scoredData.name(1:6) '_' ...
                                videoFile.scoredData.landingParams.landing_side '_' ...
                                videoFile.scoredData.landingParams.landing_location '.png'];
                        elseif i==3
                            figureName = ['acc_' videoFile.datenum '_' ...
                                videoFile.scoredData.name(1:6) '_' ...
                                videoFile.scoredData.landingParams.landing_side '_' ...
                                videoFile.scoredData.landingParams.landing_location '.png'];
                        elseif i==4
                            figureName = ['oe_' videoFile.datenum '_' ...
                                videoFile.scoredData.name(1:6) '_' ...
                                videoFile.scoredData.landingParams.landing_side '_' ...
                                videoFile.scoredData.landingParams.landing_location '.png'];
                        end
                        
                        saveas(plotHandles(i), fullfile(DirPlots_treatment, figureName) ,'png');
                    end
                end
                    
                
                
%                 keyboard;
                close(plotHandles);
            end
        end        
        
    end
end


%% Extract subset of landing trajectories that are useful to look at from the dynamics perspective
clc; close all; clear;

% Inputs
inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/videoscoring_and_tracks.mat';
load(inputFile);

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
% 
% pattern = {'spokes'};
% light = {'high'};


% selectedVideoFiles = VideoMetadata.empty;

for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        disp(['%%%%%%%%%% Treatment: ' pattern{ct_pattern} ' and ' light{ct_light} ' light condition %%%%%%%%%%%%']);
        
        relevantTreatments = treatments([treatments.videosParsed] & ...
                             strcmpi({treatments.pattern}, pattern{ct_pattern}) & ...
                             strcmpi({treatments.light}, light{ct_light}));
                         
        for ct_treatment = 1:length(relevantTreatments)
            treatment = relevantTreatments(ct_treatment);
            relevantVideos = treatment.usefulVideoFiles;
            for ct_video=1:length(relevantVideos)
                videoFile = relevantVideos(ct_video);
                track = videoFile.landingTrack;
                
                if ~isempty(track.state) && track.confirmedChecks_le
                    state_subset = track.state_LDF(track.time>=-1.5,:);
                    if any(state_subset(:,3)>0.15)
                        
                        % copy this track and videoFile information
%                         selectedVideoFiles(end+1) = videoFile;
                        
                        % extract y (less than 0.15m) when |Vy| < 0.025 m/s
                        % for the first time while a bbee approaches the
                        % surface
                        state_subset2 = state_subset(state_subset(:,3) <= 0.15,:);
                        
                        equinox_state_indx = find(state_subset2(:,6)>-0.025,1,'first');
                        if isempty(equinox_state_indx)
                            track.equinox_state_LDF = [];
                            track.stable_state_LDF = [];
                        else
                            track.equinox_state_LDF = state_subset2(equinox_state_indx,:);

                            equinox_indx_stateSubset = find(abs(state_subset(:,1)-track.equinox_state_LDF(1,1))<1e-5);
                            track.stable_state_LDF = state_subset(1:equinox_indx_stateSubset, :);
                        end
                        
                        % TO DEBUG
                        %state_subset2(abs(state_subset2(:,6))<0.025,[3 6])
                        
                        
                        
                        
                    end
%                     track.state_LDF()
                end
                
            end
        end
    end
end
save(inputFile, 'treatments');

%% Plot the dynamics perspective variables: r_ref and equinox point 
close all;
% clc; clear;

% Inputs
% inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/videoscoring_and_tracks.mat';
% load(inputFile);

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};

labels = cell(0, 1); % for x-axis of boxplots
pattern_label = {'+', 'x'}; % + is checkerboard, x is spokes
light_label = {'L', 'M', 'H'};

distance_equinox = cell(0, 1);

% 
% pattern = {'spokes'};
% light = {'high'};

    
% selectedVideoFiles = VideoMetadata.empty;
colormap = {brewermap(11,'PiYG'), ...
            brewermap(11,'RdBu')};
        
for ct_pattern = 1:length(pattern)
    
    colors = colormap{ct_pattern}(9:11,:);
    DynamicVarPlotHandle(ct_pattern) = figure; % to plot displacement with time
    subplot(3,1,1); hold on; % to plot low light condition
    subplot(3,1,2); hold on; % to plot medium light condition
    subplot(3,1,3); hold on; % to plot high light condition
    
    for ct_light = 1:length(light)
        disp(['%%%%%%%%%% Treatment: ' pattern{ct_pattern} ' and ' light{ct_light} ' light condition %%%%%%%%%%%%']);
        
        relevantTreatments = treatments([treatments.videosParsed] & ...
                             strcmpi({treatments.pattern}, pattern{ct_pattern}) & ...
                             strcmpi({treatments.light}, light{ct_light}));
                         
        relevantVideos = [relevantTreatments.usefulVideoFiles];
        
        relevantTracks = [relevantVideos.landingTrack];
        
        isUsefulTrack = false(length(relevantTracks),1);
        for ct=1:length(relevantTracks)
            x = relevantTracks(ct);
            if ~isempty(x.stable_state_LDF) && x.confirmedChecks_le && ...
               ~isempty(x.scoredLandingParams) && isequal(x.scoredLandingParams.landing_success,2) % ...
%                && strcmpi(x.scoredLandingParams.landing_location,'Tube')
                isUsefulTrack(ct) = true;
            end           
        end
        
        relevantTracks = relevantTracks(isUsefulTrack);
        
        labels{end+1, 1} = [pattern_label{ct_pattern} ' ' light_label{ct_light}];
        
        distance_equinox{end+1, 1} = arrayfun(@(x) x.equinox_state_LDF(:,3), relevantTracks);
        
        for ct_track=1:5%length(relevantTracks)
            track = relevantTracks(ct_track);
            
            figure(DynamicVarPlotHandle(ct_pattern));
            subplot(3,1,ct_light);
            
            if ct_track == 1
                plot(track.stable_state_LDF(:,3),track.stable_state_LDF(:,6)./track.stable_state_LDF(:,3) ...
                    ,'Color',colors(ct_light,:),'LineWidth',2,'DisplayName',[pattern_label{ct_pattern} ' ' light_label{ct_light}]);
            else
                plot(track.stable_state_LDF(:,3),track.stable_state_LDF(:,6)./track.stable_state_LDF(:,3) ...
                    ,'Color',colors(ct_light,:),'LineWidth',2,'HandleVisibility','off');
            end
%             keyboard;
        end
        
        
    end
end

y_equinox_figHandle = createBoxPlot(distance_equinox, labels, 'Distance from the platform at |V_{gy}| < 0.025 m/s for the 1st time');
arrayfun(@(x) setOpticalFlowFigure(x),DynamicVarPlotHandle);

%% Plot y_equinox with time within the treatment (distinguishing between treatment before)
close all;
clc; clear;

% Inputs
inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/videoscoring_and_tracks.mat';
load(inputFile);

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};

pattern_label = {'+', 'x'}; % + is checkerboard, x is spokes
light_label = {'L', 'M', 'H'};


% distance_equinox_time 
% rows are the actual treatment & columns are the prior treatment
% diagonal entries will always be empty
% 1 column & 1 row - light{1},  2 column & 2 row - light{2}, 3 column & 3 row - light{3}
% Each entry is a cell containing N arrays (time in the treatment vs...
% y_equinox), where N is the number of repititions of the...
% treatment+prior_treatment combination
distance_equinox_time = cell(length(light), length(light)); 

dynamicLearningPlots = gobjects(length(pattern), length(light)); % matrix containing figure handles
for ct_pattern = 2%1:length(pattern)
    for ct_lightTreatment=1:length(light)
        
        
        dynamicLearningPlots(ct_pattern, ct_lightTreatment) = figure; % Plotting for a particular pattern and a light condition
        hold on;
        for ct_priorTreatment=1:length(light)
            if ct_lightTreatment == ct_priorTreatment
%                 distance_equinox_time(ct_lightTreatment, ct_priorTreatment) = {};
            else
                % find y_equinox for this treatment+prior_treatment combination
                % and save it in distance_equinox_time. This is done for all
                % repititions of treatment+prior_treatment combination.
                lightTreatment = light{ct_lightTreatment};
                priorTreatment = light{ct_priorTreatment};
                patternTreatment = pattern{ct_pattern};
                
                % % Find treatment with priorTreatment
                % first find indices of treatment
                indices_treatment = find(strcmpi({treatments.pattern}, patternTreatment) & ...
                             strcmpi({treatments.light}, lightTreatment));
                indices_treatment = indices_treatment(indices_treatment > 1); % no prior treatment for treatments(1)
                                
                % find treatments with videos parsed and having priorTreatment
                treatments_subset = Lighttreatment.empty;
                for ct=1:length(indices_treatment)
                    treatment = treatments(indices_treatment(ct));
                    if treatment.videosParsed && treatment.startTime >= 093000 && ...
                            strcmpi(priorTreatment, treatments(indices_treatment(ct)-1).light) && ...
                            treatment.datenum ==  treatments(indices_treatment(ct)-1).datenum
                        treatments_subset(end+1) = treatment;                        
                    end
                end
                
                % treatments_subset contains repititions of the ...
                % treatment+priorTreatment combination
                % Compute time withtin the treatment and corresponding
                % y_equinox for each treatment+priorTreatment repitition 
                
                for ct_subset=1:length(treatments_subset)
                    treatment = treatments_subset(ct_subset);
                    relevantVideos = [treatment.usefulVideoFiles];
                    relevantTracks = [relevantVideos.landingTrack];
        
                    isUsefulTrack = false(length(relevantTracks),1);
                    for ct=1:length(relevantTracks)
                        x = relevantTracks(ct);
                        if ~isempty(x.stable_state_LDF) && x.confirmedChecks_le && ...
                           ~isempty(x.scoredLandingParams) && isequal(x.scoredLandingParams.landing_success,2) % ...
            %                && strcmpi(x.scoredLandingParams.landing_location,'Tube')
                            isUsefulTrack(ct) = true;
                        end           
                    end
                    relevantTracks = relevantTracks(isUsefulTrack);
                    if ~isempty(relevantTracks)
                        equinox_state = vertcat(relevantTracks.equinox_state_LDF);
                        datetimeSinceEpoch = datetime(equinox_state(:,1), 'convertfrom', 'posixtime', 'TimeZone', 'Europe/Zurich');
                        treatment_startDatetime = datetime([num2str(treatment.datenum) ' ' num2str(treatment.startTime, '%06.f')], ...
                            'InputFormat', 'yyyyMMdd HHmmss', 'TimeZone', 'Europe/Zurich');
                        time_hours = seconds(datetimeSinceEpoch-treatment_startDatetime)/3600;
                        distance_equinox_time{ct_lightTreatment, ct_priorTreatment}{end+1} = [time_hours equinox_state(:,3)];
                    end
                end
                
                if ~isempty(relevantTracks)
                    figure(dynamicLearningPlots(ct_pattern, ct_lightTreatment));
                    if strcmpi(priorTreatment, 'low')
                        colors = brewermap(9, '*BuGn');
                        markerShape = 'o';
                    elseif strcmpi(priorTreatment, 'medium')
                        colors = brewermap(9, '*BuPu');
                        markerShape = 's';
                    elseif strcmpi(priorTreatment, 'high')
                        colors = brewermap(9, '*GnBu');
                        markerShape = '^';
                    end
                    for ct=1:length(distance_equinox_time{ct_lightTreatment, ct_priorTreatment})
                        if ct==1
                            plot(distance_equinox_time{ct_lightTreatment, ct_priorTreatment}{ct}(:,1),distance_equinox_time{ct_lightTreatment, ct_priorTreatment}{ct}(:,2), markerShape, ...
                                'MarkerSize', 10, 'MarkerEdgeColor',colors(ct,:),'MarkerFaceColor',colors(ct,:),'LineWidth',2,'DisplayName',[pattern_label{ct_pattern} ' ' light_label{ct_priorTreatment} ' ' char(8594) ' ' pattern_label{ct_pattern} ' ' light_label{ct_lightTreatment}]);
                        else
                            plot(distance_equinox_time{ct_lightTreatment, ct_priorTreatment}{ct}(:,1),distance_equinox_time{ct_lightTreatment, ct_priorTreatment}{ct}(:,2), markerShape, ...
                                'MarkerSize', 10, 'MarkerEdgeColor',colors(ct,:),'MarkerFaceColor',colors(ct,:),'LineWidth',2,'HandleVisibility','off');
                        end
                    end
                end
                
%                 keyboard;
                
            end

        end % Data is accrued for all possible prior treatments combination
        
        ylabel('y_{equinox}', 'FontSize', 16);
        xlabel('time within the treatment (hours)', 'FontSize', 16);
        set(gca, 'FontSize', 18); grid on;
        legend(gca,'show','Location','best');
        xlim([0 1.5]);
        ylim([0 0.15]);
        
    end
end


% arrayfun(@(x) ~isempty(x.stable_state_LDF) && x.confirmedChecks_le && ...
%                            ~isempty(x.scoredLandingParams) && isequal(x.scoredLandingParams.landing_success,2), y)
%% Plot trajectory parameters for each light condition
% clc; close all; clear;

% Inputs
% inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/videoscoring_and_tracks.mat';
% load(inputFile);

% pattern = {'checkerboard', 'spokes'};
% light = {'low', 'medium', 'high'};

pattern = {'spokes'};
light = {'low'};

colormap = {brewermap(11,'PiYG'), ...
            brewermap(11,'RdBu')};


for ct_pattern = 1:length(pattern)
    
    colors = colormap{ct_pattern}(9:11,:);
    
    trajPlotHandle = figure; % to plot displacement with time
    subplot(3,1,1); 
%     hold on;
    subplot(3,1,2); 
%     hold on;
    subplot(3,1,3); 
%     hold on;

    velPlotHandle1 = figure; % to plot speed with time
    subplot(3,1,1); 
    hold on;
    subplot(3,1,2); 
    hold on;
    subplot(3,1,3); 
    hold on;

    velPlotHandle2 = figure; % to plot speed against Eucledian displacement measured from landing point
    subplot(3,1,1); 
%     hold on;
    subplot(3,1,2); 
%     hold on;
    subplot(3,1,3); 
%     hold on;

    velPlotHandle3 = figure; % to plot speed yvel against y
%     subplot(3,1,1); 
%     hold on;
%     subplot(3,1,2); 
%     hold on;
%     subplot(3,1,3); 
%     hold on;
    
    for ct_light = 1:length(light)
        disp(['%%%%%%%%%% Treatment: ' pattern{ct_pattern} ' and ' light{ct_light} ' light condition %%%%%%%%%%%%']);
        
        relevantTreatments = treatments([treatments.videosParsed] & ...
                             strcmpi({treatments.pattern}, pattern{ct_pattern}) & ...
                             strcmpi({treatments.light}, light{ct_light}));
        
        
        relevantVideos = [relevantTreatments.usefulVideoFiles];
        
        for ct_video=1:length(relevantVideos) % 17 for spokes, high
            videoFile = relevantVideos(ct_video);
            track = videoFile.landingTrack;
            
            if ~isempty(track.state_t1) && videoFile.scoredData.landingParams.landing_success == 2
                figure(trajPlotHandle);
                subplot(3,1,ct_light);
                if ct_video == 1
                    plot(track.time,track.displacement,'Color',colors(ct_light,:),'LineWidth',2,'DisplayName',['Treatment: ' pattern{ct_pattern} ' and ' light{ct_light}]);
                else
                    plot(track.time,track.displacement,'Color',colors(ct_light,:),'LineWidth',2,'HandleVisibility','off');
                end

                figure(velPlotHandle1);
                subplot(3,1,ct_light);
                if ct_video == 1
                    plot(track.time,sum((track.state(:, 5:7)).^2,2),'Color',colors(ct_light,:),'LineWidth',2,'DisplayName',['Treatment: ' pattern{ct_pattern} ' and ' light{ct_light}]);
                else
                    plot(track.time,sum((track.state(:, 5:7)).^2,2),'Color',colors(ct_light,:),'LineWidth',2,'HandleVisibility','off');
                end

                displacement = track.displacement;
                displacement(track.state(:,1)-track.state_t1(1)>1e-3) = nan;
                displacement(displacement > 0.18) = nan;
                figure(velPlotHandle2);
                subplot(3,1,ct_light);
                if ct_video == 1
                    plot(displacement, sum((track.state(:, 5:7)).^2,2),'Color',colors(ct_light,:),'LineWidth',2,'DisplayName',['Treatment: ' pattern{ct_pattern} ' and ' light{ct_light}]);
                else
                    plot(displacement, sum((track.state(:, 5:7)).^2,2),'Color',colors(ct_light,:),'LineWidth',2,'HandleVisibility','off');
                end
                xlim([0 0.18]); ylim([0 1]);
                
                figure(velPlotHandle3);
%                 subplot(3,1,ct_light);
                if ct_video == 1
                    plot(track.state(:,3), track.state(:,6), 'Color',colors(ct_light,:),'LineWidth',2,'DisplayName',['Treatment: ' pattern{ct_pattern} ' and ' light{ct_light}]);
                else
                    plot(track.state(:,3), track.state(:,6), 'Color',colors(ct_light,:),'LineWidth',2,'HandleVisibility','off');
                end
%                 xlim([0 0.18]); ylim([0 1]);
            end
            keyboard
            
        end
    end
    figure(trajPlotHandle);
    subplot(3,1,1);
    ylabel('Disp (m)', 'FontSize', 16);
    set(gca, 'FontSize', 18); grid on;
    legend(gca,'show','Location','best');
    xlim([-2 0]);
    subplot(3,1,2);
    ylabel('Disp (m)', 'FontSize', 16);
    set(gca, 'FontSize', 18); grid on;
    legend(gca,'show','Location','best');
    xlim([-2 0]);
    subplot(3,1,3);
    ylabel('Disp (m)', 'FontSize', 16);
    set(gca, 'FontSize', 18); grid on;
    legend(gca,'show','Location','best');
    xlabel('Time (s, T1 as 0)', 'FontSize', 16);
    xlim([-2 0]);

    figure(velPlotHandle1);
    

    figure(velPlotHandle2);
    subplot(3,1,1);
    ylabel('Speed (m/s)', 'FontSize', 16);
    set(gca, 'FontSize', 18); grid on;
    legend(gca,'show','Location','best');
    xlim([-0.18 0]);
    subplot(3,1,2);
    ylabel('Speed (m/s)', 'FontSize', 16);
    set(gca, 'FontSize', 18); grid on;
    legend(gca,'show','Location','best');
    xlim([-0.18 0]);
    subplot(3,1,3);
    ylabel('Speed (m/s)', 'FontSize', 16);
    set(gca, 'FontSize', 18); grid on;
    legend(gca,'show','Location','best');
    xlabel('Disp (m, T1 as 0)', 'FontSize', 16);
    xlim([-0.18 0]);
    
    figure(velPlotHandle3);
%     subplot(3,1,1);
    ylabel('v (m/s)', 'FontSize', 16);
    xlabel('y (m, T1 as 0)', 'FontSize', 16);
    set(gca, 'FontSize', 18); grid on;
%     legend(gca,'show','Location','best');
    xlim([0.25 0.42]); ylim([0 0.7]);
    
    velPlotHandle2.CurrentAxes.GridAlpha = 0.2;
%     velPlotHandle2.CurrentAxes.GridColor = colors(ct_light,:);
    grid on; 
    vline(track.state_le(3),'k','LE');
    
    keyboard;
    close all;
end




%% Functions used in this script
function figHandle = createBoxPlot(variable, labels, yxislabel)
    % for boxplots of a cell array
    col=@(x)reshape(x,numel(x),1);
    boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});

    
    figHandle = figure; hold on;
    figure(figHandle);
    boxplot2(variable, 'Labels', labels, 'OutlierSize', 0.00001);
%     hbox = gca;
    set(gca, 'FontSize', 18); grid on;
    xlabel('Treatments', 'FontSize', 18);
    ylabel(yxislabel, 'FontSize', 18);
    % ylim([0 1.5]);
    for ct = 1:length(variable)
            x=(ct+(rand(length(variable{ct, 1}),1)-0.5)/4);

            f = scatter(x(:,1),variable{ct, 1},40,'k','filled'); 
            f.MarkerFaceAlpha = 0.5;
    %         keyboard;
    end
end

function setOpticalFlowFigure(plotHandle)
    figure(plotHandle);
    subplot(3,1,1);
    ylabel('V_{gy} / y (1/s)', 'FontSize', 16);
    set(gca, 'FontSize', 18); grid on;
    legend(gca,'show','Location','best');
    xlim([0 0.2]);
    ylim([-9 1]);

    subplot(3,1,2);
    ylabel('V_{gy} / y (1/s)', 'FontSize', 16);
    set(gca, 'FontSize', 18); grid on;
    legend(gca,'show','Location','best');
    xlim([0 0.2]);
    ylim([-9 1]);

    subplot(3,1,3);
    ylabel('V_{gy} / y (1/s)', 'FontSize', 16);
    set(gca, 'FontSize', 18); grid on;
    legend(gca,'show','Location','best');
    xlim([0 0.2]);
    ylim([-9 1]);
    xlabel('y (m)', 'FontSize', 16);
end

function plot3Dand2Ddata(video_metadata)
    % video_metadata - instance of videoMetadaat class
    rootDir = '/media/reken001/Disk_08_backup/light_intensity_experiments/';
    trackingDir = 'postprocessing';
%     videosDir = 'Videos';
    
    data2d_4video = video_metadata.data2d;
    data3d_4video = video_metadata.data3d;
    data_association_4video = video_metadata.data_association;
    
    % % Plotting 2D features
    cam_names = {'Basler_22549584', 'Basler_22549585', 'Basler_22549587', 'Basler_22956425'};
    cam_ids = [0:3];

    radius = 0.5;
%     data2d_plothandle = figure;
%     
%     for i=1:length(cam_names)
%         [Image.(['X' num2str(i)]), Image.(['map' num2str(i)])] = imread(fullfile(rootDir, trackingDir, video_metadata.trackingFilenName, 'images', [cam_names{i} '.png']));
%         subplot(2,2,i); imshow([Image.(['X' num2str(i)]), Image.(['map' num2str(i)])]);
% 
%         centers = data2d_4video(data2d_4video(:,1)==cam_ids(i),5:6);
%         radii = radius*ones(size(centers,1),1);
%         viscircles(centers, radii,'Color','g');
%     %     axis square;
%     %     axis equal;
%     %     axis image;
%     end
    
    % % Plotting 3Dtrack and associated features
    ids = unique(data3d_4video(:,1));
    colors = brewermap(length(ids),'Dark2');
    for i=1:length(ids)
        disp(i);
        id = ids(i);
        dummy = data3d_4video(data3d_4video(:,1)==id,:);
        
        trajFeaturePlot = figure;
        plot3(dummy(:,4),dummy(:,5),dummy(:,6),'-o','Color',colors(i,:),'LineWidth',2,'MarkerSize',1,'DisplayName',['id: ' num2str(id)]);
        set(gca, 'FontSize', 18); grid on;
        view(0,90);
        xlim([0 1]);
        ylim([0 0.5]);
        legend(gca,'show','Location','best');
        
        % Highlighting associated 2-D features
        data2d_plothandle = figure;
        for j=1:length(cam_names)
            % read images and plot all features
            [Image.(['X' num2str(j)]), Image.(['map' num2str(j)])] = imread(fullfile(rootDir, trackingDir, video_metadata.trackingFileName, 'images', [cam_names{j} '.png']));
            subplot(2,2,j); imshow([Image.(['X' num2str(j)]), Image.(['map' num2str(j)])]);

            centers = data2d_4video(data2d_4video(:,1)==cam_ids(j),5:6);
            radii = radius*ones(size(centers,1),1);
            viscircles(centers, radii,'Color','w');

            % plot 3d-track associated features
            frame_pointID = data_association_4video(data_association_4video(:,1)==ids(i) & ...
                                                    data_association_4video(:,3)==cam_ids(j),[2 4]);
            centers = zeros(length(frame_pointID),2);

            for k=1:size(frame_pointID,1)
                centers(k,1:2) = data2d_4video(data2d_4video(:,2)==frame_pointID(k,1) & ...
                                               data2d_4video(:,1)==cam_ids(j) & ...
                                               data2d_4video(:,10)==frame_pointID(k,2),5:6);
            end
            radii = radius*ones(size(centers,1),1);
            viscircles(centers, radii,'Color',colors(i,:));
        end
        % Plotting position
        pos_plothandle = figure;
        subplot(3,1,1);
        plot(dummy(:,3), dummy(:,4),'-o');
        ylabel('X');
        subplot(3,1,2);
        plot(dummy(:,3), dummy(:,5),'-o');
        ylabel('Y');
        subplot(3,1,3);
        plot(dummy(:,3), dummy(:,6),'-o');
        ylabel('Z');

        % Plotting velocity
        vel_plothandle = figure;
        subplot(3,1,1);
        plot(dummy(:,3), dummy(:,7),'-o');
        ylabel('V_{x}');
        subplot(3,1,2);
        plot(dummy(:,3), dummy(:,8),'-o');
        ylabel('V_{y}');
        subplot(3,1,3);
        plot(dummy(:,3), dummy(:,9),'-o');
        ylabel('V_{z}');

        % Plotting position covariance
        posCov_plothandle = figure;
        subplot(3,1,1);
        plot(dummy(:,3), dummy(:,10),'-o');
        ylabel('\sigma_{xx}');
        subplot(3,1,2);
        plot(dummy(:,3), dummy(:,13),'-o');
        ylabel('\sigma_{yy}');
        subplot(3,1,3);
        plot(dummy(:,3), dummy(:,15),'-o');
        ylabel('\sigma_{zz}');

        % Plotting position covariance
        velCov_plothandle = figure;
        subplot(3,1,1);
        plot(dummy(:,3), dummy(:,16),'-o');
        ylabel('\sigma_{uu}');
        subplot(3,1,2);
        plot(dummy(:,3), dummy(:,17),'-o');
        ylabel('\sigma_{vv}');
        subplot(3,1,3);
        plot(dummy(:,3), dummy(:,18),'-o');
        ylabel('\sigma_{ww}');

        keyboard;

        close(data2d_plothandle);
        close(trajFeaturePlot);
        close(pos_plothandle);
        close(vel_plothandle);
        close(posCov_plothandle);
        close(velCov_plothandle);


    end

end