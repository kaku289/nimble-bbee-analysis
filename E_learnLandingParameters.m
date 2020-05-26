%% Learn parameters of dynamic model during landing

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

%% Select trajectories for which parameters of landing model are required to be estimated
close all;
clc; clear;

% Inputs
inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/videoscoring_and_tracks_sysID.mat';
load(inputFile);


% pattern = {'checkerboard', 'spokes'};
% light = {'low', 'medium', 'high'};

labels = cell(0, 1); % for x-axis of boxplots
pattern_label = {'+', 'x'}; % + is checkerboard, x is spokes
light_label = {'L', 'M', 'H'};

distance_equinox = cell(0, 1);

pattern = {'spokes'};
light = {'high'};
 
% selectedVideoFiles = VideoMetadata.empty;
colormap = {brewermap(11,'PiYG'), ...
            brewermap(11,'RdBu')};
        
tempPlot = figure;        
        
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
               ~isempty(x.scoredLandingParams) % && isequal(x.scoredLandingParams.landing_success,2) % ...
%                && strcmpi(x.scoredLandingParams.landing_location,'Tube')
                isUsefulTrack(ct) = true;
            end           
        end
        
        relevantTracks = relevantTracks(isUsefulTrack);
        
        labels{end+1, 1} = [pattern_label{ct_pattern} ' ' light_label{ct_light}];
        
        distance_equinox{end+1, 1} = arrayfun(@(x) x.equinox_state_LDF(:,3), relevantTracks);
        
        figure(tempPlot); hold on
        for ct_track=1:length(relevantTracks)
            track = relevantTracks(ct_track);
            plotHandles = track.plotData_LDF();
            keyboard;
            close(plotHandles);
%             plot(track.stable_state_LDF(:,3),track.stable_state_LDF(:,6)./track.stable_state_LDF(:,3) ...
%                     ,'Color',colors(ct_light,:),'LineWidth',2);
%             xlim([0 0.2]);
%             keyboard;
        end
        
%         for ct_track=1:5%length(relevantTracks)
%             track = relevantTracks(ct_track);
%             
%             figure(DynamicVarPlotHandle(ct_pattern));
%             subplot(3,1,ct_light);
%             
%             if ct_track == 1
%                 plot(track.stable_state_LDF(:,3),track.stable_state_LDF(:,6)./track.stable_state_LDF(:,3) ...
%                     ,'Color',colors(ct_light,:),'LineWidth',2,'DisplayName',[pattern_label{ct_pattern} ' ' light_label{ct_light}]);
%             else
%                 plot(track.stable_state_LDF(:,3),track.stable_state_LDF(:,6)./track.stable_state_LDF(:,3) ...
%                     ,'Color',colors(ct_light,:),'LineWidth',2,'HandleVisibility','off');
%             end
% %             keyboard;
%         end
        
        
    end
end

% y_equinox_figHandle = createBoxPlot(distance_equinox, labels, 'Distance from the platform at |V_{gy}| < 0.025 m/s for the 1st time');
% arrayfun(@(x) setOpticalFlowFigure(x),DynamicVarPlotHandle);

%% Open track_extrac_GUIDE to mark y_start and y_end for the tracks
% After marking go to the next section

%% Learn parameters 
close all;
% clc; clear;

DirPlots = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/plots/parameterEstimation';
savePlots = true;

% Inputs
% inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/videoscoring_and_tracks_sysID.mat';
% load(inputFile);

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};

labels = cell(0, 1); % for x-axis of boxplots
pattern_label = {'+', 'x'}; % + is checkerboard, x is spokes
light_label = {'L', 'M', 'H'};

distance_equinox = cell(0, 1);

% pattern = {'spokes'};
% light = {'high'};
 
% selectedVideoFiles = VideoMetadata.empty;
colormap = {brewermap(11,'PiYG'), ...
            brewermap(11,'RdBu')};

for ct_pattern = 2%1:length(pattern)
    
%     colors = colormap{ct_pattern}(9:11,:);
%     DynamicVarPlotHandle(ct_pattern) = figure; % to plot displacement with time
%     subplot(3,1,1); hold on; % to plot low light condition
%     subplot(3,1,2); hold on; % to plot medium light condition
%     subplot(3,1,3); hold on; % to plot high light condition
    
    for ct_light = 3%1:length(light)
        disp(['%%%%%%%%%% Treatment: ' pattern{ct_pattern} ' and ' light{ct_light} ' light condition %%%%%%%%%%%%']);
        
        % Create sub-directory for each treatment if it doesn't exist
        if savePlots && ~exist(fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light}]), 'dir')
            mkdir(fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light}]));
        end
        DirPlots_treatment = fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light}]);
        
        relevantTreatments = treatments([treatments.videosParsed] & ...
                             strcmpi({treatments.pattern}, pattern{ct_pattern}) & ...
                             strcmpi({treatments.light}, light{ct_light}));
                         
        relevantVideos = [relevantTreatments.usefulVideoFiles];

        % Uncomment to check how many tracks have been marked with
        % y_start and y_end using visual inspection
        relevantTracks = [relevantVideos.landingTrack];
        isUsefulTrack = false(length(relevantTracks),1);
        for ct=1:length(relevantTracks)
            x = relevantTracks(ct);
            if ~isempty(x.stable_state_LDF) && x.confirmedChecks_le && ...
               ~isempty(x.scoredLandingParams) && ~isempty(x.track_subset_sysID) && ...
               abs(x.track_subset_sysID.y_start) >= 1e-5 && ...
               abs(x.track_subset_sysID.y_end) >= 1e-5
%                && strcmpi(x.scoredLandingParams.landing_location,'Tube')
                isUsefulTrack(ct) = true;
            end           
        end
        relevantTracks = relevantTracks(isUsefulTrack);
        
% %         tempPlot = figure;
% %         for ct=1:length(relevantTracks)
% %             figure(tempPlot)
% %             track = relevantTracks(ct);
% %             % Finding iodata for system identification
% %             % for y_start_indx: find the one that is closest to the
% %             % marked point
% %             [~, y_start_indx] = min(abs(track.state_LDF(:,3)-track.track_subset_sysID.y_start));
% %             % for y_end_indx - starting from y_start_indx, find earliest time instant that is
% %             % closest (within 1e-3) to marked y_end
% %             [ y_end_indx] = find(abs(track.state_LDF(y_start_indx+1:end,3)-track.track_subset_sysID.y_end)<1e-3,1);
% %             y_end_indx = y_end_indx + y_start_indx;
% %             track.track_subset_sysID.state = track.state_LDF(y_start_indx:y_end_indx, :);
% %             
% %             y_neg_indx = find(track.state_LDF(1:y_start_indx,1) - track.state_LDF(y_start_indx,1) < -0.1, 1, 'last');
% %             track.track_subset_sysID.negTime_r = track.state_LDF(y_neg_indx:y_start_indx, :);
% % %             plot(track.track_subset_sysID.state(:,1)-track.track_subset_sysID.state(1,1), track.track_subset_sysID.state(:,6)./track.track_subset_sysID.state(:,3));
% % %             keyboard
% %         end
% %         close(tempPlot);
% %         
% %         relevantTracks = relevantTracks([2 5 8 9 11 12 14]);
% %         for ct=1:length(relevantTracks)
% %             track = relevantTracks(ct);
% %             track.plotData_Time();
% % %             outputSignalForEstimation = 'r';
% % %             [vOpt, opt_info, params, plotHandles] = landingParameterEstimation_oneOutput(track.track_subset_sysID, false, outputSignalForEstimation);
% %                 
% %         end
% %         keyboard

        
%         % % % Save moving variance plots % % %
%         for ct_video=1:length(relevantVideos)
%             videoFile = relevantVideos(ct_video);
%             track = videoFile.landingTrack;
%             
%             if ~isempty(track.track_subset_sysID) && abs(track.track_subset_sysID.y_start) >= 1e-5 && ...
%                abs(track.track_subset_sysID.y_end) >= 1e-5
% %                 plotHandles = track.plotData_Distance();
%                 plotHandles = track.plotData_Time();
% %                 figure(plotHandle);
% %                 figure('units','normalized','outerposition',[0 0 1 1]);
%                 
%                 if savePlots
%                     for i=1:length(plotHandles)
%                         % Resizing the figures
%                         plotHandles(i).Position(3) = 680;
%                         plotHandles(i).Position(4) = 545;
%                         
%                         figureName = ['Variance_' videoFile.datenum '_' ...
%                             videoFile.scoredData.name(1:6) '_' ...
%                             videoFile.scoredData.landingParams.landing_side '_' ...
%                             videoFile.scoredData.landingParams.landing_location '.png'];
%                         
%                         saveas(plotHandles(i), fullfile(DirPlots_treatment, figureName) ,'png');
%                     end
%                 end
%                 close(plotHandles);
%                 
%             end
%         end
%         % % % % % % % % % % % % 
           
           
        
        for ct_video=1:length(relevantVideos)
            videoFile = relevantVideos(ct_video);
            track = videoFile.landingTrack;
            
            testPlot = figure;
            if ~isempty(track.track_subset_sysID) && abs(track.track_subset_sysID.y_start) >= 1e-5 && ...
               abs(track.track_subset_sysID.y_end) >= 1e-5  ...
               &&  strcmpi(videoFile.datenum, '20190708') && strcmpi(videoFile.scoredData.name(1:6), '093810')
%                && strcmpi(x.scoredLandingParams.landing_location,'Tube')

                % tracks for which grey-box modeling will be performed
                disp(['Analysing, day: '  videoFile.scoredData.year ...
                    videoFile.scoredData.month  videoFile.scoredData.day ...
                    ', video: ' videoFile.scoredData.name]);
                
                % Finding iodata for system identification
                % for y_start_indx: find the one that is closest to the
                % marked point
                [~, y_start_indx] = min(abs(track.state_LDF(:,3)-track.track_subset_sysID.y_start));
                % for y_end_indx - starting from y_start_indx, find earliest time instant that is
                % closest (within 1e-3) to marked y_end           
                [ y_end_indx] = find(abs(track.state_LDF(y_start_indx+1:end,3)-track.track_subset_sysID.y_end)<1e-3,1);
                y_end_indx = y_end_indx + y_start_indx;
                track.track_subset_sysID.state = track.state_LDF(y_start_indx:y_end_indx, :);
                
                % Plotting data for visual verification
                plot(track.track_subset_sysID.state(:,1)-track.track_subset_sysID.state(1,1), track.track_subset_sysID.state(:,6)./track.track_subset_sysID.state(:,3));
%                 keyboard;
                
                % % % Estimating parameters
                % % Way ONE
                outputSignalForEstimation = 'r';
                [vOpt, opt_info, params, plotHandles] = landingParameterEstimation_oneOutput(track.track_subset_sysID, false, outputSignalForEstimation);
                
                % saving data
                track.track_subset_sysID.param_estimation(1) = outputForParameterEstimation();
                param_estimation = track.track_subset_sysID.param_estimation(1);
                param_estimation.signalForEstimation = outputSignalForEstimation; % 'r' or 'ay'
                param_estimation.param_opt = vOpt; % values of optimized parameters
                param_estimation.params = params; % values of all parameters
                param_estimation.opt_info = opt_info; % optimization info
                
                % saving plots on disk
                if savePlots
                    for i=1:length(plotHandles)
                        % Resizing the figures
                        plotHandles(i).Position(3) = 680;
                        plotHandles(i).Position(4) = 545;
                        
                        if i==1
                            figureName = ['r_' videoFile.datenum '_' ...
                                videoFile.scoredData.name(1:6) '_' ...
                                videoFile.scoredData.landingParams.landing_side '_' ...
                                videoFile.scoredData.landingParams.landing_location ...
                                '_' outputSignalForEstimation '.png'];
                        elseif i==2
                            figureName = ['traj_' videoFile.datenum '_' ...
                                videoFile.scoredData.name(1:6) '_' ...
                                videoFile.scoredData.landingParams.landing_side '_' ...
                                videoFile.scoredData.landingParams.landing_location ...
                                '_' outputSignalForEstimation '.png'];
                        else
                            error('More than 2 figures encountered. Not possible!');
                        end
                        
                        saveas(plotHandles(i), fullfile(DirPlots_treatment, figureName) ,'png');
                    end
                end
                
                
                % % Way TWO
                outputSignalForEstimation = 'ay';
                [vOpt, opt_info, params, plotHandles] = landingParameterEstimation_oneOutput(track.track_subset_sysID, false, outputSignalForEstimation);
                
                % saving data
                track.track_subset_sysID.param_estimation(2) = outputForParameterEstimation();
                param_estimation = track.track_subset_sysID.param_estimation(2);
                param_estimation.signalForEstimation = outputSignalForEstimation; % 'r' or 'ay'
                param_estimation.param_opt = vOpt; % values of optimized parameters
                param_estimation.params = params; % values of all parameters
                param_estimation.opt_info = opt_info; % optimization info
                
                % saving plots on disk
                if savePlots
                    for i=1:length(plotHandles)
                        % Resizing the figures
                        plotHandles(i).Position(3) = 680;
                        plotHandles(i).Position(4) = 545;
                        
                        if i==1
                            figureName = ['r_' videoFile.datenum '_' ...
                                videoFile.scoredData.name(1:6) '_' ...
                                videoFile.scoredData.landingParams.landing_side '_' ...
                                videoFile.scoredData.landingParams.landing_location ...
                                '_' outputSignalForEstimation '.png'];
                        elseif i==2
                            figureName = ['traj_' videoFile.datenum '_' ...
                                videoFile.scoredData.name(1:6) '_' ...
                                videoFile.scoredData.landingParams.landing_side '_' ...
                                videoFile.scoredData.landingParams.landing_location ...
                                '_' outputSignalForEstimation '.png'];
                        else
                            error('More than 2 figures encountered. Not possible!');
                        end
                        
                        saveas(plotHandles(i), fullfile(DirPlots_treatment, figureName) ,'png');
                    end
                end
                close(plotHandles);
                
            end
            close(testPlot);
            
        end
        
        
        
%         relevantTracks = [relevantVideos.landingTrack];
%         
%         isUsefulTrack = false(length(relevantTracks),1);
%         for ct=1:length(relevantTracks)
%             x = relevantTracks(ct);
%             if ~isempty(x.track_subset_sysID) && abs(x.track_subset_sysID.y_start) >= 1e-5 && ...
%                abs(x.track_subset_sysID.y_end) >= 1e-5 % ...
% %                && strcmpi(x.scoredLandingParams.landing_location,'Tube')
%                 isUsefulTrack(ct) = true;
%             end           
%         end        
%         relevantTracks = relevantTracks(isUsefulTrack); 
%         
%         if isempty(relevantTracks)
%             error('No relevantTracks found. Expand/Change your search criteria!');
%         end
%         
        
%         figure;
%         for ct_track=3%1:length(relevantTracks)
%             % extract the part of the track to be used for system
%             % identification
%             track = relevantTracks(ct_track);
%             [~, y_end_indx] = min(abs(track.state_LDF(:,3)-track.track_subset_sysID.y_end));
%             [~, y_start_indx] = min(abs(track.state_LDF(y_end_indx+1:end,3)-track.track_subset_sysID.y_start));
%             y_start_indx = y_start_indx + y_end_indx;
%             
%             track.track_subset_sysID.state = track.state_LDF(y_end_indx:y_start_indx, :);
%             
% %             open_system('landingDynamics');
% %             figure;
% %             plot(track.track_subset_sysID.state(:,1)-track.track_subset_sysID.state(1,1), track.track_subset_sysID.state(:,3));
% %             pause;
% %             plot(track.track_subset_sysID.state(:,1)-track.track_subset_sysID.state(1,1), track.track_subset_sysID.state(:,6));
% %             pause;
%             plot(track.track_subset_sysID.state(:,1)-track.track_subset_sysID.state(1,1), track.track_subset_sysID.state(:,6)./track.track_subset_sysID.state(:,3));
%             pause;
%             
%             outputSignalForEstimation = 'r';
%             [vOpt, opt_info, figureHandles] = landingParameterEstimation_oneOutput(track.track_subset_sysID, false, outputSignalForEstimation);
%             
%             % saving data 
%             track.track_subset_sysID.param_estimation(1) = outputForParameterEstimation();
%             param_estimation = track.track_subset_sysID.param_estimation(1);
%             param_estimation.signalForEstimation = outputSignalForEstimation; % 'r' or 'ay'
%             param_estimation.param_opt = vOpt; % values of optimized parameters
% %             param_estimation.params = ; % values of all parameters
%             param_estimation.opt_info = opt_info; % optimization info
%             
            
            
%             
%             figure; plot(track.track_subset_sysID.state(:,3), track.track_subset_sysID.state(:,6)./track.track_subset_sysID.state(:,3))
%             keyboard; close all;
    end
end
keyboard;
inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/videoscoring_and_tracks_sysID.mat';
save(inputFile, 'treatments');
bdclose('landingDynamics');

%% Select specific tracks 

close all;
% clc; clear;

DirPlots = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/plots/parameterEstimation';
savePlots = true;

% Inputs
% inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/videoscoring_and_tracks_sysID.mat';
% load(inputFile);

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};

labels = cell(0, 1); % for x-axis of boxplots
pattern_label = {'+', 'x'}; % + is checkerboard, x is spokes
light_label = {'L', 'M', 'H'};

distance_equinox = cell(0, 1);

% pattern = {'spokes'};
% light = {'high'};
 
% selectedVideoFiles = VideoMetadata.empty;
colormap = {brewermap(11,'PiYG'), ...
            brewermap(11,'RdBu')};

for ct_pattern = 2%1:length(pattern)
    
%     colors = colormap{ct_pattern}(9:11,:);
%     DynamicVarPlotHandle(ct_pattern) = figure; % to plot displacement with time
%     subplot(3,1,1); hold on; % to plot low light condition
%     subplot(3,1,2); hold on; % to plot medium light condition
%     subplot(3,1,3); hold on; % to plot high light condition
    
    for ct_light = 3%1:length(light)
        disp(['%%%%%%%%%% Treatment: ' pattern{ct_pattern} ' and ' light{ct_light} ' light condition %%%%%%%%%%%%']);
        
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
            
            
            if strcmpi(videoFile.datenum, '20190708') && strcmpi(videoFile.scoredData.name(1:6), '094953')
%                && strcmpi(x.scoredLandingParams.landing_location,'Tube')
                track.plotData_Distance();
                track.plotData_Time1();
                
                
                % % Save data for parameter estimation using GUI % %
                state_LDF = track.state_LDF(track.time>=-1.5,:);
                state_LDF(:,1) = state_LDF(:,1)-state_LDF(1,1);
                
                negTime = state_LDF(state_LDF(:,1)<=0.44 &  state_LDF(:,1)>=0.44-0.06,[1 3 6]);
                negTime_r = [negTime(:,1) negTime(:,3)./negTime(:,2)];
                state_LDF = state_LDF(state_LDF(:,1)>=0.44 &  state_LDF(:,1)<=0.84,:);
                negTime_r(:,1) = negTime_r(:,1)-state_LDF(1,1);
                state_LDF(:,1) = state_LDF(:,1)-state_LDF(1,1);
                negTime_r(end+1,:) = [state_LDF(1,1), state_LDF(1,6)./state_LDF(1,3)];
                
                iodata = [state_LDF(:,[1 3 6 9]) state_LDF(:,6)./state_LDF(:,3) state_LDF(:,9)./state_LDF(:,6)]; % [ time y Vgy ay Vgy/y ay/y]
%                 save(['example_track_' pattern{ct_pattern} '_' light{ct_light} '_' videoFile.datenum '_' videoFile.scoredData.name(1:6)],'track','iodata')
%                 keyboard;
                
                % % Estimate parameters using the rather "crude way" - vary
                % sensory delay in steps and compute PID parameters using
                % least squares
                % tau - sensory delay
                tau_step = 0.001; %  steps in which sensory delay is varied
                tau_init = 0;
                tau_max = 0.06;
                taus = tau_init:tau_step:tau_max;
                K = zeros(3,length(taus));
                r_ref = -5;
                
                % Plotting a(t) with a(t) computed using estimated
                % parameters
                accPlot = figure;
                
                % Storing parameters
                J_mse = zeros(length(taus), 1);
                J_mae = zeros(length(taus), 1);
                
                for ct=1:length(taus)
                    tau = taus(ct);
                    
                    % Based on a/v
%                     r_true = iodata(:,[1 6]);
                    
                    % Based on v/y
                    r_true = iodata(:,[1 5]);
                    
                    r_measured = [r_true(:,1) interp1(r_true(:,1), r_true(:,2), r_true(:,1)-tau, 'makima', NaN)];
                    
                    error = [r_measured(:,1) r_ref-r_measured(:,2)];
                    error(isnan(error(:,2)),:) = [];
                    error_sum = cumtrapz(error(:,1), error(:,2));
                    error_diff = diffxy(error(:,1), error(:,2));
                    
                    A = [error(:,2) error_sum error_diff];
                    
                    B = iodata(iodata(:,1)>=tau,4);
                    
                    assert(size(A,1) == size(B,1));
                    
                    % Computing Kp, Ki, Kd for each tau
                    K(:,ct) = A\B;
                    
                    % Plotting
                    disp(tau);
                    figure(accPlot);
                    cla;
                    plot(iodata(:,1), iodata(:,4), 'b'); hold on;
%                     plot(error(:,1), B, 'b'); hold on;
                    plot(error(:,1), A*K(:,ct), 'r');
                    
                    J_mse(ct) = sum((B - A*K(:,ct)).^2)/length(B);
                    J_mae(ct) = sum(abs(B - A*K(:,ct)))/length(B);
                    
%                     keyboard
                    
                end
                
                figure;
                plot(taus, J_mse, 'b'); hold on;
                plot(taus, J_mae, 'r');
                
                
                keyboard;
                
            end
            close all;
            
            
        end
        
        
    end
end


%% Functions used in this script
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
