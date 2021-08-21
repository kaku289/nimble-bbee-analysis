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

% % to include loSR Surrey for boxplot
% addpath('./lib/MatlabToolbox');
%% Find r* entry intervals
% Also find videos asociated with each track

% Inputs
close all; clc;
% clear;

if isunix
    inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A1_rref.mat';
    outputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A2_rrefEntry.mat';
    DirPlots = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/plots/BlindLandingTracks/rref_entry_with_time';
    % DirPlots = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/plots/BlindLandingTracks/rref_entry';
    % DirPlots = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/plots/BlindLandingTracks/rref_estimate_3dspeed';
elseif ispc
    inputFile = 'D:/light_intensity_experiments/postprocessing/BlindLandingtracks_A1_rref.mat';
    outputFile = 'D:/light_intensity_experiments/postprocessing/BlindLandingtracks_A2_rrefEntry.mat';
    DirPlots = 'D:/light_intensity_experiments/postprocessing/plots/BlindLandingTracks/rref_entry_with_time';
end


% load(inputFile);


delPreviousPlots = false; % BE CAREFUL - when set to true, all previously saved plots are deleted
savePlots = false;
savePDFs = false;

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
behaviour = {'rising','constant','sleeping'};

factors = [0.25:0.25:2.5];
factors = [1];
for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        for ct_behaviour = 2%1:length(behaviour)
            
            % Delete previous plots
            if delPreviousPlots && savePlots && exist(fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour}]), 'dir')
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
            
            
            for ct_treatment=1:length(relevantTreatments)
                treatment = relevantTreatments(ct_treatment);
                
                videoTimes = [[treatment.videosInfo.startTime]' [treatment.videosInfo.endTime]'];
                if ~isempty(treatment.landingTracks)
                    disp(['Into, day: ' num2str(treatment.datenum) ...
                          ', pattern: ' treatment.pattern ...
                          ', light: ' treatment.light ...
                          ', behaviour: ' behaviour{ct_behaviour}]);
                    
                      for ct_track=1:length(treatment.landingTracks) % for each landing track
                          track = treatment.landingTracks(ct_track);
                          for ct_excerpt=1:length(track.state_LDF) % for each track excerpt
                              excerpt = track.state_LDF(ct_excerpt);
                              
                              excerpt.find_rrefEntry();
                              
                              if savePlots
                                  for ct_factor=1:length(factors)

                                      plotHandles = excerpt.plot_rrefsEntry_with_time(factors(ct_factor));
%                                       plotHandles = excerpt.plot_rrefsEntry(factors(ct_factor));
%                                       plotHandles = excerpt.plot_rrefs_with3dspeed(factors(ct_factor));

                                      if ~isempty(plotHandles)

                                          % Resizing the figures
                                          for i=1:length(plotHandles)
                                              plotHandles(i).Position(3) = 680;
                                              plotHandles(i).Position(4) = 545;

                                              if i==1
                                                  figureName = ['fac_' num2str(factors(ct_factor),'%0.2f') '_' ...
                                                              num2str(treatment.datenum) '_' ...
                                                              num2str(treatment.startTime) '_' num2str(treatment.endTime) ...
                                                              '_obj' num2str(track.obj_id) '_track' ...
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
            end % for treatment
            
        end
    end
end
keyboard
save(outputFile, 'treatments');
keyboard;



%% Run estimation of entry dynamics
% a) data is prefiltered
% b) estimation is done for individual entry segments
% c) estimation is done for three orders of transfer functions and 10
% factors

clc; close all;
% clear;

if isunix
    inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A2_rrefEntryEstimation_everything_f1pt5_all.mat';
elseif ispc
    inputFile = 'D:/light_intensity_experiments/postprocessing/BlindLandingtracks_A2_rrefEntryEstimation_everything_f1pt5_all.mat';
end

% inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A2_rrefEntry.mat';
% load(inputFile);
treatments = treatments(1:14*8); % Taking experiments for 2 patterns * 3 lights


pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
behaviour = {'rising','constant','sleeping'};

% factors = [0.5:0.5:2.5];
factors = [1.5]; % chosen factor

data_all = struct.empty;

clear data4est data4est_lowpass estimatedData data1 data1_filtered indexes; clc;
data4est = cell(length(pattern), length(light), length(factors)); % Unmerged dataset
data4est_lowpass = cell(length(pattern), length(light), length(factors)); % Unmerged dataset prefiltered with low bandpass filter 
fc = 5; % cut-off frequency in Hz (used for prefiltering the data)

np = 1:3; % number of poles
nz = 0; % number of zeros
nPZ = struc(np, nz);
nPZ(nPZ(:,2)>nPZ(:,1),:) = [];

opt = tfestOptions('InitializeMethod','all','Display','off','SearchMethod','auto');
opt.SearchOptions.MaxIterations = 3000;
warning('off','Ident:estimation:transientDataCorrection')
warning('off','Ident:idmodel:size4isNotNX')
parpool
pause(10)
parfevalOnAll(gcp(), @warning, 0, 'off', 'Ident:estimation:transientDataCorrection')
% parfevalOnAll(gcp(), @warning, 0, 'off', 'Ident:idmodel:size4isNotNX')
% tic
for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        for ct_behaviour = 2%1:length(behaviour)
            disp(' ');
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
            
            %%%%%%%%%%%%%% Each light and pattern combination %%%%%%%%%%%
            landingTracks = [relevantTreatments.landingTracks];
            state_LDF = [landingTracks.state_LDF];
            rrefEntrySegs = [state_LDF.rrefEntrySegments];
            
            landingTracks_indx4stateLDF = arrayfun(@(x) x*ones(length(landingTracks(x).state_LDF),1),1:length(landingTracks), 'UniformOutput', false);
            landingTracks_indx4stateLDF = vertcat(landingTracks_indx4stateLDF{:});
            
            % Finding # of "tracks" (state LDFs) that contain rref entry segments
            % for each factor
            for ct_fac = 1:length(factors)
                
                factor = factors(ct_fac);
                disp(['Running factor: ' num2str(factor)]);
                
                indices = arrayfun(@(x) ~isempty(x.rrefEntrySegments(abs([x.rrefEntrySegments.factor]-factor)<1e-6).intervals),state_LDF);
                
                data_all(ct_pattern, ct_light, ct_fac).tracks_fac = state_LDF(indices);
                data_all(ct_pattern, ct_light, ct_fac).ct_pattern = ct_pattern;
                data_all(ct_pattern, ct_light, ct_fac).ct_light = ct_light;
                data_all(ct_pattern, ct_light, ct_fac).ct_factor = factor;
                data_all(ct_pattern, ct_light, ct_fac).landingTrack = landingTracks(landingTracks_indx4stateLDF(indices));
                
                % Extract data for system identification
%                 clear estimatedData;
                estimatedData.iddata = {};
                estimatedData_prefiltered.iddata = {};
                indexes = [];
                
                for ct1=1:length(data_all(ct_pattern, ct_light, ct_fac).tracks_fac)
                    state = data_all(ct_pattern, ct_light, ct_fac).tracks_fac(ct1).filteredState;
                    rrefEntrySegments_fac = data_all(ct_pattern, ct_light, ct_fac).tracks_fac(ct1).rrefEntrySegments(abs([data_all(ct_pattern, ct_light, ct_fac).tracks_fac(ct1).rrefEntrySegments.factor]-factor)<1e-6);
                    
%                     assert(length(rrefEntrySegments_fac) == 1);
                    for ct2=1:size(rrefEntrySegments_fac.intervals,1)
                        output = -state(rrefEntrySegments_fac.intervals(ct2,1):rrefEntrySegments_fac.intervals(ct2,3),6)./...
                            state(rrefEntrySegments_fac.intervals(ct2,1):rrefEntrySegments_fac.intervals(ct2,3),3);
                        input = -rrefEntrySegments_fac.rmean(ct2)*ones(length(output),1);
                        dt = mean(diff(state(rrefEntrySegments_fac.intervals(ct2,1):rrefEntrySegments_fac.intervals(ct2,3),1)));
                        
                        [num, den] = butter(2, fc/(1/dt/2),'low');
                        filt_output = filtfilt(num, den, output);
                        
                        data1 = iddata(output, input, dt);
                        data1_filtered = iddata(filt_output, input, dt);
                        
                        estimatedData.iddata{end+1} = data1;
                        estimatedData_prefiltered.iddata{end+1} = data1_filtered;
                        
                        indexes(end+1) = ct1;                        
                    end
                end
                data4est{ct_pattern, ct_light, ct_fac} = estimatedData;
                data4est_lowpass{ct_pattern, ct_light, ct_fac} = estimatedData_prefiltered;
                
                data4est{ct_pattern, ct_light, ct_fac}.track_indexes = indexes;
                data4est_lowpass{ct_pattern, ct_light, ct_fac}.track_indexes = indexes;
                
                % Run system Identification
                dummy = idtf(zeros(0,0,0));
                for ct2=1:size(nPZ,1)
                    parfor ct1=1:length(estimatedData_prefiltered.iddata)
                        % ADD initial condition as output. see tfest form
                        % in help - DO THIS BEFORE RUNNING IT AGAIN
                        dummy(:,:,ct1,ct2) = tfest(estimatedData_prefiltered.iddata{ct1}, nPZ(ct2,1), nPZ(ct2,2), opt);
                    end
                end
                data4est_lowpass{ct_pattern, ct_light, ct_fac}.tfest = dummy;
                
            end
            
           
            
        end
        
%         outputFile = ['/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A2_rrefEntryEstimation_' pattern{ct_pattern} ...
%                   '_' light{ct_light} '.mat'];
%         save(outputFile,'data4est','data4est_lowpass','data_all');      
    end
%     keyboard;
end
% toc
keyboard

outputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A2_rrefEntryEstimation_f1.mat';
save(outputFile,'data4est','data4est_lowpass','data_all');

outputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A2_rrefEntryEstimation_everything_f1.mat';
save(outputFile);
keyboard;
%%
labels = cell(0, 1); % for x-axis of boxplots
% pattern_label = {'+', 'x'}; % + is checkerboard, x is spokes
% light_label = {'L', 'M', 'H'};

% AICc = cell(length(data),1);
% labels = cell(length(data),1);


    
for ct_fac=1:length(factors)
    
    AICc = cell(size(nPZ,1),1);
%     AIC = cell(size(nPZ,1),1);
%     BIC = cell(size(nPZ,1),1);
    fitPercent = cell(size(nPZ,1),1);
    for ct_order = 1:size(nPZ,1)
        labels{ct_order, 1} = [num2str(ct_order)];
        
        for ct_pattern = 1:length(pattern)
            for ct_light = 1:length(light)
                
                n = length(data4est_lowpass{ct_pattern, ct_light, ct_fac}.track_indexes);
                AICc{ct_order, 1} = [AICc{ct_order, 1} arrayfun(@(x) data4est_lowpass{ct_pattern, ct_light, ct_fac}.tfest(:,:,x,ct_order).Report.Fit.AICc,1:n)];
%                 AIC{ct_order, 1} = [AIC{ct_order, 1} arrayfun(@(x) data4est_lowpass{ct_pattern, ct_light, ct_fac}.tfest(:,:,x,ct_order).Report.Fit.AIC,1:n)];
%                 BIC{ct_order, 1} = [BIC{ct_order, 1} arrayfun(@(x) data4est_lowpass{ct_pattern, ct_light, ct_fac}.tfest(:,:,x,ct_order).Report.Fit.BIC,1:n)];
                fitPercent{ct_order, 1} = [fitPercent{ct_order, 1} arrayfun(@(x) data4est_lowpass{ct_pattern, ct_light, ct_fac}.tfest(:,:,x,ct_order).Report.Fit.FitPercent,1:n)];
            end
        end

    end
%     createBoxPlot(BIC, labels, 'BIC')
%     createBoxPlot(AIC, labels, 'AIC')
    createBoxPlot(AICc, labels, 'AICc')
%     createBoxPlot(fitPercent, labels, 'FitPercent')
    title(['Factor = ' num2str(factors(ct_fac))], 'FontSize', 18);
end

keyboard;
for ct=1:length(data4est_unmerged)
    
    parfor ct1=1:length(data4est_unmerged{ct}.iddata)
        
        dummy(:,:,ct1) = tfest(data4est_unmerged{ct}.iddata{ct1}, 2, 0, opt);
             
    end
    data4est_unmerged{ct}.tfest = dummy;
end
toc


for ct=1%:length(data4est_unmerged)
%     data4est_unmerged{ct}.tfest(:,:,size(nPZ,1),length(data4est_unmerged{ct}.iddata)) = idtf(zeros(0,0,0));
    for ct1=1:length(data4est_unmerged{ct}.iddata)
        for ct2=1:size(nPZ,1)
            data4est_unmerged{ct}.tfest(:,:,ct2,ct1) = tfest(data4est_unmerged{ct}.iddata{ct1}, nPZ(ct2,1), nPZ(ct2,2), opt);
        end        
    end
end

%% Loading data and save plots for a specified factor
% Plots are plotted for each track and highlights rref segment(s), entry
% segment(s), filtered data, and data estimated from system identification
clc; close all;
% clear;
if isunix
    dataFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A2_rrefEntryEstimation_everything_f1.mat';
    DirPlots = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/plots/BlindLandingTracks/rref_entry_estimation2';
elseif ispc
    dataFile = 'D:/light_intensity_experiments/postprocessing/BlindLandingtracks_A2_rrefEntryEstimation_everything_f1.mat';
    DirPlots = 'D:/light_intensity_experiments/postprocessing/plots/BlindLandingTracks/rref_entry_estimation2';
end
% load(dataFile);

delPreviousPlots = true; % BE CAREFUL - when set to true, all previously saved plots are deleted
savePlots = true;
savePDFs = true;

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
chosen_fac = 1;
% behaviour = {'rising','constant','sleeping'};

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
        
        for ct1=1:length(data_all(ct_pattern,ct_light).tracks_fac)
            % create tf vector in order of rrefEntrySegments
            tfs = data4est_lowpass{ct_pattern,ct_light}.tfest(:,:,data4est_lowpass{ct_pattern,ct_light}.track_indexes == ct1, 2);
            sysiddata = data4est_lowpass{ct_pattern,ct_light}.iddata(data4est_lowpass{ct_pattern,ct_light}.track_indexes == ct1);
            
            % plot data
%             plotHandles = data_all(ct_pattern,ct_light).tracks_fac(ct1).plot_rrefsEntry_withActualFilteredEstimatedData(chosen_fac, sysiddata, tfs);
            plotHandles = data_all(ct_pattern,ct_light).tracks_fac(ct1).plot_rrefsEntry_withActualFilteredEstimatedData2(chosen_fac, sysiddata, tfs);
%             data_all(ct_pattern,ct_light).tracks_fac(ct1).plot_rrefsEntry_withSimulatedData(chosen_fac, tfs);
            
            if ~isempty(plotHandles)
                
                % Resizing the figures
                for i=1:length(plotHandles)
                    plotHandles(i).Position(3) = 680;
                    plotHandles(i).Position(4) = 545;
                    
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

%% Compare rise times from 3 following conditions:
% 1) data with tf simulation -  This is to make sure if second-order tf is a good approximation of the closed-loop control system
% 2) Actual data
% 3) Low-pass filtered data 
%
model_order_indx = 2;
labels = cell(0, 1); % for x-axis of boxplots
pattern_label = {'+', 'x'}; % + is checkerboard, x is spokes
light_label = {'L', 'M', 'H'};
diff_tr_actual_tf = cell(0, 1); % differences in rise times between actual data and second-order tf
diff_tr_actual_lowpass = cell(0, 1); % differences in rise times between actual data and second-order tf
dt = 0.00571966171264648;
rEdges = [0.5 0.7];
for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        labels{ct_pattern, ct_light} = [light_label{data_all(ct_pattern,ct_light).ct_light} '' pattern_label{data_all(ct_pattern,ct_light).ct_pattern}];
        
        n = length(data4est_lowpass{ct_pattern, ct_light, 1}.track_indexes);
        %                 AICc{ct_order, 1} = [AICc{ct_order, 1} arrayfun(@(x) data4est_lowpass{ct_pattern, ct_light, ct_fac}.tfest(:,:,x,ct_order).Report.Fit.AICc,1:n)];
        fitPercent{ct_pattern, ct_light} = arrayfun(@(x) data4est_lowpass{ct_pattern, ct_light}.tfest(:,:,x,model_order_indx).Report.Fit.FitPercent,1:n);
        
        for ct=1:n
            rref = data4est{ct_pattern, ct_light}.iddata{ct}.u(1);
            r = data4est{ct_pattern, ct_light}.iddata{ct}.y;
              
            if any(r<rEdges(1)*rref) && any(r>rEdges(2)*rref) % && fitPercent{ct_pattern, ct_light}(ct)>95
                tr_actual = diff(interp1(r, [0:length(r)-1]*dt, rEdges'*rref));
                tr_lowpass = diff(interp1(data4est_lowpass{ct_pattern, ct_light}.iddata{ct}.y, [0:length(r)-1]*dt, rEdges'*rref));
                [rsim,fit,rsim0] = compare(data4est_lowpass{ct_pattern, ct_light}.iddata{ct}, data4est_lowpass{ct_pattern, ct_light}.tfest(:,:,ct,model_order_indx));
                tr_simulated = diff(interp1(rsim.y, [0:length(rsim.y)-1]*dt, rEdges'*rref));
                diff_tr_actual_tf{ct_pattern, ct_light}(ct) = tr_actual - tr_simulated;   
                diff_tr_actual_lowpass{ct_pattern, ct_light}(ct) = tr_actual - tr_lowpass;   
            else
                diff_tr_actual_tf{ct_pattern, ct_light}(ct) = nan;  
                diff_tr_actual_lowpass{ct_pattern, ct_light}(ct) = nan;
            end
        end
    end
end

createBoxPlot({diff_tr_actual_tf{:}}', {labels{:}}', 'Diff. in rise times between simulated and actual (s)');
createBoxPlot({diff_tr_actual_lowpass{:}}', {labels{:}}', 'Diff. in rise times between lowpass and actual (s)');


%% Create fit percent boxplots and free parameters boxplots
% % Fit percent boxplots
% For each model order (treatments are pooled together for each model
% order)
for ct_fac=1:length(factors)
    
%     AICc = cell(size(nPZ,1),1);
    fitPercent = cell(size(nPZ,1),1);
    for ct_order = 1:size(nPZ,1)
        labels{ct_order, 1} = [num2str(ct_order)];
        
        for ct_pattern = 1:length(pattern)
            for ct_light = 1:length(light)
                
                n = length(data4est_lowpass{ct_pattern, ct_light, ct_fac}.track_indexes);
%                 AICc{ct_order, 1} = [AICc{ct_order, 1} arrayfun(@(x) data4est_lowpass{ct_pattern, ct_light, ct_fac}.tfest(:,:,x,ct_order).Report.Fit.AICc,1:n)];
                fitPercent{ct_order, 1} = [fitPercent{ct_order, 1} arrayfun(@(x) data4est_lowpass{ct_pattern, ct_light, ct_fac}.tfest(:,:,x,ct_order).Report.Fit.FitPercent,1:n)];
            end
        end

    end
%     createBoxPlot(AICc, labels, 'AICc')
    createBoxPlot(fitPercent, labels, 'FitPercent')
    title(['Factor = ' num2str(factors(ct_fac))], 'FontSize', 18);
end
%%
% For each treatment with a selected model order (=2)
chosen_fac = 1.5;

model_order_indx = 2;
labels = cell(0, 1); % for x-axis of boxplots
pattern_label = {'+', 'x'}; % + is checkerboard, x is spokes
light_label = {'L', 'M', 'H'};
fitPercent = cell(0, 1);
for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        labels{ct_pattern, ct_light} = [light_label{data_all(ct_pattern,ct_light).ct_light} '' pattern_label{data_all(ct_pattern,ct_light).ct_pattern}];
        
        n = length(data4est_lowpass{ct_pattern, ct_light, ct_fac}.track_indexes);
        %                 AICc{ct_order, 1} = [AICc{ct_order, 1} arrayfun(@(x) data4est_lowpass{ct_pattern, ct_light, ct_fac}.tfest(:,:,x,ct_order).Report.Fit.AICc,1:n)];
        fitPercent{ct_pattern, ct_light} = arrayfun(@(x) data4est_lowpass{ct_pattern, ct_light}.tfest(:,:,x,model_order_indx).Report.Fit.FitPercent,1:n);
    end
end
createBoxPlot({fitPercent{1,1:3} fitPercent{2,1:3}}', {labels{1,1:3} labels{2,1:3}}', 'FitPercent'); ylim([50 100])

% Free Parameters box plots (for selected model order and selected factor)
labels = cell(0, 1); % for x-axis of boxplots

fitParameters = cell(0, 1);
P1 = cell(0, 1);
P2 = cell(0, 1);
P3 = cell(0, 1);

gain = cell(0,1); % gain
zeta = cell(0,1); % damping ratio
omega = cell(0,1); % natural frequency

ystart_entry = cell(0,1); % starting y at rref entry
ystart_rref = cell(0,1); % starting y at rref entry


clear dummy;
for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        
        %Compute y_start for each rref entry segment
        ystart_indx = arrayfun(@(x) x.rrefEntrySegments(abs([x.rrefEntrySegments.factor]-1.5)<1e-5).intervals(:,1) , data_all(ct_pattern,ct_light).tracks_fac,'UniformOutput',false);
        ystart_indx = vertcat(ystart_indx{:});
        ystart_entry{ct_pattern, ct_light} = arrayfun(@(x) -data_all(ct_pattern,ct_light).tracks_fac(data4est_lowpass{ct_pattern, ct_light}.track_indexes(x)).filteredState(ystart_indx(x),3) ,1:length(ystart_indx));
        % 
        
        labels{ct_pattern, ct_light} = [light_label{data_all(ct_pattern,ct_light).ct_light} '' pattern_label{data_all(ct_pattern,ct_light).ct_pattern}];
        n = length(data4est_lowpass{ct_pattern, ct_light, ct_fac}.track_indexes);
        
        dummy = arrayfun(@(x)  data4est_lowpass{ct_pattern, ct_light}.tfest(:,:,x,model_order_indx).Report.Parameters.ParVector, 1:n, 'UniformOutput', false);
        fitParameters{ct_pattern, ct_light} = horzcat(dummy{:})';
        assert(all(fitParameters{ct_pattern, ct_light}(:,4) == 0));
        
        P1{ct_pattern, ct_light} = fitParameters{ct_pattern, ct_light}(:,1);
        P2{ct_pattern, ct_light} = fitParameters{ct_pattern, ct_light}(:,2);
        P3{ct_pattern, ct_light} = fitParameters{ct_pattern, ct_light}(:,3);
        
        gain{ct_pattern, ct_light} = P1{ct_pattern, ct_light}./P3{ct_pattern, ct_light};
        omega{ct_pattern, ct_light} = (P3{ct_pattern, ct_light}).^0.5;
        zeta{ct_pattern, ct_light} = P2{ct_pattern, ct_light}./P3{ct_pattern, ct_light}/2;
    end
end

fitParameter1_figHandle = createBoxPlot({P1{:}}', {labels{:}}', 'P1'); ylim([-Inf 1000])
fitParameter2_figHandle = createBoxPlot({P2{:}}', {labels{:}}', 'P2'); ylim([-Inf 60])
fitParameter3_figHandle = createBoxPlot({P3{:}}', {labels{:}}', 'P3'); ylim([-Inf 1000])
%% 
close all;
fitParameter1_figHandle = createBoxPlot({gain{:}}', {labels{:}}', 'gain'); ylim([-0.1 2])
fitParameter2_figHandle = createBoxPlot({omega{:}}', {labels{:}}', 'omega'); ylim([-Inf 60])
fitParameter3_figHandle = createBoxPlot({zeta{:}}', {labels{:}}', 'zeta'); ylim([-0.02 0.1])

%%
gain_figHandle = createBoxPlot({gain{:}}', {labels{:}}', 'Gain'); ylim([0 2])
omega_figHandle = createBoxPlot({omega{:}}', {labels{:}}', 'Natural frequency'); ylim([-Inf 40])
zeta_figHandle = createBoxPlot({zeta{:}}', {labels{:}}', 'Damping'); ylim([-Inf 0.1])
%%
figure;
plot(horzcat(ystart_entry{:}),vertcat(gain{:}),'.'); ylim([0 2]);
figure;
plot([ystart_entry{:}],vertcat(omega{:}),'.'); ylim([-Inf 60]);
figure;
plot([ystart_entry{:}],vertcat(zeta{:}),'.'); ylim([-Inf 0.2]);

%% Create contour plot of accleration on 
% 1) r/r* and y plane (distance to start of rref interval)
% 2) r/r* and time plane (time to start of rref interval)
pattern_label = {'+', 'x'}; % + is checkerboard, x is spokes
light_label = {'L', 'M', 'H'};

yentryend = cell(0,1); % starting y at rref entry
tentryend = cell(0,1); % starting y at rref entry
rref = cell(0,1);

state_entry = cell(0,1); %[time to rrefstart, distance to rref start, r/rref, ay]
ay_range = [];

clear dummy;
for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        
        %Compute y_start for each rref segment
        yentryend_indx = arrayfun(@(x) x.rrefEntrySegments([x.rrefEntrySegments.factor]==1).intervals(:,2) , data_all(ct_pattern,ct_light).tracks_fac,'UniformOutput',false);
        yentryend_indx = vertcat(yentryend_indx{:});
        
        yentrystart_indx = arrayfun(@(x) x.rrefEntrySegments([x.rrefEntrySegments.factor]==1).intervals(:,1) , data_all(ct_pattern,ct_light).tracks_fac,'UniformOutput',false);
        yentrystart_indx = vertcat(yentrystart_indx{:});
        
        dummy = arrayfun(@(x) x.rrefEntrySegments([x.rrefEntrySegments.factor]==1).rref , data_all(ct_pattern,ct_light).tracks_fac,'UniformOutput',false);
        rref{ct_pattern, ct_light} = vertcat(dummy{:});
        
        yentryend{ct_pattern, ct_light} = arrayfun(@(x) -data_all(ct_pattern,ct_light).tracks_fac(data4est_lowpass{ct_pattern, ct_light}.track_indexes(x)).filteredState(yentryend_indx(x),3) ,1:length(yentryend_indx));
        tentryend{ct_pattern, ct_light} = arrayfun(@(x) data_all(ct_pattern,ct_light).tracks_fac(data4est_lowpass{ct_pattern, ct_light}.track_indexes(x)).filteredState(yentryend_indx(x),1) ,1:length(yentryend_indx));
        
        dummy = arrayfun(@(x) data_all(ct_pattern,ct_light).tracks_fac(data4est_lowpass{ct_pattern, ct_light}.track_indexes(x)).filteredState(yentrystart_indx(x):yentryend_indx(x),:) ,1:length(yentryend_indx),'UniformOutput',false);
        
        state_entry{ct_pattern, ct_light} = arrayfun(@(x) [dummy{x}(:,1)-tentryend{ct_pattern,ct_light}(x)  -dummy{x}(:,3)-yentryend{ct_pattern,ct_light}(x)  dummy{x}(:,6)./dummy{x}(:,3)./rref{ct_pattern,ct_light}(x) dummy{x}(:,9)],1:length(dummy),'UniformOutput',false);
        dummy = vertcat(state_entry{ct_pattern,ct_light}{:});
        ay_range{ct_pattern, ct_light} = [min(dummy(:,4)) max(dummy(:,4))];
        % 
        
        labels{ct_pattern, ct_light} = [light_label{data_all(ct_pattern,ct_light).ct_light} '' pattern_label{data_all(ct_pattern,ct_light).ct_pattern}];
    end
end

% % per treatment plots
% for ct_pattern = 1%:length(pattern)
%     for ct_light = 1%:length(light)
%         figHandles(ct_pattern,ct_light) = figure; hold on;
%         state = state_entry{ct_pattern,ct_light};
%         for ct=1:length(state)
%             plot(state{ct}(:,2), state{ct}(:,3))
%         end
%     end
% end

% combined plot
ayrange = vertcat(ay_range{:});
ayrange = [min(ayrange(:,1)) max(ayrange(:,2))];
ayrange = [floor(ayrange(1)*10)/10 ceil(max(ayrange(2))*10)/10];
ayrange = [-5 5];

% Choose colormap and find data edges for ay in ayrange
cmap = jet(round(diff(ayrange)/1));
edges = round(linspace(ayrange(1),ayrange(2),size(cmap,1)+1),2); % # edges = # color bins + 1

close all;
figure; hold on;
colormap(cmap);

for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        
        state = state_entry{ct_pattern,ct_light};
        for ct=1:length(state)
            [data_cmap, ~] = discretize(state{ct}(:,4), edges);
            data_cmap(state{ct}(:,4)<=ayrange(1)) = 1;
            data_cmap(state{ct}(:,4)>=ayrange(2)) = size(cmap,1);
            scatter(state{ct}(:,2), state{ct}(:,3), 3, cmap(data_cmap',:),'filled','o');
            
%             plot(state{ct}(:,2), state{ct}(:,3),'b')
        end
    end
end
cmap_bar = colorbar('eastoutside');
caxis(edges([1 end]));
cmap_bar.Ticks = [edges(1) 0 edges(end)];
cmap_bar.Label.String = 'a (ms-2)';
cmap_bar.Label.FontSize = 16;
ylabel('r/r*', 'FontSize', 14);
xlabel('Distance to start of rref (m)', 'FontSize', 14);
set(gca, 'FontSize', 16);


%% Write files for statistical analysis in R
% Create dataset for delta t, mean V and mean ay in 30-90% change in r
clc;
writeFile = true;
if isunix
    r_file = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/rTransientsData_ACC_Rstudio.txt';
elseif ispc
    r_file = 'D:/Disk_08_backup/light_intensity_experiments/postprocessing/rTransientsData_ACC_Rstudio.txt';
end
data_write = [];

nBins = length(rBinValues)-1;

for ct_fac=1:length(factors)
    factor = factors(ct_fac);
    approach_number = 0;
    
    for ct_pattern = 1:length(pattern)
        for ct_light = 1:length(light)
            
            rBins = repmat(1:nBins, size(approachACC{ct_pattern, ct_light, ct_fac},1), 1);
            N = numel(rBins);
            % collecting column-wise i.e., for each rbin
            data_write = [data_write; ... 
                          rBins(:) ct_pattern*ones(N,1) ct_light*ones(N,1) ...
                          risetimesACC{ct_pattern, ct_light, ct_fac}(:) ...
                          distanceACC{ct_pattern, ct_light, ct_fac}(:) ...
                          velocityACC{ct_pattern, ct_light, ct_fac}(:) ...
                          meanVelocityACC{ct_pattern, ct_light, ct_fac}(:) ...
                          meanAccACC{ct_pattern, ct_light, ct_fac}(:) ...
                          meanRdotACC{ct_pattern, ct_light, ct_fac}(:) ...
                          rrefACC{ct_pattern, ct_light, ct_fac}(:) ...
                          approach_number+approachACC{ct_pattern, ct_light, ct_fac}(:) ...
                          landingSideACC{ct_pattern, ct_light, ct_fac}(:) ...
                          timeACC{ct_pattern, ct_light, ct_fac}(:) ...
                          dayACC{ct_pattern, ct_light, ct_fac}(:) ...
                          factor*ones(N,1)];
                      
           approach_number = approach_number + max(approachACC{ct_pattern, ct_light, ct_fac}(:));
        end
    end
end

if writeFile
    T = array2table(data_write, ...
        'VariableNames',{'rbins', 'pattern', 'light', 'delta_t', 'delta_y', 'delta_V', 'mean_V', 'mean_a', 'mean_rdot', 'rref', ...
                         'approach','landingSide','time','day','factor'});
    writetable(T,r_file);
end


%% Loading data and collecting segments with entry dynamics
clc; close all;
% clear;
% 
% inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A2_rrefEntry.mat';
% load(inputFile);
treatments = treatments(1:14*8); % Taking experiments for 2 patterns * 3 lights


pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
behaviour = {'rising','constant','sleeping'};

factors = [0.25:0.25:2.5];
chosen_fac = 1;

data = struct.empty;
% tracks_fac(length(pattern), length(light)) = filteredState_BlindLandingtrack.empty; % tracks for chosen factor

clear dummy;

for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        for ct_behaviour = 2%1:length(behaviour)
            disp(' ');
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
            
            %%%%%%%%%%%%%% Each light and pattern combination %%%%%%%%%%%
            landingTracks = [relevantTreatments.landingTracks];
            state_LDF = [landingTracks.state_LDF];
            rrefEntrySegs = [state_LDF.rrefEntrySegments];
            
            landingTracks_indx4stateLDF = arrayfun(@(x) x*ones(length(landingTracks(x).state_LDF),1),1:length(landingTracks), 'UniformOutput', false);
            landingTracks_indx4stateLDF = vertcat(landingTracks_indx4stateLDF{:});
            
%             % Display info about # of landing tracks
%             disp(['# of distinct Flydra objects (landingTracks): ' num2str(length(landingTracks))]);
%             disp(['# of state LDFs (landingTracks): ' num2str(length(state_LDF))]);
%             disp(['Size of rrefSegment vector: ' num2str(length(rrefSegs))]);
            
            % Finding # of "tracks" (state LDFs) that contain rref entry segments
            % for each factor
            N = zeros(1, length(factors));
            for ct_fac = 1:length(factors)
                factor = factors(ct_fac);
                indices = arrayfun(@(x) ~isempty(x.rrefEntrySegments(abs([x.rrefEntrySegments.factor]-factor)<1e-6).intervals),state_LDF);
                N(ct_fac) = sum(indices);
            end
            
%             % Display 
%             disp(['# of state LDFs containing rref entry segments for different factors: ' num2str(N)]);
            
            %%%%%% Collect state_LDF wise %%%%%%%%%%%%%%%%%%
            % For chosen factor, collect all tracks that do contain
            % non-empty rref entry segments
            indices = arrayfun(@(x) ~isempty(x.rrefEntrySegments(abs([x.rrefEntrySegments.factor]-chosen_fac)<1e-6).intervals), state_LDF);
            dummy.tracks_fac = state_LDF(indices);
            dummy.ct_pattern = ct_pattern;
            dummy.ct_light = ct_light;
            dummy.landingTrack = landingTracks(landingTracks_indx4stateLDF(indices));
            data = [data; dummy];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%% Collect each entry segment wise
%             for ct_treatment=1:length(relevantTreatments) % for each relevant treatment
%                 treatment = relevantTreatments(ct_treatment);
% %                 arrayfun(@(x) x.setLandingSide(),[treatment.landingTracks.state_LDF]'); % to store landing side in the rrefSegments
%                 data1 = arrayfun(@(x) x.rrefEntrySegments,[treatment.landingTracks.state_LDF]','UniformOutput',false);
%                 data1 = horzcat(data1{:});
%                 
%                 % Discard empty intervals
%                 indices = arrayfun(@(x) ~isempty(x.intervals), data1);
%                 
%                 % Save additional data
%                 data1 = data1(indices);
% %                 [data1.pattern] = deal(ct_pattern);
% %                 [data1.light] = deal(ct_light);
% %                 [data1.day] = deal(treatment.datenum);
% %                 [data1.time] = deal(treatment.startTime);
% %     
% %                 data = [data data1];
% % %                 size(data)
% % %                 keyboard;
%             end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end

%% Collect data for sytem identification
clear data4est data4est_unmerged; clc;
data4est = cell(length(data),1); % Merged for sysID
data4est_unmerged = cell(length(data),1); % Unmerged dataset

data4est_unmerged_lowpass = cell(length(data),1); % Unmerged dataset prefiltered with low bandpass filter 
fc = 5; % cut-off frequency in Hz (used for prefiltering the data)

for ct=1:length(data)
    clear estimatedData;
%     estimatedData.tfest = {};
    estimatedData.iddata = {};
    
    estimatedData_prefiltered.iddata = {};

    for ct1=1:length(data(ct).tracks_fac)
        state = data(ct).tracks_fac(ct1).filteredState;
        rrefEntrySegments_fac = data(ct).tracks_fac(ct1).rrefEntrySegments(abs([data(ct).tracks_fac(ct1).rrefEntrySegments.factor]-chosen_fac)<1e-6);
        
        assert(length(rrefEntrySegments_fac) == 1);
        for ct2=1:size(rrefEntrySegments_fac.intervals,1)
            output = -state(rrefEntrySegments_fac.intervals(ct2,1):rrefEntrySegments_fac.intervals(ct2,3),6)./...
                state(rrefEntrySegments_fac.intervals(ct2,1):rrefEntrySegments_fac.intervals(ct2,3),3);
            input = -rrefEntrySegments_fac.rmean(ct2)*ones(length(output),1);
            dt = mean(diff(state(rrefEntrySegments_fac.intervals(ct2,1):rrefEntrySegments_fac.intervals(ct2,3),1)));
            
            [num, den] = butter(2, fc/(1/dt/2),'low');
            filt_output = filtfilt(num, den, output);
            
            data1 = iddata(output, input, dt);
            data1_filtered = iddata(filt_output, input, dt);
            
            if ~isempty(data4est{ct})%exist('data4est','var')
                data4est{ct} = merge(data4est{ct}, data1);
                estimatedData.iddata{end+1} = data1;
                
                estimatedData_prefiltered.iddata{end+1} = data1_filtered;
%                 estimatedData.tfest{end+1} = idtf(zeros(0,0,0));
%                 data4est_unmerged{ct}{end+1} = data1;
            else
                data4est{ct} = data1;
                estimatedData.iddata{1} = data1;
                
                estimatedData_prefiltered.iddata{end+1} = data1_filtered;

%                 estimatedData.tfest{1} = idtf(zeros(0,0,0));
%                 data4est_unmerged{ct}{1} = data1;
            end
        end
    end
    data4est_unmerged{ct} = estimatedData;
    
    data4est_unmerged_lowpass{ct} = estimatedData_prefiltered;
end

% % Testing filter order % %
% fc = 5;
% [num, den] = butter(2, fc/(1/dt/2),'low');
% filt_output2 = filtfilt(num, den, output);
% 
% [num, den] = butter(3, fc/(1/dt/2),'low');
% filt_output3 = filtfilt(num, den, output);
% 
% [num, den] = butter(4, fc/(1/dt/2),'low');
% filt_output4 = filtfilt(num, den, output);
% 
% plot(output); hold on; plot(filt_output2,'r'); plot(filt_output3,'g'); plot(filt_output4,'k');

% % comparing filtered and actual data % % 
close all;
ct = 422;
cond = 6;
plot(data4est_unmerged{cond}.iddata{ct}); hold on;
plot(data4est_unmerged_lowpass{cond}.iddata{ct});

%% create plots for plotting low-pass filtered data with "actual" data

% find indexes of landing tracks for each rrefEntrySegment
clear indexes;
for ct=1:length(data)
    indexes = [];
    for ct1=1:length(data(ct).tracks_fac)
        state = data(ct).tracks_fac(ct1).filteredState;
        rrefEntrySegments_fac = data(ct).tracks_fac(ct1).rrefEntrySegments(abs([data(ct).tracks_fac(ct1).rrefEntrySegments.factor]-chosen_fac)<1e-6);
        
        assert(length(rrefEntrySegments_fac) == 1);
        for ct2=1:size(rrefEntrySegments_fac.intervals,1)
            indexes(end+1) = ct1;
            
            
        end
    end
    data4est_unmerged{ct}.track_indexes = indexes;
end

% create validation plot of system identification for each trackfor ct_fac=10%1:length(factors) % Create plot for each factor

DirPlots = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/plots/BlindLandingTracks/rref_entry_lowpassfiltered_vs_actual';
delPreviousPlots = false; % BE CAREFUL - when set to true, all previously saved plots are deleted
savePlots = false;
savePDFs = false;

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
behaviour = {'rising','constant','sleeping'};


chosen_fac = 1;        
for ct=1:length(data)
    % Delete previous plots
    if delPreviousPlots && savePlots && exist(fullfile(DirPlots, [pattern{data(ct).ct_pattern} '_' light{data(ct).ct_light}]), 'dir')
        rmdir(fullfile(DirPlots, [pattern{data(ct).ct_pattern} '_' light{data(ct).ct_light}]),'s');
    end
    
    % Create sub-directory for each treatment if it doesn't exist
    if savePlots && ~exist(fullfile(DirPlots, [pattern{data(ct).ct_pattern} '_' light{data(ct).ct_light}]), 'dir')
        mkdir(fullfile(DirPlots, [pattern{data(ct).ct_pattern} '_' light{data(ct).ct_light}]));
    end
    DirPlots_treatment = fullfile(DirPlots, [pattern{data(ct).ct_pattern} '_' light{data(ct).ct_light}]);
    
    for ct1=1:length(data(ct).tracks_fac)
        % create tf vector in order of rrefEntrySegments
%         tfs = data4est_unmerged{ct}.tfest(:,:,data4est_unmerged{ct}.track_indexes == ct1);
        
        % plot data
        plotHandles = data(ct).tracks_fac(ct1).plot_rrefsEntry_actual_vs_filtered(chosen_fac, fc);

        if ~isempty(plotHandles)
            
            % Resizing the figures
            for i=1:length(plotHandles)
                plotHandles(i).Position(3) = 680;
                plotHandles(i).Position(4) = 545;
                
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


%% Estimate models 
% Models to be tried
np = 1:3;
nz = 0:3;
nPZ = struc(np, nz);
nPZ(nPZ(:,2)>nPZ(:,1),:) = [];

warning('off','Ident:estimation:transientDataCorrection')

% Estimate individually on all extracted segments 
opt = tfestOptions('InitializeMethod','all','Display','off','SearchMethod','auto');
opt.SearchOptions.MaxIterations = 3000;
for ct=1%:length(data4est_unmerged)
%     data4est_unmerged{ct}.tfest(:,:,size(nPZ,1),length(data4est_unmerged{ct}.iddata)) = idtf(zeros(0,0,0));
    for ct1=1:length(data4est_unmerged{ct}.iddata)
        for ct2=1:size(nPZ,1)
            data4est_unmerged{ct}.tfest(:,:,ct2,ct1) = tfest(data4est_unmerged{ct}.iddata{ct1}, nPZ(ct2,1), nPZ(ct2,2), opt);
        end        
    end
end

figure; hold on;
clear AICc;
for ct=1%:length(data4est_unmerged)
    for ct1=1:35%length(data4est_unmerged{ct}.iddata)
        for ct2=1:size(nPZ,1)
            AICc{ct}(ct2,ct1) = aic(data4est_unmerged{ct}.tfest(:,:,ct2,ct1), 'AICc');
        end
        plot(AICc{ct}(:,ct1))
    end
end

%% Using specific model (np=2, nz=0) - unfiltered data
% Estimation on individual tracks
opt = tfestOptions('InitializeMethod','all','Display','off','SearchMethod','auto');
opt.SearchOptions.MaxIterations = 3000;
warning('off','Ident:estimation:transientDataCorrection')
clear dummy;
tic
for ct=1:length(data4est_unmerged)
    
    parfor ct1=1:length(data4est_unmerged{ct}.iddata)
        
        dummy(:,:,ct1) = tfest(data4est_unmerged{ct}.iddata{ct1}, 2, 0, opt);
             
    end
    data4est_unmerged{ct}.tfest = dummy;
end
toc

%% create plots for estimation on individual tracks
% fit percent plot
labels = cell(0, 1); % for x-axis of boxplots
pattern_label = {'+', 'x'}; % + is checkerboard, x is spokes
light_label = {'L', 'M', 'H'};

fitPercentData = cell(length(data),1);
labels = cell(length(data),1);

for ct=1:length(data)
    
    labels{ct, 1} = [light_label{data(ct).ct_light} '' pattern_label{data(ct).ct_pattern}];
    
    
    fitPercentData{ct, 1} = arrayfun(@(x)  data4est_unmerged{ct}.tfest(:,:,x).Report.Fit.FitPercent,1:length(data4est_unmerged{ct}.tfest));
    
end
fitPercentData_figHandle = createBoxPlot(fitPercentData, labels, 'Fit percent');

%% create plots for estimation on individual tracks

% find indexes of landing tracks for each rrefEntrySegment
clear indexes;
for ct=1:length(data)
    indexes = [];
    for ct1=1:length(data(ct).tracks_fac)
        state = data(ct).tracks_fac(ct1).filteredState;
        rrefEntrySegments_fac = data(ct).tracks_fac(ct1).rrefEntrySegments(abs([data(ct).tracks_fac(ct1).rrefEntrySegments.factor]-chosen_fac)<1e-6);
        
        assert(length(rrefEntrySegments_fac) == 1);
        for ct2=1:size(rrefEntrySegments_fac.intervals,1)
            indexes(end+1) = ct1;
            
            
        end
    end
    data4est_unmerged{ct}.track_indexes = indexes;
end

% create validation plot of system identification for each trackfor ct_fac=10%1:length(factors) % Create plot for each factor

DirPlots = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/plots/BlindLandingTracks/rref_entry_sysIdentified_onIndividualTracks';
delPreviousPlots = false; % BE CAREFUL - when set to true, all previously saved plots are deleted
savePlots = false;
savePDFs = false;

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
behaviour = {'rising','constant','sleeping'};


for ct=1:length(data)
    % Delete previous plots
    if delPreviousPlots && savePlots && exist(fullfile(DirPlots, [pattern{data(ct).ct_pattern} '_' light{data(ct).ct_light}]), 'dir')
        rmdir(fullfile(DirPlots, [pattern{data(ct).ct_pattern} '_' light{data(ct).ct_light}]),'s');
    end
    
    % Create sub-directory for each treatment if it doesn't exist
    if savePlots && ~exist(fullfile(DirPlots, [pattern{data(ct).ct_pattern} '_' light{data(ct).ct_light}]), 'dir')
        mkdir(fullfile(DirPlots, [pattern{data(ct).ct_pattern} '_' light{data(ct).ct_light}]));
    end
    DirPlots_treatment = fullfile(DirPlots, [pattern{data(ct).ct_pattern} '_' light{data(ct).ct_light}]);
    
    for ct1=1:length(data(ct).tracks_fac)
        % create tf vector in order of rrefEntrySegments
        tfs = data4est_unmerged{ct}.tfest(:,:,data4est_unmerged{ct}.track_indexes == ct1);
        
        % plot data
        plotHandles = data(ct).tracks_fac(ct1).plot_rrefsEntry_withSimulatedData(chosen_fac, tfs);

        if ~isempty(plotHandles)
            
            % Resizing the figures
            for i=1:length(plotHandles)
                plotHandles(i).Position(3) = 680;
                plotHandles(i).Position(4) = 545;
                
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

%% create plots for estimation on individual tracks
% Free Parameters box plots
labels = cell(0, 1); % for x-axis of boxplots
pattern_label = {'+', 'x'}; % + is checkerboard, x is spokes
light_label = {'L', 'M', 'H'};

fitParameters = cell(length(data),1);
P1 = cell(length(data),1);
P2 = cell(length(data),1);
P3 = cell(length(data),1);
labels = cell(length(data),1);
clear dummy;
for ct=1:length(data)
    
    labels{ct, 1} = [light_label{data(ct).ct_light} '' pattern_label{data(ct).ct_pattern}];
    
    dummy = arrayfun(@(x)  data4est_unmerged{ct}.tfest(:,:,x).Report.Parameters.ParVector, 1:length(data4est_unmerged{ct}.tfest), 'UniformOutput', false);
    fitParameters{ct, 1} = horzcat(dummy{:})';
    assert(all(fitParameters{ct, 1}(:,4) == 0));
    
    P1{ct,1} = fitParameters{ct, 1}(:,1);
    P2{ct,1} = fitParameters{ct, 1}(:,2);
    P3{ct,1} = fitParameters{ct, 1}(:,3);
end

fitParameter1_figHandle = createBoxPlot(P1, labels, 'P1'); ylim([-Inf 1000])
fitParameter2_figHandle = createBoxPlot(P2, labels, 'P2'); ylim([-Inf 1000])
fitParameter3_figHandle = createBoxPlot(P3, labels, 'P3'); ylim([-Inf 1000])

%% Finding average step response
t = 0:0.005:1;
clear stepResponses meanstepResponse stdStepResponse;
stepResponses = cell(length(data4est_unmerged),1);
meanStepResponse = cell(length(data4est_unmerged),1);
stdStepResponse = cell(length(data4est_unmerged),1);
clear dummy;
for ct=1:length(data4est_unmerged)
    dummy = arrayfun(@(x)  step(data4est_unmerged{ct}.tfest(:,:,x),t), 1:length(data4est_unmerged{ct}.tfest), 'UniformOutput', false);
    stepResponses{ct,1} = horzcat(dummy{:});
    meanStepResponse{ct,1} = mean(stepResponses{ct,1},2);
%     stdStepResponse{ct,1} = std(stepResponses{ct,1},2);
end

stepResponse_fig = figure; hold on;
colors = [120 120 120; 67,162,202; 245 130 46; 120 120 120; 67,162,202; 245 130 46]./255;
for ct=1:length(data)
    plot(t,meanStepResponse{ct,1}, ['-' pattern_label{data(ct).ct_pattern}],'Color',colors(ct,:));
    
%     [y,t] = step(tf_all{ct}(:,:,tf_chosen));
%     
%     plot(t,y, ['-' pattern_label{data(ct).ct_pattern}],'Color',colors(ct,:));
    
%     [y,t] = step(tf_all{ct+3}(:,:,tf_chosen));
% %     figure(fig_spoke);
%     plot(t,y,'x','Color',colors(ct,:));
    
end


%% Collective estimation on all tracks from same treatment
opt = tfestOptions('InitializeMethod','all','Display','off','SearchMethod','auto');
opt.SearchOptions.MaxIterations = 3000;
warning('off','Ident:estimation:transientDataCorrection')
np = 1:3;
nz = 0:3;
nPZ = struc(np, nz);
nPZ(nPZ(:,2)>nPZ(:,1),:) = [];
tf_all = cell(length(data4est),1);
clear dummy;
for ct=1:length(data4est)
    dummy1 = data4est{ct};
    parfor ct1=1:size(nPZ,1)
        dummy(:,:,ct1) = tfest(dummy1, nPZ(ct1,1), nPZ(ct1,2), opt);
    end
    tf_all{ct} = dummy;
end

fig_AICc = figure; hold on;
colors = [120 120 120; 67,162,202; 245 130 46]./255;
markers = {'-+','-x'};
for ct=1:length(data)
    aicc = nan(size(nPZ,1),1);
    for ct1 = 1:size(nPZ,1)
        aicc(ct1,1) = tf_all{ct}(:,:,ct1).Report.Fit.AICc;
    end
%     figure(fig_AICc);
    plot(aicc,markers{data(ct).ct_pattern}, 'Color', colors(data(ct).ct_light,:));
end

% % % together on tracks with same environmental conditions
% % opt = tfestOptions('InitializeMethod','all','Display','off','SearchMethod','auto');
% % opt.SearchOptions.MaxIterations = 3000;
% % warning('off','Ident:estimation:transientDataCorrection')
% % 
% % % collectively for all segments together
% % for ct=1:length(data)
% %     tf_all{ct} = tfest(data4est{ct}, 2, 0, opt);
% % end
%% Selecting specific model (np=2, nz=0) based on AICc variation
tf_chosen=3;

labels = cell(0, 1); % for x-axis of boxplots
pattern_label = {'+', 'x'}; % + is checkerboard, x is spokes
light_label = {'L', 'M', 'H'};

fitPercentData = cell(length(data),1);
labels = cell(length(data),1);

for ct=1:length(data)
    
    labels{ct, 1} = [light_label{data(ct).ct_light} '' pattern_label{data(ct).ct_pattern}];
    
    fitPercentData{ct, 1} = tf_all{ct}(:,:,tf_chosen).Report.Fit.FitPercent;
    
end

fitPercentData_figHandle = createBoxPlot(fitPercentData, labels, 'Fit percent');

stepResponse_fig = figure; hold on;
colors = [120 120 120; 67,162,202; 245 130 46]./255;
for ct=1:length(data)
    figure;
    step(tf_all{ct}(:,:,tf_chosen))
    
%     [y,t] = step(tf_all{ct}(:,:,tf_chosen));
%     
%     plot(t,y, ['-' pattern_label{data(ct).ct_pattern}],'Color',colors(ct,:));
    
%     [y,t] = step(tf_all{ct+3}(:,:,tf_chosen));
% %     figure(fig_spoke);
%     plot(t,y,'x','Color',colors(ct,:));
    
end

%% Create system identification plots for validation
close all;
DirPlots = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/plots/BlindLandingTracks/rref_entry_sysIdentified';
delPreviousPlots = false; % BE CAREFUL - when set to true, all previously saved plots are deleted
savePlots = true;
savePDFs = false;

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
behaviour = {'rising','constant','sleeping'};

factors = [0.25:0.25:2.5];
chosen_fac = [1];
chosen_tf = 6;
for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        for ct_behaviour = 2%1:length(behaviour)
            assert(data((ct_pattern-1)*3+ct_light).ct_pattern == ct_pattern && data((ct_pattern-1)*3+ct_light).ct_light == ct_light);            
            
            % Delete previous plots
            if delPreviousPlots && savePlots && exist(fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour}]), 'dir')
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
            
            
            for ct_treatment=1:length(relevantTreatments)
                treatment = relevantTreatments(ct_treatment);
                
                videoTimes = [[treatment.videosInfo.startTime]' [treatment.videosInfo.endTime]'];
                if ~isempty(treatment.landingTracks)
                    disp(['Into, day: ' num2str(treatment.datenum) ...
                          ', pattern: ' treatment.pattern ...
                          ', light: ' treatment.light ...
                          ', behaviour: ' behaviour{ct_behaviour}]);
                    
                      for ct_track=1:length(treatment.landingTracks) % for each landing track
                          track = treatment.landingTracks(ct_track);
                          for ct_excerpt=1:length(track.state_LDF) % for each track excerpt
                              excerpt = track.state_LDF(ct_excerpt);
                              
                              
                              
                              if savePlots
                                  for ct_factor=1:length(chosen_fac)

                                      plotHandles = excerpt.plot_rrefsEntry_withSimulatedData(chosen_fac, tf_all{(ct_pattern-1)*3+ct_light}(:,:,chosen_tf));
%                                       plotHandles = excerpt.plot_rrefs_with3dspeed(factors(ct_factor));

                                      if ~isempty(plotHandles)

                                          % Resizing the figures
                                          for i=1:length(plotHandles)
                                              plotHandles(i).Position(3) = 680;
                                              plotHandles(i).Position(4) = 545;

                                              if i==1
                                                  figureName = ['fac_' num2str(chosen_fac(ct_factor),'%0.2f') '_' ...
                                                              num2str(treatment.datenum) '_' ...
                                                              num2str(treatment.startTime) '_' num2str(treatment.endTime) ...
                                                              '_obj' num2str(track.obj_id) '_track' ...
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
            end % for treatment
            
        end
    end
end

%% Correlation of increase/decrease in r with increase/decrease in V
data4correlation = cell(length(data),1);
mean_ay = cell(length(data),1);
delta_r = cell(length(data),1);
delta_V = cell(length(data),1);
delta_y = cell(length(data),1);

for ct=1:length(data)
    
    for ct1=1:length(data(ct).tracks_fac)
        state = data(ct).tracks_fac(ct1).filteredState;
        y = -state(:,3); V = state(:,6); ay = state(:,9); r = -state(:,6)./state(:,3);
        
        rrefEntrySegments_fac = data(ct).tracks_fac(ct1).rrefEntrySegments(abs([data(ct).tracks_fac(ct1).rrefEntrySegments.factor]-chosen_fac)<1e-6);
        
        assert(length(rrefEntrySegments_fac) == 1);
        for ct2=1:size(rrefEntrySegments_fac.intervals,1)
            entryStart_indx = rrefEntrySegments_fac.intervals(ct2,1);
            entryEnd_indx = rrefEntrySegments_fac.intervals(ct2,2);
            
            if ~isempty(data4correlation{ct})
                data4correlation{ct}.data{end+1} = [y(entryStart_indx) V(entryStart_indx) r(entryStart_indx);
                                       y(entryEnd_indx)   V(entryEnd_indx)   -rrefEntrySegments_fac.rmean(ct2)];
                mean_ay{ct}(end+1) = mean(ay(entryStart_indx:entryEnd_indx));
                delta_r{ct}(end+1) = -rrefEntrySegments_fac.rmean(ct2) - r(entryStart_indx);
                delta_V{ct}(end+1) = V(entryEnd_indx) - V(entryStart_indx);
                delta_y{ct}(end+1) = y(entryEnd_indx) - y(entryStart_indx);
            else
                data4correlation{ct}.data{1} = [y(entryStart_indx) V(entryStart_indx) r(entryStart_indx);
                                       y(entryEnd_indx)   V(entryEnd_indx)   -rrefEntrySegments_fac.rmean(ct2)];
                mean_ay{ct}(1) = mean(ay(entryStart_indx:entryEnd_indx));
                delta_r{ct}(1) = -rrefEntrySegments_fac.rmean(ct2) - r(entryStart_indx);
                delta_V{ct}(1) = V(entryEnd_indx) - V(entryStart_indx);
                delta_y{ct}(1) = y(entryEnd_indx) - y(entryStart_indx);
            end
        end
    end
end
close all;
figure; hold on;
plot(horzcat(delta_r{:}), horzcat(mean_ay{:}),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
yline(0); xline(0);
ylabel('mean ay during entry (ms-2)', 'FontSize', 16);
xlabel('r* - r1 (s-1)', 'FontSize', 16);
set(gca, 'FontSize', 16); %grid on;

figure; hold on;
plot(horzcat(delta_V{:}), horzcat(mean_ay{:}),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
yline(0); xline(0);
ylabel('mean ay during entry (ms-2)', 'FontSize', 16);
xlabel('V2 - V1 (ms-1)', 'FontSize', 16);
set(gca, 'FontSize', 16); %grid on;

figure; hold on;
plot(horzcat(delta_y{:}), horzcat(delta_r{:}),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
yline(0); xline(0);
ylabel('r* - r1 (s-1)', 'FontSize', 16);
xlabel('y2 - y1 (ms-1)', 'FontSize', 16);
set(gca, 'FontSize', 16); %grid on;

sum(horzcat(mean_ay{:}) > 0 & horzcat(delta_r{:}) > 0)/sum(horzcat(delta_r{:}) > 0)
sum(horzcat(mean_ay{:}) < 0 & horzcat(delta_r{:}) < 0)/sum(horzcat(delta_r{:}) < 0)
%% ONLY PRE FILTERED DATA is used in the code below %%%%%%%
%% Using specific model (np=2, nz=0) - prefiltered data
% Estimation on individual tracks
opt = tfestOptions('InitializeMethod','all','Display','off','SearchMethod','auto');
opt.SearchOptions.MaxIterations = 3000;
warning('off','Ident:estimation:transientDataCorrection')
clear dummy;
tic
for ct=1:length(data4est_unmerged_lowpass)
    
    parfor ct1=1:length(data4est_unmerged_lowpass{ct}.iddata)
        
        dummy(:,:,ct1) = tfest(data4est_unmerged_lowpass{ct}.iddata{ct1}, 2, 0, opt);
             
    end
    data4est_unmerged_lowpass{ct}.tfest = dummy;
end
toc

%% Using specific model (np=2, nz=0) - prefiltered data
% Estimation on individual tracks
opt = tfestOptions('InitializeMethod','all','Display','off','SearchMethod','auto','InitialCondition','estimate');
opt.SearchOptions.MaxIterations = 3000;
warning('off','Ident:estimation:transientDataCorrection')
clear dummy;
tic
for ct=1:length(data4est_unmerged_lowpass)
    
    parfor ct1=1:length(data4est_unmerged_lowpass{ct}.iddata)
        
        dummy(:,:,ct1) = tfest(data4est_unmerged_lowpass{ct}.iddata{ct1}, 2, 0, opt);
             
    end
    data4est_unmerged_lowpass{ct}.tfest2 = dummy;
end
toc

%% create plots for estimation on individual tracks - prefiltered data
% fit percent plot
labels = cell(0, 1); % for x-axis of boxplots
pattern_label = {'+', 'x'}; % + is checkerboard, x is spokes
light_label = {'L', 'M', 'H'};

fitPercentData = cell(length(data),1);
labels = cell(length(data),1);

for ct=1:length(data)
    
    labels{ct, 1} = [light_label{data(ct).ct_light} '' pattern_label{data(ct).ct_pattern}];
    
    
    fitPercentData{ct, 1} = arrayfun(@(x)  data4est_unmerged_lowpass{ct}.tfest(:,:,x).Report.Fit.FitPercent,1:length(data4est_unmerged_lowpass{ct}.tfest));
    
end
fitPercentData_figHandle = createBoxPlot(fitPercentData, labels, 'Fit percent');

%% create plots for estimation on individual tracks - prefiltered data
% Free Parameters box plots
labels = cell(0, 1); % for x-axis of boxplots
pattern_label = {'+', 'x'}; % + is checkerboard, x is spokes
light_label = {'L', 'M', 'H'};

fitParameters = cell(length(data),1);
P1 = cell(length(data),1);
P2 = cell(length(data),1);
P3 = cell(length(data),1);
labels = cell(length(data),1);
clear dummy;
for ct=1:length(data)
    
    labels{ct, 1} = [light_label{data(ct).ct_light} '' pattern_label{data(ct).ct_pattern}];
    
    dummy = arrayfun(@(x)  data4est_unmerged_lowpass{ct}.tfest(:,:,x).Report.Parameters.ParVector, 1:length(data4est_unmerged_lowpass{ct}.tfest), 'UniformOutput', false);
    fitParameters{ct, 1} = horzcat(dummy{:})';
    assert(all(fitParameters{ct, 1}(:,4) == 0));
    
    P1{ct,1} = fitParameters{ct, 1}(:,1);
    P2{ct,1} = fitParameters{ct, 1}(:,2);
    P3{ct,1} = fitParameters{ct, 1}(:,3);
end
close all;
fitParameter1_figHandle = createBoxPlot(P1, labels, 'P1'); ylim([-Inf 1000])
fitParameter2_figHandle = createBoxPlot(P2, labels, 'P2'); ylim([-Inf 60])
fitParameter3_figHandle = createBoxPlot(P3, labels, 'P3'); ylim([-Inf 1000])

%% Estimation based on PRE FILTERED DATA is used in the code below %%%%%%%
%% Estimation on individual tracks based on median values of parameters obtained from previous step
opt = tfestOptions('Display','off');
opt.SearchOptions.MaxIterations = 3000;
warning('off','Ident:estimation:transientDataCorrection')
clear dummy init_sys;

tic
for ct=1:length(data4est_unmerged_lowpass)
    num = median(P1{ct});
    den = [1 median(P2{ct}) median(P3{ct})];
    init_sys(:,:,ct) = idtf(num, den);
    
    parfor ct1=1:length(data4est_unmerged_lowpass{ct}.iddata)
        
        dummy(:,:,ct1) = tfest(data4est_unmerged_lowpass{ct}.iddata{ct1}, init_sys(:,:,ct));
             
    end
    data4est_unmerged_lowpass{ct}.tfest_refined = dummy;
end
toc

%% create plots for estimation on individual tracks - prefiltered data
% fit percent plot
labels = cell(0, 1); % for x-axis of boxplots
pattern_label = {'+', 'x'}; % + is checkerboard, x is spokes
light_label = {'L', 'M', 'H'};

fitPercentData = cell(length(data),1);
labels = cell(length(data),1);

for ct=1:length(data)
    
    labels{ct, 1} = [light_label{data(ct).ct_light} '' pattern_label{data(ct).ct_pattern}];
    
    
    fitPercentData{ct, 1} = arrayfun(@(x)  data4est_unmerged_lowpass{ct}.tfest_refined(:,:,x).Report.Fit.FitPercent,1:length(data4est_unmerged_lowpass{ct}.tfest));
    
end
fitPercentData_figHandle = createBoxPlot(fitPercentData, labels, 'Fit percent');

%% create plots for estimation on individual tracks - prefiltered data
% Free Parameters box plots
labels = cell(0, 1); % for x-axis of boxplots
pattern_label = {'+', 'x'}; % + is checkerboard, x is spokes
light_label = {'L', 'M', 'H'};

fitParameters = cell(length(data),1);
P1_refined = cell(length(data),1);
P2_refined = cell(length(data),1);
P3_refined = cell(length(data),1);
labels = cell(length(data),1);
clear dummy;
for ct=1:length(data)
    
    labels{ct, 1} = [light_label{data(ct).ct_light} '' pattern_label{data(ct).ct_pattern}];
    
    dummy = arrayfun(@(x)  data4est_unmerged_lowpass{ct}.tfest_refined(:,:,x).Report.Parameters.ParVector, 1:length(data4est_unmerged_lowpass{ct}.tfest), 'UniformOutput', false);
    fitParameters{ct, 1} = horzcat(dummy{:})';
    assert(all(fitParameters{ct, 1}(:,4) == 0));
    
    P1_refined{ct,1} = fitParameters{ct, 1}(:,1);
    P2_refined{ct,1} = fitParameters{ct, 1}(:,2);
    P3_refined{ct,1} = fitParameters{ct, 1}(:,3);
end
% close all;
fitParameter1_figHandle = createBoxPlot(P1_refined, labels, 'P1'); ylim([-Inf 1000])
fitParameter2_figHandle = createBoxPlot(P2_refined, labels, 'P2'); ylim([-Inf 60])
fitParameter3_figHandle = createBoxPlot(P3_refined, labels, 'P3'); ylim([-Inf 1000])

%% Plot "rise" times in different environmental conditions

clc; close all;
% clear;
% 
% inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A2_rrefEntry.mat';
inputFile = 'D:/light_intensity_experiments/postprocessing/BlindLandingtracks_A2_rrefEntry.mat';
% load(inputFile);


treatments = treatments(1:14*8); % Taking experiments for 2 patterns * 3 lights


pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
behaviour = {'rising','constant','sleeping'};

factors = [0.25:0.25:2.5];

data_all = struct.empty;

rPercentIntervalsACC = 30:60:90; % in percent of r*
rPercentIntervalsDEC = 190:-10:110; % in percent of r*

rIntervals = 0.5:0.5:10; % in terms of r

% Performance parameters
risetimesACC = cell.empty;
distanceACC = cell.empty;
velocityACC = cell.empty;
deltarACC = cell.empty;
accACC = cell.empty;
meanDistanceACC = cell.empty;
meanVelocityACC = cell.empty;
meanAccACC = cell.empty;
meanRdotACC = cell.empty;
rrefACC = cell.empty;
approachACC = cell.empty;
landingSideACC = cell.empty;
timeACC = cell.empty;
dayACC = cell.empty;
r0ACC = cell.empty; % r at rPercentIntervalsACC(1)
rdot0ACC = cell.empty; % rdot at rPercentIntervalsACC(1)
rdot0_filteredACC = cell.empty; % rdot at rPercentIntervalsACC(1)


risetimesDEC = cell.empty;
distanceDEC = cell.empty;
velocityDEC = cell.empty;
deltarDEC = cell.empty;
accDEC = cell.empty;
meanDistanceDEC = cell.empty;
meanVelocityDEC = cell.empty;
meanAccDEC = cell.empty;
meanRdotDEC = cell.empty;
rrefDEC = cell.empty;
approachDEC = cell.empty;
landingSideDEC = cell.empty;
timeDEC = cell.empty;
dayDEC = cell.empty;

entryData = cell.empty;
fc= 5; % Hz

for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        for ct_behaviour = 2%1:length(behaviour)
            disp(' ');
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
            
            %%%%%%%%%%%%%% Each light and pattern combination %%%%%%%%%%%
            landingTracks = [relevantTreatments.landingTracks];
            state_LDF = [landingTracks.state_LDF];
            rrefEntrySegs = [state_LDF.rrefEntrySegments];
            
            landingTracks_indx4stateLDF = arrayfun(@(x) x*ones(length(landingTracks(x).state_LDF),1),1:length(landingTracks), 'UniformOutput', false);
            landingTracks_indx4stateLDF = vertcat(landingTracks_indx4stateLDF{:});
            
            % Finding # of "tracks" (state LDFs) that contain rref entry segments
            % for each factor
            for ct_fac = 1:length(factors)
                factor = factors(ct_fac);
                indices = arrayfun(@(x) ~isempty(x.rrefEntrySegments(abs([x.rrefEntrySegments.factor]-factor)<1e-6).intervals),state_LDF);
                
                data_all(ct_pattern, ct_light, ct_fac).tracks_fac = state_LDF(indices);
                data_all(ct_pattern, ct_light, ct_fac).ct_pattern = ct_pattern;
                data_all(ct_pattern, ct_light, ct_fac).ct_light = ct_light;
                data_all(ct_pattern, ct_light, ct_fac).ct_factor = factor;
                data_all(ct_pattern, ct_light, ct_fac).landingTrack = landingTracks(landingTracks_indx4stateLDF(indices));
                
                risetimesACC{ct_pattern, ct_light, ct_fac} = [];
                distanceACC{ct_pattern, ct_light, ct_fac} = [];
                velocityACC{ct_pattern, ct_light, ct_fac} = [];
                deltarACC{ct_pattern, ct_light, ct_fac} = [];
                accACC{ct_pattern, ct_light, ct_fac} = [];
                meanDistanceACC{ct_pattern, ct_light, ct_fac} = [];
                meanVelocityACC{ct_pattern, ct_light, ct_fac} = [];
                meanAccACC{ct_pattern, ct_light, ct_fac} = [];
                meanRdotACC{ct_pattern, ct_light, ct_fac} = [];
                rrefACC{ct_pattern, ct_light, ct_fac} = [];
                approachACC{ct_pattern, ct_light, ct_fac} = [];
                landingSideACC{ct_pattern, ct_light, ct_fac} = [];
                timeACC{ct_pattern, ct_light, ct_fac} = [];
                dayACC{ct_pattern, ct_light, ct_fac} = [];
                r0ACC{ct_pattern, ct_light, ct_fac} = [];
                rdot0ACC{ct_pattern, ct_light, ct_fac} = [];
                rdot0_filteredACC{ct_pattern, ct_light, ct_fac} = [];

                risetimesDEC{ct_pattern, ct_light, ct_fac} = [];
                distanceDEC{ct_pattern, ct_light, ct_fac} = [];
                velocityDEC{ct_pattern, ct_light, ct_fac} = [];
                deltarDEC{ct_pattern, ct_light, ct_fac} = [];
                accDEC{ct_pattern, ct_light, ct_fac} = [];
                meanDistanceDEC{ct_pattern, ct_light, ct_fac} = [];
                meanVelocityDEC{ct_pattern, ct_light, ct_fac} = [];
                meanAccDEC{ct_pattern, ct_light, ct_fac} = [];
                meanRdotDEC{ct_pattern, ct_light, ct_fac} = [];
                rrefDEC{ct_pattern, ct_light, ct_fac} = [];
                approachDEC{ct_pattern, ct_light, ct_fac} = [];
                landingSideDEC{ct_pattern, ct_light, ct_fac} = [];
                timeDEC{ct_pattern, ct_light, ct_fac} = [];
                dayDEC{ct_pattern, ct_light, ct_fac} = [];

                entryData{ct_fac} = [];
                
                for ct1=1:sum(indices) % for each state_LDF
                    state = data_all(ct_pattern, ct_light, ct_fac).tracks_fac(ct1).filteredState;
                    y = -state(:,3); V = state(:,6); ay = state(:,9); r = -state(:,6)./state(:,3); t = state(:,1)-state(1,1); 
                    rdot = diffxy(t, r); dt = t(2)-t(1);
                    
                    rrefEntrySegments_fac = data_all(ct_pattern, ct_light, ct_fac).tracks_fac(ct1).rrefEntrySegments(abs([data_all(ct_pattern, ct_light, ct_fac).tracks_fac(ct1).rrefEntrySegments.factor]-factor)<1e-6);
                    
                    for ct2=1:size(rrefEntrySegments_fac.intervals,1)
                        entryStart_indx = rrefEntrySegments_fac.intervals(ct2,1);
                        entryEnd_indx = rrefEntrySegments_fac.intervals(ct2,2);
                        steadyEnd_indx = rrefEntrySegments_fac.intervals(ct2,3);
                        
                        rref = -rrefEntrySegments_fac.rmean(ct2);
                        
                        % make entryEnd_indx to be one after 90%rref or
                        % until r is increasing
                        if r(entryEnd_indx) < 0.9*rref
%                             entryEnd_indx = entryEnd_indx + find(r(entryEnd_indx:steadyEnd_indx)>0.9*rref,1) - 1;
                            entryEnd_indx = entryEnd_indx + find(diff(r(entryEnd_indx:steadyEnd_indx))<0,1) - 1;
                        end
                        
                        y_part = y(entryStart_indx:entryEnd_indx); V_part = V(entryStart_indx:entryEnd_indx); ay_part = ay(entryStart_indx:entryEnd_indx);
                        rdot_part = rdot(entryStart_indx:entryEnd_indx);
%                         [num, den] = butter(2, fc/(1/dt/2),'low');
%                         r_filt = filtfilt(num, den, r(entryStart_indx:steadyEnd_indx)); rdot_filt = diffxy([0:1:length(r_filt)-1]'*dt, r_filt);
                        
                        entryData{ct_fac} = [entryData{ct_fac}; ...
                                              t(entryStart_indx:entryEnd_indx)  y(entryStart_indx:entryEnd_indx) ...
                                              V(entryStart_indx:entryEnd_indx)  ay(entryStart_indx:entryEnd_indx) ...
                                              r(entryStart_indx:entryEnd_indx)  rdot(entryStart_indx:entryEnd_indx) ...
                                              rref*ones(entryEnd_indx-entryStart_indx+1,1)];
                        
                        if r(entryStart_indx) > rref % populate risetimesDEC
                            rBinValues = rPercentIntervalsDEC*rref/100;
%                             rBinValues = rIntervals;
                            t_at_rBinValues = interp1(r(entryStart_indx:entryEnd_indx), t(entryStart_indx:entryEnd_indx), rBinValues, 'makima', nan);
                            risetimesDEC{ct_pattern, ct_light, ct_fac}(end+1,:) = diff(t_at_rBinValues);
                            
                            %                     if any(diff(t_at_rBinValues) < 0)
                            %                         keyboard;
                            %                     end
                            
                        elseif r(entryStart_indx) < rref % populate risetimesACC
                            rBinValues = rPercentIntervalsACC*rref/100;
%                             rBinValues = rIntervals;
                            
                            t_at_rBinValues = interp1(r(entryStart_indx:entryEnd_indx), t(entryStart_indx:entryEnd_indx), rBinValues, 'makima', nan);
                            y_at_rBinValues = interp1(r(entryStart_indx:entryEnd_indx), y(entryStart_indx:entryEnd_indx), rBinValues, 'makima', nan);
                            V_at_rBinValues = interp1(r(entryStart_indx:entryEnd_indx), V(entryStart_indx:entryEnd_indx), rBinValues, 'makima', nan);
                            a_at_rBinValues = interp1(r(entryStart_indx:entryEnd_indx), ay(entryStart_indx:entryEnd_indx), rBinValues, 'makima', nan);
                            rdot_at_rBinValues = interp1(t(entryStart_indx:entryEnd_indx), rdot_part,  t_at_rBinValues, 'makima', nan);
                            
                            if any(diff(t_at_rBinValues) < 0) %isnan(diff(t_at_rBinValues))
                                keyboard;
                            end
                            risetimesACC{ct_pattern, ct_light, ct_fac}(end+1,:) = diff(t_at_rBinValues);
                            distanceACC{ct_pattern, ct_light, ct_fac}(end+1,:) = diff(y_at_rBinValues);
                            velocityACC{ct_pattern, ct_light, ct_fac}(end+1,:) = diff(V_at_rBinValues);
                            deltarACC{ct_pattern, ct_light, ct_fac}(end+1,:) = diff(rBinValues) + diff(t_at_rBinValues) - diff(t_at_rBinValues);
                            accACC{ct_pattern, ct_light, ct_fac}(end+1,:) = diff(a_at_rBinValues);
                            
                            meanDistanceACC{ct_pattern, ct_light, ct_fac}(end+1,:) = arrayfun(@(x) mean(y_part(r(entryStart_indx:entryEnd_indx) > rBinValues(x) & ...
                                                                                               r(entryStart_indx:entryEnd_indx) <= rBinValues(x+1))),1:length(rBinValues)-1) ...
                                                                                     + diff(t_at_rBinValues) - diff(t_at_rBinValues); % + diff(t_at_rBinValues) - diff(t_at_rBinValues) is used to make values nan where diff(t_at_rBinValues) is also nan 
                            meanVelocityACC{ct_pattern, ct_light, ct_fac}(end+1,:) = arrayfun(@(x) mean(V_part(r(entryStart_indx:entryEnd_indx) > rBinValues(x) & ...
                                                                                               r(entryStart_indx:entryEnd_indx) <= rBinValues(x+1))),1:length(rBinValues)-1) ...
                                                                                     + diff(t_at_rBinValues) - diff(t_at_rBinValues); % + diff(t_at_rBinValues) - diff(t_at_rBinValues) is used to make values nan where diff(t_at_rBinValues) is also nan 
                            meanAccACC{ct_pattern, ct_light, ct_fac}(end+1,:) = arrayfun(@(x) mean(ay_part(r(entryStart_indx:entryEnd_indx) > rBinValues(x) & ...
                                                                                               r(entryStart_indx:entryEnd_indx) <= rBinValues(x+1))),1:length(rBinValues)-1) ...
                                                                                     + diff(t_at_rBinValues) - diff(t_at_rBinValues);
                            meanRdotACC{ct_pattern, ct_light, ct_fac}(end+1,:) = arrayfun(@(x) mean(rdot_part(r(entryStart_indx:entryEnd_indx) > rBinValues(x) & ...
                                                                                               r(entryStart_indx:entryEnd_indx) <= rBinValues(x+1))),1:length(rBinValues)-1) ...
                                                                                     + diff(t_at_rBinValues) - diff(t_at_rBinValues);
                            rrefACC{ct_pattern, ct_light, ct_fac}(end+1,:) = rref*ones(1,length(rBinValues)-1);
                            
                            % at the beginning of complete rise segment
%                             r0ACC{ct_pattern, ct_light, ct_fac}(end+1,:) = r(entryStart_indx)*ones(1,length(rBinValues)-1);
%                             rdot0ACC{ct_pattern, ct_light, ct_fac}(end+1,:) = rdot(entryStart_indx)*ones(1,length(rBinValues)-1);
                            
                            % at the beginning of each bin segment
                            r0ACC{ct_pattern, ct_light, ct_fac}(end+1,:) = rBinValues(1:end-1);
                            rdot0ACC{ct_pattern, ct_light, ct_fac}(end+1,:) = rdot_at_rBinValues(1:end-1);
                            
                            approachACC{ct_pattern, ct_light, ct_fac}(end+1,:) = ct1*ones(1,length(rBinValues)-1);
                            if strcmpi(data_all(ct_pattern, ct_light, ct_fac).tracks_fac(ct1).landingSide,"Hive")
                                landingSideACC{ct_pattern, ct_light, ct_fac}(end+1,:) = 1*ones(1,length(rBinValues)-1);
                            elseif strcmpi(data_all(ct_pattern, ct_light, ct_fac).tracks_fac(ct1).landingSide,"Feeder")
                                landingSideACC{ct_pattern, ct_light, ct_fac}(end+1,:) = 2*ones(1,length(rBinValues)-1);
                            end
                            timeACC{ct_pattern, ct_light, ct_fac}(end+1,:) = str2double(data_all(ct_pattern, ct_light, ct_fac).landingTrack(ct1).foldername(10:13))*ones(1,length(rBinValues)-1);
                            dayACC{ct_pattern, ct_light, ct_fac}(end+1,:) = data_all(ct_pattern, ct_light, ct_fac).landingTrack(ct1).datenum*ones(1,length(rBinValues)-1);
                           
                        end
                        
                    end
                    
                end

            end
        end
    end
end
%% Write files for statistical analysis in R
% writing percentage change in r file
% required columns =
clc;
writeFile = true;
removeNAN = true;
% r_file = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/absolute_rTransientsData_ACC_Rstudio.txt';
r_file = 'D:/light_intensity_experiments/postprocessing/absolute_rTransientsData_ACC_Rstudio.txt';
data_write = [];

nBins = length(rBinValues)-1;

for ct_fac=1:length(factors)
    factor = factors(ct_fac);
    approach_number = 0;
    
    for ct_pattern = 1:length(pattern)
        for ct_light = 1:length(light)
            
            rBins = repmat(1:nBins, size(approachACC{ct_pattern, ct_light, ct_fac},1), 1);
            N = numel(rBins);
            % collecting column-wise i.e., for each rbin
            data_write = [data_write; ... 
                          rBins(:) ct_pattern*ones(N,1) ct_light*ones(N,1) ...
                          risetimesACC{ct_pattern, ct_light, ct_fac}(:) ...
                          distanceACC{ct_pattern, ct_light, ct_fac}(:) ...
                          velocityACC{ct_pattern, ct_light, ct_fac}(:) ...
                          deltarACC{ct_pattern, ct_light, ct_fac}(:) ...
                          meanDistanceACC{ct_pattern, ct_light, ct_fac}(:) ...
                          meanVelocityACC{ct_pattern, ct_light, ct_fac}(:) ...
                          meanAccACC{ct_pattern, ct_light, ct_fac}(:) ...
                          meanRdotACC{ct_pattern, ct_light, ct_fac}(:) ...
                          rrefACC{ct_pattern, ct_light, ct_fac}(:) ...
                          r0ACC{ct_pattern, ct_light, ct_fac}(:) ...
                          rdot0ACC{ct_pattern, ct_light, ct_fac}(:) ...
                          approach_number+approachACC{ct_pattern, ct_light, ct_fac}(:) ...
                          landingSideACC{ct_pattern, ct_light, ct_fac}(:) ...
                          timeACC{ct_pattern, ct_light, ct_fac}(:) ...
                          dayACC{ct_pattern, ct_light, ct_fac}(:) ...
                          factor*ones(N,1)];
                      
           approach_number = approach_number + max(approachACC{ct_pattern, ct_light, ct_fac}(:));
        end
    end
end
if removeNAN
    data_write(isnan(data_write(:,4)) | data_write(:,14)<0,:) = [];
end

if writeFile
    T = array2table(data_write, ...
        'VariableNames',{'rbins', 'pattern', 'light', 'delta_t', 'delta_y', 'delta_V', 'delta_r', 'mean_y', 'mean_V', 'mean_a', 'mean_rdot', 'rref', 'r0', 'rdot0', ...
                         'approach','landingSide','time','day','factor'});
    writetable(T,r_file);
end

%% Plot rise time as a function of rdot and rref 
close all;
% For 30-45%


% For 45-60%
data_ss = data_write(data_write(:,1)==3 & abs(data_write(:,end) - 1.75) < 1e-6 & data_write(:,18) <= 20190710,:);
rref = data_ss(:,12); rdot0 = data_ss(:,14); delta_t = data_ss(:,4); light = data_ss(:,3);

% Untransformed domain
y = delta_t; x1 = rref; x2 = rdot0;
% plot3(log(rref), log(rdot0), log(delta_t),'.')
yrange = [min(y) max(y)];
% yrange = [floor(yrange(1)*1000)/1000 ceil(max(yrange(2))*1000)/1000];
yrange = [round(yrange(1)*100)/100 round(max(yrange(2))*100)/100];

% Choose colormap and find data edges for V in Vrange
cmap = jet(round(diff(yrange)/0.01));
edges = round(linspace(yrange(1),yrange(2),size(cmap,1)+1),2); % # edges = # color bins + 1

figure2d = figure;hold on;
colormap(cmap);
for ct=1:3 % for each light condition
    subplot(1,3,ct);
    y_ss = y(light==ct);
    [data_cmap, ~] = discretize(y_ss, edges);
    data_cmap(y_ss<=edges(1)) = 1; data_cmap(y_ss>=edges(end)) = size(cmap,1);
    scatter(x1(light==ct), x2(light==ct), 10, cmap(data_cmap',:),'filled','o');
    set(gca, 'FontSize', 16);
end
figure(figure2d);
cmap_bar = colorbar('eastoutside');
caxis(edges([1 end]));
cmap_bar.Ticks = [edges(1) edges(end)];
cmap_bar.Label.String = 'delta_t';
cmap_bar.Label.FontSize = 16;
ylabel('rdot0 (s-2)', 'FontSize', 14);
xlabel('rref (s-1)', 'FontSize', 14);

%% % %  log domain
close all;
chosen_fac = 1.5;
for rbin = 4%2:5
    data_ss = data_write(data_write(:,1)==rbin & abs(data_write(:,end) - chosen_fac) < 1e-6 & data_write(:,18) <= 20190710,:);
    rref = data_ss(:,12); rdot0 = data_ss(:,14); delta_t = data_ss(:,4); light = data_ss(:,3);
    y = log(delta_t); x1 = log(rref); x2 = log(rdot0);
    
    yrange = [min(y) max(y)];
    yrange = [round(yrange(1)*10)/10 round(max(yrange(2))*10)/10];

    % Choose colormap and find data edges for V in Vrange
    cmap = jet(round(diff(yrange)/0.1));
    edges = round(linspace(yrange(1),yrange(2),size(cmap,1)+1),2); % # edges = # color bins + 1

    figure2d = figure; hold on;
    colormap(cmap);
    
    % Fit from R
    if rbin==2
        y_est = @(x1,x2,light) -2.23580-0.78793*x2+0.48181*x1+0.14966*x1*x2+light; light_intercepts = [0 -0.05197 -0.09037];
    elseif rbin==3
        y_est = @(x1,x2,light) -2.57597-0.61084*x2+0.73729*x1+light; light_intercepts = [0 -0.04894 -0.06576];
    elseif rbin==4
        y_est = @(x1,x2,light) -2.52836-0.61336*x2+0.68441*x1+light; light_intercepts = [0 -0.01388 -0.04033];
    elseif rbin==5
        y_est = @(x1,x2,light) -2.62646-0.54172*x2+0.80827*x1-0.06759*x1*x2+light; light_intercepts = [0 -0.04770 -0.06759];
    end

    for ct=1:3 % for each light condition
        subplot(1,3,ct); hold on;
        caxis(edges([1 end])); fcontour(@(x1,x2) y_est(x1,x2,light_intercepts(ct)), [floor(min(x1)) ceil(max(x1)*10)/10 0 ceil(max(x2)*10)/10], 'Fill', 'on', 'LevelList', (edges(1:end-1)+edges(2:end))/2);

        y_ss = y(light==ct);
        [data_cmap, ~] = discretize(y_ss, edges);
        data_cmap(y_ss<=edges(1)) = 1; data_cmap(y_ss>=edges(end)) = size(cmap,1);
        scatter(x1(light==ct), x2(light==ct), 20, cmap(data_cmap',:),'filled','o','MarkerEdgeColor',[0 0 0]);
        set(gca, 'FontSize', 16);
        xlim([floor(min(x1)) ceil(max(x1)*10)/10]);
    %     ylim([floor(min(x2)) ceil(max(x2)*10)/10]);
        ylim([0 ceil(max(x2)*10)/10]);
        cmap_bar = colorbar('northoutside');
        caxis(edges([1 end]));
        cmap_bar.Ticks = [edges(1) edges(end)];
        cmap_bar.Label.String = 'log(delta t)';
        cmap_bar.Label.FontSize = 16;

    end
    figure(figure2d);
    subplot(1,3,1);
    ylabel('log(rdot0) (s-2)', 'FontSize', 14);
    title('low light', 'FontSize', 14);
    subplot(1,3,2);
    xlabel('log(r*) (s-1)', 'FontSize', 14);
    title('medium light', 'FontSize', 14);
    subplot(1,3,3);
    title('high light', 'FontSize', 14);

    
    figure3d = figure; hold on;
    colormap(cmap);
    % Fit from R
    y_est = @(x1,x2,light) -2.55985-0.61423*x2+0.72759*x1+light; light_intercepts = [0 -0.04157 -0.06486];
    x1_grid = floor(min(x1)):0.1:ceil(max(x1)*10)/10;
    x2_grid = 0:0.1:ceil(max(x2)*10)/10;
    [X1, X2] = meshgrid(x1_grid, x2_grid);
    for ct=1:3 % for each light condition
        subplot(1,3,ct); hold on;
        caxis(edges([1 end])); surf(X1, X2, y_est(X1, X2,light_intercepts(ct)), 'EdgeColor', 'none');
        
        
        y_ss = y(light==ct);
        [data_cmap, ~] = discretize(y_ss, edges);
        data_cmap(y_ss<=edges(1)) = 1; data_cmap(y_ss>=edges(end)) = size(cmap,1);
        scatter3(x1(light==ct), x2(light==ct), y_ss, 20, cmap(data_cmap',:),'filled','o','MarkerEdgeColor',[0 0 0]);
        set(gca, 'FontSize', 16);
        xlim([floor(min(x1)) ceil(max(x1)*10)/10]);
        %     ylim([floor(min(x2)) ceil(max(x2)*10)/10]);
        ylim([0 ceil(max(x2)*10)/10]);
        cmap_bar = colorbar('northoutside');
        caxis(edges([1 end]));
        cmap_bar.Ticks = [edges(1) edges(end)];
        cmap_bar.Label.String = 'log(delta t)';
        cmap_bar.Label.FontSize = 16;
        
    end
    figure(figure3d);
    subplot(1,3,1);
    zlabel('log(delta t) (s-2)', 'FontSize', 14);
    ylabel('log(rdot0) (s-2)', 'FontSize', 14);
    xlabel('log(r*) (s-1)', 'FontSize', 14);
    title('low light', 'FontSize', 14);
    subplot(1,3,2);
    zlabel('log(delta t) (s-2)', 'FontSize', 14);
    ylabel('log(rdot0) (s-2)', 'FontSize', 14);
    xlabel('log(r*) (s-1)', 'FontSize', 14);
    title('medium light', 'FontSize', 14);
    subplot(1,3,3);
    zlabel('log(delta t) (s-2)', 'FontSize', 14);
    ylabel('log(rdot0) (s-2)', 'FontSize', 14);
    xlabel('log(r*) (s-1)', 'FontSize', 14);
    title('high light', 'FontSize', 14);
    
    keyboard;

end

figure; colormap(cmap); caxis(edges([1 end])); fcontour(y_est, [floor(min(x1)) ceil(max(x1)*10)/10 0 ceil(max(x2)*10)/10], 'Fill', 'on', 'LevelList', (edges(1:end-1)+edges(2:end))/2); colorbar('northoutside');
figure; contourf(X1, X2, Y, edges); colormap(cmap);

Y = y_est(X1, X2);
figure;  view(2)
for ct=1:3 % for each light condition
    
end

%% Plot rise time vs meany vs meanr
mean_y = data_write(:,8);
rref = data_write(:,12);
riseTimes = data_write(:,4);
light = data_write(:,3); f = data_write(:,end);
figure;
% scatter3(mean_y,rref,riseTimes,...
%         'MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 .75 .75])
scatter3(log(mean_y),log(rref),log(riseTimes),4,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75])
view(-30,10);
xlabel('Mean y (m)', 'FontSize', 16);
ylabel('r* (1/s)', 'FontSize', 16);
zlabel('delta t (s)', 'FontSize', 16);
set(gca, 'FontSize', 18); grid on;
view([0 0]);
view([90 0]);
view([0 90]);
% xlim([0 Inf]); ylim([0 Inf]);
% writing absolute change in r file
%%
% close all;
meanRisetimes = cellfun(@(x) nanmean(x) ,risetimesACC,'UniformOutput' ,false);
medianRisetimes = cellfun(@(x) nanmedian(x) ,risetimesACC,'UniformOutput' ,false);
nValues = cellfun(@(x) sum(~isnan(x)) ,risetimesACC,'UniformOutput' ,false);

nx = size(risetimesACC{1,1,1},2);
ny = 3;
cmap = [120 120 120; 67,162,202; 245 130 46]./255;

for ct_fac=4%1:length(factors) % Create plot for each factor
    
    figure;
    subplot(2,1,1);
    ndata = max([size(risetimesACC{1,1,ct_fac},1) size(risetimesACC{1,2,ct_fac},1) size(risetimesACC{1,3,ct_fac},1)]);
    y=cellfun(@(x) [x;nan(ndata-size(x,1),nx)],risetimesACC(1,1:3,ct_fac),'un',0);
%     y=cellfun(@(x) [x;nan(ndata-size(x,1),nx)],meanVelocityACC(1,1:3,ct_fac),'un',0);
%     y=cellfun(@(x) [x;nan(ndata-size(x,1),nx)],distanceACC(1,1:3,ct_fac),'un',0);
    y = cat(3,y{:}); % Concatenating along third dimension
    
%     y = permute(y, [2 3 1]);
%     h1 = boxplot2(y); 
%     for ii = 1:ny
%         structfun(@(x) set(x(ii,:), 'color', cmap(ii,:), ...
%             'markeredgecolor', cmap(ii,:)), h1);
%     end
%     set([h1.lwhis h1.uwhis], 'linestyle', '-');
%     set(h1.out, 'marker', '.');
    
    y = permute(y, [1 2 3 ]);
% 	xlabels = arrayfun(@(x) [num2str(rIntervals(x),'%0.1f') '-' num2str(rIntervals(x+1),'%0.1f')], 1:nx,'un',0);
	xlabels = arrayfun(@(x) [num2str(rPercentIntervalsACC(x),'%2.0f') '-' num2str(rPercentIntervalsACC(x+1),'%2.0f') '%r*'], 1:nx,'un',0);
    h1 = iosr.statistics.boxPlot(y, 'sampleSize', true, 'showMean', true, 'symbolMarker', '.','x', xlabels);
    h1.lineColor{1} = cmap(1,:); h1.lineColor{2} = cmap(2,:); h1.lineColor{3} = cmap(3,:); 
    h1.meanColor{1} = cmap(1,:); h1.meanColor{2} = cmap(2,:); h1.meanColor{3} = cmap(3,:); 
    h1.medianColor{1} = cmap(1,:); h1.medianColor{2} = cmap(2,:); h1.medianColor{3} = cmap(3,:); 
    h1.symbolColor{1} = cmap(1,:); h1.symbolColor{2} = cmap(2,:); h1.symbolColor{3} = cmap(3,:); 
    h1.GroupLabels = {{'L', 'M', 'H'}};
    h1.showLegend = true;
    title('checkerboard');
    ylabel('delta t');
%     ylabel('mean V');

%     set(gca, 'FontSize', 18);
%     h1.x = ;
%     h1.boxColor{1} = cmap(1,:); h1.boxColor{2} = cmap(2,:); h1.boxColor{3} = cmap(3,:); 

    subplot(2,1,2);
    ndata = max([size(risetimesACC{2,1,ct_fac},1) size(risetimesACC{2,2,ct_fac},1) size(risetimesACC{2,3,ct_fac},1)]);
    y=cellfun(@(x) [x;nan(ndata-size(x,1),nx)],risetimesACC(2,1:3,ct_fac),'un',0);
%     y=cellfun(@(x) [x;nan(ndata-size(x,1),nx)],meanVelocityACC(2,1:3,ct_fac),'un',0);
%     y=cellfun(@(x) [x;nan(ndata-size(x,1),nx)],distanceACC(2,1:3,ct_fac),'un',0);
    y = cat(3,y{:}); % Concatenating along third dimension
    
    y = permute(y, [1 2 3 ]);
    h1 = iosr.statistics.boxPlot(y, 'sampleSize', true, 'showMean', true, 'symbolMarker', '.','x', xlabels);
    h1.lineColor{1} = cmap(1,:); h1.lineColor{2} = cmap(2,:); h1.lineColor{3} = cmap(3,:); 
    h1.meanColor{1} = cmap(1,:); h1.meanColor{2} = cmap(2,:); h1.meanColor{3} = cmap(3,:); 
    h1.medianColor{1} = cmap(1,:); h1.medianColor{2} = cmap(2,:); h1.medianColor{3} = cmap(3,:); 
    h1.symbolColor{1} = cmap(1,:); h1.symbolColor{2} = cmap(2,:); h1.symbolColor{3} = cmap(3,:); 
    ylabel('delta t');
%     ylabel('mean V');
%     set(gca, 'FontSize', 18);
    title('spoke');
end
    
%% plot rref, r and rdot 3D-plot
for ct_fac=4%1:length(factors) % Create plot for each factor
    plot3(entryData{ct_fac}(:,7), entryData{ct_fac}(:,5), entryData{ct_fac}(:,6),'.')
    xlabel('r*', 'FontSize', 16);
    ylabel('r', 'FontSize', 16);
    zlabel('rdot', 'FontSize', 16);
    set(gca, 'FontSize', 18); grid on;
    
    
    figure;
    plot(entryData{ct_fac}(:,5), entryData{ct_fac}(:,3),'.');
end

%%
close all;
meanRisetimes = cellfun(@(x) nanmean(x) ,risetimesACC,'UniformOutput' ,false);
medianRisetimes = cellfun(@(x) nanmedian(x) ,risetimesACC,'UniformOutput' ,false);
nValues = cellfun(@(x) sum(~isnan(x)) ,risetimesACC,'UniformOutput' ,false);

% meanPlot = figure;
% medianPlot = figure;

for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        n = size(risetimesACC{1,1,1},2);
        meanrisetimes = vertcat(meanRisetimes{ct_pattern, ct_light, :});
        medianrisetimes = vertcat(medianRisetimes{ct_pattern, ct_light, :});
        
        figure;
        
        for ct=1:n
            subplot(n,1,ct);
            plot(factors, meanrisetimes(:,ct),'Linewidth',2);
%             plot(factors, medianrisetimes(:,ct),'Linewidth',2);
            ylabel([num2str(rPercentIntervalsACC(ct)) '-' num2str(rPercentIntervalsACC(ct+1)) '% r*']);
            ylabel([num2str(rIntervals(ct)) '-' num2str(rIntervals(ct+1)) '% r*']);
            xlim([factors(1) factors(end)]);
            if ct==1
                title(['Pattern: ' pattern{ct_pattern} ', light: ' light{ct_light}], 'FontSize', 16);
            end
        end
        xlabel('factor f', 'FontSize', 16);
%         set(gca, 'FontSize', 16);
        
    end
end

meanRisetimes = vertcat(meanRisetimes{:});


medianRisetimes = vertcat(medianRisetimes{:});


n = vertcat(n{:});



for ct=1:length(data)
    
    risetimesDEC{ct, length(factors)} = [];
    risetimesACC{ct, length(factors)} = [];
    
    for ct1=1:length(data(ct).tracks_fac)
        state = data(ct).tracks_fac(ct1).filteredState;
        y = -state(:,3); V = state(:,6); ay = state(:,9); r = -state(:,6)./state(:,3); t = state(:,1)-state(1,1);
        
        for ct_factor=1:length(factors)
            fac = factors(ct_factor);
        
            rrefEntrySegments_fac = data(ct).tracks_fac(ct1).rrefEntrySegments(abs([data(ct).tracks_fac(ct1).rrefEntrySegments.factor]-fac)<1e-6);
        
            assert(length(rrefEntrySegments_fac) == 1);
            for ct2=1:size(rrefEntrySegments_fac.intervals,1)
                entryStart_indx = rrefEntrySegments_fac.intervals(ct2,1);
                entryEnd_indx = rrefEntrySegments_fac.intervals(ct2,2);

                rref = -rrefEntrySegments_fac.rmean(ct2);
            
            
                if r(entryStart_indx) > rref % populate risetimesDEC
                    rBinValues = rPercentIntervalsDEC*rref/100;
                    t_at_rBinValues = interp1(r(entryStart_indx:entryEnd_indx), t(entryStart_indx:entryEnd_indx), rBinValues, 'makima', nan);
                    risetimesDEC{ct, ct_factor}(end+1,:) = diff(t_at_rBinValues);

%                     if any(diff(t_at_rBinValues) < 0)
%                         keyboard;
%                     end

                elseif r(entryStart_indx) < rref % populate risetimesACC
                    rBinValues = rPercentIntervalsACC*rref/100;
                    t_at_rBinValues = interp1(r(entryStart_indx:entryEnd_indx), t(entryStart_indx:entryEnd_indx), rBinValues, 'makima', nan);
                    risetimesACC{ct, ct_factor}(end+1,:) = diff(t_at_rBinValues);

%                     if any(diff(t_at_rBinValues) < 0)
%                         keyboard;
%                     end

                end
            
            end
        end
    end
end
keyboard
meanRisetimes = cellfun(@(x) nanmean(x) ,risetimesACC,'UniformOutput' ,false);
meanRisetimes = vertcat(meanRisetimes{:});

medianRisetimes = cellfun(@(x) nanmedian(x) ,risetimesACC,'UniformOutput' ,false);
medianRisetimes = vertcat(medianRisetimes{:});

n = cellfun(@(x) sum(~isnan(x)) ,risetimesACC,'UniformOutput' ,false);
n = vertcat(n{:});

figure; histogram(risetimesACC{4}(:,8));
figure; histogram(risetimesACC{5}(:,8));
figure; histogram(risetimesACC{6}(:,8));


close all;
figure; boxplot2(risetimesACC{1}, 10:10:80); ylim([0 0.1]);
figure; boxplot2(risetimesACC{2}, [10:10:80] + 5);  ylim([0 0.1]);
figure; boxplot2(risetimesACC{3}, [10:10:80] + 10);  ylim([0 0.1]);

LowLight_xvalues_acc = 14:10:84;
LowLight_xvalues_dec = 114:10:194;

for ct_pattern = 1:2
    indexes = find([data.ct_pattern] == ct_pattern);
    
    risetimesACC_plot = risetimesACC{indexes};
    xvalues = [LowLight_xvalues_acc LowLight_xvalues_acc+1 LowLight_xvalues_acc+2];
    
%     risetimesDEC_plot = risetimesDEC{indexes};
    
    
end
%%
fig_checkerboard = figure; hold on;
fig_spoke = figure; hold on;
colors = [120 120 120; 67,162,202; 245 130 46]./255;
P = pzoptions;
for ct=1:3
    figure(fig_checkerboard);
    iopzplot(tf_all{ct});
    
    figure(fig_spoke);
    iopzplot(tf_all{ct+3});
end
figure(fig_checkerboard);
legend('L','M','H');

figure;
h = iopzplot(tf_all{ct});
compare(data4est_unmerged{3}.iddata{280},tf_all{3})
compare(data4est_unmerged{3}.iddata{280},data4est_unmerged{3}.tfest(:,:,280))
compare(data4est{1},tf_all{1})
resid(data4est,tf_all)

for ct=1:length(estimatedData.iddata)
    figure;
    resid(estimatedData.iddata{ct},tf_all)
end


% h = iopzplot(tf_all,pzoptions('ConfidenceRegionNumberSD',95));
% showConfidence(h);

% Create extended data
output = [track.track_subset_sysID(ct_excerpt).state(:,6)./track.track_subset_sysID(ct_excerpt).state(:,3)];
output = [output; output(end)*ones(20,1)];
input = track.track_subset_sysID(ct_excerpt).rref*ones(length(output),1);
input = [input; input(end)*ones(20,1)];
dummy = iddata(output,input,dt);
compare(dummy,tf_all)


% if exist('data4est','var')
%     data4est = merge(data4est, data1);
%     estimatedData.iddata{end+1} = data1;
%     %                                       estimatedData.tfest(:,:,end+1) = tfest(data1, 3, 3, Options);
% else
%     data4est = data1;
%     estimatedData.iddata{1} = data1;
%     %                                       estimatedData.tfest(:,:,1) = tfest(data1, 3, 3, Options);
% end

%%
opt = tfestOptions('InitializeMethod','all','Display','on','SearchMethod','auto');
opt.SearchOptions.MaxIterations = 3000;
warning('off','Ident:estimation:transientDataCorrection')
a = tfest(data4est{2}, 2, 0, opt)

%% Functions used in this script
function figHandle = createBoxPlot(variable, labels, yxislabel)
    % for boxplots of a cell array
    col=@(x)reshape(x,numel(x),1);
    boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});

    
    figHandle = figure; hold on;
    figure(figHandle);
    boxplot2(variable, 'Labels', labels, 'OutlierSize', 0.00001);
%     hbox = gca;
    set(gca, 'FontSize', 18); %grid on;
%     xlabel('Treatments', 'FontSize', 18);
    ylabel(yxislabel, 'FontSize', 18);
    % ylim([0 1.5]);
    for ct = 1:length(variable)
            x=(ct+(rand(length(variable{ct, 1}),1)-0.5)/4);

            f = scatter(x(:,1),variable{ct, 1},10,'k','filled'); 
            f.MarkerFaceAlpha = 0.5;
    %         keyboard;
    end
end

















