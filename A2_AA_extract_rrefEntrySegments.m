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

%% Find r* entry intervals
% Also find videos asociated with each track

% Inputs
close all; clc;
clear;

inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A1_rref.mat';
outputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A2_rrefEntry.mat';
load(inputFile);

DirPlots = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/plots/BlindLandingTracks/rref_entry';
% DirPlots = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/plots/BlindLandingTracks/rref_estimate_3dspeed';
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

                                      plotHandles = excerpt.plot_rrefsEntry(factors(ct_factor));
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
            
            % Finding # of "tracks" (state LDFs) that contain rref segments
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
for ct=1:length(data)
    clear estimatedData;
%     estimatedData.tfest = {};
    estimatedData.iddata = {};

    for ct1=1:length(data(ct).tracks_fac)
        state = data(ct).tracks_fac(ct1).filteredState;
        rrefEntrySegments_fac = data(ct).tracks_fac(ct1).rrefEntrySegments(abs([data(ct).tracks_fac(ct1).rrefEntrySegments.factor]-chosen_fac)<1e-6);
        
        assert(length(rrefEntrySegments_fac) == 1);
        for ct2=1:size(rrefEntrySegments_fac.intervals,1)
            output = -state(rrefEntrySegments_fac.intervals(ct2,1):rrefEntrySegments_fac.intervals(ct2,3),6)./...
                state(rrefEntrySegments_fac.intervals(ct2,1):rrefEntrySegments_fac.intervals(ct2,3),3);
            input = -rrefEntrySegments_fac.rmean(ct2)*ones(length(output),1);
            dt = mean(diff(state(rrefEntrySegments_fac.intervals(ct2,1):rrefEntrySegments_fac.intervals(ct2,3),1)));
            
            data1 = iddata(output, input, dt);
            if ~isempty(data4est{ct})%exist('data4est','var')
                data4est{ct} = merge(data4est{ct}, data1);
                estimatedData.iddata{end+1} = data1;
%                 estimatedData.tfest{end+1} = idtf(zeros(0,0,0));
%                 data4est_unmerged{ct}{end+1} = data1;
            else
                data4est{ct} = data1;
                estimatedData.iddata{1} = data1;
%                 estimatedData.tfest{1} = idtf(zeros(0,0,0));
%                 data4est_unmerged{ct}{1} = data1;
            end
        end
    end
    data4est_unmerged{ct} = estimatedData;
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

%% Using specific model (np=2, nz=0)
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

% create validation plot of system identification for each track
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

fitParameter1_figHandle = createBoxPlot(P1, labels, 'P1');
fitParameter2_figHandle = createBoxPlot(P2, labels, 'P2');
fitParameter3_figHandle = createBoxPlot(P3, labels, 'P3');

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

















