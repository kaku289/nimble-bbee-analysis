%% Loading data and creating plots for Introduction figure
% Plots are plotted for each track and highlights rref segment(s), entry
% segment(s), filtered data, and data estimated from system identification

% clc; clear; close all;
% clc; 
% clear all; close all;

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

%%

if isunix
%     dataFile = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/BlindLandingtracks_A2_rrefEntryEstimation_everything_f1.mat';
%     dataFile = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/BlindLandingtracks_A2_rrefEntryEstimation_everything_f1pt5_all.mat';
    DirPlots = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/plots/BlindLandingTracks/articleA3_plots';
elseif ispc
%     dataFile = 'D:/steady_wind_experiments/postprocessing/BlindLandingtracks_A2_rrefEntryEstimation_everything_f1.mat';
%     dataFile = 'D:/steady_wind_experiments/postprocessing/BlindLandingtracks_A2_rrefEntryEstimation_everything_f1pt5_all.mat';
    DirPlots = 'D:/steady_wind_experiments/postprocessing/plots/BlindLandingTracks/articleA2_plots';
end
% load(dataFile);

% For f=1 (based on track numbers)
delPreviousPlots = true; % BE CAREFUL - when set to true, all previously saved plots are deleted
savePlots = true;
savePDFs = true;

winds = unique([treatments.wind]);
behaviour = {'rising','constant','sleeping'};

factors = [0.25:0.25:2.5];
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

            % 1. For plotting lowpass filtered signal and coresponding output
            % from sysID
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

%% For f=1 (based on track identity - more general/easier than above)

if isunix
%     dataFile = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/BlindLandingtracks_A2_rrefEntryEstimation_everything_f1.mat';
%     dataFile = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/BlindLandingtracks_A2_rrefEntryEstimation_everything_f1pt5_all.mat';
    DirPlots = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/plots/BlindLandingTracks/articleA3_plots';
elseif ispc
%     dataFile = 'D:/steady_wind_experiments/postprocessing/BlindLandingtracks_A2_rrefEntryEstimation_everything_f1.mat';
%     dataFile = 'D:/steady_wind_experiments/postprocessing/BlindLandingtracks_A2_rrefEntryEstimation_everything_f1pt5_all.mat';
    DirPlots = 'D:/steady_wind_experiments/postprocessing/plots/BlindLandingTracks/articleA3_plots';
end
% load(dataFile);

delPreviousPlots = true; % BE CAREFUL - when set to true, all previously saved plots are deleted
savePlots = true;
savePDFs = true;

winds = unique([treatments.wind]);
behaviour = {'rising','constant','sleeping'};

factors = [0.25:0.25:2.5];
chosen_fac = 1;

% behaviour = {'rising','constant','sleeping'};
selected_tracks = {'20190722_123000_140000_obj5712_trackXXX_excerpt1'};

assert(length(pattern)==size(data_all,1) && length(light)==size(data_all,2));
for ct_wind=6%1:length(winds)
        % Delete previous plots
        if delPreviousPlots && savePlots && exist(fullfile(DirPlots, [num2str(winds(ct_wind)) '_' behaviour{ct_behaviour}]), 'dir')
            rmdir(fullfile(DirPlots, [num2str(winds(ct_wind)) '_' behaviour{ct_behaviour}]),'s');
        end
        
        % Create sub-directory for each treatment if it doesn't exist
        if savePlots && ~exist(fullfile(DirPlots, [num2str(winds(ct_wind)) '_' behaviour{ct_behaviour}]), 'dir')
            mkdir(fullfile(DirPlots, [num2str(winds(ct_wind)) '_' behaviour{ct_behaviour}]));
        end
        DirPlots_treatment = fullfile(DirPlots, [num2str(winds(ct_wind)) '_' behaviour{ct_behaviour}]);
        
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

            % 1. For plotting lowpass filtered signal and coresponding output
            % from sysID
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

%% To generate Fitpercent plot for sysID model output
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

clear dummy;
for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        
        %Compute y_start for each rref entry segment
        ystart_indx = arrayfun(@(x) x.rrefEntrySegments([x.rrefEntrySegments.factor]==1).intervals(:,1) , data_all(ct_pattern,ct_light).tracks_fac,'UniformOutput',false);
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

%
gain_figHandle = createBoxPlot({gain{:}}', {labels{:}}', 'Gain'); ylim([0 2])
omega_figHandle = createBoxPlot({omega{:}}', {labels{:}}', 'Natural frequency'); ylim([-Inf 40])
zeta_figHandle = createBoxPlot({zeta{:}}', {labels{:}}', 'Damping'); ylim([-Inf 0.1])


%%
% First extract rdot estimate for chosen factor using
% A2_AA_extract_rrefEntrySegments_with_rdot_ystart.m
clc;
if isunix
    DirPlots = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/plots/BlindLandingTracks/articleA3_plots';
elseif ispc
    DirPlots = 'D:/steady_wind_experiments/postprocessing/plots/BlindLandingTracks/articleA3_plots';
end

delPreviousPlots = true; % BE CAREFUL - when set to true, all previously saved plots are deleted
savePlots = true;
savePDFs = true;

winds = unique([treatments.wind]);
behaviour = {'rising','constant','sleeping'};

factors = [0.25:0.25:2.5];
chosen_fac = 1;

% behaviour = {'rising','constant','sleeping'};
selected_tracks = {'20190716_110000_123000_obj4691_track130_excerpt1',...
    '20190717_123000_140000_obj688_track18_excerpt3',...
    '20190722_123000_140000_obj5712_trackXXX_excerpt1', ...
    '20190722_123000_140000_obj3201_trackXXX_excerpt2',...
    '20190719_080000_093000_obj6358_track151_excerpt1',...
    '20190718_153000_170000_obj4121_track145_excerpt1'};

plotHandles = [];
for ct_wind=6%1:length(winds)
    
        % Delete previous plots
        if delPreviousPlots && savePlots && exist(fullfile(DirPlots, [num2str(winds(ct_wind)) '_' behaviour{ct_behaviour}]), 'dir')
            rmdir(fullfile(DirPlots, [num2str(winds(ct_wind)) '_' behaviour{ct_behaviour}]),'s');
        end
        
        % Create sub-directory for each treatment if it doesn't exist
        if savePlots && ~exist(fullfile(DirPlots, [num2str(winds(ct_wind)) '_' behaviour{ct_behaviour}]), 'dir')
            mkdir(fullfile(DirPlots, [num2str(winds(ct_wind)) '_' behaviour{ct_behaviour}]));
        end
        DirPlots_treatment = fullfile(DirPlots, [num2str(winds(ct_wind)) '_' behaviour{ct_behaviour}]);
        
        
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
        
        treatmentIndices = cell(length(relevantTreatments),1);
        for ct_treatment=1:length(relevantTreatments) % for each relevant treatment
            treatment = relevantTreatments(ct_treatment);
            dummy = ct_treatment*ones(length([treatment.landingTracks]),1);
%             dummy = arrayfun(@(x,indx) indx*ones(length(x.state_LDF),1), [treatment.landingTracks], treatmentIndices4landingTracks,'UniformOutput',false);
            treatmentIndices{ct_treatment} = dummy;
                
                
%                 hastakeoff_pertreatment = arrayfun(@(x) x.hasTakeoff(treatment.landingDiscs),[treatment.landingTracks.state_LDF]','UniformOutput',false);
% %                 hastakeoff_pertreatment = horzcat(hastakeoff_pertreatment{:});
%                 hastakeoff{ct_treatment} = hastakeoff_pertreatment;
        end
%         hastakeoff = vertcat(hastakeoff{:})
        treatmentIndices = vertcat(treatmentIndices{:});  
        
        landingTracks = [relevantTreatments.landingTracks];
        
        % Gather details of all landing tracks
        datenum = [landingTracks.datenum];
        obj_id = [landingTracks.obj_id];
        
        for ct2=1:length(selected_tracks)
            selected_track_str = selected_tracks{ct2};
            selected_track_str_parts = strsplit(selected_track_str,'_');
            
            indx = find(datenum == str2num(selected_track_str_parts{1}) & obj_id == str2num(selected_track_str_parts{4}(4:end)));
%             assert(length(indx)==1); % If not, tighten the indx selection with starting time of the treatment
            
            ct1 = indx(1)           
           
            excerpt = str2num(selected_track_str_parts{end}(end));
            hasTakeoff = landingTracks(ct1).state_LDF(excerpt).hasTakeoff(relevantTreatments(treatmentIndices(ct1)).landingDiscs)

            % For plotting rdotestimate 
%             plotHandles = landingTracks(ct1).state_LDF(excerpt).plot_rrefsEntry_with_rdotestimate(chosen_fac)
            
            % For plotting rref curves (same as Article 1) - size (680 545)
%             plotHandles = landingTracks(ct1).state_LDF(excerpt).plot_rrefs(chosen_fac);
            

            % For plotting V vs t curve
%             plotHandles = landingTracks(ct1).state_LDF(excerpt).plot_rrefsEntry_for_Vvst(chosen_fac)
            
            % For plotting V,r,A vs y curve
%             plotHandles = landingTracks(ct1).state_LDF(excerpt).plot_rrefsEntry_with_acc(chosen_fac)
            
            
            if ~isempty(plotHandles)
                
                % Resizing the figures
                for i=1:length(plotHandles)
                    plotHandles(i).Position(3) = 560;
                    plotHandles(i).Position(4) = 560;
                    
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