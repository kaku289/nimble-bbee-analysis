%% Importing previously created libraries

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


%% Find r* intervals
% Also find videos asociated with each track

% Inputs
close all; clc;
% clear;

if isunix
    inputFile = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/BlindLandingtracks_A3_LDF.mat';
    outputFile = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/BlindLandingtracks_A3_LDF_rref_sideAccData.mat';
    DirPlots = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/plots/BlindLandingTracks/rref_estimate_sideAcc';
elseif ispc
    inputFile = 'G:/steady_wind_experiments/postprocessing/BlindLandingtracks_A3_LDF.mat';
    outputFile = 'G:/steady_wind_experiments/postprocessing/BlindLandingtracks_A3_LDF_rref_sideAccData.mat';
    DirPlots = 'G:/steady_wind_experiments/postprocessing/plots/BlindLandingTracks/rref_estimate_sideAcc';
end

delPreviousPlots = false; % BE CAREFUL - when set to true, all previously saved plots are deleted
savePlots = false;
savePDFs = false;

% load(inputFile);
winds = unique([treatments.wind]);
behaviour = {'rising','constant','sleeping'};

rmse_func = @(idx, indices, data) sqrt(sum((data(indices(idx,1):indices(idx,2)) - mean(data(indices(idx,1):indices(idx,2)))).^2)/(indices(idx,2)-indices(idx,1)+1));

min_gap = 14; % the first straight line is between current point and (current point + min_gap)th point
max_gap = 49;
params = [0.53 4.22]; % [sigma_{rmean-c_{r vs y, 0-1}}, sigma_{m_{r vs y, 0-1}}]
factors = [0.25:0.25:2.5];
time_window = [min_gap max_gap];
% factors = [1];
ct_data = 0;
for ct_wind = 1:length(winds)
    
    for ct_behaviour = 2%1:length(behaviour)
        
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
        
        
        for ct_treatment=1:length(relevantTreatments)
            treatment = relevantTreatments(ct_treatment);
            
            videoTimes = [[treatment.videosInfo.startTime]' [treatment.videosInfo.endTime]'];
            if ~isempty(treatment.landingTracks)
                disp(['Into, day: ' num2str(treatment.datenum) ...
                    ', wind: ' num2str(treatment.wind) ...
                    ', behaviour: ' behaviour{ct_behaviour}]);
                
                for ct_track=1:length(treatment.landingTracks) % for each landing track
                    track = treatment.landingTracks(ct_track);
                    for ct_excerpt=1:length(track.state_LDF) % for each track excerpt
                        excerpt = track.state_LDF(ct_excerpt);
                        
                        excerpt.compute_rref(params, factors, time_window);
                        %                               excerpt.compute_rref_with3dspeed(params, factors, time_window);
                        
                        if savePlots
                            for ct_factor=1:length(factors)
                                
                                plotHandles = excerpt.plot_rrefs(factors(ct_factor));
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
keyboard
save(outputFile, 'treatments', '-v7.3', '-nocompression');
keyboard;

%% Extract rref data for statistical analysis
clc; close all;
% clear;

if isunix
    inputFile = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/BlindLandingtracks_A3_LDF_rref_sideAccData.mat';
elseif ispc
    inputFile = 'G:/steady_wind_experiments/postprocessing/BlindLandingtracks_A3_LDF_rref_sideAccData.mat';
end

% load(inputFile);

winds = unique([treatments.wind]);
behaviour = {'rising','constant','sleeping'};

data = data4rrefEstimate.empty;
hasTakeoff = []; % whether or not the landing track has takeoff dynamics in it as defined in filteredState_BlindLandingtrack.hasTakeoff(..) method
for ct_wind = 1:length(winds)
    for ct_behaviour = 2%1:length(behaviour)
        
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
        
        hasUniformHwData = arrayfun(@(x) x.hwData.hasUniformHwData,relevantTreatments);
        relevantTreatments = relevantTreatments(hasUniformHwData);
        
        
        for ct_treatment=1:length(relevantTreatments) % for each relevant treatment
            treatment = relevantTreatments(ct_treatment);
            
            arrayfun(@(x) x.setLandingSide(),[treatment.landingTracks.state_LDF]'); % to store landing side in the rrefSegments
            arrayfun(@(x) x.compute_params_basedon_3dspeed(), [treatment.landingTracks.state_LDF]');
            
            data1 = arrayfun(@(x) x.rrefSegments,[treatment.landingTracks.state_LDF]','UniformOutput',false);
            data1 = horzcat(data1{:});
            
            % compute whether or not instability follows
            arrayfun(@(x) x.compute_instabilityFollows2(),[treatment.landingTracks.state_LDF]);

            hastakeoff_pertreatment = arrayfun(@(x) x.hasTakeoff(treatment.landingDiscs)*ones(1,length(x.rrefSegments)),[treatment.landingTracks.state_LDF]','UniformOutput',false);
            hastakeoff_pertreatment = horzcat(hastakeoff_pertreatment{:});
            % Discard empty intervals
            indices = arrayfun(@(x) ~isempty(x.intervals_ti), data1);
            
            % Save additional data
            data1 = data1(indices);
%             [data1.pattern] = deal(ct_pattern);
%             [data1.light] = deal(ct_light);
            [data1.wind] = deal(winds(ct_wind));
            [data1.day] = deal(treatment.datenum);
            [data1.time] = deal(treatment.startTime);
            
            data = [data data1];
            hasTakeoff = [hasTakeoff hastakeoff_pertreatment(indices)];
            %                 size(data)
            %                 keyboard;
        end
    end
end

%% Write file for statistical analysis in R
% This file contains data for all the factors
clc;
writeFile = true;

if isunix
    r_file = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/data_all_rref_sideAcc_Rstudio.txt';
elseif ispc
    r_file = 'G:/steady_wind_experiments/postprocessing/data_all_rref_sideAcc_Rstudio.txt';
end

factors = [0.25:0.25:2.5];
data_write = [];
for ct_factor=1:length(factors)
    factor = factors(ct_factor);
    
    data_fac = data(abs([data.factor]-factor)<1e-6)';
    hasTakeoff_fac = hasTakeoff(abs([data.factor]-factor)<1e-6)';
    N = length(data_fac);
    N_fromtakeoff = sum(hasTakeoff_fac);
    % Create nominal and ordinal variables
    approach = arrayfun(@(i) i*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    side = arrayfun(@(i) data_fac(i).side*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    wind = arrayfun(@(i) data_fac(i).wind*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
%     pattern = arrayfun(@(i) data_fac(i).pattern*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
%     light = arrayfun(@(i) data_fac(i).light*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    time = arrayfun(@(i) data_fac(i).time*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    day = arrayfun(@(i) data_fac(i).day*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    
    hasTakeoff_fac = arrayfun(@(i) hasTakeoff_fac(i)*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    
    % create other variables
    y = -vertcat(data_fac.ymean_ti);
    r = -vertcat(data_fac.rref_ti);
    v = vertcat(data_fac.vmean_ti);
    delta_y = vertcat(data_fac.yrange);
    delta_t = vertcat(data_fac.fd_actual_ti);
    hasInstabilityAssociated = vertcat(data_fac.instabilityFollows);
    y_rrefEnd = vertcat(data_fac.y_rrefEnd);
    instability_deltat = vertcat(data_fac.instability_deltat);
    instability_meanv = vertcat(data_fac.instability_meanv);
    
    
    speed3d = vertcat(data_fac.speed3d_mean_ti);
%     rSpeed3d = -vertcat(data_fac.rmean_speed3d_ti);
    
%     logy = log(y);
%     logr = log(V);
    
    rrefStartState = vertcat(data_fac.rrefStartState);
    rrefMeanState = vertcat(data_fac.rrefMeanState);
    rrefDeltaState = vertcat(data_fac.rrefDeltaState);

    data_write = [data_write; ...
        vertcat(approach{:}) vertcat(side{:}) vertcat(wind{:}) ...
        vertcat(time{:}) vertcat(day{:}) y r v factor*ones(size(r,1),1) speed3d vertcat(hasTakeoff_fac{:}) ...
        delta_y delta_t hasInstabilityAssociated y_rrefEnd ...
        instability_deltat instability_meanv rrefStartState rrefMeanState rrefDeltaState(:,1:6)];
    
    N1 = sum(arrayfun(@(x) size(x.intervals_ti,1)>1, data_fac));
    disp(['# tracks: ' num2str(N) ', # data points: ' num2str(size(r,1)), ...
          ', # tracks with >1 r*: ' num2str(N1)]);
    disp(['# tracks from takeoff: ' num2str(N_fromtakeoff)]);
end

if writeFile
    T = array2table(data_write, ...
        'VariableNames',{'approach','landingSide','wind','time','day','y',...
        'r','v','threshold', 'speed3d', 'hasTakeoff','delta_y','delta_t','hia','y_rrefEnd',...
        'i_deltat','i_meanv','x0','y0','z0','u0','v0','w0','ax0','ay0','az0',...
        'xMean','yMean','zMean','uMean','vMean','wMean','axMean','ayMean','azMean',...
        'deltaX','deltaY','deltaZ','deltaU','deltaV','deltaW'});
    writetable(T,r_file);
end


%% Find r* entry intervals
% Also find videos asociated with each track

% Inputs
close all; clc;
% clear;

if isunix
    inputFile = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/BlindLandingtracks_A3_LDF_rref_sideAccData.mat';
    outputFile = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/BlindLandingtracks_A3_LDF_rref_rrefEntry_with_rdot_sideAccData.mat';
    DirPlots = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/plots/BlindLandingTracks/rref_entry_sideAcc';
elseif ispc
    inputFile = 'G:/steady_wind_experiments/postprocessing/BlindLandingtracks_A3_LDF_rref_sideAccData.mat';
    outputFile = 'G:/steady_wind_experiments/postprocessing/BlindLandingtracks_A3_LDF_rref_rrefEntry_with_rdot_sideAccData.mat';
    DirPlots = 'G:/steady_wind_experiments/postprocessing/plots/BlindLandingTracks/rref_entry_sideAcc';
end

% load(inputFile);

delPreviousPlots = false; % BE CAREFUL - when set to true, all previously saved plots are deleted
savePlots = false;
savePDFs = false;

winds = unique([treatments.wind]);
behaviour = {'rising','constant','sleeping'};

factors = [0.25:0.25:2.5];
% factors = [1];

winds_vel = [0, 0.28, 0.98, 1.77, 2.54, 3.41];
% interval_for_rdot_estimate = [0.2 0.8]; % in percentage of rref
for ct_wind = 1:length(winds)
    for ct_behaviour = 2%1:length(behaviour)
        
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
        
        hasUniformHwData = arrayfun(@(x) x.hwData.hasUniformHwData,relevantTreatments);
        relevantTreatments = relevantTreatments(hasUniformHwData);
        
        
        for ct_treatment=1:length(relevantTreatments)
            treatment = relevantTreatments(ct_treatment);
            
%             videoTimes = [[treatment.videosInfo.startTime]' [treatment.videosInfo.endTime]'];
            if ~isempty(treatment.landingTracks)
                disp(['Into, day: ' num2str(treatment.datenum) ...
                    ', wind: ' num2str(treatment.wind) ...
                    ', pattern: ' treatment.pattern ...
                    ', behaviour: ' behaviour{ct_behaviour}]);
                
                for ct_track=1:length(treatment.landingTracks) % for each landing track
                    track = treatment.landingTracks(ct_track);
                    for ct_excerpt=1:length(track.state_LDF) % for each track excerpt
                        excerpt = track.state_LDF(ct_excerpt);
                        
                        excerpt.find_rrefEntry();
                        excerpt.find_rdot_estimate_in_rrefEntry([winds_vel(winds(ct_wind)) 0 0]);
                              
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
%                                       plotHandles(i).Position(3) = 680;
%                                       plotHandles(i).Position(4) = 545;
                                        
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
keyboard
save(outputFile, 'treatments', '-v7.3', '-nocompression');
keyboard;

%% Extract rrefEntry data for statistical analysis
clc; close all;

winds = unique([treatments.wind]);
behaviour = {'rising','constant','sleeping'};

data = data4rrefEntry.empty;
hasTakeoff = []; % whether or not the landing track has takeoff dynamics in it as defined in filteredState_BlindLandingtrack.hasTakeoff(..) method

for ct_wind = 1:length(winds)
    for ct_behaviour = 2%1:length(behaviour)
        
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
        
        hasUniformHwData = arrayfun(@(x) x.hwData.hasUniformHwData,relevantTreatments);
        relevantTreatments = relevantTreatments(hasUniformHwData);
        
        
        for ct_treatment=1:length(relevantTreatments) % for each relevant treatment
            treatment = relevantTreatments(ct_treatment);
            
            arrayfun(@(x) x.setLandingSide(),[treatment.landingTracks.state_LDF]'); % to store landing side in the rrefSegments
            
            data1 = arrayfun(@(x) x.rrefEntrySegments,[treatment.landingTracks.state_LDF]','UniformOutput',false);
            data1 = horzcat(data1{:});
            
            hastakeoff_pertreatment = arrayfun(@(x) x.hasTakeoff(treatment.landingDiscs)*ones(1,length(x.rrefEntrySegments)),[treatment.landingTracks.state_LDF]','UniformOutput',false);
            hastakeoff_pertreatment = horzcat(hastakeoff_pertreatment{:});
            
            % Discard empty intervals
            indices = arrayfun(@(x) ~isempty(x.intervals), data1);
            
            % Save additional data
            data1 = data1(indices);
            [data1.wind] = deal(ct_wind);
            [data1.day] = deal(treatment.datenum);
            [data1.time] = deal(treatment.startTime);
            
            data = [data data1];
            hasTakeoff = [hasTakeoff hastakeoff_pertreatment(indices)];
            %                 size(data)
            %                 keyboard;
        end
    end
end

%% Write file for statistical analysis in R
% This file contains data for all the factors
clc;
writeFile = true;
if isunix
    r_file = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/data_all_rdot_sideAcc_Rstudio.txt';
elseif ispc
    r_file = 'G:/steady_wind_experiments/postprocessing/data_all_rdot_sideAcc_Rstudio.txt';
end
factors = [0.25:0.25:2.5];
data_write = [];
for ct_factor=1:length(factors)
    factor = factors(ct_factor);
    
    data_fac = data(abs([data.factor]-factor)<1e-6)';
    hasTakeoff_fac = hasTakeoff(abs([data.factor]-factor)<1e-6)';
    N = length(data_fac);
    
    % Create nominal and ordinal variables
    approach = arrayfun(@(i) i*ones(size(data_fac(i).intervals,1),1),1:N,'UniformOutput',false);
    side = arrayfun(@(i) data_fac(i).side*ones(size(data_fac(i).intervals,1),1),1:N,'UniformOutput',false);
    wind = arrayfun(@(i) data_fac(i).wind*ones(size(data_fac(i).intervals,1),1),1:N,'UniformOutput',false);
    time = arrayfun(@(i) data_fac(i).time*ones(size(data_fac(i).intervals,1),1),1:N,'UniformOutput',false);
    day = arrayfun(@(i) data_fac(i).day*ones(size(data_fac(i).intervals,1),1),1:N,'UniformOutput',false);
    
    hasTakeoff_fac = arrayfun(@(i) hasTakeoff_fac(i)*ones(size(data_fac(i).intervals,1),1),1:N,'UniformOutput',false);
    
    % create other variables
    y = -vertcat(data_fac.ymean_for_rdot);
    r = -vertcat(data_fac.rref);
    rdot = -vertcat(data_fac.slope_rvst);
    isRise = vertcat(data_fac.isRise);
    delta_r = vertcat(data_fac.delta_r);
    
    yStart = -vertcat(data_fac.yEntryStart);
    vStart = -vertcat(data_fac.vEntryStart);
    rStart = -vertcat(data_fac.rEntryStart);
    delta_U = vertcat(data_fac.delta_Uentry);
    delta_V = vertcat(data_fac.delta_Ventry);
    delta_W = vertcat(data_fac.delta_Wentry);
    delta_t = vertcat(data_fac.delta_tentry);
    amean = vertcat(data_fac.amean_entry);
    
    mean_Ua = vertcat(data_fac.mean_Ua);

    entryStartState = vertcat(data_fac.entryStartState);
    entryMeanState = vertcat(data_fac.entryMeanState);
    
    data_write = [data_write; ...
        vertcat(approach{:}) vertcat(side{:}) vertcat(wind{:}) ...
        vertcat(time{:}) vertcat(day{:}) y r rdot factor*ones(size(r,1),1) ...
        isRise delta_r vertcat(hasTakeoff_fac{:}) yStart vStart rStart delta_V delta_t amean mean_Ua...
        entryStartState entryMeanState delta_U delta_W];
end

if writeFile
    T = array2table(data_write, ...
        'VariableNames',{'approach','landingSide','wind','time','day','y',...
        'rref','rdot','threshold', 'isRise', 'delta_r', 'hasTakeoff', ...
        'ystart', 'vstart', 'rstart', 'deltaV', 'deltaT', 'Amean','mean_Ua',...
        'x0','y0','z0','u0','v0','w0','ax0','ay0','az0',...
        'xMean','yMean','zMean','uMean','vMean','wMean','axMean','ayMean','azMean','deltaU','deltaW'});
    writetable(T,r_file);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Code from hereon is not used for analysis
% %% Computing side and vertical movement data during transient and steady-state segments

% Inputs
close all; clc;
% clear;

if isunix
    inputFile = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/BlindLandingtracks_A3_LDF_rref_rrefEntry_with_rdot.mat';
    outputFile = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/BlindLandingtracks_A3_LDF_rref_rrefEntry_with_rdot_sideAccData.mat';
    DirPlots = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/plots/BlindLandingTracks/rref_entry_sideAcc';
elseif ispc
    inputFile = 'G:/steady_wind_experiments/postprocessing/BlindLandingtracks_A3_LDF_rref_rrefEntry_with_rdot.mat';
    outputFile = 'G:/steady_wind_experiments/postprocessing/BlindLandingtracks_A3_LDF_rref_rrefEntry_with_rdot_sideAccData.mat';
    DirPlots = 'G:/steady_wind_experiments/postprocessing/plots/BlindLandingTracks/rref_entry_sideAcc';
end

load(inputFile);

delPreviousPlots = false; % BE CAREFUL - when set to true, all previously saved plots in DirPlots are deleted
savePlots = false;
savePDFs = false;




