%% To find values of different parameters using optimization

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

%% Find r* intervals
% Also find videos asociated with each track

% Inputs
close all; clc;
clear;

inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A1.mat';
outputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A1_rref.mat';
outputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A1_rref_3dspeed.mat';
load(inputFile);

DirPlots = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/plots/BlindLandingTracks/rref_estimate';
% DirPlots = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/plots/BlindLandingTracks/rref_estimate_3dspeed';
delPreviousPlots = false; % BE CAREFUL - when set to true, all previously saved plots are deleted
savePlots = false;
savePDFs = false;

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high', 'lower'};
behaviour = {'rising','constant','sleeping'};

rmse_func = @(idx, indices, data) sqrt(sum((data(indices(idx,1):indices(idx,2)) - mean(data(indices(idx,1):indices(idx,2)))).^2)/(indices(idx,2)-indices(idx,1)+1));

min_gap = 14; % the first straight line is between current point and (current point + min_gap)th point
max_gap = 49;
params = [0.53 4.22]; % [sigma_{rmean-c_{r vs y, 0-1}}, sigma_{m_{r vs y, 0-1}}]
factors = [0.25:0.25:2.5];
time_window = [min_gap max_gap];
% factors = [1];
ct_data = 0;
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
                              
%                               excerpt.compute_rref(params, factors, time_window);
                              excerpt.compute_rref_with3dspeed(params, factors, time_window);
                              
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
end
keyboard
save(outputFile, 'treatments');
keyboard;

%% Extract rref data for statistical analysis
clc; close all;
% clear;
if isunix
    inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A1_rref.mat';
% inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A1_rref_3dspeed.mat';
elseif ispc
    inputFile = 'D:/light_intensity_experiments/postprocessing/BlindLandingtracks_A1_rref.mat';
end
% load(inputFile);
treatments = treatments(1:14*8); % Taking experiments for first 14 days
% combination

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high', 'lower'};
behaviour = {'rising','constant','sleeping'};

data = data4rrefEstimate.empty;
hasTakeoff = []; % whether or not the landing track has takeoff dynamics in it as defined in filteredState_BlindLandingtrack.hasTakeoff(..) method
for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
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
            
            
            for ct_treatment=1:length(relevantTreatments) % for each relevant treatment
                treatment = relevantTreatments(ct_treatment);
                
                arrayfun(@(x) x.setLandingSide(),[treatment.landingTracks.state_LDF]'); % to store landing side in the rrefSegments
                arrayfun(@(x) x.compute_params_basedon_3dspeed(), [treatment.landingTracks.state_LDF]');
                
                data1 = arrayfun(@(x) x.rrefSegments,[treatment.landingTracks.state_LDF]','UniformOutput',false);
                data1 = horzcat(data1{:});
                
                hastakeoff_pertreatment = arrayfun(@(x) x.hasTakeoff(treatment.landingDiscs)*ones(1,length(x.rrefSegments)),[treatment.landingTracks.state_LDF]','UniformOutput',false);
                hastakeoff_pertreatment = horzcat(hastakeoff_pertreatment{:});
                % Discard empty intervals
                indices = arrayfun(@(x) ~isempty(x.intervals_ti), data1);
                
                % Save additional data
                data1 = data1(indices);
                [data1.pattern] = deal(ct_pattern);
                [data1.light] = deal(ct_light);
                [data1.day] = deal(treatment.datenum);
                [data1.time] = deal(treatment.startTime);
    
                data = [data data1];
                hasTakeoff = [hasTakeoff hastakeoff_pertreatment(indices)];
%                 size(data)
%                 keyboard;
            end
        end
    end
end

%% Write file for statistical analysis in R
% This file contains data for all the factors
clc;
writeFile = true;
if isunix
    r_file = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/data_all_rref_Rstudio.txt';
% r_file = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/data_all_rref_Rstudio_3dspeed.txt';
elseif ispc
    r_file = 'D:/light_intensity_experiments/postprocessing/data_all_rref_Rstudio.txt';
end
factors = [0.25:0.25:2.5];
data_write = [];
for ct_factor=1:length(factors)
    factor = factors(ct_factor);
    
    data_fac = data(abs([data.factor]-factor)<1e-6)';
    hasTakeoff_fac = hasTakeoff(abs([data.factor]-factor)<1e-6)';
    N = length(data_fac);
    
    % Create nominal and ordinal variables
    approach = arrayfun(@(i) i*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    side = arrayfun(@(i) data_fac(i).side*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    pattern = arrayfun(@(i) data_fac(i).pattern*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    light = arrayfun(@(i) data_fac(i).light*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    time = arrayfun(@(i) data_fac(i).time*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    day = arrayfun(@(i) data_fac(i).day*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    
    hasTakeoff_fac = arrayfun(@(i) hasTakeoff_fac(i)*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);

    % create other variables
    y = -vertcat(data_fac.ymean_ti);
    r = -vertcat(data_fac.rref_ti);
    v = vertcat(data_fac.vmean_ti);
    
    speed3d = vertcat(data_fac.speed3d_mean_ti);
%     rSpeed3d = -vertcat(data_fac.rmean_speed3d_ti);
    
%     logy = log(y);
%     logr = log(V);
        
    data_write = [data_write; ...
        vertcat(approach{:}) vertcat(side{:}) vertcat(pattern{:}) ...
        vertcat(light{:}) vertcat(time{:}) vertcat(day{:}) y r v factor*ones(size(r,1),1) speed3d vertcat(hasTakeoff_fac{:})];
    
    N1 = sum(arrayfun(@(x) size(x.intervals_ti,1)>1, data_fac));
    disp(['# tracks: ' num2str(N) ', # data points: ' num2str(size(r,1)), ...
          ', # tracks with >1 r*: ' num2str(N1)]);
end

if writeFile
    T = array2table(data_write, ...
        'VariableNames',{'approach','landingSide','pattern','light','time','day','y','r','v','threshold', 'speed3d', 'hasTakeoff'});
    writetable(T,r_file);
end


%% Write file for statistical analysis in R
% This file contains data for all the factors, but only those tracks that
% contain >1 r* segments
clc;
writeFile = true;
r_file = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/data_multiple_rrefs_Rstudio.txt';
r_file = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/data_multiple_rrefs_Rstudio_3dspeed.txt';
data_write = [];
factors = [0.25:0.25:2.5];
for ct_factor=1:length(factors)
    factor = factors(ct_factor);
    
    data_fac = data(abs([data.factor]-factor)<1e-6)';
    
    % Choose segments that contain >1 r* segments
    has_multiple_rrefs = arrayfun(@(x) size(x.intervals_ti,1) > 1,data_fac);
    data_fac = data_fac(has_multiple_rrefs);
    
    N = length(data_fac);
    
    % Create nominal and ordinal variables
    approach = arrayfun(@(i) i*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    side = arrayfun(@(i) data_fac(i).side*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    pattern = arrayfun(@(i) data_fac(i).pattern*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    light = arrayfun(@(i) data_fac(i).light*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    time = arrayfun(@(i) data_fac(i).time*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    day = arrayfun(@(i) data_fac(i).day*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    
    % create other variables
    y = -vertcat(data_fac.ymean_ti);
    r = -vertcat(data_fac.rref_ti);
    v = vertcat(data_fac.vmean_ti);
    
%     logy = log(y);
%     logr = log(V);
    
    data_write = [data_write; ...
        vertcat(approach{:}) vertcat(side{:}) vertcat(pattern{:}) ...
        vertcat(light{:}) vertcat(time{:}) vertcat(day{:}) y r v factor*ones(size(r,1),1)];
    
    N1 = sum(arrayfun(@(x) size(x.intervals_ti,1)>1, data_fac));
    disp(['# tracks: ' num2str(N) ', # data points: ' num2str(size(r,1)), ...
          ', # tracks with >1 r*: ' num2str(N1)]);
end

if writeFile
    T = array2table(data_write, ...
        'VariableNames',{'approach','landingSide','pattern','light','time','day','y','r','v','threshold'});
    writetable(T,r_file);
end

%% Write file for statistical analysis in R
% This file contains data for all the factors, but only those tracks that
% contain exactly 1 r* segment
clc;
writeFile = true;
r_file = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/data_one_rrefs_Rstudio.txt';
data_write = [];
factors = [0.25:0.25:2.5];
for ct_factor=1:length(factors)
    factor = factors(ct_factor);
    
    data_fac = data(abs([data.factor]-factor)<1e-6)';
    
    % Choose segments that contain >1 r* segments
    has_one_rref = arrayfun(@(x) size(x.intervals_ti,1) == 1, data_fac);
    data_fac = data_fac(has_one_rref);
    
    N = length(data_fac);
    
    % Create nominal and ordinal variables
    approach = arrayfun(@(i) i*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    side = arrayfun(@(i) data_fac(i).side*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    pattern = arrayfun(@(i) data_fac(i).pattern*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    light = arrayfun(@(i) data_fac(i).light*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    time = arrayfun(@(i) data_fac(i).time*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    day = arrayfun(@(i) data_fac(i).day*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    
    % create other variables
    y = -vertcat(data_fac.ymean_ti);
    r = -vertcat(data_fac.rref_ti);
    v = vertcat(data_fac.vmean_ti);
    
%     logy = log(y);
%     logr = log(V);
    
    data_write = [data_write; ...
        vertcat(approach{:}) vertcat(side{:}) vertcat(pattern{:}) ...
        vertcat(light{:}) vertcat(time{:}) vertcat(day{:}) y r v factor*ones(size(r,1),1)];
    
    N1 = sum(arrayfun(@(x) size(x.intervals_ti,1)>1, data_fac));
    disp(['# tracks: ' num2str(N) ', # data points: ' num2str(size(r,1)), ...
          ', # tracks with >1 r*: ' num2str(N1)]);
end

if writeFile
    T = array2table(data_write, ...
        'VariableNames',{'approach','landingSide','pattern','light','time','day','y','r','v','threshold'});
    writetable(T,r_file);
end


%% Statistical analysis
writeFile = true;
r_file = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/data_rref_Rstudio.txt';
for ct_factor=4%1:length(factors)
    factor = factors(ct_factor);
    
    data_fac = data(abs([data.factor]-factor)<1e-6)';
    N = length(data_fac);
    
    % Create nominal and ordinal variables
    approach = arrayfun(@(i) i*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    approachNom = nominal(vertcat(approach{:}));
    
    side = arrayfun(@(i) data_fac(i).side*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    sideNom = nominal(vertcat(side{:})); % side
    
    pattern = arrayfun(@(i) data_fac(i).pattern*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    patternNom = nominal(vertcat(pattern{:}));
    
    light = arrayfun(@(i) data_fac(i).light*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    lightOrd = ordinal(vertcat(light{:}));
    
    time = arrayfun(@(i) data_fac(i).time*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    timeOrd = ordinal(vertcat(time{:}));
    
    day = arrayfun(@(i) data_fac(i).day*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    dayOrd = ordinal(vertcat(day{:}));
    
    % create other variables
    logy = log(-vertcat(data_fac.ymean_ti));
    logr = log(-vertcat(data_fac.rref_ti));
    
    if writeFile
        T = array2table([vertcat(approach{:}) vertcat(side{:}) vertcat(pattern{:}) ...
            vertcat(light{:}) vertcat(time{:}) vertcat(day{:}) logy logr], ...
            'VariableNames',{'approach','landingSide','pattern','light','time','day','logy','logr'});
        writetable(T,r_file);
    end
    
    A = dummyvar(patternNom); % N by 2
    B = dummyvar(lightOrd); % N by 3
    C = dummyvar(sideNom); % N by 2
    N = length(approach);
    % 
    
    % Create matrices
%     X = [ones(N, 1)    logy    s   B    C    A.*logy    B.*logy...
%         A(:,1).*logy.*B(:,1)    A(:,1).*logy.*B(:,2)    A(:,1).*logy.*B(:,3)...
%         A(:,2).*logy.*B(:,1)    A(:,2).*logy.*B(:,2)    A(:,2).*logy.*B(:,3)];
%     Z = {[ones(N, 1)], [ones(N,1)], [ones(N,1)]};
%     G = {[dayOrd], [timeOrd], [approach]};
%     lme = fitlmematrix(X,logr,Z,G)
    T = table(approachNom, sideNom, patternNom, ...
              lightOrd, timeOrd, dayOrd, logy, logr);
    
%     lme = fitlme(T, 'logr ~ logy + sideNom + lightOrd*logy + patternNom*logy + patternNom*lightOrd + sideNom*logy + (1|approachNom) + (1|timeOrd) + (1|dayOrd)');
    lme = fitlme(T, 'logr ~ logy + lightOrd*logy + patternNom*logy + patternNom*lightOrd + (1|approachNom) + (1|timeOrd) + (1|dayOrd) + (1|sideNom)');
    
end



%% Extract rref data for plotting and analysing 
clc; close all;
% clear;

% inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A1.mat';
% load(inputFile);

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
behaviour = {'rising','constant','sleeping'};

data = data4rrefEstimate.empty;
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
            
            
            for ct_treatment=1:length(relevantTreatments) % for each relevant treatment
                treatment = relevantTreatments(ct_treatment);
                
                data1 = arrayfun(@(x) x.rrefSegments,[treatment.landingTracks.state_LDF]','UniformOutput',false);
                data = [data horzcat(data1{:})];
%                 size(data)
%                 keyboard;
            end
        end
    end
end

%% Establish factor and time-window independence
clc; close all;
red_cmap = [252,187,161
252,146,114
251,106,74
239,59,44
203,24,29
165,15,21
103,0,13]./255;

data = data';
time_windows = [min_gap:max_gap]' + 1;

intercepts = nan*ones(length(factors), length(time_windows));
slopes = nan*ones(length(factors), length(time_windows));
coeffs_fac = nan*zeros(length(factors),2);
for ct_fac=1:length(factors)
    factor = factors(ct_fac);
    data_fac = data(abs([data.factor]-factor)<1e-6);
    indices = arrayfun(@(x) ~isempty(x.intervals_td), data_fac);
    data_fac = data_fac(indices);
%     size(data_fac)
    rref = arrayfun(@(x) x.rref_td,data_fac,'UniformOutput',false);
    rref = vertcat(rref{:});
    
    ymean = arrayfun(@(x) x.ymean_td,data_fac,'UniformOutput',false);
    ymean = vertcat(ymean{:});
    
    vmean = arrayfun(@(x) x.vmean_td,data_fac,'UniformOutput',false);
    vmean = vertcat(vmean{:});
    
    
    
    tws = arrayfun(@(x) diff(x.intervals_td,[],2) + 1,data_fac,'UniformOutput',false);
    tws = vertcat(tws{:});
    for ct_tw=1:length(time_windows)
        tw = time_windows(ct_tw);
        indices = tws==tw;
        coeffs = [ones(sum(indices),1) log(-ymean(indices))]\log(-rref(indices));
        intercepts(ct_fac,ct_tw) = coeffs(1);
        slopes(ct_fac,ct_tw) = coeffs(2);
        
%         figure;
%         scatter(-ymean(indices), -rref(indices),40,red_cmap(1,:),'filled','s','MarkerEdgeColor','k');
        
    end
end
[X,Y] = meshgrid(time_windows, factors);
figure;
surf(X,Y,intercepts,'FaceAlpha',0.5, 'EdgeColor', 'none')
view(90,0);
figure;
surf(X,Y,slopes,'FaceAlpha',0.5, 'EdgeColor', 'none')
view(90,0);
view(0,0);
figure;
plot(factors, intercepts,'o')
figure;
plot(factors, slopes,'o')
%% NOTES 
% 1) Can't use post-hoc test emmeans available at https://nl.mathworks.com/matlabcentral/fileexchange/71970-emmeans
%    This is because " For now, only output from fitglme can be used. 
%    Major limitation is that only interactions between categorical predictor
%    variables are accepted (not between continuous variables or categorical-continuous interactions)."
% 2) Plus, we need a function similar to emtrends in R to assess comparison
%    of difference in slopes of covariate logy in different conditions.
%    Such a function is not available in Matlab as of now.
% 3) This means R has to be used for the statistical analysis.

%%
% clc; close all;

red_cmap = [252,187,161
252,146,114
251,106,74
239,59,44
203,24,29
165,15,21
103,0,13]./255;

data = data';
coeffs = nan*zeros(length(factors),2);
coeffs1 = nan*zeros(length(factors),3);
% coeffs_robust = nan*zeros(length(factors),2);
for ct=1:length(factors)
    data_fac = data(abs([data.factor]-factors(ct))<1e-6);
    indices = arrayfun(@(x) ~isempty(x.intervals_ti), data_fac);
    data_fac = data_fac(indices);
    size(data_fac)
    
    rref = arrayfun(@(x) x.rref_ti,data_fac,'UniformOutput',false);
    rref = vertcat(rref{:});
    
    ymean = arrayfun(@(x) x.ymean_ti,data_fac,'UniformOutput',false);
    ymean = vertcat(ymean{:});
    
    vmean = arrayfun(@(x) x.vmean_ti,data_fac,'UniformOutput',false);
    vmean = vertcat(vmean{:});
    
    xmean = arrayfun(@(x) x.xmean_ti,data_fac,'UniformOutput',false);
    xmean = vertcat(xmean{:});
    
    zmean = arrayfun(@(x) x.zmean_ti,data_fac,'UniformOutput',false);
    zmean = vertcat(zmean{:});
    
    xTravelled = arrayfun(@(x) x.xTravelled_ti,data_fac,'UniformOutput',false);
    xTravelled = vertcat(xTravelled{:});
    
    yTravelled = arrayfun(@(x) x.yTravelled_ti,data_fac,'UniformOutput',false);
    yTravelled = vertcat(yTravelled{:});
    
    zTravelled = arrayfun(@(x) x.zTravelled_ti,data_fac,'UniformOutput',false);
    zTravelled = vertcat(zTravelled{:});
    
    fd_analytical = arrayfun(@(x) x.fd_analytical_ti,data_fac,'UniformOutput',false);
    fd_analytical = vertcat(fd_analytical{:});
    
    fd_actual = arrayfun(@(x) x.fd_actual_ti,data_fac,'UniformOutput',false);
    fd_actual = vertcat(fd_actual{:});
    
    nPoints = arrayfun(@(x) diff(x.intervals_ti,[],2)+1,data_fac,'UniformOutput',false);
    nPoints = vertcat(nPoints{:});
    
    if ~isempty(ymean) && ~isempty(rref)
        coeffs(ct,:) = [ones(size(ymean)) log(-ymean)]\log(-rref);
%         f = fitlm(log(-ymean), log(-rref));
%         coeffs(ct,:) = [f.Coefficients.Estimate]';
%         
%         f1 = fitlm([log(-ymean) (nPoints)], log(-rref));
%         coeffs1(ct,:) = [f1.Coefficients.Estimate]';
%         
%         f2 = fitlm([log(-ymean) (nPoints)], log(-rref), 'interactions')
%         f3=fitlm([xmean log(-ymean) zmean xTravelled yTravelled zTravelled nPoints],log(-rref))
%         f4=fitlm([log(-ymean) log(yTravelled) log(nPoints)],log(-rref))
%         f5=fitlm([log(-ymean) log(yTravelled) log(nPoints)],log(-rref), 'interactions')
%         fitlm(log(-ymean), log(vmean))
%         fitlm([log(-ymean) nPoints], log(vmean))
%         fitlm(yTravelled, -rref)
%         fitlm(log(vmean), log(-rref))
%         fitlm(yTravelled, nPoints)
%         f1 = fitlm([log(-ymean) (nPoints)], log(-rref))
%         f1 = fitlm(nPoints, -rref)
%         f1 = fitlmematrix(log(-ymean), log(-rref), nPoints, [])
%         f1 = fitlmematrix([ones(size(ymean)) log(-ymean)], log(-rref), ones(size(ymean)), nPoints)

%         f_robust = fitlm(log(-ymean), log(-rref), 'RobustOpts', 'on');
%         coeffs_robust(ct,:) = [f_robust.Coefficients.Estimate]';
    
%         figure;
%         scatter(nPoints,-rref,40,red_cmap(1,:),'filled','s','MarkerEdgeColor','k');
        
%         histogram(nPoints)

%         figure;
%         scatter(-ymean, -rref,40,red_cmap(1,:),'filled','s','MarkerEdgeColor','k');
%         figure;
%         scatter(-ymean, vmean,40,red_cmap(1,:),'filled','s','MarkerEdgeColor','k');
%         figure;
%         scatter(xmean, -rref,40,red_cmap(1,:),'filled','s','MarkerEdgeColor','k');
%         figure;
%         scatter(zmean, -rref,40,red_cmap(1,:),'filled','s','MarkerEdgeColor','k');
%         figure;
%         scatter(xTravelled, -rref,40,red_cmap(1,:),'filled','s','MarkerEdgeColor','k');
%         figure;
%         scatter(yTravelled, -rref,40,red_cmap(1,:),'filled','s','MarkerEdgeColor','k');
%         figure;
%         scatter(zTravelled, -rref,40,red_cmap(1,:),'filled','s','MarkerEdgeColor','k');
%         figure;
%         scatter(nPoints, -rref,40,red_cmap(1,:),'filled','s','MarkerEdgeColor','k');
%         figure;
%         scatter(-ymean, nPoints,40,red_cmap(1,:),'filled','s','MarkerEdgeColor','k');
%         figure;
%         scatter(yTravelled, nPoints,40,red_cmap(1,:),'filled','s','MarkerEdgeColor','k');
%         figure;
%         scatter(log(yTravelled), log(-rref),40,red_cmap(1,:),'filled','s','MarkerEdgeColor','k');
%         figure;
%         scatter(log(-ymean), log(-rref),40,red_cmap(1,:),'filled','s','MarkerEdgeColor','k');
%         figure;
%         scatter(log(nPoints), log(-rref),40,red_cmap(1,:),'filled','s','MarkerEdgeColor','k');
%         figure;
%         scatter(log(vmean), log(-rref),40,red_cmap(1,:),'filled','s','MarkerEdgeColor','k');
    
%     figure;
%     scatter(-ymean, -rref,40,red_cmap(1,:),'filled','s','MarkerEdgeColor','k');
    
%     figure;
%     scatter(log(-ymean), log(-rref),40,red_cmap(1,:),'filled','s','MarkerEdgeColor','k');
    end
end
figure;
plot(factors, coeffs)
ylim([-1 -0.5])
% 
% figure;
% plot(factors,coeffs_robust)

%% FOR DEBUGGING
a = treatments(78).landingTracks(80)
a.removeDuplicateTracks();

%%


% % Eliminate one data point that has abs(dof_analytical-dof_actual) > 1/fps;
% to_del = [];
% for ct=1:length(data)
%     for ct1=1:length(data(ct).dataPerTrackExcerpt)
%         if abs(data(ct).dataPerTrackExcerpt(ct1).dof_analytical-data(ct).dataPerTrackExcerpt(ct1).dof_actual)>1/175
%             to_del = [to_del; ct, ct1]
%         end
%     end
% end
% data(to_del(1,1)) = [];
keyboard;
% Information about the extracted data
has_lt1_rrefs = arrayfun(@(x) isempty(x.dataPerTrackExcerpt) ,data); % has_lessThanOne_rrefs
empty_data = data(has_lt1_rrefs); % should be empty

has_gt1_rrefs = arrayfun(@(x) length(x.dataPerTrackExcerpt)>1 ,data); % has_greaterThanOne_rrefs
data_multiple_rrefs = data(has_gt1_rrefs);
data_single_rref = data(~has_gt1_rrefs);

all_data4rrefEstimates = horzcat(data(:).dataPerTrackExcerpt);
rref_fitVvsy = [all_data4rrefEstimates.rref]';
rref_meanVbyy = [all_data4rrefEstimates.meanVbyy]';
Rsquared = [all_data4rrefEstimates.Rsquared]';
RMSE = arrayfun(@(x) x.model.RMSE, all_data4rrefEstimates); 
vmean = [all_data4rrefEstimates.vmean]';
ymean = [all_data4rrefEstimates.ymean]';
dof_analytical = [all_data4rrefEstimates.dof_analytical]';
dof_actual = [all_data4rrefEstimates.dof_actual]';

nPoints = arrayfun(@(x) size(x.state4rrefEstimate,1), all_data4rrefEstimates)';
xTravelled = arrayfun(@(x) abs(max(x.state4rrefEstimate(:,2)) - min(x.state4rrefEstimate(:,2))), all_data4rrefEstimates)';
yTravelled = arrayfun(@(x) abs(diff(x.state4rrefEstimate([1 end],3),1)), all_data4rrefEstimates)';
zTravelled = arrayfun(@(x) abs(max(x.state4rrefEstimate(:,4)) - min(x.state4rrefEstimate(:,4))), all_data4rrefEstimates)';
xMean = arrayfun(@(x) mean(x.state4rrefEstimate(:,2)), all_data4rrefEstimates)';
yMean = arrayfun(@(x) mean(x.state4rrefEstimate(:,3)), all_data4rrefEstimates)';
zMean = arrayfun(@(x) mean(x.state4rrefEstimate(:,4)), all_data4rrefEstimates)';

% Compare rref_meanVbyy and rref_fitVvsy
[rref_meanVbyy rref_fitVvsy]
max(abs(rref_meanVbyy-rref_fitVvsy))
mean(abs(rref_meanVbyy-rref_fitVvsy))
median(abs(rref_meanVbyy-rref_fitVvsy))
figure;
boxplot(abs(rref_meanVbyy-rref_fitVvsy), 'OutlierSize', 0.1, 'Labels', 'r*-rmean'); hold on;
f = scatter(ones(size(rref_fitVvsy)).*(1+(rand(size(rref_fitVvsy))-0.5)/10), abs(rref_meanVbyy-rref_fitVvsy), 'k', 'filled');
f.MarkerFaceAlpha = 0.5;
set(gca, 'FontSize', 16);

figure; 
max(Rsquared)
mean(Rsquared)
median(Rsquared)
min(Rsquared)
boxplot(Rsquared, 'OutlierSize', 0.1, 'Labels', 'R^2'); hold on;
f = scatter(ones(size(Rsquared)).*(1+(rand(size(Rsquared))-0.5)/10),Rsquared,'k','filled');
f.MarkerFaceAlpha = 0.5;
set(gca, 'FontSize', 16);

figure; 
max(RMSE)
mean(RMSE)
median(RMSE)
min(RMSE)
boxplot(RMSE, 'OutlierSize', 0.1, 'Labels', 'RMSE'); hold on;
f = scatter(ones(size(RMSE)).*(1+(rand(size(RMSE))-0.5)/10),RMSE,'k','filled');
f.MarkerFaceAlpha = 0.5;
set(gca, 'FontSize', 16);

figure;
histogram(Rsquared);
xlabel('R^2', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);

% Histogram of r*
figure;
histogram(rref_fitVvsy);
xlabel('Reference rate of expansion (1/s)', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);

% Insights into the selected data intervals
min(nPoints)
figure;
histogram(nPoints);
xlabel('# of data points for r* estimate', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);

figure;
histogram(yTravelled);
xlabel('y displacement for r* estimate (m)', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);

% Compare analytical and actual duration of flights for r* intervals
figure;
histogram(dof_analytical); hold on;
histogram(dof_actual);
xlabel('Duration of flight (s)', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);

[dof_analytical dof_actual]
max(abs(dof_analytical-dof_actual))
figure;
boxplot(abs(dof_analytical-dof_actual), 'OutlierSize', 0.1, 'Labels', 'dof_analytical-dof_actual'); hold on;
f = scatter(ones(size(dof_actual)).*(1+(rand(size(dof_actual))-0.5)/10), abs(dof_analytical-dof_actual), 'k', 'filled');
f.MarkerFaceAlpha = 0.5;
set(gca, 'FontSize', 16);

% % Compute extra parameters
% ystart = zeros(length(all_states),1);
% ideal_tof = zeros(length(all_states),1); % time of flight computed analytically for the interval in which r is held constant
% real_tof = zeros(length(all_states),1); % time of flight actually observed for the interval in which r is held constant
% 
% tofEstimate = false(length(all_states),1);
% for ct=1:length(all_states)
%     ystart(ct) = all_states(ct).state4rrefEstimate(1,3);
%     ideal_tof(ct) = log(abs(all_states(ct).state4rrefEstimate(end,3)/all_states(ct).state4rrefEstimate(1,3)))/rref_fitVvsy(ct);
%     real_tof(ct) = diff(all_states(ct).state4rrefEstimate([1, end],1));
%     if ~isempty(all_states(ct).state4tofEstimate)
%         state = all_states(ct).state4tofEstimate;
%         tofEstimate(ct) = true;
%         all_states(ct).distance = sum(sum(diff(state(:,2:4)).^2,2).^0.5);
%         all_states(ct).ydistance = sum(abs(diff(state(:,3))));
%         all_states(ct).displacement = sum(diff(state([1, end],2:4)).^2,2).^0.5;
%         all_states(ct).ydisplacement = abs(diff(state([1, end],3)));
%         all_states(ct).ideal_tof_ylimit = log(abs(state(end,3)/state(1,3)))/rref_fitVvsy(ct);
%         all_states(ct).real_tof_ylimit = state(end,1)-state(1,1);
%     else
%         all_states(ct).distance = [];
%         all_states(ct).ydistance = [];
%         all_states(ct).displacement = [];
%         all_states(ct).ydisplacement = [];
%         all_states(ct).ideal_tof_ylimit = [];
%         all_states(ct).real_tof_ylimit = [];
%     end
% end
% turtuosity = [all_states(:).distance]'./[all_states(:).displacement]';
% yturtuosity = [all_states(:).displacement]'./[all_states(:).ydisplacement]';
% [turtuosity  yturtuosity]
% figure;
% boxplot(turtuosity); hold on;
% f = scatter(ones(size(turtuosity)).*(1+(rand(size(turtuosity))-0.5)/10),turtuosity,'k','filled');
% 
% 
% [ideal_tof real_tof]
% delta_t = abs([ideal_tof-real_tof]);
% max(delta_t) % It implies manual picking of intervals is super good!!
% figure; 
% boxplot(delta_t, 'OutlierSize', 0.00000001, 'Labels', 'fd_{ideal}-fd_{actual} (s)'); hold on;
% f = scatter(ones(size(delta_t)).*(1+(rand(size(delta_t))-0.5)/10),delta_t,'k','filled');
% f.MarkerFaceAlpha = 0.5;
% set(gca, 'FontSize', 16);


% % Compute tof flight (ideal vs real) for ylimit distance travel
% % ideal - when bbee would have continued to fly at the computed r
% ideal_tof_ylimit = [all_states(:).ideal_tof_ylimit]';
% real_tof_ylimit = [all_states(:).real_tof_ylimit]';
% 
% [ideal_tof_ylimit./real_tof_ylimit turtuosity]
% figure; hold on;
% plot(turtuosity, real_tof_ylimit./ideal_tof_ylimit,'o')
% yline(1,'-')
% xlabel('Turtuosity', 'FontSize', 16);
% ylabel('fd_{real}/fd_{ideal}', 'FontSize', 16);
% set(gca, 'FontSize', 16);

% Create plot of ALL r vs y points
states = vertcat(all_data4rrefEstimates(:).state4rrefEstimate);
b1 = robustfit(states(:,3),states(:,6)./states(:,3));
c1 = fit(states(:,3),states(:,6)./states(:,3),'poly1');
figure;
hold on;
plot(states(:,3),states(:,6)./states(:,3),'s','MarkerSize',6,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6]);
plot(states(:,3),c1(states(:,3)),'b','LineWidth',2);
plot(states(:,3),b1(1)+b1(2)*states(:,3),'k','LineWidth',2);
ylabel('Reference rate of expansion (1/s)', 'FontSize', 16);
xlabel('y from the platform (m)', 'FontSize', 16);
set(gca, 'FontSize', 18); grid on;

% Plot mean rref vs ymean
b = robustfit(ymean,rref_fitVvsy);
% mdl = fitlm(ymean,rref,'RobustOpts','on');
[c, gof] = fit(ymean,rref_fitVvsy,'poly1');
figure;
hold on;
% plot(ymean,rref_fitVvsy,'s','MarkerSize',6,...
%     'MarkerEdgeColor','red',...
%     'MarkerFaceColor',[1 .6 .6]);
% % f = scatter(ymean,rref_fitVvsy,'s',...
% %     'MarkerEdgeColor',[215 48 39]./255,...
% %     'MarkerFaceColor',[215 48 39]./255);
% % f.MarkerFaceAlpha = 0.5;
p1 = plot(ymean,rref_fitVvsy,'s','MarkerSize',6,...
    'MarkerEdgeColor',[215 48 39]./255,...
    'MarkerFaceColor',[252,146,114]./255);
% plot(ymean,b(1)+b(2)*ymean,'k','LineWidth',2);
plot(ymean,c(ymean),'Color', [69 117 180]./255,'LineWidth',3);
% plot(ymean,mdl(ymean),'b^','LineWidth',2);
% ylabel('Mean reference rate of expansion (1/s)', 'FontSize', 16);
ylabel('Reference rate of expansion, r* (1/s)', 'FontSize', 16);
xlabel('Mean y from the platform (m)', 'FontSize', 16);
set(gca, 'FontSize', 18); grid on;

% Create plot of ALL Vgy vs y points
states = vertcat(all_data4rrefEstimates(:).state4rrefEstimate);
b1 = robustfit(states(:,3),states(:,6));
c1 = fit(states(:,3),states(:,6),'poly1');
figure;
hold on;
plot(states(:,3),states(:,6),'s','MarkerSize',6,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6]);
plot(states(:,3),c1(states(:,3)),'b','LineWidth',2);
plot(states(:,3),b1(1)+b1(2)*states(:,3),'k','LineWidth',2);
ylabel('V_{gy} (m/s)', 'FontSize', 16);
xlabel('y from the platform (m)', 'FontSize', 16);
set(gca, 'FontSize', 18); grid on;

% Plot vmean vs ymean
% b = robustfit(ymean,vmean);
c = fit(ymean,vmean,'poly1');
figure;
hold on;
plot(ymean,vmean,'s','MarkerSize',6,...
    'MarkerEdgeColor',[215 48 39]./255,...
    'MarkerFaceColor',[252,146,114]./255);
% plot(ymean,b(1)+b(2)*ymean,'r','LineWidth',2);
plot(ymean,c(ymean),'Color', [69 117 180]./255,'LineWidth',3);
ylabel('Mean V_{gy} (m/s)', 'FontSize', 16);
xlabel('Mean y from the platform (m)', 'FontSize', 16);
set(gca, 'FontSize', 18); grid on;

%% Create figures for the Article - Bbees regulate r* during landing
% Histogram of r*
figure;
histogram(-rref_fitVvsy/2, [0:0.5:5]);
xlabel('Estimated set-points of optical expansion rate, r* (rad/s)', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);

% Plot mean rref vs ymean with yTravelled colormap
y1 = linspace(min(ymean),max(ymean),100);
red_cmap = [252,187,161
252,146,114
251,106,74
239,59,44
203,24,29
165,15,21
103,0,13]./255;
[data_cmap, edges] = find_cmap(yTravelled, red_cmap);
figure;
colormap(red_cmap);
hold on;
% plot(ymean,rref_fitVvsy,'s','MarkerSize',6,...
%     'MarkerEdgeColor','red',...
%     'MarkerFaceColor',[1 .6 .6]);
f = scatter(-ymean,-rref_fitVvsy,40,red_cmap(data_cmap',:),'filled','s');
% f = scatter(-ymean,-rref_fitVvsy,'s',...
%     'MarkerEdgeColor',[215 48 39]./255,...
%     'MarkerFaceColor',[252,146,114]./255);
% f.MarkerFaceAlpha = 0.9;
cmap_bar = colorbar('eastoutside');
caxis(edges([1 end]));
cmap_bar.Ticks = edges;
cmap_bar.Label.String = 'y travelled';
cmap_bar.Label.FontSize = 16;
% p1 = plot(-ymean,-rref_fitVvsy,'s','MarkerSize',6,...
%     'MarkerEdgeColor',[215 48 39]./255,...
%     'MarkerFaceColor',[252,146,114]./255);
% plot(ymean,b(1)+b(2)*ymean,'k','LineWidth',2);

% plot(-y1,c1(-1./y1),'Color', [69 117 180]./255,'LineWidth',3);

set(gca, 'FontSize', 18);
% plot(ymean,mdl(ymean),'b^','LineWidth',2);
% ylabel('Mean reference rate of expansion (1/s)', 'FontSize', 16);
ylabel('Estimated set-point of optical expansion rate, r* (rad/s)', 'FontSize', 16);
xlabel('Mean y (m)', 'FontSize', 16);

% % Plot connecting lines
ct2plot = [24 26 33];
for ct=1:length(ct2plot)
    ct_now = ct2plot(ct);
    indx = find(arrayfun(@(x) x.dataPerTrackExcerpt(1), data) == data_multiple_rrefs(ct_now).dataPerTrackExcerpt(1));
    yTravelled_indx = arrayfun(@(x) abs(diff(x.state4rrefEstimate([1 end],3),1)), data(indx).dataPerTrackExcerpt)';
    
    % Plot straight line
    plot(-[data(indx).dataPerTrackExcerpt.ymean], -[data(indx).dataPerTrackExcerpt.rref] ...
    ,'Color', [225 225 225]./255/2,'LineWidth',1);

    % Highlight data points
    scatter(-[data(indx).dataPerTrackExcerpt.ymean], -[data(indx).dataPerTrackExcerpt.rref],40,...
    red_cmap(discretize(yTravelled_indx, edges),:),'filled','s','MarkerEdgeColor','k');

end
plot(-ymean,c(-ymean),'Color', [69 117 180]./255,'LineWidth',3);
plot(-ymean,c(-ymean),'Color', [111 172 239]./255,'LineWidth',3);

% % Plot only multiple rref points
figure;
hold on;
ct2plot = 1:length(data_multiple_rrefs); % [24 26 33]
for ct=1:length(ct2plot)
    ct_now = ct2plot(ct);
    indx = find(arrayfun(@(x) x.dataPerTrackExcerpt(1), data) == data_multiple_rrefs(ct_now).dataPerTrackExcerpt(1));
    yTravelled_indx = arrayfun(@(x) abs(diff(x.state4rrefEstimate([1 end],3),1)), data(indx).dataPerTrackExcerpt)';
    
    % Plot straight line
    plot(-[data(indx).dataPerTrackExcerpt.ymean], -[data(indx).dataPerTrackExcerpt.rref] ...
    ,'Color', [225 225 225]./255/2,'LineWidth',1);

    % Highlight data points
    scatter(-[data(indx).dataPerTrackExcerpt.ymean], -[data(indx).dataPerTrackExcerpt.rref],40,...
    red_cmap(1,:),'filled','s','MarkerEdgeColor','k');
%     scatter(-[data(indx).dataPerTrackExcerpt.ymean], -[data(indx).dataPerTrackExcerpt.rref],40,...
%     red_cmap(discretize(yTravelled_indx, edges),:),'filled','s','MarkerEdgeColor','k');
end
set(gca, 'FontSize', 18);
ylabel('Estimated set-point of optical expansion rate, r* (rad/s)', 'FontSize', 16);
xlabel('Mean y (m)', 'FontSize', 16);

% % % % % % % % % % % % % % % % indx_24 = find([data.track] == data_multiple_rrefs(24).track);
% % % % % % % % % % % % % % % % indx_26 = find([data.track] == data_multiple_rrefs(26).track);
% % % % % % % % % % % % % % % % indx_33 = find([data.track] == data_multiple_rrefs(33).track);
% % % % % % % % % % % % % % % % yTravelled_24 = arrayfun(@(x) abs(diff(x.state4rrefEstimate([1 end],3),1)), data(indx_24).dataPerTrackExcerpt)';
% % % % % % % % % % % % % % % % yTravelled_26 = arrayfun(@(x) abs(diff(x.state4rrefEstimate([1 end],3),1)), data(indx_26).dataPerTrackExcerpt)';
% % % % % % % % % % % % % % % % yTravelled_33 = arrayfun(@(x) abs(diff(x.state4rrefEstimate([1 end],3),1)), data(indx_33).dataPerTrackExcerpt)';
% % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % plot(-[data(indx_33).dataPerTrackExcerpt.ymean], -[data(indx_33).dataPerTrackExcerpt.rref] ...
% % % % % % % % % % % % % % % %     ,'Color', [150 150 150]./255/2,'LineWidth',2);
% % % % % % % % % % % % % % % % plot(-[data(indx_24).dataPerTrackExcerpt.ymean], -[data(indx_24).dataPerTrackExcerpt.rref] ...
% % % % % % % % % % % % % % % %     ,'Color', [150 150 150]./255/2,'LineWidth',2);
% % % % % % % % % % % % % % % % plot(-[data(indx_26).dataPerTrackExcerpt.ymean], -[data(indx_26).dataPerTrackExcerpt.rref] ...
% % % % % % % % % % % % % % % %     ,'Color', [150 150 150]./255/2,'LineWidth',2);
% % % % % % % % % % % % % % % % scatter(-[data(indx_33).dataPerTrackExcerpt.ymean], -[data(indx_33).dataPerTrackExcerpt.rref],40,...    
% % % % % % % % % % % % % % % %     red_cmap(discretize(yTravelled_33, edges),:),'filled','s','MarkerEdgeColor','k');
% % % % % % % % % % % % % % % % scatter(-[data(indx_24).dataPerTrackExcerpt.ymean], -[data(indx_24).dataPerTrackExcerpt.rref],40,...
% % % % % % % % % % % % % % % %     red_cmap(discretize(yTravelled_24, edges),:),'filled','s','MarkerEdgeColor','k');
% % % % % % % % % % % % % % % % scatter(-[data(indx_26).dataPerTrackExcerpt.ymean], -[data(indx_26).dataPerTrackExcerpt.rref],40,...
% % % % % % % % % % % % % % % %     red_cmap(discretize(yTravelled_26, edges),:),'filled','s','MarkerEdgeColor','k');

% Plot yTravelled vs rref
figure;
d1=fitlm(yTravelled,-rref_fitVvsy)
hold on;
plot(yTravelled,-rref_fitVvsy,'s','MarkerSize',6,...
    'MarkerEdgeColor',[215 48 39]./255,...
    'MarkerFaceColor',[252,146,114]./255);
% plot(ymean,b(1)+b(2)*ymean,'r','LineWidth',2);
ylabel('Reference rate of expansion, r* (1/s)', 'FontSize', 16);
xlabel('Distance travelled (m)', 'FontSize', 16);
set(gca, 'FontSize', 18);

% Plot mean rref vs ymean 
red_cmap = [252,187,161
252,146,114
251,106,74
239,59,44
203,24,29
165,15,21
103,0,13]./255;
d1=fitlm(log(-ymean),log(-rref_fitVvsy))
d1.Rsquared.Ordinary
modelfun = @(b,x)(exp(b(1))*x.^(b(2)));
figure;
hold on;
% plot(ymean,rref_fitVvsy,'s','MarkerSize',6,...
%     'MarkerEdgeColor','red',...
%     'MarkerFaceColor',[1 .6 .6]);
f = scatter(-ymean,-rref_fitVvsy,40,red_cmap(1,:),'filled','s','MarkerEdgeColor','k');
% f = scatter(-ymean,-rref_fitVvsy,'s',...
%     'MarkerEdgeColor',[215 48 39]./255,...
%     'MarkerFaceColor',[252,146,114]./255);
% f.MarkerFaceAlpha = 0.9;
set(gca, 'FontSize', 18);
ylabel('Estimated set-point of optical expansion rate, r* (rad/s)', 'FontSize', 16);
xlabel('Mean y (m)', 'FontSize', 16);
plot(min(-ymean):0.001:max(-ymean),modelfun(d1.Coefficients.Estimate,min(-ymean):0.001:max(-ymean)),'Color', [69 117 180]./255,'LineWidth',3);

% % Plot log(rref) vs log(ymean)
modelfun = @(b,x)(b(1)+b(2).*x);
x_vec = min(log(-ymean)):0.001:max(log(-ymean));
figure; hold on;
f = scatter(log(-ymean),log(-rref_fitVvsy),40,red_cmap(1,:),'filled','s','MarkerEdgeColor','k');
plot(x_vec,modelfun(d1.Coefficients.Estimate,x_vec),'Color', [69 117 180]./255,'LineWidth',3);
set(gca, 'FontSize', 18);
ylabel('log(r*) (rad/s)', 'FontSize', 16);
xlabel('log(y_{mean}) (m)', 'FontSize', 16);

% Plot vmean vs ymean 
d2=fitlm(log(-ymean),log(vmean))
d2.Rsquared.Ordinary
modelfun = @(b,x)(exp(b(1))*x.^(b(2)));
x_vec = min(-ymean):0.001:max(-ymean);
figure;
hold on;
f = scatter(-ymean,vmean,40,red_cmap(1,:),'filled','s','MarkerEdgeColor','k');
% plot(-ymean,vmean,'s','MarkerSize',6,...
%     'MarkerEdgeColor',[215 48 39]./255,...
%     'MarkerFaceColor',[252,146,114]./255);
plot(x_vec,modelfun(d2.Coefficients.Estimate,x_vec),'Color', [69 117 180]./255,'LineWidth',3);
ylabel('Mean V (m/s)', 'FontSize', 16);
xlabel('Mean y (m)', 'FontSize', 16);
set(gca, 'FontSize', 18);

% % Plot log(-ymean) vs log(vmean) to check if it is a straight line or
% not
modelfun = @(b,x)(b(1)+b(2).*x);
x_vec = min(log(-ymean)):0.001:max(log(-ymean));
figure; hold on;
f = scatter(log(-ymean),log(vmean),40,red_cmap(1,:),'filled','s','MarkerEdgeColor','k');
plot(x_vec,modelfun(d2.Coefficients.Estimate,x_vec),'Color', [69 117 180]./255,'LineWidth',3);
set(gca, 'FontSize', 18);
ylabel('log(V_{mean}) (m/s)', 'FontSize', 16);
xlabel('log(y_{mean}) (m)', 'FontSize', 16);

figure; hold on;
f = scatter(-ymean,vmean,20,red_cmap(1,:),'filled','o');
x_vec = min(-ymean):0.001:max(-ymean);
plot(x_vec,modelfun(d2.Coefficients.Estimate,x_vec),'Color', [69 117 180]./255,'LineWidth',3);
set(gca, 'FontSize', 18);
ylabel('V (m/s)', 'FontSize', 16);
xlabel('y (m)', 'FontSize', 16);
set(gca, 'FontSize', 18);
yticks([0:0.2:1])
ylim([0 1])
plot([0 0.3], [0 0.3],'--k','LineWidth',2);
plot([0 0.3], [0 2*0.3],'--k','LineWidth',2);
plot([0 0.3], [0 3*0.3],'--k','LineWidth',2);
plot([0 0.3], [0 4*0.3],'--k','LineWidth',2);
plot([0 0.3], [0 5*0.3],'--k','LineWidth',2);
plot([0 0.3], [0 6*0.3],'--k','LineWidth',2);
plot([0 0.3], [0 7*0.3],'--k','LineWidth',2);
% steady-state performance of the controller
figure;
d1=fitlm(-rref_fitVvsy, RMSE)
d1.Rsquared.Ordinary
x_vec = min(-rref_fitVvsy):0.01:max(-rref_fitVvsy);
modelfun = @(b,x)(b(1)+b(2).*x);
hold on;
f = scatter(-rref_fitVvsy, RMSE,40,red_cmap(1,:),'filled','s','MarkerEdgeColor','k');
plot(x_vec,modelfun(d1.Coefficients.Estimate,x_vec),'Color', [69 117 180]./255,'LineWidth',3);
set(gca, 'FontSize', 18);
ylabel('RMSE in r* (rad/s)', 'FontSize', 16);
xlabel('Estimated set-point of optical expansion rate, r* (rad/s)', 'FontSize', 16);
xticks([0:2:10])

% Checking dependency of yTravelled, xTravelled etc. on rref and vmean
% variation (???????)
d1=fitlm([xMean -yMean zMean xTravelled yTravelled zTravelled],-rref_fitVvsy)
d1=fitlm([xMean -yMean zMean],-rref_fitVvsy)
d1=fitlm([log(abs(xMean)) log(-yMean) log(abs(zMean))], log(-rref_fitVvsy))
d1=fitlm([xMean ],-rref_fitVvsy)
d1=fitlm([zMean],-rref_fitVvsy)
d1=fitlm([xTravelled],-rref_fitVvsy)
d1=fitlm([yTravelled],-rref_fitVvsy)
d1=fitlm([zTravelled],-rref_fitVvsy)
d2=fitlm(log(-ymean),log(-rref_fitVvsy))
figure;
plot(xTravelled, -rref_fitVvsy, 's')
plot(yTravelled, -rref_fitVvsy, 's')
plot(zTravelled, -rref_fitVvsy, 's')
plot(xMean, -rref_fitVvsy, 's')
plot(zMean, -rref_fitVvsy, 's')
plot(log(abs(xMean)), log(-rref_fitVvsy), 's')
plot(log(abs(zMean)), log(-rref_fitVvsy), 's')
% Trajectories' views
close all;
trajPlot = figure; hold on;
for ct=1:length(data)
    xyz = data(ct).dataTrackExcerpt.filteredState(:,[2 3 4]); % the complete trajectory
    
    if ~any(xyz(:,2)>-5e-3)
%         plot3(xyz(:,1), xyz(:,2), xyz(:,3),'Color',[0 97 205]./255); 
        plot3(xyz(:,1), -1*xyz(:,2), xyz(:,3),'Color',[117 153 242]./255); 
%         plot3(xyz(:,1), xyz(:,2), xyz(:,3),'Color',[252,187,161]./255); 
%         plot3(xyz(:,1), xyz(:,2), xyz(:,3));
    else
        disp(true)
    end
end

for ct=1:length(all_data4rrefEstimates)
    xyz4rref = all_data4rrefEstimates(ct).state4rrefEstimate(:,[2 3 4]); % trajectory used for estimating r*
%     plot3(xyz4rref(:,1), xyz4rref(:,2), xyz4rref(:,3),'Color',[215 48 39]./255);
    plot3(xyz4rref(:,1), -1*xyz4rref(:,2), xyz4rref(:,3),'Color',[239 38 38]./255, 'Linewidth', 1);
end
% % Landing Disc
% [X,Y,Z] = cylinder(treatment.landingDiscs(1).radius);
radius = treatment.landingDiscs(1).radius;

figure(trajPlot);
% h=mesh(X,Y,Z,'facecolor',[1 0 0]); % draw landing disc
plot3([-radius; radius], [0 0], [0 0], 'LineWidth', 2, 'Color', [83 83 83]./255);
plot3([0 0], [0 0], [-radius; radius], 'LineWidth', 2, 'Color', [83 83 83]./255);
zlabel('z (m)', 'FontSize', 14);
ylabel('y (m)', 'FontSize', 14);
xlabel('x (m)', 'FontSize', 14);
set(gca, 'FontSize', 16);
view(0,90);
view(-90,0);


% % % Create plot of all V vs y curves used for r*
states = vertcat(all_data4rrefEstimates(:).state4rrefEstimate);
b1 = robustfit(states(:,3),states(:,6));
c1 = fit(states(:,3),states(:,6),'poly1');
d1 = fitlm(-1*states(:,3),states(:,6),'Intercept',false);
figure;
hold on;
% plot(-1*states(:,3),states(:,6),'r');
plot(-1*states(:,3),states(:,6),'o','MarkerSize',2,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6]);
plot(-1*states(:,3),d1.Coefficients.Estimate*(-1*states(:,3)),'b','LineWidth',2);
% plot(states(:,3),b1(1)+b1(2)*states(:,3),'k','LineWidth',2);
ylabel('V_{gy} (m/s)', 'FontSize', 16);
xlabel('y from the platform (m)', 'FontSize', 16);
set(gca, 'FontSize', 18); grid on;

modelfun = @(b,x)(b(1)*x.^(b(2)+1));
% beta0 = [d1.Coefficients.Estimate;1];
beta0 = [d1.Coefficients.Estimate;0.5];
beta = nlinfit(-1*states(:,3),states(:,6),modelfun,beta0)
plot(0:0.001:max(-states(:,3)),modelfun([1.461 -0.415],0:0.001:max(-states(:,3))),'k','LineWidth',2);

% plotting log(-y) vs log(v) for all points to check if it is a straight line or not
figure;
d1=fitlm(-states(:,3),states(:,6))
d2=fitlm(log(-states(:,3)),log(states(:,6)))
d1.Rsquared.Ordinary
d2.Rsquared.Ordinary
hold on;
plot(log(-states(:,3)),log(states(:,6)),'o','MarkerSize',2,...
    'MarkerEdgeColor',[215 48 39]./255,...
    'MarkerFaceColor',[252,146,114]./255);

% Plot all trajectories together
figure;
hold on;
for ct=1:length(data)
    plot(-data(ct).dataTrackExcerpt.filteredState(:,3), data(ct).dataTrackExcerpt.filteredState(:,6),'b.')
end
xlim([0 Inf])
ylim([0 Inf])
%% Check Rsquare computations
for ct=1:length(all_data4rrefEstimates)
    state = all_data4rrefEstimates(ct).state4rrefEstimate;
    % Compute Rsquare manually
    r2_calc = calculateR2(state(:,6), -all_data4rrefEstimates(ct).rref*-state(:,3), true);
    
    % Display computed and manually calculated Rsquare values
    disp([r2_calc all_data4rrefEstimates(ct).Rsquared]);
end
%% Check RMSE computations
for ct=1:length(all_data4rrefEstimates)
    state = all_data4rrefEstimates(ct).state4rrefEstimate;
    % Compute Rsquare manually
    RMSE_calc = calculateRMSE(state(:,6), -all_data4rrefEstimates(ct).rref*-state(:,3), true);
    
    % Display computed and manually calculated SE values
    disp([ct RMSE_calc all_data4rrefEstimates(ct).model.RMSE]);
%     disp([SE_calc all_data4rrefEstimates(ct).model.SSE]);

end
%%
% Creating V vs y and V/y vs y with multiple r*

for ct=1:length(data_multiple_rrefs)
%     plotHandles = plotMultiple_rrefs(data_multiple_rrefs(ct).dataPerTrackExcerpt, data_multiple_rrefs(ct).dataTrackExcerpt);
%     plotHandles = plotMultiple_rrefs_2(data_multiple_rrefs(ct).dataPerTrackExcerpt, data_multiple_rrefs(ct).dataTrackExcerpt);
    plotHandles = plotMultiple_rrefs_3(data_multiple_rrefs(ct).dataPerTrackExcerpt, data_multiple_rrefs(ct).dataTrackExcerpt);
    pause(0.1);
%     if true && ~isempty(plotHandles)
%         % Resizing the figures
%         for i=1:length(plotHandles)
%             plotHandles(i).Position(3) = 680;
%             plotHandles(i).Position(4) = 545;
%             
% %             if i==1
% %                 figureName = ['rrefEstimate_' num2str(ct) ...
% %                     '.png'];
% %             end
% %             
% %             saveas(plotHandles(i), fullfile(DirPlots_treatment, figureName) ,'png');
%         end
%         
% %         close(plotHandles);
%     end
end

%%
% Creating V vs y and V/y vs y with single r*
for ct=1:length(data_single_rref)
    plotHandles = plotMultiple_rrefs_2(data_single_rref(ct).dataPerTrackExcerpt, data_single_rref(ct).dataTrackExcerpt);
    pause(0.1);
%     if true && ~isempty(plotHandles)
%         % Resizing the figures
%         for i=1:length(plotHandles)
%             plotHandles(i).Position(3) = 680;
%             plotHandles(i).Position(4) = 545;
%             
%             if i==1
%                 figureName = ['One_rrefEstimate_' num2str(ct) ...
%                     '.png'];
%             end
%             
%             saveas(plotHandles(i), fullfile(DirPlots_treatment, figureName) ,'png');
%         end
%         
%         close(plotHandles);
%     end
end
% Good ones for the article - 7, 11, 34, 45, ... , 269
%%
% Plotting 45 and 269 for plotting
ct_selected = [45, 269];
for ct1=1:length(ct_selected)
    ct = ct_selected(ct1);
    plotHandles = plotMultiple_rrefs_2(data_single_rref(ct).dataPerTrackExcerpt, data_single_rref(ct).dataTrackExcerpt);
end
%% Analyse min velocity during transition from one r* to another
i = 0; % transition count
for ct=1:length(data_multiple_rrefs)
    rrefs = -[data_multiple_rrefs(ct).dataPerTrackExcerpt.rref];
    [sorted_rrefs, indx] = sort(rrefs);
    for ct1=1:length(rrefs)-1 % for each transition
        i = i+1;
        
        % find min vel_y in transition
        transition(i).t_start = data_multiple_rrefs(ct).dataPerTrackExcerpt(indx(ct1)).state4rrefEstimate(end,1);
        transition(i).t_end = data_multiple_rrefs(ct).dataPerTrackExcerpt(indx(ct1+1)).state4rrefEstimate(1,1);
        assert(transition(i).t_start <= transition(i).t_end);
        transition(i).ct = ct;
        transition(i).from = rrefs(indx(ct1));
        transition(i).to = rrefs(indx(ct1+1));
        transition(i).vmin = min(data_multiple_rrefs(ct).dataTrackExcerpt.filteredState(transition(i).t_start <= data_multiple_rrefs(ct).dataTrackExcerpt.filteredState(:,1) & ...
           data_multiple_rrefs(ct).dataTrackExcerpt.filteredState(:,1) <= transition(i).t_end, 6));
    end
%     keyboard;
end
figure;
histogram([transition.vmin]);
xlabel('......', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);

%% Analysing data_single_rref for hovering phases
before_hover = false(length(data_single_rref),1);
after_hover = false(length(data_single_rref),1);
y_before_hover = zeros(length(data_single_rref),1); % this sdoes NOT work
y_after_hover = zeros(length(data_single_rref),1);
for ct=1:length(data_single_rref)
    minTime = min(data_single_rref(ct).dataPerTrackExcerpt.state4rrefEstimate(:,1));
    maxTime = max(data_single_rref(ct).dataPerTrackExcerpt.state4rrefEstimate(:,1));
    
    dummy = data_single_rref(ct).dataTrackExcerpt.filteredState;
    if any(dummy(dummy(:,1) < minTime,6)  < 0.05)
        before_hover(ct) = true;
        indx = find(dummy(dummy(:,1) < minTime,6)  < 0.05, 1, 'last');
        y_before_hover = dummy(indx, 3);
    end
    if any(dummy(dummy(:,1) > maxTime,6) < 0.05) 
        after_hover(ct) = true;
        indx = find(dummy(dummy(:,1) > maxTime,6) < 0.05, 1);
        y_after_hover = dummy(indx, 3);
    end
end
[sum(before_hover), length(data_single_rref)]
[sum(after_hover), length(data_single_rref)]
[sum(before_hover & after_hover), length(data_single_rref)]
[sum(before_hover | after_hover), length(data_single_rref)]
%% Display track with multiple r* that have a video associated with it
for ct=1:length(data_multiple_rrefs)
    if ~isempty(data_multiple_rrefs(ct).videoInfo)
        disp(ct);
    end
end
%% Selecting movie corresponding to data_multiple_rrefs(29)
% Find frames in the movie between y=-0.25 to y=0
% This y away from the approaching bumblebee (landing disc coordinate
% frames are different in code and in the article)
ct_selected=29;

% Plotting tau function
data_multiple_rrefs(ct_selected).track.plotDataLDF_Time2_forArticle(data_multiple_rrefs(ct_selected).dataTrackExcerpt, t_start);

% Finding time instant that is y=-0.25
t_start = interp1(data_multiple_rrefs(ct_selected).dataTrackExcerpt.filteredState(:,3),data_multiple_rrefs(ct_selected).dataTrackExcerpt.filteredState(:,1),-0.25);
t_end = data_multiple_rrefs(ct_selected).dataTrackExcerpt.filteredState(end,1);

% Plotting variables with time
plotHandle = data_multiple_rrefs(ct_selected).track.plotDataLDF_Time_forArticle(data_multiple_rrefs(ct_selected).dataTrackExcerpt, t_start);
figure(plotHandle);
subplot(3,1,1); grid off
ylim([0 0.3]);
yticks([0:0.1:0.3]);
subplot(3,1,2); grid off
ylim([-0.2 0.4]);
yticks([-0.2:0.2:0.4]);
subplot(3,1,3); grid off
ylim([0 3]);
yticks([0:1:3]);
% xlim([0 Inf]);
xlim([t_start-t_end 0]);

% Plotting variables with distance
plotHandle = plotMultiple_rrefs_3(data_multiple_rrefs(ct_selected).dataPerTrackExcerpt, data_multiple_rrefs(ct_selected).dataTrackExcerpt);
figure(plotHandle);
subplot(2,1,1); grid off
ylim([-0.2 0.4]);
yticks([-0.2:0.2:0.4]);
subplot(2,1,2); grid off
ylim([0 2]);
yticks([0:1:2]);

% load movie data
[videodata, timestamps] = fmf_read(fullfile(data_multiple_rrefs(ct_selected).videoInfo(1).folder, data_multiple_rrefs(ct_selected).videoInfo(1).name));

% Extract and save frames between t_start and t_end
images = videodata(:,:,t_start <= timestamps & timestamps <= t_end);

% Save video
Directory = '/home/reken001/Pulkit/MATLAB/graphs';
v = VideoWriter(fullfile(Directory, [data_multiple_rrefs(ct_selected).videoInfo(1).name(1:end-4) '.avi']),'Uncompressed AVI');
v.FrameRate = 30;
open(v);
for i=1:size(images,3)
    writeVideo(v,images(:,:,i));
end
close(v);

% Save images
Directory = '/home/reken001/Pulkit/MATLAB/graphs';
filename = fullfile(Directory, data_multiple_rrefs(ct_selected).videoInfo(1).name(1:end-4));
mkdir(filename);
for i=1:15:size(images,3)
    imwrite(images(:,:,i), fullfile(filename, [num2str(i,'%03.f') '.png']));
end
%%%%%%%%% Next movie %%%%%%%%%%%%%%%%%%%%%%
% Find frames in the movie between y=-0.25 to y=0
ct_selected=40;

% Finding time instant that is y=-0.25
t_start = interp1(data_multiple_rrefs(ct_selected).dataTrackExcerpt.filteredState(:,3),data_multiple_rrefs(ct_selected).dataTrackExcerpt.filteredState(:,1),-0.25);
t_end = data_multiple_rrefs(ct_selected).dataTrackExcerpt.filteredState(end,1);

% load movie data
[videodata, timestamps] = fmf_read(fullfile(data_multiple_rrefs(ct_selected).videoInfo(1).folder, data_multiple_rrefs(ct_selected).videoInfo(1).name));

% Extract and save frames between t_start and t_end
images = videodata(:,:,t_start <= timestamps & timestamps <= t_end);

% Save video and images
Directory = '/home/reken001/Pulkit/MATLAB/graphs';
v = VideoWriter(fullfile(Directory, [data_multiple_rrefs(ct_selected).videoInfo(1).name(1:end-4) '.avi']),'Uncompressed AVI');
v.FrameRate = 30;
open(v);
for i=1:size(images,3)
    writeVideo(v,images(:,:,i));
end
close(v);
%%
close all;
chosen = [24 26 29 33];
for ct1=1:length(chosen)
    ct = chosen(ct1);
%     plotHandles = plotMultiple_rrefs(data_multiple_rrefs(ct).dataPerTrackExcerpt, data_multiple_rrefs(ct).dataTrackExcerpt);
    plotHandles = plotMultiple_rrefs_2(data_multiple_rrefs(ct).dataPerTrackExcerpt, data_multiple_rrefs(ct).dataTrackExcerpt);
    pause(0.1);
end
%% To display results
pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
behaviour = {'rising','constant','sleeping'};

for ct_pattern = 1%1:length(pattern)
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
            
            % % REMOVE THIS CODE IN THE FINAL VERSION OF THE CODE
%             appData_filename = ['GUIDE_appData_Blindtracks_' pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour} '.mat'];
%             if exist(fullfile(DataDir, appData_filename), 'file')
%                 load(fullfile(DataDir, appData_filename));
%             end
            
            for ct_treatment=1:length(relevantTreatments)
                treatment = relevantTreatments(ct_treatment);
                if ~isempty(treatment.landingTracks)
                    disp(['Into, day: ' num2str(treatment.datenum) ...
                          ', pattern: ' treatment.pattern ...
                          ', light: ' treatment.light ...
                          ', behaviour: ' behaviour{ct_behaviour}]);
                    
                      for ct_track=1:length(treatment.landingTracks) % for each landing track
                          for ct_excerpt=1:length(treatment.landingTracks(ct_track).state_LDF) % for each track excerpt
                              
                              if ct_track >= 14
                                  error('Mid script reached, this error message was made on purpose');
                              end
                              
                              if treatment.landingTracks(ct_track).DataGUI(ct_excerpt).isExcerptUseful ... 
                                 && abs(treatment.landingTracks(ct_track).DataGUI(ct_excerpt).y_mode1(1,1)) >= 1e-3
                                    
                                  disp(['Analysing, track: ' num2str(ct_track) ...
                                        ', excerpt: ' num2str(ct_excerpt)]);
                                  track = treatment.landingTracks(ct_track);
                                  
                                  
                                  
                                  disp([track.track_subset_sysID(ct_excerpt).param_estimation.xOpt]);
                                  disp([track.track_subset_sysID(ct_excerpt).param_estimation.cost]);
%                                   keyboard;
                
                              end
                          end
                      end
                end
            end
                              
            
        end
    end
end
%% Parameter sweep
% clc; clear;
% close all; 

% Inputs
DataDir = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/';
inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks.mat';
load(inputFile);

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
behaviour = {'rising','constant','sleeping'};
% 
% pattern = {'spokes'};
% light = {'high'};
% behaviour = {'constant'};

save_rrefEstimatePlots = true;

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
            
            % % REMOVE THIS CODE IN THE FINAL VERSION OF THE CODE
%             appData_filename = ['GUIDE_appData_Blindtracks_' pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour} '.mat'];
%             if exist(fullfile(DataDir, appData_filename), 'file')
%                 load(fullfile(DataDir, appData_filename));
%             end
            
            for ct_treatment=1:length(relevantTreatments)
                treatment = relevantTreatments(ct_treatment);
                if ~isempty(treatment.landingTracks)
                    disp(['Into, day: ' num2str(treatment.datenum) ...
                          ', pattern: ' treatment.pattern ...
                          ', light: ' treatment.light ...
                          ', behaviour: ' behaviour{ct_behaviour}]);
                    
                      for ct_track=1:length(treatment.landingTracks) % for each landing track
                          for ct_excerpt=1:length(treatment.landingTracks(ct_track).state_LDF) % for each track excerpt
                              treatment.landingTracks(ct_track).track_subset_sysID(ct_excerpt) = trackForLandingModel();
                              
                              if ct_track >= 60
                                  error('Mid script reached, this error message was made on purpose');
                              end
                              
                              if treatment.landingTracks(ct_track).DataGUI(ct_excerpt).isExcerptUseful ... 
                                 && abs(treatment.landingTracks(ct_track).DataGUI(ct_excerpt).y_mode1(1,1)) >= 1e-3
                                    
                                  disp(['Analysing, track: ' num2str(ct_track) ...
                                        ', excerpt: ' num2str(ct_excerpt)]);
                                  track = treatment.landingTracks(ct_track);
                                  excerpt = treatment.landingTracks(ct_track).state_LDF(ct_excerpt);
                                  DataGUI = treatment.landingTracks(ct_track).DataGUI(ct_excerpt);
                                  
                                  % Estimate rref
                                  y_start_indx = find(abs(excerpt.filteredState(:,3)-DataGUI.y_mode1(1,1))<=5e-3 & abs(excerpt.filteredState(:,6)-DataGUI.Vgy_mode1(1,1))<=1e-1,1);
                                  y_end_indx = find(abs(excerpt.filteredState(:,3)-DataGUI.y_mode1(1,2))<=5e-3 & abs(excerpt.filteredState(:,6)-DataGUI.Vgy_mode1(1,2))<=1e-1,1);
                                  if isempty(y_start_indx) || isempty(y_end_indx)
                                      warning('start and end indx not found! :/');
                                      continue;
                                  end
                                  
                                  track.track_subset_sysID(ct_excerpt).state4rrefEstimate = excerpt.filteredState(y_start_indx:y_end_indx, :);
                                  plotHandle = track.track_subset_sysID(ct_excerpt).estimate_rref(save_rrefEstimatePlots, excerpt);
                                  
                                  if save_rrefEstimatePlots % Saving plots to the hard disk
                                      
                                  end
                                  
%                                   keyboard;
                                  
                                  % Estimate other parameters
                                  
                                  
                                  
                                  % Finding iodata for system identification
                                  % for y(1,1): find the one that is closest to the
                                  % marked point
                                  y_start_indx = find(abs(excerpt.filteredState(:,3)-DataGUI.y(1,1))<=5e-3 & abs(excerpt.filteredState(:,6)-DataGUI.Vgy(1,1))<=1e-1,1);
                                  y_end_indx = find(abs(excerpt.filteredState(:,3)-DataGUI.y(1,2))<=5e-3 & abs(excerpt.filteredState(:,6)-DataGUI.Vgy(1,2))<=1e-1,1);
                                  y_neg_indx = find(excerpt.filteredState(y_start_indx:y_end_indx,1) - excerpt.filteredState(y_start_indx,1) > 0.06, 1, 'first');
                                  if isempty(y_start_indx) || isempty(y_end_indx) || isempty(y_neg_indx)
                                      warning('start, end or neg indx not found! :/');
                                      continue;
                                  end
%                                   y_start_indx = y_start_indx + 10;
%                                   y_end_indx = y_end_indx - 10;
                                  
                                  track.track_subset_sysID(ct_excerpt).state = excerpt.filteredState(y_start_indx+y_neg_indx:y_end_indx, :);
                                  
%                                   y_neg_indx = find(excerpt.filteredState(1:y_start_indx,1) - excerpt.filteredState(y_start_indx,1) < -0.1, 1, 'last');                                  
                                  
                                  track.track_subset_sysID(ct_excerpt).negTime_state = excerpt.filteredState(y_start_indx:y_start_indx+y_neg_indx, :);
          
                                  % Plotting data for visual verification
                                  figure;
                                  plot(track.track_subset_sysID(ct_excerpt).state(:,1)-track.track_subset_sysID(ct_excerpt).state(1,1), track.track_subset_sysID(ct_excerpt).state(:,6)./track.track_subset_sysID(ct_excerpt).state(:,3),'b');
                                  hold on;
                                  plot(track.track_subset_sysID(ct_excerpt).negTime_state(:,1)-track.track_subset_sysID(ct_excerpt).state(1,1), track.track_subset_sysID(ct_excerpt).negTime_state(:,6)./track.track_subset_sysID(ct_excerpt).negTime_state(:,3),'r');
                                  xlabel('time'); ylabel('V_{gy}/y');
                                  
%                                   % Learning parameters
%                                   outputSignalForEstimation = 'ay';
                                  plotHandles = estimateLandingParameters(track.track_subset_sysID(ct_excerpt));
                                    
                                 
                                  
                                  disp([track.track_subset_sysID(ct_excerpt).param_estimation.xOpt]);
                                  disp([track.track_subset_sysID(ct_excerpt).param_estimation.cost]);
                                  pause(0.1);
%                                   keyboard;
                
                              end
                          end
                      end
                end
            end
                              
            
        end
    end
end
save(inputFile,'treatments');

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

function plotHandles = plotMultiple_rrefs(dataPerTrackExcerpt, dataTrackExcerpt)
    % Plots V vs y and V/y vs y highling change in r* within the same track
    % excerpt
    
    % dataPerTrackExcerpt - array of data4rrefEstimate instances
    % dataTrackExcerpt - 1x1 filteredState_BlindLandingtrack instance
    
    % Find track data that is not used for r* estimates, but needs to be
    % plotted.
    data4rref = vertcat(dataPerTrackExcerpt.state4rrefEstimate);
%     ymin = min(data4rref(:,3));
%     [~, indx_ymin] = min(abs(dataTrackExcerpt.filteredState(:,3)-ymin));
%     indx_start = find(dataTrackExcerpt.filteredState(1:indx_ymin,3)-ymin<-0.05,1,'last');
    
    ymin = min(data4rref(:,3));
    [~, indx_ymin] = min(abs(dataTrackExcerpt.filteredState(:,3)-ymin));
    time_ymin = dataTrackExcerpt.filteredState(indx_ymin,1);
%     % Start from 0.06s data earlier than dataTrackExcerpt.filteredState(indx_ymax,1)
    indx_start = find(dataTrackExcerpt.filteredState(:,1)-time_ymin<-0.06,1);
    if isempty(indx_start)
        indx_start = 1;
    end
    complete_state = dataTrackExcerpt.filteredState(indx_start:end,:);
    
    plotHandles = figure;
    p1 = subplot(2,1,1); hold on;
    plot(complete_state(:,3),complete_state(:,6),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
    ylabel('V_{gy} (m/s)', 'FontSize', 16);
    set(gca, 'FontSize', 16); grid on;
    title(['r* : ' num2str([dataPerTrackExcerpt.rref],3)], 'FontSize', 16);
    ylim([-0.1 p1.YLim(2)]);
%     yticks([-0.2:0.2:0.6]);
%     xticks([p1.XLim(1):0.1:0]);
    
    
    p2 = subplot(2,1,2); hold on;
    plot(complete_state(:,3),complete_state(:,6)./complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
    ylabel('Rate of expansion (1/s)', 'FontSize', 16);
    ylabel('V_{gy} / y (1/s)', 'FontSize', 16);
    xlabel('y (m)', 'FontSize', 16);
    set(gca, 'FontSize', 16); grid on;
    ylim([ceil(min(data4rref(:,6)./data4rref(:,3))-1) 1]);
%     yticks(ceil(min(data4rref(:,6)./data4rref(:,3))-1):1:1);
%     xticks(p1.XLim(1):0.1:0);
    
    linkaxes([p1, p2],'x');
    if abs(floor((min(data4rref(:,3))/0.05))*0.05 - min(data4rref(:,3))) > 0.025
        xmin = floor((min(data4rref(:,3))/0.05))*0.05;
    else
        xmin = floor((min(data4rref(:,3))/0.05))*0.05-0.05;
    end
    xlim([xmin 0]);
    
    for ct=1:length(dataPerTrackExcerpt)
        obj = dataPerTrackExcerpt(ct); 
        state_subset = obj.state4rrefEstimate;
        figure(plotHandles(1))
        subplot(2,1,1);
        plot(state_subset(:,3),state_subset(:,6),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
        plot(state_subset(:,3),obj.rref*state_subset(:,3),'LineWidth',2,'Color',[69 117 180]./255);
        plot([0 state_subset(end,3)],[0 obj.rref*state_subset(end,3)],'--','LineWidth',2,'Color',[69 117 180]./255);
        
        subplot(2,1,2);
        plot(state_subset(:,3),state_subset(:,6)./state_subset(:,3),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
        %                      plot(state_subset(:,3),c(1)+c(2)*state_subset(:,6)./state_subset(:,3),'b','LineWidth',2);
        plot(state_subset([1,end],3),[obj.rref, obj.rref],'LineWidth',2,'Color',[69 117 180]./255);
        plot([p1.XLim(1) state_subset(1,3)],[obj.rref obj.rref],'--','LineWidth',2,'Color',[69 117 180]./255);
    end
    
    
end

function plotHandles = plotMultiple_rrefs_2(dataPerTrackExcerpt, dataTrackExcerpt)
    % Plots V vs y and V/y vs y highling change in r* within the same track
    % excerpt with positive V and positive y
    
    % dataPerTrackExcerpt - array of data4rrefEstimate instances
    % dataTrackExcerpt - 1x1 filteredState_BlindLandingtrack instance
    
    % Find track data that is not used for r* estimates, but needs to be
    % plotted.
    data4rref = vertcat(dataPerTrackExcerpt.state4rrefEstimate);
%     ymin = min(data4rref(:,3));
%     [~, indx_ymin] = min(abs(dataTrackExcerpt.filteredState(:,3)-ymin));
%     indx_start = find(dataTrackExcerpt.filteredState(1:indx_ymin,3)-ymin<-0.05,1,'last');
    
    ymin = min(data4rref(:,3));
    [~, indx_ymin] = min(abs(dataTrackExcerpt.filteredState(:,3)-ymin));
    time_ymin = dataTrackExcerpt.filteredState(indx_ymin,1);
%     % Start from 0.06s data earlier than dataTrackExcerpt.filteredState(indx_ymax,1)
    indx_start = find(dataTrackExcerpt.filteredState(:,1)-time_ymin<-0.06,1);
    if isempty(indx_start)
        indx_start = 1;
    end
    complete_state = dataTrackExcerpt.filteredState(indx_start:end,:);
    
    plotHandles = figure;
    p1 = subplot(2,1,1); hold on;
    plot(-complete_state(:,3),complete_state(:,6),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
    ylabel('V (m/s)', 'FontSize', 16);
    set(gca, 'FontSize', 16); %grid on;
    title(['r* : ' num2str(-1*[dataPerTrackExcerpt.rref],3)], 'FontSize', 16);
%     ylim([-0.2 p1.YLim(2)]);
    ylim([-0.2 0.6]);
    yticks([-0.2:0.2:0.6]);
%     xticks([p1.XLim(1):0.1:0]);
    
    
    p2 = subplot(2,1,2); hold on;
%     plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
    plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
    ylabel('Rate of expansion (1/s)', 'FontSize', 16);
    ylabel('r (1/s)', 'FontSize', 16);
%     ylabel('V_{gy} / y (1/s)', 'FontSize', 16);
    xlabel('y (m)', 'FontSize', 16);
    set(gca, 'FontSize', 16); %grid on;
    ylim([0 ceil(-1*ceil(min(data4rref(:,6)./data4rref(:,3))-1))]);
    yticks(0:1:ceil(-1*ceil(min(data4rref(:,6)./data4rref(:,3))-1)));
%     xticks(p1.XLim(1):0.1:0);
    
    linkaxes([p1, p2],'x');
    if abs(floor((min(data4rref(:,3))/0.05))*0.05 - min(data4rref(:,3))) > 0.025
        xmin = floor((min(data4rref(:,3))/0.05))*0.05;
    else
        xmin = floor((min(data4rref(:,3))/0.05))*0.05-0.05;
    end
%     xlim([0 -xmin]);
    xlim([0 Inf]);
    
    for ct=1:length(dataPerTrackExcerpt)
        obj = dataPerTrackExcerpt(ct); 
        state_subset = obj.state4rrefEstimate;
        figure(plotHandles(1))
        subplot(2,1,1);
        plot(-state_subset(:,3),state_subset(:,6),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
        plot(-state_subset(:,3),-obj.rref*-state_subset(:,3),'LineWidth',2,'Color',[69 117 180]./255);
        plot([0 -state_subset(end,3)],[0 -obj.rref*-state_subset(end,3)],'--','LineWidth',2,'Color',[69 117 180]./255);
        
        subplot(2,1,2);
        plot(-state_subset(:,3),state_subset(:,6)./-state_subset(:,3),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
        %                      plot(state_subset(:,3),c(1)+c(2)*state_subset(:,6)./state_subset(:,3),'b','LineWidth',2);
        plot(-state_subset([1,end],3),[-obj.rref, -obj.rref],'LineWidth',2,'Color',[69 117 180]./255);
        plot([p1.XLim(1) -state_subset(1,3)],[-obj.rref -obj.rref],'--','LineWidth',2,'Color',[69 117 180]./255);
    end
    
end

function plotHandles = plotMultiple_rrefs_3(dataPerTrackExcerpt, dataTrackExcerpt)
    % Plots V vs y and V/y vs y highling change in r* within the same track
    % excerpt with positive V and positive y
    % Also plots acceleration in y with y
    
    % dataPerTrackExcerpt - array of data4rrefEstimate instances
    % dataTrackExcerpt - 1x1 filteredState_BlindLandingtrack instance
    
    % Find track data that is not used for r* estimates, but needs to be
    % plotted.
    data4rref = vertcat(dataPerTrackExcerpt.state4rrefEstimate);
%     ymin = min(data4rref(:,3));
%     [~, indx_ymin] = min(abs(dataTrackExcerpt.filteredState(:,3)-ymin));
%     indx_start = find(dataTrackExcerpt.filteredState(1:indx_ymin,3)-ymin<-0.05,1,'last');
    
    ymin = min(data4rref(:,3));
    [~, indx_ymin] = min(abs(dataTrackExcerpt.filteredState(:,3)-ymin));
    time_ymin = dataTrackExcerpt.filteredState(indx_ymin,1);
%     % Start from 0.06s data earlier than dataTrackExcerpt.filteredState(indx_ymax,1)
    indx_start = find(dataTrackExcerpt.filteredState(:,1)-time_ymin<-0.06,1);
    if isempty(indx_start)
        indx_start = 1;
    end
    complete_state = dataTrackExcerpt.filteredState(indx_start:end,:);
    
    plotHandles = figure;
    p1 = subplot(3,1,1); hold on;
    plot(-complete_state(:,3),complete_state(:,6),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
    ylabel('V (m/s)', 'FontSize', 16);
    set(gca, 'FontSize', 16); %grid on;
    title(['r* : ' num2str(-1*[dataPerTrackExcerpt.rref],3)], 'FontSize', 16);
%     ylim([-0.2 p1.YLim(2)]);
    ylim([-0.2 0.6]);
    yticks([-0.2:0.2:0.6]);
%     xticks([p1.XLim(1):0.1:0]);
    
    
    p2 = subplot(3,1,2); hold on;
%     plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
    plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
    ylabel('Rate of expansion (1/s)', 'FontSize', 16);
    ylabel('r (1/s)', 'FontSize', 16);
%     ylabel('V_{gy} / y (1/s)', 'FontSize', 16);
%     xlabel('y (m)', 'FontSize', 16);
    set(gca, 'FontSize', 16); %grid on;
    ylim([0 ceil(-1*ceil(min(data4rref(:,6)./data4rref(:,3))-1))]);
    yticks(0:1:ceil(-1*ceil(min(data4rref(:,6)./data4rref(:,3))-1)));
%     xticks(p1.XLim(1):0.1:0);
    
    
    if abs(floor((min(data4rref(:,3))/0.05))*0.05 - min(data4rref(:,3))) > 0.025
        xmin = floor((min(data4rref(:,3))/0.05))*0.05;
    else
        xmin = floor((min(data4rref(:,3))/0.05))*0.05-0.05;
    end
%     xlim([0 -xmin]);
    xlim([0 Inf]);
    
    p3 = subplot(3,1,3); hold on;
%     plot(-complete_state(:,3),complete_state(:,6)./-complete_state(:,3),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
    plot(-complete_state(:,3),complete_state(:,9),'.','MarkerSize',10,'MarkerFaceColor',[252,187,161]./255, 'MarkerEdgeColor',[252,187,161]./255');
    ylabel('acc (m/s^2)', 'FontSize', 16);
    xlabel('y (m)', 'FontSize', 16);
    
    linkaxes([p1, p2, p3],'x');
    
    for ct=1:length(dataPerTrackExcerpt)
        obj = dataPerTrackExcerpt(ct); 
        state_subset = obj.state4rrefEstimate;
        figure(plotHandles(1))
        subplot(3,1,1);
        plot(-state_subset(:,3),state_subset(:,6),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
        plot(-state_subset(:,3),-obj.rref*-state_subset(:,3),'LineWidth',2,'Color',[69 117 180]./255);
        plot([0 -state_subset(end,3)],[0 -obj.rref*-state_subset(end,3)],'--','LineWidth',2,'Color',[69 117 180]./255);
        
        subplot(3,1,2);
        plot(-state_subset(:,3),state_subset(:,6)./-state_subset(:,3),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
        %                      plot(state_subset(:,3),c(1)+c(2)*state_subset(:,6)./state_subset(:,3),'b','LineWidth',2);
        plot(-state_subset([1,end],3),[-obj.rref, -obj.rref],'LineWidth',2,'Color',[69 117 180]./255);
        plot([p1.XLim(1) -state_subset(1,3)],[-obj.rref -obj.rref],'--','LineWidth',2,'Color',[69 117 180]./255);
        
        subplot(3,1,3);
        plot(-state_subset(:,3),state_subset(:,9),'.','MarkerSize',12,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
        %                      plot(state_subset(:,3),c(1)+c(2)*state_subset(:,6)./state_subset(:,3),'b','LineWidth',2);
        
    end
    
end


function [data_cmap, edges] = find_cmap(data, cmap)
% Finds color values for each data point in data
% This is required because usually size(data,1) < size(cmap,1)

% Method 1 - Creating N bins of data and allocating each data point to that bin. N is size(cmap,1)-1
edges = round(linspace(min(data),max(data),size(cmap,1)+1),3); % # edges = # color bins + 1
edges(1) = floor(min(data)*1000)/1000;
edges(end) = ceil(max(data)*1000)/1000;
[data_cmap, edges] = discretize(data,edges);
end

function R2 = calculateR2(z,z_est,isRTO)
% calcuateR2 Cacluate R-squared
% R2 = calcuateR2(z,z_est) takes two inputs - The observed data x and its
% estimation z_est (may be from a regression or other model), and then
% compute the R-squared value a measure of goodness of fit. R2 = 0
% corresponds to the worst fit, whereas R2 = 1 corresponds to the best fit.
% 
% isRTO - if model is Regression through origin (RTO) or not (true or
% false)
if isRTO
    R2 = norm(z_est-mean(z))^2/(norm(z-z_est)^2 + norm(z_est-mean(z))^2);
%     R2 = 1-norm(z-z_est)^2/norm(z-mean(z))^2;
%     R2 = norm(z_est)^2/norm(z)^2; 
else
    r = z-z_est;
    normr = norm(r);
    SSE = normr.^2;
    SSR = norm(z_est-mean(z))^2;
    SST = SSE + SSR; %norm(z-mean(z))^2;
    R2 = 1 - SSE/SST;
end

% mean_z = mean(z);
% num = sum((z_est-mean_z).^2);
% den = sum((z-mean_z).^2);
% R2 = num/den;
end

function RMSE = calculateRMSE(z,z_est,isRTO)
    if isRTO
        k = 1;
    else
        k = 2;
    end
%     SSE = norm(z-z_est)^2; 
    RMSE = sqrt(norm(z-z_est)^2/(length(z)-k));
end