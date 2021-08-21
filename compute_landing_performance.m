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
%% Compute landing performance

% Inputs
close all; clc;
% clear;

if isunix
    inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A1_rref.mat';
elseif ispc
    inputFile = 'D:/light_intensity_experiments/postprocessing/BlindLandingtracks_A1_rref.mat';
end

% load(inputFile);

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
behaviour = {'rising','constant','sleeping'};

y = [0.25 0.15]; % y's between which performance parameters are calculated
clear dummy;

data = struct.empty;
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

            
            %%%%%%%%%%%%%% Each treatment %%%%%%%%%%%
            landingTracks = [relevantTreatments.landingTracks];
            state_LDF = [landingTracks.state_LDF];

            treatment_indx4landingTracks = arrayfun(@(x) x*ones(length(relevantTreatments(x).landingTracks),1),1:length(relevantTreatments), 'UniformOutput', false);
            treatment_indx4landingTracks = vertcat(treatment_indx4landingTracks{:});
            treatment_indx4stateLDF = arrayfun(@(x) treatment_indx4landingTracks(x)*ones(length(landingTracks(x).state_LDF),1),1:length(landingTracks), 'UniformOutput', false);
            treatment_indx4stateLDF = vertcat(treatment_indx4stateLDF{:});

            landingTracks_indx4stateLDF = arrayfun(@(x) x*ones(length(landingTracks(x).state_LDF),1),1:length(landingTracks), 'UniformOutput', false);
            landingTracks_indx4stateLDF = vertcat(landingTracks_indx4stateLDF{:});

            landingDiscs = {relevantTreatments(treatment_indx4stateLDF).landingDiscs};
            hastakeoff = arrayfun(@(x) state_LDF(x).hasTakeoff(landingDiscs{x}),1:length(state_LDF));

            performance_params = arrayfun(@(x) x.compute_landing_performance(y(1),y(2)), state_LDF, 'UniformOutput', false);
            performance_params = vertcat(performance_params{:});

            delta_t = performance_params(:,1);

            indices = ~isnan(delta_t);

            % store data
            dummy.tracks_fac = state_LDF(indices);
            dummy.pattern = (ct_pattern);
            dummy.light = (ct_light);
            dummy.landingTrack = landingTracks(landingTracks_indx4stateLDF(indices));
            dummy.hastakeoff = hastakeoff(indices);
            dummy.delta_t = delta_t(indices);
            data = [data; dummy];
        end
        
    end
end


%%
data_write = [];
approach = 0;
for ct=1:length(data)
    N = length(data(ct).delta_t);
    landingSide = {data(ct).tracks_fac.landingSide};
    isHive = cellfun(@(x) strcmpi(x,'hive'),landingSide);
    side = ones(length(isHive),1);
    side(~isHive) = deal(2);
    
    data_write = [data_write;
                  data(ct).delta_t data(ct).pattern*ones(N,1) data(ct).light*ones(N,1) data(ct).hastakeoff' [approach+1:approach+N]' side [data(ct).landingTrack.datenum]'];
              
    approach = approach + N;
end

writeFile = true;
if isunix
    r_file = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/data_landing_performance_Rstudio.txt';
elseif ispc
    r_file = 'D:/light_intensity_experiments/postprocessing/data_landing_performance_Rstudio.txt';
end

if writeFile
    T = array2table(data_write, ...
        'VariableNames',{'delta_t','pattern','light','hasTakeoff',...
                         'approach','landingSide','day'});
    writetable(T,r_file);
end
