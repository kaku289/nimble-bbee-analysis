%% Extract track excerpts using custom GUI

%%
clc; clear; close all;

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


%%
% Inputs
DataDir = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/';
inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A2_manual.mat';
load(inputFile);

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
behaviour = {'rising','constant','sleeping'};
% 
% pattern = {'spokes'};
% light = {'high'};
% behaviour = {'constant'};

treatments = treatments(1:14*8);
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
            
            appData_filename = ['GUIDE_appData_Blindtracks_' pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour} '_A2_manual.mat'];
            if exist(fullfile(DataDir, appData_filename), 'file')
                load(fullfile(DataDir, appData_filename));
            end
            
        end
    end
end

%% Display track numbers which have more than one rref
tracks = [relevantTreatments.landingTracks];
for ct_track = 1:length(tracks)
    for ct_excerpt = 1:length(tracks(ct_track).DataGUI)
        if ~isempty(tracks(ct_track).DataGUI)  && sum(tracks(ct_track).DataGUI(ct_excerpt).y_mode1(:,1)~=0) > 1
            disp(['track: ' num2str(ct_track) ', excerpt: ' num2str(ct_excerpt) ', rref #: ' num2str(sum(tracks(ct_track).DataGUI(ct_excerpt).y_mode1(:,1)~=0))]);
        end
    end
end

% Run 