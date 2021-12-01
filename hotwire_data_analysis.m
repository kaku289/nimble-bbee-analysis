
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
%% Find maximum change in velocity over 11 days

% Inputs
close all; clc;
% clear;

inputFile = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/BlindLandingtracks_A3_LDF_rref.mat';
outputFile = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/BlindLandingtracks_A3_LDF_rref_rrefEntry.mat';
outputFile = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/BlindLandingtracks_A3_LDF_rref_rrefEntry_with_rdot.mat';
% load(inputFile);

DirPlots = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/plots/BlindLandingTracks/rref_entry';
DirPlots = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/plots/BlindLandingTracks/rref_entry_rdotsim';
% DirPlots = '/media/reken001/Disk_08_backup/wind_intensity_experiments/postprocessing/plots/BlindLandingTracks/rref_estimate_3dspeed';
delPreviousPlots = false; % BE CAREFUL - when set to true, all previously saved plots are deleted
savePlots = false;
savePDFs = false;

winds = unique([treatments.wind]);
behaviour = {'rising','constant','sleeping'};

factors = [0.25:0.25:2.5];
% factors = [1];

winds_vel = [0, 0.28, 0.98, 1.77, 2.54, 3.41];
maxdiff = zeros(size(winds));
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
        hwData = [relevantTreatments.hwData];
        meanVoltage = [hwData.meanVoltage];
        
        maxdiff(ct_wind) = (max(meanVoltage) - min(meanVoltage))/min(meanVoltage);
    end
end