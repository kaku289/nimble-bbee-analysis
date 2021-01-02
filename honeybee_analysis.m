%% This file performs a "dirty" analysis on preliminary data sent by Dr. Baird for honeybee approaches

% to add definition of classes being used for this analysis
addpath('./lib/li_analysis');

% to include higher order accurate differentiation function
addpath('./lib/diffxy');

%% Load data
clc; close all;
clear;

dataDir = '/media/reken001/Disk_12/honeybee_experiments/landing_tracks';

% dataFiles = dir(fullfile(dataDir,'*.mat'));
% Look for dataFiles in dataDir and its subfolders
subdir = dir(dataDir);
subdir = subdir([subdir.isdir]);
subdir=subdir(~ismember({subdir.name},{'..'}));
for ct=1:length(subdir)
    filesmat = dir(fullfile(subdir(ct).folder,subdir(ct).name,'*.mat'));
    filescsv = dir(fullfile(subdir(ct).folder,subdir(ct).name,'*.csv'));
    dataFiles{ct} = [filesmat; filescsv];
end

% dataFiles = vertcat(dataFiles{:});

DirPlots = '/media/reken001/Disk_12/honeybee_experiments/plots/rref_estimate';
delPreviousPlots = true; % BE CAREFUL - when set to true, all previously saved plots are deleted
savePlots = true;
savePDFs = false;

% Delete previous plots
if delPreviousPlots && savePlots && exist(fullfile(DirPlots), 'dir')
    rmdir(fullfile(DirPlots),'s');
end

% Create directory if it doesn't exist
if savePlots && ~exist(fullfile(DirPlots), 'dir')
    mkdir(fullfile(DirPlots));
end
DirPlots_treatment = fullfile(DirPlots);

% Defining parameters for rref estimate
min_gap = 14; % the first straight line is between current point and (current point + min_gap)th point
max_gap = 99;
params = [0.53 4.22]; % [sigma_{rmean-c_{r vs y, 0-1}}, sigma_{m_{r vs y, 0-1}}]
factors = [0.25:0.25:2.5];
factors = 1;
time_window = [min_gap max_gap];

for ct_treatment = 1:length(dataFiles)
    for ct=1:length(dataFiles{ct_treatment})
        % loading track data
        if strcmpi(dataFiles{ct_treatment}(ct).name(end-2:end),'csv')
            data = readmatrix(fullfile(dataFiles{ct_treatment}(ct).folder, dataFiles{ct_treatment}(ct).name));
            dataFiles{ct_treatment}(ct).data.x = data(:,1);
            dataFiles{ct_treatment}(ct).data.y = data(:,2);
            dataFiles{ct_treatment}(ct).data.z = data(:,3);
        elseif strcmpi(dataFiles{ct_treatment}(ct).name(end-2:end),'mat')
            dataFiles{ct_treatment}(ct).data = load(fullfile(dataFiles{ct_treatment}(ct).folder, dataFiles{ct_treatment}(ct).name));
        end
        
        % creating landing track object
        landingTracks{ct_treatment}(ct) = BlindLandingtrack();
        N = length(dataFiles{ct_treatment}(ct).data.x);
        rawState = [ct*ones(N,1) [1:N]' [0:N-1]'/400 dataFiles{ct_treatment}(ct).data.x/1000 -dataFiles{ct_treatment}(ct).data.y/1000 dataFiles{ct_treatment}(ct).data.z/1000 ...
            zeros(N, 12)];% 18 columns (obj_id	frame	timestamp	x	y	z	xvel	yvel	zvel P00	P01	P02	P11	P12	P22	P33	P44	P55)
        [~,pattern,~] = fileparts(dataFiles{ct_treatment}(ct).folder);
        landingTracks{ct_treatment}(ct).rawTrack = rawState_BlindLandingtrack(rawState, pattern);
        landingTracks{ct_treatment}(ct).dt = 1/400;
        landingTracks{ct_treatment}(ct).filterRawTracks(20);
        
        landingTracks{ct_treatment}(ct).raw_LDF = landingTracks{ct_treatment}(ct).rawTrack;
        landingTracks{ct_treatment}(ct).state_LDF = landingTracks{ct_treatment}(ct).state;
        
        %     landingTracks(ct).plotDataLDF_Time_rawVSfiltered(1, []);
        %     landingTracks(ct).plotDataLDF_Distance(landingTracks(ct).state_LDF);
        
        landingTracks{ct_treatment}(ct).state_LDF.compute_rref(params, factors, time_window);
        
        if savePlots
            for ct_factor=1:length(factors)
                
                plotHandles = landingTracks{ct_treatment}(ct).state_LDF.plot_rrefs(factors(ct_factor));
                %                                       plotHandles = excerpt.plot_rrefs_with3dspeed(factors(ct_factor));
                %             plotHandles = landingTracks(ct).state_LDF.plot_states();
                if ~isempty(plotHandles)
                    
                    % Resizing the figures
                    for i=1:length(plotHandles)
                        plotHandles(i).Position(3) = 680;
                        plotHandles(i).Position(4) = 545;
                        
                        if i==1
                            figureName = ['fac_' num2str(factors(ct_factor),'%0.2f') '_' ...
                                pattern '_' (dataFiles{ct_treatment}(ct).name(1:end-4)) ...
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
end


%% create 3d trajectory figure
% dataFiles = vertcat(dataFiles{2:end});
trajFig = figure; hold on;
for ct=1:length(dataFiles)
    plot3(dataFiles(ct).data.x, dataFiles(ct).data.y, dataFiles(ct).data.z);
end
xlabel('x')
ylabel('y')
zlabel('z')
view([0 90]);
view([90 0]);