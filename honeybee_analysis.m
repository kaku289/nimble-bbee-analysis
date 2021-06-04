%% This file performs the analysis on data sent by Dr. Baird for honeybee approaches

% to add definition of classes being used for this analysis
addpath('./lib/li_analysis');

% to include higher order accurate differentiation function
addpath('./lib/diffxy');

%% Load data
clc; close all;
clear;

dataDir = '/media/reken001/Disk_12/honeybee_experiments/landing_tracks';
outputFile = '/media/reken001/Disk_12/honeybee_experiments/postprocessing/BlindLandingtracks_A4_LDF_rref.mat';

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
delPreviousPlots = false; % BE CAREFUL - when set to true, all previously saved plots are deleted
savePlots = false;
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
max_gap = 199;
% from honeybee data
params = [0.29 1.02]; % [sigma_{rmean-c_{r vs y, 0-1}}, sigma_{m_{r vs y, 0-1}}]
% from bumblebee data
% params = [0.53 4.22]; % [sigma_{rmean-c_{r vs y, 0-1}}, sigma_{m_{r vs y, 0-1}}]

factors = [0.25:0.25:2.5];
% factors = [1:0.5:2.5];
% factors = 1.5;
time_window = [min_gap max_gap];

% indexes_tracks_to_not_analyse = [18:51 148:158 167:176 189:199]; %due to y-offset 
count = 0;
for ct_treatment = 1:length(dataFiles)
    landingTracks{ct_treatment} = BlindLandingtrack.empty();
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
        
        % Discard last point as it can be erroneous
            dataFiles{ct_treatment}(ct).data.x = dataFiles{ct_treatment}(ct).data.x(1:end-1);
            dataFiles{ct_treatment}(ct).data.y = dataFiles{ct_treatment}(ct).data.y(1:end-1);
            dataFiles{ct_treatment}(ct).data.z = dataFiles{ct_treatment}(ct).data.z(1:end-1);
            
        % Make correction for y-offset (by setting 0,0,0 as the origin
        % of the reference frame (for trajectories with yend>0.005m)
        if -dataFiles{ct_treatment}(ct).data.y(end) > 0.005*1000
            dataFiles{ct_treatment}(ct).data.x = -dataFiles{ct_treatment}(ct).data.x(end,1) + dataFiles{ct_treatment}(ct).data.x;
            dataFiles{ct_treatment}(ct).data.y = -dataFiles{ct_treatment}(ct).data.y(end,1) + dataFiles{ct_treatment}(ct).data.y;
            dataFiles{ct_treatment}(ct).data.z = -dataFiles{ct_treatment}(ct).data.z(end,1) + dataFiles{ct_treatment}(ct).data.z;
        end
        
        % creating landing track object
        landingTracks{ct_treatment}(ct) = BlindLandingtrack();
        N = length(dataFiles{ct_treatment}(ct).data.x);
        rawState = [ct*ones(N,1) [1:N]' [0:N-1]'/400 dataFiles{ct_treatment}(ct).data.x/1000 -dataFiles{ct_treatment}(ct).data.y/1000 dataFiles{ct_treatment}(ct).data.z/1000 ...
            zeros(N, 12)];% 18 columns (obj_id	frame	timestamp	x	y	z	xvel	yvel	zvel P00	P01	P02	P11	P12	P22	P33	P44	P55)
        [~,pattern,~] = fileparts(dataFiles{ct_treatment}(ct).folder);
        landingTracks{ct_treatment}(ct).rawTrack = rawState_BlindLandingtrack(rawState, pattern);
        landingTracks{ct_treatment}(ct).dt = 1/175; %1/400;
        landingTracks{ct_treatment}(ct).filterRawTracks(20);
        
        landingTracks{ct_treatment}(ct).raw_LDF = landingTracks{ct_treatment}(ct).rawTrack;
        landingTracks{ct_treatment}(ct).state_LDF = landingTracks{ct_treatment}(ct).state;
        
        %     landingTracks(ct).plotDataLDF_Time_rawVSfiltered(1, []);
        %     landingTracks(ct).plotDataLDF_Distance(landingTracks(ct).state_LDF);
        
%         landingTracks{ct_treatment}(ct).state_LDF.compute_rref(params, factors, time_window);
        if all(landingTracks{ct_treatment}(ct).state_LDF.filteredState(:,3)<0.06)
            landingTracks{ct_treatment}(ct).state_LDF.compute_rref(params, factors, time_window);
            has_no_yoffset = true;
        else
%             disp(['Skipping track ' num2str(ct+length([landingTracks{1:end-1}])) 'for r* computation due to y-offset']);
            disp(['Skipping track ' num2str(ct) ' from treatment ' num2str(ct_treatment-1) ' for r* computation due to y-offset']);
            has_no_yoffset = false;
            count = count + 1;
        end
        
        % Extracting additional info from file name
        [~,basedir,~] = fileparts(dataFiles{ct_treatment}(ct).folder);
        filename = dataFiles{ct_treatment}(ct).name;
        [datenum, pattern, patternnum, beeID, setID, flightID] = BlindLandingtrack.extractInfo(basedir, filename);
        
        landingTracks{ct_treatment}(ct).filename = filename;
        landingTracks{ct_treatment}(ct).basedir = basedir;
        
        landingTracks{ct_treatment}(ct).datenum = datenum;
        landingTracks{ct_treatment}(ct).pattern = pattern;
        landingTracks{ct_treatment}(ct).patternnum = patternnum;
        landingTracks{ct_treatment}(ct).beeID = beeID;
        landingTracks{ct_treatment}(ct).setID = setID;
        landingTracks{ct_treatment}(ct).flightID = flightID;
        
        
        if savePlots && has_no_yoffset
            for ct_factor=1:length(factors)
                
                
                % Save y,V,A,r vs t curve and V,r vs y curve for Figure 1
%                 track_info = '2810051100 2 05 -b';
%                 if strcmpi(dataFiles{ct_treatment}(ct).name(1:end-4),track_info) && strcmpi(pattern, 'Random Julesz')
%                     t_start = landingTracks{ct_treatment}(ct).state_LDF.filteredState(1,1)
%                     plotHandles = BlindLandingtrack.plotDataLDF_Time3_forArticle(landingTracks{ct_treatment}(ct).state_LDF, ...
%                         t_start)
%                     figure(plotHandles);
%                     subplot(4,1,1); grid off
%                     ylim([0 0.6]);
%                     yticks([0:0.2:0.6]);
%                     subplot(4,1,2); grid off
%                     ylim([-0.2 1.5]);
%                     yticks([0:0.5:1.5]);
%                     subplot(4,1,3); grid off
%                     ylim([-4 4]);
%                     yticks([-4:4:4]);
%                     subplot(4,1,4); grid off
%                     ylim([0 8]);
%                     yticks([0:2:8]);
%                     
%                     plotHandles = landingTracks{ct_treatment}(ct).state_LDF.plot_states();
%                     figure(plotHandles);
%                     subplot(2,1,1); grid off
%                     ylim([0 1.5]);
%                     yticks([0:0.5:1.5]);
%                     subplot(2,1,2); grid off
%                     ylim([0 8]);
%                     yticks([0:2:8]);
%                 else
%                     plotHandles = [];
%                 end
                
                
                
                plotHandles = landingTracks{ct_treatment}(ct).state_LDF.plot_rrefs(factors(ct_factor));
%                 plotHandles = landingTracks{ct_treatment}(ct).state_LDF.plot_rrefs_parallax(factors(ct_factor));
                %                                       plotHandles = excerpt.plot_rrefs_with3dspeed(factors(ct_factor));
                %             plotHandles = landingTracks(ct).state_LDF.plot_states();
                if ~isempty(plotHandles)
                    
                    % Resizing the figures
                    for i=1:length(plotHandles)
                        plotHandles(i).Position(3) = 680;
                        plotHandles(i).Position(4) = 545;
%                         plotHandles(i).Position(4) = 800;
                        
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
keyboard;
save(outputFile, 'landingTracks');


%% Extract rref data for statistical analysis
clc; close all;
% clear;

inputFile = '/media/reken001/Disk_12/honeybee_experiments/postprocessing/BlindLandingtracks_A4_LDF_rref.mat';
% load(inputFile);

landing_tracks = [landingTracks{:}];

data = data4rrefEstimate.empty;
for ct_track = 1:length(landing_tracks)
    data1 = landing_tracks(ct_track).state_LDF.rrefSegments;
    
    % Save additional data
    [data1.day] = deal(landing_tracks(ct_track).datenum);
    [data1.pattern] = deal(landing_tracks(ct_track).pattern);
    [data1.patternnum] = deal(landing_tracks(ct_track).patternnum);
    [data1.setID] = deal(landing_tracks(ct_track).setID);
    [data1.beeID] = deal(landing_tracks(ct_track).beeID);
    [data1.flightID] = deal(landing_tracks(ct_track).flightID);
    
    data = [data data1];
end

% Discard empty intervals
indices = arrayfun(@(x) ~isempty(x.intervals_ti), data);
data = data(indices);

%% Write file for statistical analysis in R
% This file contains data for all the factors
clc;
writeFile = true;
r_file = '/media/reken001/Disk_12/honeybee_experiments/postprocessing/data_all_rref_Rstudio_hbeeParams.txt';
factors = [0.25:0.25:2.5];
data_write = [];
beeID = {};
for ct_factor=1:length(factors)
    factor = factors(ct_factor);
    
    data_fac = data(abs([data.factor]-factor)<1e-6)';
    N = length(data_fac);
    
    % Create nominal and ordinal variables (numeric data)
    approach = arrayfun(@(i) i*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    flightID = arrayfun(@(i) data_fac(i).flightID*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    setID = arrayfun(@(i) data_fac(i).setID*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    day = arrayfun(@(i) data_fac(i).day*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    patternnum = arrayfun(@(i) data_fac(i).patternnum*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    
    for i=1:N
         dummy = arrayfun(@(x) data_fac(i).beeID, 1:size(data_fac(i).intervals_ti,1), 'UniformOutput', false);
         beeID = [beeID; dummy'];
    end
    
%     approach = arrayfun(@(i) i*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
%     side = arrayfun(@(i) data_fac(i).side*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    
    % create other variables
    y = -vertcat(data_fac.ymean_ti);
    r = -vertcat(data_fac.rref_ti);
    v = vertcat(data_fac.vmean_ti);
    
%     speed3d = vertcat(data_fac.speed3d_mean_ti);

    data_write = [data_write; ...
        vertcat(approach{:}) vertcat(flightID{:}) vertcat(setID{:}) vertcat(day{:}) ...
        vertcat(patternnum{:}) y r v factor*ones(size(r,1),1)];
    
    N1 = sum(arrayfun(@(x) size(x.intervals_ti,1)>1, data_fac));
    disp(['# tracks: ' num2str(N) ', # data points: ' num2str(size(r,1)), ...
          ', # tracks with >1 r*: ' num2str(N1)]);
end

% % % expansioncue variable
% % data_write(:,9) = ones(length(data_write),1);
% % data_write(data_write(:,4) == 4 | data_write(:,4) == 11,9) = 0;

if writeFile
    T = array2table(data_write, ...
        'VariableNames',{'approach','flightID','setID','day','patternnum','y','r','v','threshold'});
    T.beeID = beeID;
    writetable(T,r_file);
end


%% Plot y(end) for filteredState
landing_tracks = [landingTracks{:}];
yend = arrayfun(@(x) x.state_LDF.filteredState(end,3),landing_tracks);
figure;
histogram(yend);
xlabel('y a the end', 'FontSize', 16);
ylabel('Occurences', 'FontSize', 16);
set(gca, 'FontSize', 16);

figure; histogram(yend(yend<=0.005));

yendminus1 = arrayfun(@(x) x.state_LDF.filteredState(end-1,3),landing_tracks);
figure;
histogram(yendminus1);
xlabel('y a the end but one point', 'FontSize', 16);
ylabel('Occurences', 'FontSize', 16);
set(gca, 'FontSize', 16);

figure; histogram(yendminus1(yendminus1<=0.005));

yendminus2 = arrayfun(@(x) x.state_LDF.filteredState(end-2,3),landing_tracks);
figure;
histogram(yendminus2);
xlabel('y a the end but two point', 'FontSize', 16);
ylabel('Occurences', 'FontSize', 16);
set(gca, 'FontSize', 16);

figure; histogram(yendminus2(yendminus2<=0.005));

y = abs([yend'-yendminus1' yendminus1'-yendminus2'])*1000;
y(y(:,1)>1)'


%% Plot y(end) for rawState
landing_tracks = [landingTracks{:}];
yend = arrayfun(@(x) x.rawTrack.rawState(end,5),landing_tracks);
figure;
histogram(yend);
xlabel('y at the end', 'FontSize', 16);
ylabel('Occurences', 'FontSize', 16);
set(gca, 'FontSize', 16);

figure; histogram(yend(yend<=0.005));

yendminus1 = arrayfun(@(x) x.rawTrack.rawState(end-1,5),landing_tracks);
figure;
histogram(yendminus1);
xlabel('y at the end but one point', 'FontSize', 16);
ylabel('Occurences', 'FontSize', 16);
set(gca, 'FontSize', 16);

figure; histogram(yendminus1(yendminus1<=0.005));

yendminus2 = arrayfun(@(x) x.rawTrack.rawState(end-2,5),landing_tracks);
figure;
histogram(yendminus2);
xlabel('y at the end but two point', 'FontSize', 16);
ylabel('Occurences', 'FontSize', 16);
set(gca, 'FontSize', 16);

figure; histogram(yendminus2(yendminus2<=0.005));

y = abs([yend'-yendminus1' yendminus1'-yendminus2'])*1000;
y(y(:,1)>1)'
%% Plot histogram of rref for static, expanding and contracting spirals
landing_tracks = [landingTracks{:}];
state_LDF = [landing_tracks.state_LDF];
rref_segments = [state_LDF.rrefSegments];

chosen_fac = 2.5;
rref_segments1 = rref_segments(abs([rref_segments.factor] - chosen_fac)<1e-6 & ...
                              (abs([rref_segments.patternnum] - 7)<1e-6 | ...
                               abs([rref_segments.patternnum] - 8)<1e-6 | ...
                               abs([rref_segments.patternnum] - 9)<1e-6));
rref_static = vertcat(rref_segments1.rmean_ti);

rref_segments2 = rref_segments(abs([rref_segments.factor] - chosen_fac)<1e-6 & ...
                              (abs([rref_segments.patternnum] - 12)<1e-6 | ...
                               abs([rref_segments.patternnum] - 13)<1e-6 | ...
                               abs([rref_segments.patternnum] - 14)<1e-6));
rref_expanding = vertcat(rref_segments2.rmean_ti);

rref_segments3 = rref_segments(abs([rref_segments.factor] - chosen_fac)<1e-6 & ...
                              (abs([rref_segments.patternnum] - 15)<1e-6 | ...
                               abs([rref_segments.patternnum] - 16)<1e-6 | ...
                               abs([rref_segments.patternnum] - 17)<1e-6));
rref_contracting = vertcat(rref_segments3.rmean_ti);

[mean(-rref_static) mean(-rref_expanding) mean(-rref_contracting)]
[median(-rref_static) median(-rref_expanding) median(-rref_contracting)]
% % Histogram of r*
% figure; hold on;
% histogram(-rref_static);
% histogram(-rref_expanding);
% histogram(-rref_contracting);
% xlabel('Reference rate of expansion, r* (s-1)', 'FontSize', 16);
% ylabel('Occurences', 'FontSize', 16);
% set(gca, 'FontSize', 16);

map = brewermap(3,'Set1'); 
figure;
histogram(-rref_static,'facecolor',map(1,:),'facealpha',.5,'edgecolor','none')
hold on
histogram(-rref_expanding,'facecolor',map(2,:),'facealpha',.5,'edgecolor','none')
histogram(-rref_contracting,'facecolor',map(3,:),'facealpha',.5,'edgecolor','none')
set(gca, 'FontSize', 16);
box off
axis tight
legend('Static','Expanding','Contracting','location','northeast')
legend boxoff

%% Extract dataset only for static patterns
clc; close all;
clear;

inputFile = '/media/reken001/Disk_12/honeybee_experiments/postprocessing/BlindLandingtracks_A4_LDF_rref.mat';
load(inputFile);

landing_tracks = [landingTracks{:}];
static_patternnums = [1 2 3 4 5 6 10 11 18];
landingTracks = landing_tracks(arrayfun(@(x) ismember(x,static_patternnums), [landing_tracks.patternnum]));
% landingTracks = landing_tracks([landing_tracks.patternnum] <= 11);

% Loading data and collecting segments with rref
factors = [0.25:0.25:2.5];
chosen_fac = 1.5;

data = struct.empty;
clear dummy;

state_LDF = [landingTracks.state_LDF];
rrefSegs = [state_LDF.rrefSegments];

landingTracks_indx4stateLDF = arrayfun(@(x) x*ones(length(landingTracks(x).state_LDF),1),1:length(landingTracks), 'UniformOutput', false);
landingTracks_indx4stateLDF = vertcat(landingTracks_indx4stateLDF{:});


% Finding # of "tracks" (state LDFs) that contain rref segments
% for each factor
N = zeros(1, length(factors));
for ct_fac = 1:length(factors)
    factor = factors(ct_fac);
    indices = arrayfun(@(x) ~isempty(x.rrefSegments(abs([x.rrefSegments.factor]-factor)<1e-6).intervals_ti),state_LDF);
    N(ct_fac) = sum(indices);
end

% Display
disp(['# of state LDFs containing rrefs for different factors: ' num2str(N)]);

% For chosen factor, collect all tracks that do contain
% non-empty rref segments
indices = arrayfun(@(x) ~isempty(x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).intervals_ti), state_LDF);
dummy.tracks_fac = state_LDF(indices);
dummy.landingTrack = landingTracks(landingTracks_indx4stateLDF(indices));
data = [data; dummy];

%% Create Trajectories' views with constant-r segments highlighted
close all;
trajPlot = figure; hold on;
skip_step = 1;
N = 0;
for ct=1:length(data)
    for ct1 = 1:skip_step:length(data(ct).tracks_fac)
        xyz = data(ct).tracks_fac(ct1).filteredState(:,[2 3 4]); % the complete trajectory

%         if ~any(xyz(:,2)>-5e-3)
    %         plot3(xyz(:,1), xyz(:,2), xyz(:,3),'Color',[0 97 205]./255); 
%             plot3(xyz(:,1), -1*xyz(:,2), xyz(:,3),'Color',[117 153 242]./255); 
            plot3(xyz(:,1), -1*xyz(:,2), xyz(:,3),'Color',[252,187,161]./255); 
    %         plot3(xyz(:,1), xyz(:,2), xyz(:,3));
%         else
%             disp(true)
%         end
    end
    N = N + length(1:skip_step:length(data(ct).tracks_fac));
end

for ct=1:length(data)
    for ct1 = 1:skip_step:length(data(ct).tracks_fac)
        intervals = data(ct).tracks_fac(ct1).rrefSegments(...
            abs([data(ct).tracks_fac(ct1).rrefSegments.factor]-chosen_fac)<1e-6).intervals_ti;
        for ct2 = 1:size(intervals,1)
            xyz4rref = data(ct).tracks_fac(ct1).filteredState(intervals(ct2,1):intervals(ct2,2),[2 3 4]); % trajectory used for estimating r*
            plot3(xyz4rref(:,1), -1*xyz4rref(:,2), xyz4rref(:,3),'Color',[215 48 39]./255, 'Linewidth', 1);
%             plot3(xyz4rref(:,1), -1*xyz4rref(:,2), xyz4rref(:,3),'Color',[239 38 38]./255, 'Linewidth', 1);
        end
    end
end

% % Landing Disc
% [X,Y,Z] = cylinder(treatment.landingDiscs(1).radius);
radius = 0.30; % in m

figure(trajPlot);
% h=mesh(X,Y,Z,'facecolor',[1 0 0]); % draw landing disc
plot3([-radius; radius], [0 0], [0 0], 'LineWidth', 2, 'Color', [83 83 83]./255);
plot3([0 0], [0 0], [-radius; radius], 'LineWidth', 2, 'Color', [83 83 83]./255);
zlabel('z (m)', 'FontSize', 14);
ylabel('y (m)', 'FontSize', 14);
xlabel('x (m)', 'FontSize', 14);
axis equal;
xlim([-0.3 0.3]);
xticks([-0.3:0.1:0.3]);
ylim([0 0.6]);
yticks([0:0.1:0.6]);
zlim([-0.3 0.3]);
zticks([-0.3:0.1:0.3]);
set(gca, 'FontSize', 16);
view(0,90);
view(-90,0);

%% Histogram of rref
close all;
for ct=1:length(data)
    dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).rref_ti,data(ct).tracks_fac,'UniformOutput',false);
    data(ct).rref = vertcat(dummy{:});
    
    dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).rmean_ti,data(ct).tracks_fac,'UniformOutput',false);
    data(ct).rmean = vertcat(dummy{:});

    dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).ymean_ti,data(ct).tracks_fac,'UniformOutput',false);
    data(ct).ymean = vertcat(dummy{:});
    
    dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).fd_analytical_ti,data(ct).tracks_fac,'UniformOutput',false);
    data(ct).fd_analytical_ti = vertcat(dummy{:});
    
    dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).fd_actual_ti,data(ct).tracks_fac,'UniformOutput',false);
    data(ct).fd_actual_ti = vertcat(dummy{:});
    
end
% figure;
% histogram(-vertcat(data.rref), [0:0.5:8]);
% xlabel('Estimated set-points, r* (1/s)', 'FontSize', 16);
% ylabel('Occurances', 'FontSize', 16);
% set(gca, 'FontSize', 16);

dummy = abs(vertcat(data.rref)-vertcat(data.rmean));
[mean(dummy) median(dummy) max(dummy)]

dummy = abs(vertcat(data.fd_analytical_ti)-vertcat(data.fd_actual_ti));
[mean(dummy) median(dummy) max(dummy)]

figure;
% histogram(-vertcat(data.rmean), [0:0.5:8]);
% histogram(-vertcat(data.rmean), [0:0.5:9.5]);
histfit(-vertcat(data.rmean),13,'gamma')
xlabel('Estimated set-points, r* (1/s)', 'FontSize', 16);
ylabel('Occurences', 'FontSize', 16);
set(gca, 'FontSize', 16);
dummy = -vertcat(data.rmean);
[mean(dummy) median(dummy) max(dummy) min(dummy)]
pd = fitdist(dummy,'Gamma');

%% Panel d for Figure 4 (ymean vs rref)
% plotted only for tracks containing >1 r* segments

close all; clc;
saveDir = '/home/reken001/Pulkit/graphs_temp';
cmap = [60 180 75; 245 130 48; 0 130 200]/255; % [green orange blue]


delta_rrref_posjumps = cell(length(data),1); % change in rref between two consecutive rref segments
delta_rrref_negjumps = cell(length(data),1); % change in rref between two consecutive rref segments

ratio_rrref_posjumps = cell(length(data),1); % change in rref between two consecutive rref segments
ratio_rrref_negjumps = cell(length(data),1); % change in rref between two consecutive rref segments

data_ymean = {}; data_rmean = {};
for ct=1:length(data)
    has_multiple_rrefs = arrayfun(@(x) length(x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).rref_ti)>1,data(ct).tracks_fac);
    
    N = sum(has_multiple_rrefs);
    data_ss = data(ct).tracks_fac(has_multiple_rrefs);
    
    rmean = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).rmean_ti,data_ss,'UniformOutput',false);

    ymean = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).ymean_ti,data_ss,'UniformOutput',false);
    
    data_ymean{end+1} = ymean;
    data_rmean{end+1} = rmean;
    
    landingTracks = data(ct).landingTrack(has_multiple_rrefs);
%     track_indx = arrayfun(@(x) x*ones(length(data_ss(x).rrefSegments(abs([data_ss(x).rrefSegments.factor]-chosen_fac)<1e-6).ymean_ti),1),1:N,'UniformOutput',false);
    
    jumps = cell(size(ymean));
    jump_ratios = cell(size(ymean));
    disp(['# of tracks with >1 r* segments for chosen factor: ' num2str(N)]);
    plotHandle = figure; hold on;
    for ct1=1:N
        
%         figure(plotHandle);
        % Exclude tracks wherein bbee moves back and forth between landing
        % platforms
%         if any(abs(diff(data_ss(ct1).filteredState(data_ss(ct1).filteredState(:,6) > 0,3))) > 0.05)
%             continue;
%         end
        
        y_r = sortrows([-[ymean{ct1}] -[rmean{ct1}]],1,'descend');
        
        % Finding whether the next r* is higher or lower
        jumps{ct1} = diff(y_r);
        jump_ratios{ct1} = y_r(2:end,2)./y_r(1:end-1,2);
        
        % Plot straight line
%         plot(-[ymean{ct1}], -[rmean{ct1}] ...
%         ,'Color', [225 225 225]./255/2,'LineWidth',1);
        if rem(ct1,2) == 1
            plot(-[ymean{ct1}], -[rmean{ct1}] ...
            ,'Color', [225 225 225]./255/2,'LineWidth',1);
        end
    
        
        % Highlight data points
% % %         scatter(-[ymean{ct1}], -[rmean{ct1}],5,...
% % %         'k','filled','o');
% %         
        color_indx = 1:1:length(ymean{ct1});
        color_indx(color_indx > 3) = 3;
        scatter(-[ymean{ct1}], -[rmean{ct1}],10,...
        cmap(color_indx,:),'filled','o');
    
% % %         % Plot only negative jumps
% % %         if any(jumps{ct1}(:,2)<-2)
% % %             plot(-[ymean{ct1}], -[rmean{ct1}] ...
% % %             ,'Color', [225 225 225]./255/2,'LineWidth',1);
% % %         end
% % 

        
    end
    disp(['# r*segments in tracks with >1 r* segments for chosen fac: ' num2str(length(vertcat(ymean{:})))]);
    
    dummy = vertcat(jumps{:});
    delta_rrref_posjumps{ct} = dummy(dummy(:,2)>0,2);
    delta_rrref_negjumps{ct} = dummy(dummy(:,2)<0,2);
    
    dummy1 = vertcat(jump_ratios{:});
    assert(sum((dummy1>1) == (dummy(:,2)>0)) == length(dummy1))
    ratio_rrref_posjumps{ct} = dummy1(dummy1>1);
    ratio_rrref_negjumps{ct} = dummy1(dummy1<1);
    
    disp(['# of +ve jumps in r* for chosen fac: ' num2str(sum(dummy(:,2)>0))]);
    disp(['# of -ve jumps in r* for chosen fac: ' num2str(sum(dummy(:,2)<0))]);
    
    disp(['% of # of +ve jumps in r* for chosen fac: ' num2str(sum(dummy(:,2)>0)/(sum(dummy(:,2)<0)+sum(dummy(:,2)>0)))]);
    disp(['% of # of -ve jumps in r* for chosen fac: ' num2str(sum(dummy(:,2)<0)/(sum(dummy(:,2)<0)+sum(dummy(:,2)>0)))]);
    
    disp(['Mean +ve jumps in r* for chosen fac: ' num2str(mean(dummy(dummy(:,2)>0, 2)))]);
    disp(['Mean -ve jumps in r* for chosen fac: ' num2str(mean(dummy(dummy(:,2)<0, 2)))]);
    
    plotHandle.Position(3) = 680; plotHandle.Position(4) = 545;
    
    figure(plotHandle);
    
    set(gca, 'FontSize', 16);
    ylabel('Estimated set-points, r* (1/s)', 'FontSize', 16);
    xlabel('y* (m)', 'FontSize', 16);
    
%     print(plotHandle, fullfile(saveDir,['multiple_rrefs_' pattern{data(ct).ct_pattern} '_' light{data(ct).ct_light}]), '-dpdf');

end

plotHandle1 = figure; hold on;
% cmap = [1 0 0; 0 1 0; 0 0 1];
cmap = [60 180 75; 245 130 48; 0 130 200]/255; % [green orange blue]
for ct=1:length(data_ymean)
    for ct1=1:length(data_ymean{ct})
% %         if rem(ct1,15) == 1
% %             plot(-[ymean{ct1}], -[rmean{ct1}] ...
% %                 ,'Color', [225 225 225]./255/2,'LineWidth',1);
% %         end

        color_indx = 1:1:length(data_ymean{ct}{ct1});
        color_indx(color_indx > 3) = 3;
        scatter(-[data_ymean{ct}{ct1}], -[data_rmean{ct}{ct1}],10,...
                cmap(color_indx,:),'filled','o');
    end
end
set(gca, 'FontSize', 16);
ylabel('Estimated set-points, r* (s-1)', 'FontSize', 16);
xlabel('y* (m)', 'FontSize', 16);
% Plot fit averaged over all light conditions
ymean = arrayfun(@(x) vertcat(data_ymean{x}{:}), 1:length(data_ymean), 'UniformOutput', false);
ymean = vertcat(ymean{:});
y_vec = min(-ymean):0.001:max(-ymean);
modelfun = @(b,x)(exp(b(1))*x.^(b(2)));
Coefficients = [-0.8945715 -0.7468271]; % from R
plot(y_vec,modelfun(Coefficients,y_vec),'Color', [0 0 0], 'LineWidth', 2);


[mean(vertcat(delta_rrref_posjumps{:})) std(vertcat(delta_rrref_posjumps{:}))]
length(vertcat(delta_rrref_posjumps{:}))
[mean(vertcat(delta_rrref_negjumps{:})) std(vertcat(delta_rrref_negjumps{:}))]
length(vertcat(delta_rrref_negjumps{:}))

delta_rref_per_treatment_mean = arrayfun(@(x) [mean(vertcat(delta_rrref_posjumps{x})) mean(vertcat(delta_rrref_negjumps{x}))], 1:length(delta_rrref_posjumps), 'UniformOutput', false);
delta_rref_per_treatment_mean = vertcat(delta_rref_per_treatment_mean{:})

figure;
histogram([vertcat(delta_rrref_posjumps{:}); vertcat(delta_rrref_negjumps{:})])
% histfit([vertcat(delta_rrref_posjumps{:}); vertcat(delta_rrref_negjumps{:})], [], 'Gamma')
xlabel('\Delta r*', 'FontSize', 16);
ylabel('Occurences', 'FontSize', 16);
set(gca, 'FontSize', 16);

[mean(vertcat(ratio_rrref_posjumps{:})) std(vertcat(ratio_rrref_posjumps{:}))]
[mean(vertcat(ratio_rrref_negjumps{:})) std(vertcat(ratio_rrref_negjumps{:}))]

dummy = cellfun(@(x) [mean(x) std(x)],delta_rrref_posjumps,'UniformOutput',false);
dummy = vertcat(dummy{:})
for ct=1:1
    disp([num2str(dummy(ct,1),'%0.3f') ' [' num2str(dummy(ct,2),'%0.3f') ']']);
end
dummy =cellfun(@(x) [mean(x) std(x)],delta_rrref_negjumps,'UniformOutput',false);
vertcat(dummy{:})
dummy =cellfun(@(x) [mean(x) std(x)],ratio_rrref_posjumps,'UniformOutput',false);
vertcat(dummy{:})
dummy =cellfun(@(x) [mean(x) std(x)],ratio_rrref_negjumps,'UniformOutput',false);
vertcat(dummy{:})

figure(plotHandle);
has_one_rref = arrayfun(@(x) length(x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).rref_ti)==1,data(ct).tracks_fac);

N = sum(has_one_rref);
data_ss = data(ct).tracks_fac(has_one_rref);

rmean = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).rmean_ti,data_ss,'UniformOutput',false);
ymean = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).ymean_ti,data_ss,'UniformOutput',false);
    
scatter(-[ymean{:}], -[rmean{:}],10,...
                [207 52 118]./255,'filled','o');

%% Panel d for Figure 4 (ymean vs rref)
% plotted for all tracks

close all; clc;
saveDir = '/home/reken001/Pulkit/graphs_temp';
cmap = [60 180 75; 245 130 48; 0 130 200]/255; % [green orange blue]

data_ymean = {}; data_rmean = {};
for ct=1:length(data)
    
    data_ss = data(ct).tracks_fac;
    
    rmean = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).rmean_ti,data_ss,'UniformOutput',false);

    ymean = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).ymean_ti,data_ss,'UniformOutput',false);
    
    data_ymean{end+1} = ymean;
    data_rmean{end+1} = rmean;
end

plotHandle1 = figure; hold on;
% cmap = [1 0 0; 0 1 0; 0 0 1];
cmap = [60 180 75; 245 130 48; 0 130 200]/255; % [green orange blue]
for ct=1:length(data_ymean)
    for ct1=1:length(data_ymean{ct})
% %         if rem(ct1,15) == 1
% %             plot(-[ymean{ct1}], -[rmean{ct1}] ...
% %                 ,'Color', [225 225 225]./255/2,'LineWidth',1);
% %         end

        color_indx = 1:1:length(data_ymean{ct}{ct1});
        color_indx(color_indx > 3) = 3;
        scatter(-[data_ymean{ct}{ct1}], -[data_rmean{ct}{ct1}],10,...
                cmap(color_indx,:),'filled','o');
    end
end
set(gca, 'FontSize', 16);
ylabel('Estimated set-points, r* (s-1)', 'FontSize', 16);
xlabel('y* (m)', 'FontSize', 16);
% Plot fit averaged over all light conditions
ymean = arrayfun(@(x) vertcat(data_ymean{x}{:}), 1:length(data_ymean), 'UniformOutput', false);
ymean = vertcat(ymean{:});
y_vec = 0.05:0.001:max(-ymean);
modelfun = @(b,x)(exp(b(1))*x.^(b(2)));
Coefficients = [0.78618 -0.21784]; % from R
plot(y_vec,modelfun(Coefficients,y_vec),'Color', [0 0 0], 'LineWidth', 2);


%% Plot ymean vs rmean
% Panel a
close all;

taudots = -0.21784;

intercepts = 0.78618; 
          
points_cmap = [252,187,161;
200,200,200]./255;

line_cmap = [215,48,31;
37,37,37]./255;


dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).rmean_ti,data.tracks_fac,'UniformOutput',false);
rmean= [vertcat(dummy{:})];

dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).ymean_ti,data.tracks_fac,'UniformOutput',false);
ymean = [vertcat(dummy{:})];


data2plot = [ymean rmean 2*ones(size(rmean))];
data2plot = data2plot(randperm(size(data2plot, 1)), :);
figure; hold on;
% plot lines
y_vec = min(-ymean):0.001:max(-ymean);
modelfun = @(b,x)(exp(b(1))*x.^(b(2)));

Coefficients = [mean(intercepts) mean(taudots)];
plot(y_vec,modelfun(Coefficients,y_vec),'Color', line_cmap(2,:),'LineWidth',3);

set(gca, 'FontSize', 16);
ylabel('Estimated set-points, r* (1/s)', 'FontSize', 16);
xlabel('y* (m)', 'FontSize', 16);

scatter(-data2plot(:,1),-data2plot(:,2),10,points_cmap(round(data2plot(:,3)),:),'filled','o');
% ylim([0 10])


% Panel b - log plots
modelfun = @(b,x)(b(1)+b(2).*x);
figure; hold on;
Coefficients = [mean(intercepts) mean(taudots)];
plot(log(y_vec),modelfun(Coefficients,log(y_vec)),'Color', line_cmap(2,:),'LineWidth',3);
set(gca, 'FontSize', 16);
ylabel('log(r*) (1/s)', 'FontSize', 16);
xlabel('log(y*) (m)', 'FontSize', 16);

scatter(log(-data2plot(:,1)),log(-data2plot(:,2)),10,points_cmap(round(data2plot(:,3)),:),'filled','o');



%% Plot ymean vs rmean with yrange covered with rref highlighted
close all;

taudots = -0.21784;
intercepts = 0.78618;
          
points_cmap = [252,187,161;
200,200,200]./255;

line_cmap = [215,48,31;
37,37,37]./255;


dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).rmean_ti,data.tracks_fac,'UniformOutput',false);
rmean= [vertcat(dummy{:})];

dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).ymid,data.tracks_fac,'UniformOutput',false);
ymid = [vertcat(dummy{:})];

dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).ymean_ti,data.tracks_fac,'UniformOutput',false);
ymean = [vertcat(dummy{:})];

dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).yrange,data.tracks_fac,'UniformOutput',false);
yrange = [vertcat(dummy{:})];

dummy = arrayfun(@(x) x*ones(length(data.tracks_fac(x).rrefSegments(abs([data.tracks_fac(x).rrefSegments.factor]-chosen_fac)<1e-6).yrange),1),1:length(data.tracks_fac),'UniformOutput',false);
indices = [vertcat(dummy{:})];
tracks_fac = data.tracks_fac(indices)

% data2plot = [ymean rmean 2*ones(size(rmean))];
data2plot = [ymid rmean 2*ones(size(rmean))];
% data2plot = data2plot(randperm(size(data2plot, 1)), :);
figure; hold on;
% plot lines
% y_vec = min(-ymean):0.001:max(-ymean);
% modelfun = @(b,x)(exp(b(1))*x.^(b(2)));
% 
% Coefficients = [mean(intercepts) mean(taudots)];
% plot(y_vec,modelfun(Coefficients,y_vec),'Color', line_cmap(2,:),'LineWidth',3);

set(gca, 'FontSize', 16);
ylabel('Estimated set-points, r* (1/s)', 'FontSize', 16);
xlabel('y* (m)', 'FontSize', 16);

plot_every = 2;
% scatter(-data2plot(:,1),-data2plot(:,2),10,points_cmap(round(data2plot(:,3)),:),'filled','o');
errorbar(-data2plot(1:plot_every:end,1),-data2plot(1:plot_every:end,2),yrange(1:plot_every:end)/2,'horizontal','LineStyle','none','Marker','o','MarkerSize',4 ...
    ,'MarkerFaceColor',points_cmap(2,:));
% hold on;
% scatter(-ymean,-data2plot(:,2),10,points_cmap(round(data2plot(:,3)),:),'filled','o');

% ylim([0 10])


% sum((-ymean-yrange/2)<0)
% tracks_fac_ss = tracks_fac((-ymean-yrange/2)<0)
% for ct=1:3
%     tracks_fac_ss(ct).plot_rrefs(chosen_fac);
% end

%% Plot ymean vs rmean with groups from average analysis
close all;

groups = {[6],[4,5]};
cmap = lines(length(groups));

          
points_cmap = [252,187,161;
200,200,200]./255;

line_cmap = [215,48,31;
37,37,37]./255;



dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).rmean_ti,data.tracks_fac,'UniformOutput',false);
rmean= [vertcat(dummy{:})];

dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).ymid,data.tracks_fac,'UniformOutput',false);
ymid = [vertcat(dummy{:})];

dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).ymean_ti,data.tracks_fac,'UniformOutput',false);
ymean = [vertcat(dummy{:})];

dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).yrange,data.tracks_fac,'UniformOutput',false);
yrange = [vertcat(dummy{:})];

dummy = arrayfun(@(x) x*ones(length(data.tracks_fac(x).rrefSegments(abs([data.tracks_fac(x).rrefSegments.factor]-chosen_fac)<1e-6).yrange),1),1:length(data.tracks_fac),'UniformOutput',false);
indices = [vertcat(dummy{:})];
tracks_fac = data.tracks_fac(indices)

dummy = arrayfun(@(x) data.landingTrack(x).patternnum*ones(length(data.tracks_fac(x).rrefSegments(abs([data.tracks_fac(x).rrefSegments.factor]-chosen_fac)<1e-6).yrange),1),1:length(data.tracks_fac),'UniformOutput',false);
patternnums = [vertcat(dummy{:})];

plotHandle = figure; hold on;
legend_str = {};
for ct=1:length(groups)
    for ct1=1:length(groups{ct})
        scatter(-ymean(patternnums==groups{ct}(ct1)),-rmean(patternnums==groups{ct}(ct1)),10,cmap(ct,:),'filled','o');
    end
    legend_str{end+1} = num2str(groups{ct});
end

set(gca, 'FontSize', 16);
ylabel('Estimated set-points, r* (1/s)', 'FontSize', 16);
xlabel('y* (m)', 'FontSize', 16);
legend(legend_str);


% data2plot = [ymean rmean 2*ones(size(rmean))];
data2plot = [ymid rmean 2*ones(size(rmean))];
% data2plot = data2plot(randperm(size(data2plot, 1)), :);
figure; hold on;
% plot lines
% y_vec = min(-ymean):0.001:max(-ymean);
% modelfun = @(b,x)(exp(b(1))*x.^(b(2)));
% 
% Coefficients = [mean(intercepts) mean(taudots)];
% plot(y_vec,modelfun(Coefficients,y_vec),'Color', line_cmap(2,:),'LineWidth',3);

set(gca, 'FontSize', 16);
ylabel('Estimated set-points, r* (1/s)', 'FontSize', 16);
xlabel('y* (m)', 'FontSize', 16);

plot_every = 2;
% scatter(-data2plot(:,1),-data2plot(:,2),10,points_cmap(round(data2plot(:,3)),:),'filled','o');
errorbar(-data2plot(1:plot_every:end,1),-data2plot(1:plot_every:end,2),yrange(1:plot_every:end)/2,'horizontal','LineStyle','none','Marker','o','MarkerSize',4 ...
    ,'MarkerFaceColor',points_cmap(2,:));
% hold on;
% scatter(-ymean,-data2plot(:,2),10,points_cmap(round(data2plot(:,3)),:),'filled','o');

%% Comparison of bbees and hbees
close all;
fig1 = figure; hold on;
x = 0:0.1:8;
bbee_rref_hist = makedist('Gamma','a',3.59,'b',0.65);
hbee_rref_hist = fitdist(-vertcat(data.rmean),'Gamma');
line_cmap = [215,48,31;
        37,37,37]./255;
% plot(h(1).XData, pdf(bbee_rref_hist, h(1).XData)*sum(h(1).YData*(diff(h(1).XData(1:2)))),'g','Linewidth',2);
plot(x, gampdf(x, bbee_rref_hist.a, bbee_rref_hist.b),'Color', line_cmap(1,:) ,'Linewidth',2);
plot(x, gampdf(x, hbee_rref_hist.a, hbee_rref_hist.b),'Color', line_cmap(2,:),'Linewidth',2);
xlabel('Estimated set-points, r* (1/s)', 'FontSize', 16);
ylabel('Occurences', 'FontSize', 16);
legend({'bumblebees','honeybees'})
set(gca, 'FontSize', 16);

intercepts_bbee = [-1.456806   -1.379435   -1.287747        -1.270285   -1.192914   -1.101225  -1.000403   -0.923032   -0.8313437        -0.8138813  -0.7365104  -0.6448221];
fig2 = figure; hold on;
y_vec_bbee = 0.05:0.001:0.35;
y_vec_hbee = 0.03:0.001:0.56;
modelfun = @(b,x)(exp(b(1))*x.^(b(2)));
Coefficients_bbee = [mean(intercepts_bbee) -0.872];
Coefficients_hbee = [intercepts taudots];
plot(y_vec_bbee,modelfun(Coefficients_bbee,y_vec_bbee),'Color', line_cmap(1,:),'LineWidth',3);
plot(y_vec_hbee,modelfun(Coefficients_hbee,y_vec_hbee),'Color', line_cmap(2,:),'LineWidth',3);
legend({'bumblebees','honeybees'})
set(gca, 'FontSize', 16);
ylabel('Estimated set-points, r* (1/s)', 'FontSize', 16);
xlabel('y* (m)', 'FontSize', 16);

%% Plot r* vs y*

% landing_tracks = [landingTracks{:}];
state_LDF = [landing_tracks.state_LDF];
rref_segments = [state_LDF.rrefSegments];

chosen_fac = 1.0;
rref_segments = rref_segments(abs([rref_segments.factor] - chosen_fac)<1e-6);
ymean = vertcat(rref_segments.ymean_ti); rref = vertcat(rref_segments.rmean_ti);

% Histogram of r*
figure;
histogram(-rref);
xlabel('Reference rate of expansion, r* (s-1)', 'FontSize', 16);
ylabel('Occurences', 'FontSize', 16);
set(gca, 'FontSize', 16);

% Plot mean rref vs ymean 
d1=fitlm(log(-ymean),log(-rref))
% d1.Rsquared.Ordinary
modelfun = @(b,x)(exp(b(1))*x.^(b(2)));
figure; hold on;
f = scatter(-ymean,-rref,40,[252,187,161]./255,'filled','s','MarkerEdgeColor','k');
set(gca, 'FontSize', 18);
ylabel('Estimated set-point of optical expansion rate, r* (s-1)', 'FontSize', 16);
xlabel('Distance from the platform, y* (m)', 'FontSize', 16);
plot(min(-ymean):0.001:max(-ymean),modelfun(d1.Coefficients.Estimate,min(-ymean):0.001:max(-ymean)),'Color', [69 117 180]./255,'LineWidth',3);

% % Plot log(rref) vs log(ymean)
modelfun = @(b,x)(b(1)+b(2).*x);
x_vec = min(log(-ymean)):0.001:max(log(-ymean));
figure; hold on;
f = scatter(log(-ymean),log(-rref),40,[252,187,161]./255,'filled','s','MarkerEdgeColor','k');
plot(x_vec,modelfun(d1.Coefficients.Estimate,x_vec),'Color', [69 117 180]./255,'LineWidth',3);
set(gca, 'FontSize', 18);
ylabel('log(r*) (s-1)', 'FontSize', 16);
xlabel('log(y*) (m)', 'FontSize', 16);


% for ct_treatment = 1:length(dataFiles)
%     if ~isempty(landingTracks{ct_treatment})
%         for ct=1:length(landingTracks{ct_treatment})
% 
%         end
%     end
% end


%% create 3d trajectory figure
% dataFiles = vertcat(dataFiles{2:end});
trajFig = figure; hold on;
landing_tracks = [landingTracks{:}];
for ct=1:length(landing_tracks)
    plot3(landing_tracks(ct).state_LDF.filteredState(:,2), landing_tracks(ct).state_LDF.filteredState(:,3),...
        landing_tracks(ct).state_LDF.filteredState(:,4));
end
zlabel('z (m)', 'FontSize', 14);
ylabel('y (m)', 'FontSize', 14);
xlabel('x (m)', 'FontSize', 14);
axis equal;
set(gca, 'FontSize', 16);
view([0 90]);
view([90 0]);