%% Extract track excerpts using custom GUI
% This is done to estimate parameters for rref automatic extraction
%%
% Note that some custom defined class definitions come from "legacy" code
% (code written before analysis framework of Article 1 was finalized)
clc; clear; close all;

% to generate colormaps
addpath('./lib/DrosteEffect-BrewerMap-221b913');

% to add definition of classes being used for this analysis
addpath('./lib/li_analysis');

% to include class definitions used in videoscorer_data
addpath('./lib/video_scorer/');
% rmpath('./legacy_lib/video_scorer/Source');

% to include higher order accurate differentiation function
addpath('./lib/diffxy');

% to include fmf reader
addpath('./lib/flymovieformat');

% to include hline and vline function
addpath('./lib/hline_vline');

% % to include all custom defined class definitions
% addpath('/home/reken001/Pulkit/MATLAB');

%% Create treatments and landing tracks

dataDir = '/media/reken001/Disk_12/honeybee_experiments/landing_tracks';

% Look for dataFiles in dataDir and its subfolders
subdir = dir(dataDir);
subdir = subdir([subdir.isdir]);
subdir = subdir(~ismember({subdir.name},{'..'}));

treatments = Lighttreatment.empty;
for ct_treatment=1:length(subdir)
    filesmat = dir(fullfile(subdir(ct_treatment).folder,subdir(ct_treatment).name,'*.mat'));
    filescsv = dir(fullfile(subdir(ct_treatment).folder,subdir(ct_treatment).name,'*.csv'));
    if ~isempty(filesmat) || ~isempty(filescsv)
%         [~,pattern,~] = fileparts(subdir(ct_treatment).folder);
        pattern = subdir(ct_treatment).name;
        treatments(end+1) = Lighttreatment('Nan', pattern, 'Nan', nan, nan);
        dataFiles{ct_treatment} = [filesmat; filescsv];
        
        % Initializing the landing tracks
        nObjects = length(dataFiles{ct_treatment});
        landingTracks = BlindLandingtrack.empty(nObjects, 0);
        
        % Fill landing tracks
        for ct=1:nObjects
            landingTracks(ct) = BlindLandingtrack(pattern, 'Nan', pattern, dataFiles{ct_treatment}(ct).name, ct);
            landingTracks_per_object = landingTracks(ct);
            
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
            
%             if ct==18
%                 keyboard;
%             end

            % Make correction for y-offset (by setting 0,0,0 as the origin
            % of the reference frame (for trajectories with yend>0.005m)
            if -dataFiles{ct_treatment}(ct).data.y(end) > 0.005*1000
                dataFiles{ct_treatment}(ct).data.x = -dataFiles{ct_treatment}(ct).data.x(end,1) + dataFiles{ct_treatment}(ct).data.x;
                dataFiles{ct_treatment}(ct).data.y = -dataFiles{ct_treatment}(ct).data.y(end,1) + dataFiles{ct_treatment}(ct).data.y;
                dataFiles{ct_treatment}(ct).data.z = -dataFiles{ct_treatment}(ct).data.z(end,1) + dataFiles{ct_treatment}(ct).data.z;
            end
            
            N = length(dataFiles{ct_treatment}(ct).data.x);
            t = [0:N-1]'/400;
            rawState = [ct*ones(N,1) [1:N]' t dataFiles{ct_treatment}(ct).data.x/1000 -dataFiles{ct_treatment}(ct).data.y/1000 dataFiles{ct_treatment}(ct).data.z/1000 ...
            diffxy(t, dataFiles{ct_treatment}(ct).data.x/1000) diffxy(t,-dataFiles{ct_treatment}(ct).data.y/1000) diffxy(t, dataFiles{ct_treatment}(ct).data.z/1000)...
            zeros(N, 9)];% 18 columns (obj_id	frame	timestamp	x	y	z	xvel	yvel	zvel P00	P01	P02	P11	P12	P22	P33	P44	P55)
            
            landingTracks_per_object.rawTrack(end+1) = rawState_BlindLandingtrack(rawState, 'Hive');
            landingTracks_per_object.dt = 1/400;
            landingTracks_per_object.filterRawTracks(20);   
            
            landingTracks_per_object.raw_LDF = landingTracks_per_object.rawTrack;
            landingTracks_per_object.state_LDF = landingTracks_per_object.state;
        end
        treatments(end).landingTracks = landingTracks(arrayfun(@(x) ~isempty(x.rawTrack), landingTracks));
    end
    
end


relevantTreatments = treatments;

%% Copy previous data extracted using GUI (when there was y-offset in the trajectories) (not always needed, check if you need it!)
% prevData = load('/media/reken001/Disk_12/honeybee_experiments/postprocessing/with_yoffset/BlindLandingtracks.mat');
prevData = load('/media/reken001/Disk_12/honeybee_experiments/postprocessing/BlindLandingtracks_manualrref.mat');
prevTreatments = prevData.treatments;

assert(length(prevTreatments)==length(treatments));
for ct_treatment=1:length(prevTreatments)
    prevTreatment = prevTreatments(ct_treatment);
    treatment = treatments(ct_treatment);    
    assert(length(prevTreatment.landingTracks)==length(treatment.landingTracks));
    for ct_track=1:length(treatment.landingTracks)
        treatment.landingTracks(ct_track).DataGUI = prevTreatment.landingTracks(ct_track).DataGUI;
    end
end

%% Extract already present GUI data 
outputDataDir = '/media/reken001/Disk_12/honeybee_experiments/postprocessing';
appData_filename = ['GUIDE_appData_Blindtracks_UNKNOWN.mat'];
if exist(fullfile(outputDataDir, appData_filename), 'file')
    load(fullfile(outputDataDir, appData_filename));
end
tracks_filename = ['BlindLandingtracks.mat'];
if exist(fullfile(outputDataDir, tracks_filename), 'file')
    load(fullfile(outputDataDir, tracks_filename));
    relevantTreatments = treatments;
end

%% Go below after selecting rref intervals using GUI

%% Extract data saved using GUI
% Inputs
close all;
% 
inputFile = '/media/reken001/Disk_12/honeybee_experiments/postprocessing/BlindLandingtracks_manualrref.mat';
load(inputFile);

% Plot code is not working yet. 
DirPlots = '/media/reken001/Disk_12/honeybee_experiments/plots/manual_rref_estimate';
delPreviousPlots = false;
savePlots = false;

rref_meanVbyy = [];
rref_fitVvsy = [];
Rsquared = [];
ymean = [];
vmean = [];
all_rvsy = [];
clear data;


ylimit = 0.15; % delta y for which ds/y will be calculated

ct_data = 0;

% 
% % Delete previous plots
% if delPreviousPlots && savePlots && exist(fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour}]), 'dir')
%     rmdir(fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour}]),'s');
% end
% 
% % Create sub-directory for each treatment if it doesn't exist
% if savePlots && ~exist(fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour}]), 'dir')
%     mkdir(fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour}]));
% end
% DirPlots_treatment = fullfile(DirPlots, [pattern{ct_pattern} '_' light{ct_light} '_' behaviour{ct_behaviour}]);

           
relevantTreatments = treatments;           
for ct_treatment=1:length(relevantTreatments)
    treatment = relevantTreatments(ct_treatment);
    
    if ~isempty(treatment.landingTracks)
        disp(['Into' ', pattern: ' treatment.pattern]);
        
        for ct_track=1:length(treatment.landingTracks) % for each landing track
            for ct_excerpt=1:length(treatment.landingTracks(ct_track).state_LDF) % for each track excerpt
                treatment.landingTracks(ct_track).track_subset_sysID(ct_excerpt) = trackForLandingModel();
                
                if ~isempty(treatment.landingTracks(ct_track).DataGUI) && ...
                        abs(treatment.landingTracks(ct_track).DataGUI(ct_excerpt).y_mode1(1,1)) >= 1e-3
                    
                    disp(['Analysing, track: ' num2str(ct_track) ...
                        ', excerpt: ' num2str(ct_excerpt)]);
                    track = treatment.landingTracks(ct_track);
                    excerpt = treatment.landingTracks(ct_track).state_LDF(ct_excerpt);
                    DataGUI = treatment.landingTracks(ct_track).DataGUI(ct_excerpt);
                    
                    nSegments = sum(DataGUI.y_mode1(:,1)~=0);
                    
                    time4rrefEstimate = [];
                    % Estimate rref and compute parameters
                    % Finding state4rrefEstimate
                    for ct=1:nSegments % For each identified segment
                        y_start_indx = find(abs(excerpt.filteredState(:,3)-DataGUI.y_mode1(ct,1))<=1e-4 & abs(excerpt.filteredState(:,6)-DataGUI.Vgy_mode1(ct,1))<=1e-4,1);
                        y_end_indx = find(abs(excerpt.filteredState(:,3)-DataGUI.y_mode1(ct,2))<=1e-4 & abs(excerpt.filteredState(:,6)-DataGUI.Vgy_mode1(ct,2))<=1e-4,1);
                        if isempty(y_start_indx) || isempty(y_end_indx)
                            warning('Estimate rref: start and end indx not found! :/');
                            continue;
                        end
                        
                        state4rrefEstimate = excerpt.filteredState(y_start_indx:y_end_indx, :);
                        time4rrefEstimate = [time4rrefEstimate; state4rrefEstimate(:,1)];
                        
                        if abs(max(state4rrefEstimate(:,3))-min(state4rrefEstimate(:,3)))>0.02
                            track.track_subset_sysID(ct_excerpt).data4rrefEstimate(end+1) = data4rrefEstimate();
                            track.track_subset_sysID(ct_excerpt).data4rrefEstimate(end).state4rrefEstimate = state4rrefEstimate;
                        end
                    end
                    
                    
                    % Run computation of parameters
                    if ~isempty(track.track_subset_sysID(ct_excerpt).data4rrefEstimate)
                        plotHandles = track.track_subset_sysID(ct_excerpt).estimate_rref(savePlots, excerpt);
                        
                        if savePlots && ~isempty(plotHandles)
                            % Resizing the figures
                            for i=1:length(plotHandles)
                                plotHandles(i).Position(3) = 680;
                                plotHandles(i).Position(4) = 545;
                                
                                if i==1
                                    figureName = ['rrefEstimate_' num2str(treatment.datenum) '_' ...
                                        num2str(treatment.startTime) '_' num2str(treatment.endTime) ...
                                        '_obj' num2str(track.obj_id) '_track' ...
                                        num2str(ct_track) '_excerpt' num2str(ct_excerpt) ...
                                        '.png'];
                                elseif i==2
                                    figureName = ['zoom_rrefEstimate_' num2str(treatment.datenum) '_' ...
                                        num2str(treatment.startTime) '_' num2str(treatment.endTime) ...
                                        '_obj' num2str(track.obj_id) '_track' ...
                                        num2str(ct_track) '_excerpt' num2str(ct_excerpt) ...
                                        '.png'];
                                end
                                
                                saveas(plotHandles(i), fullfile(DirPlots_treatment, figureName) ,'png');
                            end
                            
                            close(plotHandles);
                        end
                        %                           keyboard;
                        
                        
                        % Store them by reference for
                        % plotting later
                        ct_data = ct_data + 1;
                        data(ct_data).dataPerTrackExcerpt = track.track_subset_sysID(ct_excerpt).data4rrefEstimate;
                        data(ct_data).dataTrackExcerpt = excerpt;
                        data(ct_data).track = track;
                        
                        
                        %                                       if abs(data(ct_data).dataPerTrackExcerpt(1).dof_analytical-data(ct_data).dataPerTrackExcerpt(1).dof_actual)>1/175
                        %                                           error('Found it!!!');
                        %                                       end
                    end
                    
                    
                end
            end
        end
    end
end



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

ystart = arrayfun(@(x) (x.state4rrefEstimate(1,3)), all_data4rrefEstimates)';
speedmean = arrayfun(@(x) mean((sum((x.state4rrefEstimate(:,5:7)).^2, 2)).^0.5), all_data4rrefEstimates)';
r_basedon_speedmean = arrayfun(@(x) mean((sum((x.state4rrefEstimate(:,5:7)).^2, 2)).^0.5./x.state4rrefEstimate(:,3)), all_data4rrefEstimates)';

%%
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
ylabel('Estimated set-point of optical expansion rate, r* (s-1)', 'FontSize', 16);
xlabel('Distance from the platform, y* (m)', 'FontSize', 16);
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

%% Find parameters for automatic r( extraction algorithm

close all;
clc;
% rdotPlot = figure; hold on;
rdot_threshold = 5;

linearfit_Vvsy = zeros(length(all_data4rrefEstimates),2);
linearfit_rvsy = zeros(length(all_data4rrefEstimates),2);
linearfit_rvsy_firsthalf = zeros(length(all_data4rrefEstimates),2);
linearfit_rvsy_secondhalf = zeros(length(all_data4rrefEstimates),2);
linearfit_rvst = zeros(length(all_data4rrefEstimates),2);
rdot_mean = zeros(length(all_data4rrefEstimates),1);

ay_mean = zeros(length(all_data4rrefEstimates),1);
ay_mean_firsthalf = zeros(length(all_data4rrefEstimates),1);
ay_mean_secondhalf = zeros(length(all_data4rrefEstimates),1);

y_mean = zeros(length(all_data4rrefEstimates),1);
v_mean = zeros(length(all_data4rrefEstimates),1);
rdot_mean_firsthalf = zeros(length(all_data4rrefEstimates),1);
rdot_mean_secondhalf = zeros(length(all_data4rrefEstimates),1);
isrdot_lower_than_threshold = false(length(all_data4rrefEstimates),1);
for ct=1:length(all_data4rrefEstimates)
    state = all_data4rrefEstimates(ct).state4rrefEstimate;
    r = state(:,6)./state(:,3);
    y = state(:,3);
    rdot = abs(state(:,9)./state(:,3)-state(:,6).^2./state(:,3).^2);
    ay = state(:,9);
    rdot_mean(ct) = mean(rdot);
    ay_mean(ct) = mean(ay);
    y_mean(ct) = mean(state(:,3));
    v_mean(ct) = mean(state(:,6));
    isrdot_lower_than_threshold(ct) = all(abs(rdot)<rdot_threshold);
    if rem(size(state,1),2) == 0
        rdot_mean_firsthalf(ct) = mean(rdot(1:length(rdot)/2));
        rdot_mean_secondhalf(ct) = mean(rdot(length(rdot)/2+1:end));
        
        ay_mean_firsthalf(ct) = mean(ay(1:length(ay)/2));
        ay_mean_secondhalf(ct) = mean(ay(length(ay)/2+1:end));
        
        linearfit_rvsy_firsthalf(ct,:) = [ones(length(y)/2,1) y(1:length(y)/2)]\(r(1:length(r)/2));
        linearfit_rvsy_secondhalf(ct,:) = [ones(length(y)/2,1) y(length(y)/2+1:end)]\(r(length(r)/2+1:end));
    else
        rdot_mean_firsthalf(ct) = mean(rdot(1:floor(length(rdot)/2)+1));
        rdot_mean_secondhalf(ct) = mean(rdot(floor(length(rdot)/2)+1:end));
        
        ay_mean_firsthalf(ct) = mean(ay(1:floor(length(ay)/2)+1));
        ay_mean_secondhalf(ct) = mean(ay(floor(length(ay)/2)+1:end));
        
        linearfit_rvsy_firsthalf(ct,:) = [ones(floor(length(y)/2)+1,1) y(1:floor(length(y)/2)+1)]\(r(1:floor(length(r)/2)+1));
        linearfit_rvsy_secondhalf(ct,:) = [ones(floor(length(y)/2)+1,1) y(floor(length(y)/2)+1:end)]\(r(floor(length(r)/2)+1:end));
    end
    
    
    m1 = fitlm(state(:,3), state(:,6));
    m2 = fitlm(state(:,3), state(:,6)./state(:,3));
    m3 = fitlm(state(:,1)-state(1,1), state(:,6)./state(:,3));
%     m3 = fitlm(1:size(state,1), state(:,6)./state(:,3));
    
    linearfit_Vvsy(ct,:) = m1.Coefficients.Estimate;
    linearfit_rvsy(ct,:) = m2.Coefficients.Estimate;
    linearfit_rvst(ct,:) = m3.Coefficients.Estimate;
    
    
%     disp([linearfit_Vvsy(ct,2) linearfit_rvsy(ct,1) linearfit_rvst(ct,1) all_data4rrefEstimates(ct).rref])
    
%     figure(rdotPlot);
%     plot(state(:,1)-state(1,1), rdot, 'b')
end


%%% For chosen distributions when mean=0, sigma1=0.53, sigma2=4.22
ptc = [makedist('tLocationScale','mu',0,'sigma',0.53,'nu',5.02);
       makedist('tLocationScale','mu',0,'sigma',2*0.53,'nu',3.24); 
       makedist('tLocationScale','mu',0,'sigma',2*0.53,'nu',3.92);
       makedist('tLocationScale','mu',0,'sigma',4.22,'nu',2.3);
       makedist('tLocationScale','mu',0,'sigma',2*4.22,'nu',1.82);
       makedist('tLocationScale','mu',0,'sigma',2*4.22,'nu',1.64)];

% figure;
% subplot(2,1,1);
% ht = histfit([rref_meanVbyy-linearfit_rvsy(:,1)], [], 'tLocationScale');
% 
% subplot(2,1,2);
% hn = histfit([rref_meanVbyy-linearfit_rvsy(:,1)], [], 'Normal');
figure;
histogram([rref_meanVbyy-rref_fitVvsy])

figure; hold on;
h = histfit([rref_meanVbyy-linearfit_rvsy(:,1)], [], 'tLocationScale'); 
xlabel('r_{[mean, 0-1]}-c_{[0-1]}', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);
xlim([-4 4]); ylim([0 80]);
plot(h(2).XData, pdf(ptc(1), h(2).XData)*sum(h(1).YData*(diff(h(1).XData(1:2)))),'g','Linewidth',2);

figure; hold on;
h = histfit([rref_meanVbyy-linearfit_rvsy_firsthalf(:,1)], [], 'tLocationScale');
xlabel('r_{[mean, 0-1]}-c_{[0-0.5]}', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);
% xlim([-8 8]); ylim([0 80]);
plot(h(2).XData, pdf(ptc(2), h(2).XData)*sum(h(1).YData*(diff(h(1).XData(1:2)))),'g','Linewidth',2);

figure; hold on;
h = histfit([rref_meanVbyy-linearfit_rvsy_secondhalf(:,1)], [], 'tLocationScale');
xlabel('r_{[mean, 0-1]}-c_{[0.5-1]}', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);
% xlim([-8 8]); ylim([0 80]);
plot(h(2).XData, pdf(ptc(3), h(2).XData)*sum(h(1).YData*(diff(h(1).XData(1:2)))),'g','Linewidth',2);

figure; hold on;
h = histfit(linearfit_rvsy(:,2), [], 'tLocationScale');
xlabel('m_{[0-1]}', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);
% xlim([-50 50]); ylim([0 200]);
plot(h(2).XData, pdf(ptc(4), h(2).XData)*sum(h(1).YData*(diff(h(1).XData(1:2)))),'g','Linewidth',2);

figure; hold on;
h = histfit(linearfit_rvsy_firsthalf(:,2), [], 'tLocationScale');
xlabel('m_{[0-0.5]}', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);
% xlim([-100 100]); ylim([0 200]);
plot(h(2).XData, pdf(ptc(5), h(2).XData)*sum(h(1).YData*(diff(h(1).XData(1:2)))),'g','Linewidth',2);

figure; hold on;
h = histfit(linearfit_rvsy_secondhalf(:,2), [], 'tLocationScale');
xlabel('m_{[0.5-1]}', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);
% xlim([-100 100]); ylim([0 200]);
plot(h(2).XData, pdf(ptc(6), h(2).XData)*sum(h(1).YData*(diff(h(1).XData(1:2)))),'g','Linewidth',2);

%%
close all;
figure;
subplot(2,2,1);
h = histogram([rref_meanVbyy-linearfit_rvsy(:,1)]); hold on;
histfit([rref_meanVbyy-linearfit_rvsy(:,1)], [], 'tLocationScale');
xlabel('r_{mean}-c_{r vs y}', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);
% pd1 = fitdist([rref_meanVbyy-linearfit_rvsy(:,1)],'Normal')
pd1 = fitdist([rref_meanVbyy-linearfit_rvsy(:,1)],'tLocationScale')

subplot(2,2,2);
histogram([rref_meanVbyy-linearfit_rvsy_firsthalf(:,1)]);
histfit([rref_meanVbyy-linearfit_rvsy_firsthalf(:,1)], [], 'tLocationScale');
xlabel('r_{mean}-c_{r vs y, 0-0.5}', 'FontSize', 16);
set(gca, 'FontSize', 16);
% pd2 = fitdist([rref_meanVbyy-linearfit_rvsy_firsthalf(:,1)],'Normal')
pd2 = fitdist([rref_meanVbyy-linearfit_rvsy_firsthalf(:,1)],'tLocationScale')

subplot(2,2,3);
histogram([rref_meanVbyy-linearfit_rvsy_secondhalf(:,1)]);
histfit([rref_meanVbyy-linearfit_rvsy_secondhalf(:,1)], [], 'tLocationScale');
xlabel('r_{mean}-c_{r vs y, 0.5-1}', 'FontSize', 16);
set(gca, 'FontSize', 16);
% pd3 = fitdist([rref_meanVbyy-linearfit_rvsy_secondhalf(:,1)],'Normal')
pd3 = fitdist([rref_meanVbyy-linearfit_rvsy_secondhalf(:,1)],'tLocationScale')

figure;
subplot(2,2,1);
histogram(linearfit_rvsy(:,2));
histfit(linearfit_rvsy(:,2), [], 'tLocationScale');
xlabel('m_{r vs y}', 'FontSize', 16);
set(gca, 'FontSize', 16);
% pd4 = fitdist(linearfit_rvsy(:,2),'Normal')
pd4 = fitdist(linearfit_rvsy(:,2),'tLocationScale')

subplot(2,2,2);
histogram(linearfit_rvsy_firsthalf(:,2));
histfit(linearfit_rvsy_firsthalf(:,2), [], 'tLocationScale');
xlabel('m_{r vs y, 0-0.5}', 'FontSize', 16);
set(gca, 'FontSize', 16);
% pd5 = fitdist(linearfit_rvsy_firsthalf(:,2),'Normal')
pd5 = fitdist(linearfit_rvsy_firsthalf(:,2),'tLocationScale')

subplot(2,2,3);
histogram(linearfit_rvsy_secondhalf(:,2));
histfit(linearfit_rvsy_secondhalf(:,2), [], 'tLocationScale');
xlabel('m_{r vs y, 0.5-1}', 'FontSize', 16);
set(gca, 'FontSize', 16);
% pd6 = fitdist(linearfit_rvsy_secondhalf(:,2),'Normal')
pd6 = fitdist(linearfit_rvsy_secondhalf(:,2),'tLocationScale')
