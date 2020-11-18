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

%% Loading data and collecting segments with rref
clc; close all;
% clear;
% 
if isunix
    inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A1_rref.mat';
elseif ispc
    inputFile = 'D:/light_intensity_experiments/postprocessing/BlindLandingtracks_A1_rref.mat';
end
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
            rrefSegs = [state_LDF.rrefSegments];
            
            treatment_indx4stateLDF = arrayfun(@(x) x*ones(length([relevantTreatments(x).landingTracks.state_LDF]),1),1:length(relevantTreatments), 'UniformOutput', false);
            treatment_indx4stateLDF = vertcat(treatment_indx4stateLDF{:});
            
            landingTracks_indx4stateLDF = arrayfun(@(x) x*ones(length(landingTracks(x).state_LDF),1),1:length(landingTracks), 'UniformOutput', false);
            landingTracks_indx4stateLDF = vertcat(landingTracks_indx4stateLDF{:});
            
            % Display info about # of landing tracks
            disp(['# of distinct Flydra objects (landingTracks): ' num2str(length(landingTracks))]);
            disp(['# of state LDFs (landingTracks): ' num2str(length(state_LDF))]);
%             disp(['Size of rrefSegment vector: ' num2str(length(rrefSegs))]);
            
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
            dummy.ct_pattern = ct_pattern;
            dummy.ct_light = ct_light;
            dummy.landingTrack = landingTracks(landingTracks_indx4stateLDF(indices));
            dummy.landingDiscs = vertcat(relevantTreatments(treatment_indx4stateLDF(indices)).landingDiscs);
            data = [data; dummy];
            
            
%             for ct_treatment=1:length(relevantTreatments) % for each relevant treatment
%                 treatment = relevantTreatments(ct_treatment);
%                 arrayfun(@(x) x.setLandingSide(),[treatment.landingTracks.state_LDF]'); % to store landing side in the rrefSegments
%                 data1 = arrayfun(@(x) x.rrefSegments,[treatment.landingTracks.state_LDF]','UniformOutput',false);
%                 data1 = horzcat(data1{:});
%                 
%                 % Discard empty intervals
%                 indices = arrayfun(@(x) ~isempty(x.intervals_ti), data1);
%                 
%                 % Save additional data
%                 data1 = data1(indices);
%                 [data1.pattern] = deal(ct_pattern);
%                 [data1.light] = deal(ct_light);
%                 [data1.day] = deal(treatment.datenum);
%                 [data1.time] = deal(treatment.startTime);
%     
%                 data = [data data1];
% %                 size(data)
% %                 keyboard;
%             end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end

%% If landing tack started from the opposite disc
%% Plot histogram of y_start and V_start for landing approaches
% data is only for chosen factor
close all;
y_start = cell(length(data),1);
V_start = cell(length(data),1);
V3d_start = cell(length(data),1);
for ct=1:length(data)
    data(ct).hasTakeoff = false(1,length(data(ct).tracks_fac));
    for ct1=1:length(data(ct).tracks_fac)
        data(ct).hasTakeoff(ct1) = data(ct).tracks_fac(ct1).hasTakeoff(data(ct).landingDiscs(ct1,:));
        
        [y_start{ct}(ct1,1), ymin_indx] = min(data(ct).tracks_fac(ct1).filteredState(:,3));
        V_start{ct}(ct1,1) = data(ct).tracks_fac(ct1).filteredState(ymin_indx,6);
        V3d_start{ct}(ct1,1) = (sum(data(ct).tracks_fac(ct1).filteredState(ymin_indx,5:7).^2))^0.5;
    end
%     figure;
%     histogram(-vertcat(y_start{ct}));
end

V3d = vertcat(V3d_start{:});
V = vertcat(V_start{:});

%% Plots with takeoff
close all;
% starting velocity
figure;
title('Golden flights', 'FontSize', 16);
subplot(1,2,1);
histogram(V([data.hasTakeoff]));
xlabel('V_0 (m/s)', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);

subplot(1,2,2);
histogram(V3d([data.hasTakeoff]));
xlabel('V3d_0 (m/s)', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);

% Create trajectory views
% close all;
trajPlot = figure; hold on;
skip_step = 10;
N = 0;
for ct=1:length(data)
    tracks = data(ct).tracks_fac(data(ct).hasTakeoff);
    for ct1 = 1:skip_step:length(tracks)
        xyz = tracks(ct1).filteredState(:,[2 3 4]); % the complete trajectory

        plot3(xyz(:,1), -1*xyz(:,2), xyz(:,3),'Color',[252,187,161]./255);

    end
    N = N + length(1:skip_step:length(tracks));
end

for ct=1:length(data)
    tracks = data(ct).tracks_fac(data(ct).hasTakeoff);
    for ct1 = 1:skip_step:length(tracks)
        intervals = tracks(ct1).rrefSegments(...
            abs([tracks(ct1).rrefSegments.factor]-chosen_fac)<1e-6).intervals_ti;
        for ct2 = 1:size(intervals,1)
            xyz4rref = tracks(ct1).filteredState(intervals(ct2,1):intervals(ct2,2),[2 3 4]); % trajectory used for estimating r*
            plot3(xyz4rref(:,1), -1*xyz4rref(:,2), xyz4rref(:,3),'Color',[215 48 39]./255, 'Linewidth', 1);
        end
    end
end

% % Landing Disc
% [X,Y,Z] = cylinder(treatment.landingDiscs(1).radius);
radius = treatments(1).landingDiscs(1).radius;

figure(trajPlot);
% h=mesh(X,Y,Z,'facecolor',[1 0 0]); % draw landing disc
plot3([-radius; radius], [0 0], [0 0], 'LineWidth', 2, 'Color', [83 83 83]./255);
plot3([0 0], [0 0], [-radius; radius], 'LineWidth', 2, 'Color', [83 83 83]./255);
zlabel('z (m)', 'FontSize', 14);
ylabel('y (m)', 'FontSize', 14);
xlabel('x (m)', 'FontSize', 14);
title('With take-off', 'FontSize', 16);
axis equal;
xlim([-0.4 0.4]);
xticks([-0.4:0.2:0.4]);
ylim([0 0.4]);
yticks([0:0.1:0.4]);
zlim([-0.25 0.25]);
zticks([-0.25:0.25:0.25]);
set(gca, 'FontSize', 16);
view(0,90);
view(-90,0);

% % % Plot ymean vs rmean
% close all;
% from analysis in R
taudots = [-0.8443717 -0.8443717 -0.8443717;
           -0.7799524 -0.7799524 -0.7799524];
% taudot_checker = -0.847; 
% taudot = -0.780; % from analysis in R

intercepts = [-1.117776 -1.045171 -0.9404179;
              -0.9693891 -0.8967832 -0.805922];  

points_cmap = [252,187,161;
200,200,200]./255;

line_cmap = [215,48,31;
37,37,37]./255;

rmean = cell(1,length(pattern));
ymean = cell(1,length(pattern));
for ct=1:length(pattern)
    data_ss = vertcat(data([data.ct_pattern]==ct));
    
    rmean{ct} = [];
    ymean{ct} = [];
    for ct_light=1:length(data_ss)
        dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).rmean_ti,data_ss(ct_light).tracks_fac(data_ss(ct_light).hasTakeoff),'UniformOutput',false);
        rmean{ct} = [rmean{ct}; vertcat(dummy{:})];
        
        dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).ymean_ti,data_ss(ct_light).tracks_fac(data_ss(ct_light).hasTakeoff),'UniformOutput',false);
        ymean{ct} = [ymean{ct}; vertcat(dummy{:})];
    end
end
data2plot = [vertcat(ymean{:}) vertcat(rmean{:}) [ones(size(rmean{1})); 1*ones(size(rmean{2}))]];
data2plot = data2plot(randperm(size(data2plot, 1)), :);
figure; hold on;
% % % plot lines
% % y_vec = min(-ymean{1}):0.001:max(-ymean{1});
% % modelfun = @(b,x)(exp(b(1))*x.^(b(2)));
% % for ct=1:length(pattern)
% %     Coefficients = [mean(intercepts(ct,:)) mean(taudots(ct,:))];
% %     plot(y_vec,modelfun(Coefficients,y_vec),'Color', line_cmap(ct,:),'LineWidth',3);
% % end
% % legend({'checkerboard','spoke'})
set(gca, 'FontSize', 16);
ylabel('Estimated set-points, r* (1/s)', 'FontSize', 16);
xlabel('y* (m)', 'FontSize', 16);
title('With take-off', 'FontSize', 16);

scatter(-data2plot(:,1),-data2plot(:,2),10,points_cmap(round(data2plot(:,3)),:),'filled','o');
% % for ct=1:length(pattern)
% %     Coefficients = [mean(intercepts(ct,:)) mean(taudots(ct,:))];
% %     plot(y_vec,modelfun(Coefficients,y_vec),'Color', line_cmap(ct,:),'LineWidth',3);
% % end


% create log plots
figure;
scatter(log(-data2plot(:,1)),log(-data2plot(:,2)),10,points_cmap(round(data2plot(:,3)),:),'filled','o');
xlim([-4 -1]);
xticks([-4:1:-1]);
ylim([-4 3]);
yticks([-4:1:3]);
set(gca, 'FontSize', 16);
ylabel('log(r*) (1/s)', 'FontSize', 16);
xlabel('log(y*) (m)', 'FontSize', 16);
title('With take-off', 'FontSize', 16);
%% Plots without takeoff
close all;
% starting velocity
figure;
title('Sideways flights', 'FontSize', 16);
subplot(1,2,1);
histogram(V(~[data.hasTakeoff]));
xlabel('V_0 (m/s)', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);

subplot(1,2,2);
histogram(V3d(~[data.hasTakeoff]));
xlabel('V3d_0 (m/s)', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);


% Create trajectory views
% close all;
trajPlot = figure; hold on;
skip_step = 10;
N = 0;
for ct=1:length(data)
    tracks = data(ct).tracks_fac(~data(ct).hasTakeoff);
    for ct1 = 1:skip_step:length(tracks)
        xyz = tracks(ct1).filteredState(:,[2 3 4]); % the complete trajectory

        plot3(xyz(:,1), -1*xyz(:,2), xyz(:,3),'Color',[252,187,161]./255);

    end
    N = N + length(1:skip_step:length(tracks));
end

for ct=1:length(data)
    tracks = data(ct).tracks_fac(~data(ct).hasTakeoff);
    for ct1 = 1:skip_step:length(tracks)
        intervals = tracks(ct1).rrefSegments(...
            abs([tracks(ct1).rrefSegments.factor]-chosen_fac)<1e-6).intervals_ti;
        for ct2 = 1:size(intervals,1)
            xyz4rref = tracks(ct1).filteredState(intervals(ct2,1):intervals(ct2,2),[2 3 4]); % trajectory used for estimating r*
            plot3(xyz4rref(:,1), -1*xyz4rref(:,2), xyz4rref(:,3),'Color',[215 48 39]./255, 'Linewidth', 1);
        end
    end
end

% % Landing Disc
% [X,Y,Z] = cylinder(treatment.landingDiscs(1).radius);
radius = treatments(1).landingDiscs(1).radius;

figure(trajPlot);
% h=mesh(X,Y,Z,'facecolor',[1 0 0]); % draw landing disc
plot3([-radius; radius], [0 0], [0 0], 'LineWidth', 2, 'Color', [83 83 83]./255);
plot3([0 0], [0 0], [-radius; radius], 'LineWidth', 2, 'Color', [83 83 83]./255);
zlabel('z (m)', 'FontSize', 14);
ylabel('y (m)', 'FontSize', 14);
xlabel('x (m)', 'FontSize', 14);
title('Without take-off', 'FontSize', 16);
axis equal;
xlim([-0.4 0.4]);
xticks([-0.4:0.2:0.4]);
ylim([0 0.4]);
yticks([0:0.1:0.4]);
zlim([-0.25 0.25]);
zticks([-0.25:0.25:0.25]);
set(gca, 'FontSize', 16);
view(0,90);
view(-90,0);


% % % Plot ymean vs rmean
% close all;
% from analysis in R
taudots = [-0.8443717 -0.8443717 -0.8443717;
           -0.7799524 -0.7799524 -0.7799524];
% taudot_checker = -0.847; 
% taudot = -0.780; % from analysis in R

intercepts = [-1.117776 -1.045171 -0.9404179;
              -0.9693891 -0.8967832 -0.805922];  

points_cmap = [252,187,161;
200,200,200]./255;

line_cmap = [215,48,31;
37,37,37]./255;

rmean = cell(1,length(pattern));
ymean = cell(1,length(pattern));
for ct=1:length(pattern)
    data_ss = vertcat(data([data.ct_pattern]==ct));
    
    rmean{ct} = [];
    ymean{ct} = [];
    for ct_light=1:length(data_ss)
        dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).rmean_ti,data_ss(ct_light).tracks_fac(~data_ss(ct_light).hasTakeoff),'UniformOutput',false);
        rmean{ct} = [rmean{ct}; vertcat(dummy{:})];
        
        dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).ymean_ti,data_ss(ct_light).tracks_fac(~data_ss(ct_light).hasTakeoff),'UniformOutput',false);
        ymean{ct} = [ymean{ct}; vertcat(dummy{:})];
    end
end
data2plot = [vertcat(ymean{:}) vertcat(rmean{:}) [ones(size(rmean{1})); 1*ones(size(rmean{2}))]];
data2plot = data2plot(randperm(size(data2plot, 1)), :);
figure; hold on;
% % % plot lines
% % y_vec = min(-ymean{1}):0.001:max(-ymean{1});
% % modelfun = @(b,x)(exp(b(1))*x.^(b(2)));
% % for ct=1:length(pattern)
% %     Coefficients = [mean(intercepts(ct,:)) mean(taudots(ct,:))];
% %     plot(y_vec,modelfun(Coefficients,y_vec),'Color', line_cmap(ct,:),'LineWidth',3);
% % end
% % legend({'checkerboard','spoke'})
set(gca, 'FontSize', 16);
ylabel('Estimated set-points, r* (1/s)', 'FontSize', 16);
xlabel('y* (m)', 'FontSize', 16);
title('Without take-off', 'FontSize', 16);

scatter(-data2plot(:,1),-data2plot(:,2),10,points_cmap(round(data2plot(:,3)),:),'filled','o');
% % for ct=1:length(pattern)
% %     Coefficients = [mean(intercepts(ct,:)) mean(taudots(ct,:))];
% %     plot(y_vec,modelfun(Coefficients,y_vec),'Color', line_cmap(ct,:),'LineWidth',3);
% % end

% create log plots
figure;
scatter(log(-data2plot(:,1)),log(-data2plot(:,2)),10,points_cmap(round(data2plot(:,3)),:),'filled','o');
xlim([-4 -1]);
xticks([-4:1:-1]);
ylim([-4 3]);
yticks([-4:1:3]);
set(gca, 'FontSize', 16);
ylabel('log(r*) (1/s)', 'FontSize', 16);
xlabel('log(y*) (m)', 'FontSize', 16);
title('Without take-off', 'FontSize', 16);

%%
figure;
histogram(-vertcat(y_start{:}));
xlabel('y_0 (m)', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);

figure;
histogram(vertcat(V_start{:}));
xlabel('V_0 (m/s)', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);

figure;
histogram(vertcat(V3d_start{:}));
xlabel('V3d_0 (m/s)', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);

figure;
hist3([-vertcat(y_start{:}) vertcat(V_start{:})]);
xlabel('y_0 (m)', 'FontSize', 16);
ylabel('V_0 (m/s)', 'FontSize', 16);
zlabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);
view([-146 19])

figure;
hist3([-vertcat(y_start{:}) vertcat(V3d_start{:})]);
xlabel('y_0 (m)', 'FontSize', 16);
ylabel('V3d_0 (m/s)', 'FontSize', 16);
zlabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);
view([-146 19])

[max(-vertcat(y_start{:})) min(-vertcat(y_start{:}))]
[max(vertcat(V_start{:})) min(vertcat(V_start{:}))]