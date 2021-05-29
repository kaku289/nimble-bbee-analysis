%% Comparison of honeybees and bumblebees landing strategies
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

%% Load and process landings from free-flight for bumblebees

clc; close all;
% clear;
% 
inputFile = '/media/reken001/Disk_12/light_intensity_experiments/postprocessing/BlindLandingtracks_A1_rref.mat';
% inputFile = 'D:/light_intensity_experiments/postprocessing/BlindLandingtracks_A1_rref.mat';
% load(inputFile);
treatments = treatments(1:14*8); % Taking experiments for 2 patterns * 3 lights


pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
behaviour = {'rising','constant','sleeping'};

factors = [0.25:0.25:2.5];
chosen_fac = 1;

data = struct.empty;
clear dummy N N_fromTakeoff N_fromFreeflight N_stateLDFs;
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
            
            treatment_indx4landingTracks = arrayfun(@(x) x*ones(length(relevantTreatments(x).landingTracks),1),1:length(relevantTreatments), 'UniformOutput', false);
            treatment_indx4landingTracks = vertcat(treatment_indx4landingTracks{:});
            treatment_indx4stateLDF = arrayfun(@(x) treatment_indx4landingTracks(x)*ones(length(landingTracks(x).state_LDF),1),1:length(landingTracks), 'UniformOutput', false);
            treatment_indx4stateLDF = vertcat(treatment_indx4stateLDF{:});
            
            landingTracks_indx4stateLDF = arrayfun(@(x) x*ones(length(landingTracks(x).state_LDF),1),1:length(landingTracks), 'UniformOutput', false);
            landingTracks_indx4stateLDF = vertcat(landingTracks_indx4stateLDF{:});
            
            % Display info about # of landing tracks
            disp(['# of distinct Flydra objects (landingTracks): ' num2str(length(landingTracks))]);
            disp(['# of state LDFs (landingTracks): ' num2str(length(state_LDF))]);
%             disp(['Size of rrefSegment vector: ' num2str(length(rrefSegs))]);
            
            % Finding # of "tracks" (state LDFs) that contain rref segments
            % for each factor
            N{ct_light,ct_pattern} = zeros(1, length(factors));
            landingDiscs = {relevantTreatments(treatment_indx4stateLDF).landingDiscs};
            hastakeoff = arrayfun(@(x) state_LDF(x).hasTakeoff(landingDiscs{x}),1:length(state_LDF));
            N_stateLDFs{ct_light,ct_pattern} = [length(state_LDF), sum(hastakeoff), sum(~hastakeoff)];
            N_fromTakeoff{ct_light,ct_pattern} = zeros(1, length(factors));
            N_fromFreeflight{ct_light,ct_pattern} = zeros(1, length(factors));
            for ct_fac = 1:length(factors)
                factor = factors(ct_fac);
                indices = arrayfun(@(x) ~isempty(x.rrefSegments(abs([x.rrefSegments.factor]-factor)<1e-6).intervals_ti),state_LDF);
                N{ct_light,ct_pattern}(ct_fac) = sum(indices);
                N_fromTakeoff{ct_light,ct_pattern}(ct_fac) = sum(indices & hastakeoff);
                N_fromFreeflight{ct_light,ct_pattern}(ct_fac) = sum(indices & ~hastakeoff);
            end
            
            % Display 
            disp(['# of state LDFs containing rrefs for different factors: ' num2str(N{ct_light,ct_pattern})]);
            disp(['# of state LDFs containing rrefs for different factors from Takeoff: ' num2str(N_fromTakeoff{ct_light,ct_pattern})]);
            disp(['# of state LDFs containing rrefs for different factors from Freeflight: ' num2str(N_fromFreeflight{ct_light,ct_pattern})]);
            
            % For chosen factor, collect all tracks that do contain
            % non-empty rref segments
            indices = arrayfun(@(x) ~isempty(x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).intervals_ti), state_LDF);
            dummy.tracks_fac = state_LDF(indices);
            dummy.ct_pattern = ct_pattern;
            dummy.ct_light = ct_light;
            dummy.landingTrack = landingTracks(landingTracks_indx4stateLDF(indices));
            treatment_indx4stateLDF = treatment_indx4stateLDF(indices);
%             dummy.landingDiscs = {relevantTreatments(treatment_indx4stateLDF).landingDiscs};
            
            dummy.hastakeoff = hastakeoff(indices); %arrayfun(@(x) dummy.tracks_fac(x).hasTakeoff(dummy.landingDiscs{x}),1:length(dummy.tracks_fac));
            
            data = [data; dummy];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end

bbee_data = data;
for ct=1:length(bbee_data)
    hastakeoff = bbee_data(ct).hastakeoff;
    
    dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).rref_ti,bbee_data(ct).tracks_fac(~hastakeoff),'UniformOutput',false);
    bbee_data(ct).rref = vertcat(dummy{:});
    
    dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).rmean_ti,bbee_data(ct).tracks_fac(~hastakeoff),'UniformOutput',false);
    bbee_data(ct).rmean = vertcat(dummy{:});

    dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).ymean_ti,bbee_data(ct).tracks_fac(~hastakeoff),'UniformOutput',false);
    bbee_data(ct).ymean = vertcat(dummy{:});
    
end
%% Load and process landings of honeybes
clc; close all;

inputFile = '/media/reken001/Disk_12/honeybee_experiments/postprocessing/BlindLandingtracks_A4_LDF_rref.mat';
load(inputFile);

landing_tracks = [landingTracks{:}];
static_patternnums = [1 2 3 4 5 6 10 11 18];
landingTracks = landing_tracks(arrayfun(@(x) ismember(x,static_patternnums), [landing_tracks.patternnum]));

factors = [0.25:0.25:2.5];
chosen_fac = 1.5;

data = struct.empty;
clear dummy;

state_LDF = [landingTracks.state_LDF];
rrefSegs = [state_LDF.rrefSegments];

landingTracks_indx4stateLDF = arrayfun(@(x) x*ones(length(landingTracks(x).state_LDF),1),1:length(landingTracks), 'UniformOutput', false);
landingTracks_indx4stateLDF = vertcat(landingTracks_indx4stateLDF{:});

% For chosen factor, collect all tracks that do contain
% non-empty rref segments
indices = arrayfun(@(x) ~isempty(x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).intervals_ti), state_LDF);
dummy.tracks_fac = state_LDF(indices);
dummy.landingTrack = landingTracks(landingTracks_indx4stateLDF(indices));
data = [data; dummy];

hbee_data = data;

for ct=1:length(hbee_data)
    dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).rref_ti,hbee_data(ct).tracks_fac,'UniformOutput',false);
    hbee_data(ct).rref = vertcat(dummy{:});
    
    dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).rmean_ti,hbee_data(ct).tracks_fac,'UniformOutput',false);
    hbee_data(ct).rmean = vertcat(dummy{:});

    dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).ymean_ti,hbee_data(ct).tracks_fac,'UniformOutput',false);
    hbee_data(ct).ymean = vertcat(dummy{:});
    
end

dummy = abs(vertcat(hbee_data.rref)-vertcat(hbee_data.rmean));
[mean(dummy) median(dummy) max(dummy)]


%% Make comparison plots
close all;
points_cmap = [252,187,161;
        200,200,200]./255;
    
line_cmap = [215,48,31;
        37,37,37]./255;
    
map = brewermap(3,'Set1');

% Plotting set-point histogram
figure; hold on;
hbee_histfit = histfit(-vertcat(hbee_data.rmean),13,'gamma');
hbee_pd = fitdist(-vertcat(hbee_data.rmean),'Gamma');
% plot(hbee_histfit(2).XData, pdf(hbee_pd, hbee_histfit(2).XData)*sum(hbee_histfit(1).YData*(diff(hbee_histfit(1).XData(1:2)))),'g','Linewidth',0.5);
figure;
bbee_histfit = histfit(-vertcat(bbee_data.rmean),[],'gamma');
bbee_pd = fitdist(-vertcat(bbee_data.rmean),'Gamma');

figure;
hbee_hist = histogram(-vertcat(hbee_data.rmean),13,'facecolor',points_cmap(1,:),'facealpha',.5,'edgecolor','none','Normalization','pdf');
hold on;
bbe_hist = histogram(-vertcat(bbee_data.rmean),67,'facecolor',points_cmap(2,:),'facealpha',.5,'edgecolor','none','Normalization','pdf');
plot(hbee_histfit(2).XData, pdf(hbee_pd, hbee_histfit(2).XData),'color',line_cmap(1,:),'Linewidth',2);
plot(bbee_histfit(2).XData, pdf(bbee_pd, bbee_histfit(2).XData),'color',line_cmap(2,:),'Linewidth',2);

legend('honeybees','bumblebees','fontsize',16)
xlabel('Estimated set-points, r* (s-1)', 'FontSize', 16);
ylabel('Probability density function', 'FontSize', 16);
set(gca, 'FontSize', 16);


% Plotting rmean vs ymean
figure; hold on;
ymean = -hbee_data(ct).ymean;
y_vec = 0.05:0.001:max(ymean);
modelfun = @(b,x)(exp(b(1))*x.^(b(2)));
Coefficients = [0.78618 -0.21784]; % from R
plot(y_vec,modelfun(Coefficients,y_vec),'Color', line_cmap(1,:), 'LineWidth', 2);
ymean = -bbee_data(ct).ymean;
y_vec = min(ymean):0.001:max(ymean);
modelfun = @(b,x)(exp(b(1))*x.^(b(2)));
Coefficients = [-0.81634 -0.74543]; % from R (bumblebee_landing_dynamics_free-flight.R in A4 R-statistics)
plot(y_vec,modelfun(Coefficients,y_vec),'Color', line_cmap(2,:), 'LineWidth', 2);
scatter(-hbee_data(ct).ymean,-hbee_data(ct).rref,10,points_cmap(1,:),'filled','s');
scatter(-bbee_data(ct).ymean,-bbee_data(ct).rref,10,points_cmap(2,:),'filled','s');
set(gca, 'FontSize', 16);
ylabel('Estimated set-points, r* (s-1)', 'FontSize', 16);
xlabel('Mean y (m)', 'FontSize', 16);


%% CODE below is not used

%%
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


figure;
histfit(-vertcat(hbee_data.rmean),13,'gamma')
xlabel('Estimated set-points, r* (1/s)', 'FontSize', 16);
ylabel('Occurences', 'FontSize', 16);
set(gca, 'FontSize', 16);
dummy = -vertcat(hbee_data.rmean);
[mean(dummy) median(dummy) max(dummy) min(dummy)]
pd = fitdist(dummy,'Gamma');

%%
close all;

% figure;
% histogram(-vertcat(bbee_data.rref), [0:0.5:8]);
% xlabel('Estimated set-points, r* (1/s)', 'FontSize', 16);
% ylabel('Occurances', 'FontSize', 16);
% set(gca, 'FontSize', 16);

dummy = abs(vertcat(bbee_data.rref)-vertcat(bbee_data.rmean));
[mean(dummy) median(dummy) max(dummy)]

figure;
% histogram(-vertcat(bbee_data.rmean), [0:0.5:8]);
% histogram(-vertcat(bbee_data.rmean), [0:0.5:9.5]);
histfit(-vertcat(bbee_data.rmean),[],'Gamma')
xlabel('Estimated set-points, r* (s-1)', 'FontSize', 16);
ylabel('Occurences', 'FontSize', 16);
set(gca, 'FontSize', 16);
dummy = -vertcat(bbee_data.rmean);
[mean(dummy) median(dummy) max(dummy) min(dummy)]
pd = fitdist(dummy,'Gamma');

figure; hold on;

    
% intercepts (freeflight)
c = [-1.000403   -0.923032   -0.8313437;
        -0.8138813  -0.7365104  -0.6448221];
% slopes (freeflight)
m = [-0.7769429  -0.7769429  -0.7769429;
        -0.6902597  -0.6902597  -0.6902597];
    
f = scatter(-bbee_data(ct).ymean,-bbee_data(ct).rref,10,points_cmap(1,:),'filled','s');
set(gca, 'FontSize', 16);
ylabel('Estimated set-points, r* (s-1)', 'FontSize', 16);
xlabel('Mean y (m)', 'FontSize', 16);
ymean = -bbee_data(ct).ymean;
y_vec = min(ymean):0.001:max(ymean);
modelfun = @(b,x)(exp(b(1))*x.^(b(2)));
Coefficients = [-0.81634 -0.74543]; % from R (bumblebee_landing_dynamics_free-flight.R in A4 R-statistics)
plot(y_vec,modelfun(Coefficients,y_vec),'Color', line_cmap(1,:), 'LineWidth', 2);


