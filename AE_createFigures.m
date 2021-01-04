%% Create Figures for the article

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

%% Figure 1 panels


% %%%%%%%%%%%%% Find initial conditions for V at y = 0.3m for pigeons
syms t
xdot = (0.22*t-0.057);
tau = (0.719*t-0.033);
f = tau*xdot - 0.3;
t = roots(fliplr(double(coeffs(f))))
% coeff = [0.719*0.22 -0.057*0.719-0.033*0.22 0.033*0.057-0.3]

% f = @(t) (0.719*t-0.033)*(0.22*t-0.057) + 0.3; % at x = 0.3m
% 
% 
% t = fzero(f, 1.5)

% t = [-1.23, 1.54];
double(subs(xdot))
double(subs(tau))
double(subs(f))
% %%%%%%%%%%%%%%%%%%%%%% at x=0.3m, xdot = 0.3 -> tau = r = 1
saveDir = '/home/reken001/Pulkit/graphs_temp';
close all;

fig1 = figure; hold on; pbaspect ([0.9 1 1]);
fig2 = figure; hold on; pbaspect ([0.9 1 1]);
fig3 = figure; hold on; pbaspect ([0.9 1 1]);

% Plot constant-r approach (honeybee)
Vi = 0.3; % Initial velocity
yi = 0.3; % Initial distance


figure(fig1); 
y = 0:0.001:yi;
V_const_r = y./(yi/Vi);
plot(y,V_const_r,'b', 'LineWidth', 2);

figure(fig2);
y = 0.04:0.001:yi;
r_const_r = Vi/yi*ones(length(y),1);
plot(y,r_const_r,'b', 'LineWidth', 2);

figure(fig3);
plot(y,0*ones(length(y),1),'b', 'LineWidth', 2);

% Plot constant-taudot approach (pigeons) 
Vi = 0.3; % Initial velocity
yi = 0.3; % Initial distance
taudot = -0.775;

figure(fig1); 
y = 0:0.001:yi;
V_const_taudot = y.^(1+taudot)./(yi.^(1+taudot)/Vi);
plot(y,V_const_taudot,'r', 'LineWidth', 2);

figure(fig2);
y = 0.04:0.001:yi;
r_const_taudot = y.^(taudot)./(yi.^(1+taudot)/Vi);
plot(y,r_const_taudot,'r', 'LineWidth', 2);

figure(fig3);
plot(y,taudot*ones(length(y),1),'r', 'LineWidth', 2);

% Plot bumblebees approach
ymean = [0.25 0.15 0.05]; rmean = [1.15 1.8 4.25];

figure(fig1); 
for ct=1:length(ymean)
    y = ymean(ct)-0.02:0.001:ymean(ct)+0.02;
    plot(y,rmean(ct)*y,'k', 'LineWidth', 2);
end

figure(fig2);
for ct=1:length(ymean)
    y = ymean(ct)-0.02:0.001:ymean(ct)+0.02;
    plot(y,rmean(ct)*ones(length(y),1),'k', 'LineWidth', 2);
end

figure(fig3);
for ct=1:length(ymean)
    y = ymean(ct)-0.02:0.001:ymean(ct)+0.02;
    plot(y,zeros(length(y),1),'k', 'LineWidth', 3);
end


% make figure neat and tidy
figure(fig1)
xlabel('y (m)', 'FontSize', 16);
ylabel('V (m/s)', 'FontSize', 16);
set(gca, 'FontSize', 16);
legend('Zero taudot apporach', 'Negative taudot approach', 'Hybrid');
% legend('$\dot{\tau} = 0$', '$\dot{\tau} = -0.775$','Interpreter','latex');
ylim([0 0.4]);
print(fig1, fullfile(saveDir,'General_Vvsy'), '-dpdf');

figure(fig2)
xlim([0 Inf]);
xlabel('y (m)', 'FontSize', 16);
ylabel('r (1/s)', 'FontSize', 16);
set(gca, 'FontSize', 16);
% legend('$\dot{\tau} = 0$', '$\dot{\tau} = -0.75$','Interpreter','latex');
% legend('\tau dot = 0', '\tau dot = -0.75');
print(fig2, fullfile(saveDir,'General_rvsy'), '-dpdf');

% y = 0.05:0.001:yi;
% r_const_r = Vi/yi*ones(length(y),1);
% r_const_taudot1 = y.^(-0.76)./(yi.^(1-0.76)/Vi);
% r_const_taudot2 = y.^(-0.83)./(yi.^(1-0.83)/Vi);
% fig2 = figure; hold on;
% plot(y,r_const_r,'b', 'LineWidth', 2);
% plot(y,r_const_taudot1,'r', 'LineWidth', 2);
% plot(y,r_const_taudot2,'g', 'LineWidth', 2);
% xlim([0 Inf]);
% xlabel('y (m)', 'FontSize', 16);
% ylabel('r (1/s)', 'FontSize', 16);
% set(gca, 'FontSize', 16);
% legend('$\dot{\tau} = 0$', '$\dot{\tau} = -0.76$', '$\dot{\tau} = -0.83$', 'Interpreter','latex');

figure(fig3);
xlim([0 Inf]);
xlabel('y (m)', 'FontSize', 16);
ylabel('$\dot{\tau}$', 'FontSize', 16,'Interpreter','latex');
set(gca, 'FontSize', 16);
% legend('$\dot{\tau} = 0$', '$\dot{\tau} = -0.775$', 'Interpreter','latex');
ylim([-0.8 0.2])
print(fig3, fullfile(saveDir,'General_taudotvsy'), '-dpdf');

%% Figure 2 panels

clc; close all;
% clear;
% 
% inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A1_rref_videos.mat';
% load(inputFile);
treatments = treatments(1:14*8); % Taking experiments for 2 patterns * 3 lights

DirPlots = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/plots/BlindLandingTracks/mult_rref_with_videos';
delPreviousPlots = false; % BE CAREFUL - when set to true, all previously saved plots are deleted
savePlots = false;
savePDFs = false;

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
behaviour = {'rising','constant','sleeping'};

factors = [0.25:0.25:2.5];
chosen_fac = 1;

data4video = struct.empty;
% tracks_fac(length(pattern), length(light)) = filteredState_BlindLandingtrack.empty; % tracks for chosen factor

for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        for ct_behaviour = 2%1:length(behaviour)
            disp(' ');
            disp(['Pattern: ' pattern{ct_pattern} ...
                  ', light: ' light{ct_light} ...
                  ', behaviour: ' behaviour{ct_behaviour}]);
              
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
            
            clear dummy;
            for ct_treatment=1:length(relevantTreatments)
                
                treatment = relevantTreatments(ct_treatment);
                videoTimes = [[treatment.videosInfo.startTime]' [treatment.videosInfo.endTime]'];
                
                % collect time vectors for tracks with multiple r*
                landingTracks = [treatment.landingTracks];
                state_LDF = [landingTracks.state_LDF];
                has_multiple_rrefs = arrayfun(@(x) size(x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).intervals_ti, 1) > 1, state_LDF);
                state_LDF_multiple_rrefs = state_LDF(has_multiple_rrefs);
                
                landingTracks_indx4stateLDF = arrayfun(@(x) x*ones(length(landingTracks(x).state_LDF),1),1:length(landingTracks), 'UniformOutput', false);
                landingTracks_indx4stateLDF = vertcat(landingTracks_indx4stateLDF{:});
                landingTracks_indx4stateLDF = landingTracks_indx4stateLDF(has_multiple_rrefs);
                
                for ct_track = 1:length(state_LDF_multiple_rrefs)
                    x = state_LDF_multiple_rrefs(ct_track);
                    intervals_ti = x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).intervals_ti;
                    time4rrefEstimate = arrayfun(@(i) state_LDF_multiple_rrefs(ct_track).filteredState(intervals_ti(i,1):intervals_ti(i,2), 1),1:size(intervals_ti,1), 'UniformOutput', false);
                    
                    time4rrefEstimate = time4rrefEstimate{:};
                    
                    % Find if there is any video
                    if ~isempty(videoTimes)
                        isVideoUseful = min(time4rrefEstimate)>= videoTimes(:,1) & max(time4rrefEstimate) <= videoTimes(:,2);
                    else
                        isVideoUseful = false;
                    end
                    
                    if any(isVideoUseful)
                        data4video(end+1).videoInfo = treatment.videosInfo(isVideoUseful);
                        data4video(end).state_LDF = x;
                        data4video(end).landingTrack = landingTracks(landingTracks_indx4stateLDF(ct_track));
                        data4video(end).treatment = treatment;
                        
                        if savePlots
                            track = landingTracks(landingTracks_indx4stateLDF(ct_track));
                            plotHandles = x.plot_rrefs(chosen_fac);
                            
                            if ~isempty(plotHandles)
                                
                                % Resizing the figures
                                for i=1:length(plotHandles)
                                    plotHandles(i).Position(3) = 680;
                                    plotHandles(i).Position(4) = 545;
                                    
                                    if i==1
                                        figureName = ['fac_' num2str(chosen_fac,'%0.2f') '_' ...
                                            num2str(treatment.datenum) '_' ...
                                            num2str(treatment.startTime) '_' num2str(treatment.endTime) ...
                                            '_obj' num2str(track.obj_id) '.png'];
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
%                     else
%                         data4video(end+1).videoInfo = recordedVideosInformation.empty;
                    end
                    
                end
                % 
                
            end
            
            
            
        end
    end
end

% Display track with multiple r* that have a video associated with it
for ct=1:length(data4video)
    if ~isempty(data4video(ct).videoInfo)
        track = data4video(ct).landingTrack;
        disp([num2str(track.datenum) '_' track.foldername(10:15) '_obj' num2str(track.obj_id) ]);
%         disp(ct);
    else
%         disp('empty');
    end
end    

% Find video for a particular track
for ct=1:length(data4video)
    if ~isempty(data4video(ct).videoInfo) && data4video(ct).treatment.datenum == 20190704 ...
            && data4video(ct).treatment.startTime == 110000 && data4video(ct).landingTrack.obj_id == 3481
        keyboard;
        data4video(ct).videoInfo
        track = data4video(ct).landingTrack;
        disp([num2str(track.datenum) '_' track.foldername(10:15) '_obj' num2str(track.obj_id) ]);
%         disp(ct);
    else
%         disp('empty');
    end
end  

% Selected track for video
% fac_1.00_20190704_110000_123000_obj1443_track52_excerpt1
% fac_1.00_20190704_110000_123000_obj1491_track56_excerpt1
ct_selected = 13;
track = data4video(ct_selected).landingTrack;
video = data4video(ct_selected).videoInfo;
stateLDF = data4video(ct_selected).state_LDF;

% load movie data
% [videodata, timestamps] = fmf_read(fullfile(data4video(ct_selected).videoInfo.folder, data4video(ct_selected).videoInfo.name));
[videodata, timestamps] = fmf_read(fullfile(data4video(ct_selected).videoInfo.folder, strrep(data4video(ct_selected).videoInfo.name, '84.fmf', '85.fmf')));
t_start = data4video(ct_selected).state_LDF.filteredState(1,1);
t_end = data4video(ct_selected).state_LDF.filteredState(end,1);
intervals_ti = stateLDF.rrefSegments(abs([stateLDF.rrefSegments.factor]-chosen_fac)<1e-6).intervals_ti;


% Extract and save frames between t_start and t_end
images = videodata(:,:,t_start <= timestamps & timestamps <= t_end);

% Save video
Directory = '/home/reken001/Pulkit/MATLAB/graphs';
v = VideoWriter(fullfile(Directory, [data4video(ct_selected).videoInfo.name(1:end-4) '.avi']),'Uncompressed AVI');
v.FrameRate = 30;
open(v);
for i=1:size(images,3)
    writeVideo(v,images(:,:,i));
end
close(v);

% Write videos for individual rref segments
for ct=1:size(intervals_ti,1)
    images = videodata(:,:,stateLDF.filteredState(intervals_ti(ct,1),1) <= timestamps & timestamps <= stateLDF.filteredState(intervals_ti(ct,2),1));
    v = VideoWriter(fullfile(Directory, [data4video(ct_selected).videoInfo.name(1:end-4) '_rref' num2str(ct) '.avi']),'Uncompressed AVI');
    v.FrameRate = 30;
    open(v);
    for i=1:size(images,3)
        writeVideo(v,images(:,:,i));
    end
    close(v);
    
end

[~, start_frame] = min(abs(stateLDF.filteredState(1,1)-timestamps));
[~, end_frame] = min(abs(stateLDF.filteredState(end,1)-timestamps));
images = videodata(:,:,start_frame:end_frame);
% Save images
Directory = '/home/reken001/Pulkit/MATLAB/graphs';
filename = fullfile(Directory, strrep(data4video(ct_selected).videoInfo.name, '84.fmf', '85'));
% rmdir(filename);
mkdir(filename);
for i=1:15:size(images,3)
    imwrite(images(:,:,i), fullfile(filename, [num2str(i,'%03.f') '.png']));
end

%% Loading data and collecting segments with rref
% clc; close all;
% clear;
% 
% inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A1_rref.mat';
inputFile = 'D:/light_intensity_experiments/postprocessing/BlindLandingtracks_A1_rref.mat';
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


%% Create panels a and b for Figure 3 - Trajectories' views
close all;
trajPlot = figure; hold on;
skip_step = 10;
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
radius = treatments(1).landingDiscs(1).radius;

figure(trajPlot);
% h=mesh(X,Y,Z,'facecolor',[1 0 0]); % draw landing disc
plot3([-radius; radius], [0 0], [0 0], 'LineWidth', 2, 'Color', [83 83 83]./255);
plot3([0 0], [0 0], [-radius; radius], 'LineWidth', 2, 'Color', [83 83 83]./255);
zlabel('z (m)', 'FontSize', 14);
ylabel('y (m)', 'FontSize', 14);
xlabel('x (m)', 'FontSize', 14);
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


%% Panel c for Figure 3 - Histogram rref
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
histfit(-vertcat(data.rmean),[],'Gamma')
xlabel('Estimated set-points, r* (1/s)', 'FontSize', 16);
ylabel('Occurences', 'FontSize', 16);
set(gca, 'FontSize', 16);
dummy = -vertcat(data.rmean);
[mean(dummy) median(dummy) max(dummy) min(dummy)]
pd = fitdist(dummy,'Gamma');

arrayfun(@(x) length(vertcat(data(x).rmean)), 1:length(data))
% Compute pdf for each light*pattern combination
figure;
subplot(3,2,1);
% histogram(-vertcat(data(1).rmean));
histfit(-vertcat(data(1).rmean),[],'Gamma')
legend([pattern{data(1).ct_pattern} ', ' light{data(1).ct_light}]);
set(gca, 'FontSize', 16);
ylim([0 150]);
xlim([0 8]);
subplot(3,2,3);
% histogram(-vertcat(data(2).rmean));
histfit(-vertcat(data(2).rmean),[],'Gamma')
legend([pattern{data(2).ct_pattern} ', ' light{data(2).ct_light}]);
set(gca, 'FontSize', 16);
ylim([0 150]);
xlim([0 8]);
subplot(3,2,5);
% histogram(-vertcat(data(3).rmean));
histfit(-vertcat(data(3).rmean),[],'Gamma')
legend([pattern{data(3).ct_pattern} ', ' light{data(3).ct_light}]);
set(gca, 'FontSize', 16);
ylabel('Occurences', 'FontSize', 16);
ylim([0 150]);
xlim([0 8]);


subplot(3,2,2);
% histogram(-vertcat(data(4).rmean));
histfit(-vertcat(data(4).rmean),[],'Gamma')
legend([pattern{data(4).ct_pattern} ', ' light{data(4).ct_light}]);
set(gca, 'FontSize', 16);
ylim([0 150]);
xlim([0 8]);
subplot(3,2,4);
% histogram(-vertcat(data(5).rmean));
histfit(-vertcat(data(5).rmean),[],'Gamma')
legend([pattern{data(5).ct_pattern} ', ' light{data(5).ct_light}]);
set(gca, 'FontSize', 16);
ylim([0 150]);
xlim([0 8]);
subplot(3,2,6);
% histogram(-vertcat(data(6).rmean));
histfit(-vertcat(data(6).rmean),[],'Gamma')
set(gca, 'FontSize', 16);
xlabel('Estimated set-points, r* (1/s)', 'FontSize', 16);
legend([pattern{data(6).ct_pattern} ', ' light{data(6).ct_light}]);
ylim([0 150]);
xlim([0 8]);

for ct=1:length(data)
    pdf_rref(ct) = fitdist(-vertcat(data(ct).rmean),'Gamma');
    mean_rref1(ct) = mean(pdf_rref(ct));
%     mean_rref2(ct) = mean(-vertcat(data(ct).rmean));
    median_rref1(ct) = median(pdf_rref(ct));
%     median_rref2(ct) = median(-vertcat(data(ct).rmean));
end

%% Panel d for Figure 3 (delta y1* and delta y2*)
close all; clc;

dy = cell(length(data),1); % y travelled during rrefs
dY = cell(length(data),1); % total y travelled within a track
dY_samesizeas_dy = cell(length(data),1);

% for tracks containing 2 rrefs
dy_mr = cell(length(data),1); % total y travelled during rrefs
dY_mr = cell(length(data),1); % total y travelled within a track

dt_mr_actual = cell(length(data),1); % actual time between start of first and end of second rref segment
dt_mr_analytical = cell(length(data),1); % theoretical time between start of first and end of second rref segment had the bbee continued to fly at a constant rref

rmean1_mr = cell(length(data),1); % rref in first segment (out of two)
y_mean_mr = cell(length(data),1); % mean y between the start of first segment and end of second segment

% No hover in between
has_hover = cell(length(data),1); % if there is a hover transition (V<0.05m/s) in between two rref segments
has_hover_3dspeed = cell(length(data),1); % if there is a hover transition (V<0.05m/s) in between two rref segments

% for constant tau-dot strategy
tau_dot1 = cell(length(data),1); % tau-dot at the beginning of the first r* segment
taudot = -0.812; % average in our complete dataset
Vmean_const_taudot = cell(length(data),1); % mean speed had the bbees underwent constant-taudot strategy

for ct=1:length(data)
    % Find things for rref segments (delta t, delta y)
    dt_rref = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).fd_actual_ti,data(ct).tracks_fac,'UniformOutput',false);
    
    yTravelled = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).yTravelled_ti,data(ct).tracks_fac,'UniformOutput',false);
    
    % Find things for complete appraoch tracks (delta T, delta Y)
    for ct1=1:length(data(ct).tracks_fac)
        
        [~, ymax_indx] = max(-data(ct).tracks_fac(ct1).filteredState(:,3));
        [~, ymin_indx] = min(-data(ct).tracks_fac(ct1).filteredState(:,3));
        
        dT{ct}(ct1,1) = diff(data(ct).tracks_fac(ct1).filteredState([ymax_indx, end],1));
        dY{ct}(ct1,1) = diff(data(ct).tracks_fac(ct1).filteredState([ymax_indx, ymin_indx],3));
    end
    
    assert(length(dt_rref) == length(dT{ct}));
    assert(length(yTravelled) == length(dY{ct}));
    
    dt{ct} = vertcat(dt_rref{:});
    dy{ct} = vertcat(yTravelled{:});
    
    dummy = arrayfun(@(x) dT{ct}(x)*ones(length(dt_rref{x}),1),1:length(dT{ct}),'UniformOutput', false);
    dT_samesizeas_dt{ct} = vertcat(dummy{:});
     
    dummy = arrayfun(@(x) dY{ct}(x)*ones(length(yTravelled{x}),1),1:length(dY{ct}),'UniformOutput', false);
    dY_samesizeas_dy{ct} = vertcat(dummy{:});
    
    % Only for tracks containing 2 rrefs
    has_two_rrefs = arrayfun(@(x) length(x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).rref_ti)==2,data(ct).tracks_fac);
    N = sum(has_two_rrefs);
    data_ss = data(ct).tracks_fac(has_two_rrefs);
    
    for ct1=1:length(data_ss)
        indxs = data_ss(ct1).rrefSegments(abs([data_ss(ct1).rrefSegments.factor]-chosen_fac)<1e-6).intervals_ti([1 end]);
        assert(length(indxs)==2);
        
        dy_mr{ct}(ct1,1) = diff(data_ss(ct1).filteredState(indxs,3));
        y_mean_mr{ct}(ct1,1) = -mean(data_ss(ct1).filteredState(indxs,3));
        
        dt_mr_actual{ct}(ct1,1) = diff(data_ss(ct1).filteredState(indxs,1));
        ymean = data_ss(ct1).rrefSegments(abs([data_ss(ct1).rrefSegments.factor]-chosen_fac)<1e-6).ymean_ti;
        rmean = data_ss(ct1).rrefSegments(abs([data_ss(ct1).rrefSegments.factor]-chosen_fac)<1e-6).rmean_ti;
        dummy1 = sortrows([-[ymean] -[rmean]],1,'descend');
        rmean1_mr{ct}(ct1,1) = dummy1(1,2);
        if -data_ss(ct1).filteredState(indxs(2),3) < -data_ss(ct1).filteredState(indxs(1),3)
            dy_mr{ct}(ct1,1) = diff(data_ss(ct1).filteredState(indxs,3));
            dt_mr_actual{ct}(ct1,1) = diff(data_ss(ct1).filteredState(indxs,1));
            dt_mr_analytical{ct}(ct1,1) = -(log(abs(data_ss(ct1).filteredState(indxs(2),3)/data_ss(ct1).filteredState(indxs(1),3))))./rmean1_mr{ct}(ct1,1);
            
            tau_dot{ct}(ct1,1) = (data_ss(ct1).filteredState(indxs(1),3)*data_ss(ct1).filteredState(indxs(1),9))./(data_ss(ct1).filteredState(indxs(1),6))^2-1;
            
            y = -data_ss(ct1).filteredState(indxs(2),3):0.001:-data_ss(ct1).filteredState(indxs(1),3);
            yi = -data_ss(ct1).filteredState(indxs(1),3); Vi = data_ss(ct1).filteredState(indxs(1),6);
            V_const_taudot = y.^(1+taudot)./(yi.^(1+taudot)/Vi);
%             
%             figure; hold on;
%             plot(y, rmean1_mr{ct}(ct1,1)*y,'r');
%             plot(-data_ss(ct1).filteredState(indxs(1):indxs(2),3), data_ss(ct1).filteredState(indxs(1):indxs(2),6),'g');
%             plot(y, V_const_taudot, 'b');
            
            Vmean_const_taudot{ct}(ct1,1) = mean(V_const_taudot);
            if any(data_ss(ct1).filteredState(indxs(1):indxs(2),6) < 0.05)
                has_hover{ct}(ct1,1) = true;
            else
                has_hover{ct}(ct1,1) = false;
            end
            
            if any((sum(data_ss(ct1).filteredState(indxs(1):indxs(2),5:7).^2, 2)).^0.5 < 0.15)
                has_hover_3dspeed{ct}(ct1,1) = true;
            else
                has_hover_3dspeed{ct}(ct1,1) = false;
            end
        else
            has_hover{ct}(ct1,1) = false;
            has_hover_3dspeed{ct}(ct1,1) = false;
            dy_mr{ct}(ct1,1) = nan;
            dt_mr_actual{ct}(ct1,1) = diff(data_ss(ct1).filteredState(indxs,1));
            dt_mr_analytical{ct}(ct1,1) = nan;
            Vmean_const_taudot{ct}(ct1,1) = nan;
        end
        disp([dt_mr_actual{ct}(ct1,1) dt_mr_analytical{ct}(ct1,1)])
    end
    
    dY_mr{ct} = dY{ct}(has_two_rrefs);    
end

figure;
subplot(1,2,1);
% histfit(vertcat(dy{:})./vertcat(dY_samesizeas_dy{:}),[],'Gamma')
histogram(vertcat(dy{:})./vertcat(dY_samesizeas_dy{:}), 0:0.05:1);
% histogram(vertcat(dy{:})./vertcat(dY_samesizeas_dy{:}));
xlabel('\Delta y1* / \Delta y', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);
% xlim([0 1]);

subplot(1,2,2);
% histfit(vertcat(dy_mr{:})./vertcat(dY_mr{:}),[],'Gamma')
% histogram(vertcat(dy_mr{:})./vertcat(dY_mr{:}), 0:0.05:0.8);
histogram(vertcat(dy_mr{:})./vertcat(dY_mr{:}));
xlabel('\Delta y2* / \Delta y', 'FontSize', 16);
% ylabel('Occurances', 'FontSize', 16);
% ylim([0 300])
xlim([0 1])
set(gca, 'FontSize', 16);

figure;
dt_ratio = vertcat(dt_mr_actual{:})./vertcat(dt_mr_analytical{:});
histogram(dt_ratio);
% boxplot(vertcat(dt_mr_actual{:})./vertcat(dt_mr_analytical{:}))
xlabel('\Delta t / \Delta t at rref1', 'FontSize', 16);
% ylabel('Occurances', 'FontSize', 16);
% ylim([0 300])
% xlim([0 1])
set(gca, 'FontSize', 16);
dt_norm_mean = nanmean(dt_ratio)
dt_norm_med = nanmedian(dt_ratio)

figure;
Vy2 = vertcat(dy_mr{:})./vertcat(dt_mr_actual{:});
Vy1 = vertcat(dy_mr{:})./vertcat(dt_mr_analytical{:});
histogram(Vy2./Vy1);
% boxplot(vertcat(dt_mr_actual{:})./vertcat(dt_mr_analytical{:}))
xlabel('\Delta V / \Delta V at rref1', 'FontSize', 16);
% ylabel('Occurances', 'FontSize', 16);
% ylim([0 300])
% xlim([0 1])
set(gca, 'FontSize', 16);
dV_norm_mean = nanmean(Vy2./Vy1)
dV_norm_med = nanmedian(Vy2./Vy1)

figure;
histogram((vertcat(dy_mr{:})./vertcat(dt_mr_actual{:}))./(vertcat(rmean1_mr{:}).*vertcat(y_mean_mr{:})));
xlabel('\Delta V / \Delta V at rref1 using method 2', 'FontSize', 16);
% ylabel('Occurances', 'FontSize', 16);
% ylim([0 300])
% xlim([0 1])
set(gca, 'FontSize', 16);
nanmean((vertcat(dy_mr{:})./vertcat(dt_mr_actual{:}))./(vertcat(rmean1_mr{:}).*vertcat(y_mean_mr{:})))
nanmedian((vertcat(dy_mr{:})./vertcat(dt_mr_actual{:}))./(vertcat(rmean1_mr{:}).*vertcat(y_mean_mr{:})))

figure;
has_hover = vertcat(has_hover{:});
has_hover_3dspeed = vertcat(has_hover_3dspeed{:});
sum(has_hover)/length(has_hover)
sum(has_hover_3dspeed)/length(has_hover_3dspeed)
Vratio = (vertcat(dy_mr{:})./vertcat(dt_mr_actual{:}))./(vertcat(rmean1_mr{:}).*vertcat(y_mean_mr{:}));
histogram(Vratio(~has_hover));
xlabel('\Delta V / \Delta V at rref1 for no hover in between', 'FontSize', 16);
nanmean(Vratio(~has_hover))
nanmedian(Vratio(~has_hover))

figure;
dt_ratio = vertcat(dt_mr_actual{:})./vertcat(dt_mr_analytical{:});
histogram(dt_ratio(~has_hover));
% boxplot(vertcat(dt_mr_actual{:})./vertcat(dt_mr_analytical{:}))
xlabel('\Delta t / \Delta t at rref1', 'FontSize', 16);
% ylabel('Occurances', 'FontSize', 16);
% ylim([0 300])
% xlim([0 1])
set(gca, 'FontSize', 16);
dt_norm_mean = nanmean(dt_ratio(~has_hover))
dt_norm_med = nanmedian(dt_ratio(~has_hover))

% figure;
% boxplot(vertcat(tau_dot{:}))

clc;
[mean(vertcat(dy{:})), mean(vertcat(dY{:})), mean(vertcat(dy{:}))/mean(vertcat(dY{:}))]
[std(vertcat(dy{:})), std(vertcat(dY{:})), std(vertcat(dy{:}))/mean(vertcat(dY{:}))]

[nanmean(vertcat(dy_mr{:})), nanmean(vertcat(dY_mr{:})), nanmean(vertcat(dy_mr{:}))/mean(vertcat(dY_mr{:}))]
[nanstd(vertcat(dy_mr{:})), nanstd(vertcat(dY_mr{:})), nanstd(vertcat(dy_mr{:}))/mean(vertcat(dY_mr{:}))]

mean_dist = arrayfun(@(x) [mean(vertcat(dy{x})), mean(vertcat(dY{x})), mean(vertcat(dy{x}))/mean(vertcat(dY{x}))], 1:length(dy), 'UniformOutput', false);
mean_dist = vertcat(mean_dist{:})
std_dist = arrayfun(@(x) [std(vertcat(dy{x})), std(vertcat(dY{x})), std(vertcat(dy{x}))/std(vertcat(dY{x}))], 1:length(dy), 'UniformOutput', false);
std_dist = vertcat(std_dist{:})

mean_dist_mr = arrayfun(@(x) [nanmean(vertcat(dy_mr{x})), nanmean(vertcat(dY_mr{x})), nanmean(vertcat(dy_mr{x}))/mean(vertcat(dY_mr{x}))], 1:length(dy_mr), 'UniformOutput', false);
mean_dist_mr = vertcat(mean_dist_mr{:})
std_dist_mr = arrayfun(@(x) [nanstd(vertcat(dy_mr{x})), nanstd(vertcat(dY_mr{x})), nanstd(vertcat(dy_mr{x}))/nanstd(vertcat(dY_mr{x}))], 1:length(dy_mr), 'UniformOutput', false);
std_dist_mr = vertcat(std_dist_mr{:})

for ct=1:6
%     disp([num2str(mean_dist_mr(ct,1),'%0.3f') ' [' num2str(std_dist_mr(ct,1),'%0.3f') ']']);
    disp([num2str(mean_dist(ct,1),'%0.3f') ' [' num2str(std_dist(ct,1),'%0.3f') ']']);
end

clc;

V_const_r = vertcat(rmean1_mr{:}).*vertcat(y_mean_mr{:});
V_hybrid = vertcat(dy_mr{:})./vertcat(dt_mr_actual{:});
V_const_taudot = vertcat(Vmean_const_taudot{:});

disp('Vmean+-std for constant r1* strategy');
[nanmean(V_const_r) nanstd(V_const_r)]

disp('Vmean+-std for hybrid strategy');
[nanmean(V_hybrid) nanstd(V_hybrid)]

disp('Vmean+-std for constant tau-dot strategy');
[nanmean(V_const_taudot) nanstd(V_const_taudot)]

disp('Mean+-std ratio of V_hybrid/V_V_const_r')
[nanmean(V_hybrid./V_const_r) nanstd(V_hybrid./V_const_r)]

disp('Mean+-std ratio of V_hybrid/V_const_taudot')
[nanmean(V_hybrid./V_const_taudot) nanstd(V_hybrid./V_const_taudot)]

%% Panel d for Figure 3 - delta y* / delta y histogram 
close all; clc;
dt = cell(length(data),1);
dT = cell(length(data),1); % track duration
dT_samesizeas_dt = cell(length(data),1);

dy = cell(length(data),1); % y travelled during rrefs
dY = cell(length(data),1); % total y travelled within a track
dY_samesizeas_dy = cell(length(data),1);

% for multiple rrefs in a track (_mr)
dy_mr = cell(length(data),1); % total y travelled during rrefs
dY_mr = cell(length(data),1); % total y travelled within a track

for ct=1:length(data)
    % Find things for rref segments (delta t, delta y)
    dt_rref = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).fd_actual_ti,data(ct).tracks_fac,'UniformOutput',false);
    
    yTravelled = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).yTravelled_ti,data(ct).tracks_fac,'UniformOutput',false);
    
    % Find things for complete appraoch tracks (delta T, delta Y)
    for ct1=1:length(data(ct).tracks_fac)
        
        [~, ymax_indx] = max(-data(ct).tracks_fac(ct1).filteredState(:,3));
        [~, ymin_indx] = min(-data(ct).tracks_fac(ct1).filteredState(:,3));
        
        dT{ct}(ct1,1) = diff(data(ct).tracks_fac(ct1).filteredState([ymax_indx, end],1));
        dY{ct}(ct1,1) = diff(data(ct).tracks_fac(ct1).filteredState([ymax_indx, ymin_indx],3));
    end
    
    assert(length(dt_rref) == length(dT{ct}));
    assert(length(yTravelled) == length(dY{ct}));
    
    dt{ct} = vertcat(dt_rref{:});
    dy{ct} = vertcat(yTravelled{:});
    
    dummy = arrayfun(@(x) dT{ct}(x)*ones(length(dt_rref{x}),1),1:length(dT{ct}),'UniformOutput', false);
    dT_samesizeas_dt{ct} = vertcat(dummy{:});
     
    dummy = arrayfun(@(x) dY{ct}(x)*ones(length(yTravelled{x}),1),1:length(dY{ct}),'UniformOutput', false);
    dY_samesizeas_dy{ct} = vertcat(dummy{:});
    
    % Only for tracks containing multiple rrefs
    has_multiple_rrefs = arrayfun(@(x) length(x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).rref_ti)>1,data(ct).tracks_fac);
    N = sum(has_multiple_rrefs);
    data_ss = data(ct).tracks_fac(has_multiple_rrefs);
    
    dy_mr{ct} = arrayfun(@(x) sum(x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).yTravelled_ti),data_ss)';
    dY_mr{ct} = dY{ct}(has_multiple_rrefs);
    
end

figure;
histogram(vertcat(dt{:}))
xlabel('\Delta t (s)', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);

figure;
histogram(vertcat(dT{:}))
xlabel('\Delta T (s)', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);

figure;
histogram(vertcat(dt{:})./vertcat(dT_samesizeas_dt{:}))
xlabel('\Delta t / \Delta T', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);


figure;
% histogram(vertcat(dy{:}))
histfit(vertcat(dy{:}),[],'Gamma');
xlabel('\Delta y (m)', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);

figure;
histogram(vertcat(dY{:}))
xlabel('\Delta Y (m)', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);

figure;
subplot(2,1,1);
% histfit(vertcat(dy{:})./vertcat(dY_samesizeas_dy{:}),[],'Gamma')
ax1 = histogram(vertcat(dy{:})./vertcat(dY_samesizeas_dy{:}), 0:0.05:0.8);
% xlabel('\Delta y / \Delta Y', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);

subplot(2,1,2);
% histfit(vertcat(dy_mr{:})./vertcat(dY_mr{:}),[],'Gamma')
ax2 = histogram(vertcat(dy_mr{:})./vertcat(dY_mr{:}), 0:0.05:0.8);
xlabel('\Delta y / \Delta Y', 'FontSize', 16);
ylabel('Occurances', 'FontSize', 16);
ylim([0 300])
set(gca, 'FontSize', 16);

clc;

[mean(vertcat(dt{:})), mean(vertcat(dT{:})), mean(vertcat(dt{:}))/mean(vertcat(dT{:}))]

[mean(vertcat(dy{:})), mean(vertcat(dY{:})), mean(vertcat(dy{:}))/mean(vertcat(dY{:}))]
[std(vertcat(dy{:})), std(vertcat(dY{:})), std(vertcat(dy{:}))/mean(vertcat(dY{:}))]

[mean(vertcat(dy_mr{:})), mean(vertcat(dY_mr{:})), mean(vertcat(dy_mr{:}))/mean(vertcat(dY_mr{:}))]
[std(vertcat(dy_mr{:})), std(vertcat(dY_mr{:})), std(vertcat(dy_mr{:}))/mean(vertcat(dY_mr{:}))]
%% Panel d for Figure 4 (ymean vs rref) for each light*pattern combination
% plotted only for tracks containing >1 r* segments
close all; clc;
saveDir = '/home/reken001/Pulkit/graphs_temp';
red_cmap = [252,187,161
252,146,114
251,106,74
239,59,44
203,24,29
165,15,21
103,0,13]./255;
cmap = [60 180 75; 245 130 48; 0 130 200]/255; % [green orange blue]


delta_rrref_posjumps = cell(length(data),1); % change in rref between two consecutive rref segments
delta_rrref_negjumps = cell(length(data),1); % change in rref between two consecutive rref segments

ratio_rrref_posjumps = cell(length(data),1); % change in rref between two consecutive rref segments
ratio_rrref_negjumps = cell(length(data),1); % change in rref between two consecutive rref segments

data_ymean = {}; data_rmean = {};
for ct=1:length(data)
    disp([pattern{data(ct).ct_pattern} '_' light{data(ct).ct_light}]);
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
% % %         plot(-[ymean{ct1}], -[rmean{ct1}] ...
% % %         ,'Color', [225 225 225]./255/2,'LineWidth',1);
        if rem(ct1,10) == 1
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
    legend([pattern{data(ct).ct_pattern} ' ' light{data(ct).ct_light}]);
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
for ct=1:6
    disp([num2str(dummy(ct,1),'%0.3f') ' [' num2str(dummy(ct,2),'%0.3f') ']']);
end
dummy =cellfun(@(x) [mean(x) std(x)],delta_rrref_negjumps,'UniformOutput',false);
vertcat(dummy{:})
dummy =cellfun(@(x) [mean(x) std(x)],ratio_rrref_posjumps,'UniformOutput',false);
vertcat(dummy{:})
dummy =cellfun(@(x) [mean(x) std(x)],ratio_rrref_negjumps,'UniformOutput',false);
vertcat(dummy{:})
%% Figure 5 (ymean vs rmean for each pattern) average over all lights
close all;
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
        dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).rmean_ti,data_ss(ct_light).tracks_fac,'UniformOutput',false);
        rmean{ct} = [rmean{ct}; vertcat(dummy{:})];
        
        dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).ymean_ti,data_ss(ct_light).tracks_fac,'UniformOutput',false);
        ymean{ct} = [ymean{ct}; vertcat(dummy{:})];
    end
end
data2plot = [vertcat(ymean{:}) vertcat(rmean{:}) [ones(size(rmean{1})); 2*ones(size(rmean{2}))]];
data2plot = data2plot(randperm(size(data2plot, 1)), :);
figure; hold on;
% plot lines
y_vec = min(-ymean{1}):0.001:max(-ymean{1});
modelfun = @(b,x)(exp(b(1))*x.^(b(2)));
for ct=1:length(pattern)
    Coefficients = [mean(intercepts(ct,:)) mean(taudots(ct,:))];
    plot(y_vec,modelfun(Coefficients,y_vec),'Color', line_cmap(ct,:),'LineWidth',3);
end
legend({'checkerboard','spoke'})
set(gca, 'FontSize', 16);
ylabel('Estimated set-points, r* (1/s)', 'FontSize', 16);
xlabel('y* (m)', 'FontSize', 16);

scatter(-data2plot(:,1),-data2plot(:,2),10,points_cmap(round(data2plot(:,3)),:),'filled','o');
for ct=1:length(pattern)
    Coefficients = [mean(intercepts(ct,:)) mean(taudots(ct,:))];
    plot(y_vec,modelfun(Coefficients,y_vec),'Color', line_cmap(ct,:),'LineWidth',3);
end

% log plots
modelfun = @(b,x)(b(1)+b(2).*x);
figure; hold on;
for ct=1:length(pattern)
    Coefficients = [mean(intercepts(ct,:)) mean(taudots(ct,:))];
    plot(log(y_vec),modelfun(Coefficients,log(y_vec)),'Color', line_cmap(ct,:),'LineWidth',3);
end
legend({'checkerboard','spoke'})
%% Figure 5 (ymean vs rmean for each pattern) and one light condition (medium)
% Panel a
close all;
chosen_light = 2; % medium

taudots = [-0.8443719 -0.8443719 -0.8443719;
           -0.7799527 -0.7799527 -0.7799527];
% taudot_checker = -0.847; 
% taudot = -0.780; % from analysis in R

intercepts = [-1.117778 -1.045172 -0.9543105;
              -0.9693903 -0.8967847 -0.8059232]; 
          
points_cmap = [252,187,161;
200,200,200]./255;

line_cmap = [215,48,31;
37,37,37]./255;

rmean = cell(1,length(pattern));
ymean = cell(1,length(pattern));
for ct=1:length(pattern)
    data_ss = vertcat(data([data.ct_pattern]==ct & [data.ct_light]==chosen_light));
    
    rmean{ct} = [];
    ymean{ct} = [];
    for ct_light=1:length(data_ss)
        dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).rmean_ti,data_ss(ct_light).tracks_fac,'UniformOutput',false);
        rmean{ct} = [rmean{ct}; vertcat(dummy{:})];
        
        dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).ymean_ti,data_ss(ct_light).tracks_fac,'UniformOutput',false);
        ymean{ct} = [ymean{ct}; vertcat(dummy{:})];
    end
end
data2plot = [vertcat(ymean{:}) vertcat(rmean{:}) [ones(size(rmean{1})); 2*ones(size(rmean{2}))]];
data2plot = data2plot(randperm(size(data2plot, 1)), :);
figure; hold on;
% plot lines
y_vec = min(-ymean{1}):0.001:max(-ymean{1});
modelfun = @(b,x)(exp(b(1))*x.^(b(2)));
for ct=1:length(pattern)
    Coefficients = [mean(intercepts(ct,chosen_light)) mean(taudots(ct,chosen_light))];
    plot(y_vec,modelfun(Coefficients,y_vec),'Color', line_cmap(ct,:),'LineWidth',3);
end
legend({'checkerboard','spoke'})
set(gca, 'FontSize', 16);
ylabel('Estimated set-points, r* (1/s)', 'FontSize', 16);
xlabel('y* (m)', 'FontSize', 16);

scatter(-data2plot(:,1),-data2plot(:,2),10,points_cmap(round(data2plot(:,3)),:),'filled','o');
ylim([0 10])


% Panel b - log plots
modelfun = @(b,x)(b(1)+b(2).*x);
figure; hold on;
for ct=1:length(pattern)
    Coefficients = [mean(intercepts(ct,:)) mean(taudots(ct,:))];
    plot(log(y_vec),modelfun(Coefficients,log(y_vec)),'Color', line_cmap(ct,:),'LineWidth',3);
end
legend({'checkerboard','spoke'})
set(gca, 'FontSize', 16);
ylabel('log(r*) (1/s)', 'FontSize', 16);
xlabel('log(y*) (m)', 'FontSize', 16);

scatter(log(-data2plot(:,1)),log(-data2plot(:,2)),10,points_cmap(round(data2plot(:,3)),:),'filled','o');

% Panel c
% from R
means = [0.580 0.637 0.728 0.575 0.663 0.755]; % mean r* at ymean
se = [0.0345 0.0332 0.0319 0.0341 0.0323 0.0313];
ymean = exp(-1.996058); % mean(data$logy) in R
ylim([-4 3])


%% Figure supplement (ymean vs vmean for each pattern) and one light condition
% Panel a
close all;
chosen_light = 3;

taudots = [0.1554519 0.1554519 0.1554519;
           0.2198993 0.2198993 0.2198993] - 1;
% taudot_checker = -0.847; 
% taudot = -0.780; % from analysis in R

intercepts = [-1.118053 -1.045507 -0.954663;
              -0.969576 -0.8970297 -0.806186]; 
          
points_cmap = [252,187,161;
200,200,200]./255;

line_cmap = [215,48,31;
37,37,37]./255;

vmean = cell(1,length(pattern));
ymean = cell(1,length(pattern));
for ct=1:length(pattern)
    data_ss = vertcat(data([data.ct_pattern]==ct & [data.ct_light]==chosen_light));
    
    vmean{ct} = [];
    ymean{ct} = [];
    for ct_light=1:length(data_ss)
        dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).vmean_ti,data_ss(ct_light).tracks_fac,'UniformOutput',false);
        vmean{ct} = [vmean{ct}; -vertcat(dummy{:})];
        
        dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).ymean_ti,data_ss(ct_light).tracks_fac,'UniformOutput',false);
        ymean{ct} = [ymean{ct}; vertcat(dummy{:})];
    end
end
data2plot = [vertcat(ymean{:}) vertcat(vmean{:}) [ones(size(vmean{1})); 2*ones(size(vmean{2}))]];
data2plot = data2plot(randperm(size(data2plot, 1)), :);
figure; hold on;

% plot lines
modelfun = @(b,x)(exp(b(1))*x.^(b(2)));
for ct=1:length(pattern)
    y_vec = min(-ymean{ct}):0.001:max(-ymean{ct});
    Coefficients = [mean(intercepts(ct,chosen_light)) 1+mean(taudots(ct,chosen_light))];
    plot(y_vec,modelfun(Coefficients,y_vec),'Color', line_cmap(ct,:),'LineWidth',3);
    
    y_vec = 0:0.001:min(-ymean{ct});
    plot(y_vec,modelfun(Coefficients,y_vec), '-.', 'Color', line_cmap(ct,:),'LineWidth',3);

end
legend({'checkerboard',' ','spoke'})
set(gca, 'FontSize', 16);
ylabel('V* (m/s)', 'FontSize', 16);
xlabel('y* (m)', 'FontSize', 16);
ylim([0 0.8]);
scatter(-data2plot(:,1),-data2plot(:,2),10,points_cmap(round(data2plot(:,3)),:),'filled','o');


% Panel b - log plots
modelfun = @(b,x)(b(1)+(1+b(2)).*x);
figure; hold on;
for ct=1:length(pattern)
    y_vec = min(-ymean{ct}):0.001:max(-ymean{ct});
    Coefficients = [mean(intercepts(ct,:)) mean(taudots(ct,:))];
    plot(log(y_vec),modelfun(Coefficients,log(y_vec)),'Color', line_cmap(ct,:),'LineWidth',3);
end
legend({'checkerboard','spoke'})
set(gca, 'FontSize', 16);
ylabel('log(V*) (m/s)', 'FontSize', 16);
xlabel('log(y*) (m)', 'FontSize', 16);
ylim([-5 0]);
scatter(log(-data2plot(:,1)),log(-data2plot(:,2)),10,points_cmap(round(data2plot(:,3)),:),'filled','o');

% Panel c
% from R
means = [0.580 0.637 0.728 0.575 0.663 0.755]; % mean r* at ymean
se = [0.0345 0.0332 0.0319 0.0341 0.0323 0.0313];
ymean = exp(-1.996058); % mean(data$logy) in R

%% Figure Supplement (ymean vs rmean for each pattern) average over all lights
% Only with tracks containing >1 r* segments

close all;
% from analysis in R
taudots = [-0.7834318 -0.7834318 -0.7834318;
           -0.7101198 -0.7101198 -0.7101198];

intercepts = [-1.047322 -0.9905023 -0.8825408;
              -0.8899269 -0.8331076 -0.7251461];  

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
        has_multiple_rrefs = arrayfun(@(x) length(x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).rref_ti)>1,data_ss(ct_light).tracks_fac);
        
        dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).rmean_ti,[data_ss(ct_light).tracks_fac(has_multiple_rrefs)],'UniformOutput',false);
        rmean{ct} = [rmean{ct}; vertcat(dummy{:})];
        
        dummy = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-chosen_fac)<1e-6).ymean_ti,[data_ss(ct_light).tracks_fac(has_multiple_rrefs)],'UniformOutput',false);
        ymean{ct} = [ymean{ct}; vertcat(dummy{:})];
    end
end
data2plot = [vertcat(ymean{:}) vertcat(rmean{:}) [ones(size(rmean{1})); 2*ones(size(rmean{2}))]];
data2plot = data2plot(randperm(size(data2plot, 1)), :);
figure; hold on;
% plot lines
y_vec = min(-ymean{1}):0.001:max(-ymean{1});
modelfun = @(b,x)(exp(b(1))*x.^(b(2)));
for ct=1:length(pattern)
    Coefficients = [mean(intercepts(ct,:)) mean(taudots(ct,:))];
    plot(y_vec,modelfun(Coefficients,y_vec),'Color', line_cmap(ct,:),'LineWidth',3);
end
legend({'checkerboard','spoke'})
set(gca, 'FontSize', 16);
ylabel('Estimated set-points, r* (1/s)', 'FontSize', 16);
xlabel('y* (m)', 'FontSize', 16);

scatter(-data2plot(:,1),-data2plot(:,2),10,points_cmap(round(data2plot(:,3)),:),'filled','o');
for ct=1:length(pattern)
    Coefficients = [mean(intercepts(ct,:)) mean(taudots(ct,:))];
    plot(y_vec,modelfun(Coefficients,y_vec),'Color', line_cmap(ct,:),'LineWidth',3);
end

% log plots
modelfun = @(b,x)(b(1)+b(2).*x);
figure; hold on;
for ct=1:length(pattern)
    Coefficients = [mean(intercepts(ct,:)) mean(taudots(ct,:))];
    plot(log(y_vec),modelfun(Coefficients,log(y_vec)),'Color', line_cmap(ct,:),'LineWidth',3);
end
legend({'checkerboard','spoke'})

% Panel b - log plots
modelfun = @(b,x)(b(1)+(1+b(2)).*x);
figure; hold on;
for ct=1:length(pattern)
    y_vec = min(-ymean{ct}):0.001:max(-ymean{ct});
    Coefficients = [mean(intercepts(ct,:)) mean(taudots(ct,:))];
    plot(log(y_vec),modelfun(Coefficients,log(y_vec)),'Color', line_cmap(ct,:),'LineWidth',3);
end
legend({'checkerboard','spoke'})
set(gca, 'FontSize', 16);
ylabel('log(V*) (m/s)', 'FontSize', 16);
xlabel('log(y*) (m)', 'FontSize', 16);

%% ymean vs rmean for each light*pattern combination
close all;
red_cmap = [252,187,161
252,146,114
251,106,74
239,59,44
203,24,29
165,15,21
103,0,13]./255;
for ct=1:length(data)
    figure;
    f = scatter(-data(ct).ymean,-data(ct).rref,40,red_cmap(1,:),'filled','s');
    dummy = [pattern(data(ct).ct_pattern ) ', ' light(data(ct).ct_light) ];
    legend([dummy{:}]);
    set(gca, 'FontSize', 16);
    ylabel('Estimated set-points, r* (1/s)', 'FontSize', 16);
    xlabel('Mean y (m)', 'FontSize', 16);
    
end

%% Checking PERCENTAGE of +ve and -ve jumps for all factors
clear dummy;
treatments = treatments(1:14*8); % Taking experiments for 2 patterns * 3 lights


pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
behaviour = {'rising','constant','sleeping'};

factors = [0.25:0.25:2.5];
chosen_fac = 1;

N_posJumps = cell(length(factors),1);
N_negJumps = cell(length(factors),1);
% tracks_fac(length(pattern), length(light)) = filteredState_BlindLandingtrack.empty; % tracks for chosen factor

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
            
            % For chosen factor, collect all tracks that do contain
            % non-empty rref segments
            for ct_fac = 1:length(factors)
                factor = factors(ct_fac);
                indices = arrayfun(@(x) ~isempty(x.rrefSegments(abs([x.rrefSegments.factor]-factor)<1e-6).intervals_ti), state_LDF);
                dummy.tracks_fac = state_LDF(indices);
                dummy.ct_pattern = ct_pattern;
                dummy.ct_light = ct_light;
                
                has_multiple_rrefs = arrayfun(@(x) length(x.rrefSegments(abs([x.rrefSegments.factor]-factor)<1e-6).rref_ti)>1,dummy.tracks_fac);
                N = sum(has_multiple_rrefs);
                data_ss = dummy.tracks_fac(has_multiple_rrefs);
                
                rmean = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-factor)<1e-6).rmean_ti,data_ss,'UniformOutput',false);
                
                ymean = arrayfun(@(x) x.rrefSegments(abs([x.rrefSegments.factor]-factor)<1e-6).ymean_ti,data_ss,'UniformOutput',false);
                
                jumps = cell(size(ymean));
                disp(['# of tracks with >1 r* segments for ' num2str(factor) ' factor: ' num2str(N)]);
                for ct1=1:N
                    % Finding whether the next r* is higher or lower
                    jumps{ct1} = diff(sortrows([-[ymean{ct1}] -[rmean{ct1}]],1,'descend'));
                end
                disp(['# r*segments in tracks with >1 r* segments for ' num2str(factor) ' fac: ' num2str(length(vertcat(ymean{:})))]);
                
                dummy1 = vertcat(jumps{:});
                disp(['# of +ve jumps in r* for chosen fac: ' num2str(sum(dummy1(:,2)>0))]);
                disp(['# of -ve jumps in r* for chosen fac: ' num2str(sum(dummy1(:,2)<0))]);
                
                disp(['% of # of +ve jumps in r* for chosen fac: ' num2str(sum(dummy1(:,2)>0)/(sum(dummy1(:,2)<0)+sum(dummy1(:,2)>0)))]);
                disp(['% of # of -ve jumps in r* for chosen fac: ' num2str(sum(dummy1(:,2)<0)/(sum(dummy1(:,2)<0)+sum(dummy1(:,2)>0)))]);
                
                N_posJumps{ct_fac}(ct_pattern, ct_light) = sum(dummy1(:,2)>0);
                N_negJumps{ct_fac}(ct_pattern, ct_light) = sum(dummy1(:,2)<0);
            end
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end
figure; hold on;
percent_pos = arrayfun(@(i) sum(sum(N_posJumps{i}))./sum(sum(N_posJumps{i}+N_negJumps{i})),1:length(factors));
percent_neg = arrayfun(@(i) sum(sum(N_negJumps{i}))./sum(sum(N_posJumps{i}+N_negJumps{i})),1:length(factors));
plot(factors,percent_pos*100);
plot(factors,percent_neg*100);
legend('% positive jumps', '% negative jumps');
xlabel('factors');
ylabel('percentage');

%% Plot histogram of y_start and V_start for landing approaches
% data is only for chosen factor
y_start = cell(length(data),1);
V_start = cell(length(data),1);
for ct=1:length(data)
    for ct1=1:length(data(ct).tracks_fac)
        [y_start{ct}(ct1,1), ymin_indx] = min(data(ct).tracks_fac(ct1).filteredState(:,3));
        V_start{ct}(ct1,1) = data(ct).tracks_fac(ct1).filteredState(ymin_indx,6);
    end
%     figure;
%     histogram(-vertcat(y_start{ct}));
end
% figure;
% histogram(-vertcat(y_start{1:3}));
% figure;
% histogram(-vertcat(y_start{3:end}));

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
hist3([-vertcat(y_start{:}) vertcat(V_start{:})]);
xlabel('y_0 (m)', 'FontSize', 16);
ylabel('V_0 (m/s)', 'FontSize', 16);
zlabel('Occurances', 'FontSize', 16);
set(gca, 'FontSize', 16);
view([-146 19])

[max(-vertcat(y_start{:})) min(-vertcat(y_start{:}))]
[max(vertcat(V_start{:})) min(vertcat(V_start{:}))]

%% Make curve of V vs y for comaprison with Chang et. al. (2016)
close all;

% Plot all individual tracks
figure; hold on;
for ct=1:length(data)
    for ct1=1:length(data(ct).tracks_fac)
        % y is pointing inside the platform here, therefore looking at ymin
        % for maximum distance away from the platform
        [~, ymin_indx] = min(data(ct).tracks_fac(ct1).filteredState(:,3));
        
        % Plot the approach between max distance to min distance
        plot(-data(ct).tracks_fac(ct1).filteredState(ymin_indx:end,3), data(ct).tracks_fac(ct1).filteredState(ymin_indx:end,6),'Color',[189,189,189]/255)
%         plot(data(ct).tracks_fac(ct1).filteredState(:,3), data(ct).tracks_fac(ct1).filteredState(:,6),'Color',[189,189,189]/255)
    end
end
% xlim([-0.4 0])
% ylim([0 0.5])

figure; hold on;
% Plot mean track
ybins = -0.4:0.001:-0.02;
tracks_fac = [data.tracks_fac];
data_yV = arrayfun(@(x) x.filteredState(x.filteredState(:,3)>=ybins(1) & x.filteredState(:,3)<=ybins(end),[3 6]), tracks_fac, 'UniformOutput', false);
data_yV = vertcat(data_yV{:});
y = [];
V = [];
semV = []; % standard error of the means, SEM = std(data)/sqrt(length(data));
for ct=1:length(ybins)-1
    y = [y mean(ybins(ct:ct+1))];
    
    dummy = data_yV(data_yV(:,1)>=ybins(ct) & data_yV(:,1)<ybins(ct+1),2);
    V = [V mean(dummy)];
    semV = [semV std(dummy)/sqrt(length(dummy))];
end

plot(-y,V+semV,'--k');
plot(-y,V-semV,'--k');
fill([-y fliplr(-y)],[V+semV fliplr(V-semV)], 'g');
plot(-y,V,'k','Linewidth',1)
set(gca, 'FontSize', 16);
ylabel('V (m/s)', 'FontSize', 16);
xlabel('y (m)', 'FontSize', 16);
xlim([0 0.4])
ylim([0 0.3])
yticks([0:0.1:0.3]);


figure;
plot(-y, V./-y, 'k')

% % Plot V and y 3d histogram
% figure;
% hist3([-data_yV(:,1) data_yV(:,2)]);
% xlabel('y (m)', 'FontSize', 16);
% ylabel('V (m/s)', 'FontSize', 16);
% zlabel('Occurances', 'FontSize', 16);
% set(gca, 'FontSize', 16);
% view([-146 19])
%% Plot mean track for each pattern*light combination
close all; clc;
cmap = [107,174,214; 49,130,189; 8,81,156;
        116,196,118; 49,163,84; 0,109,44]./255;
cmap = [223,194,125; 191,129,45; 140,81,10;
        128,205,193; 53,151,143; 1,102,94]./255;
    
line_cmap = [252,141,89; 239,101,72; 215,48,31;
             189,189,189; 150,150,150; 115,115,115]./255;  
         
% using model with random factors (lmer)        
intercepts = [0.004601463, -0.006137268, -0.009488072;
              0.03938711, 0.01846472,   0.01432766];
slopes = [2.11 2.26 2.69; 1.47 1.99 2.28];

% simple linear regression (lm)
% intercepts = [-0.008710276, -0.01219497, -0.01865054;
%               0.02511735, 0.008219641,   0.003003411];
% slopes = [1.86 1.97 2.36; 1.26 1.71 1.97];

yrange = 0.002:0.001:0.12;
func = (@(b,x) b(1)+b(2)*x);

ybins = -0.4:0.001:-0.02;
fig1 = figure;
subplot(2,1,1); hold on;
subplot(2,1,2); hold on;

fig2 = figure;
subplot(2,1,1); hold on;
subplot(2,1,2); hold on;

yV = cell(length(data),1);
for j=1:length(data)
    tracks_fac = data(j).tracks_fac;
    data_yV = arrayfun(@(x) x.filteredState(x.filteredState(:,3)>=ybins(1) & x.filteredState(:,3)<=ybins(end),[3 6]), tracks_fac, 'UniformOutput', false);
    data_yV = vertcat(data_yV{:});
    y = [];
    V = [];
    for ct=1:length(ybins)-1
        y = [y mean(ybins(ct:ct+1))];
        V = [V mean(data_yV(data_yV(:,1)>=ybins(ct) & data_yV(:,1)<ybins(ct+1),2))];
    end
    
    figure(fig1)
    subplot(2,1,1)
    plot(-y,V,'Color',cmap(j,:), 'Linewidth', 2);
    
    subplot(2,1,2)
    plot(-y,V./-y,'Color',cmap(j,:), 'Linewidth', 2);

    Legend{j} = [pattern{data(j).ct_pattern} ', ' light{data(j).ct_light}];
    
    figure(fig2)
    subplot(2,1,data(j).ct_pattern);
    plot(yrange, func([intercepts(data(j).ct_pattern, data(j).ct_light) slopes(data(j).ct_pattern, data(j).ct_light)], yrange), 'Color', line_cmap(j,:), 'Linewidth', 2);
    plot(-y,V,'Color',cmap(j,:), 'Linewidth', 2);
end
figure(fig1);
subplot(2,1,1);
xlim([0 0.4])
ylim([0 0.3])
yticks([0:0.1:0.3]);
set(gca, 'FontSize', 16);
ylabel('V (m/s)', 'FontSize', 16);

subplot(2,1,2);
set(gca, 'FontSize', 16);
ylabel('r (1/s)', 'FontSize', 16);
xlabel('y (m)', 'FontSize', 16);
legend(Legend);




%% Make curve of V vs y for approaches starting with high speeds
close all;

% Plot all individual tracks
figure; hold on;
for ct=1:length(data)
    for ct1=1:length(data(ct).tracks_fac)
        % y is pointing inside the platform here, therefore looking at ymin
        % for maximum distance away from the platform
        [y_min, ymin_indx] = min(data(ct).tracks_fac(ct1).filteredState(:,3));
        V_min = data(ct).tracks_fac(ct1).filteredState(ymin_indx,6);
        
        if -y_min > 0.35 && V_min > 0.5
            % Plot the approach 
            plot(-data(ct).tracks_fac(ct1).filteredState(ymin_indx:end,3), data(ct).tracks_fac(ct1).filteredState(ymin_indx:end,6),'Color',[189,189,189]/255)
        end
    end
end

%% Independence with factor f
close all;
stat_coeffs = [ 0.25000000 -1.08846514 -0.86200270  0.10949501  0.06657793  0.19319169  0.05412709
                0.50000000 -1.05869342 -0.83253711  0.12561129  0.06431783  0.14973640  0.04629612
                0.75000000 -1.11607669 -0.84773447  0.14519103  0.07566884  0.15934225  0.05720759
                1.00000000 -1.11777760 -0.84437193  0.14838726  0.07260567  0.16346712  0.06441925
                1.25000000 -1.12105049 -0.84082487  0.17032351  0.08581880  0.16437322  0.07830893
                1.50000000 -1.09069319 -0.81935199  0.15186821  0.09522983  0.16607479  0.06834520
                1.75000000 -1.08058143 -0.81209318  0.12979338  0.09425142  0.16696699  0.06089656
                2.00000000 -1.07886332 -0.80658529  0.14057835  0.09354859  0.17341495  0.06461033
                2.25000000 -1.06817079 -0.79794962  0.13338053  0.08850605  0.17219498  0.06026348
                2.50000000 -1.08102648 -0.80145288  0.11972203  0.09113614  0.17120298  0.05549719]; 
                %[f     (Intercept)          log(y)        pattern2          light2          light3 log(y):pattern2

                
se=[  0.25 0.09780280 0.04621912 0.12484919 0.03397752 0.03272497 0.06050423
  0.50 0.06284643 0.02686628 0.08357760 0.02152004 0.02065386 0.03639810
  0.75 0.05245319 0.02126458 0.06672413 0.01810199 0.01724768 0.02826893
  1.00 0.04753647 0.01790679 0.06081615 0.01598175 0.01520607 0.02386188
  1.25 0.04453820 0.01580678 0.05643285 0.01479127 0.01407374 0.02100452
  1.50 0.04293834 0.01456209 0.05242353 0.01408837 0.01338564 0.01922425
  1.75 0.04275667 0.01374808 0.05140469 0.01341242 0.01271499 0.01806045
  2.00 0.03999263 0.01284847 0.04922143 0.01285597 0.01217027 0.01692365
  2.25 0.03938018 0.01213523 0.04945570 0.01242411 0.01176058 0.01602704
 2.50 0.03935380 0.01170116 0.04832551 0.01216373 0.01149490 0.01540906];

figure;
tiledlayout(6,1);

nexttile
errorbar(stat_coeffs(:,1),stat_coeffs(:,2), 2*se(:,2),'Linewidth',2)
ylim([-1.3 -0.8]);
yticks([-1.3:0.25:-0.8]);
ylabel('Intercept', 'FontSize', 16);
set(gca, 'FontSize', 16);

nexttile
errorbar(stat_coeffs(:,1),stat_coeffs(:,3), 2*se(:,3),'Linewidth',2)
ylim([-1 -0.7]);
yticks([-1:0.15:-0.7]);
ylabel('log(y)', 'FontSize', 16);
set(gca, 'FontSize', 16);

nexttile
errorbar(stat_coeffs(:,1),stat_coeffs(:,4), 2*se(:,4),'Linewidth',2)
ylim([-0.2 0.4]);
yticks([-0.2:0.3:0.4]);
ylabel('pattern2', 'FontSize', 16);
set(gca, 'FontSize', 16);

nexttile
errorbar(stat_coeffs(:,1),stat_coeffs(:,5), 2*se(:,5),'Linewidth',2)
ylim([-0.1 0.2]);
yticks([-0.1:0.15:0.2]);
ylabel('light2', 'FontSize', 16);
set(gca, 'FontSize', 16);

nexttile
errorbar(stat_coeffs(:,1),stat_coeffs(:,6), 2*se(:,6),'Linewidth',2)
ylim([0.1 0.3]);
yticks([0.1 0.2 0.3]);
ylabel('light3', 'FontSize', 16);
set(gca, 'FontSize', 16);

nexttile
errorbar(stat_coeffs(:,1),stat_coeffs(:,7), 2*se(:,7),'Linewidth',2)
ylim([-0.1 0.2]);
yticks([-0.1:0.15:0.2]);
ylabel('log(y):pattern2', 'FontSize', 16);
set(gca, 'FontSize', 16);

xlabel('factor f', 'FontSize', 16);

% %%%% Load data_write corresponding to data_all_rref_Rstudio.txt file.
% data_write can be loaded from AD_extract_rref_intervals.m
factors = unique(data_write(:,10));
clear pt;
for ct=1:length(factors)
    factor=factors(ct);
    dummy = data_write(abs(data_write(:,10)-factor)<1e-6, 8);
    pt(ct) = fitdist(dummy,'Gamma');
    median_rref(ct) = median(dummy);
    percentiles(ct,:) = prctile(dummy,[25, 50,75]);
    mean_rref(ct) = mean(dummy);
end
dummy = arrayfun(@(x) paramci(x), pt, 'UniformOutput', false);
ci = vertcat(dummy{:});
a = [pt.a]; b = [pt.b];
negpos_a = arrayfun(@(x) abs(paramci(pt(x),'Parameter','a')'-pt(x).a),1:length(pt), 'UniformOutput', false); negpos_a = vertcat(negpos_a{:});
negpos_b = arrayfun(@(x) abs(paramci(pt(x),'Parameter','b')'-pt(x).b),1:length(pt), 'UniformOutput', false); negpos_b = vertcat(negpos_b{:});
median_pt = arrayfun(@(x) median(x), pt); mean_pt = arrayfun(@(x) mean(x), pt);

figure;
tiledlayout(2,1);

nexttile
errorbar(factors,a,negpos_a(:,1),negpos_a(:,2) ,'Linewidth',2)
% ylim([-1.3 -0.8]);
% yticks([-1.3:0.25:-0.8]);
ylabel('a', 'FontSize', 16);
set(gca, 'FontSize', 16);

nexttile
errorbar(factors,b, negpos_b(:,1), negpos_b(:,2),'Linewidth',2)
% ylim([-1 -0.7]);
% yticks([-1:0.15:-0.7]);
ylabel('b', 'FontSize', 16);
set(gca, 'FontSize', 16);

figure;
tiledlayout(2,1);

nexttile
plot(factors,mean_rref ,'Linewidth',2)
ylim([0 4]);
% yticks([-1.3:0.25:-0.8]);
ylabel('mean r*', 'FontSize', 16);
set(gca, 'FontSize', 16);

nexttile
errorbar(factors,median_rref, diff(percentiles(:,1:2),[],2), diff(percentiles(:,2:3),[],2),'Linewidth',2)
ylim([0 4]);
% yticks([-1:0.15:-0.7]);
ylabel('median r*', 'FontSize', 16);
set(gca, 'FontSize', 16);
xlabel('factor f', 'FontSize', 16);
%% Displaying information about the # of tracks in each treatment
% This file contains data for all the factors
clc;
factors = unique([data.factor]);
for ct_factor=1:length(factors)
    factor = factors(ct_factor);
    
    data_fac = data(abs([data.factor]-factor)<1e-6)';
    N = length(data_fac);
    
    % Create nominal and ordinal variables
    approach = arrayfun(@(i) i*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    side = arrayfun(@(i) data_fac(i).side*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    pattern = arrayfun(@(i) data_fac(i).pattern*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    light = arrayfun(@(i) data_fac(i).light*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    time = arrayfun(@(i) data_fac(i).time*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    day = arrayfun(@(i) data_fac(i).day*ones(size(data_fac(i).intervals_ti,1),1),1:N,'UniformOutput',false);
    
    % create other variables
    logy = log(-vertcat(data_fac.ymean_ti));
    logr = log(-vertcat(data_fac.rref_ti));
    
    data_write = [data_write; ...
        vertcat(approach{:}) vertcat(side{:}) vertcat(pattern{:}) ...
        vertcat(light{:}) vertcat(time{:}) vertcat(day{:}) logy logr factor*ones(size(logr,1),1)];
    
    N1 = sum(arrayfun(@(x) size(x.intervals_ti,1)>1, data_fac));
    disp(['# tracks: ' num2str(N) ', # data points: ' num2str(size(logr,1)), ...
          ', # tracks with >1 r*: ' num2str(N1)]);
end

if writeFile
    T = array2table(data_write, ...
        'VariableNames',{'approach','landingSide','pattern','light','time','day','logy','logr','threshold'});
    writetable(T,r_file);
end

%% FUNCTIONS USED IN THIS SCRIPT
function [data_cmap, edges] = find_cmap_dataedges(data, cmap)
% Finds color values for each data point in data
% This is required because usually size(data,1) > size(cmap,1)

% Method 1 - Creating N bins of data and allocating each data point to that bin. N is size(cmap,1)-1
edges = round(linspace(min(data),max(data),size(cmap,1)+1),3); % # edges = # color bins + 1
edges(1) = floor(min(data)*1000)/1000;
edges(end) = ceil(max(data)*1000)/1000;
[data_cmap, edges] = discretize(data,edges);
end

function [data_cmap, edges] = find_cmap(data, cmap)
% Finds color values for each data point in data
% This is required because usually size(data,1) > size(cmap,1)

% Method 1 - Creating N bins of data and allocating each data point to that bin. N is size(cmap,1)-1
edges = round(linspace(min(data),max(data),size(cmap,1)+1),3); % # edges = # color bins + 1
edges(1) = floor(min(data)*1000)/1000;
edges(end) = ceil(max(data)*1000)/1000;
[data_cmap, edges] = discretize(data,edges);
end