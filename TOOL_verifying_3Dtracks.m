%% 
clc; clear; close all;
addpath('/home/reken001/Pulkit/lib/flymovieformat');
addpath('/home/reken001/Pulkit/lib/DrosteEffect-BrewerMap-221b913');

%% video file
file = '/media/reken001/Disk_09/unsteady_wind_experiments/Videos/2019/07/28/091348-105845-Basler_22549584.fmf';
file = '/media/reken001/Disk_09/unsteady_wind_experiments/Videos/2019/07/28/092136-187644-Basler_22549584.fmf';
[video_data, timestamps] = fmf_read(file);

%% tracking file
path = '/media/reken001/Disk_09/unsteady_wind_experiments/test/postprocessing/20190723_110007.mainbrain.offline/';
path = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/20190723_110007.mainbrain.offline/';
% path = '/media/reken001/Disk_09/unsteady_wind_experiments/test/Tracking/20190728_080818.mainbrain/';
file3d = 'kalman_estimates.csv';
file2d = 'data2d_distorted.csv';
data_association = 'data_association.csv';

Data3d = importdata([path file3d]);
data3d = Data3d.data;
clear Data3d;

Data2d = importdata([path file2d]);
data2d = Data2d.data;
clear Data2d;

Data_association = importdata([path data_association]);
data_association = Data_association.data;
clear Data_association;

%%
data3d_4video = data3d(data3d(:,3)>=timestamps(1) & data3d(:,3)<=timestamps(end), [1:9]);
data2d_4video = data2d(data2d(:,3)>=timestamps(1) & data2d(:,3)<=timestamps(end), [1:6]);
% Plotting 2d feature data over images
cam_names = {'Basler_22549584', 'Basler_22549585', 'Basler_22549587', 'Basler_22956425'};
cam_ids = [0:3];

radius = 0.5;
data2d_plothandle = figure;
for i=1:length(cam_names)
    [Image.(['X' num2str(i)]), Image.(['map' num2str(i)])] = imread(fullfile(path, 'images', [cam_names{i} '.png']));
    subplot(2,2,i); imshow([Image.(['X' num2str(i)]), Image.(['map' num2str(i)])]);
    
    centers = data2d_4video(data2d_4video(:,1)==cam_ids(i),5:6);
    radii = radius*ones(size(centers,1),1);
    viscircles(centers, radii,'Color','g');
%     axis square;
%     axis equal;
%     axis image;
end

%% Plotting the tracks
trajFeaturePlot = figure;
set(gca, 'FontSize', 18);

% delete data from data3d_video that is out of the wind tunnel/camera-views

% Method 1
% unrealistic_tracks = data3d_4video(data3d_4video(:,4)>=1.2 | data3d_4video(:,4)<=-0.2 | ...
%                                    data3d_4video(:,5)>=0.5 | data3d_4video(:,5)<=0.0 | ...
%                                    data3d_4video(:,6)>=0.5 | data3d_4video(:,6)<=0.0);
% data3d_4video(ismember(data3d_4video(:,1),unique(unrealistic_tracks(:,1))),:) = [];

% Method 2
[video_data, timestamps] = fmf_read(file);
data3d_4video = data3d(data3d(:,3)>=timestamps(1) & data3d(:,3)<=timestamps(end), [1:9]);
data2d_4video = data2d(data2d(:,3)>=timestamps(1) & data2d(:,3)<=timestamps(end), [1:6]);
% data3d_4video(data3d_4video(:,4)>=1.2 | data3d_4video(:,4)<=-0.2 | ...
%                                    data3d_4video(:,5)>=0.5 | data3d_4video(:,5)<=-0.0 | ...
%                                    data3d_4video(:,6)>=0.5 | data3d_4video(:,6)<=-0.0,:) = [];

% 
ids = unique(data3d_4video(:,1));
hold on;
colors = brewermap(length(ids),'Dark2');
for i=1:length(ids)
    id = ids(i);
    dummy = data3d_4video(data3d_4video(:,1)==id,:);
    plot3(dummy(:,4),dummy(:,5),dummy(:,6),'o','Color',colors(i,:),'LineWidth',2,'MarkerSize',1,'DisplayName',['id: ' num2str(id)]);
end
figure(trajFeaturePlot);
xlabel('x (m)', 'FontSize', 16);
ylabel('y (m)', 'FontSize', 16);
zlabel('z (m)', 'FontSize', 16);
set(gca, 'FontSize', 18); grid on;
legend(gca,'show','Location','best');
xlim([0 1]);
ylim([0 0.5]);

%% Plotting the 3-d tracks alongwith associated features in 2-D cameras views that made the 3d-track 

set(gca, 'FontSize', 18);

% delete data from data3d_video that is out of the wind tunnel/camera-views

% Method 1
% unrealistic_tracks = data3d_4video(data3d_4video(:,4)>=1.2 | data3d_4video(:,4)<=-0.2 | ...
%                                    data3d_4video(:,5)>=0.5 | data3d_4video(:,5)<=0.0 | ...
%                                    data3d_4video(:,6)>=0.5 | data3d_4video(:,6)<=0.0);
% data3d_4video(ismember(data3d_4video(:,1),unique(unrealistic_tracks(:,1))),:) = [];

% Method 2
[video_data, timestamps] = fmf_read(file);
data3d_4video = data3d(data3d(:,3)>=timestamps(1) & data3d(:,3)<=timestamps(end), [1:9]);
data2d_4video = data2d(data2d(:,3)>=timestamps(1) & data2d(:,3)<=timestamps(end), [1:10]);
data_association_4video = data_association(min(data2d_4video(:,2)) <= data_association(:,2) & ...
                                           data_association(:,2) <= max(data2d_4video(:,2)),:);
data3d_4video(data3d_4video(:,4)>=1.2 | data3d_4video(:,4)<=-0.2 | ...
                                   data3d_4video(:,5)>=0.5 | data3d_4video(:,5)<=-0.0 | ...
                                   data3d_4video(:,6)>=0.5 | data3d_4video(:,6)<=-0.0,:) = [];

% 
ids = unique(data3d_4video(:,1));
hold on;
colors = brewermap(length(ids),'Dark2');
for i=1:length(ids)
    id = ids(i);
    dummy = data3d_4video(data3d_4video(:,1)==id,:);
    
    trajFeaturePlot = figure;
    plot3(dummy(:,4),dummy(:,5),dummy(:,6),'-o','Color',colors(i,:),'LineWidth',2,'MarkerSize',1,'DisplayName',['id: ' num2str(id)]);
    set(gca, 'FontSize', 18); grid on;
    view(0,90);
    xlim([0 1]);
    ylim([0 0.5]);
    legend(gca,'show','Location','best');
    % Highlighting associated 2-D features
    data2d_plothandle = figure;
    for j=1:length(cam_names)
        % read images and plot all features
        [Image.(['X' num2str(j)]), Image.(['map' num2str(j)])] = imread(fullfile(path, 'images', [cam_names{j} '.png']));
        subplot(2,2,j); imshow([Image.(['X' num2str(j)]), Image.(['map' num2str(j)])]);

        centers = data2d_4video(data2d_4video(:,1)==cam_ids(j),5:6);
        radii = radius*ones(size(centers,1),1);
        viscircles(centers, radii,'Color','w');
        
        % plot 3d-track associated features
        frame_pointID = data_association_4video(data_association_4video(:,1)==ids(i) & ...
                                                data_association_4video(:,3)==cam_ids(j),[2 4]);
        centers = zeros(length(frame_pointID),2);
        
        for k=1:size(frame_pointID,1)
            centers(k,1:2) = data2d_4video(data2d_4video(:,2)==frame_pointID(k,1) & ...
                                           data2d_4video(:,1)==cam_ids(j) & ...
                                           data2d_4video(:,10)==frame_pointID(k,2),5:6);
        end
        radii = radius*ones(size(centers,1),1);
        viscircles(centers, radii,'Color',colors(i,:));
    end
    pause;
    
    close(data2d_plothandle);
    close(trajFeaturePlot);
    
    
end
