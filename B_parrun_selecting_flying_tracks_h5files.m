%% To select (from h5 files) and write flying tracks (in h5_flying_kalman_estimates.csv files)
clc; clear; close all;

addpath('./lib/flymovieformat');
addpath('./lib/DrosteEffect-BrewerMap-221b913');

%% tracking file
if isunix
    path = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/';
elseif ispc
    path = 'D:/steady_wind_experiments/postprocessing/';
end

% paths/locations in h5 files
data3d_path = '/kalman_estimates';
data2d_path = '/data2d_distorted';

% data_association_path = '/data_association';

data3d_outputpath = '/flying_kalman_estimates';
file3d_output = 'h5_flying_kalman_estimates.csv';

% h5_pattern = {'20190705*.mainbrain.offline.h5', ...
%               '20190708*.mainbrain.offline.h5', ...
%               '20190709*.mainbrain.offline.h5', ...
%               '20190710*.mainbrain.offline.h5'};

h5_pattern = {'*.mainbrain.offline.h5'};

% Finding list (cell array) of required mainbrain.offline folders
% folders = GetFolders(path,{'.mainbrain.offline'});
% folders = {'20190706_123007.mainbrain.offline'};

% folders = dir(fullfile(path, h5_pattern));
folders = cellfun(@(x) dir(fullfile(path, x)), h5_pattern, 'UniformOutput', false);
folders = vertcat(folders{:});

% Read h5 files created offline, select flying tracks and write new data to
% h5_flying_kalman_estimates.csv
for ct_folder = 1:length(folders)
    disp(['Opening ' folders(ct_folder).name]);
    
    tic;
    
    % Read the data file
    data3d_original = h5read(fullfile(folders(ct_folder).folder, folders(ct_folder).name), data3d_path);
    
    data3d = data3d_original; % structure is copied (not from a handle class)
    
    nTracks = max(data3d_original.obj_id)+1; % number of tracks/objects in the file; contains 0 as obj_id 
    trackStatus = zeros(nTracks,1); % 0 - keep it, 1 - delete it, 2 - contains points outside physical domain
    
    % Apply selection filters
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    step = nTracks-1; % Code is NOT yet compatible with any other step
    if rem(nTracks-1, step) == 0
        ct = [0:step:(nTracks-1)];
    else
        ct = [[0:step:(nTracks-1)], (nTracks-1)];
    end
%     ct(1) = 1;
    
    tic
    for i=1:(length(ct)-1)
        % reducing communication overhead
%         M = 0;
        
        obj_ids = data3d.obj_id;
        xyz = [data3d.x data3d.y data3d.z];
        parfor ct_track=ct(i):1:ct(i+1)
            track_flags = obj_ids==ct_track;
            track_flags_indices = find(track_flags);
            
            track = xyz(track_flags,:);
            trackLength = length(track);
            
            in = IsInsideBox(track, [-0.1, 1.1, 0.0, 0.48, 0.0, 0.48]);
            out = ~in;
            if sum(out)/trackLength >= 0.33
                trackStatus(ct_track+1) = 1;
                continue;
            elseif sum(out) >= 1
                track(out,:) = [];
                trackLength = length(track);
                
%                 track_flags(track_flags_indices(out)) = [];
%                 track_flags_indices(out) = [];
                trackStatus(ct_track+1) = 2;
            end
            
            nearLeftWall = track(:,2) <= 0.06;
            nearRightWall = track(:,2) >= 0.42;
            if sum(nearLeftWall)/trackLength >= 0.95 || ...
               sum(nearRightWall)/trackLength >= 0.95
                trackStatus(ct_track+1) = 1;
                continue;
            end
        
            nearTop = track(:,3) >= 0.43;
            nearBottom = track(:,3) <= 0.05;
            if sum(nearTop)/length(nearTop) >= 0.999 || sum(nearBottom)/length(nearBottom) >= 0.999
                trackStatus(ct_track+1) = 1;
                continue;
            end
            
%             if (length(track) <= 150 && trackStatus(ct_track+1)==2) || length(track) <= 50
%                 trackStatus(ct_track+1) = 1;
%                 continue;
%             end
            if length(track) <= 50
                trackStatus(ct_track+1) = 1;
                continue;
            end
        end
        % update data3d variable
        tracks2del = find(trackStatus==1);
        obj_ids2del = tracks2del-1;
        rows_to_remove = find(ismember(data3d.obj_id, obj_ids2del) | ~NanIsInsideBox([data3d.x data3d.y data3d.z], [-0.1, 1.1, 0.0, 0.48, 0.0, 0.48]));
        data3d = structfun(@(x) (removerows(x, 'ind', rows_to_remove)), data3d, 'UniformOutput', false);
        
       
%         keyboard;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Write flying_kalman_estimates
    % Write flying_kalman_estimates.csv
    writecell(fieldnames(data3d)', fullfile(folders(ct_folder).folder, folders(ct_folder).name(1:end-3), file3d_output));
    data2write = [double(data3d.obj_id)  double(data3d.frame) (data3d.timestamp) ...
                                                                                                    data3d.x data3d.y data3d.z ...
                                                                                                    data3d.xvel data3d.yvel data3d.zvel ...
                                                                                                    data3d.P00 data3d.P01 data3d.P02 ...
                                                                                                    data3d.P11 data3d.P12 data3d.P22 ...
                                                                                                    data3d.P33 data3d.P44 data3d.P55];
    dlmwrite(fullfile(folders(ct_folder).folder, folders(ct_folder).name(1:end-3), file3d_output), data2write , '-append', 'precision', 20);
    
%     fileattrib(fullfile(folders(ct_folder).folder, folders(ct_folder).name),'+w');
%     h5create(fullfile(folders(ct_folder).folder, folders(ct_folder).name), [data3d_outputpath '/obj_id'], [length(flying_kalman_estimates.obj_id) 1], 'Datatype', 'uint64');
%     h5write(fullfile(folders(ct_folder).folder, folders(ct_folder).name), [data3d_outputpath '/obj_id'], flying_kalman_estimates.obj_id);
    
    
    toc
    % display done message
    disp(['Finished with selecting tracks for ' folders(ct_folder).name]);
end

%% GO BELOW only for check
data3d = readmatrix(fullfile(path,folders{1},file3d),'Range',2);
data2d = readmatrix(fullfile(path,folders{1},file2d),'Range',2);

data_association = readmatrix(fullfile(path,folders{1},data_association),'Range',2);

%% video file
file = '/media/reken001/Disk_08_backup/light_intensity_experiments/preprocessing/2019/07/06/123533-52287-Basler_22549584.fmf';

[video_data, timestamps] = fmf_read(file);


%% Plotting 2d feature data over images
[~, timestamps] = fmf_read(file);
data2d_4video = data2d(data2d(:,3)>=timestamps(1) & data2d(:,3)<=timestamps(end), [1:6]);
path = ['/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/' folders{1}];

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

%% Display 2d-orientation over the video data for a particular track
[video_data, timestamps] = fmf_read(file);
data3d_4video = data3d(data3d(:,3)>=timestamps(1) & data3d(:,3)<=timestamps(end), :);
data2d_4video = data2d(data2d(:,3)>=timestamps(1) & data2d(:,3)<=timestamps(end), [1:10]);
data_association_4video = data_association(min(data2d_4video(:,2)) <= data_association(:,2) & ...
                                           data_association(:,2) <= max(data2d_4video(:,2)),:);
                                      

nImages = 10; % # of randomnly selected images for which orientation will be plotted
trackID = 816;
frames = sort(datasample(1:length(timestamps),nImages,'Replace',false));
track_4video = sortrows(data3d_4video(data3d_4video(:,1)==trackID,:), [3,2]);

camID = 0; % for 84 camera
for ct=1:length(frames)
    orientationImage = figure;
    imshow(video_data(:,:,frames(ct))); hold on;
    
    % Identify center and 2d-slope 
    timestamp = timestamps(frames(ct));
%     [minabsdiff, indx] = min(abs(track_4video(:,3)-timestamp));
%     minabsdiff
%     frame = track_4video(indx,2);
    if timestamp >= track_4video(1,3) && timestamp <= track_4video(end,3)
        frame = track_4video(abs(track_4video(:,3)-timestamp)<=1e-5,2); % multiple frames
        frame = frame(1);
        pointID = data_association_4video(data_association_4video(:,1)==trackID & ...
                                          data_association_4video(:,2)==frame & ...
                                          data_association_4video(:,3)==camID,[4]);
        if ~isempty(pointID)
    %         point2D = data2d_4video(data2d_4video(:,2)==frame & ...
    %                            data2d_4video(:,1)==camID & ...
    %                            data2d_4video(:,10)==pointID(1),5:9);
    %                                        
    %     % 
    %         drawline(point2D(1:2), point2D(4), orientationImage);

    %     [minabsdiff, indx] = min(abs(data2d_4video(:,3)-timestamp));
    %     minabsdiff
    %     frame = track_4video(indx,2);
            point2D = data2d_4video(data2d_4video(:,3)==timestamp & ...
                               data2d_4video(:,1)==camID & ...
                               data2d_4video(:,10)==pointID(1),5:9);

        % 
            drawline(point2D(1:2), point2D(4), orientationImage);

            keyboard;
        end
    end
    close(orientationImage);
end

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
data3d_4video = data3d(data3d(:,3)>=timestamps(1) & data3d(:,3)<=timestamps(end), :);
data2d_4video = data2d(data2d(:,3)>=timestamps(1) & data2d(:,3)<=timestamps(end), [1:10]);
data_association_4video = data_association(min(data2d_4video(:,2)) <= data_association(:,2) & ...
                                           data_association(:,2) <= max(data2d_4video(:,2)),:);
% data3d_4video(data3d_4video(:,4)>=1.2 | data3d_4video(:,4)<=-0.2 | ...
%                                    data3d_4video(:,5)>=0.5 | data3d_4video(:,5)<=-0.0 | ...
%                                    data3d_4video(:,6)>=0.5 | data3d_4video(:,6)<=-0.0,:) = [];

% 
ids = unique(data3d_4video(:,1));
nTracks = length(ids);
for ct_track = 1:nTracks
        track_flags = data3d_4video(:,1)==ids(ct_track);
        track = data3d_4video(track_flags,:);
        
        % If more than 33% of the track is outside physical dimensions of
        % wind-tunnel, delete it
        in = IsInsideBox(track(:,4:6), [-0.1, 1.1, 0.0, 0.48, 0.0, 0.48]);
        out = ~in;
        if sum(out)/length(out) >= 0.33
            data3d_4video(track_flags,:) = [];
            disp(['1. Deleting track ' num2str(ids(ct_track))]);
            continue; % start next iteration
        elseif sum(out) >= 1
            % Delete points that lie outside physical dimensions of the tunnel
%             track(out,:) = [];
            track_flags_indices = find(track_flags);
%             track_flags(track_flags_indices(out)) = [];
            
            data3d_4video(track_flags_indices(out),:) = [];
            
            track_flags = data3d_4video(:,1)==ids(ct_track);
            track = data3d_4video(track_flags,:);
            
            disp(['1. Removed points from track ' num2str(ids(ct_track))]);
        end

        % If more than 75% of the track is within 6cms from side walls, delete
        % it
        nearLeftWall = track(:,5) <= 0.06;
        nearRightWall = track(:,5) >= 0.42;
        if sum(nearLeftWall)/length(nearLeftWall) >= 0.75 || ...
           sum(nearRightWall)/length(nearRightWall) >= 0.75
            data3d_4video(data3d_4video(:,1)==ids(ct_track),:) = [];
            display(['2. Deleting track ' num2str(ids(ct_track))]);
            continue;
        end
        
        % If more than 90% of the track is within 5cms from top and botttom
        % walls, delete it
        nearTop = track(:,6) >= 0.43;
        nearBottom = track(:,6) <= 0.05;
        if sum(nearTop)/length(nearTop) >= 0.90 || sum(nearBottom)/length(nearBottom) >= 0.90
            data3d_4video(data3d_4video(:,1)==ids(ct_track),:) = [];
            display(['3. Deleting track ' num2str(ids(ct_track))]);
            continue;
        end

        % If length of the track is less than 150 points delete it
        if length(track) <= 100
            data3d_4video(data3d_4video(:,1)==ids(ct_track),:) = [];
            display(['4. Deleting track ' num2str(ids(ct_track))]);
            continue;
        end

        % based on velocity covariance
        % based on position covariance


end
    

ids = unique(data3d_4video(:,1));
colors = brewermap(length(ids),'Dark2');
for i=1:length(ids)
    disp(i);
    id = ids(i);
    dummy = data3d_4video(data3d_4video(:,1)==id,:);
    if length(dummy)< 200
        disp('length of the track is shorter than 200 points. Not displaying it');
        continue;
    end
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
    % Plotting position
    pos_plothandle = figure;
    subplot(3,1,1);
    plot(dummy(:,3), dummy(:,4),'-o');
    ylabel('X');
    subplot(3,1,2);
    plot(dummy(:,3), dummy(:,5),'-o');
    ylabel('Y');
    subplot(3,1,3);
    plot(dummy(:,3), dummy(:,6),'-o');
    ylabel('Z');
    
    % Plotting velocity
    vel_plothandle = figure;
    subplot(3,1,1);
    plot(dummy(:,3), dummy(:,7),'-o');
    ylabel('V_{x}');
    subplot(3,1,2);
    plot(dummy(:,3), dummy(:,8),'-o');
    ylabel('V_{y}');
    subplot(3,1,3);
    plot(dummy(:,3), dummy(:,9),'-o');
    ylabel('V_{z}');
    
    % Plotting position covariance
    posCov_plothandle = figure;
    subplot(3,1,1);
    plot(dummy(:,3), dummy(:,10),'-o');
    ylabel('\sigma_{xx}');
    subplot(3,1,2);
    plot(dummy(:,3), dummy(:,13),'-o');
    ylabel('\sigma_{yy}');
    subplot(3,1,3);
    plot(dummy(:,3), dummy(:,15),'-o');
    ylabel('\sigma_{zz}');
    
    % Plotting position covariance
    velCov_plothandle = figure;
    subplot(3,1,1);
    plot(dummy(:,3), dummy(:,16),'-o');
    ylabel('\sigma_{uu}');
    subplot(3,1,2);
    plot(dummy(:,3), dummy(:,17),'-o');
    ylabel('\sigma_{vv}');
    subplot(3,1,3);
    plot(dummy(:,3), dummy(:,18),'-o');
    ylabel('\sigma_{ww}');
    
    keyboard;
    
    close(data2d_plothandle);
    close(trajFeaturePlot);
    close(pos_plothandle);
    close(vel_plothandle);
    close(posCov_plothandle);
    close(velCov_plothandle);
    
    
end
%% function used in the script
function drawline(center, slope_deg, figure_handle)
    figure(figure_handle); hold on;
    viscircles(center, 1,'Color','red');
    
    
    x = [center(1)+20*cos(slope_deg), center(1)-20*cos(slope_deg)];
    y = [center(2)+20*sin(slope_deg), center(2)-20*sin(slope_deg)];
    line(x,y,'Color','green','LineWidth', 2);
    
    x = [center(1)+20*cos(-slope_deg), center(1)-20*cos(-slope_deg)];
    y = [center(2)+20*sin(-slope_deg), center(2)-20*sin(-slope_deg)];
    line(x,y,'Color','red','LineWidth', 2);
%     
%     x = [center(1)+20*cos(-slope_deg+pi/2), center(1)-20*cos(-slope_deg+pi/2)];
%     y = [center(2)+20*sin(-slope_deg+pi/2), center(2)-20*sin(-slope_deg+pi/2)];
%     line(x,y,'Color','white','LineWidth', 2);
%     
%     x = [center(1)+20*cos(slope_deg+pi/2), center(1)-20*cos(slope_deg+pi/2)];
%     y = [center(2)+20*sin(slope_deg+pi/2), center(2)-20*sin(slope_deg+pi/2)];
%     line(x,y,'Color','magenta','LineWidth', 2);
    

end
function in = IsInsideBox(data,box)
    % nans are counted outside box
    % data = [X Y Z]
    % box = [xmin xmax ymin ymax zmin zmax]
    in = [(data(:,1) >= box(1) & data(:,1) <= box(2)) & ...
          (data(:,2) >= box(3) & data(:,2) <= box(4)) & ...
          (data(:,3) >= box(5) & data(:,3) <= box(6))];        
end
function in = NanIsInsideBox(data,box)
    % Nans are counted inside box
    % data = [X Y Z]
    % box = [xmin xmax ymin ymax zmin zmax]
    in = [(data(:,1) >= box(1) & data(:,1) <= box(2)) & ...
          (data(:,2) >= box(3) & data(:,2) <= box(4)) & ...
          (data(:,3) >= box(5) & data(:,3) <= box(6))]; 
    in(isnan(data(:,1)) | isnan(data(:,2)) | isnan(data(:,1))) = true;
end

% Functions used in this script file

function [folders] = GetFolders(direc,pattern)
    % Get a list of all folders in a directory containing patterns in their names.
    
    % direc - path where to look
    % pattern - cell array containing string patterns that must be in the
    % names of the folders

    all    = dir(direc);
    folders = all([all(:).isdir]);
    names    = {folders.name};

    % Get a logical vector that discards folders that aren't required
    flags = ~strcmp(names, '.') & ...
              ~strcmp(names, '..') & ~strcmp(cellfun(@(x){x(1)},names),'~');%~cellfun('isempty',strfind(A,B));
    for ct=1:length(pattern)
       flags = flags & contains(names,pattern{ct});
    end

    % Extract only those that contain pattern
    folders = names(flags);
end