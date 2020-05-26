
cam_names = {'Basler_22549584', 'Basler_22549585', 'Basler_22549587', 'Basler_22956425'};
cam_ids = [0:3];

radius = 0.5;
data2d_plothandle = figure;
for i=1:length(cam_names)
    [Image.(['X' num2str(i)]), Image.(['map' num2str(i)])] = imread(fullfile(rootDir, trackingDir, treatmentTrackingFiles(indx).name, 'images', [cam_names{i} '.png']));
    subplot(2,2,i); imshow([Image.(['X' num2str(i)]), Image.(['map' num2str(i)])]);
    
    centers = data2d_4video(data2d_4video(:,1)==cam_ids(i),5:6);
    radii = radius*ones(size(centers,1),1);
    viscircles(centers, radii,'Color','g');
%     axis square;
%     axis equal;
%     axis image;
end


%%

ids = unique(data3d_4video(:,1));
colors = brewermap(length(ids),'Dark2');
for i=1:length(ids)
    disp(i);
    id = ids(i);
    dummy = data3d_4video(data3d_4video(:,1)==id,:);
%     if length(dummy)< 200
%         disp('length of the track is shorter than 200 points. Not displaying it');
%         continue;
%     end
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
        [Image.(['X' num2str(j)]), Image.(['map' num2str(j)])] = imread(fullfile(rootDir, trackingDir, treatmentTrackingFiles(indx).name, 'images', [cam_names{j} '.png']));
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

%% Rough work

data3d = h5read(fullfile(rootDir, trackingDir, '20190708_140014.mainbrain.offline.h5'), '/kalman_estimates');
data3d = [double(data3d.obj_id)  double(data3d.frame) (data3d.timestamp) ...
                                                                                                    data3d.x data3d.y data3d.z ...
                                                                                                    data3d.xvel data3d.yvel data3d.zvel ...
                                                                                                    data3d.P00 data3d.P01 data3d.P02 ...
                                                                                                    data3d.P11 data3d.P12 data3d.P22 ...
                                                                                                    data3d.P33 data3d.P44 data3d.P55];
                                                                                                
 data3d_4video = data3d(...
                        data3d(:,3)>=videoFile.timestamps(1) & ...
                        data3d(:,3)<=videoFile.timestamps(end), :);
                    
                    
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