%% Plotting trajectories and mean approach for ALL tracks
%

% Extract all approaches
close all; clc;
% clear;



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

%%
inputFile = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/BlindLandingtracks_A3_LDF.mat';
inputFile = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/BlindLanding_and_aborted_tracks_A3_LDF.mat';
% load(inputFile);

winds = unique([treatments.wind]);
behaviour = {'rising','constant','sleeping'};
data_all = struct.empty;

Nperday = zeros(length(winds),1); % no. of state LDFs per day in different winds
Nperday_aborted = zeros(length(winds),1); % no. of aborted state LDFs per day in different winds
NrelevantTreatments = zeros(length(winds),1);
for ct_wind = 1:length(winds)
        for ct_behaviour = 2%1:length(behaviour)
            clear dummy;
            
            disp(['Wind: ' num2str(winds(ct_wind)) ...
                  ', behaviour: ' behaviour{ct_behaviour}]);
              
            % Selecting relevant treatments
            if strcmpi(behaviour{ct_behaviour}, 'rising')
                relevantTreatments = treatments(rem(1:length(treatments), 8)==1);
            elseif strcmpi(behaviour{ct_behaviour}, 'constant')
                relevantTreatments = treatments( [treatments.wind] == winds(ct_wind) & ...
                                     rem(1:length(treatments), 8)>1 & ...
                                     rem(1:length(treatments), 8)<8);
                hasUniformHwData = arrayfun(@(x) x.hwData.hasUniformHwData,relevantTreatments);
                relevantTreatments = relevantTreatments(hasUniformHwData);
%                 landingTracks = [relevantTreatments.landingTracks];
%                                  
%                 % # of landing tracks
%                 disp(['# of distinct Flydra objects (landingTracks): ' num2str(length(landingTracks))]);
            elseif strcmpi(behaviour{ct_behaviour}, 'sleeping')
                relevantTreatments = treatments(rem(1:length(treatments), 8)==0);
            else
                error('What other treatments did you perform dude?')
            end
            hasUniformHwData = arrayfun(@(x) x.hwData.hasUniformHwData,relevantTreatments);
            relevantTreatments = relevantTreatments(hasUniformHwData);
            
            landingTracks = [relevantTreatments.landingTracks];
            abortedLandingTracks = [relevantTreatments.abortedLandingTracks];
            
            startTimes = arrayfun(@(x) (x.startTime)*ones(length([x.landingTracks.state_LDF]), 1), ...
                relevantTreatments, 'UniformOutput', false);
            startTimes = vertcat(startTimes{:});
            
            days = arrayfun(@(x) (x.datenum)*ones(length([x.landingTracks.state_LDF]), 1), ...
                relevantTreatments, 'UniformOutput', false);
            days = vertcat(days{:});
            
            % # of landing tracks
            disp(['# of distinct Flydra objects (landingTracks): ' num2str(length(landingTracks))]);
            
            % # of distinct state LDFs
            disp(['# of distinct state LDFs: ' num2str(length([landingTracks.state_LDF]))]);
            
            
            dummy.state_LDF = [landingTracks.state_LDF];
%             dummy.pattern = (ct_pattern)*ones(1, length(dummy.state_LDF));
%             dummy.light = (ct_light)*ones(1, length(dummy.state_LDF));
            dummy.wind = (ct_wind)*ones(1, length(dummy.state_LDF));
            dummy.day = days';
            dummy.time = startTimes';
            
            data_all = [data_all; dummy];
            
            Nperday(ct_wind) = length([landingTracks.state_LDF]);
            Nperday_aborted(ct_wind) = length([abortedLandingTracks.state_LDF]);
            NrelevantTreatments = length(relevantTreatments);
        end

end
% keyboard;
figure;
% subplot(2,1,1);
bar(winds, Nperday./NrelevantTreatments);
ylabel('Avg. no. of landings per day', 'FontSize', 14);
xlabel('Wind conditions', 'FontSize', 14);
set(gca, 'FontSize', 16);
ylim([0 500]);
figure;
% subplot(2,1,2);
bar(winds, Nperday_aborted./NrelevantTreatments);
ylabel('Avg. no. of aborted landings per day', 'FontSize', 14);
xlabel('Wind conditions', 'FontSize', 14);
set(gca, 'FontSize', 16);
ylim([0 500]);
figure;
bar(winds, (Nperday_aborted)./(Nperday_aborted+Nperday)*100);
ylabel('Avg. % of aborted landings per day', 'FontSize', 14);
xlabel('Wind conditions', 'FontSize', 14);
set(gca, 'FontSize', 16);

% %  Write file for analysis in R for 10005 tracks
writeFile = false;
r_file = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/data_all_trajs_Rstudio.txt';
data_write = [];
approach_no = 0;
yrange = [-0.12 -0.02];
if writeFile
    for ct=1:length(data_all)
        for ct1=1:1:length(data_all(ct).state_LDF)
            dummy_y = data_all(ct).state_LDF(ct1).filteredState(:,3); % the complete trajectory
%             [~, ymin_indx] = min(data_all(ct).state_LDF(ct1).filteredState(:,3));
% 
%             % create other variables
%             y = -data_all(ct).state_LDF(ct1).filteredState(ymin_indx:end, 3);
%             v = data_all(ct).state_LDF(ct1).filteredState(ymin_indx:end, 6);
            
            y = -data_all(ct).state_LDF(ct1).filteredState(dummy_y >= yrange(1) & dummy_y <=yrange(2), 3);
            v = data_all(ct).state_LDF(ct1).filteredState(dummy_y >= yrange(1) & dummy_y <=yrange(2), 6);
            r = v./y;

            N = length(r);
            % Create nominal and ordinal variables

            data_all(ct).state_LDF(ct1).setLandingSide(); % to store landing side in the rrefSegments

            approach_no = approach_no + 1;
            approach = approach_no*ones(N, 1);
            if strcmpi(data_all(ct).state_LDF(ct1).landingSide, 'Hive')
                side = 1*ones(N, 1);
            elseif strcmpi(data_all(ct).state_LDF(ct1).landingSide, 'Feeder')
                side = 2*ones(N, 1);
            end

            pattern = data_all(ct).pattern(ct1)*ones(N, 1);
            light = data_all(ct).light(ct1)*ones(N, 1);
            time = data_all(ct).time(ct1)*ones(N, 1);
            day = data_all(ct).day(ct1)*ones(N, 1);



            %     logy = log(y);
            %     logr = log(V);

            data_write = [data_write; ...
                approach side pattern light time day y r v];

        end
    end
end

if writeFile
    T = array2table(data_write, ...
        'VariableNames',{'approach','landingSide','pattern','light','time','day','y','r','v'});
    writetable(T,r_file);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
close all;

% Plot all individual tracks (every 70th track)
% Find Vrange
skip_step = 150;
Vrange = [];
for ct=1:length(data_all)
    for ct1=1:skip_step:length(data_all(ct).state_LDF)
        xyz = data_all(ct).state_LDF(ct1).filteredState(:,[2 3 4]); % the complete trajectory
        [~, ymin_indx] = min(xyz(:,2));
        
        Vrange(end+1,:) = [min(data_all(ct).state_LDF(ct1).filteredState(ymin_indx:end,6)) max(data_all(ct).state_LDF(ct1).filteredState(ymin_indx:end,6))]; 
    end
end
Vrange = [min(Vrange(:,1)) max(Vrange(:,2))];
Vrange = [floor(Vrange(1)*10)/10 ceil(max(Vrange(2))*10)/10];

% Choose colormap and find data edges for V in Vrange
cmap = jet(round(diff(Vrange)/0.05));
edges = round(linspace(Vrange(1),Vrange(2),size(cmap,1)+1),2); % # edges = # color bins + 1

% [data_cmap, edges] = discretize(data,edges);
% 
% [data_cmap, edges] = find_cmap(yTravelled, red_cmap);
% f = scatter(-ymean,-rref_fitVvsy,40,red_cmap(data_cmap',:),'filled','s');


trajPlot = figure; hold on;
colormap(cmap);
tracks_plotted = 0;
for ct=1:length(data_all)
    for ct1=1:skip_step:length(data_all(ct).state_LDF)
        % y is pointing inside the platform here, therefore looking at ymin
        % for maximum distance away from the platform
        xyz = data_all(ct).state_LDF(ct1).filteredState(:,[2 3 4]); % the complete trajectory
        [~, ymin_indx] = min(xyz(:,2));
        
        [data_cmap, ~] = discretize(data_all(ct).state_LDF(ct1).filteredState(ymin_indx:end,6), edges);
        scatter3(xyz(ymin_indx:end,1), -1*xyz(ymin_indx:end,2), xyz(ymin_indx:end,3), 2, cmap(data_cmap',:),'filled','o');
%         plot3(xyz(ymin_indx:end,1), -1*xyz(ymin_indx:end,2), xyz(ymin_indx:end,3),'Color',[117 153 242]./255); 
        
        % Use discretize to find correspondece within colormap as per V
        
        % Plot the approach between max distance to min distance
%         plot(-data_all(ct).state_LDF(ct1).filteredState(ymin_indx:end,3), data_all(ct).state_LDF(ct1).filteredState(ymin_indx:end,6),'Color',[189,189,189]/255)
%         plot(data(ct).tracks_fac(ct1).filteredState(:,3), data(ct).tracks_fac(ct1).filteredState(:,6),'Color',[189,189,189]/255)

        tracks_plotted = tracks_plotted + 1;
    end
end
% % Landing Disc
% [X,Y,Z] = cylinder(treatment.landingDiscs(1).radius);
radius = treatments(1).landingDiscs(1).radius;

figure(trajPlot);

cmap_bar = colorbar('eastoutside');
caxis(edges([1 end]));
cmap_bar.Ticks = [edges(1) 0 edges(end)];
cmap_bar.Label.String = 'V';
cmap_bar.Label.FontSize = 16;

% h=mesh(X,Y,Z,'facecolor',[1 0 0]); % draw landing disc
plot3([-radius; radius], [0 0], [0 0], 'LineWidth', 2, 'Color', [83 83 83]./255);
plot3([0 0], [0 0], [-radius; radius], 'LineWidth', 2, 'Color', [83 83 83]./255);
zlabel('z (m)', 'FontSize', 14);
ylabel('y (m)', 'FontSize', 14);
xlabel('x (m)', 'FontSize', 14);
axis equal;
xlim([-0.4 0.4]);
xticks([-0.4:0.2:0.4]);
ylim([0 0.32]);
yticks([0:0.1:0.4]);
zlim([-0.25 0.25]);
zticks([-0.25:0.25:0.25]);
set(gca, 'FontSize', 16);

view(0,90);
view(-90,0);
% Collect tracks from ymax to end

for ct=1:length(data_all)
    data_all(ct).filteredStates_all = struct.empty;
    for ct1=1:1:length(data_all(ct).state_LDF)
        % y is pointing inside the platform here, therefore looking at ymin
        % for maximum distance away from the platform
        xyz = data_all(ct).state_LDF(ct1).filteredState(:,[2 3 4]); % the complete trajectory
        [~, ymin_indx] = min(data_all(ct).state_LDF(ct1).filteredState(:,3));
        data_all(ct).filteredStates_all(ct1).state = data_all(ct).state_LDF(ct1).filteredState(ymin_indx:end,:);
    end
end
% Plot mean track
ybins = -0.4:0.005:-0.01;
tracks = [data_all.filteredStates_all];
data_xyzuvw = arrayfun(@(x) x.state(x.state(:,3)>=ybins(1) & x.state(:,3)<=ybins(end),[2:7]), tracks, 'UniformOutput', false);
data_xyzuvw = vertcat(data_xyzuvw{:});
xyzuvw = [];
y = [];
sem_xyzuvw = []; % standard error of the means, SEM = std(data)/sqrt(length(data));
r = [];
sem_r = [];
for ct=1:length(ybins)-1
    dummy = data_xyzuvw(data_xyzuvw(:,2)>=ybins(ct) & data_xyzuvw(:,2)<ybins(ct+1), [1 2 3 4 5 6]);
    y = [y; mean(ybins(ct:ct+1))];
    xyzuvw = [xyzuvw; mean(dummy)];
    r = [r; mean(dummy(:,5)./dummy(:,2))];
    sem_xyzuvw = [sem_xyzuvw; std(dummy)/sqrt(size(dummy,1))];
%     sem_xyzuvw = [sem_xyzuvw; std(dummy)];
    sem_r = [sem_r; std(dummy(:,5)./dummy(:,2))/sqrt(size(dummy,1))];
%     sem_r = [sem_r; std(dummy(:,5)./dummy(:,2))];
end

figure(trajPlot);

% plot(-y,V+semV,'--k');
% plot(-y,V-semV,'--k');
% fill([-y fliplr(-y)],[V+semV fliplr(V-semV)], 'g');
plot3(xyzuvw(:,1),-xyzuvw(:,2),xyzuvw(:,3),'k','Linewidth',3)
set(gca, 'FontSize', 16);


figure; hold on;
% plot(-y,xyzuvw(:,4)+sem_xyzuvw(:,4),'--k');
% plot(-y,xyzuvw(:,4)-sem_xyzuvw(:,4),'--k');
fill([-y; flipud(-y)],[xyzuvw(:,4)+sem_xyzuvw(:,4); flipud(xyzuvw(:,4)-sem_xyzuvw(:,4))], 'g');
plot(-y,xyzuvw(:,4),'k','Linewidth',1);
ylabel('U (m/s)', 'FontSize', 16);
xlabel('y (m)', 'FontSize', 16);
set(gca, 'FontSize', 16);
ylim([-0.15 0.15]);
yticks([-0.15:0.15:0.15]);
xlim([0 0.32]);
% xticks([0:0.08:0.32]);

figure; hold on;
yrange = [0.02 0.12]; 
rref_R = 2.34; % from  analysis in R

subplot(2,1,1); hold on;
% plot(-y,xyzuvw(:,5)+sem_xyzuvw(:,5),'--k');
% plot(-y,xyzuvw(:,5)-sem_xyzuvw(:,5),'--k');
fill([-y; flipud(-y)],[xyzuvw(:,5)+sem_xyzuvw(:,5); flipud(xyzuvw(:,5)-sem_xyzuvw(:,5))],[252,187,161]./255, 'EdgeColor', [252,187,161]./255);
plot(-y,xyzuvw(:,5),'Color',[252,187,161]./255,'Linewidth',1);
ylabel('V (m/s)', 'FontSize', 16);
xlabel('y (m)', 'FontSize', 16);
set(gca, 'FontSize', 16);
% ylim([0 0.3]);
yticks([0:0.1:0.3]);
% xlim([0 0.32]);
y_in_yrange = -y(-y>=yrange(1) & -y<=yrange(2));
plot(y_in_yrange, xyzuvw(-y>=yrange(1) & -y<=yrange(2), 5),'.','MarkerSize',10,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
plot(y_in_yrange, rref_R*y_in_yrange,'LineWidth',2,'Color',[69 117 180]./255);
plot([0 y_in_yrange(end)],[0 rref_R*y_in_yrange(end)],'--','LineWidth',2,'Color',[69 117 180]./255);
title(['r* : ' num2str(rref_R,3)], 'FontSize', 16);

subplot(2,1,2); hold on;
a = fill([-y; flipud(-y)],[-r+sem_r; flipud(-r-sem_r)],[252,187,161]./255, 'EdgeColor', [252,187,161]./255);
plot(-y,-r,'Color',[252,187,161]./255,'Linewidth',1);
ylabel('r (s-1)', 'FontSize', 16);
xlabel('y (m)', 'FontSize', 16);
set(gca, 'FontSize', 16);
% ylim([0 0.3]);
% yticks([0:0.1:0.3]);
% xlim([0 0.32]);
plot(y_in_yrange, -r(-y>=yrange(1) & -y<=yrange(2)),'.','MarkerSize',10,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
plot(y_in_yrange, rref_R*ones(length(y_in_yrange),1),'LineWidth',2,'Color',[69 117 180]./255);
plot([0 y_in_yrange(end)],[rref_R rref_R],'--','LineWidth',2,'Color',[69 117 180]./255);


% xticks([0:0.08:0.32]);

figure; hold on;
% plot(-y,xyzuvw(:,6)+sem_xyzuvw(:,6),'--k');
% plot(-y,xyzuvw(:,6)-sem_xyzuvw(:,6),'--k');
fill([-y; flipud(-y)],[xyzuvw(:,6)+sem_xyzuvw(:,6); flipud(xyzuvw(:,6)-sem_xyzuvw(:,6))], 'g');
plot(-y,xyzuvw(:,6),'k','Linewidth',1);
ylabel('W (m/s)', 'FontSize', 16);
xlabel('y (m)', 'FontSize', 16);
set(gca, 'FontSize', 16);
ylim([-0.15 0.15]);
yticks([-0.15:0.15:0.15]);
xlim([0 0.32]);
xticks([0:0.08:0.32]);

figure; hold on;
plot(-y,xyzuvw(:,4)./xyzuvw(:,5),'k','Linewidth',1);
ylabel('U/V', 'FontSize', 16);
xlabel('y (m)', 'FontSize', 16);
set(gca, 'FontSize', 16);
xlim([0 0.32]);
xticks([0:0.08:0.32]);

figure; hold on;
plot(-y,xyzuvw(:,6)./xyzuvw(:,5),'k','Linewidth',1);
ylabel('W/V', 'FontSize', 16);
xlabel('y (m)', 'FontSize', 16);
set(gca, 'FontSize', 16);
xlim([0 0.32]);
xticks([0:0.08:0.32]);

xlim([0 0.4])
ylim([0 0.3])
yticks([0:0.1:0.3]);
view(0,90);
view(-90,0);

figure;
plot(-y, xyzuvw(:,4)./-y, 'k')

% % Plot V and y 3d histogram
% figure;
% hist3([-data_yV(:,1) data_yV(:,2)]);
% xlabel('y (m)', 'FontSize', 16);
% ylabel('V (m/s)', 'FontSize', 16);
% zlabel('Occurances', 'FontSize', 16);
% set(gca, 'FontSize', 16);
% view([-146 19])


%% %%%%%%%%%%% Plot average approach for different wind conditions
close all;
ybins = -0.4:0.005:-0.01;

cmap = jet(6);
cmap([2 5],:) = [];
colors = [cmap; 1 0 1; 1 0 0];

fig1 = figure; hold on;
winds = unique([data_all.wind]);
yrange = [0.02 0.12];
for ct_wind=1:length(data_all)
    tracks = arrayfun(@(x) data_all(x).filteredStates_all(data_all(x).wind == ct_wind),1:length(data_all), 'UniformOutput', false);
    tracks = [tracks{:}]; %[data_all([data_all.light] == ct_light).filteredStates_all];
    data_xyzuvw = arrayfun(@(x) x.state(x.state(:,3)>=ybins(1) & x.state(:,3)<=ybins(end),[2:7]), tracks, 'UniformOutput', false);
    data_xyzuvw = vertcat(data_xyzuvw{:});
    xyzuvw = [];
    y = [];
    sem_xyzuvw = []; % standard error of the means, SEM = std(data)/sqrt(length(data));
    r = [];
    sem_r = [];
    for ct=1:length(ybins)-1
        dummy = data_xyzuvw(data_xyzuvw(:,2)>=ybins(ct) & data_xyzuvw(:,2)<ybins(ct+1), [1 2 3 4 5 6]);
        y = [y; mean(ybins(ct:ct+1))];
        xyzuvw = [xyzuvw; mean(dummy)];
        r = [r; mean(dummy(:,5)./dummy(:,2))];
        sem_xyzuvw = [sem_xyzuvw; std(dummy)/sqrt(size(dummy,1))];
        %     sem_xyzuvw = [sem_xyzuvw; std(dummy)];
        sem_r = [sem_r; std(dummy(:,5)./dummy(:,2))/sqrt(size(dummy,1))];
        %     sem_r = [sem_r; std(dummy(:,5)./dummy(:,2))];
    end
    
    
    
    
    subplot(2,1,1); hold on;
    % plot(-y,xyzuvw(:,5)+sem_xyzuvw(:,5),'--k');
    % plot(-y,xyzuvw(:,5)-sem_xyzuvw(:,5),'--k');
    fill([-y; flipud(-y)],[xyzuvw(:,5)+sem_xyzuvw(:,5); flipud(xyzuvw(:,5)-sem_xyzuvw(:,5))], colors(ct_wind,:), 'EdgeColor', colors(ct_wind,:));
%     plot(-y,xyzuvw(:,5),'Color',[252,187,161]./255,'Linewidth',1);
    
%     y_in_yrange = -y(-y>=yrange(1) & -y<=yrange(2));
%     plot(y_in_yrange, xyzuvw(-y>=yrange(1) & -y<=yrange(2), 5),'.','MarkerSize',10,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
%     plot(y_in_yrange, rref_R*y_in_yrange,'LineWidth',2,'Color',[69 117 180]./255);
%     plot([0 y_in_yrange(end)],[0 rref_R*y_in_yrange(end)],'--','LineWidth',2,'Color',[69 117 180]./255);
%     title(['r* : ' num2str(rref_R,3)], 'FontSize', 16);
    
    subplot(2,1,2); hold on;
    a = fill([-y; flipud(-y)],[-r+sem_r; flipud(-r-sem_r)], colors(ct_wind,:), 'EdgeColor', colors(ct_wind,:));
%     plot(-y,-r,'Color',[252,187,161]./255,'Linewidth',1);
    
    % ylim([0 0.3]);
    % yticks([0:0.1:0.3]);
    xlim([0 0.32]);
%     plot(y_in_yrange, -r(-y>=yrange(1) & -y<=yrange(2)),'.','MarkerSize',10,'MarkerFaceColor',[215 48 39]./255, 'MarkerEdgeColor',[215 48 39]./255);
%     plot(y_in_yrange, rref_R*ones(length(y_in_yrange),1),'LineWidth',2,'Color',[69 117 180]./255);
%     plot([0 y_in_yrange(end)],[rref_R rref_R],'--','LineWidth',2,'Color',[69 117 180]./255);


end

figure(fig1);
subplot(2,1,1);
ylabel('V (ms-1)', 'FontSize', 16);
set(gca, 'FontSize', 16);
% ylim([0 0.3]);
% yticks([0:0.1:0.3]);
% xlim([0 0.32]);
% xline(0.04); xline(0.11);
subplot(2,1,2);
ylabel('r (s-1)', 'FontSize', 16);
xlabel('y (m)', 'FontSize', 16);
set(gca, 'FontSize', 16);
xline(0.04); xline(0.11);
legend(arrayfun(@(x) num2str(x),winds,'UniformOutput',false));