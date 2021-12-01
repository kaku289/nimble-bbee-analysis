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

% to include boxplot2
addpath('./lib/boxplot2-pkg/boxplot2');
addpath('./lib/boxplot2-pkg/minmax');

% % to include loSR Surrey for boxplot
% addpath('./lib/MatlabToolbox');
%% Find r* entry intervals
% Also find videos asociated with each track

% Inputs
close all; clc;
% clear;

inputFile = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/BlindLandingtracks_A3_LDF_rref.mat';
outputFile = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/BlindLandingtracks_A3_LDF_rref_rrefEntry.mat';
outputFile = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/BlindLandingtracks_A3_LDF_rref_rrefEntry_with_rdot.mat';
% load(inputFile);

DirPlots = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/plots/BlindLandingTracks/rref_entry';
DirPlots = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/plots/BlindLandingTracks/rref_entry_rdotsim';
% DirPlots = '/media/reken001/Disk_08_backup/wind_intensity_experiments/postprocessing/plots/BlindLandingTracks/rref_estimate_3dspeed';
delPreviousPlots = false; % BE CAREFUL - when set to true, all previously saved plots are deleted
savePlots = false;
savePDFs = false;

winds = unique([treatments.wind]);
behaviour = {'rising','constant','sleeping'};

factors = [0.25:0.25:2.5];
% factors = [1];

winds_vel = [0, 0.28, 0.98, 1.77, 2.54, 3.41];
% interval_for_rdot_estimate = [0.2 0.8]; % in percentage of rref
for ct_wind = 1:length(winds)
    for ct_behaviour = 2%1:length(behaviour)
        
        % Delete previous plots
        if delPreviousPlots && savePlots && exist(fullfile(DirPlots, [num2str(winds(ct_wind)) '_' behaviour{ct_behaviour}]), 'dir')
            rmdir(fullfile(DirPlots, [num2str(winds(ct_wind)) '_' behaviour{ct_behaviour}]),'s');
        end
        
        % Create sub-directory for each treatment if it doesn't exist
        if savePlots && ~exist(fullfile(DirPlots, [num2str(winds(ct_wind)) '_' behaviour{ct_behaviour}]), 'dir')
            mkdir(fullfile(DirPlots, [num2str(winds(ct_wind)) '_' behaviour{ct_behaviour}]));
        end
        DirPlots_treatment = fullfile(DirPlots, [num2str(winds(ct_wind)) '_' behaviour{ct_behaviour}]);
        
        % Selecting relevant treatments
        if strcmpi(behaviour{ct_behaviour}, 'rising')
            relevantTreatments = treatments(rem(1:length(treatments), 8)==1);
        elseif strcmpi(behaviour{ct_behaviour}, 'constant')
            relevantTreatments = treatments( [treatments.wind] == winds(ct_wind) & ...
                rem(1:length(treatments), 8)>1 & ...
                rem(1:length(treatments), 8)<8);
        elseif strcmpi(behaviour{ct_behaviour}, 'sleeping')
            relevantTreatments = treatments(rem(1:length(treatments), 8)==0);
        else
            error('What other treatments did you perform dude?')
        end
        
        hasUniformHwData = arrayfun(@(x) x.hwData.hasUniformHwData,relevantTreatments);
        relevantTreatments = relevantTreatments(hasUniformHwData);
        
        
        for ct_treatment=1:length(relevantTreatments)
            treatment = relevantTreatments(ct_treatment);
            
%             videoTimes = [[treatment.videosInfo.startTime]' [treatment.videosInfo.endTime]'];
            if ~isempty(treatment.landingTracks)
                disp(['Into, day: ' num2str(treatment.datenum) ...
                    ', wind: ' num2str(treatment.wind) ...
                    ', pattern: ' treatment.pattern ...
                    ', wind: ' treatment.wind ...
                    ', behaviour: ' behaviour{ct_behaviour}]);
                
                for ct_track=1:length(treatment.landingTracks) % for each landing track
                    track = treatment.landingTracks(ct_track);
                    for ct_excerpt=1:length(track.state_LDF) % for each track excerpt
                        excerpt = track.state_LDF(ct_excerpt);
                        
                        excerpt.find_rrefEntry();
                        excerpt.find_rdot_estimate_in_rrefEntry([winds_vel(winds(ct_wind)) 0 0]);
                              
                        if savePlots
                            for ct_factor=1:length(factors)
                                plotHandles = excerpt.plot_rdotSimulation_with_actualdata(factors(ct_factor));
%                                       plotHandles = excerpt.plot_rrefsEntry_with_rdotestimate(factors(ct_factor));

%                                       plotHandles = excerpt.plot_rrefsEntry_with_time(factors(ct_factor));
%                                       plotHandles = excerpt.plot_rrefsEntry(factors(ct_factor));
%                                       plotHandles = excerpt.plot_rrefs_with3dspeed(factors(ct_factor));

                                
                                if ~isempty(plotHandles)
                                    
                                    % Resizing the figures
                                    for i=1:length(plotHandles)
                                        % for plot_rdotSimulation_with_actualdata
                                        plotHandles(i).Position(3) = 885;
                                        plotHandles(i).Position(4) = 820;
                                              
                                        % for rest
%                                       plotHandles(i).Position(3) = 680;
%                                       plotHandles(i).Position(4) = 545;
                                        
                                        if i==1
                                            figureName = ['fac_' num2str(factors(ct_factor),'%0.2f') '_' ...
                                                num2str(treatment.datenum) '_' ...
                                                num2str(treatment.startTime) '_' num2str(treatment.endTime) ...
                                                '_obj' num2str(track.obj_id) '_track' ...
                                                num2str(ct_track) '_excerpt' num2str(ct_excerpt) ...
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
                        
                    end % for excerpt
                end % for track
            end
        end % for treatment
        
    end
end
keyboard
save(outputFile, 'treatments', '-v7.3', '-nocompression');
keyboard;

%% Plot Rsquare and difference between y_distance travelled
% if isunix
%     outputFile = '/media/reken001/Disk_08_backup/wind_intensity_experiments/postprocessing/BlindLandingtracks_A2_rrefEntry_with_rdot.mat';
% elseif ispc
%     outputFile = 'D:/wind_intensity_experiments/postprocessing/BlindLandingtracks_A2_rrefEntry_with_rdot.mat';
% end
% load(outputFile);

labels = cell(0, 1); % for x-axis of boxplots
wind_label = {'0', '1', '2', '3', '4', '5'}; 
Rsquare = cell(0,1);

ymean = cell(0,1);
rdot = cell(0,1);
rref = cell(0,1);
isRise = cell(0,1);
delta_y_analytical = cell(0,1);
delta_y_actual = cell(0,1);
diff_delta_y = cell(0,1);

factors = 1;

winds = unique([treatments.wind]);
behaviour = {'rising','constant','sleeping'};
ct_behaviour = 2;
for ct_fac=1:length(factors) 
    for ct_wind = 1:length(winds)
        labels{ct_wind} = [wind_label{ct_wind}];
            
        % Selecting relevant treatments
        if strcmpi(behaviour{ct_behaviour}, 'rising')
            relevantTreatments = treatments(rem(1:length(treatments), 8)==1);
        elseif strcmpi(behaviour{ct_behaviour}, 'constant')
            relevantTreatments = treatments( [treatments.wind] == winds(ct_wind) & ...
                rem(1:length(treatments), 8)>1 & ...
                rem(1:length(treatments), 8)<8);
        elseif strcmpi(behaviour{ct_behaviour}, 'sleeping')
            relevantTreatments = treatments(rem(1:length(treatments), 8)==0);
        else
            error('What other treatments did you perform dude?')
        end
        
        hasUniformHwData = arrayfun(@(x) x.hwData.hasUniformHwData,relevantTreatments);
        relevantTreatments = relevantTreatments(hasUniformHwData);
        
        landingTracks = [relevantTreatments.landingTracks];
        state_LDF = [landingTracks.state_LDF];
        rrefEntrySegments = [state_LDF.rrefEntrySegments];
        rrefEntrySegments = rrefEntrySegments(abs([rrefEntrySegments.factor]-factors(ct_fac))<1e-6);
        
        R2 = arrayfun(@(x) x.Rsquare_rvst, rrefEntrySegments, 'UniformOutput', false);
        Rsquare{ct_wind} = vertcat(R2{:});
        
        
        ymean_dum = arrayfun(@(x) x.ymean_for_rdot, rrefEntrySegments, 'UniformOutput', false);
        ymean{ct_wind} = vertcat(ymean_dum{:});
        
        rdot_dum = arrayfun(@(x) x.slope_rvst, rrefEntrySegments, 'UniformOutput', false);
        rdot{ct_wind} = vertcat(rdot_dum{:});
        
        rref_dum = arrayfun(@(x) x.rref, rrefEntrySegments, 'UniformOutput', false);
        rref{ct_wind} = vertcat(rref_dum{:});
        
        isRise_dum = arrayfun(@(x) x.isRise, rrefEntrySegments, 'UniformOutput', false);
        isRise{ct_wind} = vertcat(isRise_dum{:});
        
        delta_y_analytical_dum = arrayfun(@(x) x.delta_y_analytical, rrefEntrySegments, 'UniformOutput', false);
        delta_y_analytical{ct_wind} = vertcat(delta_y_analytical_dum{:});
        
        delta_y_actual_dum = arrayfun(@(x) x.delta_y_actual, rrefEntrySegments, 'UniformOutput', false);
        delta_y_actual{ct_wind} = vertcat(delta_y_actual_dum{:});
        
        diff_delta_y_dum{ct_wind} = delta_y_analytical{ct_wind} - delta_y_actual{ct_wind};
        
    end
    %     figure;
    %     subplot(1,2,1);
    createBoxPlot({Rsquare{:}}', {labels{:}}', 'R2')
    title(['Factor = ' num2str(factors(ct_fac))], 'FontSize', 18);
    ylim([0.5 1]);
    yticks([0.5:0.1:1]);
    
%     figure;
%     subplot(1,2,2);
    createBoxPlot({diff_delta_y_dum{:}}', {labels{:}}', 'Diff delta_y')
    title(['Factor = ' num2str(factors(ct_fac))], 'FontSize', 18);
end

% For article
quantile(vertcat(Rsquare{:}), [0.25 0.5 0.75])
quantile(vertcat(diff_delta_y_dum{:}), [0.25 0.5 0.75])
% Plot rdot vs ymean
figure;
plot(vertcat(ymean{:}), vertcat(rdot{:}),'o');


figure;
plot(vertcat(rref{:}), vertcat(rdot{:}),'o');


%% Extract rrefEntry data for statistical analysis
clc; close all;

winds = unique([treatments.wind]);
behaviour = {'rising','constant','sleeping'};

data = data4rrefEntry.empty;
hasTakeoff = []; % whether or not the landing track has takeoff dynamics in it as defined in filteredState_BlindLandingtrack.hasTakeoff(..) method

for ct_wind = 1:length(winds)
    for ct_behaviour = 2%1:length(behaviour)
        
        % Selecting relevant treatments
        if strcmpi(behaviour{ct_behaviour}, 'rising')
            relevantTreatments = treatments(rem(1:length(treatments), 8)==1);
        elseif strcmpi(behaviour{ct_behaviour}, 'constant')
            relevantTreatments = treatments( [treatments.wind] == winds(ct_wind) & ...
                rem(1:length(treatments), 8)>1 & ...
                rem(1:length(treatments), 8)<8);
        elseif strcmpi(behaviour{ct_behaviour}, 'sleeping')
            relevantTreatments = treatments(rem(1:length(treatments), 8)==0);
        else
            error('What other treatments did you perform dude?')
        end
        
        hasUniformHwData = arrayfun(@(x) x.hwData.hasUniformHwData,relevantTreatments);
        relevantTreatments = relevantTreatments(hasUniformHwData);
        
        
        for ct_treatment=1:length(relevantTreatments) % for each relevant treatment
            treatment = relevantTreatments(ct_treatment);
            
            arrayfun(@(x) x.setLandingSide(),[treatment.landingTracks.state_LDF]'); % to store landing side in the rrefSegments
            
            data1 = arrayfun(@(x) x.rrefEntrySegments,[treatment.landingTracks.state_LDF]','UniformOutput',false);
            data1 = horzcat(data1{:});
            
            hastakeoff_pertreatment = arrayfun(@(x) x.hasTakeoff(treatment.landingDiscs)*ones(1,length(x.rrefEntrySegments)),[treatment.landingTracks.state_LDF]','UniformOutput',false);
            hastakeoff_pertreatment = horzcat(hastakeoff_pertreatment{:});
            
            % Discard empty intervals
            indices = arrayfun(@(x) ~isempty(x.intervals), data1);
            
            % Save additional data
            data1 = data1(indices);
            [data1.wind] = deal(ct_wind);
            [data1.day] = deal(treatment.datenum);
            [data1.time] = deal(treatment.startTime);
            
            data = [data data1];
            hasTakeoff = [hasTakeoff hastakeoff_pertreatment(indices)];
            %                 size(data)
            %                 keyboard;
        end
    end
end

%% Write file for statistical analysis in R
% This file contains data for all the factors
clc;
writeFile = true;
if isunix
    r_file = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/data_all_rdot_Rstudio.txt';
elseif ispc
    r_file = 'D:/steady_wind_experiments/postprocessing/data_all_rdot_Rstudio.txt';
end
factors = [0.25:0.25:2.5];
data_write = [];
for ct_factor=1:length(factors)
    factor = factors(ct_factor);
    
    data_fac = data(abs([data.factor]-factor)<1e-6)';
    hasTakeoff_fac = hasTakeoff(abs([data.factor]-factor)<1e-6)';
    N = length(data_fac);
    
    % Create nominal and ordinal variables
    approach = arrayfun(@(i) i*ones(size(data_fac(i).intervals,1),1),1:N,'UniformOutput',false);
    side = arrayfun(@(i) data_fac(i).side*ones(size(data_fac(i).intervals,1),1),1:N,'UniformOutput',false);
    wind = arrayfun(@(i) data_fac(i).wind*ones(size(data_fac(i).intervals,1),1),1:N,'UniformOutput',false);
    time = arrayfun(@(i) data_fac(i).time*ones(size(data_fac(i).intervals,1),1),1:N,'UniformOutput',false);
    day = arrayfun(@(i) data_fac(i).day*ones(size(data_fac(i).intervals,1),1),1:N,'UniformOutput',false);
    
    hasTakeoff_fac = arrayfun(@(i) hasTakeoff_fac(i)*ones(size(data_fac(i).intervals,1),1),1:N,'UniformOutput',false);
    
    % create other variables
    y = -vertcat(data_fac.ymean_for_rdot);
    r = -vertcat(data_fac.rref);
    rdot = -vertcat(data_fac.slope_rvst);
    isRise = vertcat(data_fac.isRise);
    delta_r = vertcat(data_fac.delta_r);
    
    yStart = -vertcat(data_fac.yEntryStart);
    vStart = -vertcat(data_fac.vEntryStart);
    rStart = -vertcat(data_fac.rEntryStart);
    delta_V = vertcat(data_fac.delta_Ventry);
    delta_t = vertcat(data_fac.delta_tentry);
    amean = vertcat(data_fac.amean_entry);
    
    mean_Ua = vertcat(data_fac.mean_Ua);
    
    data_write = [data_write; ...
        vertcat(approach{:}) vertcat(side{:}) vertcat(wind{:}) ...
        vertcat(time{:}) vertcat(day{:}) y r rdot factor*ones(size(r,1),1) isRise delta_r vertcat(hasTakeoff_fac{:}) yStart vStart rStart delta_V delta_t amean mean_Ua];
end

if writeFile
    T = array2table(data_write, ...
        'VariableNames',{'approach','landingSide','wind','time','day','y','rref','rdot','threshold', 'isRise', 'delta_r', 'hasTakeoff', 'ystart', 'vstart', 'rstart', 'deltaV', 'deltaT', 'Amean','mean_Ua'});
    writetable(T,r_file);
end

%% Acceleration vs airspeed plot for article-3 discussion
% Data is obtained from R script: 0.5821531 0.6591997 1.1471 1.862773 2.572135 3.448812
Amean = [1.327746 1.394981 1.470067 1.576800 1.709370 1.967488];
mean_Ua = [0.5821531 0.6591997 1.1471 1.862773 2.572135 3.448812];
theUltimateFIT = fitlm(mean_Ua, Amean);
plot(theUltimateFIT)
xlabel('median airspeed in each wind condition, Ua (ms-1)', 'FontSize', 16);
ylabel('mean acceleration (ms-2)', 'FontSize', 16);
set(gca, 'FontSize', 16);
%% Plot histogram of rdote
chosen_fac = 1;
dayCol = 5;
rdotCol = 8;
factorCol = 9;
isriseCol = 10;
deltaVCol = 16;
deltatCol = 17;
% ameanCol = 18;
ameanCol = 20;

data_write(:,20) = data_write(:,deltaVCol)./data_write(:,deltatCol);
dvdtCol = size(data_write,2);
% data_ss = data_write(abs(data_write(:,factorCol) - chosen_fac) < 1e-6, :);
data_ss = data_write(data_write(:,isriseCol) == 1 & abs(data_write(:,factorCol) - chosen_fac) < 1e-6, :);

figure;
% histogram(data_ss(:,rdotCol));
histfit(data_ss(:,rdotCol),[],'Gamma')
xlim([0 50]);
% histogram(-vertcat(data.rmean), [0:0.5:8]);
% histogram(-vertcat(data.rmean), [0:0.5:9.5]);

xlabel('Estimated optical-expansion-rate, rdote (s-2)', 'FontSize', 16);
ylabel('Occurences', 'FontSize', 16);
set(gca, 'FontSize', 16);
pd = fitdist(data_ss(:,rdotCol),'Gamma');
median(data_ss(:,rdotCol))

%% Plot statistical model from R
close all;
% R model
approachCol = 1;
landingSideCol = 2;
windCol = 3;
timeCol = 4;
dayCol = 5;
yCol = 13;
rrefCol = 7;
rdotCol = 8;
factorCol = 9;
isriseCol = 10;
deltarCol = 11;
deltaVCol = 16;
deltatCol = 17;
hasTakeoffCol = 12;
% ameanCol = 17;
% data_write(:,end+1) = data_write(:,rrefCol)-data_write(:,deltarCol);
% r0Col = size(data_write,2);

% rdotCol = ameanCol; %dvdtCol;
% get subset of data

chosen_fac = 1;
data_ss = data_write(abs(data_write(:,factorCol) - chosen_fac) < 1e-6, :);
data_ss = data_write(data_write(:,isriseCol) == 1 & abs(data_write(:,factorCol) - chosen_fac) < 1e-6, :);
data_ss = data_write(data_write(:,isriseCol) == 1 & data_write(:,hasTakeoffCol) == 0 & abs(data_write(:,factorCol) - chosen_fac) < 1e-6, :);
% data_ss = data_write(data_write(:,isriseCol) == 1 & data_write(:,hasTakeoffCol) == 1 & abs(data_write(:,factorCol) - chosen_fac) < 1e-6 & data_write(:,dayCol) <= 20190710, :);
% data_ss = data_write(data_write(:,isriseCol) == 1 & data_write(:,hasTakeoffCol) == 0 & data_write(:,ameanCol) >0 & abs(data_write(:,factorCol) - chosen_fac) < 1e-6 & data_write(:,dayCol) <= 20190710, :);

% Factor = 1, with ystart as y and inclduing hasTakeoff as a factor
% From free-fwind landing
coeffs = [1.46603342 -0.28862172  0.11087894  0.23470481  0.03418062  0.24288514 -0.36214815 -0.06926833]; % [intercept logy wind4 wind6 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
% % From take-off landing
% coeffs = [1.46603342 -0.28862172  0.11087894  0.23470481-0.03537827  0.03418062  0.24288514 -0.36214815 -0.06926833]; % [intercept logy wind4 wind6 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]

chosen_winds = [1 4 6];

% Create basic plots
figure;
subplot(1,2,1)
plot(log(data_ss(:, yCol)), log(data_ss(:,rdotCol)),'.');
set(gca, 'FontSize', 14);
ylabel('log(rdot) (s-2)', 'FontSize', 14);
xlabel('log(y) (m)', 'FontSize', 14);
subplot(1,2,2)
plot((data_ss(:, yCol)), (data_ss(:,rdotCol)),'.');
set(gca, 'FontSize', 14);
ylabel('rdot (s-2)', 'FontSize', 14);
xlabel('y (m)', 'FontSize', 14);

figure;
subplot(1,2,1)
plot(log(data_ss(:, deltarCol)), log(data_ss(:,rdotCol)),'.');
set(gca, 'FontSize', 14);
ylabel('log(rdot) (s-2)', 'FontSize', 14);
xlabel('log(deltar) (s-1)', 'FontSize', 14);
subplot(1,2,2)
plot((data_ss(:, deltarCol)), (data_ss(:,rdotCol)),'.');
set(gca, 'FontSize', 14);
ylabel('rdot (s-2)', 'FontSize', 14);
xlabel('deltar (s-1)', 'FontSize', 14);

figure;
subplot(1,2,1)
plot(log(data_ss(:, rrefCol)), log(data_ss(:,rdotCol)),'.');
set(gca, 'FontSize', 14);
ylabel('log(rdot) (s-2)', 'FontSize', 14);
xlabel('log(rref) (s-1)', 'FontSize', 14);
subplot(1,2,2)
plot((data_ss(:, rrefCol)), (data_ss(:,rdotCol)),'.');
set(gca, 'FontSize', 14);
ylabel('rdot (s-2)', 'FontSize', 14);
xlabel('rref (s-1)', 'FontSize', 14);

% figure;
% subplot(1,2,1)
% plot(log(data_ss(:, r0Col)), log(data_ss(:,rdotCol)),'.');
% set(gca, 'FontSize', 14);
% ylabel('log(rdot) (s-2)', 'FontSize', 14);
% xlabel('log(r0) (s-1)', 'FontSize', 14);
% subplot(1,2,2)
% plot((data_ss(:, r0Col)), (data_ss(:,rdotCol)),'.');
% set(gca, 'FontSize', 14);
% ylabel('rdot (s-2)', 'FontSize', 14);
% xlabel('r0 (s-1)', 'FontSize', 14);

%% Plot rdot as a function of y, delta_r and "fixed" rref (FINAL - TO BE USED)
% One panel for the main figures (y on continuous scale)

% deltar as x1, y as x2 and rref as x3
close all;
rdot = data_ss(:,rdotCol); y = data_ss(:,yCol); wind = data_ss(:,windCol); deltar = data_ss(:,deltarCol); rref = data_ss(:,rrefCol);

% values for which lines will be plotted
% y_percentile_5_50_95 = quantile(y, [0.05, 0.5, 0.95]);
% rref_percentile_5_50_95 = quantile(rref, [0.05, 0.5, 0.95]);
% 
% Untransformed domain
res = rdot; x1 = deltar; x2 = y; x3 = rref; x4 = wind;

x3_centers = quantile(x3,[0.15 0.5 0.85]);% quantile(x3,[0.25 0.5 0.75]); % values for which rref data is binned
binwidth = 1; % data within [x3_center(1)-binwidth/2 x3_center(1)+binwidth/2] is considered to be at x3_center(1)

x2_forplot = [];
for ct=1:3 % for each wind condition
    for ct1=1:length(x3_centers)  % for each x3
        dummy = x2(data_ss(:,windCol) == ct & (data_ss(:,rrefCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,rrefCol) <= x3_centers(ct1)+binwidth/2));
        x2_forplot = [x2_forplot; dummy]; 
    end
end
% x2range = [min(x2_forplot) max(x2_forplot)];
% x2range = [round(x2range(1)*10)/10 round(max(x2range(2))*10)/10];
x2range = [0.05 0.45]; % y range

% cmap=brewermap(round(diff(x2range)/0.01),'Set1'); 
% cmap = jet(round(diff(x2range)/0.01));  
cmap = copper(round(diff(x2range)/0.01));  

% cmap = [0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.4660, 0.6740, 0.1880];
% % % % % % cmap = [253,208,162; 253,174,107; 253,141,60; 241,105,19; 217,72,1; 166,54,3]./255;
edges = round(linspace(x2range(1),x2range(2),size(cmap,1)+1),2); % # edges = # color bins + 1

% For plotting continuous lines
% x2quantiles = quantile(x2_forplot,[0.05 0.5 0.95]) 
x2quantiles = quantile(data_ss(:,yCol),[0.05 0.5 0.95]); % x2 values for which lines are plotted
x2quantiles = [0.1 0.2 0.3 0.4]; %quantile(data_ss(:,yCol),[0.05 0.5 0.95]) % x2 values for which lines are plotted
% deltar_values_cont = 0.2:0.01:4; % deltar values for which lines are plotted

figure2d = figure; hold on;
colormap(cmap);

figure2d_log = figure; hold on;
colormap(cmap);

for dum=3%1:3 % for each wind condition
    ct = chosen_winds(dum);
    for ct1=2%1:length(x3_centers)  % for each x3
        figure(figure2d); hold on;
%         subplot(3, length(x3_centers), (ct-1)*length(x3_centers)+ct1); hold on;
        data_ss2 = data_ss(data_ss(:,windCol) == ct & (data_ss(:,rrefCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,rrefCol) <= x3_centers(ct1)+binwidth/2),:);
        y_ss2 = data_ss2(:,yCol);
%         deltar_ss2 = data_ss2(:,deltarCol);
%         [min(deltar_ss2) max(deltar_ss2)]
        [data_cmap, ~] = discretize(y_ss2, edges);
        data_cmap(y_ss2<=edges(1)) = 1; data_cmap(y_ss2>=edges(end)) = size(cmap,1);
        scatter(data_ss2(:,deltarCol), data_ss2(:,rdotCol), 10, cmap(data_cmap',:),'filled','o');
        
        % Plot continous lines
%         x2quantiles = quantile(deltar_ss2,[0.05 0.5 0.95]) % deltar values for which lines are plotted
        deltar_values_cont = min(data_ss2(:,deltarCol)):0.01:max(data_ss2(:,deltarCol)); % deltar values for which lines are plotted
        dum = [];
        for ct2=1:length(x2quantiles)
            rdot_response = coeffs(1) + coeffs(2)*log(x2quantiles(ct2)) + coeffs(5)*log(deltar_values_cont) + ...
                coeffs(6)*log(x3_centers(ct1)) + coeffs(7)*log(deltar_values_cont)*log(x2quantiles(ct2)) + ...
                + coeffs(8)*log(x3_centers(ct1))*log(deltar_values_cont);
            if ct == 4
                % coeffs = [intercept logy wind2 wind3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
                rdot_response = rdot_response + coeffs(3);
            elseif ct == 6
                rdot_response = rdot_response + coeffs(4);
            end
            linecolor = discretize(x2quantiles(ct2), edges);
            dum(ct2) = linecolor;
            plot(deltar_values_cont, exp(rdot_response),'Color',cmap(linecolor,:),'LineWidth', 2);
        end
        display(dum)
        
        set(gca, 'FontSize', 16);
        ylim([0, 35]);
        yticks([0:5:35]);
        if ct1==1
            xlim([0 2.5]);
            xticks([0:0.5:2.5]);
        elseif ct1==2
            xlim([0 3.3]);
            xticks([0:0.5:3]);
        elseif ct1==3
            xlim([0 4]);
            xticks([0:0.5:4]);
        end
%         figure2d.Position = [-1313 46 895 820];

        
        figure(figure2d_log); hold on;
%         subplot(3, length(x3_centers), (ct-1)*length(x3_centers)+ct1); hold on;
        data_ss2 = data_ss(data_ss(:,windCol) == ct & (data_ss(:,rrefCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,rrefCol) <= x3_centers(ct1)+binwidth/2),:);
        y_ss2 = data_ss2(:,yCol);
        [data_cmap, ~] = discretize(y_ss2, edges);
        data_cmap(y_ss2<=edges(1)) = 1; data_cmap(y_ss2>=edges(end)) = size(cmap,1);
        scatter(log(data_ss2(:,deltarCol)), log(data_ss2(:,rdotCol)), 10, cmap(data_cmap',:),'filled','o');
        
        % Plot continous lines
        deltar_values_cont = min(data_ss2(:,deltarCol)):0.01:max(data_ss2(:,deltarCol)); % deltar values for which lines are plotted
        for ct2=1:length(x2quantiles)
            rdot_response = coeffs(1) + coeffs(2)*log(x2quantiles(ct2)) + coeffs(5)*log(deltar_values_cont) + ...
                coeffs(6)*log(x3_centers(ct1)) + coeffs(7)*log(deltar_values_cont)*log(x2quantiles(ct2)) + ...
                + coeffs(8)*log(x3_centers(ct1))*log(deltar_values_cont);
            if ct == 2
                % coeffs = [intercept logy wind2 wind3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
                rdot_response = rdot_response + coeffs(3);
            elseif ct == 3
                rdot_response = rdot_response + coeffs(4);
            end
            linecolor = discretize(x2quantiles(ct2), edges);
            plot(log(deltar_values_cont), rdot_response,'Color',cmap(linecolor,:),'LineWidth', 2);
        end
        
        set(gca, 'FontSize', 16);
        ylim([1, 4]);
        yticks([1:1:4]);
        xlim([-1.5 1.5]);
        xticks([-1.5:0.5:1.5]);
%         figure2d_log.Position = [-1313 46 895 820];

    end
end
% figure(figure2d);
% cmap_bar = colorbar('eastoutside');
% caxis(edges([1 end]));
% cmap_bar.Ticks = [edges(1) edges(end)];
% cmap_bar.Label.String = 'delta_r';
% cmap_bar.Label.FontSize = 16;


figure;
colormap(cmap);
cmap_bar = colorbar('eastoutside');
caxis(edges([1 end]));
cmap_bar.Ticks = [edges(1) x2quantiles edges(end)];
cmap_bar.Label.String = 'y';
cmap_bar.Label.FontSize = 16;
ylabel('rdot (s-2)', 'FontSize', 14);
xlabel('deltar (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);

%% Plot rdot as a function of delta_r, rref and "fixed" y (TO BE USED)
% deltar as x1, rref as x2 and y as x3
close all;
rdot = data_ss(:,rdotCol); y = data_ss(:,yCol); wind = data_ss(:,windCol); deltar = data_ss(:,deltarCol); rref = data_ss(:,rrefCol);

% values for which lines will be plotted
y_percentile_5_50_95 = quantile(y, [0.05, 0.5, 0.95]);
rref_percentile_5_50_95 = quantile(rref, [0.05, 0.5, 0.95]);

% Untransformed domain
res = rdot; x1 = deltar; x2 = rref; x3 = y; x4 = wind;

x3_centers = quantile(x3,[0.25 0.5 0.75]); % values for which rref data is binned
binwidth = 0.05; % data within [x3_center(1)-binwidth/2 x3_center(1)+binwidth/2] is considered to be at x3_center(1)
% x3_centers = [0.1 0.2 0.3];
% binwidth = 0.1; % data within [x3_center(1)-binwidth/2 x3_center(1)+binwidth/2] is considered to be at x3_center(1)

x2_forplot = [];
for dum=1:3 % for each wind condition
    ct = chosen_winds(dum);
    for ct1=1:length(x3_centers)  % for each x3
        dummy = x2(data_ss(:,windCol) == ct & (data_ss(:,yCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,yCol) <= x3_centers(ct1)+binwidth/2));
        x2_forplot = [x2_forplot; dummy]; 
    end
end
% x2range = rref_percentile_5_50_95([1 end]);
% x2range = [round(x2range(1)*10)/10 round(max(x2range(2))*10)/10];
x2range = [1 5.5];

% cmap=brewermap(round(diff(x2range)/0.5),'Set1'); 
% cmap = jet(round(diff(x2range)/0.1));  
% cmap = [0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.4660, 0.6740, 0.1880];
cmap = copper(round(diff(x2range)/0.15));  

edges = round(linspace(x2range(1),x2range(2),size(cmap,1)+1),2); % # edges = # color bins + 1

% For plotting continuous lines
% x2quantiles = quantile(x2_forplot,[0.05 0.5 0.95]) 
x2quantiles = quantile(data_ss(:,rrefCol),[0.05 0.5 0.95]) % x2 values for which lines are plotted
% deltar_values_cont = 0.2:0.01:4; % deltar values for which lines are plotted

figure2d = figure; hold on;
colormap(cmap);

figure2d_log = figure; hold on;
colormap(cmap);
for dummy_ct=1%1:3 % for each wind condition
    ct=chosen_winds(dummy_ct);
    for ct1=2%1:length(x3_centers)  % for each x3
        figure(figure2d);
%         subplot(3, length(x3_centers), (dummy_ct-1)*length(x3_centers)+ct1); hold on;
        data_ss2 = data_ss(data_ss(:,windCol) == ct & (data_ss(:,yCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,yCol) <= x3_centers(ct1)+binwidth/2),:);
        rref_ss2 = data_ss2(:,rrefCol);
%         deltar_ss2 = data_ss2(:,deltarCol);
%         [min(deltar_ss2) max(deltar_ss2)]
        [data_cmap, ~] = discretize(rref_ss2, edges);
        data_cmap(rref_ss2<=edges(1)) = 1; data_cmap(rref_ss2>=edges(end)) = size(cmap,1);
        scatter(data_ss2(:,deltarCol), data_ss2(:,rdotCol), 10, cmap(data_cmap',:),'filled','o');
        
        % Plot continous lines
%         x2quantiles = quantile(deltar_ss2,[0.05 0.5 0.95]) % deltar values for which lines are plotted
        deltar_values_cont = min(data_ss2(:,deltarCol)):0.01:max(data_ss2(:,deltarCol)); % deltar values for which lines are plotted
        dum = [];
        for ct2=1:length(x2quantiles)
            rdot_response = coeffs(1) + coeffs(2)*log(x3_centers(ct1)) + coeffs(5)*log(deltar_values_cont) + ...
                coeffs(6)*log(x2quantiles(ct2)) + coeffs(7)*log(deltar_values_cont)*log(x3_centers(ct1)) + ...
                + coeffs(8)*log(x2quantiles(ct2))*log(deltar_values_cont);
            if dummy_ct == 2
                % coeffs = [intercept logy wind2 wind3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
                rdot_response = rdot_response + coeffs(3);
            elseif dummy_ct == 3
                rdot_response = rdot_response + coeffs(4);
            end
            linecolor = discretize(x2quantiles(ct2), edges);
            dum(ct2) = linecolor;
            plot(deltar_values_cont, exp(rdot_response),'Color',cmap(linecolor,:),'LineWidth', 2);
        end
        display(dum)
        
        set(gca, 'FontSize', 16);
        ylim([0, 35]);
        yticks([0:5:35]);
        if ct1==1
            xlim([0 5]);
            xticks([0:1:5]);
        elseif ct1==2
            xlim([0 6]);
            xticks([0:1:6]);
        elseif ct1==3
            xlim([0 6]);
            xticks([0:1:6]);
        end
%         figure2d.Position = [-1313 46 895 820];
        
        figure(figure2d_log);
%         subplot(3, length(x3_centers), (dummy_ct-1)*length(x3_centers)+ct1); hold on;
        data_ss2 = data_ss(data_ss(:,windCol) == ct & (data_ss(:,yCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,yCol) <= x3_centers(ct1)+binwidth/2),:);
        rref_ss2 = data_ss2(:,rrefCol);
        [data_cmap, ~] = discretize(rref_ss2, edges);
        data_cmap(rref_ss2<=edges(1)) = 1; data_cmap(rref_ss2>=edges(end)) = size(cmap,1);
        scatter(log(data_ss2(:,deltarCol)), log(data_ss2(:,rdotCol)), 10, cmap(data_cmap',:),'filled','o');
        
        % Plot continous lines
        deltar_values_cont = min(data_ss2(:,deltarCol)):0.01:max(data_ss2(:,deltarCol)); % deltar values for which lines are plotted
        for ct2=1:length(x2quantiles)
            rdot_response = coeffs(1) + coeffs(2)*log(x3_centers(ct1)) + coeffs(5)*log(deltar_values_cont) + ...
                coeffs(6)*log(x2quantiles(ct2)) + coeffs(7)*log(deltar_values_cont)*log(x3_centers(ct1)) + ...
                + coeffs(8)*log(x2quantiles(ct2))*log(deltar_values_cont);
            if dummy_ct == 2
                % coeffs = [intercept logy wind2 wind3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
                rdot_response = rdot_response + coeffs(3);
            elseif dummy_ct == 3
                rdot_response = rdot_response + coeffs(4);
            end
            linecolor = discretize(x2quantiles(ct2), edges);
            plot(log(deltar_values_cont), rdot_response,'Color',cmap(linecolor,:),'LineWidth', 2);
        end
        
        set(gca, 'FontSize', 16);
        ylim([0, 4]);
        yticks([0:1:4]);
        xlim([-1.5 2]);
        xticks([-1.5:0.5:2]);
%         figure2d_log.Position = [-1313 46 895 820];
    end
end
% figure(figure2d);
% cmap_bar = colorbar('eastoutside');
% caxis(edges([1 end]));
% cmap_bar.Ticks = [edges(1) edges(end)];
% cmap_bar.Label.String = 'delta_r';
% cmap_bar.Label.FontSize = 16;


figure;
colormap(cmap);
cmap_bar = colorbar('eastoutside');
caxis(edges([1 end]));
cmap_bar.Ticks = [edges(1) x2quantiles edges(end)];
% cmap_bar.Ticks = [edges(1) edges(end)];
cmap_bar.Label.String = 'rref';
cmap_bar.Label.FontSize = 16;
ylabel('rdot (s-2)', 'FontSize', 14);
xlabel('deltar (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);

%% Plot rdot as a function of rref, delta_r and "fixed" y (NOT TO BE USED)
% deltar as x2, rref as x1 and y as x3
% There is some issure - need to be fixed

close all;
rdot = data_ss(:,rdotCol); y = data_ss(:,yCol); wind = data_ss(:,windCol); deltar = data_ss(:,deltarCol); rref = data_ss(:,rrefCol);

% values for which lines will be plotted
y_percentile_5_50_95 = quantile(y, [0.05, 0.5, 0.95]);
rref_percentile_5_50_95 = quantile(rref, [0.05, 0.5, 0.95]);
deltar_percentile_5_50_95 = quantile(deltar, [0.05, 0.5, 0.95]);

% Untransformed domain
res = rdot; x2 = deltar; x1 = rref; x3 = y; x4 = wind;

x3_centers = quantile(x3,[0.25 0.5 0.75]); % values for which y data is binned
binwidth = 0.05; % data within [x3_center(1)-binwidth/2 x3_center(1)+binwidth/2] is considered to be at x3_center(1)
% x3_centers = [0.1 0.2 0.3];
% binwidth = 0.1; % data within [x3_center(1)-binwidth/2 x3_center(1)+binwidth/2] is considered to be at x3_center(1)

x2_forplot = [];
for dum=1:3 % for each wind condition
    ct = chosen_winds(dum);
    for ct1=1:length(x3_centers)  % for each x3
        dummy = x2(data_ss(:,windCol) == ct & (data_ss(:,yCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,yCol) <= x3_centers(ct1)+binwidth/2));
        x2_forplot = [x2_forplot; dummy]; 
    end
end
% x2range = rref_percentile_5_50_95([1 end]);
% x2range = [round(x2range(1)*10)/10 round(max(x2range(2))*10)/10];
x2range = [1 3];

% cmap=brewermap(round(diff(x2range)/0.5),'Set1'); 
% cmap = jet(round(diff(x2range)/0.1));  
% cmap = [0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.4660, 0.6740, 0.1880];
cmap = copper(round(diff(x2range)/0.1));  

edges = round(linspace(x2range(1),x2range(2),size(cmap,1)+1),2); % # edges = # color bins + 1

% For plotting continuous lines
% x2quantiles = quantile(x2_forplot,[0.05 0.5 0.95]) 
% x2quantiles = quantile(data_ss(:,deltarCol),[0.25 0.5 0.75]) % x2 values for which lines are plotted
x2quantiles = [1.5 2 2.5]; % x2 values for which lines are plotted
% deltar_values_cont = 0.2:0.01:4; % deltar values for which lines are plotted

figure2d = figure; hold on;
colormap(cmap);

figure2d_log = figure; hold on;
colormap(cmap);
for dummy_ct=1:3 % for each wind condition
    ct=chosen_winds(dummy_ct);
    for ct1=1:length(x3_centers)  % for each x3
        figure(figure2d);
        subplot(3, length(x3_centers), (dummy_ct-1)*length(x3_centers)+ct1); hold on;
        data_ss2 = data_ss(data_ss(:,windCol) == ct & (data_ss(:,yCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,yCol) <= x3_centers(ct1)+binwidth/2),:);
%         rref_ss2 = data_ss2(:,rrefCol);
        deltar_ss2 = data_ss2(:,deltarCol);
%         [min(deltar_ss2) max(deltar_ss2)]
        [data_cmap, ~] = discretize(deltar_ss2, edges);
        data_cmap(deltar_ss2<=edges(1)) = 1; data_cmap(deltar_ss2>=edges(end)) = size(cmap,1);
        scatter(data_ss2(:,rrefCol), data_ss2(:,rdotCol), 10, cmap(data_cmap',:),'filled','o');
        
        % Plot continous lines
%         x2quantiles = quantile(deltar_ss2,[0.05 0.5 0.95]) % deltar values for which lines are plotted
%         deltar_values_cont = min(data_ss2(:,deltarCol)):0.01:max(data_ss2(:,deltarCol)); % deltar values for which lines are plotted
        rref_values_cont = min(data_ss2(:,rrefCol)):0.01:max(data_ss2(:,rrefCol)); % rref values for which lines are plotted
        dum = [];
        for ct2=1:length(x2quantiles)
            rdot_response = coeffs(1) + coeffs(2)*log(x3_centers(ct1)) + coeffs(5)*log(x2quantiles(ct2) ) + ...
                coeffs(6)*log(rref_values_cont) + coeffs(7)*log(x2quantiles(ct2))*log(x3_centers(ct1)) + ...
                + coeffs(8)*log(x2quantiles(ct2))*log(rref_values_cont);
            if dummy_ct == 2
                % coeffs = [intercept logy wind2 wind3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
                rdot_response = rdot_response + coeffs(3);
            elseif dummy_ct == 3
                rdot_response = rdot_response + coeffs(4);
            end
            linecolor = discretize(x2quantiles(ct2), edges);
            dum(ct2) = linecolor;
            plot(rref_values_cont, exp(rdot_response),'Color',cmap(linecolor,:),'LineWidth', 2);
        end
        display(dum)
        
        set(gca, 'FontSize', 16);
        ylim([0, 35]);
        yticks([0:5:35]);
%         if ct1==1
%             xlim([0 5]);
%             xticks([0:1:5]);
%         elseif ct1==2
%             xlim([0 6]);
%             xticks([0:1:6]);
%         elseif ct1==3
%             xlim([0 6]);
%             xticks([0:1:6]);
%         end
        figure2d.Position = [-1313 46 895 820];
        
        figure(figure2d_log);
        subplot(3, length(x3_centers), (dummy_ct-1)*length(x3_centers)+ct1); hold on;
        data_ss2 = data_ss(data_ss(:,windCol) == ct & (data_ss(:,yCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,yCol) <= x3_centers(ct1)+binwidth/2),:);
        deltar_ss2 = data_ss2(:,deltarCol);
        [data_cmap, ~] = discretize(deltar_ss2, edges);
        data_cmap(deltar_ss2<=edges(1)) = 1; data_cmap(deltar_ss2>=edges(end)) = size(cmap,1);
        scatter(log(data_ss2(:,rrefCol)), log(data_ss2(:,rdotCol)), 10, cmap(data_cmap',:),'filled','o');
        
        % Plot continous lines
        rref_values_cont = min(data_ss2(:,rrefCol)):0.01:max(data_ss2(:,rrefCol)); % rref values for which lines are plotted
        for ct2=1:length(x2quantiles)
            rdot_response = coeffs(1) + coeffs(2)*log(x3_centers(ct1)) + coeffs(5)*log(x2quantiles(ct2) ) + ...
                coeffs(6)*log(rref_values_cont) + coeffs(7)*log(x2quantiles(ct2) )*log(x3_centers(ct1)) + ...
                + coeffs(8)*log(x2quantiles(ct2))*log(rref_values_cont);
            if dummy_ct == 2
                % coeffs = [intercept logy wind2 wind3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
                rdot_response = rdot_response + coeffs(3);
            elseif dummy_ct == 3
                rdot_response = rdot_response + coeffs(4);
            end
            linecolor = discretize(x2quantiles(ct2), edges);
            plot(log(rref_values_cont), rdot_response,'Color',cmap(linecolor,:),'LineWidth', 2);
        end
        
        set(gca, 'FontSize', 16);
        ylim([0, 4]);
        yticks([0:1:4]);
        xlim([-1.5 2]);
        xticks([-1.5:0.5:2]);
        figure2d_log.Position = [-1313 46 895 820];
    end
end
% figure(figure2d);
% cmap_bar = colorbar('eastoutside');
% caxis(edges([1 end]));
% cmap_bar.Ticks = [edges(1) edges(end)];
% cmap_bar.Label.String = 'delta_r';
% cmap_bar.Label.FontSize = 16;


figure;
colormap(cmap);
cmap_bar = colorbar('eastoutside');
caxis(edges([1 end]));
cmap_bar.Ticks = [edges(1) edges(end)];
cmap_bar.Label.String = 'deltar';
cmap_bar.Label.FontSize = 16;
ylabel('rdot (s-2)', 'FontSize', 14);
xlabel('rref (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);

%% Plot rdot as a function of y, delta_r and "fixed" rref (FINAL - TO BE USED) - Supplement figure
% One panel for the main figures (y on continuous scale)

% deltar as x1, y as x2 and rref as x3
close all;
rdot = data_ss(:,rdotCol); y = data_ss(:,yCol); wind = data_ss(:,windCol); deltar = data_ss(:,deltarCol); rref = data_ss(:,rrefCol);

% values for which lines will be plotted
% y_percentile_5_50_95 = quantile(y, [0.05, 0.5, 0.95]);
% rref_percentile_5_50_95 = quantile(rref, [0.05, 0.5, 0.95]);
% 
% Untransformed domain
res = rdot; x1 = deltar; x2 = y; x3 = rref; x4 = wind;

x3_centers = quantile(x3,[0.15 0.5 0.85]);% quantile(x3,[0.25 0.5 0.75]); % values for which rref data is binned
binwidth = 1; % data within [x3_center(1)-binwidth/2 x3_center(1)+binwidth/2] is considered to be at x3_center(1)

x2_forplot = [];
for ct=1:3 % for each wind condition
    for ct1=1:length(x3_centers)  % for each x3
        dummy = x2(data_ss(:,windCol) == ct & (data_ss(:,rrefCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,rrefCol) <= x3_centers(ct1)+binwidth/2));
        x2_forplot = [x2_forplot; dummy]; 
    end
end
% x2range = [min(x2_forplot) max(x2_forplot)];
% x2range = [round(x2range(1)*10)/10 round(max(x2range(2))*10)/10];
x2range = [0.05 0.35]; % y range

% cmap=brewermap(round(diff(x2range)/0.01),'Set1'); 
% cmap = jet(round(diff(x2range)/0.01));  
cmap = copper(round(diff(x2range)/0.01));  

% cmap = [0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.4660, 0.6740, 0.1880];
% % % % % % cmap = [253,208,162; 253,174,107; 253,141,60; 241,105,19; 217,72,1; 166,54,3]./255;
edges = round(linspace(x2range(1),x2range(2),size(cmap,1)+1),2); % # edges = # color bins + 1

% For plotting continuous lines
% x2quantiles = quantile(x2_forplot,[0.05 0.5 0.95]) 
x2quantiles = quantile(data_ss(:,yCol),[0.05 0.5 0.95]); % x2 values for which lines are plotted
x2quantiles = [0.1 0.2 0.3]; %quantile(data_ss(:,yCol),[0.05 0.5 0.95]) % x2 values for which lines are plotted
% deltar_values_cont = 0.2:0.01:4; % deltar values for which lines are plotted

figure2d = figure; hold on;
colormap(cmap);

figure2d_log = figure; hold on;
colormap(cmap);
for dummy_ct=1:3 % for each wind condition
    ct = chosen_winds(dummy_ct);
    for ct1=1:length(x3_centers)  % for each x3
        figure(figure2d); hold on;
        subplot(3, length(x3_centers), (dummy_ct-1)*length(x3_centers)+ct1); hold on;
        data_ss2 = data_ss(data_ss(:,windCol) == ct & (data_ss(:,rrefCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,rrefCol) <= x3_centers(ct1)+binwidth/2),:);
        y_ss2 = data_ss2(:,yCol);
%         deltar_ss2 = data_ss2(:,deltarCol);
%         [min(deltar_ss2) max(deltar_ss2)]
        [data_cmap, ~] = discretize(y_ss2, edges);
        data_cmap(y_ss2<=edges(1)) = 1; data_cmap(y_ss2>=edges(end)) = size(cmap,1);
        scatter(data_ss2(:,deltarCol), data_ss2(:,rdotCol), 10, cmap(data_cmap',:),'filled','o');
        
        % Plot continous lines
%         x2quantiles = quantile(deltar_ss2,[0.05 0.5 0.95]) % deltar values for which lines are plotted
        deltar_values_cont = min(data_ss2(:,deltarCol)):0.01:max(data_ss2(:,deltarCol)); % deltar values for which lines are plotted
        dum = [];
        for ct2=1:length(x2quantiles)
            rdot_response = coeffs(1) + coeffs(2)*log(x2quantiles(ct2)) + coeffs(5)*log(deltar_values_cont) + ...
                coeffs(6)*log(x3_centers(ct1)) + coeffs(7)*log(deltar_values_cont)*log(x2quantiles(ct2)) + ...
                + coeffs(8)*log(x3_centers(ct1))*log(deltar_values_cont);
            if ct == 2
                % coeffs = [intercept logy wind2 wind3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
                rdot_response = rdot_response + coeffs(3);
            elseif ct == 3
                rdot_response = rdot_response + coeffs(4);
            end
            linecolor = discretize(x2quantiles(ct2), edges);
            dum(ct2) = linecolor;
            plot(deltar_values_cont, exp(rdot_response),'Color',cmap(linecolor,:),'LineWidth', 2);
        end
        display(dum)
        
        set(gca, 'FontSize', 16);
        ylim([0, 35]);
        yticks([0:5:35]);
        if ct1==1
            xlim([0 2.5]);
            xticks([0:0.5:2.5]);
        elseif ct1==2
            xlim([0 3.3]);
            xticks([0:0.5:3]);
        elseif ct1==3
            xlim([0 4]);
            xticks([0:0.5:4]);
        end
%         figure2d.Position = [-1313 46 895 820];

        
        figure(figure2d_log); hold on;
        subplot(3, length(x3_centers), (dummy_ct-1)*length(x3_centers)+ct1); hold on;
        data_ss2 = data_ss(data_ss(:,windCol) == ct & (data_ss(:,rrefCol) >= x3_centers(ct1)-binwidth/2 & data_ss(:,rrefCol) <= x3_centers(ct1)+binwidth/2),:);
        y_ss2 = data_ss2(:,yCol);
        [data_cmap, ~] = discretize(y_ss2, edges);
        data_cmap(y_ss2<=edges(1)) = 1; data_cmap(y_ss2>=edges(end)) = size(cmap,1);
        scatter(log(data_ss2(:,deltarCol)), log(data_ss2(:,rdotCol)), 10, cmap(data_cmap',:),'filled','o');
        
        % Plot continous lines
        deltar_values_cont = min(data_ss2(:,deltarCol)):0.01:max(data_ss2(:,deltarCol)); % deltar values for which lines are plotted
        for ct2=1:length(x2quantiles)
            rdot_response = coeffs(1) + coeffs(2)*log(x2quantiles(ct2)) + coeffs(5)*log(deltar_values_cont) + ...
                coeffs(6)*log(x3_centers(ct1)) + coeffs(7)*log(deltar_values_cont)*log(x2quantiles(ct2)) + ...
                + coeffs(8)*log(x3_centers(ct1))*log(deltar_values_cont);
            if ct == 2
                % coeffs = [intercept logy wind2 wind3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
                rdot_response = rdot_response + coeffs(3);
            elseif ct == 3
                rdot_response = rdot_response + coeffs(4);
            end
            linecolor = discretize(x2quantiles(ct2), edges);
            plot(log(deltar_values_cont), rdot_response,'Color',cmap(linecolor,:),'LineWidth', 2);
        end
        
        set(gca, 'FontSize', 16);
        ylim([1, 4]);
        yticks([1:1:4]);
        xlim([-1.5 1.5]);
        xticks([-1.5:0.5:1.5]);
%         figure2d_log.Position = [-1313 46 895 820];

    end
end
% figure(figure2d);
% cmap_bar = colorbar('eastoutside');
% caxis(edges([1 end]));
% cmap_bar.Ticks = [edges(1) edges(end)];
% cmap_bar.Label.String = 'delta_r';
% cmap_bar.Label.FontSize = 16;


figure;
colormap(cmap);
cmap_bar = colorbar('eastoutside');
caxis(edges([1 end]));
cmap_bar.Ticks = [edges(1) x2quantiles edges(end)];
cmap_bar.Label.String = 'y';
cmap_bar.Label.FontSize = 16;
ylabel('rdot (s-2)', 'FontSize', 14);
xlabel('deltar (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);

%% Plot rdot as a function of rref (FINAL TO BE USED) - Supplement figure
% rref as x1
close all;
rdot = data_ss(:,rdotCol); y = data_ss(:,yCol); wind = data_ss(:,windCol); deltar = data_ss(:,deltarCol); rref = data_ss(:,rrefCol);

% values for which lines will be plotted
y_chosen = 0.21; %quantile(y, 0.5); % mean(y)
y_bin = 0.06; %mean(y) std(y)
deltar_chosen = 1.68; quantile(deltar, 0.5); %mean(delta_r) std(delta_r)
deltar_bin = 0.8;

% Untransformed domain
res = rdot; x1 = rref;

% x3_centers = quantile(x3,[0.25 0.5 0.75]); % values for which rref data is binned
% binwidth = 0.05; % data within [x3_center(1)-binwidth/2 x3_center(1)+binwidth/2] is considered to be at x3_center(1)

x1_forplot = [];
for ct=1:3 % for each wind condition
    dummy = x1(data_ss(:,windCol) == ct & (data_ss(:,yCol) >= y_chosen-y_bin/2 & data_ss(:,yCol) <= y_chosen+y_bin/2) & ...
                                           (data_ss(:,deltarCol) >= deltar_chosen-deltar_bin/2 & data_ss(:,deltarCol) <= deltar_chosen+deltar_bin/2));
    x1_forplot = [x1_forplot; dummy];
end
x1range = [min(x1_forplot) max(x1_forplot)];
x1range = [floor(x1range(1)*10)/10 round(max(x1range(2))*10)/10];
rref_values_cont = x1range(1):0.1:x1range(2);
% x2range = [1 5.5];

% cmap=brewermap(round(diff(x2range)/0.5),'Set1'); 
% cmap = jet(round(diff(x2range)/0.1));  
cmap = [0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.4660, 0.6740, 0.1880];
edges = round(linspace(x2range(1),x2range(2),size(cmap,1)+1),2); % # edges = # color bins + 1

% For plotting continuous lines
% x2quantiles = quantile(x2_forplot,[0.05 0.5 0.95]) 
% x2quantiles = quantile(data_ss(:,rrefCol),[0.05 0.5 0.95]) % x2 values for which lines are plotted
 % deltar values for which lines are plotted

figure2d = figure; hold on;
colormap(cmap);

figure2d_log = figure; hold on;
colormap(cmap);
for dummy_ct=1:3 % for each wind condition
        ct = chosen_winds(dummy_ct);
        figure(figure2d); hold on;
        subplot(3, 1, dummy_ct); hold on;
        data_ss2 = data_ss(data_ss(:,windCol) == ct & (data_ss(:,yCol) >= y_chosen-y_bin/2 & data_ss(:,yCol) <= y_chosen+y_bin/2) & ...
                                           (data_ss(:,deltarCol) >= deltar_chosen-deltar_bin/2 & data_ss(:,deltarCol) <= deltar_chosen+deltar_bin/2),:);
%         rref_ss2 = data_ss2(:,rrefCol);
%         [data_cmap, ~] = discretize(rref_ss2, edges);
%         data_cmap(rref_ss2<=edges(1)) = 1; data_cmap(rref_ss2>=edges(end)) = size(cmap,1);
        scatter(data_ss2(:,rrefCol), data_ss2(:,rdotCol), 10,'filled','o');
        
        % Plot continous lines
        
        rdot_response = coeffs(1) + coeffs(2)*log(y_chosen) + coeffs(5)*log(deltar_chosen) + ...
            coeffs(6)*log(rref_values_cont) + coeffs(7)*log(y_chosen)*log(deltar_chosen) + ...
            + coeffs(8)*log(deltar_chosen)*log(rref_values_cont);
        if ct == 2
            % coeffs = [intercept logy wind2 wind3 log(deltar) log(rref) logy:log(deltar) log(deltar):log(rref)]
            rdot_response = rdot_response + coeffs(3);
        elseif ct == 3
            rdot_response = rdot_response + coeffs(4);
        end
%         linecolor = discretize(x2quantiles(ct2), edges);
        plot(rref_values_cont, exp(rdot_response),'Color',[0 0 0],'LineWidth', 2);
        
        set(gca, 'FontSize', 16);
        ylim([5, 20]);
        yticks([5:5:20]);
%         if ct1==1
            xlim([1.5 5]);
            xticks([1:1:5]);
%         elseif ct1==2
%             xlim([0 6]);
%             xticks([0:1:6]);
%         elseif ct1==3
%             xlim([0 6]);
%             xticks([0:1:6]);
%         end
%         figure2d.Position = [326    42   335   954];
        
        figure(figure2d_log); hold on;
        subplot(3, 1, dummy_ct); hold on;
        data_ss2 = data_ss(data_ss(:,windCol) == ct & (data_ss(:,yCol) >= y_chosen-y_bin/2 & data_ss(:,yCol) <= y_chosen+y_bin/2) & ...
                                           (data_ss(:,deltarCol) >= deltar_chosen-deltar_bin/2 & data_ss(:,deltarCol) <= deltar_chosen+deltar_bin/2),:);
        scatter(log(data_ss2(:,rrefCol)), log(data_ss2(:,rdotCol)), 10,'filled','o');
        
        % Plot continous lines
        plot(log(rref_values_cont), rdot_response,'Color',[0 0 0],'LineWidth', 2);
        
        
        set(gca, 'FontSize', 16);
        ylim([1.5 3]);
        yticks([1.5:0.5:3]);
        xlim([0.5 1.75]);
        xticks([0.5:0.5:1.5]);
%         figure2d_log.Position = [667    34   335   961];
end
% figure(figure2d);
% cmap_bar = colorbar('eastoutside');
% caxis(edges([1 end]));
% cmap_bar.Ticks = [edges(1) edges(end)];
% cmap_bar.Label.String = 'delta_r';
% cmap_bar.Label.FontSize = 16;


figure;
% colormap(cmap);
% cmap_bar = colorbar('eastoutside');
% caxis(edges([1 end]));
% cmap_bar.Ticks = [edges(1) edges(end)];
% cmap_bar.Label.String = 'rref';
% cmap_bar.Label.FontSize = 16;
ylabel('rdot (s-2)', 'FontSize', 14);
xlabel('rref (s-1)', 'FontSize', 14);
set(gca, 'FontSize', 16);

%% Functions used in this script
function figHandle = createBoxPlot(variable, labels, yxislabel)
    % for boxplots of a cell array
    col=@(x)reshape(x,numel(x),1);
    boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});

    
    figHandle = figure; hold on;
    figure(figHandle);
    boxplot2(variable, 'Labels', labels, 'OutlierSize', 0.00001);
%     hbox = gca;
    set(gca, 'FontSize', 18); grid on;
%     xlabel('Treatments', 'FontSize', 18);
    ylabel(yxislabel, 'FontSize', 18);
    % ylim([0 1.5]);
    for ct = 1:length(variable)
            x=(ct+(rand(length(variable{ct, 1}),1)-0.5)/4);

            f = scatter(x(:,1),variable{ct, 1},40,'k','filled'); 
            f.MarkerFaceAlpha = 0.5;
    %         keyboard;
    end
end

















