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
% DirPlots = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/plots/BlindLandingTracks/rref_estimate_3dspeed';
delPreviousPlots = false; % BE CAREFUL - when set to true, all previously saved plots are deleted
savePlots = false;
savePDFs = false;

winds = unique([treatments.wind]);
behaviour = {'rising','constant','sleeping'};

factors = [0.25:0.25:2.5];
% factors = [1];

interval_for_rdot_estimate = [0.2 0.8]; % in percentage of rref
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
                    ', light: ' treatment.light ...
                    ', behaviour: ' behaviour{ct_behaviour}]);
                
                for ct_track=1:length(treatment.landingTracks) % for each landing track
                    track = treatment.landingTracks(ct_track);
                    for ct_excerpt=1:length(track.state_LDF) % for each track excerpt
                        excerpt = track.state_LDF(ct_excerpt);
                        
                        excerpt.find_rrefEntry();
                        excerpt.find_rdot_estimate_in_rrefEntry(interval_for_rdot_estimate);
                              
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
%     outputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A2_rrefEntry_with_rdot.mat';
% elseif ispc
%     outputFile = 'D:/light_intensity_experiments/postprocessing/BlindLandingtracks_A2_rrefEntry_with_rdot.mat';
% end
% load(outputFile);

labels = cell(0, 1); % for x-axis of boxplots
wind_label = {'0', '1', '2', '3', '4', '5'}; % + is checkerboard, x is spokes

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
    
    data_write = [data_write; ...
        vertcat(approach{:}) vertcat(side{:}) vertcat(wind{:}) ...
        vertcat(time{:}) vertcat(day{:}) y r rdot factor*ones(size(r,1),1) isRise delta_r vertcat(hasTakeoff_fac{:}) yStart vStart rStart delta_V delta_t amean];
end

if writeFile
    T = array2table(data_write, ...
        'VariableNames',{'approach','landingSide','wind','time','day','y','rref','rdot','threshold', 'isRise', 'delta_r', 'hasTakeoff', 'ystart', 'vstart', 'rstart', 'deltaV', 'deltaT', 'Amean'});
    writetable(T,r_file);
end

%% Plot histogram of rdote
chosen_fac = 1;
dayCol = 5;
rdotCol = 8;
factorCol = 9;
isriseCol = 10;
data_ss = data_write(data_write(:,isriseCol) == 1 & abs(data_write(:,factorCol) - chosen_fac) < 1e-6, :);

figure;
% histogram(data_ss(:,rdotCol));
histfit(data_ss(:,rdotCol),[],'Gamma')
xlim([0 40]);
% histogram(-vertcat(data.rmean), [0:0.5:8]);
% histogram(-vertcat(data.rmean), [0:0.5:9.5]);

xlabel('Estimated optical-expansion-rate, rdote (s-2)', 'FontSize', 16);
ylabel('Occurences', 'FontSize', 16);
set(gca, 'FontSize', 16);
pd = fitdist(data_ss(:,rdotCol),'Gamma');
median(data_ss(:,rdotCol))

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

















