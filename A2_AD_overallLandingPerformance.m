%% This file plots landing performance parameters (time-to-land and distance travelled)

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


%% 

% Extract all 10005 approaches
close all;
clc; 
% clear;
% inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A1_rref.mat';
inputFile = 'D:/light_intensity_experiments/postprocessing/BlindLandingtracks_A1_rref.mat';
% load(inputFile);
treatments = treatments(1:14*8);

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
behaviour = {'rising','constant','sleeping'};
data_all = struct.empty;
for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        for ct_behaviour = 2%1:length(behaviour)
            clear dummy;

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
            landingTracks = [relevantTreatments.landingTracks];
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
            dummy.pattern = (ct_pattern)*ones(1, length(dummy.state_LDF));
            dummy.light = (ct_light)*ones(1, length(dummy.state_LDF));
            dummy.day = days';
            dummy.time = startTimes';
            
            data_all = [data_all; dummy];
  
        end
    end
end

%%
close all;

y_start = [-0.30:0.06:-0.18];
y_end = -0.02;
dt = cell.empty;
ds = cell.empty;

labels = cell(0, 1); % for x-axis of boxplots
pattern_label = {'+', 'x'}; % + is checkerboard, x is spokes
light_label = {'L', 'M', 'H'};

yEndObs = cell.empty;
for ct_y=1:length(y_start)
    for ct=1:length(data_all)
        labels{ct} = [light_label{data_all(ct).light(end)} '' pattern_label{data_all(ct).pattern(end)}];
        
        for ct1=1:1:length(data_all(ct).state_LDF)
            
            
            txyz = data_all(ct).state_LDF(ct1).filteredState(:,[1 2 3 4]); % the complete trajectory
            indx_start = find(txyz(:,3)<y_start(ct_y),1,'last');
            indx_end = find(txyz(:,3)<y_end,1,'last');
            
            if ~isempty(indx_start) && ~isempty(indx_end) && indx_end>indx_start && indx_end~=size(txyz,1)
                dt{ct_y}{ct}(ct1) = diff(interp1(txyz(indx_start:indx_end+1,3), txyz(indx_start:indx_end+1,1), [y_start(ct_y); y_end]));
                
                dummy = interp1(txyz(indx_start:indx_end+1,3), txyz(indx_start:indx_end+1,[2 3 4]), [y_start(ct_y); y_end]);
                xyz = [dummy(1,:); txyz(indx_start+1:indx_end-1, [2 3 4]); dummy(2,:)];
                ds{ct_y}{ct}(ct1) = sum((sum((diff(xyz)).^2,2)).^0.5);
                
%                 if ds{ct_y}{ct}(ct1)<1e-4
%                     keyboard;
%                 end
            else
                dt{ct_y}{ct}(ct1) = nan;
                ds{ct_y}{ct}(ct1) = nan;
            end
            
%             yEndObs{ct}(ct1) = txyz(end,3);
        end
    end
end

% yEndObs_figHandle = createBoxPlot(yEndObs, labels, 'y end observed');
for ct=2%1:length(y_start)
    figHandle = createBoxPlot(dt{ct}, labels, ['dt from y=' num2str(y_start(ct),'%0.2f') ' to y=' num2str(y_end,'%0.2f') 'm']);
    figHandle = createBoxPlot(ds{ct}, labels, ['ds from y=' num2str(y_start(ct),'%0.2f') ' to y=' num2str(y_end,'%0.2f') 'm']);
end

%% Functions used in this script
function figHandle = createBoxPlot(variable, labels, yxislabel)
    % for boxplots of a cell array
    col=@(x)reshape(x,numel(x),1);
    boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});

    
    figHandle = figure; hold on;
    figure(figHandle);
    boxplot2(variable, 'Labels', labels, 'OutlierSize', 0.00001);
%     hbox = gca;
    set(gca, 'FontSize', 18); %grid on;
%     xlabel('Treatments', 'FontSize', 18);
    ylabel(yxislabel, 'FontSize', 18);
    % ylim([0 1.5]);
    for ct = 1:length(variable)
            x=(ct+(rand(length(variable{ct}),1)-0.5)/4);

            f = scatter(x(:,1),variable{ct},10,'k','filled'); 
            f.MarkerFaceAlpha = 0.5;
    %         keyboard;
    end
end

