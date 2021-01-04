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
%% Find correlation between accleration and r*-r(t) signal
clc; close all;
% clear;
% 
% inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A2_rrefEntry.mat';
inputFile = 'D:/light_intensity_experiments/postprocessing/BlindLandingtracks_A2_rrefEntry.mat';
% load(inputFile);


treatments = treatments(1:14*8); % Taking experiments for 2 patterns * 3 lights


pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
behaviour = {'rising','constant','sleeping'};

labels = cell(0, 1); % for x-axis of boxplots
pattern_label = {'+', 'x'}; % + is checkerboard, x is spokes
light_label = {'L', 'M', 'H'};
delay = cell(length(pattern),length(light));
labels = cell(length(pattern),length(light));

factors = [0.25:0.25:2.5];

entryData = cell.empty;

for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        for ct_behaviour = 2%1:length(behaviour)
            disp(' ');
            disp(['Pattern: ' pattern{ct_pattern} ...
                  ', light: ' light{ct_light} ...
                  ', behaviour: ' behaviour{ct_behaviour}]);
              
            labels{ct_pattern, ct_light} = [light_label{ct_light} '' pattern_label{ct_pattern}];
            delay{ct_pattern, ct_light} = [];
            
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
            rrefEntrySegs = [state_LDF.rrefEntrySegments];
            
            landingTracks_indx4stateLDF = arrayfun(@(x) x*ones(length(landingTracks(x).state_LDF),1),1:length(landingTracks), 'UniformOutput', false);
            landingTracks_indx4stateLDF = vertcat(landingTracks_indx4stateLDF{:});
            
            % Finding # of "tracks" (state LDFs) that contain rref entry segments
            % for each factor
            for ct_fac = 4%1:length(factors)
                factor = factors(ct_fac);
                indices = arrayfun(@(x) ~isempty(x.rrefEntrySegments(abs([x.rrefEntrySegments.factor]-factor)<1e-6).intervals),state_LDF);
                
                tracks_fac = state_LDF(indices);
                
                entryData{ct_fac} = {};
                
                for ct1=1:sum(indices) % for each state_LDF
                    state = tracks_fac(ct1).filteredState;
                    y = -state(:,3); V = state(:,6); ay = state(:,9); r = -state(:,6)./state(:,3); t = state(:,1)-state(1,1);
                    rdot = diffxy(t, r);
                    
                    rrefEntrySegments_fac = tracks_fac(ct1).rrefEntrySegments(abs([tracks_fac(ct1).rrefEntrySegments.factor]-factor)<1e-6);
                    
                    for ct2=1:size(rrefEntrySegments_fac.intervals,1)
                        entryStart_indx = rrefEntrySegments_fac.intervals(ct2,1);
                        entryEnd_indx = rrefEntrySegments_fac.intervals(ct2,2);
                        
                        rref = -rrefEntrySegments_fac.rmean(ct2);
                        y_part = y(entryStart_indx:entryEnd_indx); V_part = V(entryStart_indx:entryEnd_indx); r_part = r(entryStart_indx:entryEnd_indx);
                        ay_part = ay(entryStart_indx:entryEnd_indx); rdot_part = rdot(entryStart_indx:entryEnd_indx);
                        
                        entryData{ct_fac}{end+1} = [...
                                              t(entryStart_indx:entryEnd_indx)  y(entryStart_indx:entryEnd_indx) ...
                                              V(entryStart_indx:entryEnd_indx)  ay(entryStart_indx:entryEnd_indx) ...
                                              r(entryStart_indx:entryEnd_indx)  rdot(entryStart_indx:entryEnd_indx) ...
                                              rref*ones(entryEnd_indx-entryStart_indx+1,1)];
                                          
                        [rdiff, lags] = xcorr(ay_part, r_part-rref);
                        [~,idx] = max((rdiff));
                        delay{ct_pattern, ct_light}(1, end+1) = lags(idx);
                        
%                         figureHandle = figure;
%                         stem(lags,rdiff)
%                         
%                         pause(2)
%                         close(figureHandle);
                        
                    end
                    
                end

            end
        end
    end
end
%%
delay_pos = cellfun(@(x) x(x>0), delay, 'UniformOutput', false);
delay_neg = cellfun(@(x) x(x<0), delay, 'UniformOutput', false);

delay_figHandle = createBoxPlot({delay{:}}, {labels{:}}, 'Fit percent');

delay_figHandle1 = createBoxPlot({delay_pos{:}}, {labels{:}}, 'Fit percent');
delay_figHandle2 = createBoxPlot({delay_neg{:}}, {labels{:}}, 'Fit percent');
%%
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
