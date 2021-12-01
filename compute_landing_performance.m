%% Compute landing performance
clc; close all;
% clear;
% 
inputFile = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/BlindLandingtracks_A3_LDF_rref.mat';
% load(inputFile);

winds = unique([treatments.wind]);
behaviour = {'rising','constant','sleeping'};

y = [0.25 0.05]; % y's between which performance parameters are calculated
clear dummy;

data = struct.empty;
for ct_wind = 1:length(winds)
    
    for ct_behaviour = 2%1:length(behaviour)
        
        disp([' wind: ' num2str(winds(ct_wind)) ...
            ', behaviour: ' behaviour{ct_behaviour}]);
        
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
        
        % Only analyse treatments that have uniform wind measurements
        % avaliable throughout their course of running
        hasUniformHwData = arrayfun(@(x) x.hwData.hasUniformHwData,relevantTreatments);
        relevantTreatments = relevantTreatments(hasUniformHwData);
        
        %%%%%%%%%%%%%% Each wind condition %%%%%%%%%%%
        landingTracks = [relevantTreatments.landingTracks];
        state_LDF = [landingTracks.state_LDF];
        
        treatment_indx4landingTracks = arrayfun(@(x) x*ones(length(relevantTreatments(x).landingTracks),1),1:length(relevantTreatments), 'UniformOutput', false);
        treatment_indx4landingTracks = vertcat(treatment_indx4landingTracks{:});
        treatment_indx4stateLDF = arrayfun(@(x) treatment_indx4landingTracks(x)*ones(length(landingTracks(x).state_LDF),1),1:length(landingTracks), 'UniformOutput', false);
        treatment_indx4stateLDF = vertcat(treatment_indx4stateLDF{:});
        
        startTimes = [relevantTreatments(treatment_indx4stateLDF).startTime];
        
        landingTracks_indx4stateLDF = arrayfun(@(x) x*ones(length(landingTracks(x).state_LDF),1),1:length(landingTracks), 'UniformOutput', false);
        landingTracks_indx4stateLDF = vertcat(landingTracks_indx4stateLDF{:});
        
        landingDiscs = {relevantTreatments(treatment_indx4stateLDF).landingDiscs};
        hastakeoff = arrayfun(@(x) state_LDF(x).hasTakeoff(landingDiscs{x}),1:length(state_LDF));
        
        performance_params = arrayfun(@(x) x.compute_landing_performance(y(1),y(2)), state_LDF, 'UniformOutput', false);
        performance_params = vertcat(performance_params{:});
        
        delta_t = performance_params(:,1);
        
        indices = ~isnan(delta_t);
        
        % store data
        dummy.tracks_fac = state_LDF(indices);
        dummy.wind = winds(ct_wind);
        dummy.landingTrack = landingTracks(landingTracks_indx4stateLDF(indices));
        dummy.hastakeoff = hastakeoff(indices);
        dummy.delta_t = delta_t(indices);
        dummy.startTimes = startTimes(indices);
        data = [data; dummy];
    end
end

cmap = jet(6);
cmap([2 5],:) = [];
colors = [cmap; 1 0 1; 1 0 0];
figure;
histogram(data(1).delta_t,'facecolor',colors(1,:),'facealpha',.5,'edgecolor','none','Normalization','pdf');
hold on;
histogram(data(6).delta_t,'facecolor',colors(6,:),'facealpha',.5,'edgecolor','none','Normalization','pdf');
legend('1','6','fontsize',16)

%%
data_write = [];
for ct=1:length(winds)
    N = length(data(ct).delta_t);
    landingSide = {data(ct).tracks_fac.landingSide};
    isHive = cellfun(@(x) strcmpi(x,'hive'),landingSide);
    side = ones(length(isHive),1);
    side(~isHive) = deal(2);
    
    data_write = [data_write;
                  data(ct).delta_t data(ct).wind*ones(N,1) data(ct).hastakeoff' data(ct).startTimes' side [data(ct).landingTrack.datenum]'];
end

writeFile = true;
r_file = '/media/reken001/Disk_07/steady_wind_experiments/postprocessing/data_landing_performance_Rstudio.txt';

if writeFile
    T = array2table(data_write, ...
        'VariableNames',{'delta_t','wind','hasTakeoff',...
                         'time','landingSide','day'});
    writetable(T,r_file);
end
