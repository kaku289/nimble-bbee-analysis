%% Miscellaenous code used for extracting information required during writing the manuscript

%% Find how many tracks from spoke, high condition were manually analyzed
% Last track analysed manually
% 
% BlindLandingtrack with properties:
% 
%               foldername: '20190704_110002.mainbrain.offline'
%                  datenum: 20190704
%                  pattern: 'spokes'
%                    light: 'high'
%                   obj_id: 4196
%                 rawTrack: [1×1 rawState_BlindLandingtrack]
%                    state: [1×1 filteredState_BlindLandingtrack]
%           dayTimeInstant: ''
%     treatmentTimeInstant: ''
%                     time: []
%             displacement: []
%                     vbyz: NaN
%                       dt: 0.0057
%                state_LDF: [1×1 filteredState_BlindLandingtrack]
%                  DataGUI: [1×1 DataGUI_BlindTracks]
%        equinox_state_LDF: []
%         stable_state_LDF: []
%       track_subset_sysID: [1×1 trackForLandingModel]

last_track.datenum = 20190704;
last_track.pattern = 'spokes';
last_track.light = 'high';
last_track.obj_id = 4196;

% Collect all landing tracks in spoke,high condition
close all;
clc;
% inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A1_rref.mat';
% load(inputFile);
treatments = treatments(1:14*8);

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
behaviour = {'rising','constant','sleeping'};
for ct_pattern = 2%1:length(pattern)
    for ct_light = 3%1:length(light)
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
            
            indx = [landingTracks.datenum] == last_track.datenum & [landingTracks.obj_id] == last_track.obj_id;
            disp(['# of landing tracks manually analyzed: ' num2str(find(indx))]);
            track = landingTracks(indx);
            
  
        end
    end
end

%% Create data files to be uploaded with manuscript
clc; close all;
clear;
% 
inputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/BlindLandingtracks_A1_rref.mat';
load(inputFile);
treatments = treatments(1:14*8); % Taking experiments for 2 patterns * 3 lights

outputFile = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing/landingTracks_forUpload.mat';

pattern = {'checkerboard', 'spokes'};
light = {'low', 'medium', 'high'};
behaviour = {'rising','constant','sleeping'};

factors = [0.25:0.25:2.5];
chosen_fac = 1;

clear dummy data;
data = struct.empty;

for ct_pattern = 1:length(pattern)
    for ct_light = 1:length(light)
        dummy = struct.empty;
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
            for ct_treatment = 1:length(relevantTreatments)
                landingTracks = [relevantTreatments(ct_treatment).landingTracks];
                state_LDF = [landingTracks.state_LDF];
                
                for ct_track=1:length(state_LDF)
                    dummy(end+1).state = state_LDF(ct_track).filteredState(:,1:4);
                    dummy(end).state(:,3) = -dummy(end).state(:,3);
                    dummy(end).date = relevantTreatments(ct_treatment).datenum;
                    dummy(end).pattern = relevantTreatments(ct_treatment).pattern;
                    dummy(end).light = relevantTreatments(ct_treatment).light;
                    dummy(end).treatmentStartTime = relevantTreatments(ct_treatment).startTime;
                    dummy(end).treatmentEndTime = relevantTreatments(ct_treatment).endTime;
                end
            end
            
            disp(size(dummy));
            
            data = [data dummy];
            
        end
    end
end
save(outputFile, 'data');
