%% Miscellaenous code used for extracting information required during writing the manuscript A1

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