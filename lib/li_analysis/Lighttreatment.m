classdef Lighttreatment < handle
    properties (Access = public)
        datenum = nan;
        pattern = '';
        light = '';
        startTime = nan;
        endTime = nan;
        
        usefulVideoFiles 
        trackingFiles % a file struct
        tempRH 
        calib = FlydraCalibration.empty;
        
        landingDiscs = LandingDisc.empty; 
        
        landingTracks = BlindLandingtrack.empty; % to store blind landing tracks
        
        videosParsed = false; % If the videos has been scored and tracks have been "projected onto them"
        
        videosInfo = recordedVideosInformation.empty; % to store information about the recorded videos
    end
     methods
         function obj = Lighttreatment(datenum, pattern, light, startTime, endTime)
             % datenum - YYYYMMDD
             % pattern - checkerboard or spokes or grayscale
             % light - low/medium/high/lower
             % startTime - HHMMSS
             % endTime - HHMMSS
             % centers - [X1 Y1 Z1; X2 Y2 Z2] - where 1 and 2 refers to two
             % different sides (hive and feeder)
             % sides - {'Hive', 'Feeder'}
             
             obj.datenum = datenum;
             obj.pattern = pattern;
             obj.light = light;
             obj.startTime = startTime;
             obj.endTime = endTime;
             
             obj.usefulVideoFiles = VideoMetadata.empty;

         end
         
         function addLandingDisc(obj, cam_ids_and_points2D, side)
             if ~isempty(obj.landingDiscs) && any(strcmpi({obj.landingDiscs.side}, side))
                 disp('Landing disc already exists corresponding to this side. Updating the data2d...');
                 obj.landingDiscs(strcmpi({obj.landingDiscs.side}, side)) = LandingDisc(obj.datenum, obj.pattern, cam_ids_and_points2D, side);
             else
                obj.landingDiscs(end+1) = LandingDisc(obj.datenum, obj.pattern, cam_ids_and_points2D, side);
             end
         end
     end
end