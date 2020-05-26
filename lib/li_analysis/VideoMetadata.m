classdef VideoMetadata < handle
    properties (Access = public)
        scoredData = videoParameters.empty;
        data2d = []; % 2d data corresponding to the video
        data3d = []; % 3d data corresponding to the video and not the whole 3D track
        data_association = []; % to find association of 2d and 3d points
        landingTrack;
%         hasLandingTrack = false;
        
        pattern = '';
        light = '';
        trackingFileName = '';
        
        datenum = ''; % YYYYMMDD
        
    end
    methods
        function obj = VideoMetadata(data2d, data3d, data_association, pattern, light, confirmedChecks, track, videoParams, trackingFileName)
            obj.data2d = data2d;
            obj.data3d = data3d;
            obj.data_association = data_association;
            obj.scoredData = videoParameters.empty;
            
            obj.pattern = pattern;
            obj.light = light;
            obj.trackingFileName = trackingFileName;
            
            obj.datenum = [videoParams.year videoParams.month videoParams.day];
            
%             if ~isempty(track) 
%                 obj.hasLandingTrack = true;
            obj.landingTrack = Landingtrack(pattern, light, confirmedChecks, track, videoParams);
%             else
%                 obj.landingTrack = Landingtrack.empty;
%             end

            
        end
    end
end