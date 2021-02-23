classdef videoParameters < matlab.mixin.Copyable
  properties
      isVideoUseful = false;
      % True - if video contains a landing (not a false trigger)
      % False - if video is a false trigger
      
      year = NaN;
      month = NaN;
      day = NaN;
      hour = NaN;
      minutes = NaN;
      seconds = NaN;
      filenum = NaN; % file identifier (HHMMSS in name)
      
      path;
      name;
      landingParams = landingParameters();
      
      hasMultipleLandings = 0;
      
      analyzeVideo = false; %if the video file needs to be analyzed or not
      isScored = false; % if the video file has been scored or not
      
      timestamps = [];
      
  end
    methods
        function obj = videoParameters()
            obj.analyzeVideo = false;
            obj.hasMultipleLandings = 0;
            obj.isScored = 0;
            
            obj.landingParams = landingParameters();
        end
    end
    methods(Access = protected)
        % Override copyElement method:
        function cpObj = copyElement(obj)
            % Make a shallow copy of all properties
            cpObj = copyElement@matlab.mixin.Copyable(obj);
            
            % Make a deep copy of the DeepCp object
            cpObj.landingParams = copy(obj.landingParams);
        end
    end
end