classdef landingParameters < matlab.mixin.Copyable
  properties
      legextension_x = NaN
      legextension_y = NaN
      legextension_framenumber = NaN
      
      % When first body appendage touches the platform
      touchdown_x = NaN
      touchdown_y = NaN
      touchdown_framenumber = NaN
      
      % Which body appendage touched the platform
      touchdown_bodyappendage = NaN
      
      % When all legs touch the platform
      touchdown2_x = NaN
      touchdown2_y = NaN
      touchdown2_framenumber = NaN
      
      landing_success = 2
      % 2 - Direct landing
      % 1 - Not direct, but eventually landed after some tries
      % 0.5 - Platform is touched and bbee flew away 
      % 0 - Got bounced off or collided
      % -1 - Platform is not touched and bbee flew away
      
      landing_side = 'Hive'
      % 0 - Hive
      % 1 - Feeder
      
      landing_location = 'Tube'
      % 0 - Tube
      % 1 - Disc
      
      isGoldenFlight = 0;      
  end
    methods
        function obj = landingParameters()
            
        end
    end
    methods(Access = protected)
        
    end
    
end