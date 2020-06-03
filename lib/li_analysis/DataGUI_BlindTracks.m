classdef DataGUI_BlindTracks < handle
    properties (Access = public)
        
        % % Start and end of the U shape in the plots of Vgy vs y is
        % selected. Max. 5 U shaped per track are allowed
        isExcerptUseful = false;
        y = zeros(5,2);
        Vgy = zeros(5,2);
        
        % Mode 1 - to estimate r*
        y_mode1 = zeros(5,2);
        Vgy_mode1 = zeros(5,2);
        
        % Mode 2 - to select track segments with rise transients
        t_mode2 = zeros(5,2);
        y_mode2 = zeros(5,2);
        Vgy_mode2 = zeros(5,2);
        
%         state = []; % storing state to be used for optimization (for positive times)
%         negTime_r = []; % To store state for negative time
%         
%         % values set after optimization is complete
%         param_estimation = outputForParameterEstimation.empty;
    end
     methods
         function obj = DataGUI_BlindTracks()
             
         end
     end
end