classdef outputForParameterEstimation < handle
    properties (Access = public)
        signalForEstimation = 'ay'; % 'r' or 'ay'
        x0 = []; % Initial values of the parameters
        xOpt = []; % Optimized values of the parameters for that x0
        cost = [];
        eflag = [];
        output = [];
        
    end
     methods
         function obj = outputForParameterEstimation()
             
         end
     end
end