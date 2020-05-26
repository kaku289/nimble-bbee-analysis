classdef data4rrefEstimate < handle
    properties (Access = public)
        
        % state for rref estimation
        state4rrefEstimate = []; % NX10 vector - the whole shebang [ time x y z vx vy vz ax ay az]
        rref = nan; % obtained from linear fit of V vs y curve
        meanVbyy = nan; % obtained by taking mean of V/y
        vmean = nan; % Mean of Vgy in state4rrefEstimate
        ymean = nan; % Mean of y in state4rrefEstimate
        model; % linear fit model
        
        dof_analytical = nan; % Duration of flight calculated analytically 
        dof_actual = nan; % Actual duration of flight
        Rsquared = nan;
        
    end
     methods
         function obj = data4rrefEstimate()
             
         end
         
     end
end