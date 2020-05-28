classdef rrefEstimateInterval < Interval
    properties (Access = public)
%         params ; % sigmas used for thresholding to identify constant r segments
        factor ; % another factor used for thresholding to identify constant r segments
        
        xTravelled ;
        yTravelled ;
        zTravelled ;
        xmean ;
        zmean ;
        fd_analytical ; % analytically computed flight duration 
        fd_actual ; % actual flight duration
        slope_rvsy ; % N by 1
        const_rvsy ; % N by 1

        rref ; % obtained from linear fit with zero intercept of V vs y curve
        rmean ; % obtained by taking mean of V/y
        vmean ; % Mean of Vgy in state4rrefEstimate
        ymean ; % Mean of y in state4rrefEstimate
        
    end
     methods
         function obj = rrefEstimateInterval()
             % Calls superclass (Interval) constructor
             
         end
         
     end
end