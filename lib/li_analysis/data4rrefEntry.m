classdef data4rrefEntry < handle
    properties (Access = public)
        
        factor ; % another factor used for thresholding to identify constant r segments
        
        intervals; % N by 3, where N is the number of non-overlapping time-window independent intervals
        rref; % N nby 1 - linear fit of V vs y with zero intercept 
        vmean; % N by 1
        ymean; % N by 1
        rmean; % N by 1 - mean of r
%         models; % N by 1 iddtf models
        
        xTravelled;
        yTravelled;
        zTravelled;
        xmean;
        zmean;
        
        fd_analytical_ti; % analytically computed flight duration
        fd_actual_ti; % actual flight duration
        
        %%%%% Parameters for statistical analysis %%%%%
        day = ''; % date on which the track was obtained
        pattern = '';
        patternnum;
        setID;
        beeID;
        flightID;
        
    end
     methods
         function obj = data4rrefEntry()
             
         end
         
     end
end