classdef data4rrefEntry < handle
    properties (Access = public)
        
        factor ; % another factor used for thresholding to identify constant r segments
        
        intervals; % N by 2, where N is the number of non-overlapping time-window independent intervals
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
        day; % day of the experiment
        time; % time during the day (measured in blocks)
        pattern; % 1 (checkerboard) or 2 (spokes)
        light; % 1 (low), 2 (medium), 3 (high)
        side; % 1 (hive), 2 (food source)
        
    end
     methods
         function obj = data4rrefEntry()
             
         end
         
     end
end