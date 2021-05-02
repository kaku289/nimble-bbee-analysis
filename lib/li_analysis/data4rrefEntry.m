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
        day; % day of the experiment
        time; % time during the day (measured in blocks)
        patternnum;
        setID;
        beeID;
        flightID;
        
        slope_rvst ; % N by 1
        const_rvst ; % N by 1
        Rsquare_rvst; % N by 1
        
        isRise; % N y 1 (1 if the entry segment is rising, 0 if falling)
        ymean_for_rdot; % N by 1
        delta_r; % N by 1
        delta_y_analytical; % N by 1
        delta_y_actual; % N by 1
        
        yEntryStart;
        delta_Ventry;
        delta_tentry;
        amean_entry;
        
        %%%%% Parameters for statistical analysis %%%%%
        
        
        
    end
     methods
         function obj = data4rrefEntry()
             
         end
         
     end
end