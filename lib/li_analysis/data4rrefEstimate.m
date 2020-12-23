classdef data4rrefEstimate < handle
    properties (Access = public)
        params ; % sigmas used for thresholding to identify constant r segments
        factor ; % another factor used for thresholding to identify constant r segments
        intervalArray = Interval.empty ; % Interval array vector containing all intervals that satisfy the criterion
       
        intervals_ti ; % N by 2, where N is the number of non-overlapping time-window independent intervals
        rref_ti ; % N nby 1 - linear fit of V vs y with zero intercept 
        vmean_ti ; % N by 1
        ymean_ti ; % N by 1
        rmean_ti; % N by 1 - mean of r
        slope_rvsy_ti ; % N by 1
        const_rvsy_ti ; % N by 1
        
        xTravelled_ti ;
        yTravelled_ti ;
        zTravelled_ti ;
        xmean_ti ;
        zmean_ti ;
        fd_analytical_ti ; % analytically computed flight duration 
        fd_actual_ti ; % actual flight duration
        
        speed3d_mean_ti;
        rmean_speed3d_ti;
        
        intervals_td ; % M by 2, where M is the number of intervals for all time-windows. For each time-window, these intervals are non-overlapping, but across time-windows they can be overlapping
        rref_td ; % N nby 1
        vmean_td ; % N by 1
        ymean_td ; % N by 1
%         rmean_td ; % N by 1
%         slope_rvst_td ; % N by 1
%         const_rvst_td ; % N by 1
        
        %%%%% Parameters for statistical analysis %%%%%
        day; % day of the experiment
        time; % time during the day (measured in blocks)
        pattern; % 1 (checkerboard) or 2 (spokes)
        light; % 1 (low), 2 (medium), 3 (high)
        side; % 1 (hive), 2 (food source)
        wind; % 1 to 6

        
        % state for rref estimation - for legacy code
        state4rrefEstimate = []; % NX10 vector - the whole shebang [time x y z vx vy vz ax ay az]
        rref ; % obtained from linear fit of V vs y curve
        meanVbyy ; % obtained by taking mean of V/y
        vmean ; % Mean of Vgy in state4rrefEstimate
        ymean ; % Mean of y in state4rrefEstimate
        model; % linear fit model
        
        dof_analytical ; % Duration of flight calculated analytically 
        dof_actual ; % Actual duration of flight
        Rsquared ;
        
    end
     methods
         function obj = data4rrefEstimate()
             
         end
         
     end
end