classdef LandingDisc < handle
    properties (Access = public)
        datenum = nan; % YYYYMMDD
        pattern = ''; % checkerboard or spokes
        center = []; % [X Y Z] in m
        radius = 0.181/2; % in m
        side = ''; % Hive or Feeder
        
        cam_ids_and_points2D = {};% cell array {NX2 size}
                                 % 1st column contains name of the camera
                                 % 2nd column is 1X2 array containing image points in pixels
        
    end
     methods
         function obj = LandingDisc(datenum, pattern, cam_ids_and_points2D, side)
             obj.datenum = datenum;
             obj.pattern = pattern;
             
             obj.cam_ids_and_points2D = cam_ids_and_points2D;
             
             obj.side = side;
             
         end
         
         function computeCenter(obj, calib)
             % calib - FlydraCalibration instance
             obj.center = calib.find3d(obj.cam_ids_and_points2D);
         end
     end
end