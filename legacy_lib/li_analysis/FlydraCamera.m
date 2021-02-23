classdef FlydraCamera < handle
    properties (Access = public)
        fc1 ;
        fc2 ;
        cc1 ;
        cc2 ;
        k1 ;
        k2 ;
        k3 = 0;
        p1 ;
        p2 ;
        alpha_c ;
        
        pmat ;
        pmat_inv ;
        
        resolution ;
        
        cam_id = '';
        
        cam_center ;
        
    end
     methods
         function obj = FlydraCamera(cameraStruct)
             
             for i=1:length(cameraStruct) % for each children in single_camera_calibration
                 if strcmpi(cameraStruct(i).Name,'cam_id')
                     obj.cam_id = cameraStruct(i).Children.Data;
                 elseif strcmpi(cameraStruct(i).Name,'calibration_matrix')
                     obj.pmat = str2num(cameraStruct(i).Children.Data);
                     obj.pmat_inv = pinv(obj.pmat);
                     obj.pmat2cam_center();
                 elseif strcmpi(cameraStruct(i).Name,'resolution')
                     obj.resolution = str2num(cameraStruct(i).Children.Data);
                 elseif strcmpi(cameraStruct(i).Name,'non_linear_parameters')
                     for j=1:length(cameraStruct(i).Children)
                         obj.([cameraStruct(i).Children(j).Name]) = str2double(cameraStruct(i).Children(j).Children.Data);
                     end
                 end
                     
                     
                 
             end
             
             
         end
         
         function pmat2cam_center(obj)
             % pmat - 3 X 4 matrix
             % Hartley & Zisserman (2003) p. 163
             
             X = det(obj.pmat(:, [2 3 4]));
             Y = -det(obj.pmat(:, [1 3 4]));
             Z = det(obj.pmat(:, [1 2 4]));
             T = -det(obj.pmat(:, [1 2 3]));
             
             obj.cam_center = [X/T Y/T Z/T]';
         end
     end
     
     methods
     end
end