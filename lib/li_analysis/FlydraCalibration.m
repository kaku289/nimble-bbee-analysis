classdef FlydraCalibration < handle
    properties (Access = public)
        
        name = '';
        cam = FlydraCamera.empty;
        
    end
    methods
         function obj = FlydraCalibration(xmlfile)
             
             a = parseXML(xmlfile);
             obj.name = a.Name;
             
             for ct_cam = 1:length(a.Children) % for each camera
                 obj.cam(ct_cam) = FlydraCamera(a.Children(ct_cam).Children);
                 
             end
             
         end
         
         function X = find3d(obj, cam_ids_and_points2D)
             % This function computes the 3D point for given 2D points
             
             % cam_ids_and_point2D - cell array {NX2 size}
             % 1st column contains name of the camera
             % 2nd column is 1X2 array containing image points in pixels
             A = [];
             for ct=1:size(cam_ids_and_points2D,1) % for each image point
                 camera = obj.cam(strcmpi({obj.cam.cam_id}, cam_ids_and_points2D{ct,1}));
                 pmat = camera.pmat;
                 xd = cam_ids_and_points2D{ct,2}(1);
                 yd = cam_ids_and_points2D{ct,2}(2);
                 if ~isempty(xd) && ~isempty(yd) && ~ismissing(xd) && ~ismissing(yd)
                    [xl, yl] = FlydraCalibration.undistort(camera, xd, yd);
                    A(end+1,:) = [xl*pmat(3,:) - pmat(1,:)];
                    A(end+1,:) = [yl*pmat(3,:) - pmat(2,:)];                 
                 end
             end
             [~,~,V] = svd(A);
             X = V(1:end-1,end)./V(end,end);             
         end
     end
     
     methods (Static=true, Access=public)
         function xy = findUndistorted2Dfeature(cam, XYZ)
             % This function computes 2D undistorted points for a 3D point
             % XYZ in a camera cam
             
             % cam - instance of FlydraCamera class
             % XYZ - 3 by N matrix
             % xy - 2 by N matrix containing undistorted coordinates for
             % each 3D point
             
             assert(size(XYZ,1)==3, '3D point matrix must have 3 rows');
             xy = cam.pmat*[XYZ; ones(1, size(XYZ, 2))];
             xy = xy(1:2,:)./xy(3,:);
         end
         
         function [xd, yd] = distort(cam, xl, yl)
             % cam - instance of FlydraCamera class
             % xl, yl - undistorted points in pixels (each 1 by N size)
             % xd, yd - 1 by N size each
             
             x = ( xl - cam.cc1 ) / cam.fc1;
             y = ( yl - cam.cc2 ) / cam.fc2;

             r_2 = x.*x + y.*y;
             r_4 = r_2.^2;
             r_6 = r_2.*r_4;

             term1 = cam.k1*r_2 + cam.k2*r_4 + cam.k3*r_6;

             xd = x + x.*term1 + (2*cam.p1*x.*y + cam.p2*(r_2+2*x.^2));
             yd = y + y.*term1 + (cam.p1*(r_2+2*y.^2) + 2*cam.p2*x.*y);

             xd = (cam.fc1)*xd + (cam.cc1);
             yd = (cam.fc2)*yd + (cam.cc2);
             
         end
         
         function [xl, yl] = undistort(cam, xd, yd)
             % cam - instance of FlydraCamera class
             % xd, yd - distorted points in pixels (each 1 by N size)
             % xd, yd - 1 by N size each
             
             n_iterations = 5;
             
             xd = ( xd - cam.cc1 )/ cam.fc1;
             yd = ( yd - cam.cc2 )/ cam.fc2;

             xd = xd - cam.alpha_c*yd;

             x = xd;
             y = yd;

             for i=1:n_iterations
                 r_2 = x.*x + y.*y;
                 k_radial = 1.0 + r_2.*(cam.k1 + r_2.*(cam.k2 + r_2*cam.k3));
                 delta_x = 2.0*(cam.p1).*x.*y + (cam.p2).*(r_2 + 2.0.*x.*x);
                 delta_y = (cam.p1).*(r_2 + 2.0.*y.*y)+2.0*(cam.p2).*x.*y;
                 x = (xd-delta_x)./k_radial;
                 y = (yd-delta_y)./k_radial;
             end

            xl = (cam.fc1).*x + (cam.fc1*cam.alpha_c).*y + (cam.cc1);
            yl = (cam.fc2).*y + (cam.cc2);
             
         end
         
         
     end
     
end