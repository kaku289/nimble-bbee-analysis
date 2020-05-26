classdef rawState_BlindLandingtrack < handle
    properties (Access = public)
        rawState = []; % 18 columns (obj_id	frame	timestamp	x	y	z	xvel	yvel	zvel P00	P01	P02	P11	P12	P22	P33	P44	P55)
        landingSide = ''; % Hive or Feeder
    end


    methods
        function obj = rawState_BlindLandingtrack(rawState, landingSide)
            assert(size(rawState,2) == 18);
            obj.rawState = rawState;
            obj.landingSide = landingSide;
        end
        
        function instance = createNewInstanceFromSubset(obj,indx1,indx2)
            % instance - instance of rawState_BlindLandingtrack created
            % from the obj.rawState(indx1:indx2,:)
            instance = rawState_BlindLandingtrack(obj.rawState(indx1:indx2,:), obj.landingSide);
        end
    end
end
        