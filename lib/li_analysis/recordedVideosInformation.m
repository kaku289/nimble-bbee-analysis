classdef recordedVideosInformation < handle
    properties (Access = public)
        startTime = nan;
        endTime = nan;
        timestamps = [];
        nFrames = nan;
        name = '';
        folder = '';
        cam_id = nan; % 84, 85, 87 or 25
    end
     methods
         function obj = recordedVideosInformation()
             obj.startTime = nan;
             obj.endTime = nan;
             obj.timestamps = [];
             obj.nFrames = nan;
             obj.name = '';
             obj.folder = '';
             obj.cam_id = nan; % 84, 85, 87 or 25
         end
         function obj = extractInformation(obj, dirStruct)
             obj.name  = dirStruct.name;
             obj.folder = dirStruct.folder;
             obj.cam_id = str2double(dirStruct.name(end-5:end-4));
             
             [~, obj.timestamps] = fmf_read(fullfile(obj.folder, obj.name));
             obj.nFrames = length(obj.timestamps);
             obj.startTime = obj.timestamps(1);
             obj.endTime = obj.timestamps(end);
             
         end
         
     end
end