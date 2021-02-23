classdef appData_BlindtracksGUI < handle
    properties (Access = public)
        currentTrack = 1;
        dataDir = '/media/reken001/Disk_08_backup/light_intensity_experiments/postprocessing';
        filename = '';
        mode = 1; % 1(to select r*) or 2 (to select "entering" into r* excerpt
        part = 1; % to select multiple (mode 1 or 2) segments within the same track (it can be max upto 5)
    end
     methods
         function obj = appData_BlindtracksGUI()
             
         end
     end
end