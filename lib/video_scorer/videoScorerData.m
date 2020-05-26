classdef videoScorerData < handle
    properties
%         gui;
        moviesdir = ''; % Directory containing the video files
        
        videoFiles = videoParameters.empty; % Array containing names of the video files
        year;
        month;
        day;
        
        doesDataFileExist = false;
        datadir = '';  % directory where data is saved or needs to be saved
        
        currentFileIndex = 0;
       
    end
    methods 
%         function r = videoScorerGUI(obj, GUI)
%             obj.gui = GUI;
%         end

    end
end