classdef hotwireData < handle
    properties (Access = public)
        meanVoltage;
        stdVoltage;
        meanVoltage_persec = {};
        stdVoltage_persec = {};
        
        meanVel;
        stdVel;
        meanVel_persec = {};
        stdVel_persec = {};
        
        hwFiles = struct.empty;
        recordingDate = {};
        recordingTime = {};
        fs;
        
        hasUniformHwData;
    end
    methods
        function obj = hotwireData(hotwireFiles)
            % hotwireFiles - a dir structure
            
            obj.hwFiles = hotwireFiles;
            
            if ~isempty(hotwireFiles)
                obj.findParameters();
            end
            
            
        end
        
        function findParameters(obj)
            
            for ct=1:length(obj.hwFiles)
                obj.recordingDate{ct} = str2double(obj.hwFiles(ct).name(1:8));
                obj.recordingTime{ct} = str2double(obj.hwFiles(ct).name(10:15));
                obj.fs = str2double(obj.hwFiles(ct).name(17:20));
                data{ct} = readmatrix(fullfile(obj.hwFiles(ct).folder, obj.hwFiles(ct).name));
                obj.meanVoltage_persec{ct} = arrayfun(@(x)  mean(data{ct}((x-1)*obj.fs+1:x*obj.fs,1)), 1:size(data{ct},1)/obj.fs);
                obj.stdVoltage_persec{ct} = arrayfun(@(x)  std(data{ct}((x-1)*obj.fs+1:x*obj.fs,1)), 1:size(data{ct},1)/obj.fs);
                
                
            end
            data = vertcat(data{:});
            obj.meanVoltage = mean(data(:,1));
            obj.stdVoltage = std(data(:,1));
        end
        
        
    end
end