%%
clc; clear; close all;

%%
pathTracking = '/media/reken001/Disk_08_backup/light_intensity_experiments/Tracking/';
trackingFolders = GetFolders(pathTracking,{'.mainbrain'});

pathCalib = '/media/reken001/Disk_08_backup/light_intensity_experiments/Calibration/_xml/';
calibFiles = GetFilesOnly(pathCalib, {'_2.xml'});
calibDates = cellfun(@(x) x(9:16), calibFiles, 'UniformOutput', false);%{calibFiles(:)}(9:14);

for ct_folder = 1:length(trackingFolders)
%     disp(trackingFolders{ct_folder});
    
    date = trackingFolders{ct_folder}(1:8);
    
    % move original file to original_calibration.xml
    if ~isfile(fullfile(pathTracking, trackingFolders{ct_folder}, 'original_calibration.xml'))
        movefile(fullfile(pathTracking, trackingFolders{ct_folder}, 'calibration.xml'), fullfile(pathTracking, trackingFolders{ct_folder}, 'original_calibration.xml'));
    end
    
    % find closest calibration file
    dateNumbers = datenum(changeDateFormat(calibDates));
    diffDateNumbers = abs(dateNumbers-datenum(changeDateFormat({date})));
    indx = find( diffDateNumbers == min(diffDateNumbers), 1, 'last');
    
    disp(['Copying calibration from ' calibFiles{indx(end)} ' to ' trackingFolders{ct_folder}]);
    % copy the closest calibration file as calibration.xml
    copyfile(fullfile(pathCalib, calibFiles{indx(end)}), fullfile(pathTracking, trackingFolders{ct_folder}, 'calibration.xml'));
end
%% Function used in this script file
function outputDates = changeDateFormat(dates)
    % converts dates (in the format 'YYYYMMDD') to outputDates (in the
    % format 'YYYY-MM-DD')
    % dates - cell structure input ({'YYYYMMDD'} {'YYYYMMDD'})
    outputDates = cell(size(dates));
    for ct=1:length(dates)
        outputDates{ct} = [dates{ct}(1:4) '-' dates{ct}(5:6) '-' dates{ct}(7:8)];
    end
end

function [files] = GetFilesOnly(direc,filePattern)
    % Get a list of all files in a directory containing filePatterns in their names.

    % This function will not provide files starting with ~ (e.g. temporary
    % excel files)
    % direc - directory in which list of files is required
    % filePattern - cell structure containing patterns that uniquely identify all files to be listed
    % files - cell structure output containing list of file names that
    % contain ALL of the patterns listed in filePattern

    allFiles    = dir(direc);
    names    = {allFiles.name};

    % Get a logical vector that tells which is a file.
    filesFlags = ~[allFiles.isdir] & ~strcmp(names, '.') & ...
              ~strcmp(names, '..') & ~strcmp(cellfun(@(x){x(1)},names),'~');%~cellfun('isempty',strfind(A,B));
    for ct=1:length(filePattern)
       filesFlags = filesFlags & contains(names,filePattern{ct});
    end

    % Extract only those contain filePattern
    files = names(filesFlags);
end

function [folders] = GetFolders(direc,pattern)
    % Get a list of all folders in a directory containing patterns in their names.
    
    % direc - path where to look
    % pattern - cell array containing string patterns that must be in the
    % names of the folders

    all    = dir(direc);
    folders = all([all(:).isdir]);
    names    = {folders.name};

    % Get a logical vector that discards folders that aren't required
    flags = ~strcmp(names, '.') & ...
              ~strcmp(names, '..') & ~strcmp(cellfun(@(x){x(1)},names),'~');%~cellfun('isempty',strfind(A,B));
    for ct=1:length(pattern)
       flags = flags & contains(names,pattern{ct});
    end

    % Extract only those that contain pattern
    folders = names(flags);
end