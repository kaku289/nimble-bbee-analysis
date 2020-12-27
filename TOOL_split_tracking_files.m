clc; clear; close all;

%%
path = '/media/reken001/Disk_09/unsteady_wind_experiments/Tracking';
outputPath = '/media/reken001/Disk_09/unsteady_wind_experiments/postprocessing';

% deleteOldfiles = true;
% if deleteOldfiles 
%     delete(fullfile(outputPath, '*mainbrain'))
% end

file3d = 'kalman_estimates.csv';
file2d = 'data2d_distorted.csv';
data_association = 'data_association.csv';

datetimeAMS = @(unix_time) subsasgn( datetime(unix_time, 'ConvertFrom', 'posixtime' , 'TimeZone', 'UTC'), struct('type', '.', 'subs', {'TimeZone'}), 'Europe/Amsterdam');
            
files = dir(fullfile(path, '*mainbrain'));
opts = delimitedTextImportOptions;
opts.VariableNames = {'camn', 'frame', 'timestamp', 'cam_received_timestamp', 'x', 'y', 'area', 'slope', 'eccentricity', 'frame_pt_idx', 'cur_val', 'mean_val', 'sumsqf_val'};
opts.SelectedVariableNames = {'frame'};
for ct=1:length(files)
    disp(['Reading ' files(ct).name]);
    data2d = readmatrix(fullfile(files(ct).folder, files(ct).name, file2d));
    frame = data2d(:,2);
    sync_triggers = find(diff(frame)<-100);
    
    if isempty(sync_triggers)
        disp(['  No sync triggers in this file. Copying as it is ...']);
        mkdir(fullfile(outputPath, files(ct).name));
        insidefiles = dir(fullfile(files(ct).folder, files(ct).name));
        insidefiles = insidefiles(~ismember({insidefiles.name},{'.','..'}));
        for ct2=1:length(insidefiles)
            if ~strcmpi(insidefiles(ct2).name,file3d) && ~strcmpi(insidefiles(ct2).name,data_association) && ~insidefiles(ct2).isdir
                copyfile(fullfile(insidefiles(ct2).folder, insidefiles(ct2).name), fullfile(outputPath, files(ct).name));
            elseif insidefiles(ct2).isdir && strcmpi(insidefiles(ct2).name,'images')
                mkdir(fullfile(outputPath, files(ct).name, 'images'));
                copyfile(fullfile(insidefiles(ct2).folder, insidefiles(ct2).name), fullfile(outputPath, files(ct).name, 'images'));
            end
        end
    else
        sync_triggers = [0; sync_triggers; size(frame,1)];
        for ct1=1:length(sync_triggers)-1
            syncStart = datetimeAMS(data2d(sync_triggers(ct1)+1,4));
            syncEnd = datetimeAMS(data2d(sync_triggers(ct1+1),4));
            if syncEnd-syncStart > duration(0,1,0) % copy only if duration between two sync triggers is more than 1 minute
            
                % make directory(or file in flydra terminology)
%                 unix_time = data2d(sync_triggers(ct1)+1,4);
%                 folderName = [datestr(datetimeAMS(unix_time), 'yyyymmdd_hhMMss') '.mainbrain'];
                folderName = [datestr(syncStart, 'yyyymmdd_hhMMss') '.mainbrain'];
                disp(['  Creating ' folderName]);
                mkdir(fullfile(outputPath, folderName));
                % write data2d_distorted.csv file
                writecell(opts.VariableNames, fullfile(outputPath, folderName, file2d));
                dlmwrite(fullfile(outputPath, folderName, file2d), data2d(sync_triggers(ct1)+1:sync_triggers(ct1+1),:), '-append', 'precision', 20);
                % Copy other files
                insidefiles = dir(fullfile(files(ct).folder, files(ct).name));
                insidefiles = insidefiles(~ismember({insidefiles.name},{'.','..'}));
                for ct2=1:length(insidefiles)
                    if ~strcmpi(insidefiles(ct2).name,file3d) && ~strcmpi(insidefiles(ct2).name,file2d) && ~strcmpi(insidefiles(ct2).name,data_association) && ~insidefiles(ct2).isdir
                        copyfile(fullfile(insidefiles(ct2).folder, insidefiles(ct2).name), fullfile(outputPath, folderName));
                    elseif insidefiles(ct2).isdir && strcmpi(insidefiles(ct2).name,'images')
                        mkdir(fullfile(outputPath, folderName, 'images'));
                        copyfile(fullfile(insidefiles(ct2).folder, insidefiles(ct2).name), fullfile(outputPath, folderName, 'images'));
                    end
                end
            
            end
        end
    end 
   
end