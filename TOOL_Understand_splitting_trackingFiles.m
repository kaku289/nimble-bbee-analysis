clc; clear; close all;
addpath('/home/reken001/Pulkit/lib/flymovieformat');
addpath('/home/reken001/Pulkit/lib/DrosteEffect-BrewerMap-221b913');

% %% video file
% file = '/media/reken001/Disk_09/unsteady_wind_experiments/Videos/2019/07/28/091348-105845-Basler_22549584.fmf';
% file = '/media/reken001/Disk_09/unsteady_wind_experiments/Videos/2019/07/28/092136-187644-Basler_22549584.fmf';
% [video_data, timestamps] = fmf_read(file);

file3d = 'kalman_estimates.csv';
file2d = 'data2d_distorted.csv';
data_association = 'data_association.csv';

%% tracking file without any sync triggers (e.g., from steady wind experiments)
path = '/media/reken001/Disk_09/unsteady_wind_experiments/test/Tracking/20190715_110317.mainbrain/';
Data2d = readmatrix([path file2d]);

%% tracking file with sync triggers (e.g., from oscillating wind experiments)
path = '/media/reken001/Disk_09/unsteady_wind_experiments/test/Tracking/20190727_074310.mainbrain/';
Data2d_syncTrigs = readmatrix([path file2d]);

%% Plot
figure; hold on;
plot(Data2d(:,2),'b.');
plot(Data2d_syncTrigs(:,2),'r.');

%%
Data3d = importdata([path file3d]);
data3d = Data3d.data;
clear Data3d;

Data2d = importdata([path file2d]);
data2d = Data2d.data;
clear Data2d;

Data_association = importdata([path data_association]);
data_association = Data_association.data;
clear Data_association;