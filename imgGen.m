clc;clear;close all;

% data = csvread ('D:\OneDrive - Universiti Malaysia Pahang\Atik_Home\Data Files\fyp-velocity signals\Twist_3.csv');
% y = data(1:2500,:);

load('D:\OneDrive - ump.edu.my\Atik_Home\Data Files\nEEEMD\Data\130.mat');
y = X130_DE_time(1:6400,:);

dmin = min(y);
dmax = max(y);
scaled_data = (y - dmin)./(dmax - dmin) * 2 - 1;

scaled = reshape(y,[80,80]);
imshow(scaled);