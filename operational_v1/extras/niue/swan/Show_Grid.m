clear all
close all
filegrid=['fort' '.14'];%% Grid file
[fem,elebnd]=read_adcirc_mesh(filegrid);
fsz = 12; % default font size
bgc = [1 1 1]; % default background color
figure
hold on
trisurf(fem.e,fem.x,fem.y,fem.z);