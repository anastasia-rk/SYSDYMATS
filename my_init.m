clc; clear; close all;
%% Figure environment setup
set(0,'DefaultFigureWindowStyle','docked');
set(0,'defaultfigurecolor',[1 1 1]);
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultColorbarTickLabelInterpreter','latex');
set(0, 'DefaultAxesFontWeight', 'normal','DefaultAxesFontSize', 16);

%% Custom colormap
my_map = cutsom_colormap(256);

