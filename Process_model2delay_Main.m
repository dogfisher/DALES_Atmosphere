clear,clc
close all;
fclose all;
%Script for generate time series atmospheric interferogram
%Written by Fengming Hu in TJU
%Please change the input_path to the data folder in your computer
input_path='E:/TUD_Work_Data/data2/';
%===========
%Step 1: merge all nc files and calculate refactivity
%output files are N_xyz.nc (Cordinates) and N_epoch*.nc (Refractivity) files
% in folder './N_delay_out/'
process_ncfiles(input_path);
N_path=[input_path,'N_delay_out/'];

%===========
%Step 2: Set the orbit and ifg parameters
%If this step is running already, just load the SAR_geo.mat file in folder
%'./N_delay_out/'
Orb_setting(N_path);
load([N_path,'SAR_geo.mat']);

%===========
%Step 3: Resample the refractivity model to slice
Resample_slice(N_path,orb_slc);
ReN_path=[N_path,'Slant_delay_out/'];

%===========
%Step 4: Calculate the slant delay for all epoch
Slant_Delay_Calc(ReN_path,orb_slc);

%===========
%Step 5: Calculate the atmospheric interferograms; Single master or daisy
%chain and output results in figures and videos
atmos_interf_ts(ReN_path,orb_slc);
%Output files are saved in folder 'ReN_path/TS_phase_single and ReN_path/TS_phase_daisy_chain'
%Total_delay results (*.fig, *.tif and *.avi) are saved in folder 'ReN_path'
