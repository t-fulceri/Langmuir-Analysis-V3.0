clear all
close all

addpath('..\Functions\');

%Access grid scan directory

%Define the source folder (where the raw measurement data is taken)
%source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2018_08_21_Plasma_Profile_01kHz_10V_1_rep_CurrentMonitor_CompleteProfile_ECRH_Only_DeamplifiedBy200\'
source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_01_16_ECRH_Langmuir_CompleteProfile_10kHz_10V_DeamplifiedBy200\oscilloscope\'
%Set the number of repetition per position
rep_per_point = 1;
%Raster step [mm]
raster_step = 10;

%Define the source folder for the LC vacuum effect (self-inductance of the measuring circuit)
%LC_vacuum_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2018_08_21_Vacuum_Langmuir_CurrentMonitor_Central_Position_01kHz_10V_DeamplifiedBy200\';
LC_vacuum_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_01_28_Vacuum_Effects\010khz\';

%Define the targed folder (where the data will be saved)
% target_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\Reconstructed_Profiles\2018_08_21_Plasma_Profile_01kHz_10V_1_rep_CurrentMonitor_CompleteProfile_ECRH_Only_DeamplifiedBy200\';
% voltage_deamplification = 200;
target_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\Reconstructed_Profiles\2019_01_16_ECRH_Langmuir_CompleteProfile_10kHz_10V_DeamplifiedBy200\';
voltage_deamplification = 1;

iv_timeseries_cell_array = buildIVTimeseriesCellArray( source_folder, rep_per_point, raster_step, LC_vacuum_source_folder, voltage_deamplification);

%Save all the iv_timeseries
save([target_folder,'iv_timeseries_cell_array.mat'],'iv_timeseries_cell_array');
