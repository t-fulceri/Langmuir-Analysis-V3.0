clear all;
close all;

addpath('..\Functions\');


source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2018_08_09_Vacuum_LC measurements_CurrentMonitor_050kHz_4V\'
filename_cell_array = cell(1,1);
filename_cell_array{1} = 'meas_0000.h5';
%filename_cell_array{2} = 'meas_0001.h5';
% filename_cell_array{3} = 'meas_0307.h5';
% filename_cell_array{4} = 'meas_0308.h5';
% filename_cell_array{5} = 'meas_0309.h5';

[ fitted_vac_lc_effect ] = extractLC_VacuumEffect( source_folder, filename_cell_array );