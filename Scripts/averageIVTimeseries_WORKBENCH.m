%This script is to be used with stationary ECRH Plasma measurement data:
%it averages an IV timeseries over time giving a single averaged IV-curve
%as output

clear
close all

addpath('..\Functions\');

source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\Reconstructed_Profiles\2018_08_17_Plasma_Profile_500Hz_10V_1_rep_CurrentMonitor_CompleteProfile_ECRH_Only_DeamplifiedBy200\';

%Load the timeseries cell array data structure
load([source_folder, 'iv_timeseries_cell_array']);

%Select iv_timeseries data structure to be time-averaged
iv_timeseries = iv_timeseries_cell_array{11,9};

%Perform the avearge over all the cycles

iv_time_average.rise.V = cat(1,iv_timeseries.rise(:).V_axis);
iv_time_average.rise.I = cat(1,iv_timeseries.rise(:).I_of_V);

iv_time_average.fall.V = cat(1,iv_timeseries.fall(:).V_axis);
iv_time_average.fall.I = cat(1,iv_timeseries.fall(:).I_of_V);

%Reconstruct time-averaged IV-curve

rec_curve_out = reconstructCurve(iv_time_average.rise.V,iv_time_average.rise.I);
iv_time_average.rise.V_axis = rec_curve_out.x;
iv_time_average.rise.I_of_V = rec_curve_out.y_of_x;

rec_curve_out = reconstructCurve(iv_time_average.fall.V,iv_time_average.fall.I);
iv_time_average.fall.V_axis = rec_curve_out.x;
iv_time_average.fall.I_of_V = rec_curve_out.y_of_x;

%Plot time averaged IV-characteristic
plot(iv_time_average.rise.V_axis,iv_time_average.rise.I_of_V)
hold on
plot(iv_time_average.fall.V_axis,iv_time_average.fall.I_of_V)
