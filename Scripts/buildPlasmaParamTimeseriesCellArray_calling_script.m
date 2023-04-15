
addpath('..\Functions\');

%Source folder for the plasma parameters timeseries list
source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\Reconstructed_Profiles\2018_08_21_Plasma_Profile_01kHz_10V_1_rep_CurrentMonitor_CompleteProfile_ECRH_Only_DeamplifiedBy200\';
NX = 23;
NY = 18;
NR = 1;

load([source_folder,'iv_timeseries_cell_array.mat']);
[ plasma_param_timeseries_cell_array ] = buildPlasmaParamTimerseriesCellArray( [source_folder,'plasma_parameters_timeseries_list\'], NX, NY, NR);