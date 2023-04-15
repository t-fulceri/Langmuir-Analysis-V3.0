addpath('..\Functions\');

%Access grid scan directory

%Define the source folder (where the raw measurement data is taken)
%source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\Reconstructed_Profiles\2018_08_21_Plasma_Profile_01kHz_10V_1_rep_CurrentMonitor_CompleteProfile_ECRH_Only_DeamplifiedBy200\';
source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\Reconstructed_Profiles\2019_01_16_ECRH_Langmuir_CompleteProfile_10kHz_10V_DeamplifiedBy200\';
temp = load([source_folder,'iv_timeseries_cell_array']);
iv_timeseries_cell_array = temp.iv_timeseries_cell_array;

%Define the targed folder (where the data will be saved)
%target_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\Reconstructed_Profiles\2018_08_21_Plasma_Profile_01kHz_10V_1_rep_CurrentMonitor_CompleteProfile_ECRH_Only_DeamplifiedBy200\';
target_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\Reconstructed_Profiles\2019_01_16_ECRH_Langmuir_CompleteProfile_10kHz_10V_DeamplifiedBy200\';

[ iv_time_average_cell_array ] = averageIVTimeseriesCellArray( iv_timeseries_cell_array );

%Save all the iv_timeseries
save([target_folder,'iv_time_average_cell_array'],'iv_time_average_cell_array');