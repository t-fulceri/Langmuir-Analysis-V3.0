clear all
close all

addpath('..\Functions\');

%Access grid scan directory

% %Plasma measurement folder
% plasma_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2018_08_24_Plasma_Profile_100kHz_10V_1_rep_CompleteProfile_ECRH_On_E-Guns_On_X-Drive_On_1500V_Rec-Drive_On_1000V\'
% %Plasma measurement settings
% settings.voltage_sign = +1;
% settings.current_sign = -1;
% settings.voltage_deamplification = 200;
% 
% rep_per_point = 1;
% 
% %Read the vacuum measurement
% vacuum_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2018_08_22_Vacuum_Langmuir_CurrentMonitor_Central_Position_100kHz_10V_DeamplifiedBy200\'
% filename = 'meas_0000.h5';
% settings.voltage_sign = +1;
% settings.current_sign = -1;
% settings.voltage_deamplification = 200;
% v_and_i_timetrace_vacuum = readHD5Measurement(vacuum_source_folder,filename,settings);
% 
% %Define the targed folder (where the data will be saved)
% target_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\Reconstructed_Profiles\2018_08_24_Plasma_Profile_100kHz_10V_1_rep_CompleteProfile_ECRH_On_E-Guns_On_X-Drive_On_1500V_Rec-Drive_On_1000V\';


%Plasma measurement folder
plasma_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_01_16_ECRH_Langmuir_CompleteProfile_10kHz_10V_DeamplifiedBy200\oscilloscope\'
%Plasma measurement settings
settings.voltage_sign = +1;
settings.current_sign = -1;
settings.voltage_deamplification = 1;

rep_per_point = 1;

%Read the vacuum measurement
vacuum_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_01_28_Vacuum_Effects\010khz\'
filename = 'meas_0000.h5';
settings.voltage_sign = +1;
settings.current_sign = -1;
settings.voltage_deamplification = 1;
v_and_i_timetrace_vacuum = readHD5Measurement(vacuum_source_folder,filename,settings);

%Define the targed folder (where the data will be saved)
target_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\Reconstructed_Profiles\2019_01_16_ECRH_Langmuir_CompleteProfile_10kHz_10V_DeamplifiedBy200\';
map_file = 'map.map'
full_file_path = strcat(plasma_source_folder,map_file);

%Raster step
raster_step = 10;

%Read the map file
B = tdfread(full_file_path);
ind = B.x0x25_idx;
x_pos = B.x;
y_pos = B.y;
ind = ind(:,1);
x_pos = x_pos(:,1);
y_pos = y_pos(:,1);
%Number of points

%CORRECTION
x_pos_mod = 10*floor(x_pos/10);
x_pos = x_pos_mod;
%----------

NX = (max(x_pos) - min(x_pos))/raster_step + 1;
NY = (max(y_pos) - min(y_pos))/raster_step + 1;

%Organize data in multidimensional array
%Measurement number 2Dspace*1Drepetition array
meas_number_STR = NaN(NX,NY,rep_per_point);

for k_x = 1:NX
    for k_y = 1:NY
        for k_r = 1:rep_per_point
            %file_name = strcat('meas_',num2str(k+(repetition-1),'%04.f'),'.h5');
            %file_names_STR(k_x,k_y,repetition) = file_name;
            k = k_x + (k_y-1)*NX;
            meas_number_STR(k_x,k_y,k_r) = rep_per_point*(k-1)+(k_r-1);
        end
    end
end

%Select sub-grid to process
NX_begin = 1;
NX_end = NX;
nX = NX_end - NX_begin + 1;

NY_begin = 1;
NY_end = NY;
nY = NY_end - NY_begin + 1;

rep_begin = 1;
rep_end = 5;
%Self correct:
if rep_end > rep_per_point
    rep_end = rep_per_point;
end
nr = rep_end - rep_begin +1;


%Multidimensional cell array to store all the iv_timeseries
iv_timeseries_cell_array = cell(NX,NY,rep_per_point);

% %Multidimensional cell array to store all the plasma_param_timeseries_data
% plasma_param_timeseries_cell_array = cell(NX,NY,rep_per_point);

fprintf('Start of plasma parameter extraction for the whole profile.\n');
for k_x = 1:nX
    for k_y = 1:nY
        for k_r = 1:nr
            fprintf(['X_index = ',num2str(k_x),'\n']);
            fprintf(['Y_index = ',num2str(k_y),'\n']);
            fprintf(['Repetition_index = ',num2str(k_r),'\n\n']);
            file_number = meas_number_STR(k_x+NX_begin-1,k_y+NY_begin-1,k_r + rep_begin-1);
            filename = strcat('meas_',num2str(file_number,'%04.f'),'.h5')
            try
                v_and_i_timetrace_plasma = readHD5Measurement(plasma_source_folder,filename,settings);
                [ return_status_read, iv_timeseries ] = extractIVTimeseries(v_and_i_timetrace_plasma,v_and_i_timetrace_vacuum);
                %Save the first iv_timeseries which can be extracted
                if return_status_read == 0
                    iv_timeseries_cell_array{k_x,k_y,k_r} = iv_timeseries;
                end
            catch err_read
                fprintf('IV Timeseries Read failed: jumping to next interation\n\n');
                continue
            end
            

        end
    end
end

%Save all the iv_timeseries
save([target_folder,'iv_timeseries_cell_array.mat'],'iv_timeseries_cell_array');
%load([target_folder,'iv_timeseries_cell_array.mat']);

%Prepare the folder for the plasma parameters timeseries
mkdir([target_folder,'plasma_parameters_timeseries_list']);

fprintf('Start of plasma parameter extraction for the whole profile.\n');
for k_x = 1:nX
    for k_y = 1:nY
        for k_r = 1:nr
            fprintf(['X_index = ',num2str(k_x),'\n']);
            fprintf(['Y_index = ',num2str(k_y),'\n']);
            fprintf(['Repetition_index = ',num2str(k_r),'\n\n']);
            try
                [ return_status_extract, plasma_param_timeseries ] = extractPlasmaParamIVTimeseries(iv_timeseries_cell_array{k_x,k_y,k_r});
                %Prepare the file name
                file_number = meas_number_STR(k_x+NX_begin-1,k_y+NY_begin-1,k_r + rep_begin-1);
                filename = strcat('plasma_param_timeseries_',num2str(file_number,'%04.f'),'.mat')
                save([target_folder,'plasma_parameters_timeseries_list\',filename],'plasma_param_timeseries');
                fprintf('Plasma Parameter Timeseries saved.\n');
            catch err_extract
                fprintf('Plasma Parameter extraction failed: jumping to next interation...\n');
                continue
            end
        end
    end
end
