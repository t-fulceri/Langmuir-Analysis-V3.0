addpath('..\Functions\');

plasma_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2018_08_24_Plasma_Profile_100kHz_10V_1_rep_CompleteProfile_ECRH_On_E-Guns_On_X-Drive_On_1500V_Rec-Drive_On_1000V\'
filename = 'meas_0290.h5';
settings.voltage_sign = +1;
settings.current_sign = -1;
settings.voltage_deamplification = 200;
v_and_i_timetrace = readHD5Measurement(plasma_source_folder,filename,settings);

% vacuum_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2018_08_22_Vacuum_Langmuir_CurrentMonitor_Central_Position_100kHz_10V_DeamplifiedBy200\'
% filename = 'meas_0000.h5';
% settings.voltage_sign = +1;
% settings.current_sign = -1;
% settings.voltage_deamplification = 200;
% v_and_i_timetrace = readHD5Measurement(vacuum_source_folder,filename,settings);

% vacuum_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2018_08_27_Vacuum_Effects\f_100khz\'
% filename = 'meas_0000.h5';
% settings.voltage_sign = +1;
% settings.current_sign = -1;
% settings.voltage_deamplification = 200;
% v_and_i_timetrace = readHD5Measurement(vacuum_source_folder,filename,settings);

% vacuum_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2018_08_21_Plasma_Profile_01kHz_10V_1_rep_CurrentMonitor_CompleteProfile_ECRH_Only_DeamplifiedBy200\'
% filename = 'meas_0290.h5';
% settings.voltage_sign = +1;
% settings.current_sign = -1;
% settings.voltage_deamplification = 200;
% v_and_i_timetrace = readHD5Measurement(vacuum_source_folder,filename,settings);

V = v_and_i_timetrace.V_timetrace;
I = v_and_i_timetrace.I_timetrace;
ts = v_and_i_timetrace.time_step;
time_axis = (0:(length(V)-1)).*ts;

% % conversion to ms
% time_axis = 1e3.*time_axis;
% t_ax_string = 't [ms]';

%conversion to µs
time_axis = 1e6.*time_axis;
t_ax_string = 't [µs]';

close all
figure
set(gcf,'Renderer','painters');
movegui(gcf,'center');

subplot(2,1,1);
plot(time_axis,V);
set(gca,'fontsize',20);
title('Bias voltage over time');
xlabel(t_ax_string);
ylabel('V [V]');

subplot(2,1,2);
plot(time_axis,I);
set(gca,'fontsize',20);
title('Probe current over time');
xlabel(t_ax_string);
ylabel('I [A]');

set(gcf, 'Position', [100, 100, 900, 900]);