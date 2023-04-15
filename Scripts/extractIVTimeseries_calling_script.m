
addpath('..\Functions\');

%Read the plasma measurement
plasma_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2018_08_24_Plasma_Profile_100kHz_10V_1_rep_CompleteProfile_ECRH_On_E-Guns_On_X-Drive_On_1500V_Rec-Drive_On_1000V\'
filename = 'meas_0290.h5';
settings.voltage_sign = +1;
settings.current_sign = -1;
settings.voltage_deamplification = 200;
v_and_i_timetrace_plasma = readHD5Measurement(plasma_source_folder,filename,settings);

%Read the vacuum measurement
vacuum_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2018_08_22_Vacuum_Langmuir_CurrentMonitor_Central_Position_100kHz_10V_DeamplifiedBy200\'
filename = 'meas_0000.h5';
settings.voltage_sign = +1;
settings.current_sign = -1;
settings.voltage_deamplification = 200;
v_and_i_timetrace_vacuum = readHD5Measurement(vacuum_source_folder,filename,settings);


% %Read the plasma measurement
% plasma_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_01_16_ECRH_Langmuir_CentralPosition_10kHz_10V_DeamplifiedBy200\'
% settings.voltage_sign = +1;
% settings.current_sign = -1;
% settings.voltage_deamplification = 200;
% 
% %Number of repetition to be averaged
% N = 10;
% 
% %Array of measurements
% vi_tt_plas_arr = cell(1,N);
% for n = 1:N
%     filename = strcat('meas_',num2str(n-1,'%04.f'),'.h5');
%     vi_tt_plas_arr{n}  = readHD5Measurement(plasma_source_folder,filename,settings);
% end
% 
% %Repetition-averaged measurement
% vi_tt = struct;
% vi_tt.time_step = vi_tt_plas_arr{1}.time_step;
% vi_tt.date = vi_tt_plas_arr{1}.date;
% vi_tt.time = vi_tt_plas_arr{1}.time;
% vi_tt.V_ADC_step = vi_tt_plas_arr{1}.V_ADC_step;
% vi_tt.I_ADC_step = vi_tt_plas_arr{1}.I_ADC_step;
% vi_tt.sampling_freq = vi_tt_plas_arr{1}.sampling_freq;
% vi_tt.V_timetrace = zeros(size(vi_tt_plas_arr{1}.V_timetrace));
% vi_tt.I_timetrace = zeros(size(vi_tt_plas_arr{1}.I_timetrace));
% 
% %Averaging of signals----
% for n = 1:N
%     vi_tt.V_timetrace = vi_tt.V_timetrace + vi_tt_plas_arr{n}.V_timetrace;
%     vi_tt.I_timetrace = vi_tt.I_timetrace + vi_tt_plas_arr{n}.I_timetrace;
% end
% vi_tt.V_timetrace = vi_tt.V_timetrace/N;
% vi_tt.I_timetrace = vi_tt.I_timetrace/N;
% %------------------------
% 
% time_axis = vi_tt.time_step.*(0:length(vi_tt.V_timetrace)-1);
% I_timetrace = vi_tt.I_timetrace;
% V_timetrace = vi_tt.V_timetrace;
% 
% 
% %Correct for CM-effect-----------------------------------------------------
% [xData, yData] = prepareCurveData( time_axis, I_timetrace );
% 
% % Set up fittype and options.
% ft = fittype( 'B*exp(-t/tau) + A*sin(omega*t+phi)', 'independent', 't', 'dependent', 'I' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.Lower = [0.1 0.08 62000 0 0.0005];
% opts.StartPoint = [0.17 0.1 62800 0.970592781760616 0.001];
% opts.Upper = [0.2 0.12 63000 6.28 0.002];
% 
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
% A = fitresult.A;
% omega = fitresult.omega;
% phi = fitresult.phi;
% %Initial voltage offset:
% B = fitresult.B;
% %Decay constant:
% tau = fitresult.tau;
% I_time_decay = B*exp(-time_axis./tau);
% I_model = I_time_decay + A*sin(omega.*time_axis + phi);
% %Correct the signal (remove the decay-over-time effect)
% I_timetrace_corrected = I_timetrace - I_time_decay' + B;
% %Do FFT operations
% signal = I_timetrace_corrected;
% Fs = vi_tt.sampling_freq;
% %--------------------------------------------------------------------------



[ return_status, iv_timeseries ] = extractIVTimeseries( v_and_i_timetrace_plasma, v_and_i_timetrace_vacuum );
%[ return_status, iv_timeseries ] = extractIVTimeseries(vi_tt);

nc_selected = [6,9,12,18,21,24];

time_axis = 1e6.*(iv_timeseries.total_time_axis_rise + iv_timeseries.total_time_axis_fall)/2;

close all
n_sub_plot = 1;
for nc = nc_selected
    
    V_axis_rise = iv_timeseries.rise(nc).V_axis;
    I_of_V_rise = iv_timeseries.rise(nc).I_of_V;

    V_axis_fall = iv_timeseries.fall(nc).V_axis;
    I_of_V_fall = iv_timeseries.fall(nc).I_of_V;
    
    delta_I = iv_timeseries.I_ADC_step;
    %subplot(2,3,n_sub_plot);
    figure
    set(gcf,'Renderer','painters');
    errorbar(V_axis_rise,I_of_V_rise,ones(length(V_axis_rise),1).*delta_I,'or');
    hold on;
    errorbar(V_axis_fall,I_of_V_fall,ones(length(V_axis_fall),1).*delta_I,'ok');
    min_length = min([length(I_of_V_rise),length(I_of_V_fall)]);
    %I_of_V_average = (I_of_V_rise(1:min_length) + I_of_V_fall(1:min_length))./2;
    %plot(V_axis_rise(1:min_length),I_of_V_average,'-b');
    hold off;
    title(['IV characteristic at time = ',num2str(time_axis(nc)),' µs']);
%     xlim([-80,+80]);
%     ylim([-5e-3,2e-2]);
    xlabel('Bias voltage [V]');
    ylabel('Probe current [A]');
    set(gca,'fontsize',20);
    set(gcf, 'Position', [100, 100, 900, 600]);
    
    n_sub_plot = n_sub_plot + 1;

end