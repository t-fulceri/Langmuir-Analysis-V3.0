%This script analyzes a Langmuir measurement of the ECRH plasma in the
%central spot

%Tiziano Fulceri 2019-01-25

%Optimized to work with langmuir measurements with the following parameters:
% f = 10 kHz
% V_max = 100 V
% time interval = 10 ms
% Probe position = (148, 115) (Langmuir probe at the center)
% ECRH plasma ON

clear all
close all

addpath('..\Functions\');

%Read the plasma measurement
plasma_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_01_16_ECRH_Langmuir_CentralPosition_10kHz_10V_DeamplifiedBy200\'
settings.voltage_sign = +1;
settings.current_sign = -1;
settings.voltage_deamplification = 1;

%Number of repetition to be averaged
N = 10;

%Array of measurements
vi_tt_plas_arr = cell(1,N);
for n = 1:N
    filename = strcat('meas_',num2str(n-1,'%04.f'),'.h5');
    vi_tt_plas_arr{n}  = readHD5Measurement(plasma_source_folder,filename,settings);
end

%Repetition-averaged measurement
vi_tt = struct;
vi_tt.time_step = vi_tt_plas_arr{1}.time_step;
vi_tt.date = vi_tt_plas_arr{1}.date;
vi_tt.time = vi_tt_plas_arr{1}.time;
vi_tt.V_ADC_step = vi_tt_plas_arr{1}.V_ADC_step;
vi_tt.I_ADC_step = vi_tt_plas_arr{1}.I_ADC_step;
vi_tt.sampling_freq = vi_tt_plas_arr{1}.sampling_freq;
vi_tt.V_timetrace = zeros(size(vi_tt_plas_arr{1}.V_timetrace));
vi_tt.I_timetrace = zeros(size(vi_tt_plas_arr{1}.I_timetrace));

%Averaging of signals----
for n = 1:N
    vi_tt.V_timetrace = vi_tt.V_timetrace + vi_tt_plas_arr{n}.V_timetrace;
    vi_tt.I_timetrace = vi_tt.I_timetrace + vi_tt_plas_arr{n}.I_timetrace;
end
vi_tt.V_timetrace = vi_tt.V_timetrace/N;
vi_tt.I_timetrace = vi_tt.I_timetrace/N;
%------------------------

time_axis = vi_tt.time_step.*(0:length(vi_tt.V_timetrace)-1);
I_timetrace = vi_tt.I_timetrace;
V_timetrace = vi_tt.V_timetrace;


%Correct for CM-effect-----------------------------------------------------
[xData, yData] = prepareCurveData( time_axis, I_timetrace );

% Set up fittype and options.
ft = fittype( 'B*exp(-t/tau) + A*sin(omega*t+phi)', 'independent', 't', 'dependent', 'I' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0.1 0.08 62000 0 0.0005];
opts.StartPoint = [0.17 0.1 62800 0.970592781760616 0.001];
opts.Upper = [0.2 0.12 63000 6.28 0.002];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
A = fitresult.A;
omega = fitresult.omega;
phi = fitresult.phi;
%Initial voltage offset:
B = fitresult.B;
%Decay constant:
tau = fitresult.tau;
I_time_decay = B*exp(-time_axis./tau);
I_model = I_time_decay + A*sin(omega.*time_axis + phi);
%Correct the signal (remove the decay-over-time effect)
I_timetrace_corrected = I_timetrace - I_time_decay' + B;
%--------------------------------------------------------------------------

vi_tt_corrected = vi_tt;
vi_tt_corrected.I_timetrace = I_timetrace_corrected;

[ return_status, iv_timeseries ] = extractIVTimeseries(vi_tt_corrected)
[return_status, plasma_param_timeseries] = extractPlasmaParamIVTimeseries(iv_timeseries)

T_e = plasma_param_timeseries.mean.T_e;
I_sat_i = plasma_param_timeseries.mean.I_sat_i;
V_plasma = plasma_param_timeseries.mean.V_plasma;

n_max = 20;
T_e_mean = mean(T_e(1:n_max))
I_sat_i_mean = mean(I_sat_i(1:n_max))
V_plasma_mean = mean(V_plasma(1:n_max))


