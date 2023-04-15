%Thi script extracts the parameters which model the current monitor
%"capacitor charging"-like effect

%Tiziano Fulceri 2019-01-19

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
settings.voltage_deamplification = 200;

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

% Plot fit with data.
%h = plot( fitresult, xData, yData );
%legend( h, 'I_timetrace vs. time_axis', 'untitled fit 1', 'Location', 'NorthEast' );
%grid on

A = fitresult.A;
omega = fitresult.omega;
phi = fitresult.phi;
%Initial voltage offset:
B = fitresult.B
%Decay constant:
tau = fitresult.tau

I_time_decay = B*exp(-time_axis./tau);

I_model = I_time_decay + A*sin(omega.*time_axis + phi);

figure
set(gcf,'Renderer','painters');
subplot(2,1,1);
hold on
plot(time_axis.*1e6,I_timetrace,'.');
plot(time_axis.*1e6,I_model);
plot(time_axis.*1e6,I_time_decay./B);
title('Current timetrace');
xlabel('Time (탎)');
ylabel('Current (A)');
hold off

%Correct the signal (remove the decay-over-time effect)
I_timetrace_corrected = I_timetrace - I_time_decay' + B;

subplot(2,1,2);
set(gcf,'Renderer','painters');
movegui(gcf,'northwest');
plot(time_axis.*1e6,I_timetrace_corrected);
title('CM-Corrected current timetrace');
xlabel('Time (탎)');
ylabel('Current (A)');

figure
set(gcf,'Renderer','painters');
subplot(2,1,1);
plot(time_axis.*1e6,V_timetrace);
title('Voltage timetrace');
xlabel('Time (탎)');
ylabel('Voltage (V)');

subplot(2,1,2);
set(gcf,'Renderer','painters');
movegui(gcf,'northwest');
plot(time_axis.*1e6,I_timetrace_corrected);
title('CM-Corrected current timetrace');
xlabel('Time (탎)');
ylabel('Current (A)');


figure
set(gcf,'Renderer','painters');
movegui(gcf,'northeast');
plot(time_axis.*1e6,I_timetrace_corrected);
title('CM-Corrected current timetrace');
xlabel('Time (탎)');
ylabel('Current (A)');
ylim([0 0.3]);

figure
set(gcf,'Renderer','painters');
movegui(gcf,'southeast');
plot(time_axis.*1e6,I_timetrace_corrected);
title('CM-Corrected current timetrace');
xlabel('Time (탎)');
ylabel('Current (A)');
ylim([-1e-2 0]);



%Do FFT operations
signal = I_timetrace_corrected;
Fs = vi_tt.sampling_freq;

% Length of signal
L = length(signal);             

%Prepare the array of frequencies:

%Total number of frequencies
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
% Use only positive frequencies as reference
freq_pos = Fs/2*linspace(0,1,NFFT/2+1); 
freq_neg = fliplr(-freq_pos(2:end-1));
freq_sym = [freq_pos,freq_neg];

%Perform the Fast Fourier Transform using the array of frequencies
fft_signal = fft(signal,NFFT)/L;

%Consider only positive frequencies
fft_positive = fft_signal(1:NFFT/2+1);

%For every +/- frequency pair, sum their amplitudes
double_abs_fft_signal = 2*abs(fft_positive);

%Look for maximum amplitude in the spectrum (exclude DC component)
[M,M_ind] = max(double_abs_fft_signal(2:end));
%Find corresponding frequency
main_freq = freq_pos(M_ind);
%Find corresponding phase
main_phase = angle(fft_positive(M_ind));

%PLOTTING
% Plot single-sided amplitude distribution in frequency space
f_fft_signal = figure;
set(gcf,'Renderer','painters');
movegui(f_fft_signal,'north');
ax_fft_signal = axes;
double_abs_fft_signal_normalized = double_abs_fft_signal./M;
loglog(ax_fft_signal,freq_sym,2*abs(fft_signal)./M,'-');
%loglog(ax_fft_signal,freq_pos,double_abs_fft_signal_normalized,'.-');
title(ax_fft_signal,'Normalized amplitude distribution (AKA spectrum)');
xlabel(ax_fft_signal,'Frequency [Hz]');
ylabel(ax_fft_signal,'Normalized Amplitude');  

% % Plot single-sided phase distribution in frequency space
% f_phase = figure;
% set(gcf,'Renderer','painters');
% movegui(f_phase,'northeast');
% ax_phase = axes;
% semilogx(ax_phase,freq,angle(fft_positive),'.-');
% title(ax_phase,'Phase distribution');
% xlabel(ax_phase,'Frequency [Hz]');
% ylabel(ax_phase,'Phase');


%Remove frequencies lower than the peak frequency
%Leave DC component intact
fft_signal_after_HP = fft_signal;
fft_signal_after_HP(2:M_ind-1) = 0;
fft_signal_after_HP(end-M_ind+2:end) = 0;

%Consider only positive frequencies
%fft_positive_after_HP = fft_signal_after_HP(1:NFFT/2+1);

%For every +/- frequency pair, sum their amplitudes
%double_abs_fft_signal_after_HP = 2*abs(fft_positive_after_HP);

% %Look for maximum amplitude in the spectrum (exclude DC component)
% [M,M_ind] = max(double_abs_fft_signal(2:end));
% %Find corresponding frequency
% main_freq = freq(M_ind);
% %Find corresponding phase
% main_phase = angle(fft_positive(M_ind));

% Plot single-sided amplitude distribution in frequency space
f_fft_signal_after_HP = figure;
set(gcf,'Renderer','painters');
movegui(f_fft_signal_after_HP,'south');
ax_fft_signal_after_HP = axes;
loglog(ax_fft_signal_after_HP,freq_sym,2*abs(fft_signal_after_HP)./M,'.-');
%plot(ax_fft_signal_after_HP,abs(fft_signal_after_HP),'.-');
title(ax_fft_signal_after_HP,'Normalized amplitude distribution (AKA spectrum)');
xlabel(ax_fft_signal_after_HP,'Frequency [Hz]');
ylabel(ax_fft_signal_after_HP,'Normalized Amplitude');

% % Plot single-sided amplitude distribution in frequency space
% f_fft_signal_after_HP = figure;
% set(gcf,'Renderer','painters');
% movegui(f_fft_signal_after_HP,'south');
% ax_fft_signal_after_HP = axes;
% double_abs_fft_signal_normalized_after_HP = double_abs_fft_signal_after_HP./M;
% loglog(ax_fft_signal_after_HP,freq_pos,double_abs_fft_signal_normalized_after_HP,'.-');
% title(ax_fft_signal_after_HP,'Normalized amplitude distribution (AKA spectrum)');
% xlabel(ax_fft_signal_after_HP,'Frequency [Hz]');
% ylabel(ax_fft_signal_after_HP,'Normalized Amplitude');

%return to time-domain
signal_after_HP = ifft(fft_signal_after_HP,L);

% % Plot HP filtered signal
% f_signal_after_HP = figure;
% set(gcf,'Renderer','painters');
% movegui(f_signal_after_HP,'south');
% ax_signal_after_HP = axes;
% plot(time_axis.*1e6,signal_after_HP);
% title(ax_signal_after_HP,'High-Pass-filtered corrected timetrace');
% xlabel(ax_signal_after_HP,'Time (탎)');
% ylabel(ax_signal_after_HP,'Current (A)');
