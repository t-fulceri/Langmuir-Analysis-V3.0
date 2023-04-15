close all;
clear all;

source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2018_08_21_Plasma_Profile_01kHz_10V_1_rep_CurrentMonitor_CompleteProfile_ECRH_Only_DeamplifiedBy200\';
filename = 'meas_0109.h5';
voltage_sign = +1;
current_sign = -1;
voltage_deamplification = 200;

full_file_path = strcat(source_folder,filename);
%Check metadata
metadata = h5info(full_file_path);
%Extract timestep from metadata [s]
time_step = metadata.Groups.Attributes(2,1).Value(1);

%Extract date and time from metadata
date = metadata.Datasets.Attributes(1,1).Value;
time = metadata.Datasets.Attributes(2,1).Value;

%Read data
data = h5read(full_file_path,'/dataset');
%DAQ delay in s
%DAQ_delay = 20e-3;
%Transpose
data = data';
%Extract voltage
voltage = data(:,1);
%Create time axis
time_axis = time_step.*(1:length(voltage));
%Optional inversion of voltage signal (depends on the data source)
voltage = voltage_deamplification.*voltage_sign.*voltage;
%Extract current
current = data(:,2);
%Optional inversion of current signal (depends on the data source)
current = current_sign.*current;

figure
set(gcf,'Renderer','painters');
movegui(gcf,'west');
subplot(2,1,1);
plot(time_axis,voltage);
subplot(2,1,2);
plot(time_axis,current);

figure
set(gcf,'Renderer','painters');
movegui(gcf,'east');
plot(voltage,current,'.');

signal = current;
Fs = 1/time_step;

% Length of signal
L = length(signal);             

%Prepare the array of frequencies
NFFT = 2^nextpow2(L); % Next power of 2 from length of y

% Use only positive frequencies as reference
freq = Fs/2*linspace(0,1,NFFT/2+1);

%Perform the Fast Fourier Transform using the array of frequencies
fft_signal = fft(signal,NFFT)/L;

%Consider only positive frequencies
fft_positive = fft_signal(1:NFFT/2+1);

%For every +/- frequency pair, sum their amplitudes
double_abs_fft_signal = 2*abs(fft_positive);

%Look for maximum amplitude in the spectrum
[M,M_ind] = max(double_abs_fft_signal);
%Find corresponding frequency
main_freq = freq(M_ind);
%Find corresponding phase
main_phase = angle(fft_positive(M_ind));

%PLOTTING (JUST FOR DEBUG)
% Plot single-sided amplitude distribution in frequency space
f_fft_signal = figure;
movegui(f_fft_signal,'north');
set(gcf,'Renderer','painters');
ax_fft_signal = axes;
double_abs_fft_signal_normalized = double_abs_fft_signal./M;
loglog(ax_fft_signal,freq,double_abs_fft_signal_normalized,'.-');
title(ax_fft_signal,'Normalized amplitude distribution (AKA spectrum)');
xlabel(ax_fft_signal,'Frequency [Hz]');
ylabel(ax_fft_signal,'Normalized Amplitude');  

% Plot single-sided phase distribution in frequency space
f_phase = figure;
movegui(f_phase,'northeast');
set(gcf,'Renderer','painters');
ax_phase = axes;
semilogx(ax_phase,freq,angle(fft_positive),'.-');
title(ax_phase,'Phase distribution');
xlabel(ax_phase,'Frequency [Hz]');
ylabel(ax_phase,'Phase');

% Plot original signal plus extracted main component
f_supimp = figure;
movegui(f_supimp,'center');
set(gcf,'Renderer','painters');
ax_supimp = axes;
hold on;
plot(signal);
time_step = (1/Fs);
t = time_step*linspace(0,L-1,L);
main_component = M*cos(2*pi*main_freq*t + main_phase);
plot(main_component);
hold off;
title(ax_supimp,'Original signal superimposed with main Fourier component');
xlabel(ax_supimp,'Time [timesteps]');
ylabel(ax_supimp,'Phase');