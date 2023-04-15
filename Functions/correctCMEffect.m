function [ v_and_i_timetrace_CM_corrected ] = correctCMEffect( v_and_i_timetrace )
%EXTRACTCMEFFECT Removes the "capacitor-discharging" effect from the current timetrace

time_axis = v_and_i_timetrace.*(0:length(v_and_i_timetrace.V_timetrace)-1);
I_timetrace = v_and_i_timetrace.I_timetrace;
V_timetrace = v_and_i_timetrace.V_timetrace;


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
B = fitresult.B
%Decay constant:
tau = fitresult.tau

I_time_decay = B*exp(-time_axis./tau);

I_model = I_time_decay + A*sin(omega.*time_axis + phi);

%Correct the signal (remove the decay-over-time effect)
I_timetrace_corrected = I_timetrace - I_time_decay' + B;

iv_timeseries_CM_corrected = iv_timeseries;
iv_timeseries_CM_corrected.I_timetrace = I_timetrace_corrected;

end

