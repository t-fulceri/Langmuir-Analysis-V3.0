

addpath('..\Functions\');

%Read the plasma measurement
plasma_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2018_08_24_Plasma_Profile_100kHz_10V_1_rep_CompleteProfile_ECRH_On_E-Guns_On_X-Drive_On_1500V_Rec-Drive_On_1000V\'
filename = 'meas_0198.h5';
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

%[ return_status, iv_timeseries ] = extractIVTimeseries( v_and_i_timetrace_plasma );
[ return_status, iv_timeseries ] = extractIVTimeseries( v_and_i_timetrace_plasma, v_and_i_timetrace_vacuum );

nc = 24;

close all;
rise_fall_select = 1;
switch rise_fall_select
    case 1
        V_axis = iv_timeseries.rise(nc).V_axis;
        I_of_V = iv_timeseries.rise(nc).I_of_V;
        time = iv_timeseries.total_time_axis_rise(nc);
    case 2
        V_axis = iv_timeseries.fall(nc).V_axis;
        I_of_V = iv_timeseries.fall(nc).I_of_V;
        time = iv_timeseries.total_time_axis_fall(nc);
    case 3
        V_axis = cat(1,iv_timeseries.rise(nc).V_axis,iv_timeseries.fall(nc).V_axis);
        [V_axis,v_order] = sort(V_axis);
        I_of_V = cat(1,iv_timeseries.rise(nc).I_of_V,iv_timeseries.fall(nc).I_of_V);
        I_of_V = I_of_V(v_order);
        time = (iv_timeseries.total_time_axis_rise(nc)+iv_timeseries.total_time_axis_fall(nc))./2;
end
%Convert to microseconds
time = time*1e6;
delta_I = iv_timeseries.I_ADC_step;

L = length(V_axis);
n_sgolay = 5;
frame_length = round(L/4);
if mod(frame_length,2) == 0
    frame_length = frame_length + 1;
end

%Apply Savitzky�Golay filter
I_of_V_sgolay = sgolayfilt(I_of_V,n_sgolay,frame_length);
fprintf('Savitzky�Golay filter applied...\n');

%Plasma detector: abort if no plasma is present
if (max(I_of_V_sgolay)/abs(min(I_of_V_sgolay)) < 3 && max(abs(I_of_V_sgolay)) < 2*delta_I)
    fprintf('No plasma detected: abort parameter extraction\n');
    return
end

% %Focus on a meaningful range of indices (exclude margins)
% margin_percent = 0.1;
% margin = round(margin_percent*length(V_axis));
% margin_range = margin:length(V_axis)-margin;

%Polynomial approximation of I =I(V) and its derivatives
n_poly = 5

%Polynomial approximations of I(V), I'(V), I''(V)
p_I = polyfit(V_axis,I_of_V_sgolay,n_poly);
I = polyval(p_I,V_axis);

% OLD METHOD: Polynomial fit of SG-filtered data as basis for derivatives
p_DI = polyder(p_I);
DI = polyval(p_DI,V_axis);
%----------------
DI_num = diff(I);
%----------------

p_D2I = polyder(p_DI);
D2I = polyval(p_D2I,V_axis);

%Shift the I = I(V) curve upward to avoid negative values in the logarithm
if min(I_of_V_sgolay) < 0
    shift = -(1.1)*min(I_of_V_sgolay);
else
    shift = 0;
end

%Polynomial approximation of log(I+I_shift)
n_poly_log = 5

p_logI = polyfit(V_axis,log(I_of_V_sgolay + shift),n_poly_log);
logI = polyval(p_logI,V_axis);

p_DlogI = polyder(p_logI);
DlogI = polyval(p_DlogI,V_axis);

p_D2logI = polyder(p_DlogI);
D2logI = polyval(p_D2logI,V_axis);

% %Extract floating potential (zero-crossing V of I=I(V)) (first approximation)
% [ ~, V_float_ind ] = min(abs(I_of_V_sgolay));
% if ismember(V_float_ind,1:length(V_axis))
%     V_float = V_axis(V_float_ind)
% else
%     fprintf('Floating potential is out of range\n');
%     V_float = nan;
% end

%Extract plasma potential (a root of D2I)
r = roots(p_D2I);
r = r(r > min(V_axis) & r < max(V_axis) & r >= 0);

%NEEDS TO BE ADDED IN FUNCTION---------------------------------
if isempty(r)
    fprintf('Plasma potential is out of range\n');
    return
end
%--------------------------------------------------------------

V_plasma = min(r);
[~,V_plasma_ind] = min(abs(V_axis-V_plasma));

% %Adjust plasma potential value
% p_D3I = polyder(p_D2I);
% D3I = polyval(p_D3I,V_axis);
% %Extract plasma potential (a root of D3I)
% r = roots(p_D3I);
% 
% if isempty(r)
%     fprintf('Plasma potential adjustment not possible.\n');
% else
%     V_plasma = (V_plasma + max(r))/2;
%     [~,V_plasma_ind] = min(abs(V_axis-V_plasma));
% end


%Look for peaks of DLogI
[pks,locs] = findpeaks(DlogI);
% locs = locs(locs > margin_range(1) & locs < margin_range(end));
% pks = DlogI(locs);
if ~isempty(pks)
    mp = pks(1);
else
    fprintf('Derivative of log(I+I_shift) has no peaks: parameter extraction aborted\n');
    return
end

fitresult_list = cell(length(locs),1);

%NEW FIT RANGE: from beginning of V_axis to plasma potential
fit_range_right_end = V_plasma_ind;
% fit_range = intersect(margin_range,1:fit_range_right_end);
fit_range = 1:fit_range_right_end;
if isempty(fit_range)
    fprintf('Fit range does not contain any indices \n');
    return
end

V_tf = V_axis(fit_range);
I_tf = I_of_V_sgolay(fit_range);


[xData, yData] = prepareCurveData( V_tf, I_tf );

ft = fittype( 'I_0 + a*V + b*exp(c*V)', 'independent', 'V', 'dependent', 'I' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf 0 0.001 0.01];
opts.MaxFunEvals = 600;
opts.MaxIter = 600;
opts.Robust = 'Bisquare';
opts.StartPoint = [0 0.0001 0.1 mp];
opts.TolFun = 0.0001;
opts.TolX = 0.0001;
opts.Upper = [0 Inf 1 Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts )
cv = coeffvalues(fitresult);

I_0 = cv(1);
a = cv(2);
b = cv(3);
c = cv(4);

I_of_V_fitted = I_0 + a.*V_axis + b.*exp(c.*V_axis);

%Calculate floating potential using fitted model
[ ~, V_float_ind ] = min(abs(I_of_V_fitted));
if ismember(V_float_ind,1:length(V_axis))
    V_float = V_axis(V_float_ind)
else
    fprintf('Floating potential is out of range\n');
    V_float = nan;
end




I_i = I_0 + a.*V_axis;


T_fit = 1/c


%Adjust plot-range
[~,V_plasma_ind] = min(abs(V_axis-V_plasma));
plot_range = (1:V_plasma_ind);

%Adjust plasma potential value
p_D3I = polyder(p_D2I);
D3I = polyval(p_D3I,V_axis);
%Extract plasma potential (a root of D3I)
r = roots(p_D3I);

if isempty(r)
    fprintf('Plasma potential adjustment not possible.\n');
else
    V_plasma = (V_plasma + max(r))/2;
    [~,V_plasma_ind] = min(abs(V_axis-V_plasma));
end



figure
set(gcf,'Renderer','painters');
ax_IV = axes;
vec = ones(length(V_axis),1);
errorbar(ax_IV,V_axis,I_of_V,delta_I.*vec,'ko');
hold(ax_IV,'on');
plot(ax_IV,V_axis,I_of_V_sgolay,'b-','LineWidth',2);
switch rise_fall_select
    case 1
        title(ax_IV,['IV-characteristic at time = ',num2str(time),' �s (Rising edge)']);
    case 2
        title(ax_IV,['IV-characteristic at time = ',num2str(time),' �s (Falling edge)']);
end
xlabel(ax_IV,'Voltage [V]');
ylabel(ax_IV,'Current [A]');
xlim_temp = xlim;


plot(ax_IV,V_axis(plot_range),I_of_V_fitted(plot_range),'r-','LineWidth',2);
plot(ax_IV,V_axis,I_i,'r--','LineWidth',2);
hold(ax_IV,'off');

%Draw a vertical line at floating potential position
axes(ax_IV);
x_bound = [V_float V_float];
y_bound = get(ax_IV,'YLim');
line(x_bound,y_bound,'Color',[0 1 0],'LineWidth',2);

%Draw a vertical line at plasma potential position
axes(ax_IV);
x_bound = [V_plasma V_plasma];
y_bound = get(ax_IV,'YLim');
line(x_bound,y_bound,'Color',[1 0 1],'LineWidth',2);

% h_leg = legend(ax_IV,'Data points','Savitzky-Golay filtered data','Fit with model I(V) = I_0 + aV + be^{cV}','Ion current I_i(V) = I_0 + aV','Floating potential','Plasma potential');
% set(h_leg,'FontSize',10);

V_ion_saturation = V_float - 100;
I_sat_i_fit = I_0 + a*V_ion_saturation



figure
set(gcf,'Renderer','painters');
movegui('northeast');

ax_logI = subplot(3,1,1);
plot(ax_logI,V_axis,logI);
title('log(I+I_{shift})');
xlim(xlim_temp);

ax_DlogI = subplot(3,1,2);
plot(ax_DlogI,V_axis,DlogI);
title('Dlog(I+I_{shift})');

ax_D2logI = subplot(3,1,3);
plot(ax_D2logI,V_axis,D2logI);
title('D2log(I+I_{shift})');


figure
set(gcf,'Renderer','painters');
movegui('southeast');

ax_I = subplot(3,1,1);
plot(ax_I,V_axis,I);
title('I = I(V)');
xlim(xlim_temp);

ax_DI = subplot(3,1,2);
plot(ax_DI,V_axis,DI);
title('I'' = I''(V)');

ax_D2I = subplot(3,1,3);
plot(ax_D2I,V_axis,D2I);
title('I'''' = I''''(V)');
x_bound = get(ax_IV,'XLim'); 
y_bound = [0 0];
line(x_bound,y_bound,'Color',[1 0 0]);
