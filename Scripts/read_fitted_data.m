close all
clear all

addpath('..\Functions\');
source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\Reconstructed_Profiles\2018_08_22_Plasma_Profile_100kHz_10V_5_rep_CompleteProfile_X-Drive_On_1500V\'
load([source_folder,'plasma_param_timeseries_cell_array.mat']);
load([source_folder,'iv_timeseries_cell_array.mat']);
iv_timeseries = iv_timeseries_cell_array{1,1,1};

nX = 32;
nY = 18;
nr = 5;

%FLoating potential 2Dspace*1Dtime*1Drepetition array
V_float_XYTR_rise = NaN(nX,nY,100,nr);
V_float_XYTR_fall = NaN(nX,nY,100,nr);
%Plasma potential 2Dspace*1Dtime*1Drepetition array
V_plasma_XYTR_rise = NaN(nX,nY,100,nr);
V_plasma_XYTR_fall = NaN(nX,nY,100,nr);
%Electron temperature 2Dspace*1Dtime*1Drepetition array
T_e_XYTR_rise = NaN(nX,nY,100,nr);
T_e_XYTR_fall = NaN(nX,nY,100,nr);
%Ion saturation current 2Dspace*1Dtime*1Drepetition array
I_sat_i_XYTR_rise = NaN(nX,nY,100,nr);
I_sat_i_XYTR_fall = NaN(nX,nY,100,nr);

for k_x = 1:nX
    for k_y = 1:nY
        for k_r = 1:nr
            k_x_disp = k_x
            k_y_disp = k_y            
            k_r_disp = k_r
            
            plasma_param_timeseries = plasma_param_timeseries_cell_array{k_x,k_y,k_r};
            
            V_float_XYTR_rise(k_x,k_y,1:iv_timeseries.NC_tot,k_r) = plasma_param_timeseries.rise.V_float;
            V_float_XYTR_fall(k_x,k_y,1:iv_timeseries.NC_tot,k_r) = plasma_param_timeseries.fall.V_float;
            
            V_plasma_XYTR_rise(k_x,k_y,1:iv_timeseries.NC_tot,k_r) = plasma_param_timeseries.rise.V_plasma;
            V_plasma_XYTR_fall(k_x,k_y,1:iv_timeseries.NC_tot,k_r) = plasma_param_timeseries.fall.V_plasma; 
            
            T_e_XYTR_rise(k_x,k_y,1:iv_timeseries.NC_tot,k_r) = plasma_param_timeseries.rise.T_e;
            T_e_XYTR_fall(k_x,k_y,1:iv_timeseries.NC_tot,k_r) = plasma_param_timeseries.fall.T_e;
            
            I_sat_i_XYTR_rise(k_x,k_y,1:iv_timeseries.NC_tot,k_r) = plasma_param_timeseries.rise.I_sat_i;
            I_sat_i_XYTR_fall(k_x,k_y,1:iv_timeseries.NC_tot,k_r) = plasma_param_timeseries.fall.I_sat_i;
        end
    end
end