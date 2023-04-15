
% %FLoating potential 2Dspace*1Dtime*1Drepetition array
% V_float_XYTR_rise = NaN(nX,nY,100,nr);
% V_float_XYTR_fall = NaN(nX,nY,100,nr);
% %Plasma potential 2Dspace*1Dtime*1Drepetition array
% V_plasma_XYTR_rise = NaN(nX,nY,100,nr);
% V_plasma_XYTR_fall = NaN(nX,nY,100,nr);
% %Electron temperature 2Dspace*1Dtime*1Drepetition array
% T_e_XYTR_rise = NaN(nX,nY,100,nr);
% T_e_XYTR_fall = NaN(nX,nY,100,nr);
% %Ion saturation current 2Dspace*1Dtime*1Drepetition array
% I_sat_i_XYTR_rise = NaN(nX,nY,100,nr);
% I_sat_i_XYTR_fall = NaN(nX,nY,100,nr);


%             V_float_XYTR_rise(k_x,k_y,1:iv_timeseries.NC_tot,k_r) = plasma_param_timeseries.rise.V_float;
%             V_float_XYTR_fall(k_x,k_y,1:iv_timeseries.NC_tot,k_r) = plasma_param_timeseries.fall.V_float;
%             
%             V_plasma_XYTR_rise(k_x,k_y,1:iv_timeseries.NC_tot,k_r) = plasma_param_timeseries.rise.V_plasma;
%             V_plasma_XYTR_fall(k_x,k_y,1:iv_timeseries.NC_tot,k_r) = plasma_param_timeseries.fall.V_plasma; 
%             
%             T_e_XYTR_rise(k_x,k_y,1:iv_timeseries.NC_tot,k_r) = plasma_param_timeseries.rise.T_e;
%             T_e_XYTR_fall(k_x,k_y,1:iv_timeseries.NC_tot,k_r) = plasma_param_timeseries.fall.T_e;
%             
%             I_sat_i_XYTR_rise(k_x,k_y,1:iv_timeseries.NC_tot,k_r) = plasma_param_timeseries.rise.I_sat_i;
%             I_sat_i_XYTR_fall(k_x,k_y,1:iv_timeseries.NC_tot,k_r) = plasma_param_timeseries.fall.I_sat_i;



% %Save all the extracted plasma information
% save([target_folder,'plasma_param_timeseries_cell_array.mat'],'plasma_param_timeseries_cell_array');

% V_float_XYT_rise = nanmean(V_float_XYTR_rise,4);
% save([save_folder,'V_float_XYT_rise.mat'],'V_float_XYT_rise');
% V_float_XYT_fall = nanmean(V_float_XYTR_fall,4);
% save([save_folder,'V_float_XYT_fall.mat'],'V_float_XYT_fall');
% 
% V_plasma_XYT_rise = nanmean(V_plasma_XYTR_rise,4);
% save([save_folder,'V_plasma_XYT_rise.mat'],'V_plasma_XYT_rise');
% V_plasma_XYT_fall = nanmean(V_plasma_XYTR_fall,4);
% save([save_folder,'V_plasma_XYT_fall.mat'],'V_plasma_XYT_fall');
% 
% T_e_XYT_rise = nanmean(T_e_XYTR_rise,4);
% save([save_folder,'T_e_XYT_rise.mat'],'T_e_XYT_rise');
% T_e_XYT_fall = nanmean(T_e_XYTR_fall,4);
% save([save_folder,'T_e_XYT_fall.mat'],'T_e_XYT_fall');
% 
% I_sat_i_XYT_rise = nanmean(I_sat_i_XYTR_rise,4);
% save([save_folder,'I_sat_i_XYT_rise.mat'],'I_sat_i_XYT_rise');
% I_sat_i_XYT_fall = nanmean(I_sat_i_XYTR_fall,4);
% save([save_folder,'I_sat_i_XYT_fall.mat'],'I_sat_i_XYT_fall');