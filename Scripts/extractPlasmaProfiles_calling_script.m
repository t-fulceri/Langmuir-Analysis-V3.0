% %Source folder for the plasma parameters timeseries list
% source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\Reconstructed_Profiles\2018_08_24_Plasma_Profile_100kHz_10V_1_rep_CompleteProfile_ECRH_On_E-Guns_On_X-Drive_On_1500V_Rec-Drive_On_1000V\plasma_parameters_timeseries_list\';
% 
% [ plasma_param_timeseries_cell_array ] = buildPlasmaParamTimerseriesCellArray( source_folder, 23, 18, 1);

[ T_e_timeseries_XYR, I_sat_i_timeseries_XYR, V_plasma_timeseries_XYR ] = extractPlasmaParamProfiles( plasma_param_timeseries_cell_array );

NX = 23;
NY = 18;
NT = 49;

T_e_XYT = nan(NX,NY,NT);

for k_x = 1:NX
    for k_y = 1:NY
        for k_t = 1:NT
            T_e_timeseries = T_e_timeseries_XYR{k_x,k_y};
            T_e_XYT(k_x,k_y,k_t) = T_e_timeseries(k_t);
            if T_e_XYT(k_x,k_y,k_t) > 40
                T_e_XYT(k_x,k_y,k_t) = nan;
            end
        end
    end
end

close all;
%contourf(0:220,0:170,(T_e_XYT(:,:,20))');
contourf((T_e_XYT(:,:,20))');
colormap(gca,hot);
colorbar(gca);