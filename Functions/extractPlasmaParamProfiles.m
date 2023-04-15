function [ T_e_timeseries_XYR, I_sat_i_timeseries_XYR, V_plasma_timeseries_XYR ] = extractPlasmaParamProfiles( plasma_param_timeseries_cell_array )
%EXTRACTPLASMAPARAMPROFILES Extracts (2+1)D profiles of (T_e, I_sat_i, V_p) from a plasma_param_timeseries_cell_array structure

    fprintf('Executing extractPlasmaParamProfiles.\n');
    
    s = size(plasma_param_timeseries_cell_array);
    
    switch length(s)
        case 2
            NX = s(1);
            NY = s(2);
            NR = 1;
        case 3
            NX = s(1);
            NY = s(2);
            NR = s(3);
        otherwise
            fprintf('Plasma Cell Array has the wrong size\n');
            fprintf('Plasma Cell Array must be a NXxNYxNR or NXxNY cell array!\n\n');
            return
    end
    
    %Prepare the output structures
    if NR > 1
        T_e_timeseries_XYR = cell(NX,NY,NR);
        I_sat_i_timeseries_XYR = cell(NX,NY,NR);
        V_plasma_timeseries_XYR = cell(NX,NY,NR);
    else
        T_e_timeseries_XYR = cell(NX,NY);
        I_sat_i_timeseries_XYR = cell(NX,NY);
        V_plasma_timeseries_XYR = cell(NX,NY);
    end
    
    for k_x = 1:NX
        for k_y = 1:NY
            for k_r = 1:NR
                fprintf(['X_index = ',num2str(k_x),'\n']);
                fprintf(['Y_index = ',num2str(k_y),'\n']);
                fprintf(['Repetition_index = ',num2str(k_r),'\n\n']);
                try
                    if NR > 1
                        plasma_param_timeseries = plasma_param_timeseries_cell_array{k_x,k_y,k_r};
                        T_e_timeseries_XYR{k_x,k_y,k_r} = plasma_param_timeseries.mean.T_e;
                        I_sat_i_timeseries_XYR{k_x,k_y,k_r} = plasma_param_timeseries.mean.I_sat_i;
                        V_plasma_timeseries_XYR{k_x,k_y,k_r} = plasma_param_timeseries.mean.V_plasma;
                    else
                        plasma_param_timeseries = plasma_param_timeseries_cell_array{k_x,k_y};
                        T_e_timeseries_XYR{k_x,k_y} = plasma_param_timeseries.mean.T_e;
                        I_sat_i_timeseries_XYR{k_x,k_y} = plasma_param_timeseries.mean.I_sat_i;
                        V_plasma_timeseries_XYR{k_x,k_y} = plasma_param_timeseries.mean.V_plasma;
                    end
                    

                    
                catch err_extract
                        T_e_timeseries_XYR{k_x,k_y} = [];
                        I_sat_i_timeseries_XYR{k_x,k_y} = [];
                        V_plasma_timeseries_XYR{k_x,k_y} = [];
                    fprintf('Loading of plasma parameter timeseries failed...\n\n');
                    rethrow(err_extract);
                    
                end
            end
        end
    end
    
    
    fprintf('extractPlasmaParamProfiles executed successfully.\n\n\n');

end