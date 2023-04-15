function [ plasma_param_timeseries_cell_array ] = buildPlasmaParamTimerseriesCellArray( source_folder, NX, NY, NR)
%BUILDPLASMAPARAMTIMERSERIESCELLARRAY builds a plasma parameter timeseries
%cell array starting from individual plasma parameter in the source folder

    fprintf('Executing buildPlasmaParamTimerseriesCellArray.\n');
    
    %Measurement number 2Dspace*1Drepetition array
    meas_number_STR = NaN(NX,NY,NR);

    for k_x = 1:NX
        for k_y = 1:NY
            for k_r = 1:NR
                %file_name = strcat('meas_',num2str(k+(repetition-1),'%04.f'),'.h5');
                %file_names_STR(k_x,k_y,repetition) = file_name;
                k = k_x + (k_y-1)*NX;
                meas_number_STR(k_x,k_y,k_r) = NR*(k-1)+(k_r-1);
            end
        end
    end
    
    %Prepare the output plasma_param_timeseries_cell_array structure
    plasma_param_timeseries_cell_array = cell(NX,NY,NR);
    
    for k_x = 1:NX
        for k_y = 1:NY
            for k_r = 1:NR
                fprintf(['X_index = ',num2str(k_x),'\n']);
                fprintf(['Y_index = ',num2str(k_y),'\n']);
                fprintf(['Repetition_index = ',num2str(k_r),'\n\n']);
                try
                    %Produce the right file number
                    file_number = meas_number_STR(k_x,k_y,k_r);
                    filename = strcat('plasma_param_timeseries_',num2str(file_number,'%04.f'),'.mat')
                    temp_struct = load([source_folder,filename]);
                    plasma_param_timeseries = temp_struct.plasma_param_timeseries;
                    plasma_param_timeseries_cell_array{k_x,k_y,k_r} = plasma_param_timeseries;
                    fprintf('Plasma Parameter Timeseries loaded into cell array.\n\n');
                    
                    %save([target_folder,'plasma_parameters_timeseries_list\',filename],'plasma_param_timeseries');
                    %fprintf('Plasma Parameter Timeseries saved.\n');
                catch err_extract
                    plasma_param_timeseries_cell_array{NX,NY,NR} = [];
                    fprintf('Loading of plasma parameter timeseries failed: Jumping to next iteration.\n\n');
                    continue
                    
                end
            end
        end
    end
    
    
    fprintf('buildPlasmaParamTimerseriesCellArray executed successfully.\n\n\n');

end

