function [ iv_timeseries_cell_array ] = buildIVTimeseriesCellArray( source_folder, rep_per_point, raster_step, LC_vacuum_source_folder, voltage_deamplification)
%BUILDIVTIMESERIESCELLARRAY Takes raw data and organizes it into a (NX x Ny x rep_per_point) multidimensional cell array, with each cell containing one iv_timeseries structure 
%   Detailed explanation goes here

    if nargin == 4
        voltage_deamplification = 1;
    end

    %Take note of the LC_VacuumEffect
    filename_cell_array = cell(1,1);
    filename_cell_array{1} = 'meas_0000.h5';
    [ fitted_vac_lc_effect ] = extractLC_VacuumEffect( LC_vacuum_source_folder, filename_cell_array, voltage_deamplification);

    %Read the map file
    map_file = 'map.map'
    full_file_path = strcat(source_folder,map_file);
    B = tdfread(full_file_path);
    ind = B.x0x25_idx;
    x_pos = B.x;
    y_pos = B.y;
    ind = ind(:,1);
    x_pos = x_pos(:,1);
    y_pos = y_pos(:,1);

    %Number of points

    %CORRECTION IN CASE OF DEFECTIVE (NON-LUBRICATED) LINEAR MOTORS
    x_pos_mod = raster_step*floor(x_pos/raster_step);
    x_pos = x_pos_mod;
    %----------

    NX = (max(x_pos) - min(x_pos))/raster_step + 1;
    NY = (max(y_pos) - min(y_pos))/raster_step + 1;

    %Organize data in multidimensional array
    %Measurement number 2Dspace*1Drepetition array
    meas_number_STR = NaN(NX,NY,rep_per_point);

    for k_x = 1:NX
        for k_y = 1:NY
            for k_r = 1:rep_per_point
                %file_name = strcat('meas_',num2str(k+(repetition-1),'%04.f'),'.h5');
                %file_names_STR(k_x,k_y,repetition) = file_name;
                k = k_x + (k_y-1)*NX;
                meas_number_STR(k_x,k_y,k_r) = rep_per_point*(k-1)+(k_r-1);
            end
        end
    end

    %Multidimensional cell array to store all the iv_timeseries
    iv_timeseries_cell_array = cell(NX,NY,rep_per_point);

    for k_x = 1:NX
        for k_y = 1:NY
            for k_r = 1:rep_per_point
                file_number = meas_number_STR(k_x,k_y,k_r);
                filename = strcat('meas_',num2str(file_number,'%04.f'),'.h5')
                try
                    settings.voltage_sign = +1;
                    settings.current_sign = -1; %(IV-Curve convention: invert current axis);    
                    settings.voltage_deamplification = 1;
                    [ return_status_read, iv_timeseries ] = readIVTimeseries( source_folder, filename, settings);
                    %Save the first iv_timeseries which can be extracted
                    if return_status_read == 0
                        iv_timeseries_cell_array{k_x,k_y,k_r} = iv_timeseries;
                    end
                catch err_read
                    fprintf('IV Timeseries Read failed: jumping to next interation');
                    continue
                end
            end
        end
    end

end

