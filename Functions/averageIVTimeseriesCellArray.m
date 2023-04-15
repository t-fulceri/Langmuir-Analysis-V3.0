function [ iv_time_average_cell_array ] = averageIVTimeseriesCellArray( iv_timeseries_cell_array )
%AVERAGEIVTIMESERIESCELLARRAY Takes the time_average of a cell array of iv_timeseries and produces a cell array of time averages as outputs
%   Detailed explanation goes here

    %Take the size of the input cell array, normally it must be a 2D cell
    %array (1 repetition measurements) or a 3D cell array (multiple
    %repetition maeasurements)
    s = size(iv_timeseries_cell_array);    
    NX = s(1);
    NY = s(2);
    if length(s) == 2
        rep_per_point = 1;
    else
        if length(s) == 3;
            rep_per_point = s(3);
        end
    end
    
    if rep_per_point > 1
        %Fill in the new cell array
        iv_time_average_cell_array = cell(NX,NY,rep_per_point);
        for k_x = 1:NX
            for k_y = 1:NY
                for k_r = 1:rep_per_point
                    [ iv_time_average ] = averageIVTimeseries( iv_timeseries_cell_array{k_x,k_y,k_r} );
                    iv_time_average_cell_array{k_x,k_y,k_r} = iv_time_average;
                end
            end
        end
    else
        %Fill in the new cell array
        iv_time_average_cell_array = cell(NX,NY);
        for k_x = 1:NX
            for k_y = 1:NY

                [ iv_time_average ] = averageIVTimeseries( iv_timeseries_cell_array{k_x,k_y} );
                iv_time_average_cell_array{k_x,k_y} = iv_time_average;

            end
        end
    end

end

