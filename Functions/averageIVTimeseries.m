function [ iv_time_average ] = averageIVTimeseries( iv_timeseries )
%AVERAGEIVTIMESERIES Takes the time_average of an iv_timeseries under the assumption of constant plasma parameters
%   Detailed explanation goes here

    %Perform the average over all the cycles

    iv_time_average.rise.V = cat(1,iv_timeseries.rise(:).V_axis);
    iv_time_average.rise.I = cat(1,iv_timeseries.rise(:).I_of_V);

    iv_time_average.fall.V = cat(1,iv_timeseries.fall(:).V_axis);
    iv_time_average.fall.I = cat(1,iv_timeseries.fall(:).I_of_V);

    %Reconstruct time-averaged IV-curve

    rec_curve_out = reconstructCurve(iv_time_average.rise.V,iv_time_average.rise.I);
    iv_time_average.rise.V_axis = rec_curve_out.x;
    iv_time_average.rise.I_of_V = rec_curve_out.y_of_x;

    rec_curve_out = reconstructCurve(iv_time_average.fall.V,iv_time_average.fall.I);
    iv_time_average.fall.V_axis = rec_curve_out.x;
    iv_time_average.fall.I_of_V = rec_curve_out.y_of_x;

end

