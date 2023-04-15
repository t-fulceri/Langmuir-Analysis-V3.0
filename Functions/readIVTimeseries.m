function [ return_status, iv_timeseries ] = readIVTimeseries( source_folder, filename, settings)%vac_lc_effect, voltage_deamplification )
%readIVTimeseries reads the voltage/current data from one measurement file and returns a timeseries of I = I(V) curves, each corresponding to a voltage cycle

    fprintf('Runnning readIVTimeseries\n');

    [ v_and_i_timetrace_plasma ] = readHD5Measurement( source_folder, filename, settings);
    [ return_status, iv_timeseries ] = extractIVTimeseries(v_and_i_timetrace_plasma);

    return_status = 0;

end

