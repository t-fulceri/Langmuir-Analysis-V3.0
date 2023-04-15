function [ fitted_vac_lc_effect ] = extractLC_VacuumEffect( source_folder, filename_cell_array, voltage_deamplification )
%extractLC_VacuumEffect gathers the LC signal of the Langmuir measuring
%circuit from a series of measurements

    fprintf('Runnning extractLC_VacuumEffect\n');
    
    if nargin == 2
        voltage_deamplification = 1;
    end

    voltage_sign = +1;
    current_sign = -1; %(IV-Curve convention: invert current axis);

    vac_lc_meas_list.N_files = length(filename_cell_array);
    fprintf(['Number of files: ',num2str(vac_lc_meas_list.N_files),'\n']);

    fprintf(['Folder path: ',source_folder,'\n']);
    fprintf('File names:\n');
    for nf = 1:vac_lc_meas_list.N_files
        fprintf([filename_cell_array{nf},'\n']);
    end

    full_file_path = strcat(source_folder,filename_cell_array{nf});
    %Check metadata
    metadata = h5info(full_file_path);
    %Extract timestep from metadata [s]
    time_step = metadata.Groups.Attributes(2,1).Value(1);
%     time_step_array(nf) = time_step;
    %Extract date and time from metadata
%     iv_timeseries.date_array(nf) = metadata.Datasets.Attributes(1,1).Value;
%     iv_timeseries.time_array(nf) = metadata.Datasets.Attributes(2,1).Value;
    
    %Read data
    data = h5read(full_file_path,'/dataset');
    %Transpose
    data = data';
    %Extract voltage
    voltage = data(:,1);
    %Optional inversion of voltage signal (depends on the data source)
    voltage = voltage_deamplification.*voltage_sign.*voltage;
    %Extract current
    current = data(:,2);
    %Optional inversion of current signal (depends on the data source)
    current = current_sign.*current;

    %Calculate sampling frequency [Hz]
    sampling_freq = 1/time_step;

    reference = current;
    signal = current;

    %Keep note of important values
    L = length(signal);
    time_step = 1/sampling_freq;
    t = time_step*linspace(0,L-1,L);

    %Some starting parameters
    [ main_freq, main_phase ] = extractMainFreqAndPhase( reference, sampling_freq );

    pulsation_start = main_freq*2*pi*time_step;

    [xData, yData] = prepareCurveData( [], reference );

    % Set up fittype and options.
    ft = fittype( 'sin1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [max(reference) pulsation_start 2.5];

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    cv = coeffvalues(fitresult);
    A = cv(1);
    omega = cv(2)/time_step;
    phi = cv(3);

    fitted_vac_lc_effect = A*sin(omega*t + phi);
    fitted_vac_lc_effect = fitted_vac_lc_effect';


end

