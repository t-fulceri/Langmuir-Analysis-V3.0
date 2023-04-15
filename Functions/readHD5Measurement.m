function [ v_and_i_timetrace ] = readHD5Measurement( source_folder, filename, settings)
%READHD5MEASUREMENT %readIVTimeseries reads the voltage/current data from a HD5 measurement file and returns the v_and_i_timetraces data strucure

    fprintf(['readHD5Measurement on file \n',filename]);
    fprintf('Settings:\n');
    voltage_sign = settings.voltage_sign%+1;
    current_sign = settings.current_sign%-1; %(IV-Curve convention: invert current axis);    
    voltage_deamplification = settings.voltage_deamplification

    v_and_i_timetrace = struct;

    fprintf(['Folder path: ',source_folder,'\n']);
    fprintf('File name:\n');
    fprintf(filename,'\n');

    %Collect all the relevant information to form the I_cycles and V_cycles
    %frome each measurement file

    full_file_path = strcat(source_folder,filename);
    %Check metadata
    metadata = h5info(full_file_path);
    %Extract timestep from metadata [s]
    v_and_i_timetrace.time_step = metadata.Groups.Attributes(2,1).Value(1);
    %Extract date and time from metadata
    v_and_i_timetrace.date = metadata.Datasets.Attributes(1,1).Value;
    v_and_i_timetrace.time = metadata.Datasets.Attributes(2,1).Value;

    %Read data
    data = h5read(full_file_path,'/dataset');
    %DAQ delay in s
    %DAQ_delay = 20e-3;
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
    
    %Take note of V(t)[V] and I(t)[A] timetraces
    v_and_i_timetrace.V_timetrace = voltage;
    v_and_i_timetrace.I_timetrace = current;

    %Extract uncertainty coming from ADC discretization;
    v_and_i_timetrace.V_ADC_step = signalToADCStep(voltage);
    v_and_i_timetrace.I_ADC_step = signalToADCStep(current);

    %Calculate sampling frequency [Hz]
    v_and_i_timetrace.sampling_freq = 1/v_and_i_timetrace.time_step;
end

