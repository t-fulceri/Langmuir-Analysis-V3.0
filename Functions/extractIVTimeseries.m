function [ return_status, iv_timeseries ] = extractIVTimeseries( v_and_i_timetrace_plasma, v_and_i_timetrace_spurious )
%EXTRACTIVTIMESERIES Extracts a time-series of IV-characteristics from V(t) and I(t), optionally correcting for spurious effects

    fprintf('Runnning extractIVTimeseries\n');
    
    %Optional correction for vacuum effect
    if nargin == 2
        v_and_i_timetrace_plasma.I_timetrace = v_and_i_timetrace_plasma.I_timetrace - v_and_i_timetrace_spurious.I_timetrace;
        fprintf('Timetrace of current I(t) corrected for spurious effects\n');
    end
    
    iv_timeseries = struct;
    
    iv_timeseries.date = v_and_i_timetrace_plasma.date;
    iv_timeseries.time = v_and_i_timetrace_plasma.time;
    iv_timeseries.V_timetrace = v_and_i_timetrace_plasma.V_timetrace;
    iv_timeseries.I_timetrace = v_and_i_timetrace_plasma.I_timetrace;
    iv_timeseries.V_ADC_step = v_and_i_timetrace_plasma.V_ADC_step;
    iv_timeseries.I_ADC_step = v_and_i_timetrace_plasma.I_ADC_step;
    iv_timeseries.time_step = v_and_i_timetrace_plasma.time_step;
    iv_timeseries.sampling_freq = v_and_i_timetrace_plasma.sampling_freq;

    %Split the timeseries signal into chunks of equal duration
    %(each corresponding to one voltage cycle)
    voltage = iv_timeseries.V_timetrace;
    current = iv_timeseries.I_timetrace;
    sampling_freq = iv_timeseries.sampling_freq;
    time_step = iv_timeseries.time_step;
    [V_cycles,iv_timeseries.initial_time_steps_offset] =  timeseriesSplitSinefit(voltage,voltage,sampling_freq);
    I_cycles =  timeseriesSplitSinefit(voltage,current,sampling_freq);

    %Total number of cycles
    iv_timeseries.NC_tot = size(V_cycles,2);
    %Time steps per cycle
    iv_timeseries.time_steps_per_cycle = size(V_cycles,1);



    %Produce the time axis
    %time_axis = time_step*ts_per_cycle*conversion_factor*(linspace(1,NC_tot,NC_tot) - 1/2);
    %time_axis = time_axis + DAQ_delay*conversion_factor;


    %total_time_axis = (1/sampling_freq)*conversion_factor*linspace(0,length(voltage)-1,length(voltage));
    %total_time_axis = total_time_axis + DAQ_delay*conversion_factor;

    iv_timeseries.cycle_length = size(V_cycles,1);
    hl = round(iv_timeseries.cycle_length/2);

    iv_timeseries.V_cycles_rise = V_cycles(1:hl,:);
    iv_timeseries.I_cycles_rise = I_cycles(1:hl,:);

    iv_timeseries.V_cycles_fall = V_cycles(hl+1:end,:);
    iv_timeseries.I_cycles_fall = I_cycles(hl+1:end,:);
    

    iv_timeseries.time_step_in_seconds = time_step;
    iv_timeseries.sampling_freq = 1/time_step;
    %iv_timeseries.time_steps_per_cycle = time_steps_per_cycle_array(1);
    iv_timeseries.time_uncertainty = round(iv_timeseries.time_steps_per_cycle/4)*time_step;
        
    %iv_timeseries.initial_time_steps_offset = iv_timeseries.initial_time_steps_offset_array(1);


    %Total number of time steps
    iv_timeseries.total_time_steps = iv_timeseries.initial_time_steps_offset + iv_timeseries.time_steps_per_cycle*iv_timeseries.NC_tot;
    fprintf('Total number of time steps extracted\n');
    %Total time axis for rising edge
    iv_timeseries.total_time_axis_rise = iv_timeseries.time_step_in_seconds*(iv_timeseries.initial_time_steps_offset + round(iv_timeseries.time_steps_per_cycle/4) + (0:iv_timeseries.NC_tot-1)*iv_timeseries.time_steps_per_cycle);
    fprintf('Total time axis for rising edge extracted\n');
    %Total time axis for falling edge
    iv_timeseries.total_time_axis_fall = iv_timeseries.time_step_in_seconds*(iv_timeseries.initial_time_steps_offset + round(iv_timeseries.time_steps_per_cycle*(3/4)) + (0:iv_timeseries.NC_tot-1)*iv_timeseries.time_steps_per_cycle);
    fprintf('Total time axis for falling edge extracted\n');
    
    %Stack together the datapoints and data uncertainties which correspond to
    %the same cycle
    iv_cycles_stack = struct('V_cycles_rise',[],'I_cycles_rise',[],'V_cycles_fall',[],'I_cycles_fall',[]);

    iv_cycles_stack.V_cycles_rise = cat(1,iv_cycles_stack.V_cycles_rise,iv_timeseries.V_cycles_rise);
    iv_cycles_stack.I_cycles_rise = cat(1,iv_cycles_stack.I_cycles_rise,iv_timeseries.I_cycles_rise);
    iv_cycles_stack.V_cycles_fall = cat(1,iv_cycles_stack.V_cycles_fall,iv_timeseries.V_cycles_fall);
    iv_cycles_stack.I_cycles_fall = cat(1,iv_cycles_stack.I_cycles_fall,iv_timeseries.I_cycles_fall);

    fprintf('Datapoints have been stacked together\n');

    %Sort the V(t) and I(t) signals in every cycle following the V-order
    for nc = 1:iv_timeseries.NC_tot
        [iv_cycles_stack.V_cycles_rise(:,nc),sort_ind] = sort(iv_cycles_stack.V_cycles_rise(:,nc));
        iv_cycles_stack.I_cycles_rise(:,nc) = iv_cycles_stack.I_cycles_rise(sort_ind,nc);

        [iv_cycles_stack.V_cycles_fall(:,nc),sort_ind] = sort(iv_cycles_stack.V_cycles_fall(:,nc));
        iv_cycles_stack.I_cycles_fall(:,nc) = iv_cycles_stack.I_cycles_fall(sort_ind,nc);
    end
    fprintf('V(t) and I(t) signals sorted according to V-order\n');


    for nc = 1:iv_timeseries.NC_tot

        rec_curve_out = reconstructCurve(iv_cycles_stack.V_cycles_rise(:,nc),iv_cycles_stack.I_cycles_rise(:,nc));
        iv_timeseries.rise(nc).V_axis = rec_curve_out.x;
        iv_timeseries.rise(nc).I_of_V = rec_curve_out.y_of_x;

        rec_curve_out = reconstructCurve(iv_cycles_stack.V_cycles_fall(:,nc),iv_cycles_stack.I_cycles_fall(:,nc));
        iv_timeseries.fall(nc).V_axis = rec_curve_out.x;
        iv_timeseries.fall(nc).I_of_V = rec_curve_out.y_of_x;

    end
    fprintf('IV-timeseries ready\n');

    return_status = 0;

end

