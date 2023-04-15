%Source folder for the plasma parameters timeseries list
source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\Reconstructed_Profiles\2018_08_21_Plasma_Profile_01kHz_10V_1_rep_CurrentMonitor_CompleteProfile_ECRH_Only_DeamplifiedBy200\';
load([source_folder,'iv_timeseries_cell_array.mat']);
load([source_folder,'plasma_param_timeseries_cell_array.mat']);

%Take note of the size of the multidimensional array
s = size(iv_timeseries_cell_array);
NX = s(1);
NY = s(2);
if length(s) == 3
    NR = s(3);
else
    NR = 1;
end

%Need the total number of cycles, pick it from the first iv_timeseries
NC = iv_timeseries_cell_array{1,1,1}.NC_tot;

%Prepare 4 arrays of dimension [NX,NY,NR,NC,2] (one for each plasma parameter)
T_e_XYRCH      = nan(NX,NY,NR,NC,2);
I_sat_i_XYRCH  = nan(NX,NY,NR,NC,2);
V_plasma_XYRCH = nan(NX,NY,NR,NC,2);
V_float_XYRCH  = nan(NX,NY,NR,NC,2);

%Establish boundaries
T_e_min = 0;
T_e_max = 20;

I_sat_i_min = -0.2;
I_sat_i_max = 0;



%Fill the arrays with the values
for k_x = 1:NX
    for k_y = 1:NY
        for k_r = 1:NR
            fprintf(['X_index = ',num2str(k_x),'\n']);
            fprintf(['Y_index = ',num2str(k_y),'\n']);
            fprintf(['Repetition_index = ',num2str(k_r),'\n\n']);  
            try
                fprintf('Retrieving plasma parameter values...\n');
                if length(s) == 3
                    pp_ts = plasma_param_timeseries_cell_array{k_x,k_y,k_r};
                else
                    pp_ts = plasma_param_timeseries_cell_array{k_x,k_y};
                end
                for k_c = 1:NC

                    for k_rf = 1:2
                        if k_rf == 1
                            T_e      = pp_ts.rise.T_e(k_c);
                            I_sat_i  = pp_ts.rise.I_sat_i(k_c);
                            V_plasma = pp_ts.rise.V_plasma(k_c);
                            V_float  = pp_ts.rise.V_float(k_c);
                        else
                            T_e      = pp_ts.fall.T_e(k_c);
                            I_sat_i  = pp_ts.fall.I_sat_i(k_c);
                            V_plasma = pp_ts.fall.V_plasma(k_c);
                            V_float  = pp_ts.fall.V_float(k_c);
                        end
                    end
                    
                    T_ok = true;%(T_e >= T_e_min && T_e <= T_e_max);
                    I_sat_i_ok = true;%(I_sat_i >= I_sat_i_min && I_sat_i <= I_sat_i_max);
                    V_plasma_ok = (V_plasma >= 0 && isreal(V_plasma));
                    
                    if(T_ok && I_sat_i_ok && V_plasma_ok)
                        T_e_XYRCH(k_x,k_y,k_r,k_c,k_rf)      = T_e;
                        I_sat_i_XYRCH(k_x,k_y,k_r,k_c,k_rf)  = I_sat_i;
                        V_plasma_XYRCH(k_x,k_y,k_r,k_c,k_rf) = V_plasma;
                        V_float_XYRCH(k_x,k_y,k_r,k_c,k_rf)  = V_float;
                    end
                    
                end
                fprintf('Done.\n\n');
            catch
                fprintf('Failed.\n');
                fprintf('Jumping to next iteration...\n\n');
            end
            
        end
    end
end

%Average over repetitions
T_e_XYCH      = nanmean(T_e_XYRCH,3);
I_sat_i_XYCH  = nanmean(I_sat_i_XYRCH,3);
V_plasma_XYCH = nanmean(V_plasma_XYRCH,3);
V_float_XYCH  = nanmean(V_float_XYRCH,3);

%Average over rise/fall
T_e_XYC_temp      = nanmean(T_e_XYCH,5);
I_sat_i_XYC_temp  = nanmean(I_sat_i_XYCH,5);
V_plasma_XYC_temp = nanmean(V_plasma_XYCH,5);
V_float_XYC_temp  = nanmean(V_float_XYCH,5);

T_e_XYC      = nan(NX,NY,NC);
I_sat_i_XYC  = nan(NX,NY,NC);
V_plasma_XYC = nan(NX,NY,NC);
V_float_XYC  = nan(NX,NY,NC);

for k_x = 1:NX
    for k_y = 1:NY
        for k_c = 1:NC
            T_e_XYC(k_x,k_y,k_c) = T_e_XYC_temp(k_x,k_y,1,k_c);
            I_sat_i_XYC(k_x,k_y,k_c) = I_sat_i_XYC_temp(k_x,k_y,1,k_c);
            V_plasma_XYC(k_x,k_y,k_c) = V_plasma_XYC_temp(k_x,k_y,1,k_c);
            V_float_XYC(k_x,k_y,k_c) = V_float_XYC_temp(k_x,k_y,1,k_c);
        end
    end
end




close all
nc =25;

% %Smooth the data a little bit...
space_smoothing = [3,3,1];
time_smoothing = [1,1,3];
spacetime_smoothing = 3;
% T_e_XYC = smooth3(T_e_XYC,'box',space_smoothing);
% I_sat_i_XYC = smooth3(I_sat_i_XYC,'box',space_smoothing);
V_plasma_XYC = smooth3(V_plasma_XYC,'box',space_smoothing);

%try a plot
figure
movegui(gcf,'northwest');
set(gcf,'Renderer','painters');
[~,h] = contourf(T_e_XYC(:,:,nc)');
set(h,'LineColor','none');
colormap(gca,hot);
h = colorbar(gca);
ylabel(h, 'T_e [eV]');

%try a plot
figure
movegui(gcf,'center');
set(gcf,'Renderer','painters');
[~,h] = contourf(-I_sat_i_XYC(:,:,nc)');
set(h,'LineColor','none');
colormap(gca,parula);
h = colorbar(gca);
ylabel(h, 'I_sat_i [A]');

%try a plot
figure
movegui(gcf,'northeast');
set(gcf,'Renderer','painters');
[~,h] = contourf(V_plasma_XYC(:,:,nc)',15);
%set(h,'LineColor','none');
colormap(gca,bone);
h = colorbar(gca);
ylabel(h, 'V_{plasma} [V]');

