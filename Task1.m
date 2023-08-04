%% Cicunomics Challenge - Dimitrios Vogiatzis

%% 1st Task
clear;
clc;
%load the .mat file after cleaning and seperating each value in a single cell and calculating timesteps and cumulative time
load CIR001Char.mat

CSVcalc.P=CSVdata.Voltage(:).*CSVdata.Current(:);                   
CSVcalc.E=CSVcalc.P(:).*CSVtime.TIMESTEP(:);            %[W x sec or J]
CSVcalc.SE_kJ=-sum(CSVcalc.E)/1000;                     %[kJ]
CSVcalc.SE_kWh=CSVcalc.SE_kJ/3600;                      %[kWh]

%Assuming an OCV of 3.7V 
CSVcalc.V_nominal= 3.7;                                 %[V]
CSVcalc.Q=CSVcalc.SE_kJ/CSVcalc.V_nominal*1000;         %[As]
CSVcalc.CapacityAh= CSVcalc.Q*3600;                     %[Ah]

%% 2nd Task

% Seperation of discharge voltages and pauses from raw data
 CSVcalc.Voc_step=CSVdata.Voltage(1);
 a=0;
 b=0;
 time_dis=0;
 time_term=0;
 
for i=1:length(CSVdata.Voltage)
    if (CSVdata.Index(i)=='Discharge')
        a=a+1;
        time_dis=time_dis+CSVtime.TIMESTEP(i);
        CSVcalc.VoltageDischarge(a,1)=CSVdata.Voltage(i);
        CSVcalc.CurrentDischarge(a,1)=CSVdata.Current(i);
        CSVcalc.SoCDischarge(a,1)=CSVdata.SoC(i);
        CSVcalc.VoltageDischarge_time(a,1)=time_dis;
    end
     if (CSVdata.Index(i)=='Pause') && CSVdata.Voltage(i)>=3.38
        b=b+1;
        time_term=time_term+CSVtime.TIMESTEP(i);
        CSVcalc.OCV(b,1)=CSVdata.Voltage(i);
        CSVcalc.OCV_time(b,1)=time_term;
     end
end
clear a b time_dis time_term

%%Use the current interuption technique to make ECM

% Find the indices of the interruption points
CIT.interruption_indices = find(diff(CSVcalc.CurrentDischarge) == 0);

% Initialize an array to store the voltage drops
CIT.voltage_drops = zeros(length(CIT.interruption_indices), 1);

% Calculate the voltage drops for each interruption
for i = 1:length(CIT.interruption_indices)
    % Find the indices before and after the interruption
    CIT.index_before = CIT.interruption_indices(i);
    CIT.index_after = CIT.index_before + 1;

    % Get the voltages before and after the interruption
    CIT.V_before = CSVcalc.VoltageDischarge(CIT.index_before);
    CIT.V_after = CSVcalc.VoltageDischarge(CIT.index_after);

    % Calculate the voltage drop
    CIT.voltage_drops(i) = CIT.V_before - CIT.V_after;
end

% Display the voltage drops
figure;
plot(CIT.voltage_drops);
ylabel('voltage drops');
xlabel('inderruption indices');

% Initialize an array to store the time constants
CIT.time_constants = zeros(length(CIT.interruption_indices), 1);

% Calculate t_inter for each interruption
CIT.t_inter = zeros(length(CIT.interruption_indices), 1);

% Calculate the time constants for each interruption
for i = 1:length(CIT.interruption_indices)
    
    % Find the index corresponding to the interruption
    CIT.index_inter = CIT.interruption_indices(i);

    % Get the time value at the interruption index
    CIT.t_inter(i) = CSVcalc.VoltageDischarge_time(CIT.index_inter);
    
    % Find the indices before and after the interruption
    CIT.index_before = CIT.interruption_indices(i);
    CIT.index_after = CIT.index_before + 1;

    % Get the voltages before and after the interruption
    CIT.V_initial = CSVcalc.VoltageDischarge(CIT.index_before);
    CIT.V_final = CSVcalc.VoltageDischarge(CIT.index_after);

    % Find the time at which the voltage starts to drop during the interruption
    % Adjusted to handle cases where the voltage does not drop during the interruption
    CIT.voltage_drop_indices = find(CSVcalc.VoltageDischarge(CIT.index_before:end) ~= CIT.V_initial);
    if isempty(CIT.voltage_drop_indices)
        continue;  % Skip this interruption if no voltage drop is detected
    end
    CIT.t_drop_index = CIT.voltage_drop_indices(1) - 1;
    CIT.t_drop = CSVcalc.VoltageDischarge_time(CIT.index_before + CIT.t_drop_index);

    % Calculate the time constant
    CIT.time_constants(i) = (CIT.t_drop-CIT.t_inter(i)) / log(CIT.V_initial / CIT.V_final);
end

% Exclude 'inf' values from time_constants
CIT.time_constants_valid = CIT.time_constants(isfinite(CIT.time_constants));

% Filter out negative time constant values
CIT.time_constants_valid_pos = CIT.time_constants_valid(CIT.time_constants_valid >= 0);

% Display the time constants
figure;
plot(CIT.time_constants_valid_pos);
ylabel('time constant');
xlabel('index interruptions');

% Calculate the average time constant
CIT.t_avg = mean(CIT.time_constants_valid_pos);

% Calculate the resistance R1
CIT.R1= CIT.t_avg/(CSVcalc.CapacityAh*3600);

% Calculate the polarization resistance R2
CIT.R2 = CIT.R1 ./ (CSVcalc.SoCDischarge/100 .* (1 - CSVcalc.SoCDischarge/100));

% Calculate the double-layer capacitance C_dl
CIT.C_dl = CIT.t_avg ./ CIT.R2;

% Display the polarization resistance and double-layer capacitance
figure;
plot(CIT.R2);
xlabel('index interruptions');
ylabel('Polarization Resistance');

figure;
plot(CIT.C_dl);
xlabel('index interruptions');
ylabel('Double-layer capacitance');


% Filter SoC and voltage data based on SoC values greater than or equal to 10%
filtered_indices = (CSVcalc.SoCDischarge) >= 10;
SoC_filtered = CSVcalc.SoCDischarge(filtered_indices);
V_filtered = CSVcalc.VoltageDischarge(filtered_indices);

% Sort the SoC and voltage data in ascending order of SoC
[SoC_sorted, sort_idx] = sort(SoC_filtered);
V_sorted = V_filtered(sort_idx);

% Plot voltage versus SoC
figure;
plot(SoC_sorted, V_sorted);
xlim([10 100]);
xlabel('State of Charge (SoC)');
ylabel('Voltage (V)');

% Average the OCV values
OCV_avg = mean(CSVcalc.OCV);

% Simulate the voltage response for different discharge pulses (According to ECM)
ECM.simulated_voltage = zeros(size(CSVcalc.VoltageDischarge_time));  % Initialize an array to store the simulated voltage
ECM.dI_dt = zeros(size(CSVcalc.VoltageDischarge_time));              % Initialize an array to store the ECM.dI_dt

for i = 2:length(CSVcalc.VoltageDischarge_time)
    
    % Calculate the simulated voltage at each time point using the equivalent circuit model equation
    ECM.dI_dt(i) = (CSVcalc.CurrentDischarge(i) - CSVcalc.CurrentDischarge(i-1)) / (CSVcalc.VoltageDischarge_time(i) - CSVcalc.VoltageDischarge_time(i-1));             % Approximate derivative of I with respect to t
    ECM.simulated_voltage(i) = OCV_avg - CIT.R1 * CSVcalc.CurrentDischarge(i) - CIT.R2(i) * CSVcalc.CurrentDischarge(i) * CSVcalc.SoCDischarge(i)/100 + CIT.C_dl(i) * ECM.dI_dt(i);

end

% Compare the simulated voltage with the experimental data
figure;
plot(CSVcalc.VoltageDischarge_time, CSVcalc.VoltageDischarge, 'b', CSVcalc.VoltageDischarge_time, ECM.simulated_voltage, 'r');
ylim([0 6]);
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('Experimental Data', 'Simulated Voltage');
