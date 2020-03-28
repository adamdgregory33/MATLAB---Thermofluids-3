% Answer to Q3

% Engine performance model to automatise the gas turbine thermodynamic
% cycle analysis of a turbofan engine

% Outputs the total pressure and temperature of the flow at each point of the system
% labeled as:
% 0: Before intake                      5: Before low pressure turbine
% 1: Before low pressure compressor     6: Before jet pipe   
% 2: Before high pressure compressor    7: Before exhaust 
% 3: Before combustor                   8: Exhaust gases  
% 4: Before high pressure turbine       9: Bypass exhaust gases.

% These results are displayed in the command window.

% Author : Liam Bradbury
% SID : 201008735

clear;
clc;

% Input known variables
massflow_fuel = 0; %Mass flow rate of fuel (0 if not assigned)

% Properties of air
cp_cold = 1005;
R_cold = 287;

cp_hot = 1150;
R_hot = 287.3;

gamma_cold = 1.4;
gamma_hot = 1.33;

LHV = 43.1*10^6; % Lower Heating Value (J/kg) Calorific value

P0 = 28.1384*10^3; % Ambient static pressure of the surroundings
T0 = 225.815; % Ambient static temperature of the surroundings
M0 = 0.8; % Cruising Speed
PRin = 0.96; % Pressure ratio of the intake
TRin = 1; % Temperature ratio of the intake

PRlpc = 1.437; % Pressure ratio of the low pressure compressor
LPC_efficiency = 0.93; % Efficiency of the low pressure compressor
BPR = 6.3; % Bypass pressure ratio

PRhpc = 27.68042; % Pressure ratio of the high pressure compressor
HPC_iso_efficiency = 0.91; % Isentropic efficiency of the high pressure compressor

TET = 1050; % Turbine Entry Temperature
comb_efficiency = 0.99; % Combustion efficiency
ploss_comb = 0.02; % Pressure loss in the combustor

HPT_iso_efficiency = 0.94; % Isentropic efficiency of the high pressure turbine
HPT_mech_efficiency = 0.99; % Mechanical efficiency of the high pressure turbine


LPT_iso_efficiency = 0.95; % Isentropic efficiency of the low pressure turbine
LPT_mech_efficiency = 0.99; % Mechanical efficiency of the low pressure turbine

ploss_bypass = 0.1; % Pressure loss in the bypass

ploss_jetpipe = 0; % Pressure loss in the jetpipe
tloss_jetpipe = 0; % Temperature loss in the jetpipe

massflow_air = 420.24; % Relative mass flow of air at cruising speed and altitude

% #####################################
% #########Intake stage (0-1)##########
% #####################################

% Calculate the total pressure at point 0. (see function description)
Pt0 = ps2pt(P0, gamma_cold, M0);
% Calculate the total temperature at point 0. (see function description)
Tt0 = ts2tt(T0, Pt0, P0, gamma_cold);

% Calculate the total pressure at point 1 from the pressure ratio given
Pt1 = PRin*Pt0;
% Calculate the total temperature at point 1 from the temperature ratio
% given
Tt1 = TRin*Tt0;

% Convert mach number to velocity
V0 = M0*(gamma_cold*R_cold*T0)^(1/2);

% ###################################################
% ######### Low Pressure Compressor stage (1-2)######
% ###################################################

% Calculate the total pressure at point 2 from the pressure ratio of the
% low pressure compressor given.
Pt2 = PRlpc*Pt1;

% Calculate the total temperature at point 2 from isentropic relationships
% and the isentropic efficiency of the low pressure compressor given.
Tt2 = tcompression(Tt1, Pt1, Pt2, gamma_cold, LPC_efficiency);

% ###################################################
% ######### High Pressure Compressor stage (2-3)######
% ###################################################

% Calculate the total pressure at point 3 from the pressure ratio of the
% high pressure compressor given.
Pt3 = PRhpc*Pt2;

% Calculate the total temperature at point 3 from isentropic relationships
% and the isentropic efficiency of the high pressure compressor given.
Tt3 = tcompression(Tt2, Pt2, Pt3, gamma_cold, HPC_iso_efficiency);

% #####################################
% #########Combustor stage (3-4)#######
% #####################################

PRcc = 1-ploss_comb;

% The Turbine Entry Temperature is the exit temperature of the combustor
% Tt4=TET

% Calculate the total pressure at point 4 from the pressure loss of in the
% combustor given.
Pt4 = PRcc*Pt3;

% Air to Fuel Ratio derived from an energy balance. (Q_dot =
% massflow_fuel*LHV*effiency_of_Combustor)
AFR = ((LHV*comb_efficiency)/(cp_hot*(TET-Tt3)))-1;

% Calculates the mass flow rate within the core as the flow had split
% between the bypass and the core after the Low Pressure Compressor.
massflow_core = massflow_air/(1+BPR);

% If no mass flow rate of the fuel is given, calculate its value from the
% AFR
if massflow_fuel == 0
    massflow_fuel = massflow_core/AFR;
end

% ###################################################
% ######### High Pressure Turbine stage (4-5)########
% ###################################################

% Calculate the total temperature at point 5 from an energy balance of the
% high pressure spool. The energy produced by the high pressure turbine
% must be the energy required by the high pressure compressor.
Tt5 = TET - ((AFR*cp_cold*(Tt3-Tt2))/((AFR+1)*cp_hot*HPT_mech_efficiency));

% Total pressure at point 5 is calculated for isentropic equations, and the
% isentropic efficiency coefficient given.
Pt5 = pexpansion(Pt4, TET, Tt5, gamma_hot, HPT_iso_efficiency);

% ###################################################
% ######### Low Pressure Turbine stage (5-6)########
% ###################################################

% The same principles are applied here as in the high pressure turbine.
% The low pressure turbine must supply the energy required for the low
% pressure compressor.
Tt6 = Tt5 - ((1+BPR)*AFR*cp_cold*(Tt2-Tt1))/((AFR+1)*cp_hot*LPT_mech_efficiency);
Pt6 = pexpansion(Pt5, Tt5, Tt6, gamma_hot, LPT_iso_efficiency);

% ######################################
% ######### Jet pipe stage (6-7)########
% ######################################

Pt7 = (1-ploss_jetpipe)*Pt6;
Tt7 = (1-tloss_jetpipe)*Tt6;

% ######################################
% ######### Core nozzle stage (7-8)#####
% ######################################

% There are no pressure or temperature losses given so we assumed that
% there aren't any.
Pt8 = Pt7;
Tt8 = Tt7;

% Critical pressure ratio at which the flow is choked for gamma_hot = 1.333
PRcrit = 1.852;

% Critical temperature ratio at which the flow is choked for gamma_hot =
% 1.333.
TRcrit = 1.167;

% The ideal pressure and temperature ratio is calculated
PRideal = Pt8/P0;
TRideal = Tt8/Tt0;

if PRideal >= PRcrit || TRideal >= TRcrit
    M8 = 1;
    P8 = Pt8/PRcrit;
    T8 = Tt8/TRcrit;
    V8 = M8*(gamma_hot*R_hot*T8)^(1/2);
    a8_mcore = (1+1/AFR)*(R_hot*T8)/(P8*V8);
    FN_hot = ((1+1/AFR)*V8-V0+a8_mcore*(P8-P0))*massflow_core;
    
else
    % If the flow is not choked, we assume the flow is fully expanded
    % Therefore, Pe = P0
    massflow_8 = massflow_core+massflow_fuel;
    P8 = P0;
    T8 = Tt8/((Pt8/P0)^((gamma_hot-1)/gamma_hot));
    V8 = (2*cp_hot*(Tt8-T8))^(1/2);
    
    % FN_hot is the Net Uninstalled Thrust produced by the hot nozzle (ie
    % the core nozzle).
    FN_hot=massflow_8*V8-massflow_core*V0;
end


% ######################################
% ######### Bypass nozzle stage (2-9)###
% ######################################

% No temperature ratio was given at the bypass so we assume there are no
% losses between the low pressure compressor and the bypass nozzle.
Tt9 = Tt2;

% We know the flow is choked.
T9 = Tt9/TRcrit;
V9 = (gamma_cold*R_cold*T9)^(1/2);
Pt9 = (1-ploss_bypass)*PRlpc*Pt1;
P9 = Pt9/PRcrit;

% Calculates mass flow rate of the bypass.
massflow_9 = BPR*massflow_core;

% FN_cold is the net uninstalled thrust produced by the cold nozzle (ie the
% bypass nozzle).
a9_mcore = (massflow_9*R_cold*T9/(P9*V9))/massflow_core;
FN_cold = massflow_core*(BPR*(V9-V0)+a9_mcore*(P9-P0));

% NET UNINSTALLED THRUST
FN_total = FN_hot+FN_cold;

% THRUST SPECIFIC FUEL CONSUMPTION
TSFC = massflow_fuel/FN_total;

fprintf('The total pressure at point 0 is: %i Pa \n', Pt0)
fprintf('The total temperature at point 0 is: %i K \n', Tt0)
fprintf('The total pressure at point 1 is: %i Pa \n', Pt1)
fprintf('The total temperature at point 1 is: %i K \n', Tt1)
fprintf('The total pressure at point 2 is: %i Pa \n', Pt2)
fprintf('The total temperature at point 2 is: %i K \n', Tt2)
fprintf('The total pressure at point 3 is: %i Pa \n', Pt3)
fprintf('The total temperature at point 3 is: %i K \n', Tt3)
fprintf('The total pressure at point 4 is: %i Pa \n', Pt4)
fprintf('The total temperature at point 4 is: %i K \n', TET)
fprintf('The total pressure at point 5 is: %i Pa \n', Pt5)
fprintf('The total temperature at point 5 is: %i K \n', Tt5)
fprintf('The total pressure at point 6 is: %i Pa \n', Pt6)
fprintf('The total temperature at point 6 is: %i K \n', Tt6)
fprintf('The total pressure at point 7 is: %i Pa \n', Pt7)
fprintf('The total temperature at point 7 is: %i K \n', Tt7)
fprintf('The total pressure at point 8 is: %i Pa \n', Pt8)
fprintf('The total temperature at point 8 is: %i K \n', Tt8)
fprintf('The total pressure at point 9 is: %i Pa \n', Pt9)
fprintf('The total temperature at point 9 is: %i K \n', Tt9)
fprintf('\nThe net uninstalled thrust is: %i \n', FN_total)
fprintf('The thrust specific fuel consumption is: %i \n', TSFC)