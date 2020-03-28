function ANS = PRhpc2FN(PRhpc)
% Input known variables
massflow_fuel = 0;
cp_cold = 1005;
R_cold = 287;

cp_hot = 1150;
R_hot = 287.3;

gamma_cold = 1.4;
gamma_hot = 1.33;

LHV = 43.1*10^6;

P0 = 28.1384*10^3;
T0 = 225.815;
M0 = 0.8;
PRin = 0.96;
TRin = 1;

PRlpc = 1.437;
LPC_efficiency = 0.93;
BPR = 6.3;

HPC_efficiency = 0.91;

TET=1050;
comb_efficiency = 0.99;
ploss_comb = 0.02;

HPT_efficiency = 0.94;

LPT_efficiency = 0.95;

ploss_bypass = 0.1;

ploss_jetpipe = 0;
tloss_jetpipe = 0;

massflow_air = 420.24;

% Intake stage (0-1)

Pt0 = ps2pt(P0, gamma_cold, M0);
Tt0 = ts2tt(T0, Pt0, P0, gamma_cold);

Pt1 = PRin*Pt0;
Tt1 = TRin*Tt0;

V0 = M0*(gamma_cold*R_cold*T0)^(1/2);

% LP compressor stage (1-2)

Pt2 = PRlpc*Pt1;
Tt2 = tcompression(Tt1, Pt1, Pt2, gamma_cold, LPC_efficiency);

%HP compressor stage (2-3)

Pt3 = PRhpc*Pt2;
Tt3 = tcompression(Tt2, Pt2, Pt3, gamma_cold, HPC_efficiency);

% Combustor stage (3-4)

PRcc = 1-ploss_comb;
Pt4 = PRcc*Pt3;

AFR = ((LHV*comb_efficiency)/(cp_hot*(TET-Tt3)))-1;
massflow_core = massflow_air/(1+BPR);

if massflow_fuel == 0
    massflow_fuel = massflow_core/AFR;
end

%HP turbine stage (4-5)

Tt5 = TET - ((AFR*cp_cold*(Tt3-Tt2))/((AFR+1)*cp_hot));
Pt5 = pexpansion(Pt4, TET, Tt5, gamma_hot, HPT_efficiency);

% LP Turbine stage (5-6)

Tt6 = Tt5 - ((1+BPR)*AFR*cp_cold*(Tt2-Tt1))/((AFR+1)*cp_hot);
Pt6 = pexpansion(Pt5, Tt5, Tt6, gamma_hot, LPT_efficiency);

%Jet pipe stage (6-7)

Pt7 = (1-ploss_jetpipe)*Pt6;
Tt7 = (1-tloss_jetpipe)*Tt6;

% Core nozzle stage (7-8)

Pt8 = Pt7;
Tt8 = Tt7;

PRcrit = 1.852;
TRcrit = 1.167;

PRideal = Pt8/P0;

if PRideal >= PRcrit
    M8 = 1;
    P8 = Pt8/PRcrit;
    T8 = Tt8/TRcrit;
    V8 = M8*(gamma_hot*R_hot*T8)^(1/2);
    a8_mcore = (1+1/AFR)*(R_hot*T8)/(P8*V8);
    FN_hot = ((1+1/AFR)*V8-V0+a8_mcore*(P8-P0))*massflow_core;
else
    % If the flow is not choked, we assume the flow is fully expanded
    % Therefore, Pe = P0
    P8 = P0;
    T8 = Tt8/((Pt8/P0)^((gamma_hot-1)/gamma_hot));
    V8 = (2*cp_hot*(Tt8-T8))^(1/2);
    FN_hot=massflow_core*V8-massflow_air*V0;
end



% Bypass nozzle stage (2-9)

Tt9 = Tt2;
T9 = Tt9/TRcrit;
V9 = (gamma_cold*R_cold*T9)^(1/2);
Pt9 = (1-ploss_bypass)*PRlpc*Pt1;
P9 = Pt9/PRcrit;

a9_mcore = BPR*(R_hot*T8)/(P8*V8);
FN_cold = massflow_core*(BPR*(V9-V0)+a9_mcore*(P9-P0));

% NET UNINSTALLED THRUST

FN_total = FN_hot+FN_cold;
TSFC = massflow_fuel/FN_total;
ANS = [FN_total; TSFC];
end