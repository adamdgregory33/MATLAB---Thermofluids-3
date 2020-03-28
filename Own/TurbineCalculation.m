function [ThrustTotal, TSFC] = TurbineCalculation(TETInput)

yCold = 1.4;%Cold air properties - capacity ratio, heat capacity and gas constant
CpCold = 1005;
RCold = 287;

yHot = 1.333;%Hot Air properties
CpHot= 1150;
RHot = 287.3;

BPR = 6.64;%Bypass ratio
n = 1.14;%Normalised spool speed

M0 = 0.8;%Mach number of flight
PRintake = 0.96;%Intake pressure recovery

PRCritCold = ((yCold+1)./2).^(yCold./(yCold-1));%Critial pressure ratio for choking cold air
PRlp = PRCritCold;%Low pressure compressor pressure ratio
LPCEff =0.93;%Low pressure compressor isentropic efficiency

HPCEff = 0.91;%High pressure compressor isentropic efficiency
PRhp = 26.7 .* (n .^ 1.22);%High pressure compressor pressure ratio for an increase and decrease of
%10% about the provided value of n

PRcomb = 0.98;%Pressure ratio combustor
CombEff=0.99;%combuster efficiency

LHV = 43.1 .* 10.^6;%Lower Heating Value of fuel

HPTEff = 0.94;%high pressure turbine isentropic efficiency
HPSpoolEff = 0.99;%High pressure spool efficiency

LPTEff = 0.95;%Low pressure turbine isentropic efficiency
LPSpoolEff = 0.99;%low pressure spol efficiency

PRJet = 1;%Jet pressure recovery
PRCore = 1;%Core pressure recovery
PRBypass = 0.9;%Bypass pressure recovery

TET = TETInput;%Turbine entry temperature

%Initial Stage
P0= 100000.* (0.2858 + ((0.2858-0.2650) .* (17./50)));
T0=226.5 + (223.3-226.5).*(17./50);
V0 = M0 .* sqrt(yCold.*RCold.*T0);

%Total Values
Pt0 = P0 * (1+ ((M0^2)*((yCold-1)/2)))^(yCold/(yCold-1));
Tt0 = TotalTempIdeal(T0,Pt0,P0,yCold);

%Intake stage
Pt1 = PRintake * Pt0;
Tt1 = Tt0;

%LPCompressor Stage
Pt2 = PRlp .* Pt1;
Tt2_ = TotalTempIdeal(Tt1,Pt2,Pt1,yCold);
Tt2 = (Tt2_ - Tt1)./LPCEff + Tt1;

%HPCompressor Stage
Pt3 = PRhp .* Pt2;
Tt3_ = TotalTempIdeal(Tt2,Pt2,Pt1,yCold);
Tt3 = (Tt3_ - Tt2)./HPCEff + Tt2;

%Combuster Stage
Tt4= TET;
Pt4 = PRcomb * Pt3;

AFR = ((LHV.*CombEff)./(CpHot.*(Tt4-Tt3)))-1;

mAirIntake = n .* 408;%Mass flow rate of air at intake
mCore = mAirIntake./(1+BPR);
mFuel = mCore./AFR;

%HPTurbine stage
Tt5 = Tt4 - ((AFR .* CpCold .* (Tt3 - Tt2).*HPSpoolEff)./((AFR+1).*CpHot));
Tt5_ = Tt4 - ((Tt4-Tt5)./HPTEff);

Pt5 = Pt4 * (Tt5_./Tt4).^(yHot./(yHot-1));

%LPTurbine stage
Tt6 = Tt5 - (((BPR+1).*AFR.*CpCold.*(Tt2-Tt1).*LPSpoolEff)./((AFR+1).*CpHot));
Tt6_ = Tt5 - ((Tt5-Tt6)./LPTEff);

Pt6 = Pt5 .* (Tt6_./Tt5).^(yHot./(yHot-1));

%Jet pipe stage
Pt7 = PRJet .* Pt6;
Tt7 = Tt6;

%Core Nozzle stage
Pt8 = Pt7;
Tt8 = Tt7;

PRCritHot = (((yHot)+1 )/2)^(yHot/(yHot-1));
TRCritHot = (yHot+1)/2;

PRIdealHot = Pt8./P0;
TRIdealHot = Tt8./Tt0;

if PRIdealHot >= PRCritHot
    P8 = P0;
    T8 = Tt8./((Pt8./P0).^((yHot-1)./yHot));
    V8 = (2.*CpHot.*(Tt8-T8)).^(1/2);
    
else
    P8 = Pt8./PRCritHot;
    T8 = Tt8./TRCritHot;
    V8 = sqrt(yHot.*RHot.*T8);
    
end

A8 = ((mCore).*RHot.*T8)./(P8.*V8);
ThrustHot = (mCore).*V8 - mCore.*V0+(A8.*(P8-P0));

%ThrustHot = (mCore + mFuel).*V8 - mCore.*V0;


%Bypass Nozzle stage - known CHOKED FLOW
Tt9 = Tt2;
Pt9 = PRBypass * Pt2;

TRCritCold = (yCold+1)/2;

P9 = Pt9./PRCritCold;
T9 = Tt9./TRCritCold;
V9 = (yCold*RCold.*T9).^(1/2);

mBypass = mCore.*BPR;

A9 = ((mBypass).*RCold.*T9)./(P9.*V9);
ThrustCold = (mBypass).*V9 - mBypass.*V0+(A9.*(P9-P0));

%FINAL VALUES
ThrustTotal = ThrustHot + ThrustCold;
TSFC = mFuel./ThrustTotal;



end