%TET_TSFC_THRUST function calculates the TSFC and net uninstalled thrust of
%a turbofan given fixed conditions indicated below, taking Turbine entry
%temperature, high pressure compressor ratio as input parameters.

%SID: 201029633 - Adam Gregory

%To use the program, vary the input parameters below, and set the inital P0
%and T0 for the flight conditions 

function [ThrustTotal,TSFC] = TET_TSFC_THRUST(TET,PRhp,n)

%-------------------------Input Parameters -----------------------------
yCold = 1.4;%Cold air properties - capacity ratio, heat capacity and gas constant
CpCold = 1005;
RCold = 287;

yHot = 1.333;%Hot Air properties
CpHot= 1150;
RHot = 287.3;

M0 = 0.8;%Mach number of flight
BPR = 6.64;%Bypass ratio

PRintake = 0.96;%Intake pressure recovery

PRCritCold = ((yCold+1)./2).^(yCold./(yCold-1));%Critial pressure ratio for choking cold air
PRlp = PRCritCold;%Low pressure compressor pressure ratio
LPCEff =0.93;%Low pressure compressor isentropic efficiency

HPCEff = 0.91;%High pressure compressor isentropic efficiency

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

%-------------------------------------------------------------------------
%---------------------- Process Calculations------------------------------

% 0 - Initial Stage - Using ISA standard values for elevation of 9.67km
P0= 100000.* (0.2858 + ((0.2858-0.2650) .* (17./50)));
T0=226.5 + (223.3-226.5).*(17./50);
V0 = M0 .* sqrt(yCold.*RCold.*T0);

%Total Values - calculate the total pressure and temperature for cruise
%conditions
Pt0 = P0 * (1+ ((M0^2)*((yCold-1)/2)))^(yCold/(yCold-1));
Tt0 = TotalTempIdeal(T0,Pt0,P0,yCold);

% 1 - Intake stage
Pt1 = PRintake * Pt0;
Tt1 = Tt0;

% 2 - Low Pressure Compressor Stage
Pt2 = PRlp .* Pt1;
Tt2_ = TotalTempIdeal(Tt1,Pt2,Pt1,yCold);
Tt2 = (Tt2_ - Tt1)./LPCEff + Tt1;

% 3 - High pressure Compressor Stage
Pt3 = PRhp .* Pt2;
Tt3_ = TotalTempIdeal(Tt2,Pt3,Pt2,yCold);
Tt3 = (Tt3_ - Tt2)./HPCEff + Tt2;

% 4 - Combuster Stage
Tt4= TET;
Pt4 = PRcomb * Pt3;

AFR = ((LHV.*CombEff)./(CpHot.*(Tt4-Tt3)))-1;

mAirIntake = n .* 408;%Mass flow rate of air at intake
mCore = mAirIntake./(1+BPR);
mFuel = mCore./AFR;

% 5 - High Pressure Turbine stage
Tt5 = Tt4 - ((AFR .* CpCold .* (Tt3 - Tt2).*HPSpoolEff)./((AFR+1).*CpHot));
Tt5_ = Tt4 - ((Tt4-Tt5)./HPTEff);

Pt5 = Pt4 * (Tt5_./Tt4).^(yHot./(yHot-1));

% 6 - Low Pressure Turbine stage
Tt6 = Tt5 - (((BPR+1).*AFR.*CpCold.*(Tt2-Tt1).*LPSpoolEff)./((AFR+1).*CpHot));
Tt6_ = Tt5 - ((Tt5-Tt6)./LPTEff);

Pt6 = Pt5 .* (Tt6_./Tt5).^(yHot./(yHot-1));

% 7 - Jet pipe stage
Pt7 = PRJet .* Pt6;
Tt7 = Tt6;

% 8 - Core Nozzle stage
Pt8 = Pt7;
Tt8 = Tt7;

PRCritHot = (((yHot)+1 )/2)^(yHot/(yHot-1)); %Critical pressure ratio for choking to occur
TRCritHot = (yHot+1)/2;%Temperature ratio when chokin occurs

PRIdealHot = Pt8./P0;

   

 P8(PRIdealHot < PRCritHot) = P0;
 T8(PRIdealHot < PRCritHot) = Tt8(PRIdealHot < PRCritHot)./((Pt8(PRIdealHot < PRCritHot)./P0).^((yHot-1)./yHot));
 V8(PRIdealHot < PRCritHot) = (2.*CpHot.*(Tt8(PRIdealHot < PRCritHot)-T8(PRIdealHot < PRCritHot))).^(1/2);
    
 P8(PRIdealHot >= PRCritHot) = Pt8(PRIdealHot >= PRCritHot)./PRCritHot;
 T8(PRIdealHot >= PRCritHot) = Tt8(PRIdealHot >= PRCritHot)./TRCritHot;
 V8(PRIdealHot >= PRCritHot) = (yHot.*RHot.*T8(PRIdealHot >= PRCritHot)).^(1/2);

V8 = real(V8); %Removes imaginary and thus invalid values

A8 = ((mCore).*RHot.*T8)./(P8.*V8);
ThrustHot = (mCore + mFuel).*V8 - mCore.*V0 + (A8.*(P8-P0));

% 9 - Bypass Nozzle stage - As this flow is taken to be choked (and
% utilises the critical pressure ratio), no conditions need to be taken
% into account, and pressure and temperature simply need to be divided by
% critical ratios
Tt9 = Tt2;
Pt9 = PRBypass * Pt2;

TRCritCold = (yCold+1)/2;

P9 = Pt9./PRCritCold;
T9 = Tt9./TRCritCold;
V9 = (yCold*RCold.*T9).^(1/2);

mBypass = mCore.*BPR;

A9 = ((mBypass).*RCold.*T9)./(P9.*V9);
ThrustCold = (mBypass).*V9 - mBypass.*V0+(A9.*(P9-P0));

%------------------------------FINAL VALUES-------------------------------
ThrustTotal = ThrustHot + ThrustCold;
TSFC = mFuel./ThrustTotal;
end