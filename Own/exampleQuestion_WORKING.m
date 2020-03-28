clc
clear
%TurbineAnalysis.m takes inputs

h = 9670; 
BPR = 0.8;
n = 1.14;
M0 = 0.8;

yCold = 1.4;
CpCold = 1005;
RCold = 287;

yHot = 1.333;
CpHot=1150;
RHot = 287.3;

mAirIntake = n * 408;
PRintake = 1;


PRCritCold = ((yCold+1)/2)^(yCold/(yCold-1));

PRlp = 2.5; %Calculated low pressure compressor pressure ratio to be limit of choking
LPCEff =0.84;

PRhp = 3.5;
HPCEff = 0.86;

PRcomb = 0.92;
CombEff=0.98;

LHV = 43.1 * 10^6;

HPTEff = 0.9;
HPSpoolEff = 1;

LPTEff = 0.9;
LPSpoolEff = 1;

PRJet = 1;
PRCore = 1;
PRBypass = 0.95;

TET = 1250;

%Calculating pressures for the altitude
P0=101300;
T0=288.2;


V0 = M0 * sqrt(yCold*RCold*T0);


%Total pressure and temperature at intake
Pt0 = P0 * (1+ ((M0^2)*((yCold-1)/2)))^(yCold/(yCold-1));
Tt0 = T0 *((Pt0/P0)^((yCold-1)/yCold));

%Intake stage
Pt1 = PRintake * Pt0;
Tt1 = Tt0;

%LPCompressor Stage
Pt2 = PRlp .* Pt1;
Tt2_ = TotalTempIdeal(Tt1,Pt2,Pt1,yCold);
Tt2 = (Tt2_ - Tt1)./LPCEff + Tt1;

%HPCompressor Stage
Pt3 = PRhp .* Pt2;
Tt3_ = TotalTempIdeal(Tt2,Pt3,Pt2,yCold);
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

indexChoked = find(PRIdealHot >= PRCritHot);
indexNotChoked = find(PRIdealHot < PRCritHot);

P8(indexNotChoked) = P0;
T8(indexNotChoked) = Tt8(indexNotChoked)./((Pt8(indexNotChoked)./P0).^((yHot-1)./yHot));
V8(indexNotChoked) = (2.*CpHot.*(Tt8(indexNotChoked)-T8(indexNotChoked))).^(1/2);


P8(indexChoked) = Pt8(indexChoked)./PRCritHot;
T8(indexChoked) = Tt8(indexChoked)./TRCritHot;
V8(indexChoked) = (yHot.*RHot.*T8(indexChoked)).^(1/2);

V8 = real(V8);

A8 = ((mCore).*RHot.*T8)./(P8.*V8);
% ThrustHot = (mCore).*V8 - mCore.*V0+(A8.*(P8-P0));

ThrustHot = (mCore + mFuel).*V8 - mCore.*V0 + (A8.*(P8-P0));


%Bypass Nozzle stage - KNOWN CHOKED FLOW
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

ThrustHot./mCore
ThrustCold./mCore