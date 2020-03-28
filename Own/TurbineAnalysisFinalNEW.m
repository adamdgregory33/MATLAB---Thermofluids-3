clc
clear
close all

% --------------------------- Input Parameters -----------------------------

M0 = 0.8;%Mach number of aircraft

yCold = 1.4;%Cold air properties - capacity ratio, heat capacity and gas constant
CpCold = 1005;
RCold = 287;

yHot = 1.333;%Hot Air properties
CpHot= 1150;
RHot = 287.3;

PRintake = 0.96;%Intake pressure recovery

PRCritCold = ((yCold+1)./2).^(yCold./(yCold-1));%Critial pressure ratio for choking cold air
PRlp = PRCritCold;%Low pressure compressor pressure ratio - it is known to be CHOKED
LPCEff =0.93;%Low pressure compressor isentropic efficiency

HPCEff = 0.91;%High pressure compressor isentropic efficiency
BPR = 6.64;%Bypass ratio
n = 1.14;%Normalised spool speed
PRhp = [0.9.*26.7 .* (n .^ 1.22),26.7 .* (n .^ 1.22), 1.1*26.7 .* (n .^ 1.22)];%High pressure compressor pressure ratio for an increase and decrease of
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

TET = 940:.1:3000;%Turbine entry temperature

%Calculating pressures for the altitude
P0= 100000.* (0.2858 + ((0.2858-0.2650) .* (17./50)));
T0=226.5 + (223.3-226.5).*(17./50);
V0 = M0 .* sqrt(yCold.*RCold.*T0);

%Total pressure and temperature at intake
Pt0 = P0 .* (1+ ((M0.^2).*((yCold-1)./2))).^(yCold./(yCold-1));
Tt0 = T0 .*((Pt0./P0).^((yCold-1)./yCold));

%Intake stage
Pt1 = PRintake .* Pt0;
Tt1 = Tt0;

%Low pressure compressor stage
Pt2 = PRlp .* Pt1;
Tt2_ = Tt1.* (Pt2./Pt1).^((yCold-1)./yCold);
Tt2 = (Tt2_ - Tt1)./LPCEff + Tt1;

%High pressure compressor stage
Pt3 = PRhp .*Pt2;
Tt3_ = Tt2.* (Pt3./Pt2).^((yCold-1)./yCold);
Tt3 = (Tt3_ - Tt2)./HPCEff + Tt2;

%Combustor section
Tt4 = TET;
Pt4 = PRcomb .* Pt3;

AFR = [((LHV .* CombEff)./(CpHot.*(Tt4-Tt3(1))))-1,
    ((LHV .* CombEff)./(CpHot.*(Tt4-Tt3(2))))-1,
    ((LHV .* CombEff)./(CpHot.*(Tt4-Tt3(3))))-1];

mAirIntake = n .* 408;%Mass flow rate of air at intake
mCore = mAirIntake./(1+BPR);%mass flow of air rate into core of engine
mBypass = mCore.*BPR;%mass flow rate of air into engine bypass
mFuel = mCore./AFR(1,:);%Mass flow rate of the fuel into engine

%High pressure turbine
Tt5 = [Tt4 - ((AFR(1,:) .* CpCold .* (Tt3(1) - Tt2).*HPSpoolEff)./((AFR(1,:)+1).*CpHot)),
    Tt4 - ((AFR(2,:) .* CpCold .* (Tt3(2) - Tt2).*HPSpoolEff)./((AFR(2,:)+1).*CpHot)),
    Tt4 - ((AFR(3,:) .* CpCold .* (Tt3(3) - Tt2).*HPSpoolEff)./((AFR(3,:)+1).*CpHot))];
Tt5_ = Tt4 - ((Tt4-Tt5)./HPTEff);

Pt5 = [Pt4(1) .* (Tt5_(1,:)./Tt4).^(yHot/(yHot-1)),
       Pt4(2) .* (Tt5_(2,:)./Tt4).^(yHot/(yHot-1)),
       Pt4(3) .* (Tt5_(3,:)./Tt4).^(yHot/(yHot-1))];

%Low pressure turbine
Tt6 = Tt5 - (((BPR+1).*AFR.*CpCold.*(Tt2-Tt1).*LPSpoolEff)./((AFR+1).*CpHot));
Tt6_ = Tt5 - ((Tt5-Tt6)./LPTEff);

Pt6 = Pt5 .* (Tt6_./Tt5).^(yHot./(yHot-1));

%Jet pipe section
Tt7 = Tt6;
Pt7 = Pt6;

%Core Jet Nozzle 
Pt8 = Pt7;
Tt8 = Tt7;

PRCritHot = (((yHot)+1 )./2).^(yHot./(yHot-1));
TRCritHot = (yHot+1)./2;

PRIdealHot = Pt8./P0;

%FOr CHoked Scenarios
P8 = Pt8./PRCritHot;
T8 = Tt8./TRCritHot;
V8 = sqrt(yHot.*RHot.*T8);


%For non choked scenarios
P8(PRIdealHot < PRCritHot) = P0;
T8(PRIdealHot < PRCritHot) = Tt8(PRIdealHot < PRCritHot)./ ((Pt8(PRIdealHot < PRCritHot)./P8(PRIdealHot < PRCritHot)).^((yHot-1)./yHot));
V8(PRIdealHot < PRCritHot) = sqrt(2.*CpHot.*(Tt8(PRIdealHot < PRCritHot)-T8(PRIdealHot < PRCritHot)));

V8 = real(V8);

% A8overMcore = (1+(1./AFR)).*((RHot.*T8)./(P8.*V8));
% ThrustHot = (1+(1./AFR)).*V8 - V0+A8overMcore.*(P8-P0);
A8 = ((mCore + mFuel).*RHot.*T8)./(P8.*V8);
ThrustHot = (mCore+mFuel).*V8 - mCore.*V0+(A8.*(P8-P0));


%Bypass Nozzle
Tt9 = Tt2;
Pt9 = PRBypass .* Pt2;

TRCritCold = (yCold+1)./2;

P9 = Pt9./PRCritCold;
T9 = Tt9./TRCritCold;
V9 = sqrt(yCold.*RCold.*T9);

% A9overMcore = BPR .* ((RCold.*T9)./(P9.*V9));
% ThrustCold = (BPR .*(V9-V0))+A9overMcore.*(P9-P0);
A9 = ((mBypass).*RCold.*T9)./(P9.*V9);
ThrustCold = (mBypass).*V9 - mBypass.*V0+(A9.*(P9-P0));

%TotalValues
ThrustTotal = ThrustHot + ThrustCold;
TSFC = mFuel./ThrustTotal;

figure(1)
hold on
title('Plot of Thrust against TET');
plot(TET,ThrustTotal(1,:))
plot(TET,ThrustTotal(2,:))
plot(TET,ThrustTotal(3,:))
lgd = legend(num2str(round(PRhp(1),1)),num2str(round(PRhp(2),1)),num2str(round(PRhp(3),1)));
title(lgd,{'High Pressure Compressor'  'Pressure Ratio'})
axis([700 3000 0 200000])
xlabel('TET');
ylabel('Thrust (N)');


figure(2)
hold on
title('Plot of TSFC against TET');

plot(TET,TSFC(1,:))
plot(TET,TSFC(2,:))
plot(TET,TSFC(3,:))
axis([700 3000 0 0.00005])
lgd = legend(num2str(round(PRhp(1),1)),num2str(round(PRhp(2),1)),num2str(round(PRhp(3),1)));
title(lgd,{'High Pressure Compressor'  'Pressure Ratio'})
xlabel('TET');
ylabel('TSFC');



TET(TSFC == min(TSFC(1,:)))
TET(TSFC == min(TSFC(2,:)))
TET(TSFC == min(TSFC(3,:)))

