%Q4_Calculated calculates the value of TSFC and Net uninstalled thrust for
%a non mixing turbofan engine at specified TET and HPC pressure ratios

%SID: 201029633 - Adam Gregory

%To alter paramters more finely, adjust the function TET_TSFC_THRUST
%accordingly, as therea are many assumed paramters within the function that
%can be changed.

clc
clear
close all

n=1.14;%normalised spool speed
TET = 1050:.1:2000;%Turbine entry temperature
PRhpc = [0.9*26.7 .* (n .^ 1.22),26.7 .* (n .^ 1.22),1.1*26.7 .* (n .^ 1.22)];%10% increase and decrease in the HPC pressure ratio

[ThrustTotal(1,:),TSFC(1,:)] = TET_TSFC_THRUST(TET,PRhpc(1),n);
[ThrustTotal(2,:),TSFC(2,:)] = TET_TSFC_THRUST(TET,PRhpc(2),n);
[ThrustTotal(3,:),TSFC(3,:)] = TET_TSFC_THRUST(TET,PRhpc(3),n);

%------------------------------FINAL VALUES-------------------------------

minTSFC = [TET(TSFC(1,:) == nanmin(TSFC(1,:))),nanmin(TSFC(1,:));
 TET(TSFC(2,:) == nanmin(TSFC(2,:))),nanmin(TSFC(2,:));
 TET(TSFC(3,:) == nanmin(TSFC(3,:))),nanmin(TSFC(3,:))];

figure(1)
hold on
title('Plot of Thrust against TET');
plot(TET,ThrustTotal(1,:),'blue')
plot(TET,ThrustTotal(2,:),'red')
plot(TET,ThrustTotal(3,:),'green')
lgd = legend(num2str(round(PRhpc(1),1)),num2str(round(PRhpc(2),1)),num2str(round(PRhpc(3),1)));
title(lgd,{'High Pressure Compressor'  'Pressure Ratio'})
axis([min(TET) max(TET)+50 min(ThrustTotal(1,:))*0.9 max(ThrustTotal(1,:))*1.1])
xlabel('TET (K)');
ylabel('Thrust (N)');


figure(2)
hold on
title('Plot of TSFC against TET');
plot(TET,TSFC(1,:),'blue')
plot(TET,TSFC(2,:),'red')
plot(TET,TSFC(3,:),'green')
plot(minTSFC(1,1),minTSFC(1,2),'xb')
plot(minTSFC(2,1),minTSFC(2,2),'xr')
plot(minTSFC(3,1),minTSFC(3,2),'xg')
labels = {"\leftarrow"+num2str(minTSFC(1,1))+","+num2str(round(minTSFC(1,2),3,'significant')),"\leftarrow"+num2str(minTSFC(2,1))+","+num2str(round(minTSFC(2,2),3,'significant')),"\leftarrow"+num2str(minTSFC(3,1))+","+num2str(round(minTSFC(3,2),3,'significant'))};
text(minTSFC(:,1),minTSFC(:,2),labels,'HorizontalAlignment','left')
axis([min(TET) max(TET)+50 0.9*min(TSFC(1,:)) 1.1*max(TSFC(1,:))])
lgd = legend(num2str(round(PRhpc(1),1)),num2str(round(PRhpc(2),1)),num2str(round(PRhpc(3),1)));
title(lgd,{'High Pressure Compressor'  'Pressure Ratio'})
xlabel('TET (K)');
ylabel('TSFC (Kg.s/N)');
