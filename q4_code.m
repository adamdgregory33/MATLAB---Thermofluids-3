% Answer to Q4
% Engine performance model to automatise the gas turbine thermodynamic
% cycle analysis of a turbofan engine

% Author : Liam Bradbury
% SID : 201008735

clear;
clc;


% Number of TET values to search for within the given bounds
increments =10000; 

% HPC Values to graph Thrust vs TET
PRhpc_initial = 27.68042;
PRhpc_lb = 27.68042-0.1*27.68042;
PRhpc_ub = 27.68042+0.1*27.68042;

% Upper and Lower Bound of TET to output Net uninstalled Thrust and TSFC
TET_lb = 1050; % Lower bound
TET_ub = 2000; % Upper bound

% Initiate search variables and matrices
TET = linspace(TET_lb, TET_ub, increments); 
ANS_TET_initial = zeros(increments, 2);
ANS_TET_lb = zeros(increments, 2);
ANS_TET_ub = zeros(increments, 2);

% Calculates the Thrust force and TSFC for all the values of matrix TET
for i=1:increments
    ANS_TET_initial(i,:)=TET2FN(PRhpc_initial ,TET(i));
    ANS_TET_lb(i,:)=TET2FN(PRhpc_lb ,TET(i));
    ANS_TET_ub(i,:)=TET2FN(PRhpc_ub ,TET(i));
end

[min_TSFC_initial, TSFC1] = min(ANS_TET_initial(:,2));
[min_TSFC_lb, TSFC2] = min(ANS_TET_lb(:,2));
[min_TSFC_ub, TSFC3] = min(ANS_TET_ub(:,2));

TSFC1 = TET(TSFC1);
TSFC2 = TET(TSFC2);
TSFC3 = TET(TSFC3);

TSFC_points = [TSFC1, TSFC2, TSFC3];
mins_TSFC = [min_TSFC_initial, min_TSFC_lb, min_TSFC_ub];

% Plot of the TET versus the total uninstalled thrust.
figure(1)
h = plot(TET, ANS_TET_lb(:,1), 'r', TET, ANS_TET_initial(:,1), 'b', ...
    TET, ANS_TET_ub(:,1), 'k');
xlabel('TET (K)')
ylabel('Net Uninstalled Thrust (N)')
title('Net Uninstalled Thrust vs TET')
legend('HPC=24.9124','HPC=27.68042','HPC=30.4485', 'location', 'southeast')
grid on
grid minor

% Plot of the TET versus the thrust specific fuel consumption
figure(2)
f = plot(TET, ANS_TET_lb(:,2), 'r', TET, ANS_TET_initial(:,2), 'b', ...
    TET, ANS_TET_ub(:,2), 'k',TSFC1, min_TSFC_initial, 'bx', ...
    TSFC2, min_TSFC_lb, 'rx', TSFC3, min_TSFC_ub, 'kx');
strValues = strtrim(cellstr(num2str([TSFC_points(:) mins_TSFC(:)],'(%d,%d)')));
text(TSFC_points,mins_TSFC,strValues,'VerticalAlignment','bottom','FontSize', 8);
xlabel('TET (K)')
ylabel('TSFC (kg/N)')
title('Thrust Specific Fuel Consumption VS TET')
legend('HPC=24.9124','HPC=27.68042','HPC=30.4485', 'location', 'southeast')
grid on
grid minor
