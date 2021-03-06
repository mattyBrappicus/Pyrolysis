clear
%clc
%Species key:
    % 1 = c2h4
    % 2 = hcl
    % 3 = vinylCl
    % 4 = 1,1,2-trichloroethane
    % 5 = h2
    % 6 = cl2
    % 7 = 1,2-dichloroethane
    % 8 = c4h6
    % 9 = c2h2
    % 10 = c2h2cl2

%total heats of reaction, calclated in PyrolysisHRxn spreadsheet
%H_tot = H_rxn,std + deltaCp, T=773K
H_tot1 = 92.90; % units of kJ/mol
H_tot2 = 77.63;
H_tot3 = 186.0;
H_tot4 = 167.2;
H_tot5 = -173.5;
H_tot6 = 62.32;
H_tot = [H_tot1, H_tot2, H_tot3, H_tot4, H_tot5, H_tot6]; 

%Average heat capacities (Cpbar), Temp range=298K-773K
Cp_1 = 0.0640; % units of kJ/mol/K
Cp_2 = 0.0295; 
Cp_3 = 0.0215;
Cp_4 = 0.1175;
Cp_5 = 0.0293;
Cp_6 = 0.0359;
Cp_7 = 0.1063;
Cp_8 = 0.1231;
Cp_9 = 0.0547;
Cp_10 = 0.0651;
Cp_tot = [Cp_1 Cp_2 Cp_3 Cp_4 Cp_5 Cp_6 Cp_7 Cp_8 Cp_9 Cp_10];
sumCp = sum(Cp_tot);
 

%Main Properties, note that some of these values were taken from Aspen HYSYS
T0 = 491; %                 units of K, =dew point of stream S8 in HYSYS
P0 = 3000; %                units of kPa
D = 0.05; %               units of m; diameter of tube
L = 9;  %                 units of m
N = 1100; %                 number of tubes
Ac = (pi*((D^2)/4)); %      units of m^2
mu = 1.591*10^-5; %         units of kg/m/s
rho0 = 63; %              units of kg/m^3
V_r = (pi*((D^2)/4))*L; %   units of m^3

%Coolant Properties
U = 0.3; % units of kJ/(m^2*K*s)
Tc0 =560; % units of K, boiling point is 530 K
flowC = 10010/3600; % units of kg/s

%Initial molar flowrates from starting material balance
 % units of mol/s
F1_0 = 0.3732/3600; %   1 = c2h4
F2_0 = 2.2637/3600; %   2 = hcl
F3_0 = 0.0001/3600; %   3 = vinylCl
F4_0 = 3.7100/3600; %   4 = 1,1,2-trichloroethane
F5_0 = 0/3600; %        5 = h2
F6_0 = 124.5069/3600; % 6 = cl2
F7_0 = 1801.084/3600; % 7 = 1,2-dichloroethane
F8_0 = 0/3600; %        8 = c4h6
F9_0 = 0/3600; %        9 = c2h2
F10_0 = 0/3600; %       10 = c2h2cl2
F0 = [F1_0 F2_0 F3_0 F4_0 F5_0 F6_0 F7_0 F8_0 F9_0 F10_0];     
Ftotal_0 = sum(F0);

%Pressure drop initial parameters
MW = [0.02805, 0.03646, 0.06250, 0.1334, 0.00202, 0.0709, 0.09896, 0.05409, 0.026038, 0.09694]; %kg/mol
Q0 = sum(MW.*F0)/rho0; % units of m^3/s
Re0 = 4*rho0*Q0/pi/D/mu;
if Re0 < 10000
    frick0 = 16/Re0;
else
    frick0 = 1/(12.96*(log10(6.9/Re0))^2);
end
Beta0 = -32*frick0*rho0*Q0^2/Ac/pi^2/D^5; % units of kPa/m^3

%Logic
numElements = 50; % number of solver iterations
dv = V_r/numElements;
vspan = linspace(0, V_r, numElements);
y0 = [F1_0 F2_0 F3_0 F4_0 F5_0 F6_0 F7_0 F8_0 F9_0 F10_0 T0 P0 Tc0]; % load dependent variables
handleranon = @(v,y) handler(v,y,H_tot,Cp_tot,L,D,Ac,U,flowC,Ftotal_0,T0,P0,rho0,MW,mu); % use handler fxn
[ v, ysoln ] = ode15s(handleranon,vspan,y0);
conv = zeros(numElements,1);
for i = 1:numElements
    conv(i) = (1-ysoln(i,7)/ysoln(1,7));
end
disp(conv(numElements))
finalvcl = ysoln(numElements,3)/(3006*.3476)*3600;

%disp('Final Conversion: '+ num2str(conv(numElements)))
plotdata(v, ysoln, conv);
