function outvar = handler(v,y,H_tot,Cp_tot,L,D,Delta,Ac,U,flowC,Ftotal_0,T0,P0,rho0)
%Species indices key:
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
    
% Unload variables
F1 = y(1); % units of mol/s
F2 = y(2);
F3 = y(3);
F4 = y(4);
F5 = y(5);
F6 = y(6);
F7 = y(7);
F8 = y(8);
F9 = y(9);
F10 = y(10);
T = y(11); % units of K
P = y(12); % units of kPa

Tc = y(13); % units of K
Ftotal = (F1 + F2 + F3 + F4 + F5 + F6 + F7 + F8 + F9 + F10);
As = D * pi * L; % units of m^2
Do = D + 2*.0036; % units of m
Cpc = ((Tc - 273) * .0029 + 1.5041 + 273)/1000; % units of kJ/(kg*K)

% Calculate partial pressures for each species 
pp = [0,0,0,0,0,0,0,0,0,0];
for i = 1:length(pp)
    pp(i) = y(i)/Ftotal*P;
    %c(i) = (y(i)/Ftotal)*(T0/T)*P; % units of kPa Old Method
end

% Rate expressions from paper
R_kinetics = 8.3144621; % units of J/(mol*K)

% a's are the pre-exponential factors from the Lakshmanan paper
% Laksh's units for 2nd order constants are m^3/kmol/s. Constants here are for m^3/mol/s.
a1f = 10^13.6; 
a1r = 0.3*10^6; 
a2f = 0.5*10^14; 
a2r = 0.37*10^6; 
a3f = 10^13; %there is no reverse reaction here
a4f = 10^11; 
a4r = 0.8*10^6; 
a5f = 0.15*10^6; 
a5r = 0.5*10^13; 
a6f = 0.2*10^14; %note that a6 constants are for Laksh's 18th pyrolysis rxn
a6r = 0.3*10^6; 

% E's are acivation energies from Lakshmanan.
E1f = -58000; % unsure of units
E1r = -44000; 
E2f = -69000; 
E2r = -40000; 
E3f = -72000; 
E4f = -82000; % note typo in Laksh for this Ea
E4r = -38000; 
E5f = -32000; 
E5r = -73000; 
E6f = -58000; 
E6r = -44000; 

% Equilibrium rate constants
k(1) = (a1f * exp(E1f/(1.987*T)))/(a1r * exp(E1r/(1.987*T))); % units=mol/m^3
k(2) = (a2f * exp(E2f/(1.987*T)))/(a2r * exp(E2r/(1.987*T))); % units=mol/m^3
k(3) = (a3f * exp(E3f/(1.987*T))); %                            units=mol/m^3
k(4) = (a4f * exp(E4f/(1.987*T)))/(a4r * exp(E4r/(1.987*T))); % units=mol/m^3
k(5) = (a5f * exp(E5f/(1.987*T)))/(a5r * exp(E5r/(1.987*T))); % units=m^3/mol
k(6) = (a6f * exp(E6f/(1.987*T)))/(a6r * exp(E6r/(1.987*T))); % units=mol/m^3

% Rate equations (pp/RT = concentration in mol/m^3)
r1 = k(1) * pp(7) / R_kinetics / T;
r2 = k(2) * pp(3) / R_kinetics / T;
r3 = k(3) * pp(7) / R_kinetics / T;
r4 = k(4) * pp(1) / R_kinetics / T;
r5 = k(5) * pp(9) * pp(1) / R_kinetics / T;
r6 = k(6) * pp(4) / R_kinetics / T;

% Set up differential equations
outvar(1) = r3 - 1*r5;
outvar(2) = r1 + r2 + r6;
outvar(3) = r1 - 1*r2;
outvar(4) = -1*r6;
outvar(5) = r4;
outvar(6) = r3;
outvar(7) = -1*r1 - 1*r3;
outvar(8) = r5;
outvar(9) = r4 - 1*r5;
outvar(10) = r6;

% Tube side differential equation
term1 = r1 * H_tot(1);
term1 = term1 + r2 * H_tot(2);
term1 = term1 + r3 * H_tot(3);
term1 = term1 + r4 * H_tot(4);
term1 = term1 + r5 * H_tot(5);
term1 = term1 + r6 * H_tot(6);

term2 = U * 4 / D * (T - Tc); % units of kJ/(m^3*s)

% Calculate mol fractions
molFracs = [0,0,0,0,0,0,0,0,0,0];
for i = 1:length(molFracs)
    molFracs(i) = y(i)/Ftotal;
end

term3 = 0;
for i = 1:length(molFracs)
    term3 = term3 + y(i)*Cp_tot(i); %Cp in units of kJ/mol K
end

outvar(11) = (-term1 - term2) / term3; % units of K/m^3

% Ergun differetial equation
rho = rho0 * (P/P0) * (Ftotal_0/Ftotal) * (T0/T); % units of kg/m^3
outvar(12) = -Delta * (1/(Ac*rho)); % units of kPa/m^3


% Shell side differential equations
outvar(13) = -(4 * U * (T - Tc)) / (flowC * Cpc * Do); % units of K/m^3

outvar = outvar';

end
