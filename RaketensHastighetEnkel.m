function [Vraket] = RaketensHastighetEnkel(Liter,Tryck,Massa)
%% Beräknar raketens hastighet under initiala fasen när vattnet pressas ut ur flaskan
clc
close all 

%% Konstanter
r = 0.01025;    %radie av utlåpshålet  [m]
cd = 1;         % Discharge coefficient, från 8 cm diameter till 2 cm, höftning 
Dw = 998;       % Densitet vatten  [kg/m^3]
Patm = 100600;  % Air pressure i psi 100500  [N/m^2]
Pb = Tryck*Patm;    % Övertryck i flaskan  [N/m^2]
V = Liter/1000;  % Volymen av vatten i raketen   [m^3]
%V = 0.0005;        % Volymen av vatten i raketen   [m^3]
Mw = Dw*V;       % Massan av vatten        [kg] 
%Mb = 0.107;     % Raketens tomma massa, 0.107 kg    [kg]   
g = 9.82;       % Gravitations konstanten [m/s2]
Mtot = (Massa + Mw)/2;  % Totala massan raket(medelvärde)  [kg]

deltaP = Pb-Patm  % Skillnad i trycket mellan atmosfär och flaska
A = pi*r^2 % Tvärsnittsarea flaskans mynning 
MwR = A*cd*(2*Dw*deltaP)^0.5   % Massan flow rate
Vw = MwR/(Dw*A);         % Utloppshastighet vatten
Thrust = MwR*Vw;           % Thrust
ThrustM = Thrust - (Mtot*g);     % medel thrust
a = ThrustM/Mtot         % Acceleration
t = Mw/MwR         % tiden då allt vatten är ute
Vraket = a*t     
