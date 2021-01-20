function [Vb, time] = RaketensHastighet
% Beräknar raketens hastighet under initiala fasen när vattnet pressas ut
% ur flaskan

%% Tibert tipsade om att lägga på en buffert i 
%% framförallt luftdriftsfasen på kanske 10
clc
close all
clear

Patm = 1.2754; % kg/m^3
u = 0.1; % Friktionskonstanten från rampen
Angle = 45; % Startvinkel i grader
g = 9.82; % m/s^2
Mb = 0.107; % Raketens tomma massa, 0.1 kg
Ma = 0.001225*1.5; % Luftens massa, 1,225 g/l * 1,5 liter  
Mw = 1.565; % 1,565 liter vatten = 1,565 kg 
Pw = 997; % Vattnets densitet kg/m3
Pa0 = 2.5; % 2,5 Bar som en höftning
V = 0.001564; % Flaskans volym i kubikmeter
Vw0 = 0.28 * V; % Vattnets volym, 28% fylld flaska, 
FlaskaRadie = 0.0881/2; % Flaskans ytterdiameter
Sb = 2*pi*FlaskaRadie; % Cross sectional area of the rocket, 0.0881 m
S = 0.0205; % Innerdiameter flaskans mynning
Cd = 0.345; % Dragcoefficient
Ips = 1.4; % För monatomic ideal gas (wikipedia) 

Fm = u*(Mb + Ma + Mw)*g*cosd(Angle); % Friktion raket-ramp
Wx = (Mb + Ma + Mw)*g*sind(Angle); % Vik 
h = 0.01; %Steglängd
x = 0;
y = 0;
time = 0;
xMatrix = [];
yMatrix = [];

% Aktuellt luftmotstånd vid viss lufthastighet
Fd = @(xVb) 0.5*Patm*Cd*Sb*xVb^2; 
% MwT - Hastigheten massan vatten förändras (18)
MwPrim = @(xVwn) -Pw*xVwn*S;
% Vb - raketens förändringshastighet (19) - integraldelen 
VbPrim = @(xVwn, Pa) Sb*(Pa - Patm) - ((((Sb-S)^2)/(2*Sb))*(Pw*xVwn^2));
% Vwn - vattnets förändringshastighet (20) - integraldelen
VwnPrim = @(xMw) ((Sb^2-S^2)/xMw) + ((Sb - S)^2)/(Mb+Ma);
% Vpa - Lufttryckets förändringshastighet
VPaPrim = @(xVwn, Pa) ((-Ips*Pa^(Ips+1/Ips))/((V-Vw0)*Pa0^(1/Ips)))*xVwn*S;

% Sätt alla hastigheter till 0
MwT = 0; % Massans förändring acceleration kg/s^2
VbT = 0; % Raketens acceleration v/s^2
Vwn = 0; % Vattnets acceleration l/s^2
VPa = 0; % Tryckets förändring acceleration bar
Vb = 0; % Raketens hastighet v = m/s

while Mw > 0 % Kör loopen tills massan vatten i flaskan är noll
    
    %% Lite outputs
    
    % Vwn i delar
    Vwn1 = (((Sb^2)*(VPa-Patm))/S)*((1/Mw)+(1/(Mb+Ma)))
    Vwn2 = (((Vwn^2)*Pw)/(2*S))*VwnPrim(MwT)
    Vwn3 = (Sb/(S*(Mb+Ma)))*(Fd(Vb)+Fm+Wx)
    VwnTot = Vwn1 - Vwn2 - Vwn3
    
    
    % Vbt i delar
    Vbt1 = (1/(Mb+Ma))*(Fd(Vb)+Fm+Wx)
    Vbt2 = ((1/(Mb+Ma))*VbPrim(Vwn, VPa))
    VbtTot = Vbt1 + Vbt2
    
    % VPa i delar
    VPa1 = VPaPrim(Vwn, VPa)
    
    
    
    % MwT i delar
    MwT1 = MwPrim(Vwn)
    
    %% Skarpa differentialekvationer med eulers stegmetod Yk+1 = Yk + f(Yk) * h
    
    % Vwn Vattnets hastighet ut ur flaskan (20) 
    Vwn = Vwn + ((((Sb^2)*(VPa-Patm))/S)*((1/Mw)+(1/(Mb+Ma)))-(((Vwn^2)*Pw)/(2*S))*VwnPrim(Mw) - (Sb/(S*(Mb+Ma)))*(Fd(Vb)+Fm+Wx)) * h
    
    
    % VbT Raketens hastighet (19)
    VbT = VbT + ((1/(Mb+Ma))*(Fd(Vb)+Fm+Wx) + (1/(Mb+Ma))*VbPrim(Vwn, VPa)) * h
    
    % VpA Lufttryckets förändringshastighet
    VPa = VPa + VPaPrim(Vwn, VPa) * h
    
    % MwT Massan vattens förändringshastighet
    MwT = MwT + MwPrim(Vwn) * h;
    
    
    
    
    % Vattnets minskning för varje h
    Mw = Mw - MwT*h^2;
    disp("MwT")
    disp(MwT)
    
    Vb = VbT*time;
    
    % Raketens förflyttning i x-led för varje h, s = s + v*t
    x = x + Vb*time;
    disp("Mw: ")
    disp(Mw)
    disp("Vb: ")
    disp(Vb)
    
    xMatrix = [xMatrix, x];
    %yMatrix = [yMatrix, y];
    
    time = time + h
end


plot(xMatrix, Vb)
ylabel("Hastighet")
xlabel("Distance (m)")

end



