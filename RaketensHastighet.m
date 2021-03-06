function [Vb, time] = RaketensHastighet
% Beräknar raketens hastighet under initiala fasen när vattnet pressas ut
% ur flaskan

%% Tibert tipsade om att lägga på en buffert i 
%% framförallt luftdriftsfasen på kanske 10
clc
close all
clear
format short
%% Konstanter
Patm = 1.0248; % 1005 mbar = 101300 N/m2
u = 0.1; % Friktionskonstanten från rampen
Angle = 45; % Startvinkel i grader
g = 9.82; % m/s^2
Mb = 0.133; % Raketens tomma massa, 0.107 kg
Ma = 0.001225*1.5; % Luftens massa, 1,225 g/l * 1,5 liter  
Pw = 997; % Vattnets densitet kg/m3
Pa0 = 7*Patm; % 2 Bar som en höftning, 20394 kg/m2
Vol = 0.001564; % Flaskans volym i kubikmeter
Vw0 = 0.5; % Vattnets volym, 33% fylld flaska, 0,5 l 
FlaskaRadie = 0.0881/2; % Flaskans ytterdiameter 0.0881/2
Sb = pi*FlaskaRadie^2; % Tvärsnittsarea flaskan pi*r^2
S = pi*(0.0205/2)^2; % Tvärsnittsarea flaskans mynning pi*r^2 
Cd = 0.36; % Dragcoefficient 0.345
Ips = 1.4; % För monatomic ideal gas (wikipedia) 

h = 0.01; %Steglängd
x = 0;
y = 0;
time = 0;
xMatrix = [];
yMatrix = [];

%% Integraler
% Aktuellt luftmotstånd vid viss lufthastighet
Fd = @(xVb) 0.5*Patm*Cd*Sb*xVb^2; 
% MwT - Hastigheten massan vatten förändras (18)
MwPrim = @(xVwn) -Pw*xVwn*S;
% Vb - raketens förändringshastighet (19) - integraldelen 
VbPrim = @(xVwn, Pa) Sb*(Pa - Patm) - ((((Sb-S)^2)/(2*Sb))*(Pw*xVwn^2));
% Vwn - vattnets förändringshastighet (20) - integraldelen
VwnPrim = @(xMw) ((Sb^2-S^2)/xMw) + ((Sb - S)^2)/(Mb+Ma);
% Vpa - Lufttryckets förändringshastighet
VPaPrim = @(xVwn, Pa) ((-Ips*Pa^((Ips+1)/Ips))/((Vol-Vw0)*Pa0^(1/Ips)))*xVwn*S;

Mw = 0.5; % 0,5 liter vatten = 0,5 kg 
Pa = Pa0; % Absoluta trycket i flaskan, börjar på Pa0 = 2 bar i vårt fall

% Sätt alla hastigheter till 0
MwT = 0; % Massans förändring acceleration kg/s^2
VbT = 0; % Raketens acceleration v/s^2
Vwn = 0; % Vattnets acceleration l/s^2
VPa = 0; % Tryckets förändring acceleration bar
Vb = 0; % Raketens hastighet v = m/s
yVb = [];
Fm = u*(Mb + Ma + Mw)*g*cosd(Angle); % Friktion raket-ramp
Wx = (Mb + Ma + Mw)*g*sind(Angle); % Vikt mot ramp 

while time < 0.07 % Kör loopen tills massan vatten i flaskan är noll
    
    % Hastigheter per steg, t ex Vwn = Wn * h, 630 m/s * 0.01 = 6,3 m 
    
    %% Lite test outputs
    
    % Vwn i delar
    Vwn1 = (((Sb^2)*(Pa-Patm))/S)*((1/Mw)+(1/(Mb+Ma)))
    Vwn2 = (((Vwn*h^2)*Pw)/(2*S))*VwnPrim(Mw)
    Vwn3 = (Sb/(S*(Mb+Ma)))*(Fd(Vb)+Fm+Wx)
    VwnSpeed = Vwn1 - Vwn2 - Vwn3
    VwnSum = VwnSpeed*h
    
    disp("Fd: ")
    disp(Fd(Vb))
    
    % Vbt i delar
    Vbt1 = (1/(Mb+Ma))*(Fd(Vb)+Fm+Wx)
    Vbt2 = ((1/(Mb+Ma))*VbPrim(Vwn, Pa))
    VbtTot = -Vbt1 + Vbt2
    
    % VPa i delar
    VPa1 = VPaPrim(Vwn*h, Pa)
    
    % MwT i delar
    MwT1 = MwPrim(Vwn*h)
    
    %% Skarpa differentialekvationer med eulers stegmetod Yk+1 = Yk + f(Yk) * h
    
    % Vwn Vattnets hastighet ut ur flaskan (20) 
    VwnNew = Vwn + ((((Sb^2)*(Pa-Patm))/S)*((1/Mw)+(1/(Mb+Ma)))-(((Vwn^2)*Pw)/(2*S))*VwnPrim(Mw) - (Sb/(S*(Mb+Ma)))*(Fd(Vb)+Fm+Wx))*h
    
    % VpA Lufttryckets förändringshastighet
    VPaNew = VPa + VPaPrim(Vwn, Pa)*h
    
    % VbT Raketens hastighet (19)
    VbTNew = VbT + (( -(1/(Mb+Ma))*(Fd(Vb)+Fm+Wx) + (1/(Mb+Ma))*VbPrim(Vwn, Pa)))*h
    
    % MwT Massan vattens förändringshastighet
    MwTNew = MwT + MwPrim(Vwn)*h
    

    % Vattnets minskning för varje h
    Mw = Mw - MwTNew*h;
    disp("MwT")
    disp(MwT)
    
    % Pa, Trycket i flaskans minskning för varje h
    Pa = Pa - VPaNew*h;
    disp("Pa")
    disp(Pa)

    Vb = VbT*time;
    
    Vwn = VwnNew;
    VPa = VPaNew;
    VbT = VbTNew;
    MwT = MwTNew;
    
    yVb(end+1) = Vb;
    
    time = time + h;

    disp("Vb: ")
    disp(Vb)
    disp("Time: ")
    disp(time)
    
    
end



%plot(time, yVb)
%ylabel("Flaskans Hastighet (Vb)")
%xlabel("Tid (s)")

end



