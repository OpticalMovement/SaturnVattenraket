function [V, time] = RaketensHastighet
% Beräknar raketens hastighet under initiala fasen när vattnet pressas ut
% ur flaskan

%% Tibert tipsade om att lägga på en buffert i 
%% framförallt luftdriftsfasen på kanske 10


Patm = 1.2754; % kg/m^3
u = 0; % Friktionskonstanten från rampen
Angle = 45; % Startvinkel i grader
g = 9.82; % m/s^2
Mb = 0.107; % Raketens tomma massa, 0.1 kg
Ma = 0.001225*1.5; % Luftens massa, 1,225 g/l * 1,5 liter  
Mw = 1; % 1 liter vatten = 1.000 kg 
Pw = 997; % Vattnets densitet kg/m3
FlaskaRadie = 0.0881/2; % Flaskans ytterdiameter
Sb = 2*pi*FlaskaRadie; % Cross sectional area of the rocket, 0.0881 m
S = 0.0205; % Innerdiameter flaskans mynning
Cd = 0.345; % Dragcoefficient
Ips = 5/3; % För monatomic ideal gas (wikipedia 

Fm = u*(Mb + Ma + Mw)*g*cosd(Angle); % Friktion raket-ramp
Wx = (Mb + Ma + Mw)*g*sind(Angle); % Vik 
h = 0.01; %Steglängd
x = 0;
y = 0;
time = 0;
xMatrix = [];
yMatrix = [];

% Aktuellt luftmotstånd vid viss lufthastighet
Fd = @(Vb) 0.5*Patm*Cd*Sb*Vb^2; 
% MwT - Hastigheten massan vatten förändras (18)
MwPrim = @(Vwn) -Pw*Vwn*S;
% Vb - raketens förändringshastighet (19) - integraldelen 
VbPrim = @(Vwn) Sb(Pa - Patm) - ((((Sb-S)^2)/(2*Sb))*(Pw*Vwn^2));
% Vwn - vattnets förändringshastighet (20) - integraldelen
VwnPrim = @(Mw) ((Sb^2-S^2)/Mw) + ((Sb - S)^2)/(Mb+Ma);
% Vpa - Lufttryckets förändringshastighet
VPaPrim = @(Vwn) ((-Ips*Pa^(Ips+1/Ips))/((V-Vw0)*Pa0^(1/Ips)))*Vwn*S;

% Sätt alla hastigheter till 0
MwT = 0;
VbT = 0;
Vwn = 0;
VPa = 0;

while y > -0.1
    
    % Massan vattens förändringshastighet
    MwT = MwT + MwPrim(Vwn) * h;
    % Raketens hastighet (19)
    VbT = VbT - (1/(Mb+Ma))*((Fd+Fm+Wx) + (1/(Mb+Ma)))*VbPrim(Vwn) * h;
    % Vattnets hastighet ut ur flaskan (20)
    Vwn = Vwn + (((Sb^2)*(Pa-Patm))/S)*((1/Mw)+(1/(Mb+Ma))-((Vwn^2)*Pw)/(2*S))*VwnPrim(MwT) * h;
    % Lufttryckets förändringshastighet
    VPa = VPa + VwnPrim * h;
    
    
    x = x + Vb * h;
    
    xMatrix = [xMatrix, x];
    yMatrix = [yMatrix, y];
    
    time = time + h;
end

plot(xMatrix, yMatrix)
ylabel("Altitude (m)")
xlabel("Distance (m)")

end

