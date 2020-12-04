function [V, time] = RaketensHastighet
% Beräknar raketens hastighet under initiala fasen när vattnet pressas ut
% ur flaskan


Patm = 1.2754; % kg/m^3
u = 0; % Friktionskonstanten från rampen
Angle = 45; % Startvinkel i grader
g = 9.82; % m/s^2
Mb = 0.107; % Raketens tomma massa, 0.1 kg
Ma = 0.001225*1.5; % Luftens massa, 1,225 g/l * 1,5 liter  
Mw = 1; % 1 liter vatten = 1.000 kg 
FlaskaRadie = 0.0881/2; % Flaskans ytterdiameter
Sb = 2*pi*FlaskaRadie; % Cross sectional area of the rocket, 0.0881 m
S = 0.0205; % Innerdiameter flaskans mynning
Cd = 0.345; % Dragcoefficient
Fd = 0.5*Patm*Cd*Sb*Vb^2; % Luftmotstånd
Fm = u(Mb + Ma + Mw)*g*cosd(Angle); % Friktion raket-ramp
Wx = (Mb + Ma + Mw)*g*sind(Angle); % Vik 
h = 0.01; %Steglängd
x = 0;
y = 0;
time = 0;
xMatrix = [];
yMatrix = [];


% Vb - raketens hastighet (19) - integraldelen 
VbPrim = @(Vwn) Sb(Pa - Patm) - ((((Sb-S)^2)/(2*Sb))*(Pw*Vwn^2));
% MwT - Hastigheten massan vatten förändras (18)
MwPrim = @(Vwn) -Pw*Vwn*S;
% Vwn - vattnets hastighet (20) - integraldelen
VwnPrim = @(Mw) ((Sb^2-S^2)/Mw) + ((Sb - S)^2)/(Mb+Ma);

VPaPrim = @(Vwn) ((-Y*Pa^(Y+1/Y))/(V-Vw0)*Pa0^(1/Y))*Vwn*S;


while y > -0.1
    % Massan vattens förändringshastighet
    MwT = MwT + MwPrim(Vwn) * h;
    % Raketens hastighet (19)
    VbT = VbT -(1/(Mb+Ma))*(Fd+Fm+Wx)+(1/(Mb+Ma))*VbPrim(Vwn) * h;
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

