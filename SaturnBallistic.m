function [x] = SaturnBallistic

%% Konstanter
atmDens = 1.0248; % 1005 mBar = 1.0248 kg/m^3 [kg/m^3]
Cd = 0.29; % Magic
A = 0.0045; % Tvärsnittsarea flaskans kon [m^2]
Mb = 0.196; % Raketens tomma massa, 0.196 kg    [kg]  
g = 9.82; % Gravitations konstanten [m/s2]
L = 0.55; % Mängden vatten, 33% fylld flaska, 0,5 l
Bar = 7; 

% Starting conditions
Vinkel = 44; % degrees to horizontal
V = RaketensHastighetEnkel(L, Bar, Mb); % m/s
%V = 60;
x = 10; % m %% Adderar 10 m för luftens framdrivning
y = 0; % m
time = 0; % s

stepLength = 0.001; %Seconds
Vy = sind(Vinkel)*V;
Vx = cosd(Vinkel)*V;
xMatrix = [];
yMatrix = [];
tMatrix = [];

%VxPrim = @(Vx) (-atmDens .* dragCoefficient .* frontalArea .* Vx .* sqrt(Vx.^2 + Vy.^2)) ./ (2 .* emptyMass);

%VyPrim = @(Vy) -g - ((atmDens .* dragCoefficient .* frontalArea .* Vy .* sqrt(Vx.^2 + Vy.^2)) ./ (2 .* emptyMass));


while y >= 0 
    VxNew = Vx + ((-(atmDens*Cd*A)/(2*Mb)*Vx*sqrt(Vx^2+Vy^2)))*stepLength;
    VyNew = Vy + (-g-((atmDens*Cd*A)/(2*Mb)*Vy*sqrt(Vx^2+Vy^2)))*stepLength;
    Vx = VxNew;
    Vy = VyNew;
    
    x = x + Vx * stepLength;
    y = y + Vy * stepLength;
    
    xMatrix = [xMatrix, x];
    yMatrix = [yMatrix, y];
    tMatrix = [tMatrix, time];
    
    time = time + stepLength;
end

y
x
time

plot(xMatrix, yMatrix)
ylabel("Altitude (m)")
xlabel("Distance (m)")