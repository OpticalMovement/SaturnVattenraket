clear all, close all, clc

% Constants
atmDens = 1.2754; %kg/m^3
dragCoefficient = 0.345; %Magic
frontalArea = 0.005; %m^2
emptyMass = 0.15; %kg
g = 9.82; %m/s^2

% Starting conditions
angle = 45; % degrees to horizontal
V = 60; % m/s
x = 0; % m
y = 0; % m
time = 0; % s

stepLength = 0.001; %Seconds
Vy = sind(angle) .* V;
Vx = cosd(angle) .* V;
xMatrix = [];
yMatrix = [];
tMatrix = [];

%VxPrim = @(Vx) (-atmDens .* dragCoefficient .* frontalArea .* Vx .* sqrt(Vx.^2 + Vy.^2)) ./ (2 .* emptyMass);

%VyPrim = @(Vy) -g - ((atmDens .* dragCoefficient .* frontalArea .* Vy .* sqrt(Vx.^2 + Vy.^2)) ./ (2 .* emptyMass));


while y > 0 | time < 0.1
    VxNew = Vx + ((-atmDens .* dragCoefficient .* frontalArea .* Vx .* sqrt(Vx.^2 + Vy.^2)) ./ (2 .* emptyMass)) * stepLength;
    VyNew = Vy + (-g - ((atmDens .* dragCoefficient .* frontalArea .* Vy .* sqrt(Vx.^2 + Vy.^2)) ./ (2 .* emptyMass))) * stepLength;
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