atmDens = 1.2754; %kg/m^3
dragCoefficient = 0.345; %Magic
frontalArea = 0.01; %m^2
emptyMass = 0.1; %kg
g = 9.82; %m/s^2


stepLength = 0.001; %Seconds
x = 0;
Vx = 20;
y = 0;
Vy = 20;
time = 0;
xMatrix = [];
yMatrix = [];

VxPrim = @(Vx) (-atmDens .* dragCoefficient .* frontalArea .* Vx .* sqrt(Vx.^2 + Vy.^2)) ./ (2 .* emptyMass);

VyPrim = @(Vy) (-g - (atmDens .* dragCoefficient .* frontalArea .* Vy .* sqrt(Vx.^2 + Vy.^2)) ./ (2 .* emptyMass));


while y > 0 | time < 0.1
    Vx = Vx + VxPrim(Vx) * stepLength;
    Vy = Vy + VyPrim(Vy) * stepLength;
    
    x = x + Vx * stepLength;
    y = y + Vy * stepLength;
    
    xMatrix = [xMatrix, x];
    yMatrix = [yMatrix, y];
    
    time = time + stepLength;
end
y
x
time

plot(xMatrix, yMatrix)
ylabel("Altitude (m)")
xlabel("Distance (m)")