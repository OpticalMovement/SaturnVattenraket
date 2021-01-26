clear all, close all, clc
r = 0.01;      %radie av utlåpshålet  [m]
cd = 1;     %Discharge coefficient
p = 998;    %Densitet   [kg/m^3]
pf = 101300; %Air pressure   [N/m^2]
pi = 3*pf;   %Atmospheric Pressure   [N/m^2]
V = 0.001;  %volymen av vatten i raketen   [L]
mh20 = p*V;  %massan av vattent        [kg] 
ma = 0.100;      %Massan av raketen    [kg]   
g = 9.82;    %gravitations konstanten
mave = (ma + (mh20/2));   %mass of rocket   [kg]
o = 45;      %vinkeln         [grader]
Vml = V*1000;
Dc = 0.3; %drag coefficient


deltP = pi-pf;  %skillnad mellan tryken
A = pi*r^2;     %arean av utlåpshålet
m = A*cd*(2*p*deltP)^0.5;       %mass flow rate
v = m/(p*A);         %exit velocity
ft = m*v;           %thrust
f = ft - (mave*g);     %medel thrust
a = f/mave;         %acceleration
t = mh20/m;         %tiden då allt vatten är ute
vbottle = a*t;     
R = (vbottle^2 * sind(2*o))/g  %distansen den flyger utan luft 
D = 1-Dc;   %drag factor
Rd = R*D;  %drag factor 

Rf = Rd - (0.00015 * (Vml - 700)^2) + 75  %distansen den flyger med en drag factor

