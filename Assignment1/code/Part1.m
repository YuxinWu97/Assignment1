% Electron Modelling
global C
C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665; %metres (32.1740 ft) per s
np = 500;
noe = 20;
r2 = randi(360,noe,1);
nx = 200;
ny = 100;
x = randi(200,noe,1);
y = randi(100,noe,1);

% Thermal Velocity
vth = sqrt((C.kb * 300)/(C.m_0 * 0.26));
vx = vth * cos(r2) ;
vy = vth * sin(r2);
% Mean free path
MFP = vth * 0.2 * 10^-12;

colrArray = rand(noe,1);


for n = 1:np
    temp = y >= ny ;
    temp1 = y < ny ;
    temp = temp * -1;
    th = temp + temp1;    
    temp2 = y <= 0;
    temp3 = y > 0;
    temp2 = temp2 * -1;
    tl = temp2 + temp3;
    vy = vy .* th;
    vy = vy .* tl;
    tempx1 = x <= 200;    
    x = x .* tempx1;
    tempx2 = x < -0.1;
    tempx2 = tempx2 * 200;
    tempxFinal = x + tempx2;
    x = tempxFinal;
    dx = vx * (1/200000);
    dy = vy * (1/200000);
    x = x + dx;
    y = y + dy;
    vsq = (vy).^2 + (vx).^2 ;
    average = mean(vsq);
    CT = (average *(0.26)* C.m_0)/(C.kb);
   
    figure(1)
    plot(n , CT,'.r')
    title("Temperature");
    xlabel("Simulation Number")
    ylabel("Temperature(K)")
    axis([0 np 0 500]);
    hold on    
    
    
    figure(2)
    scatter(x,y,3,colrArray);
    axis([0 200 0 100]);
    xlabel("x position (nm)")
    ylabel("y position (nm)")
    title("MFP:" + MFP);
    hold on;
    
    
    pause(0.01);  
end


