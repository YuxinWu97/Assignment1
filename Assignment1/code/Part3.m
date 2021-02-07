% Enhancements
clear
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

np = 300;
noe = 100;
r2 = randi(360,noe,1); 

colorArray = rand(noe,1);

nx = 200;
ny = 100;
x =  randi(200,noe,1);
y =  randi(100,noe,1);
vth = sqrt((C.kb * 300)/(C.m_0 * 0.26));
vx = vth * cos(r2) ;
vy = vth * sin(r2);

for pos = 1: noe
    xpos = x(pos);
    if (xpos < 120 && xpos > 80)
        if (y(pos) < 40)
            xpos = xpos + 50;
            x(pos) = xpos;
        elseif(y(pos) > 60)
            xpos = xpos - 50;
            x(pos) = xpos;
        else
        end
    end
end

MFP = vth * 0.2 * 10^-12;

pScat = 1 - exp((-3 * 10^-16)/(0.2 * 10^-12));

tMatrix = zeros(noe);

for t = 1:np
    
    
    vxc = vx;
    vyc = vy;
    [n,m] = size(vx);
    [n1,m1] = size(vy);
    idx = randperm(n);
    rdomx = vx;
    rdomx(idx,1)= vx (:,1) ;
    idy = randperm(n1);
    rdomy = vy;
    rdomy(idy,1) = vy(:,1);
    %Scattered
    rScatter= rand(noe,1);
    tempScatter = rScatter < pScat;
    rdomx = tempScatter .* rdomx; 
    rdomy = tempScatter .* rdomy;   
    %Not scattered
    notScatter = rScatter >= pScat;
    vx = vx .* notScatter;
    vy = vy .* notScatter; 
    vx = vx + rdomx;
    vy = vy + rdomy;
    xc = x; 
    yc = y ;
    temp = y >= ny;
    temp1 = y < ny; 
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
    tLR1s = ( x > 80 & x < 120) & y < 40;
    tLR0s = tLR1s == 0;
    tLR1s = -1 * tLR1s;
    f = tLR1s + tLR0s;
    vx = vx .* f;
    tLR1s = ( x > 80 & x < 120) & (y < 41 & y >= 40);
    tLR0s = tLR1s == 0;
    tLR1s = -1 * tLR1s;
    f = tLR1s + tLR0s;
    vy = vy .* f;
    tUR1s = ( x > 80 & x < 120) & y > 60;
    tUR0s = tUR1s == 0;
    tUR1s = -1 * tUR1s;   
    f = tUR1s + tUR0s;
    vx = vx .* f;
    tUR1s = ( x > 80 & x < 120) & (y >59 & y < 60);
    tUR0s = tUR1s == 0;
    tUR1s = -1 * tUR1s;
    vy = vy .* f;
    dx = vx * (1/200000);
    dy = vy * (1/200000);
    x = x + dx;
    y = y + dy;
    vsq = (vy).^2 + (vx).^2 ;
    average = mean(vsq);
    tMatrix = ((vsq * 0.26 *  C.m_0)/C.kb);
    figure (3);
    CT = (average *(0.26)* C.m_0)/(C.kb);
    plot(t , CT,'.r');
    title('Temperature plot')
    xlabel('time');
    ylabel('temperature(K)');
    axis([0 300 200 400]);
    hold on;  
    figure(4); 
    scatter (x, y , 3 ,colorArray);
    axis([0 200 0 100]);
    rectangle('Position',[80 0 40 40]);
    rectangle('Position',[80 60 40 40]);
    xlabel("x");
    ylabel("y");
    hold on;
    title ("Semiconductor Temperature:" + CT);
    %pause(0.001)
    figure (3);
    hold on;
end

scatter(x,y,'r.'); 
hold on;

Elecpos = [x,y];
D = hist3(Elecpos(:,1:2),'Nbins',[20,10]);
figure (6);
surf(D);
title('Electron Density plot');
shading interp;