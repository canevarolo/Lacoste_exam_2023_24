% Esercitazione 10, esercizio 1
% Simone Canevarolo
% S269893
% 15 aprile 2024

clear all
close all
clc

ll = 10e-2; % lunghezza, m
alfa = 0.1e-4; % diffusività, m^2/s
T0 = 300; % temperatura, K

Tini = @(x) T0+50*sin(pi*x/ll);

Tan = @(x,t) T0+50*sin(pi*x/ll).*exp(-(pi/ll)^2*alfa*t);

tfin = 100; % tempo finale, s

dt = 0.1; % s
dx = 0.1e-2; % m

xx = (0:dx:ll)';
Nx = length(xx);

ii = 1;
Tm = Tini(xx);

tt = (0:dt:tfin);
Nt = length(tt);
graftt = [20,40,60,80,100];
tinput = round(graftt/dt)+1;
count = 1;

figure(1)
plot(xx,Tm,'LineWidth',2,'color','g');
Tmax = max(Tm)
hold on

for ii = 2:Nt

    aa = alfa*dt/dx^2;

    sub_diag = -aa*ones(Nx,1);
    main_diag = (1+2*aa)*ones(Nx,1);
    sup_diag = -aa*ones(Nx,1);

 %    Band = [[sub_diag;0], main_diag, [0;sup_diag]];
    Band = [sub_diag, main_diag, sup_diag];

    AA = spdiags(Band,-1:1,Nx,Nx);

    bb = Tm;

    AA(1,1) = 1;
    AA(1,2) = 0;
    bb(1) = T0;

    AA(end,end-1) = 0;
    AA(end,end) = 1;
    bb(end) = T0;

    TT = AA\bb;

if ii == tinput(count)

    plot(xx,Tm)
    hold on
    Tanalitica = Tan(xx,tt(ii));
    plot(xx,Tanalitica,'ob')
    hold on
    count = count+1;

end

    Tm = TT;

end

ylabel('Temperatura (K)')
xlabel('Lunghezza (m)')
title('Evoluzione della temperatura in funzione del tempo')
