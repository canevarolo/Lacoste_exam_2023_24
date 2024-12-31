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

dtvett = [0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20]; % s
dxvett = [0.005 0.01]; % m




% tt = (0:dt:tfin);
% Nt = length(tt);
% graftt = [20,40,60,80,100];
% tinput = round(graftt/dt)+1;
% count = 1;

% figure(1)
% plot(xx,Tm,'LineWidth',2,'color','g');
% Tmax = max(Tm)
% hold on
    
for kk = 1:length(dxvett)

dx = dxvett(kk);
xx = (0:dx:ll)';
Nx = length(xx);

Tanalitica = Tan(xx,tfin);

    for jj = 1:length(dtvett)
    
        dt = dtvett(jj);
        tt = (0:dt:tfin);
        Nt = length(tt);
    
        Tm = Tini(xx);
    
        for ii = 2:Nt
        
            aa = alfa*dt/dx^2;
        
            sub_diag = -aa*ones(Nx,1);
            main_diag = (1+2*aa)*ones(Nx,1);
            sup_diag = -aa*ones(Nx,1);
        
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
            Tm = TT;
    
        end
    
    err(jj) = norm(TT-Tanalitica)/norm(Tanalitica-T0);
    
    end

    figure(1)
    loglog(dtvett,err,'linewidth',2)
    hold on
    grid on

    pause(0.1)

end

xlabel('Precisione temporale (s)')
ylabel('Errore relativo')
title('Errore relativo a t=100 s')