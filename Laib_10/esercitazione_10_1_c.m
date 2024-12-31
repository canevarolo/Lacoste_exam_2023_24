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

dtvett = [0.05 0.1]; % s
dxvett = [0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05];

for kk = 1:length(dtvett)

    dt = dtvett(kk);
    tt = (0:dt:tfin);
    Nt = length(tt);

    for jj = 1:length(dxvett)
        
        dx = dxvett(jj);
        xx = (0:dx:ll)';
        Nx = length(xx);
        
        Tanalitica = Tan(xx,tfin);
        
        Tm = Tini(xx);

        for ii = 2:Nt % importante: è Nt NON Nx !!!
        
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
loglog(dxvett,err,'linewidth',2)
hold on
grid on

end

xlabel('Precisione spaziale (m)')
ylabel('Errore relativo')
title('Errore relativo in spazio')