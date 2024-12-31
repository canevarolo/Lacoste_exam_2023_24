% Esercitazione 10, esercizio 2
% Simone Canevarolo
% S269893
% 16 aprile 2024

clear all
close all
clc

II = 50; % corrente, A
ll = 1; % lunghezza pezzo, m
roel = 1.75e-8; % resistività elettrica, Ohm*m
alfa = 11e-5; % diffusività, m^2*s
cond = 350; % conducibilità, W/m/K
dd = 4e-3; % diametro interno, m
Abase = pi*(dd/2)^2;
Vol = Abase*ll;

T0 = 300; % K

bordosx = load('Tleft.dat');
tempo = bordosx(:,1);
Tleft = bordosx(:,2);
tmax = tempo(end);

dx = 0.01; % lunghezza, m
xx = (0:dx:ll)';
Nx = length(xx);

dt = 1;
tt = (0:dt:tmax);
Nt = length(tt);

Tm = T0*ones(Nx,1);

kk = 1;

Tmid = 300*ones(length(tempo),1);

for ii = 2:Nt

    cp = cond/roel/alfa;
    qv = II^2*roel*ll/Abase/Vol;
    pos = ll/2/dx+1;

    sub_diag = -alfa*dt/dx^2*ones(Nx,1);
    main_diag = 1+2*alfa*dt/dx^2*ones(Nx,1);
    sup_diag = -alfa*dt/dx^2*ones(Nx,1);
    
    Band = [sub_diag,main_diag,sup_diag];
    
    AA = spdiags(Band,-1:1,Nx,Nx);
    
    % bb = Tm+qv*dt/roel/cp;
%     bb = Tm+qv*dt*alfa/cond;


    if tt(ii) == tempo(kk+1)

        kk = kk+1;

    end


    AA(1,1) = 1;
    AA(1,2) = 0;
    Tm(1) = Tleft(kk);

    AA(end,end) = 1;
    AA(end,end-1) = 0;
    Tm(end) = T0;

    % bb = Tm+qv*dt/roel/cp;
    bb = Tm+qv*dt*alfa/cond;

    TT = AA\bb;
    Tmid(ii) = TT(pos);

    Tm = TT;

end

figure(1)
plot(tt,Tmid-T0)

