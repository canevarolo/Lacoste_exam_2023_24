% Esercitazione 11, esercizio 1
% Simone Canevarolo
% S269893
% 13/01/2024

% Una grande fornace industriale è sostenuta da una colonna di mattoni, di lato 1 m × 1 m.
% Durante l’operazione in condizioni stazionarie, l’installazione del forno è tale per cui tre
% facce della colonna sono mantenute a 500 K, mentre la faccia rimanente è esposta a una
% corrente di aria a 300 K, con un coefficiente di scambio termico h = 10 W/m2K. Determinare
% la mappa di temperatura sulla sezione della colonna e la potenza dispersa nell’aria, per unità
% di lunghezza della colonna.

clear all
close all
clc

%%

Taria = 300; % temperatura, [K]
Timp = 500; % temperatura imposta, [K]
ll = 1; % lunghezza, [m]

ro = 1920; % densità, kg/m^3
cond = 0.72; % conducibilità termica, W/K/m
cp = 835; % calore specifico, J/kg/K
hh = 10; % coefficiente di scambio termico, W/m^2/K

dx = 5e-2; % m
dy = dx;

xx = (0:dx:ll)';
yy = xx;

Nx = length(xx);
Ny = length(yy);
Ntot = Nx*Ny;

[xmat,ymat] = meshgrid(xx,yy);

% Costruisco matrice AA

AA = sparse([],[],[],Ntot,Ntot,5*Ntot);

% creo il cornicione della matrice
for ii=2:Nx-1
    for jj=2:Ny-1

        kk=Nx*(jj-1)+ii; % dichiaro la variabile kk
        AA(kk,kk-Nx) = 1/dy^2;
        AA(kk,kk-1) = 1/dx^2;
        AA(kk,kk) = -2/dx^2-2/dy^2;
        AA(kk,kk+1) = 1/dx^2;
        AA(kk,kk+Nx) = 1/dy^2;

    end
end

bb = zeros(Ntot,1);

%% Imposto le condizioni al contorno
% Inserisco il lato con convezione (Robin) a Nord, nei restanti ho temperatura
% imposta (Dirichlet)

% Ovest, temperatura imposta
kk = 1:Nx;
AA(kk,:) = 0;
AA(kk,kk) = eye(length(kk));
bb(kk) = Timp;

% Est, temperatura imposta
kk = Nx*(Ny-1)+1:Ntot;
AA(kk,:) = 0;
AA(kk,kk) = eye(length(kk));
bb(kk) = Timp;

% Sud, temperatura imposta
kk = Nx:Nx:Ntot;
AA(kk,:) = 0;
AA(kk,kk) = eye(length(kk));
bb(kk) = Timp;

% Nord, convezione (Robin)
ii=1;
AA(Nx+1:Nx:Ntot-Nx,:)=0; % lato nord dal secondo punto al penultimo

for jj=2:Ny-1

    kk=ii+(jj-1)*Nx;
    AA(kk,kk-Nx)=cond/2/dy^2;
    AA(kk,kk)=-(cond/dx^2+cond/dy^2+hh/dx);
    AA(kk,kk+1)=cond/dx^2;
    AA(kk,kk+Nx)=cond/2/dy^2;
    bb(kk)=-hh/cond/dx*Taria;

end

TT=AA\bb;
TTT=reshape(TT,Nx,Ny);

surf(xmat',ymat',TTT,'FaceColor','interp')
colorbar
xlabel('x (m)')
ylabel('y (m)')
zlabel('T (K)')
title('Profilo di temperatura sulla sezione')
