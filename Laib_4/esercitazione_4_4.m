% Esercitazione 4, esercizio 4
% Simone Canevarolo
% S269893
% 10/11/2023

clear all
close all
clc

%%

ll = 0.05; % m lunghezza
kk = 5; % W/m*K conducibilità termica
qv = 5e5; % W/m^3 generazione volumetrica
T0 = 300; % K temperatura imposta sulla facciata B

Tf = 273; % K temperatura
hf = 100; % W/(m^2*K) coefficiente di scambio termico

dx = 1e-4; % m passo
xx = (0:dx:ll)'; % vettore lunghezza
NN = length(xx); % scalare

%%

sub_diag = ones(NN,1);
main_diag = -2*ones(NN,1);
sup_diag = sub_diag;

Band = [[sub_diag(2:end);0],main_diag,[0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,NN,NN);

bb = -qv*(dx^2)/kk*ones(NN,1);

%%

% Robin

AA(1,1) = -1-hf*dx/kk;
AA(1,2) = 1;
bb(1) = -hf*dx*Tf/kk;

% Dirichlet

AA(end,end-1) = 0;
AA(end,end) = 1;
bb(end) = T0;

TT = AA\bb;

Tmax = max(TT)

%%

figure(1)
plot(xx,TT,'linewidth',2)
title('Andamento della temperatura nello spessore')
xlabel('spessore')
ylabel('Temperatura')

