% Esercitazione 9, esercizio 2
% Simone Canevarolo
% S269893
% 23 gennaio 2024

clear all
close all
clc

%%

Tti = 600+273; % K
Toil = 15+273; % K
dd = 10e-2; % m
ll = 20e-2; % m
raggio = dd/2; % m

hh = 100; % coefficiente di scambio termico tra barra e olio, W/m^2/K
kti = 21.9; % W/m/K
koil = 0.1; % W/m/K
roti = 4500; % kg/m^3
rooil = 900: % kg/m^3
cpti = 520; % J/kg/K
cpoil = 7200; % J/kg/K

Atot = 2*pi*(raggio)^2+2*pi*raggio+ll; % area totale, m^2
VV = pi*(raggio)^2*ll; % volume, m^3

% Verifico il numero di Biot (<0.01)
dc = Atot/VV;
Bi = hh*dc/kk

xx = (0:dl:ll)';
NN = length(xx);

%%

iter = 1;
maxii = 1000;
toll = 1e-6;
err = toll+1;

dt = 1; % tempo, s
tt = 0; % s

% Avendo coefficienti costanti, calcolo la matrice dei coefficienti al di
% fuori del ciclo for

    sub_diag = 
    main_diag = 
    sup_diag = 

    Band = 

    AA = 

while err>toll && iter<maxii





end
