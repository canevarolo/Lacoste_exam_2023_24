% esercitazione 4, esercizio 5
% Simone Canevarolo
% S269893
% 11/11/2023

clear all
close all
clc

%%

Tw = 400; % K temperatura acqua
Din = 20e-3; % m diametro interno
Dout = 24e-3; % m diametro esterno

kk = 35e-3; % W/(m*K) conducibilità
hin = 100; % W/(m^2*K) coefficiente di scambio termico interno
hout = 10; % W/(m^2*K) coefficiente di scambio termico esterno
Ta = 300; % K temperatura aria

rin = Din/2; % m raggio interno
rout = Dout/2; % m raggio esterno

dr = 1e-5; % m
rr = (0:dr:(rout-rin))';

NN = length(rr);

%%

sub_diag = 1-dr/2*rr;
main_diag = -2*ones(NN,1);
sup_diag = 1+dr/2*rr;

Band = [[sub_diag(2:end);0],main_diag,[0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,NN,NN);

bb = zeros(NN,1);

%%

% Robin

AA(1,1) = -1-hin*dr/kk;
AA(1,2) = 1;
bb(1) = -Tw*hin*dr/kk;

% Robin

AA(end,end-1) = 1;
AA(end,end) = -1-hout*dr/kk;
bb(end) = -Ta*hout*dr/kk;

TT = AA\bb;

%%

figure(1)
plot(rr,TT,'linewidth',2)
title('Andamento temperatura nel raggio')
xlabel('Raggio [m]')
ylabel('Temperatura [K]')