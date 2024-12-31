% Esercitazione 4, esercizio 3
% Simone Canevarolo
% S269893
% 27/10/2023

clear all
close all
clc

%% 

DD = 1e-2; % m diametro
LL = 4; % m lunghezza
II = 1e3; % A corrente
Tb = 77; % K temperatura
% estremità a Tb costanti, Dirichlet
hh = 500; % W/(m^2*K)
rocu = 1.75e-8; % omega*m resistività elettrica
kk = 350; % W/(m*K) coonducibilità termica

As = pi*(DD/2).^2; % m^2 area sezione (non SCAMBIO)
Alat = 2*pi*(DD/2)*LL; % m^2 area laterale
Vol = As*LL; % m^3 volume

dx = 1e-3; % m delta lunghezza
xx = (0:dx:LL)';

NN = length(xx);

sub_diag = ones(NN,1);
main_diag = (-2-hh*Alat/Vol)*ones(NN,1);
sup_diag = sub_diag;

Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,NN,NN);

bb = (-(rocu*II^2)/(As^2)-hh*As/Vol*Tb)*(dx^2/kk)*ones(NN,1);

%%

% Dirichlet
AA(1,1) = 1;
AA(1,2) = 0;
bb(1) = Tb;

% Dirichlet
AA(end,end-1) = 0;
AA(end,end) = 1;
bb(end) = Tb;

TT = AA\bb;

% non serve ma verifico
% Tmax=max(TT)

%%

figure(1)
title('Andamento della temperatura nello spessore')
plot(xx,TT,'LineWidth',2)