% Esame del 11 febbraio 2022
% Simone Canevarolo
% S269893
% 29/06/2024

clear all
close all
clc

diam = 10e-2; % m
raggio = diam/2; % m
lungh = 4; % m
lunghout = 1; % m
ss = 21e-3; % m
Taria = 283; % K
hout = 25; % W/m^2/K
hin = 500; % W/m^2/K

roaria = 1; % kg/m^3
ross = 140; % kg/m^3
kss = 0.3; % W/m/K
cpss = 400; % J/kg/K

dl = 1e-4; %m
ll = (0:dl:lungh)';
Nl = length(ll);

dr = 1e-4; 
rr = (0:dr:raggio)';
Nr = length(rr);

Areabase = pi*raggio^2; % m^2
As = 2*pi*raggio*lungh; % m^2
Vol = Areabase*lungh; % m^3

G0 = 70; % N*m^3/h
% GG è la portata
Tfumi = @(G) -1.5*(G+200)+273;
T0 = Tfumi(70);

dt = 1; % s
uu = G0/Areabase; % m/s

aa = kss*dt/dl^2/ross/cpss;
bb = hin*As*dt/Vol/ross/cpss;
cc = ross*cpss/dt;
dd = uu*dt/dl;
ee = hin*As*dt/Vol/ross/cpss;

subsub = zeros(2*Nl,1);
sub = zeros(2*Nl,1);
main = zeros(2*Nl,1);
sup = zeros(2*Nl,1);
supsup = zeros(2*Nl,1);

subsub(1:2:end) = 0;
subsub(2:2:end) = -dd;

sub(1:2:end) = -aa;
sub(2:2:end) = (1+bb+dd-ee);

main(1:2:end) = (1+2*aa-bb+ee);
main(2:2:end) = 0;

sup(1:2:end) = -aa;
sup(2:2:end) = 0;

supsup(1:2:end) = 0;
supsup(2:2:end) = 0;

Band = [subsub, sub, main, sup, supsup];

AA = spdiags(Band,-2:2, 2*Nl, 2*Nl);

bb = ones(2*Nl,1);

% Neumann adiabatico

AA(1,:) = 0;
AA(1,1) = 1;
AA(1,3) = -1;
bb(1) = 0;

% Dirichlet

AA(2,:) = 0;
AA(2,2) = 1;
bb(2) = T0; 

% Neumann adiabatico

AA(end-1,:) = 0;
AA(end-1,end-1) = 1;
AA(end-1,end-3) = -1;
bb(end-1) = 0;

% ???

AA(end,:) = 0;
AA(end,end) = uu/dl;
AA(end,end-2) = -uu/dl;
bb(end) = 0;


TT = AA\bb;

figure(1)
plot(ll,TT(1:2:end))
hold on
plot(ll,TT(2:2:end))

% RIVEDERE ERRORE DI MATRICE SINGOLARE

