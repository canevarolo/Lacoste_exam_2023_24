% Esercitazione 8, esercizio 6
% Simone Canevarolo
% S269893
% 04/01/2024

clear all
close all
clc

%%

dout = 33.7e-3; % diametro tubo [m]
rout = dout/2;
ll = 4; % lunghezza tubo [m]
spess = 2.9e-3; % spessore tubo [m]
rin = rout-spess;
kkss = 13; 

dr = 1e-6;
rr = (rin:dr:rout)';
NN = length(rr);

pp = 5; % pressione nei tubi [bar]

Test = 100+273; % temperatura esterna [K]
Tmedia = 60+273; % temperatura media [K]

Pr = 0.7; % numero di Prandl
Re = 1e5; % numero di Reynolds

kelio = 2.682e-3*(1+1.123e-3*pp)*Tmedia.^(0.71*(1-2e-4*pp));

hh = 0.023*Re^(0.8)*Pr^(0.7)*kelio/dout;

%%

sub_diag = 1-dr./(2.*rr);
main_diag = -2*ones(NN,1);
sup_diag = 1+dr./(2.*rr);

Band = [[sub_diag(2:end);0],main_diag,[0;sup_diag(1:end-1)]];

AA = spdiags(Band, -1:1, NN, NN);

bb = zeros(NN,1);

%%

AA(end,end) = 1;
AA(end,end-1) = 0;
bb(end) = Test;

AA(1,1) = kkss/dr+hh;
AA(1,2) = -kkss/dr;
bb(1) = hh*Tmedia;

TT = AA\bb;

figure(1)
plot(rr,TT-273)

Tmin = min(TT)-273
Tmax = max(TT)-273