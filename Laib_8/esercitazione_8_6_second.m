% Esercitazione 8, esercizio 6
% Simone Canevarolo
% S269893
% 12/02/2024

clear all
close all
clc

ll = 4; % lunghezza, m
spess = 2.9e-3; % spessore, m
dout = 33.7e-3; % diametro tubo, m
rout = dout/2; % raggio esterno, m
rin = rout-spess; % raggio interno, m
kkss = 13; % conducibilità termica acciaio, W/m*K

% ll >> dd,ss
% Uso l'ipotesi del tubo infinito, posso lavorare in una sola dimensione

TinHe = 20+273; % temperatura di ingresso elio, K
ToutHe = 100+273; % temperatura di uscita elio, K

TmedHe = (TinHe+ToutHe)/2; % temperatura media elio, K
% potendo considerare il tubo infinito, lavoro nel punto a Temperatura
% media che trovo come media della temperatura di ingresso

Tout = 100+273; % temperatura imposta esterno, K

Pr = 0.7; % numero di Prandl elio
Re = 1e5; % numero di Reynolds elio
pp = 5; % pressione elio, bar
% Considero 1 bar = 1e5 Pa

% kk = @(T) 2.682e-3*(1+1.123e-3.*pp).*T.^(0.71*(1-2e-4.*pp)); % conducibilità termica, W/m/K
kelio = 2.682e-3*(1+1.123e-3.*pp).*TmedHe.^(0.71*(1-2e-4.*pp)); % conducibilità termica elio, W/m/K

% T = 200:400;
% figure(2)
% plot(T,kk(T))

% Correlazione di Diltus-Boelter
% Numero di Nusselt
Nu = 0.023*Re^(0.8)*Pr^(0.7);
% Coefficiente di scambio termico
hh = Nu*kelio/dout;

dr = 1e-6;
rr = (rin:dr:rout)';
Nx= length(rr);

%%

sub_diag = 1-dr./(2.*rr);
main_diag = -2*ones(Nx,1);
sup_diag = 1+dr./(2.*rr);

Band = [[sub_diag(2:end);0],main_diag,[0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,Nx,Nx);

bb = zeros(Nx,1);

%%

AA(end,end) = 1;
AA(end,end-1) = 0;
bb(end) = Tout;

AA(1,1) = kkss/dr+hh;
AA(1,2) = -kkss/dr;
bb(1) = hh*TmedHe;

TT = AA\bb;

%%

figure(1)
plot(rr,TT-273)

Tmin = min(TT)
Tmax = max(TT)
