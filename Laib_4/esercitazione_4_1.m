% In una parete infinita di spessore δ = 50 mm (v. Fig. 1) e conducibilità termica
% k = 5 W/(m K) è presente una generazione volumetrica di calore q’’’ = 500 kW/m3.
% La superficie A è adiabatica, mentre la superficie B è mantenuta ad una temperatura costante
% T0 = 300 K. Assumendo che il problema sia 1D lungo x, calcolare analiticamente e
% graficare la distribuzione di temperatura all’interno della parete. Ricalcolare quindi la
% distribuzione di temperatura utilizzando il metodo delle differenze finite (DF) e
% confrontare la soluzione numerica con quella analitica, usando stili di linea diversi,
% spiegati in una legenda. Verificare, infine, la conservazione dell’energia. 

% Esercitazione 4, esercizio 1
% Simone Canevarolo
% S269893
% 26/10/2023


clear all
close all
clc

%%

ll = 0.05; % m lunghezza
kk = 5; % W/m*K conducibilità termica
qv = 5e5; % W/m^3 generazione volumetrica
T0 = 300; % K temperatura imposta sulla facciata B

dx = 1e-4; % m passo
xx = (0:dx:ll)'; % vettore lunghezza
NN = length(xx); % scalare

%%

sub_diag = ones(NN,1);
main_diag = -2*ones(NN,1);
sup_diag = sub_diag;

Band = [[sup_diag(2:end);0], main_diag, [0;sub_diag(1:end-1)]];

AA = spdiags(Band,-1:1,NN,NN);

bb = -qv*dx^2/kk*ones(NN,1);

%% 

% Adiabatica, uso Neumann
AA(1,1) = -1;
AA(1,2) = 1;
bb(1) = 0;

% Temperatura imposta, uso Dirichlet
AA(end,end-1) = 0;
AA(end,end) = 1;
bb(end) = T0;

TT = AA\bb;

%%

figure(1)
plot(xx*1e2,TT,'LineWidth',2)
xlabel('Lunghezza [cm]')
ylabel('Temperatura [K]')
title('Andamento nello spessore del pezzo')

