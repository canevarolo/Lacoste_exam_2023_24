% Esercitazione 8, esercizio 1
% Simone Canevarolo
% S269893
% 8/12/2023

clear 
close all
clc

%% funzione

NN = 1000; % numero molecole
TT = 300; % temperatura [K]
pp = 3.5e7; % pressione [Pa]
kk = 1.3806503e-23; % costante di Boltzmann, J/K
% per CO2
aa = 0.401; % Pa*m^6
bb = 42.7e-6; % m^3

vdw = @(V) (pp+aa*(NN./V).^2).*(V-NN*bb)-kk*NN*TT;

%% plot funzione

figure
fplot(vdw,[0.01 0.1])
hold on
plot([0.01 0.1],[0 0], '--k')
xlabel('Volume (m^3)')
ylabel('vdw function')

ylim([-1e7 1e7])

%% Punto iniziale

% Uso valore vicino allo zero

V0 = 0.04;

%% Trovo f0

V1 = fzero(vdw,V0);

% fzero con opzioni

tic
options=optimset('MaxIter',1e5,'TolX',1e-12); 
[V2,res2,~,output]=fzero(vdw,V0,options);
time2=toc;

% output risultati

iter2=output.iterations;
met=output.algorithm;

% post processing
fprintf('\nFZERO (funzione interna MATLAB):\n')
fprintf('Volume: %3.7f m^3\n',V2)
fprintf('iterazioni: %i , residuo: %.3e,  , tempo richiesto: %.4e e utilizzando il metodo "%s"\n',iter2,norm(res2),time2,met);
