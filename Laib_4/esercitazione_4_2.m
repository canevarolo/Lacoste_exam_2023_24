

% Esercitazione 4, esercizio 2
% Simone Canevarolo
% S269893
% 27/10/2023

clear all
close all
clc

%%

kk = 0.15; % W/(K*m) conducibilità termica
Din = 20e-3; % m diametro interno
Dout = 50e-3; % m diametro esterno
rin=Din/2;
rout=Dout/2;
Tw = 15+273.15; % K temperatura
hh = 1e3; % W/(m^2*K) coefficiente di scambio termico

qq = 2e3; % W/m^2 flusso termico uniforme

dr = 1e-4;
rr = (rin:dr:rout)';

NN = length(rr);

%%

sub_diag = 1-dr./(2.*rr);
main_diag = -2*ones(NN,1);
sup_diag = 1+dr./(2.*rr);

Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,NN,NN);

bb = zeros(NN,1);

%% 

% Robin
AA(1,1) = -1+hh*dr/kk;
AA(1,2) =-1;
bb(1) = -Tw*dr*hh/kk;

% Neumann
AA(end,end-1) = -1;
AA(end,end) = 1;
bb(end) = qq*dr/kk;

TT = AA\bb;
% Tmax=max(TT)

%% 

figure(1)
plot((rr)*1e2-1,TT,'LineWidth',2)
xlabel('Lunghezza [cm]')
ylabel('Temperatura [K]')
title('Andamento nello spessore del pezzo')

%% Conservazione dell'energia

flow_in = qv*;
flow_out = abs(kk/dr*(TT(end-1)-TT(end)));
err1 = abs(flow_out - flow_in)/abs(flow_in);

fprintf('Relative mistake is %.5f with method 1\n',err1)

flow_out_2 = abs(kk/dx*(TT(end-1)-TT(end))+qv*dx/2);
err2 = abs(flow_out_2 - flow_in)/abs(flow_in);

fprintf('Relative mistake is %.5f with method 2\n',err2)