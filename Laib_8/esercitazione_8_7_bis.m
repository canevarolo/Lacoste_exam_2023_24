% Esercitazione 8, esercizio 6
% Simone Canevarolo
% S269893
% 09/01/2024

clear all
close all
clc

%%

ll = 2e-2; % spessore [m]
kk = 0.10; % conducibilità termica [W/(m*K)]
Tamb = 25+273; % temperatura ambiente [K]
qf = 75; % flusso termico [W/m^2]
sigma = 5.67e-8; % costante di Boltzmann

dx = 1e-5; % [m]
xx = (0:dx:ll)';
NN = length(xx);

%%

Tguess = Tamb*ones(NN,1);
ii = 0;
iimax = 1e3;
toll = 1e-6;
err = toll + 1;

while err>toll && ii<iimax

ii = ii+1;

sub_diag = ones(NN,1);
main_diag = -2*ones(NN,1);
sup_diag = sub_diag;

Band = [[sub_diag(2:end);0],main_diag,[0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,NN,NN);

bb = zeros(NN,1);

% Dirichlet
AA(1,1) = -1;
AA(1,2) = 1;
bb(1) = -qf*dx/kk;

hirr = sigma*(Tguess(end)^2+Tamb^2)*(Tguess(end)+Tamb);

% irraggiamento a parete dx con punto fisso
AA(end,end-1) = -1;
AA(end,end) = 1+dx/kk*hirr;
bb(end) = Tamb*dx/kk*hirr;


TT = AA\bb;

err = norm(TT-Tguess)/norm(TT-Tamb);

Tguess = TT;

end

figure(1)
plot(xx*1e2,TT-273)
title('Andamento della temperatura lungo la parete')
xlabel('Spessore [cm]')
ylabel('Temperatura [°C]')
grid on