% Esercitazione 8, esercizio slides con iterazioni di punto 
% Simone Canevarolo
% S269893
% 11/12/2023

clear all
close all
clc

%%

lungh = 1; % lunghezza [m]
Tco2 = 300; % temperatura [°C]
hh = 50; % coefficiente di scambio termico [W/(m^2*K)]
ll = 0.5*lungh;

qv = @(x) 8e3.*cos(pi.*x./ll); % generazione volumetrica di calore [W/m^3]
kk = @(T) 100./(11.8+0.0238.*T)+8.775e-11.*T.^3; % conducibilità termica [W/(m*K)]

dx = 1e-4;
xx = (0:dx:ll)';
NN = length(xx);

%%

Tguess = Tco2*ones(NN,1);
toll = 1e-4;
err = toll+1;
ii = 1;
imax = 1e3;

while err>toll && ii<imax

    ii = ii+1;

sub_diag = [kk(Tguess(2:end-1))-0.25*(kk(Tguess(1:end-2)))+0.25*kk(Tguess(3:end));0;0];
main_diag = [-2*kk(Tguess)];
sup_diag = [0;0;kk(Tguess(2:end-1))+0.25*(kk(Tguess(1:end-2)))-0.25*kk(Tguess(3:end))];

Band = [sub_diag, main_diag, sup_diag];

AA = spdiags(Band,-1:1,NN,NN);

bb = -qv(xx).*dx^2;

% Dirichlet

AA(1,1) = -1;
AA(1,2) = 1;
bb(1) = 0;

% Robin

AA(end,end-1) = +kk(Tguess(end-1))./dx;
AA(end,end) = -kk(Tguess(end))./dx-hh;
bb(end) = -hh*Tco2;

TT = AA\bb;

err = norm(TT-Tguess)/norm(TT-Tco2);

Tguess = TT;

end

figure(1)
plot(xx*1e3,TT,'linewidth',2)
xlabel('Spessore [mm]')
ylabel('Temperatura [°C]')
title('Profilo della temperatura dal centro del pezzo verso esterno')