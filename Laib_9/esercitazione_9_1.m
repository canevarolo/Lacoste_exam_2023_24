% Esercitazione 9, esercizio 1
% Simone Canevarolo
% S269893
% 13/01/2024

clear all
close all
clc

%%

Tin = 450+273; % temperatura iniziale del pezzo, [K]
Tw = 15+273; % temperatura acqua, [K]
hh = 100; % coefficiente di scambio termico, [W/(m^2*K)]

% misure geometriche
Dmagg = 6e-1; % diametro maggiore, [m]
Dmin = 4e-1; % diametro minore, [m]
spessmagg = 3e-1; % spessore pezzo maggiore, [m]
spessmin = 1e-1; % spessore pezzo cilindro centrale, [m]
spesslam = 1e-2; % spessore lamina, [m]

% Ricavo per comodità le seguenti lunghezze
rmagg = Dmagg/2; % raggio maggiore, [m]
rmin = Dmin/2; % raggio minore, [m]

% Trascrivo i dati del documento
kkall = 237; % conduttività termica alluminio, [W/(m*K)]
kkacc = 14.2; % conduttività termica acciaio, [W/(m*K)]
kkh2o = 0.6; % conduttività termica acqua, [W/(m*K)]

roall = 2700; % densità alluminio, [kg/m^3]
roacc = 7978; % densità acciaio, [kg/m^3]
roh2o = 1e3; % densità acqua, [kg/m^3]

ccall = 897; % calore specifico alluminio, [J/(kg*K)]
ccacc = 480; % calore specifico acciaio, [J/(kg*K)]
cch2o = 4186; % calore specifico acqua, [J/(kg*K)]

%%

% Al fine di considerare il corpo come un punto, trovo le misure
% equivalenti e i volumi 

Vacc = 2*(pi*spesslam*(rmagg^2-rmin^2)); % volume acciaio, [m^3]
Vall = pi*(spessmagg*rmin^2+spessmin*(rmagg^2-rmin^2)); % volume alluminio, [m^3]
Vtot = Vacc+Vall;

Atot = pi*(rmagg)^2+spessmagg*2*pi*rmagg;

keq = (Vacc*kkacc+Vall*kkall)/Vtot;
roeq = (Vacc*roacc+Vall*roall)/Vtot;
cpeq = (ccall*roall*Vall+ccacc*roacc*Vacc)/(roeq*Vtot);

% Verifica del numero di Biot (<0.01)

dc = Vtot/Atot;
Bi = hh*dc/keq

tau = roeq*cpeq*Vtot/hh/Atot;

%%

toll = 1e-4;
err = toll+1;
ii = 1;
maxii = 1000;

dt = 1; % tempo, [s]
tmax = 1000; % tempo, [s]
tt = (0:dt:tmax);
ttlast = 0; % s

TT = Tin;
Tguess = Tin;

while err>toll && ii<maxii
     
    ii = ii+1;

    TT(ii) = Tguess+(Tw-Tguess)*dt/tau; 

    tt(ii) = tt(ii-1)+dt;

    err = abs(TT(ii)-Tguess)/abs(TT(ii));

end

figure(1)
plot(tt,TT)

