% Esercitazione 8, esercizio 3 - Newton
% Simone Canevarolo
% S269893
% 1/02/2024

clear all
close all
clc


Le = 1; % lunghezza pezzo, m
ll = Le/2; % lunghezza sezione pezzo per simmetria, m
Tco2 = 300; % temperatura, °C
hh = 50; % coefficiente di scambio termico, W/m^2/K

dx = 1e-3;
xx = (-ll/2:dx:0)';
Nx = length(xx);

qv = @(x) 8e3*cos(pi*x./Le); % generazione interna di calore, W/m^3
figure(1)
plot(xx,qv(xx))

kk = @(T) 100./(11.8+0.0238*T)+8.775e-11*T.^3; % conducibilità termica, W/m/K

%% Risoluzione punto fisso (non conservativo, caso 1)

ii = 1;
maxii = 1000;
toll = 1e-5;
err = toll+1;

Tguess = Tco2*ones(Nx,1);

while err>toll && ii<maxii

    ii = ii+1;

sub_diag = ones(Nx,1);
main_diag = -2*ones(Nx,1);
sup_diag = ones(Nx,1);

Band = [[sub_diag(2:end);0],main_diag,[0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,Nx,Nx);

bb = - qv(xx)*dx^2./kk(Tguess);

% Condizioni al contorno

AA(1,1) = hh+kk(Tguess(1))./dx;
AA(1,2) = -kk(Tguess(1))/dx;
bb(1) = hh*Tco2;

AA(end,end-1) = 1;
AA(end,end) = -1;
bb(end) = 0;

TT = AA\bb;

Tguess = TT;

end

figure(2)
plot(xx,TT)


%% Risoluzione punto fisso (conservativo, caso 5)

ii = 1;
maxii = 10000;
toll = 1e-5;
err = toll+1;

Tp = Tco2*ones(Nx,1);

while err>toll && ii<maxii
    ii = ii+1;

    Tim1 = Tguess(1:end-2);
    Ti = Tguess(2:end-1);
    Tip1 = Tguess(3:end);

sub_diag = -kk(Tip1)/4+kk(Tim1)./4+kk(Ti);
main_diag = -2*Tguess;
sup_diag = kk(Tip1)/4-kk(Tim1)/4+kk(Ti);

Band = [[sub_diag;0;0],main_diag,[0;0;sup_diag]];

AA = spdiags(Band,-1:1,Nx,Nx);

bb = - qv(xx)*dx^2;

% Condizioni al contorno

AA(1,1) = hh+kk(Tguess(1))./dx;
AA(1,2) = -kk(Tguess(1))/dx;
bb(1) = hh*Tco2;

AA(end,end-1) = -1;
AA(end,end) = 1;
bb(end) = 0;

TT = AA\bb;
err=norm(TT-Tguess)/norm(Tguess-Tco2);
Tguess = TT;

end

figure(3)
plot(xx,TT)

% Conservazione energia

 qin = trapz(xx,qv(xx)) %W/m2
 qout = hh*(TT(1)-Tco2)  %W/m2

 qin-qout

  if abs(qin-qout)>1e-2 %attenzione a questo numero, dipende dalla griglia!
     disp('Non conservo l''energia!')
 end

 % Non conservo l'energia nel mio caso


 %% NEWTON

ii = 1;
maxii = 10000;
toll = 1e-5;
err = toll+1;

qvec = qv(xx);
Tguess = zeros(Nx,1);

TTT = Told;

while err>toll && ii<maxii

    ii = ii+1;

    % Creo vettore funzione (colonna)
    fun = zeros(length(TTT),1);

    



end







