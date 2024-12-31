% Esercitazione 8, esercizio 6
% Simone Canevarolo
% S269893
% 26/12/2023

% La prima parete di un reattore a fusione nucleare di tipo Tokamak è
% sottoposta a enormi carichi termici provenienti dal plasma, e deve
% essere refrigerate utilizzando acqua a 70 °C. Si approssimi la prima
% parete come una lastra composta di Cu-Be (spessore 2 cm,
% conducibilità termica 400 W/(m K)), e la si consideri da un lato
% irraggiata da un plasma alla temperatura di 2000 K e dall’altro
% refrigerata dall’acqua, con un coefficiente di scambio termico tra
% l’acqua e la parete di ~10 kW/(m2 K). Si calcoli il profilo di temperatura
% nella prima parete, valutando a stazionario il carico termico radiativo
% sulla parete e verificando che la superficie esposta al plasma non superi
% la temperatura limite di 550K.

clear all
close all
clc

%%

ll = 2e-2; % lunghezza, m
Th2o = 343; % temperatura dell'acqua, K
kk = 400; % condubibilità termica, W/(m*K)
Tpl = 2e3; % temperatura del plasma, K
hh = 1e4; % coefficiente di scambio termico acqua-parete, W/(m^2*K)
Tlim = 550; % temperatura limite parete, K
sigma = 5.67e-8; % Costante di Stefan-Boltzmann [W m^-2 K^-4]

dx = 1e-4;
xx = (0:dx:ll)';
NN = length(xx);

%%

    sub_diag = ones(NN,1);
    main_diag = -2*ones(NN,1);
    sup_diag = sub_diag;

    Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];

    AA = spdiags(Band,-1:1,NN,NN);

    bb = zeros(NN,1);

    AA(end,end-1) = kk/dx;
    AA(end,end) = hh-kk/dx;
    bb(end) = -hh*Th2o;

    %% Creazione della matrice



Tguess = Th2o*ones(NN,1);
ii = 0;
iimax = 1000;
toll = 0.001;
err = toll+1;

while err>toll && ii<iimax

    ii = ii+1;


    AA(1,1) = 
    AA(1,2) = 
    bb(1) = 

    TT = AA\bb;

    err = 

end


    TT = AA\bb;

figure(1)
plot(xx*1e2,TT)
title('Andamento della temperatura lungo la parete')
xlabel('Spessore [cm]')
ylabel('Temperatura [°C]')





% Condizioni al contorno
T(1) = 100; % Temperatura iniziale a un'estremità
T(Nx) = epsilon * sigma * (T_ambiente^4 - T(Nx)^4) * dx / (alpha * dx + epsilon * sigma * dx); 
% Condizione al contorno di irraggiamento all'altra estremità

% Iterazione nel tempo
for step = 1:num_steps
    % ... Algoritmo numerico per la diffusione del calore ...
    
    % Aggiornamento della temperatura con la condizione al contorno di irraggiamento
    T(Nx) = epsilon * sigma * (T_ambiente^4 - T(Nx)^4) * dt / (alpha + epsilon * sigma * dx);
end
