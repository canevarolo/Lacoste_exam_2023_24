% Esercitazione 8, esercizio slides con iterazioni di punto fisso
% Simone Canevarolo
% S269893
% 26/12/2023

clear all
close all
clc

%%

% lungh = 1; % lunghezza [m]
% Tco2 = 300; % temperatura [°C]
% hh = 50; % coefficiente di scambio termico [W/(m^2*K)]
% ll = 0.5*lungh;
% 
% qv = @(x) 8e3.*cos(pi.*x./ll); % generazione volumetrica di calore [W/m^3]
% kk = @(T) 100./(11.8+0.0238.*T)+8.775e-11.*T.^3; % conducibilità termica [W/(m*K)]
% 
% dx = 1e-4;
% xx = (0:dx:ll)';
% NN = length(xx);

%%

function temperatura_barretta
    % Parametri del problema
    D = 0.2; % diametro in metri
    L = 1; % lunghezza in metri
    TCO2 = 300; % temperatura del CO2 in gradi Celsius
    h = 50; % coefficiente di scambio termico in W/m^2/K
    qv_max = 80000; % potenza volumetrica massima in W/m^3
    rho = 10000; % densità in g/cm^3
    cp = 300; % calore specifico in J/kg/K

    % Discretizzazione spaziale
    Nz = 100; % numero di nodi lungo z
    dz = L / Nz; % passo di discretizzazione
    z = linspace(-L/2, L/2, Nz); % coordinate z

    % Inizializzazione della temperatura
    T = TCO2 * ones(1, Nz); % inizializza la temperatura a quella del CO2

    % Metodo di Newton
    max_iter = 100;
    tol = 1e-6;

    for iter = 1:max_iter
        % Calcolo della matrice Jacobiana e del vettore dei residui
        [J, R] = calcola_jacobiana_residui(T, dz, h, qv_max, z);

        % Risoluzione del sistema lineare
        delta_T = J \ (-R');

        % Aggiornamento della temperatura
        T = T + delta_T';

        % Verifica della convergenza
        if max(abs(delta_T)) < tol
            break;
        end
    end

    % Plot della distribuzione di temperatura
    figure;
    plot(z, T);
    title('Distribuzione di temperatura nella barretta di combustibile');
    xlabel('Coordinate assiali (m)');
    ylabel('Temperatura (°C)');
end

function [J, R] = calcola_jacobiana_residui(T, dz, h, qv_max, z)
    % Calcola la matrice Jacobiana e il vettore dei residui

    Nz = length(T);
    J = zeros(Nz, Nz);
    R = zeros(1, Nz);

    for i = 2:Nz-1
        % Calcola la conducibilità termica al nodo i
        k = @(T) 100 / (11.8 + 0.0238*T) + 8.775e-11 * T^3;

        % Termini centrali per la derivata seconda
        alpha = k(T(i-1)) / dz;
        beta = - (k(T(i-1)) + k(T(i+1))) / dz^2;
        gamma = k(T(i+1)) / dz;

        % Calcola il termine sorgente
        qv = qv_max * cos(pi * z(i) / L);
        qv_next = qv_max * cos(pi * z(i+1) / L);
        qv_prev = qv_max * cos(pi * z(i-1) / L);
        q = (qv_next - qv_prev) / (2 * dz);

        % Calcola il termine di scambio termico
        qh = h * (T(i) - TCO2);

        % Calcola il residuo
        R(i) = alpha * T(i-1) + beta * T(i) + gamma * T(i+1) + q - qh;

        % Calcola la derivata del residuo rispetto a T(i-1), T(i) e T(i+1)
        dR_dT_prev = alpha;
        dR_dT = beta;
        dR_dT_next = gamma;

        % Assembla la matrice Jacobiana
        J(i, i-1) = dR_dT_prev;
        J(i, i) = dR_dT;
        J(i, i+1) = dR_dT_next;
    end

    % Condizioni al contorno
    J(1, 1) = 1; % Temperatura al centro della barretta
    J(Nz, Nz) = 1; % Temperatura alla superficie esterna
    R(1) = 0; % Condizione al contorno alla base
    R(Nz) = 0; % Condizione al contorno alla superficie esterna

end