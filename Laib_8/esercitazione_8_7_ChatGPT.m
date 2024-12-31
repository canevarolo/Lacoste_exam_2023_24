% Dati del problema
kk = 0.10; % Conducibilità termica [W/m/K]
thickness = 0.02; % Spessore della parete [m]
epsilon = 0.9; % Emissività del materiale
sigma = 5.67e-8; % Costante di Stefan-Boltzmann [W m^-2 K^-4]
T_ambient = 25; % Temperatura ambiente [°C]
flux = 75; % Flusso termico attraverso la parete [W/m^2]

% Parametri per il metodo delle iterazioni a punto fisso
max_iterations = 1000;
tolerance = 1e-6;
err = tolerance+1;

% Inizializzazione
T_guess = T_ambient; % Temperatura iniziale (può essere scelta arbitrariamente)

% Metodo delle iterazioni a punto fisso
while iteration<max_iterations && err>tolerance
    
    % Calcolo del flusso termico radiativo
    Q_rad = epsilon * sigma * (T_guess(end)^4 - T_ambient^4);
    
    % Calcolo del flusso termico conduttivo
    Q_cond = -kk * (T_guess(end) - T_ambient) / thickness;
    
    % Bilancio termico
    residual = Q_rad - Q_cond;
    
    % Aggiornamento della temperatura utilizzando il metodo di
    % Newton-Raphson (delle tangenti)
    T_guess(end + 1) = T_guess(end) - residual / (4 * epsilon * sigma * T_guess(end)^3);
    
    err = norm(TT-T_guess)/norm(TT-Tco2);
end

% Visualizzazione del risultato
xx = linspace(0, thickness, length(T_guess));
plot(xx, T_guess);
xlabel('Spessore della parete [m]');
ylabel('Temperatura [°C]');
title('Profilo di temperatura attraverso la parete');
grid on