% Esercitazione 11, esercizio 3
% Simone Canevarolo
% S269893
% 27/01/2024

clear all
close all
clc

%%

ll = 4e-1; % m
altezza = 6e-1; % m
cond = 1.5; % W/m/K
Tfissa = 200; % °C
hh = 50; % W/m^2/K
Taria = 30; % °C

dx = 1e-2;
dy = dx;

% lato verticale
xx = (0:dx:altezza)';

% lato orizzontale
yy = (0:dy:ll)';

Nx = length(xx);
Ny = length(yy);
Ntot = Nx*Ny;

[xmat,ymat] = meshgrid(xx,yy);

%%

AA = sparse([],[],[],Ntot,Ntot,5*Ntot);

for ii = 2:Nx-1
    for jj = 2:Ny-1

        kk = (jj-1)*Nx+ii;
        AA(kk,kk-Nx) = cond/dy^2;
        AA(kk,kk-1) = cond/dx^2;
        AA(kk,kk) = -2*cond*(1/dx^2+1/dy^2);
        AA(kk,kk+1) = cond/dx^2;
        AA(kk,kk+Nx) = cond/dy^2;

    end
end

bb = zeros(Ntot,1);

%%

% Neumann adiabatico lato ovest

jj = 1;
for ii = 2:Nx-1

    kk = ii;
    AA(kk,kk-1) = 1/dx^2;
    AA(kk,kk) = -2*(1/dx^2+1/dy^2);
    AA(kk,kk+1) = 1/dx^2;
    AA(kk,kk+Nx) = 2*(1/dy^2);

end

% Dirichlet lato nord

ii = 1;
for jj = 1:Ny

    kk = (jj-1)*Nx+ii;
    AA(kk,kk) = 1;
    bb(kk) = Tfissa;

end

% Robin lato est 

jj = Ny;
for ii = 2:Nx-1

    kk = Ntot-Nx+ii;
    AA(kk,kk-Nx) = 2*1/dy^2;
    AA(kk,kk-1) = 1/dx^2;
    AA(kk,kk) = -2*(1/dx^2+1/dy^2+hh/cond/dy);
    AA(kk,kk+1) = 1/dx^2;
    bb(kk) = -2*hh*Taria/dy/cond;

end

% Dirichlet lato sud

ii = Nx;
for jj = 1:Ny

    kk = (jj)*Nx;
    AA(kk,kk) = 1;
    bb(kk) = Tfissa;

end


TT=AA\bb;
TTT=reshape(TT,Nx,Ny);

surf(xmat',ymat',TTT,'FaceColor','interp')
colorbar
xlabel('x (m)')
ylabel('y (m)')
zlabel('T (K)')
title('Profilo di temperatura sulla sezione')