% load test.mat
clear; clc;
addpath('func/');
addpath('ui/');
addpath('src/');

%Perturbaciones
% Flujo volumetrivo F0
t_sf0 = 10;  %s
f0i = 5;     %gpm
f0f = 10;     %gpm

%Concentracion Ca0
t_ca0 = 12;  %s
Ca0i = 0.8;     %gpm
Ca0f = 1.6;     %gpm

% Entradas
T0 = 80;
Tj0 = 77;
entrada = [t_sf0, f0i, f0f, t_ca0, ...
           Ca0i, Ca0f, T0, Tj0];

%Condiciones iniciales
inicial = [2; ...  % h0  ft
          0.8; ... % Ca0 reactor lbmol/ft3
          0.0; ... % Cb0 reactor lbmol/ft3
          80; ...  % T0  °F
          77];     % Tj0 °F

% Data
% Estructura
D = 5.0;                            %ft  Diametro del cilindro
Ej = 0.5;                           %ft  Espesor de la chaqueta
h = 6.5;                            %ft  Altura del cilindro
hj = 5.0;                           %ft  Altura de la chaqueta
est = [D,Ej,h,hj];

% Prop. térmicas
% Transf. calor
mu = 1000;                           % Btu/h°F ft2 Coeficiente de transferencia calor
rho = 50;                            % lbm/ft3     Densidad de la solución
cp = 1;                              % Btu/lbm °F  Calor específico de la solución
rhoj = 62.3;                         % lbm/ft3     Densidad del fluido de enfriamiento
cpj = 1;                             % Btu/lbm °F  Calor específico del fluido de enfriamiento
gw = 0.78;                           % ADM         Gravedad específica solción
gj = 0.98;                           % ADM         Gravedad específica del fluido de enfriamiento
%Reaccion
E =  1500;                           % Btu/lbmol   Energia de activacion de la reacción
lambda = -60000;                     % Btu/lbmol   Calor de reacción
k0 = 10.28;                         % s^-1 C
% Const. físicas
gc = 32.2;                           % lbm ft/ lbf s2^2 factor de conversión
g = 32.2;                            % ft/s2 Constante de la gravedad
R = 1.985875;                        % Btu L/mol°R Constante de los gases ideales
termo = [mu,rho,cp,lambda,...
         rhoj,cpj,E,R,gw,...
         gj,gc,g,k0];

% Control
cv = 5.28;                           % gmp/psi^.5  Coeficiente válvula V001
cvj = 9.86;                          % gmp/psi^.5  Coeficiente válvula V002
dPj = 2;                             % psi         Caída de presión en la válvula V002

kp1 = -35;
kd1 = -0.1;
ki1 = -0.4;
N1 = 100;

kp2 = -50;
kd2 = 0.05;
ki2 = 0;
N2 = 100;

% set point
Tspi = 80.0;                         % °F Valor inicial setpoint temperatura
Tspf = 85.0;                         % °F Valor final setpoint temperatura
tsp_t = 4000;                        % s  Tiempo de cambio escalón temp.

hspi = 2.5;                          % ft Valor inicial setpoint altura
hspf = 3.0;                          % ft Valor final setpoint altura
tsp_h = 3600;                        % s  Tiempo de cambio escalón altura

control = [cv, cvj, dPj, kp1, kd1, ...
           ki1, N1, kp2, kd2, ki2, ...
           N2, Tspi, Tspf, tsp_t, ...
           hspi, hspf, tsp_h];

% Configuración
t_ini = 0;                          % s Tiempo inicial de simulacion
t_fin = 7200;                       % s Tiempo final de simulacion
sps = 0.03  ;                         % Tamaño de paso
cfg = [t_ini, t_fin, sps];

% % Structure parameters calculation
% A = (pi()*(D)^2)/4;                  % ft2 Area trans. del reactor                    
% Atc = pi()*D*hj;                     % ft2 Area trans. de calor
% vj = (((pi()*(D+2*Ej)^2)/4)-A)*hj;   % ft3 Volumen de la chaqueta
% 
% % Constantes
% k0 = 10.28;                         % s^-1 C
% k1 = 1/(448.83*A);                  % ft^-2
% k2 = (mu*Atc)/(A*rho*cp);           % 
% k3 = (lambda/(rho*cp));
% k4 = 1/vj;
% k5 = (mu*Atc)/(rhoj*cpj*vj);
% k6 = (E/R);
% k7 = cv*sqrt((rho*g)/(144*gw*gc));  
% k8 = sqrt(dPj/gj)*cvj;              
% 
% constantes = [k0, k1, k2, k3, k4, k5, k6, k7, k8];


[tiempo, cont, resultado] = icstrex(est,control, termo, inicial, entrada, cfg);