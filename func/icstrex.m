function [tiempo, cont, resultado] = icstrex(est, control, termo, inicial, entrada, cfg)

% Estructura
D = est(1);                            %ft  Diametro del cilindro
Ej = est(2);                           %ft  Espesor de la chaqueta
h = est(3);                            %ft  Altura del cilindro
hj = est(4);                           %ft  Altura de la chaqueta

% Control
cv = control(1);                       % gmp/psi^.5  Coeficiente válvula V001
cvj = control(2);                      % gmp/psi^.5  Coeficiente válvula V002
dPj = control(3);                      % psi         Caída de presión en la válvula V002
kp1 = control(4);
kd1 = control(5);
ki1 = control(6);
N1 = control(7);
kp2 = control(8);
kd2 = control(9);
ki2 = control(10);
N2 = control(11);
% set point
Tspi = control(12);                     % °F Valor inicial setpoint temperatura
Tspf = control(13);                     % °F Valor final setpoint temperatura
tsp_t = control(14);                    % s  Tiempo de cambio escalón temp.
hspi = control(15);                     % ft Valor inicial setpoint altura
hspf = control(16);                     % ft Valor final setpoint altura
tsp_h = control(17);                    % s  Tiempo de cambio escalón altura

% Prop. térmicas
mu = termo(1);                           % Btu/h°F ft2 Coeficiente de transferencia calor
rho = termo(2);                          % lbm/ft3     Densidad de la solución
cp = termo(3);                           % Btu/lbm °F  Calor específico de la solución
lambda = termo(4);                       % Btu/lbmol   Calor de reacción
rhoj = termo(5);                         % lbm/ft3     Densidad del fluido de enfriamiento
cpj =termo(6);                           % Btu/lbm °F  Calor específico del fluido de enfriamiento
E =  termo(7);                           % Btu/lbmol   Energia de activacion de la reacción
R = termo(8);                            % Btu L/mol°R Constante de los gases ideales
gw = termo(9);                           % ADM         Gravedad específica solción
gj = termo(10);                          % ADM         Gravedad específica del fluido de enfriamiento
% Const. físicas
gc = termo(11);                          % lbm ft/ lbf s2^2 factor de conversión
g = termo(12);                           % ft/s2 Constante de la gravedad


% Entradas
%Perturbaciones
% Flujo volumetrivo F0
t_sf0 = entrada(1);                      % s
f0i = entrada(2);                        % gpm
f0f = entrada(3);                        % gpm
%Concentracion Ca0
t_ca0 = entrada(4);                      % s
Ca0i = entrada(5);                       % gpm
Ca0f = entrada(6);                       % gpm
% Entradas fijas
T0 = entrada(7);                         % °F
Tj0 = entrada(8);                        % °F

% Configuración
t_ini = cfg(1);                          % s Tiempo inicial de simulacion
t_fin = cfg(2);                            % s Tiempo final de simulacion
sps = cfg(3);                              % Tamaño de paso


% CALCULO DE CONSTANTES
% Structure parameters calculation
A = (pi()*(D)^2)/4;                  % ft2 Area trans. del reactor                    
Atc = pi()*D*hj;                     % ft2 Area trans. de calor
vj = (((pi()*(D+2*Ej)^2)/4)-A)*hj;   % ft3 Volumen de la chaqueta

% Constantes
k0 = termo(13);                     % s^-1 C
k1 = 1/(448.83*A);                  % ft^-2
k2 = (mu*Atc)/(A*rho*cp);           % 
k3 = (lambda/(rho*cp));
k4 = 1/vj;
k5 = (mu*Atc)/(rhoj*cpj*vj);
k6 = (E/R);
k7 = cv*sqrt((rho*g)/(144*gw*gc));  
k8 = sqrt(dPj/gj)*cvj;              

constantes = [k0, k1, k2, k3, k4, k5, k6, k7, k8];


save test.mat;

simout = sim('cstrpid.slx','SimulationMode','normal',...
             'AbsTol','1e-5',...
             'SaveState','on',...
             'StateSaveName','xout',...
             'SaveOutput','on',...
             'OutputSaveName','yout',...
             'SaveFormat', 'StructureWithTime',...
             'SrcWorkspace','current',...
             'StartTime', num2str(t_ini),...
             'StopTime', num2str(t_fin),...
             'FixedStep',num2str(sps));
 
output = simout.yout;

tiempo = output.time;
resultado = output.signals(1).values;
cont = output.signals(2).values;      

end

