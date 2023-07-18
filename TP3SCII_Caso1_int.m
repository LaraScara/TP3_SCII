% Caso de estudio 1 - Motor CC
close all; clc;

% Importo datos
ruta_archivo = 'C:\Program Files\MATLAB\R2017b\bin\Curvas_Medidas_Motor_2023.xls';
hoja_trabajo = 'Hoja1';
mediciones = xlsread(ruta_archivo, hoja_trabajo);
% Columnas: tiempo , ángulo , velocidad angular , corriente , tension , torque

% Gráficos de curvas con los datos extraidos de la tabla
% figure % Gráfico del ángulo
% plot(mediciones(:,1),mediciones(:,2));
% title('Ángulo del motor'); legend('?(t)','Location','northeast');
% xlabel('Tiempo [s]'); ylabel('Ángulo [rad]');
% grid on; axis([0 10 -1 3]);
% 
% figure % Gráfico de la velocidad angular
% plot(mediciones(:,1),mediciones(:,3));
% title('Velocidad angular'); legend('wr(t)','Location','northeast');
% xlabel('Tiempo [s]'); ylabel('Velocidad angular [rad/s]');
% grid on; axis([0 10 -6 8]);
% 
% figure % Gráfico de la corriente de armadura
% plot(mediciones(:,1),mediciones(:,4)); 
% title('Corriente de armadura'); legend('ia(t)','Location','northeast');
% ylabel('Corriente [A]'); xlabel('Tiempo [s]');
% grid on; axis([0 10 -1 2]);
% 
% figure % Gráfico de la tension de entrada
% plot(mediciones(:,1),mediciones(:,5)); 
% title('Tensión de armadura'); legend('va(t)','Location','northeast');
% xlabel('Tiempo [s]'); ylabel('Voltaje [V]');
% grid on; axis([0 10 -1 2]);
% 
% figure % Gráfico del torque
% plot(mediciones(:,1),mediciones(:,6)); 
% title('Torque de carga'); legend('Tl(t)','Location','northeast');
% xlabel('Tiempo [s]'); ylabel('Torque [Nm]');
% grid on; axis([0 10 0 0.12]);

% Parámetros
Laa = 0.56*10^-3;
J = 0.0019;
Ra = 1.35;
Bm = 0.000792;
Ki = 0.1;
Km = 0.1;
Tl = 0.1035;

% Matrices ; X = [ia ; w ; tita ];
A = [-Ra/Laa -Km/Laa 0; Ki/J -Bm/J 0; 0 1 0];
B = [1/Laa; 0; 0];
C = [0 1 0; 0 0 1];
D = [0];

% Controlabilidad y observabilidad
Co = ctrb(A, B);
rank(Co) % = 3 por ende es controlable
Ob = obsv(A, C);
rank(Ob) % = 3 por ende es observable

% Implementación de funciones a usar
Ts = 1e-4; %Tiempo de muestreo
T = 10;
T_switch = 5;
deltat = Ts/10;
Kmax = T/Ts;
pasos = round(T/deltat);
t = 0:deltat:(T+T);
ta = 0:deltat:(T+0.3*T-deltat);
ref = (pi/2)*square(2*pi*ta/(2*T_switch)); % Función de referencia que varia entre pi/2 y -pi/2
fTl = -(Tl/2)*square(2*pi*ta/(2*T_switch))+Tl/2;% Función de torque que varia entre 0 y 0.1035

% Discretización del sistema
sys = ss(A,B,C,D);
sys_d = c2d(sys,Ts,'zoh');
A = sys_d.a;
B = sys_d.b;
C = sys_d.c;

% Matrices ampliadas debido al integrador
AA = [A,zeros(3,1);-C(2,:)*A, 1];
BB = [B;-C(2,:)*B];
CC = [C(2,:) 0];

% Controlabilidad sistema discreto ampliado
Co = ctrb(AA, BB);
rank(Co) % = 4 por ende es controlable

% LQR ; QQ = [ia ; w; tita; Ki];
QQ = diag([400 800 1 1/10]);  
RR = 80;
KK = dlqr(AA,BB,QQ,RR);
K = KK(1:3);
Ki = -KK(4);

% Condiciones iniciales
Ci = [0 0 0 0];
x = zeros(4,pasos);
x(1,1) = Ci(1);
x(2,1) = Ci(2);
x(3,1) = Ci(3);
x(4,1) = Ci(4);
x_ts = x((1:3),1);
v_ts = x(4,1);
z = 1;

for i=1:1:Kmax+1
    x_k = x_ts;
    v_k = v_ts;
    u = -K(1:3)*x_k(1:3)+Ki*v_k;
    
    for j=1:1:Ts/deltat 
        ua(z) = u;
        x1_p = -Ra*x(1,z)/Laa-Km*x(2,z)/Laa+u/Laa;
        x2_p = Ki*x(1,z)/J-Bm*x(2,z)/J-fTl(z)/J;
        x3_p = x(2,z);
        x_p_actual = [x1_p; x2_p; x3_p];
        
        x((1:3),z+1) = x((1:3),z)+deltat*x_p_actual;
        z = z+1;
    end
    v_ts = v_ts+ref(z)-C(2,:)*x_ts;
    x_ts = x((1:3),z);
end

% Gráficas
figure; hold on; grid on;
plot(t(1:length(ua)),ua,'r');
xlim([0 T]);
title('Acción de control');
xlabel('Tiempo [s]');
ylabel('Voltaje [V]');

figure; hold on; grid on;
plot(t(1:length(x(1,:))),x(1,:),'r');
plot(mediciones(:,1),mediciones(:,4),'b')
title('Corriente de armadura');
legend('ia(t)','referencia','Location','northeast');
ylabel('Corriente [A]'); xlabel('Tiempo [s]');
xlim([0 T]);

figure; hold on; grid on;
plot(t(1:length(x(3,:))),x(3,:),'r'); 
plot(t(1:length(ref)),ref,'k');
title('Ángulo del motor');
xlabel('Tiempo [s]'); ylabel('Ángulo [rad]');
legend('?(t)','referencia','Location','northeast');
xlim([0 T]);

figure; hold on; grid on;
plot(t(1:length(x(2,:))),x(2,:),'r');
plot(mediciones(:,1),mediciones(:,3),'b')
title('Velocidad angular');
legend('wr(t)','referecia','Location','northeast');
xlabel('Tiempo [s]'); ylabel('Velocidad angular [rad/s]');
xlim([0 T]);

disp('Terminado');