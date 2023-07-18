% Caso de estudio 1 - Motor CC con observador y zona muerta
% clc; clear ; close all;

% Parámetros
Laa = 0.56*10^-3;
J = 0.0019;
Ra = 1.35;
Bm = 0.000792;
Ki = 0.1;
Km = 0.1;
Tl = 0.1035;
Ts = 1e-4;

% Matrices ; X = [ia ; w ; tita ];
A = [-Ra/Laa -Km/Laa 0; Ki/J -Bm/J 0; 0 1 0];
B = [1/Laa; 0; 0];
C = [0 1 0; 0 0 1];
D = [0];

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

% Matrices del sistema observador
C_o = B';
A_o = A';
B_o = C';

% LQR
QQ = diag([400 800 1 1/10]);  
RR = 80;
KK = dlqr(AA,BB,QQ,RR);
K = KK(1:3);
Ki = -KK(4);

Q_o = 1*diag([1/10 1/10 10]);    
R_o = [200 0; 0 300];
K_o = dlqr(A_o,B_o,Q_o,R_o);

% Implementación de funciones a usar
Ts = 1e-4; %Tiempo de muestreo
T = 10;
T_switch = 5;
deltat = Ts/10;
Kmax = T/Ts;
pasos = round(T/deltat);
t = 0:deltat:(T+T);
ta = 0:deltat:(T+0.3*T-deltat);
ref = (pi/2)*square(2*pi*ta/(2*T_switch));
fTl = -(Tl/2)*square(2*pi*ta/(2*T_switch))+Tl/2;

% Condiciones iniciales
Ci = [0 0 0 0];
x = zeros(4,pasos);
x(1,1) = Ci(1);
x(2,1) = Ci(2);
x(3,1) = Ci(3);
x(4,1) = Ci(4);
x_hat(1,1) = Ci(1);
x_hat(2,1) = Ci(2);
x_hat(3,1) = Ci(3);
x_ts = x((1:3),1);
v_ts = x(4,1);
z = 1;

% Iteración
for i=1:1:Kmax+1
    x_k = x_ts;
    v_k = v_ts;
    %u = -K(1:3)*x_k(1:3)+Ki*v_k; %Sin observador
    u = -K(1:3)*x_hat(1:3)+Ki*v_k; %Con observador
    
   % Alinealidad
   alin = 3;
    if abs(u)<alin
        u = 0;
    else
        u = sign(u)*(abs(u)-alin);
    end
    uu = u;
    ys = C*x(1:3,z);
    for j=1:1:Ts/deltat 
        ua(z) = uu;
        x1_p = -Ra*x(1,z)/Laa-Km*x(2,z)/Laa+u/Laa;
        x2_p = Ki*x(1,z)/J-Bm*x(2,z)/J-fTl(z)/J;
        x3_p = x(2,z);
        x_p_actual = [x1_p; x2_p; x3_p];
        
        x((1:3),z+1) = x((1:3),z)+deltat*x_p_actual;
        z = z+1;
    end
    yhat = C*x_hat;
    e = ys-yhat;
    x_hat = A*x_hat+B*u+K_o'*e;
    v_ts = v_ts+ref(z)-C(2,:)*x_ts;
    x_ts = x((1:3),z);
end

% Gráficas
figure(1); hold on; grid on;
plot(t(1:length(ua)),ua);
title('Acción de control');
zm = ones(2,length(ua));
zm(1,:) = zm(1,:)*alin;
zm(2,:) = zm(2,:)*-alin;
plot(t(1:length(ua)),zm,'k--');
xlabel('Tiempo [s]');
ylabel('Voltaje [V]');
legend('con observador','sin observador','zona muerta');
xlim([0 T]);

figure(2); hold on; grid on;
plot(t(1:length(x(1,:))),x(1,:));
title('Corriente de armadura');
ylabel('Corriente [A]'); xlabel('Tiempo [s]');
legend('con observador','sin observador');
xlim([0 T]);

figure(3); hold on; grid on;
plot(t(1:length(x(3,:))),x(3,:));  
plot(t(1:length(ref)),ref,'k');
title('Ángulo del motor');
xlabel('Tiempo [s]'); ylabel('Ángulo [rad]');
legend('con observador','referencia','sin observador');
xlim([0 T]);
%ylim([-2.5 2])

figure(4); hold on; grid on;
plot(t(1:length(x(2,:))),x(2,:));
title('Velocidad angular');
xlabel('Tiempo [s]'); ylabel('Velocidad angular [rad/s]');
legend('con observador','sin observador');
xlim([0 T]);

figure(5); hold on; grid on;
plot(x(3,:), x(2,:));
title('Plano de fase');
xlabel('Ángulo [rad]'); ylabel('Velocidad angular [rad/s]');
legend('con observador','sin observador');
xlim([-3 3]);

disp('Terminado');