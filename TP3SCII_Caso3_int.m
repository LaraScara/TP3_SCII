% Caso de estudio 3 - Péndulo
close all; clear all; clc;

% Parámetros
m = 0.1; 
F = 0.1; 
l = 1.6; 
g = 9.8; 
M = 1.5; 
tm = 1e-4;
tf = 15;

% Matrices ; X = [delta delta_p phi phi_p]
A = [0 1 0 0 ; 0 (-F/M) (-m*g/M) 0 ; 0 0 0 1 ; 0 (-F/(l*M)) (((M+m)*-g)/(l*M)) 0];
B = [0 ; (1/M) ; 0 ; (1/(l*M))];
C = [1 0 0 0; 0 0 1 0];
D = 0;

% Controlabilidad y observabilidad
Co = ctrb(A, B);
rank(Co) % = 4 por ende es controlable
Ob = obsv(A, C);
rank(Ob) % = 4 por ende es observable

% Discretización del sistema
sys = ss(A,B,C,D);
sys_d = c2d(sys,tm,'zoh');
Ad = sys_d.a; 
Bd = sys_d.b;
Cd = sys_d.c;

% Matrices ampliadas debido al integrador
Aa = [Ad zeros(4,1);-C(1,:)*Ad 1];
Ba = [Bd ; -C(1,:)*Bd];
Ca = [1 0];

% LQR
% Primer controlador para la masa inicial
Q1 = diag([1 6000 70000 1 .0001]);   
R1 = 500;
Ka1 = dlqr(Aa,Ba,Q1,R1);
K_i1 = -Ka1(5); 
K1 = Ka1(1:4);
% Segundo controlador para la masa 10 veces nayor
Q2 = diag([1 1000 50000 1 0.0001]);   
R2 = 150;
Ka2 = dlqr(Aa,Ba,Q2,R2);
K_i2 = -Ka2(5); 
K2 = Ka2(1:4);

% Implementación de funciones a usar
KMAX = tf/tm;
dt = tm/10; 
t = 0:dt:tf;
Vh = tm/dt;
u = [];
i = 1;
referencia1 = 10;  % posición delta de referencia 1
referencia2 = 0;   % posición delta de referencia 2
n = round(tf/dt);
ref = zeros(1,n);  
masa = zeros(1,n); 
for j=1:1:n+1
    if (j < (n)/2 )
        ref(j) = referencia1;
        masa(j) = m;
    else if(j >= (n)/2 )
        ref(j) = referencia2;
        masa(j) = m*10;
        end
    end
end

% Condidiones iniciales
Xop = [0; 0; pi; 0];
Ve = 0;
u_k(1) = 0; 
X = zeros(4,round(tf/dt)); 
X(:,1) = [0; 0; pi; 0];      
x_int = [0 0 pi 0]; 

% Iteración
for ki=1:KMAX
    %Dependiendo de la masa, se elige el controlador
    if masa(i)<0.2
        K = K1;
        K_i = K_i1;
    else
        K = K2;
        K_i = K_i2;
    end
 
    u_k(ki) = -K*(x_int'-Xop)+K_i*Ve;
    Y = C*(x_int'-Xop);
 
    for kii=1:Vh
        X_a = X(:,i)'; % X = [delta delta_p phi phi_p]
        u(i) = u_k(ki);
        Xp_1 = X_a(2);
        Xp_2 = -F*X_a(2)/M-masa(i)*g*(X_a(3)-pi)/M+u(i)/M;
        Xp_3 = X_a(4);
        Xp_4 = -F*X_a(2)/(M*l)-g*(M+masa(i))*(X_a(3)-pi)/(M*l)-u(i)/(M*l);
        Xp_a = [Xp_1 , Xp_2 , Xp_3, Xp_4];
        X(:,i+1) = X_a+ dt*Xp_a;
                
        i = i+1;
    end
 
    Ve = Ve+ref(i)-C(1,:)*x_int';
    x_int = X(:,i)';
end
u(i) = u_k(ki); 

% Gráficas
figure;hold on;
plot(t,X(1,:));grid on; hold on;
title('Desplazamiento del carro');
xlabel('Tiempo [s]');
ylabel('Posición [m]');
legend('(t)','referencia');
plot(t,ref,'k');

figure;
plot(t,X(3,:));grid on;
title('Ángulo del péndulo');
xlabel('Tiempo [s]');
ylabel('Ángulo [rad]');
legend('?(t)');
ylim([0 5]);

figure;
plot(t,u);grid on;
title('Acción de control');
xlabel('Tiempo [s]');
legend('u(t)');
ylim ([-30 30]);

disp('Terminado');
