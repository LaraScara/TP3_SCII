% Caso de estudio 3 - P�ndulo
% close all; clear all; clc;

% Par�metros
m = 0.1; 
F = 0.1; 
l = 1.6; 
g = 9.8; 
M = 1.5; 
tm = 1e-4;
tf = 15;

% Matrices ; X = [delta delta_p phi phi_p]
A = [0 1 0 0 ; 0 (-F/M) (-m*g/M) 0 ; 0 0 0 1 ; 0 (-F/(l*M)) (((M+m)*-g)/(l*M)) 0];
B = [0 ; (1/M) ; 0 ; (-1/(l*M))];
C = [1 0 0 0; 0 0 1 0];
D = 0;

% Controlabilidad y observabilidad
Co = ctrb(A, B);
rank(Co) % = 4 por ende es controlable
Ob = obsv(A, C);
rank(Ob) % = 4 por ende es observable

% Discretizaci�n del sistema
sys = ss(A,B,C,D);
sys_d = c2d(sys,tm,'zoh');
Ad = sys_d.a; 
Bd = sys_d.b;
Cd = sys_d.c;

% Matrices ampliadas debido al integrador
Aa = [Ad zeros(4,1);-C(1,:)*Ad 1];
Ba = [Bd ; -C(1,:)*Bd];
Ca = [1 0];

% Matrices del sistema observador
Ao = Ad'; 
Bo = Cd';
Co = Bd';

% LQR
% Primer controlador para la masa inicial y su observador 
Q1 = diag([1 6000 70000 1 .0001]);   
R1 = 500;
Ka1 = dlqr(Aa,Ba,Q1,R1);
K_i1 = -Ka1(5); 
K1 = Ka1(1:4);
Qo1 = (diag([0.1 5000 2000 0.001])); 
Ro1 = (diag([4000 300]));
Ko1 = dlqr(Ao,Bo,Qo1,Ro1);
% Segundo controlador para la masa 10 veces mayor y su observador
Q2 = diag([1 1000 50000 1 0.0001]);   
R2 = 150;
Ka2 = dlqr(Aa,Ba,Q2,R2);
K_i2 = -Ka2(5); 
K2 = Ka2(1:4);
Qo2 = (diag([1 100 300 0.1])); 
Ro2 = (diag([100 10]));
Ko2 = dlqr(Ao,Bo,Qo2,Ro2);

% Implementaci�n de funciones a usar
KMAX = tf/tm;
dt = tm/10; 
t = 0:dt:tf;
Vh = tm/dt;
u = [];
i = 1;
referencia1 = 10;  % posici�n delta de referencia 1
referencia2 = 0;   % posici�n delta de referencia 2
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
err = 0; 
Xhat = [0; 0; pi; 0];

% Iteraci�n
for ki=1:KMAX
    %Dependiendo de la masa, se elige el controlador
    if masa(i)<0.2
        K = K1;
        K_i = K_i1;
        Ko = Ko1;
    else
        K = K2;
        K_i = K_i2;
        Ko = Ko2;
    end
 
 %   u_k(ki) = -K*(x_int'-Xop)+K_i*Ve;   % sin observador
    u_k(ki) = -K*(Xhat-Xop)+K_i*Ve;     % con observador
    Y = C*(x_int'-Xop);
    
    % Alinealidad
    alin = 0.3;
    if abs(u_k(ki))<alin
        u_k(ki) = 0;
    else
        u_k(ki) = sign(u_k(ki))*(abs(u_k(ki))-alin);
    end
 
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
    Yhat = C*(Xhat-Xop);
    err = Y-Yhat;
    Xhat = u(i-1)*Bd+Ko'*err+Ad*(Xhat-Xop)+Xop;
    Ve = Ve+ref(i)-C(1,:)*x_int';
    x_int = X(:,i)';
end
u(i) = u_k(ki); 

color = 'b';

% Gr�ficas
figure(1);hold on;
plot(t,X(1,:),color);grid on; hold on;
title('Desplazamiento del carro');
xlabel('Tiempo [s]');
ylabel('Posici�n [m]');
legend('con observador', 'sin observador', 'referencia');
plot(t,ref,'k');

figure(2); hold on;
plot(t,X(3,:),color);grid on;
title('�ngulo del p�ndulo');
xlabel('Tiempo [s]');
ylabel('�ngulo [rad]');
legend('con observador', 'sin observador');
ylim([0 5]);

figure(3); hold on;
plot(t,u,color);grid on;
title('Acci�n de control');
xlabel('Tiempo [s]');
legend('con observador', 'sin observador');
ylim ([-20 40]);

disp('Terminado');