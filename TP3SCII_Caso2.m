% Caso de estudio 2 - Avi鏮
clc; clear ; close all;

% Parametros
a = 0.07;
b = 5;
c = 150;
w = 9;

% Matrices   ;   X = [alfa phi phi_p h]
A = [-a a 0 0 ; 0 0 1 0 ; w^2 -w^2 0 0 ; c 0 0 0];
B = [0 ; 0 ; w^2*b ; 0];
C = [0 0 0 1];
D = 0;

% Polos del sistema
p1 = -15+15i;
p2 = -15-15i;
p3 = -0.5+0.5i;
p4 = -0.5-0.5i;
p = [p1 p2 p3 p4];

% Obtencion de K 
K = place(A, B, p); % con asignaci鏮 de polos

Q = diag([1 400000 1/10 1]);
R = 700000;
K = lqr(A,B,Q,R);   % con LQR

% Ganancia de referencia
G = -inv(C*inv(A-B*K)*B);

t_etapa = 1e-3;
tF = 70;
t=0:t_etapa:tF;
tiempo = round(tF/t_etapa);

%ent = 0*t+500; %referencia 1
ent = 0*t-500; %referencia 2

% Condidiones iniciales
x0 = [0 0 0 0];
alfa(1) = 0;
fi(1) = 0;
fi_p(1) = 0;
%h(1) = -100; %altura inicial 1
h(1) = 100; %altura inicial 2
u(1) = 0;

for i=1:1:tiempo
    x = [alfa(i); fi(i); fi_p(i); h(i)];
    u(i+1) = -K*x+ent(i)*G;

    % Ecuaciones diferenciales
    alfa_p = a*(fi(i)-alfa(i));
    fi_pp = -w^2 * (fi(i)-alfa(i)-b*u(i));
    h_p = c*alfa(i);

    alfa(i+1) = alfa(i) + t_etapa*alfa_p;
    fi_p(i+1) = fi_p(i) + t_etapa*fi_pp;
    fi(i+1) = fi(i) + t_etapa*fi_p(i);
    h(i+1) = h(i) + t_etapa*h_p;
end

figure; hold on; grid on;
plot(t,h);
plot(t,ent,'k');
title('Altura del avi鏮');
xlabel('Tiempo [s]');
ylabel('Altura [m]');

figure; hold on; grid on;
plot(t,alfa);
title('聲gulo con la horizontal');
xlabel('Tiempo [s]');
ylabel('聲gulo [rad]');

figure; hold on; grid on;
plot(t,fi);
title('聲gulo de cabeceo');
xlabel('Tiempo [s]');
ylabel('聲gulo [rad]');

figure; hold on; grid on;
plot(t,u);
title('Acci鏮 de control');
xlabel('Tiempo [s]');
ylabel('');

disp('Terminado')