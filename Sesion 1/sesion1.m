% SESION 1

% ORGANIZACION DE LA INTERFAZ DE USUARIO
%  1. Desktop --> Desktop Layout --> Default
%  2. Maximizar la interfaz
%  3. Anclar las ventanas de edición y gráfica

% Sesión 1: SEÑALES DISCRETAS


% Sesión 1 del laboratorio de TDSÑ

%% Apartado 1.1

% Eje "temporal" n
n=-10:20;

% Impulso unidad
d=double(n==0);

% Gráfica
stem(n,d,'k','MarkerfaceColor','k'), title('Impulso unidad \delta[n]')
xlabel('n')
axis([-10.5 20.5 -0.1 1.1])

%% Apartado 1.2
% Escalón unidad
u=double(n>=0);
stem(n,u,'k','MarkerfaceColor','k'), title('Escalón unidad u[n]')
xlabel('n')
axis([-10.5 20.5 -0.1 1.1])

%% Apartado 1.3
% Coseno
x=cos(2*pi/11*n);
stem(n,x,'k','MarkerfaceColor','k'), title('cos[\omega_0n]')
xlabel('n')
axis([-10.5 20.5 -1.1 1.1])

%% Apartado 1.4
% Coseno representado con plot y con grid
w0=2*pi/11;
x=cos(w0*n);
plot(n,x,'k','MarkerfaceColor','k'), title('cos[\omega_0n]')
grid
xlabel('n')
axis([-10 20 -1.1 1.1])

%% Apartado 1.5
% Coseno amortiguado
n=0:29;
x=cos(2*pi/11*n).*exp(-0.1*n);
plot(n,x,'k','MarkerfaceColor','k'), title('e^{-\alphan}cos[\omega_0n]')
xlabel('n')
axis([0 29 -1.1 1.1]), grid

%% Apartado 1.6
% Señal aleatoria con distribución uniforme entre -1 y 2
L=100; n=0:L-1;
x=3*rand(L,1)-1;
plot(n,x,'k','MarkerfaceColor','k'), title('r_u[n]')
xlabel('n')
axis([0 L-1 -1.1 2.1]), grid

%% Apartado 1.7
% Señal aleatoria con distribución normal, media 1 y desviación típica 2
L=100; n=0:L-1;
x=2*randn(L,1)+1;
plot(n,x,'k','MarkerfaceColor','k'), title('r_g[n]')
xlabel('n')
axis([0 L-1 1-3*2 1+3*2]), grid

% ENTRADA/SALIDA DE DATOS

%% Apartado 2.1
% Usando el terminal
m=input('Indique la masa del objeto: ');
v=input('Indique la velocidad del objeto: ');
disp(['La energía cinética del objeto vale ' num2str(1/2*m*v^2) ' julios'])
fprintf('La energia cinética del objeto vale %8.3f julios\n',1/2*m*v^2)

%% Apartado 2.2
% Usando el sistema de archivos
load a
n=0:length(a)-1;
plot(n,a), title('Forma de onda de /a/'), xlabel('n'), grid
a2=a.^2;
save a2 a2

%% Apartado 2.3
% Gráficamente
[Ind Amp]=ginput(2);
Fs=8000; Ts=1/Fs;
F1=1/((Ind(2)-Ind(1))*Ts);
fprintf('Primera resonancia aproximada %5.1f Hz\n',F1)

% EJERCICIOS PARA LOS ALUMNOS

%% Apartado 3.1
%Dibujar la señal x[n]=sen[((w0n)/pi)n] para w0=2pi/7 en el intervalo de
%índices -20 a 20
n=-20:20;
x=sin(2*pi/7*n)./(pi*n);
x(n==0)=2/7;
figure('Name','Apartado 3.1','NumberTitle','off');
stem(n,x,'k','MarkerfaceColor','k')
title('x(n)=sen(2\pi/7n)/\pin')
xlabel('n')
axis([n(1)-0.5 n(end)+0.5 min(x)-0.05 max(x)+0.05])
grid

%% Apartado 3.2
% Calcular y dibujar la respuesta impulsiva h[n]=a^n u[n] desde n=-5 a n=25
n=-5:25;
u=double(n>=0);
t=0.85.^n;
h=t.*u;
figure('Name','Apartado 3.2','NumberTitle','off');
stem(n,h,'k','MarkerfaceColor','k')
title('h(n)=0.85^n*u(n)');
axis([n(1)-0.5 n(end)+0.5 min(h)-0.1 max(h)+0.1])
grid
xlabel('n')
%% Apartado 3.3
%Cálculo de una convolución
n=0:19;
l=length(n);
x=double(n>=5 & n<=10);
n=0:29;
l_long=length(n);
h=0.85.^n;

y=conv(x,h);
l_2=length(y)-1;
m=0:l_2;
stem(m,y,'k','MarkerfaceColor','k')
title('Convolución lineal')
xlabel('n')
axis([-0.5 l+l_long-1 min(y)-0.55 max(y)+0.5])
grid


%% Apartado 3.4
% Cálculo de una serie de Fourier discreta
N=4; % Período
n=0:3*N-1; % abarca 3 períodos
w=2*pi/N;  % factor común
sumatorio=0;
    for k=0:3
        sumatorio=sumatorio + exp(1i*(pi/2)*k.*n);
    end

p=1/N*sumatorio;



%Si no quieres que te aparezca un warning puedes añadir real(p) dentro de
%stem. El resultado no cambia añadas eso o no.
stem(n,p,'k','MarkerfaceColor','k')
title('Tren de deltas periódico de periodo N=4')
ylabel('p(n)'),
xlabel('n'),grid on
delta=0.5;
axis([n(1)-delta n(end)+delta min(real(p))-0.1 max(real(p))+0.1])
