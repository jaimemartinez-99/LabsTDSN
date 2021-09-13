
%   Práctica 6: La DFT y Análisis Espectral de Sinusoides

%             ENRIQUE GONZALEZ MACÍAS
%        JAVIER LOPEZ INIESTA DÍAZ DEL CAMPO

%                   GRUPO 32.1
 
%% 1. MUESTREO DE LA TRANSFORMADA DE FOURIER

% 1.1  Calcule y dibuje el módulo y fase de las muestras X(e^jw) en el
%      periodo[0,2pi), para L=N=49

clear 

L=49;
N=49; 
w=2*pi*(0:N-1)/N;
X1 = (exp(-1j*w*((L-1)/2))).*(sin (w*L/2)./sin(w/2));
X1(1)=L;

figure('Name','Muestreo de la transformada de fourier para N = L','NumberTitle','off');
subplot(2,1,1)
    plot(w,10*log(abs(X1)));
    title('|X(2\pik/L)|');
    grid on;axis tight; ylabel('dB');xlabel('\omega');
subplot(2,1,2)
    plot(w,angle(X1));
    title ('< X(2\pik/L)');
    grid;axis tight;ylabel('rad');xlabel('\omega');
    
%% 1.2 Idem, pero tomando L=49 y N=300

L=49;
N=300; 
w=2*pi*(0:N-1)/N;
X2 = (exp(-1j*w*((L-1)/2))).*(sin (w*L/2)./sin(w/2));
X2(1)=L;

figure('Name','Muestreo de la transformada de fourier para N > L','NumberTitle','off');
subplot(2,1,1)
    plot(w,10*log(abs(X2)));
    title('|X(2\pik/L)|');
    grid on;axis tight; ylabel('dB');xlabel('\omega');
subplot(2,1,2)
    plot(w,angle(X2));
    title ('< X(2\pik/L)');
    grid;axis tight;ylabel('rad');xlabel('\omega');
    
    
%% 2. CALCULO DE LAS MUESTRAS DE UNA TRANSFORMADA DE FOURIER CON LA DFT

% 2.1 Repita los calculos del apartado 1.1 utilizando fft(...)
clear 

L=49;
N=49; 
w=2*pi*(0:N-1)/N;
signal=ones(1,L);
X3=fft(signal);% Muestrea con 49 muestras

figure('Name','Muestreo de la transformada de Fourier con la DFT para N = L','NumberTitle','off');
subplot(2,1,1)
    plot(w,10*log(abs(X3)));
    title('|X(2\pik/L)|');
    grid on;axis tight; ylabel('dB');xlabel('\omega');
subplot(2,1,2)
    plot(w,angle(X3));
    title ('< X(2\pik/L)');
    grid;axis tight;ylabel('rad');xlabel('\omega');
    
%% 2.2 Repita el apartado 1.2 usando el comando fft(x,N)
 clear

 L=49;
 N=300; 
 signal=ones(1,L);
 X4=fft(signal,N);% Muestreamos con 300 muestras 
 w=2*pi*(0:N-1)/N;
  
figure('Name','Muestreo de la transformada de Fourier con la DFT para N > L','NumberTitle','off');
subplot(2,1,1)
    plot(w,10*log(abs(X4)));
    title('|X(2\pik/L)|');
    grid on;axis tight; ylabel('dB');xlabel('\omega');
subplot(2,1,2)
    plot(w,angle(X4));
    title ('< X(2\pik/L)');
    grid;axis tight;ylabel('rad');xlabel('\omega');
     
%% 3. El comando MATLAB fftshift(...)

% 3.1  Repita los calculos del apartado 2.2 pero trabajando en el intervalo
%      de frecuencias (-pi, pi)

 w=2*pi*(0:N-1)/N -pi;
 X5=fftshift(X4);
 
figure('Name','Comando fftshift(...)','NumberTitle','off');
subplot(2,1,1)
    plot(w,10*log(abs(X5)));
    title('|X(2\pik/L)|');
    grid on;axis tight; ylabel('dB');xlabel('\omega');
subplot(2,1,2)
    plot(w,angle(X5));
    title ('< X(2\pik/L)');
    grid;axis tight;ylabel('rad');xlabel('\omega');
    
%% 4. Análisis espectral (caso determinista) con la DFT

clear

fs=20000;
L=100;
N=10000;
k=((0:N-1)/N)*(fs/2);
load('x1.mat');

% 4.1 Ventana Rectangular
xw1= x'.*(rectwin(L))';
X6=fft(xw1,N);% Muestreamos con 300 muestras 
  
% 4.2 Ventana de Hamming
xw2= x'.*(hamming(L))';
X7=fft(xw2,N);% Muestreamos con 300 muestras 

% Analisis espectral de x1.mat
figure('Name','Analisis espectral de x1.mat','NumberTitle','off');
subplot(2,1,1)
    plot(k,10*log(abs(X6)));
    title('x1: Ventana rectangular');
    grid on;axis tight; ylabel('dB');xlabel('k');
    [F1,A1]=ginput(2);
fprintf('\n Representando nuestra función vemos que en %0.0f Hz tiene una A de %0.0f dB\n', F1,A1)
subplot(2,1,2)
    plot(k,10*log(abs(X7)));
    title('x2: Ventana de Hamming');
    grid on;axis tight; ylabel('dB');xlabel('k');
    [F2,A2]=ginput(2)
fprintf('\n Representando nuestra función vemos que en %0.0f Hz tiene una A de %0.0f dB\n', F2,A2)

%% Ejercicio 5
clear

fs=1000;
N=10000;
load('x2.mat');
L=length(x);
k2=((0:N-1)/N)*(fs/2);
% ejeF=fs/M*(0:M-1)-fs/2;

% 5.1 Ventana Rectangular
xw3= x'.*(rectwin(L))';
X8=fft(xw3,N);% Muestreamos con 300 muestras 
  
% 5.2 Ventana de Hamming
xw4= x'.*(hamming(L))';
X9=fft(xw4,N);% Muestreamos con 300 muestras 

% Analisis espectral de x2.mat
figure('Name','Analisis espectral de x2.mat','NumberTitle','off');
subplot(2,1,1)
    plot(k2,10*log(abs(X8)));
    title('x2: Ventana Rectangular');
    grid;axis tight; ylabel('dB');xlabel('k')
subplot(2,1,2)
    plot(k2,10*log(abs(X9)));
    title('x2: Ventana de Hamming');
    grid;axis tight; ylabel('dB');xlabel('k')
    
% Representandolo vemos que P=2 y Q=1

[F3,A3]=ginput(3); % Para saber los pares A y F de los 2 cosenos y de la exponencial

fprintf('\n Representando nuestra función vemos que en %0.0f Hz tiene una A de %0.0f dB\n', F3,A3)