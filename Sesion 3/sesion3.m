% Sesión 3: Muestreo

% Autor: Juan Santos
% Modificaciones: Gustavo Cuevas del Río

%% Muestreo de sinusoides

% Apartado 1.1

% Datos
Frecuencias=[200 850 1850 3800]; A=1;
Tv=10e-3;
Fs=8000; Ts=1/Fs;

% Rejilla temporal para el muestreo
t=0:Ts:Tv-Ts;

% Cálculo y dibujo de las sinusoides
n=1;
for F0=Frecuencias
    subplot(2,2,n)
    plot(t,A*cos(2*pi*F0*t))
    title(['cos(2\pi' num2str(F0,4) 't)']), xlabel('t'), grid
    n=n+1;
end
pause

% Apartado 1.2
% Datos
Frecuencias=[200 850 1850 5800];
Fs=6000; Ts=1/Fs;

% Rejilla temporal para el muestreo
t=0:Ts:Tv-Ts;

% Cálculo y dibujo de las sinusoides
n=1;
for F0=Frecuencias
    subplot(2,2,n)
    plot(t,A*cos(2*pi*F0*t))
    title(['cos(2\pi' num2str(F0,4) 't)']), xlabel('t'), grid
    n=n+1;
end
pause

% Apartado 1.3

% Cálculo y dibujo de las sinusoides
% ¿Por qué aparece un efecto parecido a una modulación de amplitud?
    % 
% ¿Por qué hay un cambio de signo de una de las señales en los apartados 1.2 y 1.3?
    % Porque en el 1.2 usa un coseno y en 1.3 usa un seno.
n=1;
for F0=Frecuencias
    subplot(2,2,n)
    plot(t,A*sin(2*pi*F0*t)),
    title(['sin(2\pi' num2str(F0,4) 't)']), xlabel('t'), grid
    n=n+1;
end

subplot

%% Interpolación

% Apartado 2
load tds
%wavplay(int16(tds), 24000)
p=audioplayer(int16(tds),24000);
play(p);
pause

% Apartado 2.1
L=3;
tdsL = zeros(L*length(tds), 1);
tdsL(1:L:length(tdsL)) = tds;
%wavplay(int16(tdsL), 24000)
p=audioplayer(int16(tdsL),24000);
play(p);

N=30;
Ncomienzo=1250; 
%n=L*Ncomienzo:L*Ncomienzo+N-1; %si el indexado de vectores comenzase en n=0
n=L*Ncomienzo+1:L*Ncomienzo+N; %Pero en MATLAB los índices comienzan en n=1
stem(n,tdsL(n)), title('Señal expandida por un factor L=3'), xlabel('n'), grid
pause

% Apartado 2.2
% Explique la operación del tercer comando MATLAB
    % Donde el seno vale NaN, la funcion isnan nos devuelve un 1. Es
    % equivalente a un filtro paso bajo.
nh = -100:100;
h = L*sin(pi/L*nh)./(pi*nh);
h(isnan(h)) = 1;
plot(nh,h), title('Respuesta impulsiva del filtro interpolador'), xlabel('n'), grid
pause
stem(-20:20,h( (nh>=-20) & (nh<=20) ))
title('Respuesta impulsiva del filtro interpolador en torno a n=0'), xlabel('n'), grid
axis([-20.5 20.5 -0.3 1.1])
pause

% Apartado 2.3
tdsLi=filter(h,1,tdsL);
%wavplay(int16(tdsLi),24000)
p=audioplayer(int16(tdsLi),24000);
play(p);


% Visualización de un segmento de 30 puntos de todas las señales
Npuntos=30;
Ncomienzo=1250;
%n=L*Ncomienzo:L*Ncomienzo+Npuntos-1; %si el indexado de vectores comenzase en n=0
n=L*Ncomienzo+1:L*Ncomienzo+Npuntos; %Pero en MATLAB los índices comienzan en n=1
nnn=L*Ncomienzo+1:L:L*Ncomienzo+Npuntos;
nn=Ncomienzo+1:Ncomienzo+length(nnn);
[hmax,nmax]=max(h);
N=abs(nh(nmax)-nh(1));
plot(nnn,tds(nn),'b')
hold
stem(n,tdsLi(n+N),'g')
stem(n,tdsL(n),'r')
hold
grid, xlabel('n')
legend('tds(t)','tdsLi[n]','tdsL[n]','Location','NorthWest')

%¿Cuánto vale ese retardo?
    %Ejecutando la sección podemos observar que N=100. La semilongitud del
    %filtro, ya que son 201 muestras

%% Diezmado
%¿Cuánto vale M?
    % Vale 3 porque la Ts que tenemos es 1/12Fo y la frecuencia máxima es
    % Fs= 1/4Fo; entonces el factor M que está multiplicando a 12F0 tiene
    % que ser 3.
    
% Apartado 3.1
% Datos
F0=1e3;
T1=-1; T2=1;

% Datos derivados
Ts=1/(12*F0);
t=T1:Ts:T2;

% Cálculo y dibujo de la señal
x=(sin(2*pi*F0*t)./(pi*t)).^2;
x(isnan(x))=4*F0^2;
n=round(t/Ts);
stem(-30:30,x( (n>=-30) & (n<=30))), title('(sin(2*\piF_0nT_s)/\pinT_s)^2'), xlabel('n'), grid
axis([-30.5 30.5 -0.05*4*F0^2 1.05*4*F0^2])
pause

% Apartado 3.2
[X w]=freqz(x,1,'whole');
plot(w/pi,abs(X)), title('|X(\omega)|'), xlabel('\omega (x \pi rad/muestra)'), grid
axis tight
pause

%Mida el máximo de|X(w)| y verifique que su valor es el correcto.
 
% Apartado 3.3
L=3;
xe=zeros(L*length(x),1);
xe(1:L:length(xe))=x;
[Xe w]=freqz(xe,1,'whole');
plot(w/pi,abs(Xe)), title('|Xe(\omega)|'), xlabel('\omega (x \pi rad/muestra)'), grid
axis tight
pause

%Se ha producido “aliasing” al comprimir laseñal x[n]? ¿Por qué?
    % Sí. Porque usa el valor M=2 cuando tendría que usar M=3. Esto hace
    % que no se cumpla Nyquist
    
% Apartado 3.4
M=2;
xd=x(1:M:length(x));
[Xd w]=freqz(xd,1,'whole');
plot(w/pi,abs(Xd)), title('|Xd(\omega)|'), xlabel('\omega (x \pi rad/muestra'), grid
axis tight


%% Cuantificación

% Apartado 4.1

% Datos
Fs=8000; Ts=1/Fs;
N=8;
Xm= 2^(16-1);% El fondo de escala es el mayor nivel que se puede reproducir
             % de la señal original que está cuantificada con 16 bits.
load tds

d = 2*Xm/2^N; % éste es delta, el tamaño del escalón para N bits y Xm
tdsQ = d*round(tds/d); %este es el cuantificador por redondeo al n. m. c.
e= tdsQ-tds;  % la señal de error
n=0:length(tds)-1;
subplot(3,1,1), plot(n,tds), title('Señal tds[n]'), grid, xlabel('n')
axis tight
subplot(3,1,2), plot(n,tdsQ), title('Señal cuantificada tdsQ[n]'), grid, xlabel('n')
axis tight
subplot(3,1,3), plot(n,e), title('Señal e[n]'), grid, xlabel('n')
axis([0 n(end) -d/2*1.1 d/2*1.1])
subplot

%wavplay([int16(tds)' zeros(1,4000) int16(tdsQ)' zeros(1,4000) int16(e)'],Fs)
p=audioplayer([int16(tds)' zeros(1,4000) int16(tdsQ)' zeros(1,4000) int16(e)'],Fs);
play(p)
pause

% Apartado 4.2
%SNRglobal=10*log10(sum(tds.^2)/sum(e.^2));
SNRglobal=10*log10((tds'*tds)/(e'*e));
fprintf('SNRglobal = %5.2f dB\n',SNRglobal)
pause

%Explique el origen de la discrepancia.
    %Porque falta el término de adaptación
    
% Apartado 4.3
% Datos
Tv=20e-3;

Lv=Tv/Ts; Lsol= 0.4*Lv; Ld=Lv-Lsol;

% Cálculo de la energía localizada de la señal x
Ploctds=[];
tdsl=slocal(tds,Lv,Ld,0);
while 1
    tdsl=slocal(tds,Lv,Ld,1);
    if isempty(tdsl)
        break,   % Se alcanzó el final de la señal tds
    end
    Ploctds(end+1)=sum(tdsl.^2)/Lv;
end

% Cálculo de la energía localizada del ruido e(n)
Ploce=[];
el=slocal(e,Lv,Ld,0);
while 1
    el=slocal(e,Lv,Ld,1); % segmento local de ruido
    if isempty(el)
        break,   % Se alcanzó el final de la señal tds
    end
    Ploce(end+1)=sum(el.^2)/Lv;
end

SNRseg=10*log10(Ploctds./Ploce);

t=(0:length(tds)-1)*Ts;
subplot(211), plot(t,tds), title('Señal original tds'), xlabel('t'), grid
axis tight
tSNRseg=(0:length(SNRseg)-1)*Ld*Ts;
subplot(212), plot(tSNRseg,SNRseg,'g',[tSNRseg(1) tSNRseg(end)],[SNRglobal SNRglobal],'r')
title('SNR segmental en verde y SNRglobal en rojo'), xlabel('t'), grid
axis tight
subplot
pause

%Comente los resultados
    %Deberíamos obtener una distribución uniforme.
    
%¿A qué se debe que la SNR segmental alcance valores tan altos?
    %Porque hay un solapamiento del 40%

% Apartado 4.4
M = 50;
[h, b] = hist(e, M);
bar(b, h/(length(e)*d/M));
hold on
plot(b, (1/d)*ones(size(b)),'g')
title('Diagrama de barras del histograma de amplitudes y densidad de distribución de amplitudes teórica')
xlabel('e')
hold off
pause
%Explique la necesidad del factor de normalización del histograma 1/(length(e)*d/M).
    %El área de la función ha de valer 1. Por lo tanto, queremos que el diagrama de barras nos de resultados en torno a
    %+1/delta

% Apartado 4.5
ree=xcorr(e)/length(e);
[Mx Ind]=max(ree);
m=-30:30;
subplot(2,1,1)
stem(m,ree(m+Ind))
title('\phi_{ee}[m]'), xlabel('m'), grid

%Estimación del espectro de potencia:
%[Pee w]=periodogram(e);
%Pee=pi*Pee; % Suprimimos las unidades que introduce MATLAB (units of power per radians per sample)

%Estimación alternativa del espectro de potencia:
%(considerando que la autocorrelación es una secuencia determinista)
[Pee w]=freqz(ree,1,'whole');%Transformada de Fourier de la autocorrelación
Pee=abs(Pee); 

Pote=10*log10(d^2/12);
subplot(2,1,2)
plot(w/pi,10*log10(Pee),'g',[w(1) w(end)]/pi,[Pote Pote],'r')
title('Estimación de la densidad espectral de potencia de e[n] en verde, valor teórico en rojo')
xlabel('\omega (x \pi rad/muestra)'), ylabel('dB'), grid
subplot
axis tight

% EJERCICIOS ADICIONALES

%% Apartado 2.1
%
% Datos
errormax=1e-3;
V=10;
Tmin=0; Tmax=5;


%% Apartado 2.2
%
% Datos
cardio100;



%% Apartado 2.3

% 2.3.1
[x Fs]=audioread('Invention8.wav');
