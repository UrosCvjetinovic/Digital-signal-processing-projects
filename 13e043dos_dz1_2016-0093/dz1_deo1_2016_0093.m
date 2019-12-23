%%-----1.1
clear all
clc
close all

f1 = 1e3;       % Frekvencija f1
f2 = 1.24e3;    % Frekvencija f2
f3 = 1.3e3;     % Frekvencija f3
fs = 10e3;      % Frekvencija diskretizovanja
Ts = 1/fs;

fnzd=gcd(gcd(f1,f2),f3);      % nzd(f1,f2,f3)
N=fs/fnzd;                    % trazeni broj odabiraka

t = 0:Ts:(N-1)*Ts;            % 
% signal x(t) mozemo jedino diskretno prikazati u matlabu 
x=5*cos(2*pi*f1*t)+1000*cos(2*pi*f2*t)+10*cos(2*pi*f3*t);
 

%%
 figure(11), 
 plot(t,x)
 title('Vremenski oblik signala');
 xlabel('t[s]');
 
%DFT po def
    X1=zeros(1,N);
    n=0:N-1;    
    Ws=exp(-1j*2*pi*n/N);   %twidle
    for k = 1:N
        X1(k) = sum(x.*(Ws.^(k-1)));
    end
   
deltaF=fs/length(X1);   %frekvencijska rezolucija
 
%jednostrani spektar
X1p=X1(1:N/2);
k1=0:deltaF:deltaF*(N/2-1);

% realni i imaginarni deo
figure(12),
subplot(2,1,1);
stem(k1,real(X1p));
title('Realni deo');
ylabel('Re(X[k])');
xlabel('k');

subplot(2,1,2);
plot(k1,imag(X1p));
title('Imaginarni deo');
ylabel('Img(X[k])');
xlabel('k');

figure(13),stem(k1,abs(X1p));
title('Spektar diskretizovang signala zumiran');
ylabel('|X[k]|');
xlabel('k');


%%
%-----------1.2

w1=triang(N)';    %trougaounu funkciju

x2=x.*w1;

%DFT po def
    X2=zeros(1,N);
    n=0:N-1;    
    Ws=exp(-1j*2*pi*n/N);   %twidle
    for k = 1:N
        X2(k) = sum(x2.*(Ws.^(k-1)));
    end

X2p=X2(1:N/2);
k2 = 0:deltaF:deltaF*(N/2-1); 
 % realni i imaginarni deo
figure(21),
subplot(2,1,1);
stem(k2,real(X2p));
title('Realni deo');
ylabel('Re(X[k])');
xlabel('k');

subplot(2,1,2);
stem(k2,imag(X2p));
title('Imaginarni deo');
ylabel('Img(X[k])');
xlabel('k');

%%
%------------1.3

%Signal se odabira u 1.3N tacaka 
N3=1.3*N;
deltaF3=fs/N3;
t3=0:Ts:(N3*Ts-Ts);
x3=5*cos(2*pi*f1*t3)+1000*cos(2*pi*f2*t3)+10*cos(2*pi*f3*t3);
n3=0:N3-1;

k3 = 0:deltaF3:(floor(N3/2)-1)*deltaF3;
X3=zeros(1,N3);
 
 Ws1=exp(-1j*2*pi/N3*n3);
 for k = 1:N3
     X3(k) = sum(x3.*(Ws1.^(k-1)));
 end

figure(31),plot( k3 ,20*log10(abs(X3(1:floor(N3/2)))));
title('Amplitudska karakteristika x u 1.3*N tacaka ')
xlabel('k');
ylabel('|X|[dB]');
    

%%
    % a) pravougaona prozorska funkcija

 wa=ones(1,N3);
 y3a=x3.*wa;
 Y3a=zeros(1,N3);

 for k = 1:N3
     Y3a(k) = sum(y3a.*(Ws1.^(k-1)));
 end
Y3ap=Y3a(1:N3/2); 
 
figure(321);
subplot(2,2,1);
plot(k3 ,20*log10(abs(Y3ap)));
title('Pravougaona');
ylabel('|X[k]|[dB]');
xlabel('k');
grid on;

    % b) trougaona prozorska funkcija
    
 wb=triang(N3)';
 y3b=x3.*wb;
 Y3b=zeros(1,N3);

 for k = 1:N3
     Y3b(k) = sum(y3b.*(Ws1.^(k-1)));
 end

 Y3bp=Y3b(1:N3/2); 
 
subplot(2,2,2);
plot(k3,20*log10(abs(Y3bp)));
title('Trougaona');
ylabel('|X[k]|[dB]');
xlabel('k');
grid on;

    % c) Hanova prozorska funkcija
 wc=hann(N3)'; 
 y3c=x3.*wc;
 Y3c=zeros(1,N3);

 for k = 1:N3
     Y3c(k) = sum(y3c.*(Ws1.^(k-1)));
 end

Y3cp=Y3c(1:N3/2); 
 
subplot(2,2,3);
plot(k3,20*log10(abs(Y3cp)));
title('Hanova');
ylabel('|X[k]|[dB]');
xlabel('k');
grid on;


    % d) Hemingova prozorska funkcija
wd=hamming(N3)';
y3d=x3.*wd;
Y3d=zeros(1,N3);

 for k = 1:N3
     Y3d(k) = sum(y3d.*(Ws1.^(k-1)));
 end

 Y3dp=Y3d(1:N3/2); 
 
subplot(2,2,4);
plot(k3,20*log10(abs(Y3dp)));
title('Hemingova');
ylabel('|X[k]|[dB]');
xlabel('k');
grid on;

    % e) Blekmenova prozorska funkcija

 we=blackman(N3)'; 
 y3e=x3.*we;
 Y3e=zeros(1,N3);
 
 for k = 1:N3
     Y3e(k) = sum(y3e.*(Ws1.^(k-1)));
 end

 Y3ep=Y3e(1:N3/2); 
 
figure(322);
subplot(2,2,1);
plot(k3,20*log10(abs(Y3ep)));
title('Blekmenova');
ylabel('|X[k]|[dB]');
xlabel('k');
grid on;


    %pod f),g) i h) kajzerova prozorska funkcija za beta=3,7,10
kai1=kaiser(N3,3)';
kai2=kaiser(N3,7)';
kai3=kaiser(N3,10)';
    
    y3f=x3.*kai1;
    y3g=x3.*kai2;
    y3h=x3.*kai3;
    Y3f=zeros(1,N3);
    Y3g=zeros(1,N3);
    Y3h=zeros(1,N3);
 Ws1=exp(-1j*2*pi/N3*n3);
 for k = 1:N3
     Y3f(k) = sum(y3f.*(Ws1.^(k-1)));
     Y3g(k) = sum(y3g.*(Ws1.^(k-1)));
     Y3h(k) = sum(y3h.*(Ws1.^(k-1)));
 end

 Y3fp=Y3f(1:N3/2); 
 Y3gp=Y3g(1:N3/2); 
 Y3hp=Y3h(1:N3/2); 
 
subplot(2,2,2);
plot(k3,20*log10(abs(Y3fp))); 
title('Kajzerova Beta=3');
ylabel('|X[k]|[dB]');
xlabel('k');
grid on;

subplot(2,2,3);
plot(k3,20*log10(abs(Y3gp)));
title('Kajzerova Beta=7');
ylabel('|X[k]|[dB]');
xlabel('k');
grid on;

subplot(2,2,4);
plot(k3,20*log10(abs(Y3hp)));
title('Kajzerova Beta=10');
ylabel('|X[k]|[dB]');
xlabel('k');
grid on;

%%
%-------------1.4

 f1=1000;
 f2=1240;
 f3=1300;
 N4=1000;
 fs=10000;
 deltaF4=fs/N4*8;
 T4s=1/fs;
 t4=0:T4s:(N4-1)*T4s;
 x4=5*cos(2*pi*f1*t4)+1000*cos(2*pi*f2*t4)+10*cos(2*pi*f3*t4);    
 k4 = deltaF4:deltaF4:fs/2;
 
 %Uzimanje svakog osmog
 xd4 = downsample(x4,8);
 
 X4=fft(xd4);
 X4p=X4(1:floor(length(X4)/2));
 
 figure(4);
 subplot (211);
 stem(k1,20*log10(abs(X1p)));
 title('Amp.Spek. datog signal');
 ylabel('|X[k]|[dB]');
 xlabel('k');
 
 subplot (212);
 stem(k4,20*log10(abs(X4p)));
 title('Amp.Spek. dobijenog odabiranjem datog');
 ylabel('|X[k]|[dB]');
 xlabel('k');
 
%%
%------------------1.5
t0 = 0; tk = 5;
f5s=8000;
T5s=1/f5s;

c = t0:T5s:tk;
xchirp=chirp(c,t0,tk,f5s/2,'linear');

figure(15),
spectrogram(xchirp,'yaxis');
audio51=audioplayer(xchirp,f5s);
% play(audio);

%%
%--------------------------------1.6

xchirpd2= downsample(xchirp,2);
audio61=audioplayer(xchirpd2,f5s/2);
 %play(audio61);

xchirpd5= downsample(xchirp,5);
audio62=audioplayer(xchirpd5,f5s/5);
 %play(audio62);

%%
%--------------------------------------1.7
figure(171);
subplot(311);
s71=spectrogram(xchirp,f5s);
spectrogram(xchirp,[], [] , f5s/2,f5s,'yaxis');

title('xchirp');
ylabel('');
xlabel('');

subplot(312);
s72=spectrogram(xchirpd2,f5s);
spectrogram(xchirpd2,[], [] , f5s/2/2,f5s/2,'yaxis');
title('xchirp,d2');
xlabel('');

subplot(313);
s73=spectrogram(xchirpd5,f5s);
spectrogram(xchirpd5,[], [] , f5s/2/5,f5s/5,'yaxis');
title('xchirp,d5');
ylabel('');
xlabel('Time');


% figure(172);
% subplot(311);
% s71=spectrogram(xchirp,f5s);
% spectrogram(xchirp,'yaxis');
% 
% title('xchirp');
% ylabel('');
% xlabel('');
% 
% subplot(312);
% s72=spectrogram(xchirpd2,f5s);
% spectrogram(xchirpd2,'yaxis');
% title('xchirp,d2');
% xlabel('');
% 
% subplot(313);
% s73=spectrogram(xchirpd5,f5s);
% spectrogram(xchirpd5,'yaxis');
% title('xchirp,d5');
% ylabel('');
% xlabel('Sample');
