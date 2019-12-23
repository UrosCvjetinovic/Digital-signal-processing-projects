close all
clear all
clc
%----------------------------------2.3

%crtanje spektra izracunatog preko nase funkcije fft
f1=1000;
f2=1240;
f3=1300;
f3s=10000;
T3s=1/f3s;

fnzd=gcd(gcd(f1,f2),f3);  % trazenje nzda
N=f3s/fnzd;
n=0:T3s:(N-1)*T3s;

x=5*cos(2*pi*f1*n)+1000*cos(2*pi*f2*n)+10*cos(2*pi*f3*n);
 
 X31=fft_radix_2(x);
 deltaF31=f3s/length(X31);  % rezolucija
 X31p=X31(1:length(X31)/2);
 k31=0:deltaF31:(length(X31)/2-1)*deltaF31;
 
 figure(231),
     subplot(211);
     stem(k31,(real(X31p)));
     title('Realni deo fft radix 2');
     ylabel('Real(X[k])');
     xlabel('k');

     subplot(212);
     stem(k31,(imag(X31p)));
     title('Imaginarni deo fft radix 2');
     ylabel('Img(X[k])');
     xlabel('k');
 
 X32=fft(x);
 X32p=X32(1:length(X32)/2);
 deltaF32=f3s/length(X32);
 k32=0:deltaF32:(length(X32)/2-1)*deltaF32;
 
 figure(232);
     subplot(211);
     stem(k32,real(X32p));  
     ylabel('Real(X[k])');
     xlabel('k');
     title('Realni deo fft'); 

     subplot(212);
     stem(k32,imag(X32p)); 
     ylabel('Img(X[k])');
     xlabel('k');
     title('Imaginarni deo fft');                       
 %Zakljucujemo da zavisi koliko tacaka uzimamo za fft (dodavanje nula)
 %ako uzmemo isti broj tacaka kao za nas fft_radix_2 algoritam   
 %spektri ce biti isti                        
  %%
  %-----------------------------2.4
  
  x41=ifft_radix_2(X31); 
  n41=0:length(x41)-1;
  figure(3)
      n4=0:length(x)-1;
      subplot(211);
      plot(n4,x);
      title('Dati niz x');
      xlabel('t[s]');
      subplot(212);
      plot(n41,x41);
      title('ifft radix 2(X)');
      xlabel('t[s]');
figure(33333)  
      plot(n41,x41);
      title('ifft radix 2(X)');
      xlabel('t[s]');
      
  x42=ifft(X31);
  n42=0:length(x42)-1;
  figure(5),
      plot(n42,real(x42));
      title('ifft(X)');
      ylabel(' x = ifft(X) ');
      xlabel('t[s]');
  
  
  %%
  %--------------------------------2.5
  
 [x5,fs] = audioread('chopin.wav');
 tstart=1;
 tend=1.25;
 elapsed= tstart - tend;
 Tst=30/length(x5');
 
 x51= x5(tstart/Tst:tend/Tst)';
 t1 = 0:Tst:(length(x51)-1)*Tst;
  
sound(x51);
figure(251);
    plot(0:Tst:(length(x5)-1)*Tst,x5);
    title('Izdvojen signal');
    xlim([1 1.5]);
    xlabel('t[s]');
%%
%Merenje vremena izvrsavanja naseg algoritma fft_radix_2
 tic
    Xc51=fft_radix_2(x51);
 time1=toc;
%%

%Merenje vremena izvrsavanja DFT-a
tic
     Xc52=zeros(1,length(x51));
     n=0:length(x51)-1;
          Ws=exp(-1j*2*pi/length(x51)*n);
          for k = 1:length(x51)
              Xc52(k) = sum(x51.*(Ws.^(k-1)));
          end
time2=toc;

 figure(252);
     plot(0:length(Xc51)-1,20*log10(abs(Xc51)));
     title('Sprektar sinala dobijen nasom funckijom fft radix 2');
     ylabel('|X[k]|[dB]');
     xlabel('k');
%% 
 figure(253);
     plot(0:length(Xc52)-1,20*log10(abs(Xc52)));
     title('Sprektar signala dobijen DFT-om po definiciji');
     ylabel('|X[k]|[dB]');
     xlabel('k');





