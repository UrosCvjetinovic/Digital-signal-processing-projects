clear all
close all
clc

params = struct('RoundingMethod', 'Ceiling', 'OverflowAction','Saturate',...
    'SIGNAL_BITLENGTH',21,'SIGNAL_FRAC', 10);
N = 1024;
F = 1/128;
W = 2*pi*F;
n = -N/2:1:N/2-1;

x = exp(1i*W*n);

X = abs(FI_FFT_radix_2(x,params));
k = (-N/2:1:N/2-1)*1/N*2*pi;
Y = [X(N/2+1:N) X(1:N/2)];
subplot(2,1,1);
    stem(k,Y),axis([-pi pi 0 1.1]);
        xlabel('W[rad]');
        ylabel('|X(W)|');
        title('Spektar signala sa FI FFT radix 2');
        
X1 = abs(fft(x,1024));
k1 = (-N/2:1:N/2-1)*1/N*2*pi;
Y1 = [X(N/2+1:N) X(1:N/2)];
subplot(2,1,2);
    stem(k,Y1),axis([-pi pi 0 1.1]);
        xlabel('W[rad]');
        ylabel('|X(W)|');
        title('Spektar signala sa fft');
        
%% -------------DEO 3. tacka 4.
    spurious = zeros(1,10);
    for i = 1:10
         params = struct('RoundingMethod', 'CEILING', 'OverflowAction','Saturate',...
    'SIGNAL_BITLENGTH',30,'SIGNAL_FRAC', 10+i);  
         X = abs(FI_FFT_radix_2(x,params));
         sorted = sort(X);
         sorted = double(sorted);
         spurious(i) = sorted(N)/sorted(N-1);
    end
    k = 11:20;
    spurious = 20*log(spurious);
    figure
        stem(k,spurious);
            xlabel('Broj bita u frakciji (broj bita za predstavu je 30 ) ');
            ylabel('SFDR');


 