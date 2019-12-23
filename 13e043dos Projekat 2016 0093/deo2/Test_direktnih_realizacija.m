clear all
close all
clc
%%
%---------Pravljenje slozene sinusoide i band pass filtra
fs = 44100;
F1 = 0.3;
F2 = 0.1;
F3 = 0.4;
N = 500;
n = 0:1:N-1;
x = 0.4*sin(2*pi*F1*n)+0.2*sin(2*pi*F2*n)+0.1*sin(2*pi*F3*n);
[b,a] = bandpass_filter( fs,[7800 16200],[8000 16000], 60, 0.05 );

%%

delta = [1 zeros(1,10000)];
odziv = filter(b, a, delta);
%fvtool(b,a);
%%
%---------Originalni signal
figure
subplot(311)
    X = abs(fft(x,length(x)));
    n = (0:1:length(x)-1)*1/N;
    stem(n,X), axis([0 1/2 0 100]);
    title('Spektar originalnog signala');

%---------Filtriranje signala pomocu IIR_direct
a1 = a(2:end);%Skrati se jer se podrazumeva da je prvi element jedan
y = IIR_direct_I(b,a1,x);
subplot(312)
    Y=abs(fft(y,length(y)));
    stem(n,Y),axis([0 1/2 0 100]);
    title('Spektar IIR direct I');


%---Fixed point filtriranje
wlength = 24; % moze sa 23
flength = 16;

% Podrazumevano je da se radi zaokruzivanje i da se radi zasicenje u
% slucaju prekoracenja
FixedPointAttributes=fimath ( 'RoundingMethod', 'Floor', 'OverflowAction', 'Saturate') ;
fi_params = struct('FILTER_COEFITIENTS_BITLENGTH',  wlength, 'FILTER_COEFITIENTS_FRAC', flength, ...
                   'SIGNAL_BITLENGTH',              wlength, 'SIGNAL_FRAC',             flength);

               
FI_b = fi( b , true , fi_params.FILTER_COEFITIENTS_BITLENGTH ,...
            fi_params.FILTER_COEFITIENTS_FRAC, FixedPointAttributes);
FI_a = fi( a1 , true , fi_params.FILTER_COEFITIENTS_BITLENGTH ,...
            fi_params.FILTER_COEFITIENTS_FRAC, FixedPointAttributes);
FI_x = fi( x , true , fi_params.SIGNAL_BITLENGTH ,...
            fi_params.SIGNAL_FRAC, FixedPointAttributes);
                
% Zasicenje smo radili samo ulaznom signalu i koeficijentima, u filtru
% zasicenje zahteva dodatni hardver, pa nakon zaokruzivanja, vracamo 'Wrap'        
FixedPointAttributes.OverflowAction = 'Wrap';
FI_b = fi( FI_b , true , fi_params.FILTER_COEFITIENTS_BITLENGTH ,...
            fi_params.FILTER_COEFITIENTS_FRAC, FixedPointAttributes);
FI_a = fi( FI_a , true , fi_params.FILTER_COEFITIENTS_BITLENGTH ,...
            fi_params.FILTER_COEFITIENTS_FRAC, FixedPointAttributes);
FI_x = fi( FI_x , true , fi_params.SIGNAL_BITLENGTH ,...
            fi_params.SIGNAL_FRAC, FixedPointAttributes);

            
y1=FI_IIR_direct(FI_b,FI_a,FI_x); 
Y=abs(fft(y1,length(y1)));
%n=[0:1:length(x)-1];
subplot(313)
    stem(n,Y),axis([0 1/2 0 100]);
    title('Spektar FI IIR direct I');
    xlabel('F / Fs');

FI_a=[1 FI_a];
FI_a=double(FI_a);
FI_b=double(FI_b);
Nfreqz = 200000;                                   
[hd,Wd]=freqz(FI_b,FI_a,Nfreqz);
F= Wd/(2*pi);  
Hd=abs(hd); 
figure              
    plot(F,20*log10(Hd), 'LineWidth', 2);
        title('Amplitudska karakteristika bandpass filtra');
        xlabel('Ucestanost (Hz)');
        ylabel('|H(z)|')
 hold on
    %fvtool(FI_b,FI_a);
    Nfreqz = 200000;                                   
    [hd,Wd]=freqz(b,a,Nfreqz);
    F= Wd/(2*pi);  
    Hd=abs(hd); 
    plot(F,20*log10(Hd), 'LineWidth', 2);
        title('Amplitudska karakteristika bandpass filtra');
        xlabel('Ucestanost (Hz)');
        ylabel('|H(z)|');

        
%% --------------------- Deo 2 tacka 4        
        
nule = roots(b);
polovi = roots(a);
figure
    [hz, hp, ht] = zplane(nule,polovi); 
    set(findobj(hz, 'Type', 'line'), 'LineWidth', 1, 'MarkerSize',5);
    set(findobj(hp, 'Type', 'line'), 'LineWidth', 1.5, 'MarkerSize',10);
    set(findobj(ht, 'Type', 'line'), 'LineWidth', 2);
    title('Raspored nula i polova za double preciznost')
    xlabel('Re(z)');
    ylabel('Im(z)');

nule = roots(FI_b);
polovi = roots(FI_a);
figure
    [hz, hp, ht] = zplane(nule,polovi); 
    set(findobj(hz, 'Type', 'line'), 'LineWidth', 1, 'MarkerSize',5);
    set(findobj(hp, 'Type', 'line'), 'LineWidth', 1.5, 'MarkerSize',10);
    set(findobj(ht, 'Type', 'line'), 'LineWidth', 2);
    title('Raspored nula i polova za fixed point preciznost')
    xlabel('Re(z)');
    ylabel('Im(z)');

    
    
 %% --------- Deo 2 Tacka 5
    
figure
    subplot(411)
        t = [0:1:length(x)-1]*1/fs;
        plot(t,x),axis([0 N/2*1/fs -0.6 0.6]);
        title('Originalni signal'); 
    subplot(412)
        plot(t,y),axis([0 N/2*1/fs -0.6 0.6]);
        title('Filtriran signal sa double preciznoscu');
    subplot(413)
        plot(t,y1),axis([0 N/2*1/fs -0.6 0.6]);
        title('Filtriran signal sa konacnom preciznoscu'); 
    subplot(414)
        deltay = y1-y;
        plot(t,deltay);
        title('Razlika filtriranih signala');
        
        