%% ---- 1.1. 

% Pogledati funkciju: [h] = lowpass_filter_kaiser(Wc, Bt, Aa, Ap) 
% FIR NF filtar pomocu Kajzerove prozorske funkcije
% funkcija ima mogucnost graficke provere.
% Frekvencijske karakteristike racunate u 16384 tacaka.

%% ---- 1.2. 


% Pogledati funkciju:  [y] = dos_resample(x, up_down, r)
% 


%% ---- 1.3. 

% Pogledati funkciju:  [ y ] = IIR_equalizer( x, fs, style)
% 

%% ---- 1.4. 

clc
% NAPOMENA: ugasiti proveru unutar funkcija, iscrtavace se preko 100
%           grafika inace
fs25 = 44100;
Ts25 = 1 / fs25;
N = 22000;
n = 0:Ts25:Ts25*(N-1);
figure;
    imp = [ 1, zeros(1,N-1)];
    stem(n,imp); title('Impulsni signal'); xlabel('Vreme [s]');
figure; %POP, ROCK
    yp = IIR_equalizer(imp,fs25,'POP');
    subplot(2,1,1); stem(n,yp);  
    title('Impulsni odziv [POP]');
    xlabel('Vreme [s]');
    yr = IIR_equalizer(imp,fs25,'ROCK');
    subplot(2,1,2); stem(n,yr);  
    title('Impulsni odziv [ROCK]');
    xlabel('Vreme [s]');
figure; %DANCE, CUSTOM
    yd = IIR_equalizer(imp,fs25,'DANCE');
    subplot(2,1,1); stem(n,yd); 
    title('Impulsni odziv [DANCE]');
    xlabel('Vreme [s]');
    yc = IIR_equalizer(imp,fs25,'CUSTOM');
    subplot(2,1,2); stem(n,yc);  
    title('Impulsni odziv [CUSTOM]');
    xlabel('Vreme [s]');
    
%% ---- 1.5. 

clc
close all;
NFFT = N;
f = fs25/2*linspace(0,1,floor(NFFT/2));  % single-sided positive frequency
figure; %POP
    Yp = fft(yp,NFFT);
    n = 1:floor(length(Yp)/2);
    Ypp=Yp(n);
    subplot(211);
        semilogx(f,20*log10(abs(Ypp)));
        grid on
        title('Amplitudska karakteristika impulsnog odziva [POP]');
        xlabel('f[Hz]'); ylabel('|H(f)|[dB]');
    subplot(212);
        semilogx(f,angle(Ypp));
        grid on
        title('Fazna karakteristika impulsnog odziva [POP]');
        xlabel('f[Hz]'); ylabel('Phase[rad]');
figure; %ROCK
    Yr = fft(yr,NFFT);
    n = 1:length(Yr)/2;
    Yrp=Yr(n);
    subplot(211);
        semilogx(f,20*log10(abs(Yrp)));
        grid on
        title('Amplitudska karakteristika impulsnog odziva [ROCK]');
        xlabel('f[Hz]'); ylabel('|H(f)|[dB]');
    subplot(212);
        semilogx(f,angle(Yrp));
        grid on
        title('Fazna karakteristika impulsnog odziva [ROCK]');
        xlabel('f[Hz]'); ylabel('Phase[rad]');
figure; %DANCE
    Yd = fft(yd,NFFT);
    n = 1:length(Yd)/2;
    Ydp=Yd(n);
    subplot(211);
        semilogx(f,20*log10(abs(Ydp)));
        grid on
        title('Amplitudska karakteristika impulsnog odziva [DANCE]');
        xlabel('f[Hz]'); ylabel('|H(f)|[dB]');
    subplot(212);
        semilogx(f,angle(Ydp));
        grid on
        title('Fazna karakteristika impulsnog odziva [DANCE]');
        xlabel('f[Hz]'); ylabel('Phase[rad]');
figure; %CUSTOM
    Yc = fft(yc,NFFT);
    n = 1:length(Yc)/2;
    Ycp=Yc(n);
    subplot(211);
        semilogx(f,20*log10(abs(Ycp)));
        grid on
        title('Amplitudska karakteristika impulsnog odziva [CUSTOM]');
        xlabel('f[Hz]'); ylabel('|H(f)|[dB]');
    subplot(212);
        semilogx(f,angle(Ycp));
        grid on
        title('Fazna karakteristika impulsnog odziva [CUSTOM]');
        xlabel('f[Hz]'); ylabel('Phase[rad]');
        
%% ---- 1.6. 

[ mysong, fms] = audioread('my_song.wav');
equalized_sound = IIR_equalizer(mysong(:,1)',fms,'CUSTOM');
audiowrite('equalized_sound.wav',equalized_sound,44100);
    %%
    sound(mysong(:,1),fms)
    pause
    sound(equalized_sound,fms)

%% ---- 1.7. 

% Za crtanje grafika u ovoj tacki je iskoriscena provera unutar funkcija
% za projektovanje filtara. U svakoj funkciji bandpass_filtar,
% highpass_filtar, lowpass_filtar je promenljivoj provera dodeljena 1
% i nakon toga je pokrenuta sledeci kod (Isti filtri se projektuju da je
% bilo u pitanju ROCK, DANCE ili CUSTOM bilo pozvano)
% NAPOMENA: provera iscrtava tri grafika ili subplot zavisno kako je
% dodeljena vrednost indvorsubp
clc
close all
fs28 = 44100;
Ts28 = 1 / fs28;
N = 500;
n = 0:Ts28:Ts28*(N-1);
x = [ 1, zeros(1,N-1)];
y = IIR_equalizer(x, fs28,'POP');



%% ---- 1.8.  Nalaze se u resampling_2016_0093.m 
%% ---- 1.9. 