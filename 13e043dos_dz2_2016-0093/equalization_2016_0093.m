clear all
close all
clc
%%  -------2.1.------ Drugi deo, prva tacka
% Pogledati funkciju: [b, a] = bandpass_filter(fs, fa, fp, Aa, Ap);
% Primenjena aproksimacija je elipticka, i koriscenja je bilinearna
% transformacija za diskretizaciju

%%  -------2.2.------ Drugi deo, druga tacka
% Pogledati funkciju: [b, a] = highpass_filter(fs, fa, fp, Aa, Ap);
% Primenjena aproksimacija je elipticka, i koriscenja je bilinearna
% transformacija za diskretizaciju

%%  -------2.3.------ Drugi deo, treca tacka
fs = 44100;
PZ = 100;
Aa = 60;
Ap = 0.05;
% Filtriramo ih sa datim filtrom 9 i 10, kojim se menjanjem frek odabiranja
% dobijaju i filtri od 1 do 8
[bhp,ahp] = highpass_filter(fs,16000-PZ,16000,Aa,Ap);
fvtool(bhp,ahp);
[bbp,abp] = bandpass_filter(fs,[8000-PZ 16000+PZ], [8000 16000], Aa, Ap);
fvtool(bbp,abp);

%%  -------2.4.------ Drugi deo, cetvrta tacka
% Pogledati funkciju: [ y ] = IIR_equalizer( x, fs, style);
% Od ulaznog signala x je kreirano 10 signala, od kojih njih 9 se filtrira
% istim filtrom, na takav nacin sto se prvo smanjuje ucestanost odabiranja
% tog signala onoliko koliko odgovara odnosu propusnog opsega sa propusnim 
% opsegom projektovanog filtra 

%%  -------2.5.------ Drugi deo, peta tacka
clc
% NAPOMENA: ugasiti proveru unutar funkcija, iscrtava ce se preko 100
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

%%  -------2.6.------ Drugi deo, sesta tacka
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

%%  -------2.7.------ Drugi deo, sedma tacka
[ mysong, fms] = audioread('my_song.wav');
equalized_sound = IIR_equalizer(mysong(:,1),fms,'CUSTOM');
audiowrite('equalized_sound.wav',equalized_sound,44100);
%%
sound(mysong(:,1),fms)
pause
sound(equalized_sound,fms)
%%  -------2.8.------ Drugi deo, osma tacka
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

