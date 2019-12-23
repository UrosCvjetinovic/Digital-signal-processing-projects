clear all
close all
clc
%%  -------1.1.------ Prvi deo, prva tacka

[ main, fsm] = audioread('..\dz2_signali\allplusA.wav');
[ background, fsb] = audioread('..\dz2_signali\noteA.wav');
figure(1);
    Tsm = 1/fsm;
    tm = 0:Tsm:(length(main)-1)*Tsm;
    plot(tm, main);
    title('Vremenski oblik signala ,,allplusA"');
    xlabel('Vreme [s]');

figure(2);
    Tsb = 1/fsb;
    tb = 0:Tsb:(length(background)-1)*Tsb;
    plot(tb, background);
    title('Vremenski oblik signala ,,notaA"');
    xlabel('Vreme [s]');


%%  ---------1.2.--------Prvi deo, druga tacka
% Pogledati funkciju: [b, a] = bandstop_filter_Cheb2(fs, fa, fp, Aa, Ap);
% Primenjena aproksimacija je Inverzna Cebisevljeva, i koriscenja je
% bilinearna transformacija za diskretizaciju

%%  ---------1.3.--------Prvi deo, treca tacka

nfft = 4096; 
window_width = nfft;
overlap_num = 3/4*window_width;

ws = hamming(window_width);

%racunanje spektrograma
[B,frequencies,times] = spectrogram(main, ws, overlap_num, nfft, fsm);
B_dB = 20*log10(abs(B)); %u dB

figure(13);
% prikaz spektrograma
subplot(1,2,1);
imagesc(times, frequencies(1:end/4), B_dB(1:end/4,:));
axis('xy');
xlabel('Vreme [s]');
ylabel('Ucestanost [Hz]');
title('Spektrogram signala allplusA');


%racunanje spektrograma
[B,frequencies,times] = spectrogram(background, ws, overlap_num, nfft, fsb);
B_dB = 20*log10(abs(B)); %u dB

% prikaz spektrograma
subplot(1,2,2);
imagesc(times, frequencies(1:end/4), B_dB(1:end/4,:));
axis('xy');
xlabel('Vreme [s]');
ylabel('Ucestanost [Hz]');
title('Spektrogram signala notaA');

%% ------1.4---------------Prvi deo, cetvrta tacka


mb = 1; % main/background signal
if( mb == 1)
    signal = main;
    fs = fsm;
else
    signal = background;   % Mogucnost provere kako niz filtara utice na notu A
    fs = fsb;
end
filterdSignal = signal;
TopFrek = 4500; % frekvencija do koje cemo filtrirati
fa01 = 210; fa02 = 230; INC  = 220; %opsezi na kojima se nalaze otisci note
%--mi podesavamo, uslov je max PZ = 100
        PZ  = 100; Aa   = 60; Ap   = 1;
%         PZ  = 40; Aa   = 30; Ap   = 0.05;
%         [ PZ Aa Ap]

 [b,a] = bandstop_filter_Cheb2(fs,[fa01 fa02], [fa01-PZ fa02+PZ], Aa, Ap);
 filtersb = {b};
 filtersa = {a};
 filterdSignal = filter(b,a,filterdSignal);
 for i = 1:floor((TopFrek-fa01)/INC)
    c = INC*i;
    [b,a] = bandstop_filter_Cheb2(fs, [fa01+c fa02+c], [fa01-PZ+c fa02+PZ+c], Aa ,Ap);
    filterdSignal = filter(b,a,filterdSignal);
    filtersb(i+1) = {b};
    filtersa(i+1) = {a};
 end
%--racunanje spektrograma
[B,frequencies,times] = spectrogram(filterdSignal, ws, overlap_num, nfft, fs);
B_dB = 20*log10(abs(B)); %u dB

%--prikaz spektrograma
figure;
imagesc(times, frequencies(1:end/4), B_dB(1:end/4,:));
axis('xy');
xlabel('Vreme [s]');
ylabel('Ucestanost [Hz]');
if(mb == 1)
    title('Spektrogram signala all' );
else
    title('Spektrogram signala none' );
end
%%
% Zvucna prvera;
if(mb == 1)
    sound(main,fsm);
    pause
    sound(filterdSignal,fsm);
    audiowrite('all.wav',filterdSignal,fsm);
else
     sound(background,fsb);
     pause
     sound(filterdSignal,fsb);
     audiowrite('none.wav',filterdSignal, fsb);
end
%% -------1.5.----------------Prvi deo, peta tacka
%NAPOMENA: UGASITI provera U FUNKCIJI bandstop_filter_Cheb2
%          inace ce iscrtavati mnogo grafika
 for i = 0:floor((TopFrek-fa01)/INC)
        if(i == 0)
            close(figure(51))
            close(figure(52))
        end
            c = INC*i;
            fp1 = fa01-PZ+c;
            fa1 = fa01+c;
            fa2 = fa02+c;
            fp2 = fa02+PZ+c;
            delta = INC - 2*PZ -fa2 + fa1; %ako koristimo PZ=100,delta je 0
            
            bd = cell2mat(filtersb(i+1));
            ad = cell2mat(filtersa(i+1));
            Nfreqz = 25000;
            [hd,Wd]=freqz(bd,ad,Nfreqz);
            Hd=abs(hd);
%---Amplitudska karakteristika digitalnog filtra u logaritamskoj razmeri
            f = Wd/(2*pi)*fs;
      figure(51);
            plot(f,20*log10(Hd), 'LineWidth', 2);
            title('Amplitudska karakteristika u logaritamskoj razmeri');
            xlabel('Ucestanost (Hz)');
            ylabel('|H(z)|');

            %crtanje gabarita
        hold on
            x_1 = [fp1-delta fp1]; y_1 = [-Ap -Ap];
            x_2 = [fp1 fp1]; y_2 = [-Ap -Aa];

            x_3 = [fp2 fp2+delta]; y_3 = [-Ap -Ap];
            x_4 = [fp2 fp2]; y_4 = [-Ap -Aa];

            x_5 = [fa1 fa2]; y_5 = [-Aa -Aa];
            x_6 = [fa1 fa1]; y_6 = [0 -Aa];
            x_7 = [fa2 fa2]; y_7 = [0 -Aa];

            plot(x_1,y_1,'r',x_2,y_2,'r',x_3,y_3,'r',x_4,y_4,'r',x_5,y_5,'r',x_6,y_6,'r',x_7,y_7,'r','LineWidth',1.5);

%---Amplitudska karakteristika digitalnog filtra u linearnoj razmeri
            f = Wd/(2*pi)*fs;
      figure(52);
            plot(f,Hd, 'LineWidth', 2);
            title('Amplitudska karakteristika u linearnoj razmeri');
            xlabel('Ucestanost (Hz)');
            ylabel('|H(z)|');

        %Crtanje gabarita
            delta_a = 10^(-0.05*Aa);
            delta_p = 1 - 10^(-0.05*Ap);
        hold on
            x_1 = [fp1-delta fp1]; y_1 = [1-delta_p 1-delta_p];
            x_2 = [fp1 fp1]; y_2 = [1-delta_p delta_a];

            x_3 = [fp2 fp2+delta]; y_3 = [1-delta_p 1-delta_p];
            x_4 = [fp2 fp2]; y_4 = [1-delta_p delta_a];

            x_5 = [fa1 fa2]; y_5 = [delta_a delta_a];
            x_6 = [fa1 fa1]; y_6 = [1 delta_a];
            x_7 = [fa2 fa2]; y_7 = [1 delta_a];

            plot(x_1,y_1,'r',x_2,y_2,'r',x_3,y_3,'r',x_4,y_4,'r',...
                x_5,y_5,'r',x_6,y_6,'r',x_7,y_7,'r','LineWidth',1.5);
 end
 hold off
 