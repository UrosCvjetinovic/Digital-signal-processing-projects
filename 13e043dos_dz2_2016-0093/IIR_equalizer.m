function [ y ] = IIR_equalizer( x, fs, style)
%IIR_equalizer zvucna ekvalizacija
%   x ulazni signal
%   fs ucestanost odabiranja signala
%   style tip ekvalizacije 'POP', 'ROCK', 'DANCE','CUSTOM'

% oznake blp,bbp,bhp znace b low pass-a, b band pass-a...

% provera=1 --->crtaju se amplitudske karakteristike filtra
provera = 0;   

type_of_style = {'POP', 'ROCK', 'DANCE','CUSTOM'};
rock  = [5, 3.75, 3, 1.5, -0.5, -1.5, 1, 2.5, 3.75, 4.5];
pop   =	[-1.5, -1, 0, 1.5, 4, 4, 2, 0, -1, -1.5];
dance = [4, 7, 5, 0, 2, 4, 5, 4.5, 3.5, 0];
custom = [ 15, 5, -10, 0, -10, 10.5, 5, 8, -40,-40];
preamp = 0; 
PZ = 100;
Aa = 60;
Ap = 0.05;

switch (find(ismember(type_of_style,style)))
    case 1
        tos = pop;
    case 2
        tos = rock;
    case 3
        tos = dance;
    case 4
        tos = custom;
    otherwise
        error('Pogresno pozvana funkcija (stil nije naveden pravilno)');
end
tos = tos + preamp;
db2num = 10.^(tos/20);

%----NADALJE: indeks celije i odgovara (11-i) filtru u zadatku i=1..10

% Ulazni signal
for i = 1:10
    X(i) = {x};
end

% Prolaze kroz NF filtar
pz = [47.5, 80, 100, 120, 120, 120, 120, 120, 120];
Xf(1) = X(1);%{cell2mat(X(1))}; %    Xf(1) == filter 10 
for i = 2:10
    f = 16000/(2^(i-2));
    [blp,alp] = lowpass_filter(fs,f+pz(11-i),f,Aa,Ap);
    Xf(i) = {filter(blp,alp,cell2mat(X(i)))};
end

% Ulazi su decimirani signal
Xfd(1) = Xf(1);
for i = 2:10
    Xfd(i) = {resample(cell2mat(Xf(i)),1,2^(i-2))};
end

f = 16000;
% Filtriramo ih sa datim filtrom 9 i 10, koji se menjanjem frek odabiranja
% dobijaju i filtri od 1 do 8
[bhp,ahp] = highpass_filter(fs,f-PZ,f,Aa,Ap);% filtar 10
[bbp,abp] = bandpass_filter(fs,[f/2-PZ f+PZ], [f/2 f], Aa, Ap);%filtar 9

 Y = {filter(bhp,ahp,cell2mat(Xfd(1)))};
for i = 2:10
    Y(i) = {filter(bbp,abp,cell2mat(Xfd(i)))};
end

% Izlaze dovodimo do iste frekvencije odabiranja
Z(1) = {cell2mat(Y(1))*db2num(10)};
Z(2) = {cell2mat(Y(2))*db2num(9)};
for i = 3:10
    Z(i) = {resample(cell2mat(Y(i)), 2^(i-2), 1)*db2num(11-i)};
end

% Izlaze saberemo
y = cell2mat(Z(1))';
for i = 2:10
    p = cell2mat(Z(i))';
    y = y + p(1:length(y)); %poslednjih par iteracija ima 8 odabiraka vise
end


if(provera == 1)
%%---SPEKTROGRAMI
        nfft = 4096; 
        window_width = nfft;
        overlap_num = 3/4*window_width;
        ws = hamming(window_width);

    %Racunanje spektrograma
        [B,frequencies,times] = spectrogram(x, ws, overlap_num, nfft, fs);
        B_dB = 20*log10(abs(B)); %u dB

    %Prikaz spektrograma
    subplot(2,1,1);
        imagesc(times, frequencies(1:end/4), B_dB(1:end/4,:));
        axis('xy');
        xlabel('Vreme [s]');
        ylabel('Ucestanost [Hz]');
        title('Spektrogram signala na ulazu u equalize');

    %Racunanje spektrograma
        [B,frequencies,times] = spectrogram(y, ws, overlap_num, nfft, fs);
        B_dB = 20*log10(abs(B)); %u dB

    %Prikaz spektrograma
    subplot(2,2,1);
        imagesc(times, frequencies(1:end/4), B_dB(1:end/4,:));
        axis('xy');
        xlabel('Vreme [s]');
        ylabel('Ucestanost [Hz]');
        title('Spektrogram signala na izlazu equalizera');
end

end

