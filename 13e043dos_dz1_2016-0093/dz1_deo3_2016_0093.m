%%
clear all
clc
close all
%%
%--------------------------3.2

tic
    Rour=block_correlation(received_sequence, pn_sequence,4095);
ourtime=toc;
%%
%--------------------------3.3
tic
    Rtheir=xcorr(received_sequence,pn_sequence);
theirtime=toc;
figure(31);
    plot(1:length(Rour),Rour);
    title('block correlation');
    xlabel('Sample');

figure(32);
    plot(1:length(Rtheir),Rtheir);
    title('xcorr cross correlation');
    xlabel('Sample');



%%
%-------------------------3.4

[M,I] = max(Rour);
maximum=M;
index=I;

%%
%------------------------3.5
rawinput=received_sequence(index+1:end);
input=decode_sound(rawinput');

 %%
 audio = audioplayer(input,8000);

 audiowrite('govor.wav',input,8000);
% play(audio)
 
figure(351);

    s=spectrogram(input);
    spectrogram(input,'yaxis')
    title('Decoded input');
    xlabel('Sample');

figure(352);

    s1=spectrogram(received_sequence);
    spectrogram(received_sequence,'yaxis');
    title('Raw input');
    xlabel('Sample');
%%
%pseudo slucajna sekvenca
figure(222);
    stem(1:length(pn_sequence),pn_sequence);
    title('Pseudo random sequence');
    xlabel('Sample');
    axis([100 200 -2 2]);
    
figure(223);
    mid = floor(length(received_sequence)/2):floor(length(received_sequence)/2)+200;
    plot(mid,received_sequence(mid));
    title('Raw input signal');
    xlabel('Sample');

%%
%kros korelacija dva signala
figure(224);
    plot(1:length(Rour),Rour);
    title('Cross-Correlation');
    xlabel('Sample');




