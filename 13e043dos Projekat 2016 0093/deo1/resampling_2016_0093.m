clear all
close all
clc

N = 300;
n = 0:N-1;
fs = 44100;
x = sin(2*pi*0.1.*n) + sin(2*pi*0.6*n) + sin(2*pi*0.5*n);

%%slucaj kada je p > q

y = dos_resample_rat(x,5,3);
t1 = n/fs;
t2 = [0:length(y)-1] / (5/3*fs);
    figure
        stem(t1,x,'LineWidth',2);
    hold on
        stem(t2,y,'LineWidth',1);
            title('p=5 q=3');

%% slucaj kada je p < q

y1 = dos_resample_rat(x,3,5);
t2 = [0:length(y1)-1] / (3/5*fs);
figure
        stem(t1,x,'LineWidth',2);
    hold on
        stem(t2,y1,'LineWidth',1);
            title('p=3 q=5');



