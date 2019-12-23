function [ y ] = lowpass_Parks_McClellan( Wc,Bt,Aa,Ap)
%lowpass_Parks_McClellan Projektuje FIR NF filtar koriscenjem Parks
% MekKlenov postupak
% Wc granicna ucestanost
% Aa slabljenje u nepropusnom
% Ap slabljenje u propusnom

% provera=1 --->crtaju se amplitudske karakteristike filtra
provera = 0;

if ((Ap > Aa) || (Bt > Wc) || (Wc >= pi))  
    error('Error: Greska pri pozivanju lowpass_Parks_McClellan');
end

Nfreqz = 25000;

delta_p = ( 10^(0.05*Ap) - 1 )/( 10^(0.05*Ap) + 1 );	
delta_a = 10^(-0.05*Aa);
Wp = Wc - Bt/2; % Wp = fp/ (fs/2) *pi
Wa = Wc + Bt/2;

%---izracunavanje potrebnog reda filtra
D = (0.005309*log10(delta_p)*log10(delta_p)+0.07114*log10(delta_p)-0.4761)...
    *log10(delta_a);
D = D - (0.00266*log10(delta_p)*log10(delta_p)+0.5941*log10(delta_p)+0.4278);
f = 11.01217+0.51244*(log10(delta_p)-log10(delta_a));
M = 2*pi*D/Bt-f*Bt/(2*pi) + 1;
M = ceil(M);		% duzina impulsnog odziva
N = M-1;			% red filtra

%---projektovanje filtra na osnovu zadatih parametara
Fp_n = Wp/pi;	% normalizacija ucestanosti (fp/(fs/2))
Fa_n = Wa/pi;
Hd = [1   1    0   0];	% zeljena amplitudska karakteristika
F = [ 0 Fp_n Fa_n  1];	% vektor normalizovanih ucestanosti za koje se zadaje Hd(w)
y = firpm(N,F,Hd);	

[H,w] = freqz(y,1,Nfreqz);
% w = w*fs/2/pi;
Hfa = abs(H);
Hfp = unwrap(angle(H));


ip = ceil((Nfreqz*Wp)/(pi))+1;
ia = floor((Nfreqz*Wa)/(pi))+1;
Ha = Hfa(ia:end);
Hp = Hfa(1:ip);

if((max(Ha)<= delta_a) && (min(Hp) >= delta_p ))
    if provera==1,disp('Gabariti su zadovoljeni');end
else
    while(true)
        N=N+1;
        M=M+1;
        y = firpm(N,F,Hd);
        [H,w]=freqz(y,1,Nfreqz); 
        Hfa=abs(H);
        Ha = Hfa(ia:end);
        Hp = Hfa(1:ip);
        if((max(Ha) <= delta_a) && (min(Hp) >= delta_p ))
            break;
        end
    end

end

indeks = 0 : length(y)-1;

if(provera == 1) 
    disp('Potreban red filtra (N) je:');
    disp(N);
    figure
        subplot(221)
            stem(indeks,y),title('Impulsni odziv')
        subplot(222)
            plot(w,Hfa),title('Amplitudska karakteristika')
                xlabel('Relativna kruzna ucestanost ');
                ylabel('|H(z)|');
        subplot(223)
            plot(w,20*log10(Hfa)),title('Amplitudska karakteristika u dB')
                xlabel('Relativna kruzna ucestanost');
                ylabel('|H(z)| dB');
        subplot(224)
            plot(w,Hfp),title('Fazna karakteristika')
                xlabel('Relativna kruzna ucestanost');
                ylabel('|H(z)| dB');
end

end



