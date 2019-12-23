function [ b, a ] = bandpass_filter( fs, fa, fp, Aa, Ap )
%bandpass_filter
%   projektuje IIR filtar propusnik opsega
% fs ucestanost odabiranja
% fa (dva elem) opseg nepropusnog
% fp (dva elem) granice propusnog
% Aa slabljenje u nepropusnom
% Ap slabljenje u propusnom
% a polinom u imeniocu
% b polinom u brojiocu
% uslov je da je fa(1) < fp(1) < fp(2) < fa(2)

% provera=1 --->crtaju se amplitudske karakteristike filtra
provera = 0;
% indvorsubp ----> crta grafike pojedinacno ili na subplotu
indvorsubp = 0;

if (( fa(1) > fp(1)) || (fp(1) > fp(2)) || (fp(2) > fa(2))...
        || (fa(1) > fa(2)) || (Ap > Aa)|| (Ap < 0) || (Aa<0))
    error('Pogresno pozvano');
end
%%
    Aag = Aa;           % Samo u slucaju da se da preveliko slabljenje i filtar bude nestabilan, smanjujemo slabljenje u propusnom
    Aatemp = Aa;
    Aptemp = Ap;
    meddled = 0;
    fp1 = fp(1);
    fp2 = fp(2);
    fa1 = fa(1);
    fa2 = fa(2);
    wp1 = 2 * pi * fp1;         % gr. ucestanost propusnog opsega
    wa1 = 2 * pi * fa1;         % gr. ucestanost nepropusnog opsega
    wa2 = 2 * pi * fa2;
    wp2 = 2 * pi * fp2;
    T= 1/fs;
    Nfreqz = 25000;
    
    
f0 = sqrt(fp1*fp2);
%---predistorzija ucestanosti za bilinearnu transformaciju
wa1Pred = 2/T*tan(wa1/fs/2);
wp1Pred = 2/T*tan(wp1/fs/2);
wp2Pred = 2/T*tan(wp2/fs/2);
wa2Pred = 2/T*tan(wa2/fs/2);

%---gabariti normalizovanog prototipa
W0 = sqrt(wp1Pred*wp2Pred);
B = wp2Pred-wp1Pred;

if (wa2Pred > W0^2/wa1Pred)
    wa2Pred = W0^2/wa1Pred;
elseif (wa1Pred < W0^2/wa2Pred)
    wa1Pred = W0^2/wa2Pred;
end
%---gabariti normalizovanog prototipa
wpN=1;
waN=1/B * (wa2Pred - wa1Pred);


%%
while(1)   
%---Odredjivanje reda NF filtra
        k=sqrt(1-(wpN/waN)^2);
        D=(10^(0.1*Aatemp)-1)/(10^(0.1*Aptemp)-1);
        q0=(1/2)*((1-sqrt(k))/(1+sqrt(k)));
        q=q0+2*q0^5+15*q0^9+15*q0^13;
        N_f=ceil(real(log10(16*D)/(log10(1/q))));
        if(N_f == 0) 
            N_f = 1;
        end
%---Projektovanje Analognog filtra
        [za,pa,ka]=ellipap(N_f,Aptemp,Aatemp);
        baN=ka*poly(za);
        aaN=poly(pa);
%---Transformacija normalizovanog prototipa u NO
        %-------------------NF -> PO
        [ba,aa]=lp2bp(baN,aaN,W0,B);
%---Diskretizacija
        [bd,ad]=bilinear(ba,aa,fs);

        [hd,Wd]=freqz(bd,ad,Nfreqz);
        Hd=abs(hd);
        fd = (Wd*fs)/(2*pi);
%%

%provera gabarita
    NO_OK = 0;  %nepropusni opseg ok (kada je NO_OK = 1)
    PO_OK = 0;  %propusni opseg ok (kada je PO_OK = 1)
    df=(fs/2)/Nfreqz;
    
    ia1=ceil(fa1/df)+1;
    ip1=floor(fp1/df)+1;
    ip2=ceil(fp2/df)+1;
    ia2=floor(fa2/df)+1;
    Ha=[Hd(1:ia1)' Hd(ia2:end)'];%Amplitudska Karakteristika Neprop Opsega
    Hp=Hd(ip1:ip2);   %Amplitudska Karakteristika Propusnog Opsega
    
    maxbs = max(20*log10(Ha));
    minbp = min(20*log10(Hp));
    
    if((maxbs > -Ap) || (Aptemp <= 0))    % u nepropusnom iznad -Ap
             % error('Suvise ostri gabariti');ILI SREDITI NA SLEDECI NACIN
             % MENJAMO ULAZNE GABARITE (u proveri crtamo i njih i pocetne)
             Aag = Aag - 0.002*Aa; 
             Aatemp = Aag;
             Aptemp = Ap;
             meddled = 1;
    end
    if (Aag < 0.5*Aa)
        % ogranicavamo se da slabljenje u neprop. opsegu
        % mozemo promeniti do pola pocetne vrednosti
        error('Previše strogi gabariti');
    end
    if(maxbs > -Aag)
        Aatemp = Aatemp + 0.002*Aa;
    else
        NO_OK=1;
    end
    if(minbp < -Ap)
        if(Aptemp > 0.01*Ap)
            Aptemp=Aptemp - 0.01*Ap;
        end
    else
        PO_OK=1;
    end   
    if((NO_OK==1)&&(PO_OK==1))
        if((maxbs < -Aa) && (minbp > -Ap))
            Aag = Aa;
            meddled = 0;
        end
        break
    end
end  
b = bd;
a = ad;
if (provera == 1)
%-----------------GRAFICKA PROVERA
            f = Wd/(2*pi)*fs;
            figure;
%---Amplitudska karakteristika digitalnog filtra u linearnoj razmeri
         if(indvorsubp == 1)
         else   subplot(3,1,1);
         end
            plot(f,Hd, 'LineWidth', 2);
            grid on
            title('Amplitudska karakteristika u linearnoj razmeri');
            xlabel('Ucestanost (Hz)');
            ylabel('|H(z)|');

        % crtanje DATIH gabarita
            delta_a = 10^(-0.05*Aa);
            delta_p = 1 - 10^(-0.05*Ap);
            delta_ag = 10^(-0.05*Aag);
            
        hold on
            x_1 = [0 fa1];      y_1 = [delta_a delta_a];
            x_2 = [fa1 fa1];	y_2 = [delta_a 1-delta_p];

            x_3 = [fa2 fs/2];   y_3 = [delta_a delta_a];
            x_4 = [fa2 fa2];    y_4 = [delta_a 1-delta_p];

            x_1_g = [fa1*1e-3 fa1];	y_1_g = [delta_ag delta_ag];
            x_2_g = [fa2 fs/2];	    y_2_g = [delta_ag delta_ag];
            x_5 = [fp1 fp2];	y_5 = [1-delta_p 1-delta_p];
            x_6 = [fp1 fp1];	y_6 = [delta_a 1-delta_p];
            x_7 = [fp2 fp2];	y_7 = [delta_a 1-delta_p];
            
        % crtanje POOSTRENIH gabarita
            delta_a = 10^(-0.05*Aatemp);
            delta_p = 1 - 10^(-0.05*Aptemp);
            
            x_1v = [0 fa1];      y_1v = [delta_a delta_a];
            x_2v = [fa1 fa1];	y_2v = [delta_a 1];

            x_3v = [fa2 fs/2];   y_3v = [delta_a delta_a];
            x_4v = [fa2 fa2];    y_4v = [delta_a 1];

            x_5v = [fp1 fp2];	y_5v = [1-delta_p 1-delta_p];
            x_6v = [fp1 fp1];	y_6v = [delta_a 1-delta_p];
            x_7v = [fp2 fp2];	y_7v = [delta_a 1-delta_p];
            
            % crvene - gabariti koji su prosledjeni funkciji
            % zelena - poostreni gabariti
            % magenta(ljubicasta) - slabljenje koje smo promenili da bi
            %               napravili filtar iako su previse ostri gabariti
            plot(x_1,y_1,'r',x_2,y_2,'r',x_3,y_3,'r',x_4,y_4,'r',...
                x_5,y_5,'r',x_6,y_6,'r',x_7,y_7,'r','LineWidth',1.5);
            plot(x_1v,y_1v,':g',x_2v,y_2v,':g',x_3v,y_3v,':g',...
                x_4v,y_4v,':g',x_5v,y_5v,':g',x_6v,y_6v,':g',...
                x_7v,y_7v,':g','LineWidth',1.5);
            plot(x_1_g,y_1_g,'m', x_2_g,y_2_g,'m');
        hold off
            
%Amplitudska karakteristika digitalnog filtra u logaritamskoj razmeri
         if(indvorsubp == 1)
            figure;
         else
            subplot(3,1,2);
         end
            plot(f,20*log10(Hd), 'LineWidth', 2);
            grid on
            title('Amplitudska karakteristika u logaritamskoj razmeri');
            xlabel('Ucestanost (Hz)');
            ylabel('|H(z)| dB');

        % crtanje DATIH gabarita
        hold on
            x_1 = [0 fa1];      y_1 = [-Aa -Aa];
            x_2 = [fa1 fa1];	y_2 = [-Aa -Ap];

            x_3 = [fa2 fs/2];   y_3 = [-Aa -Aa];
            x_4 = [fa2 fa2];    y_4 = [-Aa -Ap];

            x_1_g = [fa1*1e-3 fa1];	y_1_g = [-Aag -Aag];
            x_2_g = [fa2 fs/2];	    y_2_g = [-Aag -Aag];
            x_5 = [fp1 fp2];	y_5 = [-Ap -Ap];
            x_6 = [fp1 fp1];	y_6 = [-Aa -Ap];
            x_7 = [fp2 fp2];	y_7 = [-Aa -Ap];
            
        % crtanje POOSTRENIH gabarita
            x_1v = [0 fa1];    y_1v = [-Aatemp -Aatemp];
            x_2v = [fa1 fa1]; y_2v = [-Aatemp -Aptemp];

            x_3v = [fa2 fs/2]; y_3v = [-Aatemp -Aatemp];
            x_4v = [fa2 fa2]; y_4v = [-Aatemp -Aptemp];

            x_5v = [fp1 fp2]; y_5v = [-Aptemp -Aptemp];
            x_6v = [fp1 fp1]; y_6v = [-Aatemp -Aptemp];
            x_7v = [fp2 fp2]; y_7v = [-Aatemp -Aptemp];
            
            % crvene - gabariti koji su prosledjeni funkciji
            % zelena - poostreni gabariti
            % magenta(ljubicasta) - slabljenje koje smo promenili da bi
            %               napravili filtar iako su previse ostri gabariti
            plot(x_1,y_1,'r',x_2,y_2,'r',x_3,y_3,'r',x_4,y_4,'r',...
                x_5,y_5,'r',x_6,y_6,'r',x_7,y_7,'r','LineWidth',1.5);
            plot(x_1v,y_1v,':g',x_2v,y_2v,':g',x_3v,y_3v,':g',...
                x_4v,y_4v,':g',x_5v,y_5v,':g',x_6v,y_6v,':g',...
                x_7v,y_7v,':g','LineWidth',1.5);
            plot(x_1_g,y_1_g,'m', x_2_g,y_2_g,'m');
        hold off

%Amplitudska karakteristika pomocu semilogx
         if(indvorsubp == 1)
            figure;
         else
            subplot(3,1,3);
         end
            semilogx(f,20*log10(Hd), 'LineWidth', 2);
            grid on
            title('Amplitudska karakteristika semilogx');
            xlabel('Ucestanost (Hz)');
            ylabel('|H(z)| dB');

        % crtanje DATIH gabarita
        hold on
            x_1 = [fa1*1e-3 fa1];      y_1 = [-Aa -Aa];
            x_2 = [fa1 fa1];	y_2 = [-Aa -Ap];

            x_3 = [fa2 fs/2];   y_3 = [-Aa -Aa];
            x_4 = [fa2 fa2];    y_4 = [-Aa -Ap];

            x_1_g = [fa1*1e-3 fa1];	y_1_g = [-Aag -Aag];
            x_2_g = [fa2 fs/2];	    y_2_g = [-Aag -Aag];
            x_5 = [fp1 fp2];	y_5 = [-Ap -Ap];
            x_6 = [fp1 fp1];	y_6 = [-Aa -Ap];
            x_7 = [fp2 fp2];	y_7 = [-Aa -Ap];
            
            
        % crtanje POOSTRENIH gabarita
            x_1v = [fa1*1e-1 fa1];    y_1v = [-Aatemp -Aatemp];
            x_2v = [fa1 fa1];         y_2v = [-Aatemp -Aptemp];

            x_3v = [fa2 fs/2];         y_3v = [-Aatemp -Aatemp];
            x_4v = [fa2 fa2];        y_4v = [-Aatemp -Aptemp];

            x_5v = [fp1 fp2]; y_5v = [-Aptemp -Aptemp];
            x_6v = [fp1 fp1]; y_6v = [-Aatemp -Aptemp];
            x_7v = [fp2 fp2]; y_7v = [-Aatemp -Aptemp];
            
            % crvene - gabariti koji su prosledjeni funkciji
            % zelena - poostreni gabariti
            % magenta(ljubicasta) - slabljenje koje smo promenili da bi
            %               napravili filtar iako su previse ostri gabariti
            plot(x_1,y_1,'r',x_2,y_2,'r',x_3,y_3,'r',x_4,y_4,'r',...
                x_5,y_5,'r',x_6,y_6,'r',x_7,y_7,'r','LineWidth',1.5);
            plot(x_1v,y_1v,':g',x_2v,y_2v,':g',x_3v,y_3v,':g',...
                x_4v,y_4v,':g',x_5v,y_5v,':g',x_6v,y_6v,':g',...
                x_7v,y_7v,':g','LineWidth',1.5);
            plot(x_1_g,y_1_g,'m', x_2_g,y_2_g,'m');
        hold off
end

if (meddled == 1)      % Ispis u slucaju da smo korigovali nepropusni opseg
    str_out = sprintf('Gabariti su korigovani Aa=%d ->> %d',Aa,Aag);
    display(str_out)
end
end
