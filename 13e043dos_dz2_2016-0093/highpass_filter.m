function [ b, a ] = highpass_filter( fs, fa, fp, Aa, Ap )
%highpass_filter
%   projektuje IIR filtar propusnik visokih ucestanosti
% fs ucestanost odabiranja
% fa opseg nepropusnog
% fp granica propusnog
% Aa slabljenje u nepropusnom
% Ap slabljenje u propusnom
% a polinom u imeniocu
% b polinom u brojiocu
% uslov je da je fa(1) < fp(1)

% provera=1 --->crtaju se amplitudske karakteristike filtra
provera = 0;    
% indvorsubp ----> crta grafike pojedinacno ili na subplotu
indvorsubp = 0;

if (( fa > fp) || (Ap > Aa)|| (Ap < 0) || (Aa<0))
    error('Pogresno pozvano');
end
%%
    Aag = Aa;
    Aatemp = Aa;
    Aptemp = Ap;
    wp = 2 * pi * fp;         % gr. ucestanost propusnog opsega
    wa = 2 * pi * fa;         % gr. ucestanost nepropusnog opsega
    T= 1/fs;
    Nfreqz = 25000;
    meddled = 0;
%---predistorzija ucestanosti za bilinearnu transformaciju
wpPred = 2/T*tan(wp/fs/2);
waPred = 2/T*tan(wa/fs/2);
%---gabariti normalizovanog prototipa
 wpN = 1;
 a = wpN * wpPred;
 waN = a / waPred;

%%
while(1)   
% ---Odredjivanje reda NF filtra
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
        %Transformacija normalizovanog prototipa u NO
%---NF -> VF
        [ba,aa]=lp2hp(baN,aaN,wpPred);
%---Diskretizacija
        [bd,ad]=bilinear(ba,aa,fs);

        [hd,Wd]=freqz(bd,ad,Nfreqz);
        Hd=abs(hd);
        fd = (Wd*fs)/(2*pi);
%%

%---provera gabarita
    NO_OK = 0;  %nepropusni opseg ok (kada je NO_OK = 1)
    PO_OK = 0;  %propusni opseg ok (kada je PO_OK = 1)
    df=(fs/2)/Nfreqz;

    ia=ceil(fa/df)+1;
    ip=floor(fp/df)+1;
    Ha=Hd(1:ia);               %Amplitudska Karakteristika Nepropusnog Opsega
    Hp=Hd(ip:end);             %Amplitudska Karakteristika Propusnog Opsega
    
    maxbs = max(20*log10(Ha));
    minbp = min(20*log10(Hp));
    
    if((maxbs > -Ap) || (Aptemp <= 0))    % u nepropusnom iznad -Ap
             % error('Suvise ostri gabariti');ILI SREDITI NA SLEDECI NACIN
             % MENJAMO ULAZNE GABARITE (u proveri crtamo i njih i pocetne)
             Aag = Aag - 0.01*Aa; 
             Aatemp = Aag;
             Aptemp = Ap;
             meddled = 1;
    end
    if (Aag < 0.5*Aa)
        % ogranicavamo se da slabljenje u neprop. opsegu
        % mozemo promeniti do pola pocetne vrednosti
        error('Previše strogi gabariti');
    end
    if(maxbs >-Aag)
        Aatemp=Aatemp + 0.002*Aa;
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
%---Amplitudska karakteristika digitalnog filtra u linearnoj razmeri
            figure;
         if(indvorsubp == 1)
         else
            subplot(3,1,1);
         end
            plot(f,Hd, 'LineWidth', 2);
            grid on
            title('Amplitudska karakteristika u linearnoj razmeri');
            xlabel('Ucestanost (Hz)');
            ylabel('|H(z)|');

        % crtanje DATIH gabarita
            delta_a = 10^(-0.05*Aa);
            delta_ag = 10^(-0.05*Aag);
            delta_p = 1 - 10^(-0.05*Ap);
        hold on
            x_1_g = [0 fa];    y_1_g = [delta_ag delta_ag];
            x_1 = [0 fa];      y_1 = [delta_a delta_a];
            x_2 = [fa fa];	y_2 = [delta_a 1-delta_p];
            
            x_3 = [fp fs/2];	y_3 = [1-delta_p 1-delta_p];
            x_4 = [fp fp];	y_4= [delta_a 1-delta_p];
            
        % crtanje POOSTRENIH gabarita
            delta_a = 10^(-0.05*Aatemp);
            delta_p = 1 - 10^(-0.05*Aptemp);
            
            x_1v = [0 fa];      y_1v = [delta_a delta_a];
            x_2v = [fa fa];	y_2v = [delta_a 10^(-0.05*Aa)];
            
            x_3v = [fp fs/2];	y_3v = [1-delta_p 1-delta_p];
            x_4v = [fp fp];	y_4v= [10^(-0.05*Ap) 1-delta_p];
            
            % crvene - gabariti koji su prosledjeni funkciji
            % zelena - poostreni gabariti
            % magenta(ljubicasta) - slabljenje koje smo promenili da bi
            %               napravili filtar iako su previse ostri gabariti
            plot(x_1,y_1,'r',x_2,y_2,'r',x_3,y_3,'r',x_4,y_4,...
                'r','LineWidth',1.5);
            plot(x_1v,y_1v,':g',x_2v,y_2v,':g',x_3v,y_3v,':g',...
                x_4v,y_4v,'LineWidth',1.5);
            plot(x_1_g,y_1_g,'m');
        hold off
            
%---Amplitudska karakteristika digitalnog filtra u logaritamskoj razmeri
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

        hold on
        % crtanje DATIH gabarita
            x_1_g = [0 fa];    y_1_g = [-Aag -Aag];
            x_1 = [0 fa];       y_1 = [-Aa -Aa];
            x_2 = [fa fa];      y_2 = [-Aa -Ap];

            x_3 = [fp fs/2];	y_3 = [-Ap -Ap];
            x_4 = [fp fp];      y_4 = [-Aa -Ap];
        % crtanje POOSTRENIH gabarita
            x_1v = [0 fa];       y_1v = [-Aatemp -Aatemp];
            x_2v = [fa fa];      y_2v = [-Aa -Aatemp];

            x_3v = [fp fs/2];	y_3v = [-Aptemp -Aptemp];
            x_4v = [fp fp];      y_4v = [-Ap -Aptemp];
            
           
            % crvene - gabariti koji su prosledjeni funkciji
            % zelena - poostreni gabariti
            % magenta(ljubicasta) - slabljenje koje smo promenili da bi
            %               napravili filtar iako su previse ostri gabariti
            plot(x_1,y_1,'r',x_2,y_2,'r',x_3,y_3,'r',x_4,y_4,...
                'r','LineWidth',1.5);
            plot(x_1v,y_1v,':g',x_2v,y_2v,':g',x_3v,y_3v,':g',...
                x_4v,y_4v,'LineWidth',1.5);
            plot(x_1_g,y_1_g,'m');
        hold off
            
%---Amplitudska karakteristika pomocu semilogx
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

        %crtanje DATIH gabarita
        hold on

            x_1_g = [fa*1e-2 fa];    y_1_g = [-Aag -Aag];
            x_1 = [fa*1e-2 fa];     y_1 = [-Aa -Aa];
            x_2 = [fa fa];          y_2 = [-Aa -Ap];

            x_3 = [fp fs/2];        y_3 = [-Ap -Ap];
            x_4 = [fp fp];          y_4 = [-Aa -Ap];
        % crtanje POOSTRENIH gabarita
            x_1v = [fa*1e-2 fa];    y_1v = [-Aatemp -Aatemp];
            x_2v = [fa fa];         y_2v = [-Aatemp -Aa];

            x_3v = [fp fs/2];       y_3v = [-Aptemp -Aptemp];
            x_4v = [fp fp];         y_4v = [-Ap -Aptemp];
            
           
            % crvene - gabariti koji su prosledjeni funkciji
            % zelena - poostreni gabariti
            % magenta(ljubicasta) - slabljenje koje smo promenili da bi
            %               napravili filtar iako su previse ostri gabariti
            plot(x_1,y_1,'r',x_2,y_2,'r',x_3,y_3,'r',x_4,y_4,...
                'r','LineWidth',1.5);
            plot(x_1v,y_1v,':g',x_2v,y_2v,':g',x_3v,y_3v,':g',...
                x_4v,y_4v,'LineWidth',1.5);
            plot(x_1_g,y_1_g,'m');
        hold off

end

if (meddled == 1)      % Ispis u slucaju da smo korigovali nepropusni opseg
    str_out = sprintf('Gabariti su korigovani Aa=%d ->> %d',Aa,Aag);
    display(str_out)
end

end
