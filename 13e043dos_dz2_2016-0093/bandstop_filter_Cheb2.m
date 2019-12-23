function [ b, a ] = bandstop_filter_Cheb2( fs, fa, fp, Aa, Ap )
%bandstop_filter_Cheb2 
%   projektuje IIR filtar nepropusnik opsega ucestanosti koriscenjem 
%   inverznog Cebisejevog filtra za analogni prototip.
% fs ucestanost odabiranja
% fa (dva elem) opseg nepropusnog
% fp (dva elem) granice propusnog
% Aa slabljenje u nepropusnom
% Ap slabljenje u propusnom
% a polinom u imeniocu
% b polinom u brojiocu
% uslov je da je fp(1) < fa(1) < fa(2) < fp(2), Aa > Ap
 
% provera=1 --->crtaju se amplitudske karakteristike filtra
provera = 0;
% indvorsubp ----> crta grafike pojedinacno ili na subplotu
indvorsubp = 0;    

if ((fp(1) > fa(1)) || (fa(1) > fa(2)) || (fa(2) > fp(2)) ...
        || (fp(1) >fp(2)) || (Ap > Aa) || (Ap < 0) || (Aa<0))
    error('Pogresno pozvana funkcija');
end
%%           
    % Samo u slucaju da se zada preveliko slabljenje i filtar bude vec u 
    % prvoj iteracijji nestabilan, smanjujemo slabljenje u nepropusnom
    Aag = Aa;
    % Promenljiva meddled ima vrednost 1 ako smo korigovali zadate gabarite
    meddled = 0;
    Aatemp = Aag;
    Aptemp = Ap;
    fp1 = fp(1);
    fp2 = fp(2);
    fa1 = fa(1);
    fa2 = fa(2);
    wp1 = 2 * pi * fp1;         % gr. kruzna ucestanost propusnog opsega
    wa1 = 2 * pi * fa1;         % gr. kruzna ucestanost nepropusnog opsega
    wa2 = 2 * pi * fa2;
    wp2 = 2 * pi * fp2;
    T= 1/fs;
    Nfreqz = 25000;
    
    
%predistorzija ucestanosti za bilinearnu transformaciju
    wp1Pred = 2/T*tan(wp1/fs/2);
    wa1Pred = 2/T*tan(wa1/fs/2);
    wa2Pred = 2/T*tan(wa2/fs/2);
    wp2Pred = 2/T*tan(wp2/fs/2);

W0 = sqrt(wa1Pred*wa2Pred);
B = wa2Pred - wa1Pred; 
if (wp1Pred < W0 ^ 2 / wp2Pred)
    wp1Pred = W0 ^ 2 / wp2Pred;
end
%gabariti normalizovanog prototipa
    waN = 1;
    wpN = (B * wp1Pred) / (W0^2 - wp1Pred^2);

%%
while(1)   
%---Odredjivanje reda NF filtra
        % faktor selektivnosti
        k = wpN / waN;
        % faktor diskriminacije 
        D = (10^(0.1*Aatemp)-1)/(10^(0.1*Aptemp)-1) ;
        % red filtra
        N_f = ceil(real(acosh(sqrt(D))/acosh(1/k)));
        % jer ceil(0) = 0, ovim izbegavamo zaokruzivanje na 0
        if(N_f == 0)
            N_f = 1;
        end
%---Kreiranje analognog NF prototipa
        [za,pa,ka]=cheb2ap(N_f,Aatemp);
        baN= ka * poly(za);
        aaN= poly(pa);
        
%---Transformacija normalizovanog prototipa u NO
        %---NF -> NO
        [ba,aa]=lp2bs(baN,aaN,W0,B);
%---Diskretizacija preko bilinearne transformacije
        [bd,ad]=bilinear(ba,aa,fs);

        [hd,Wd]=freqz(bd,ad,Nfreqz);
        Hd=abs(hd);
        fd = (Wd*fs)/(2*pi);
%%

%---provera gabarita
    NO_OK = 0;  %nepropusni opseg ok (kada je NO_OK = 1)
    PO_OK = 0;  %propusni opseg ok (kada je PO_OK = 1)
    df=(fs/2)/Nfreqz;

    ip1=ceil(fp1/df)+1;
    ia1=floor(fa1/df)+1;
    ia2=ceil(fa2/df)+1;
    ip2=floor(fp2/df)+1;
    Ha=Hd(ia1:ia2);              %Amplitudska Karakteristika Neprop. Opsega
    Hp=[Hd(1:ip1)',Hd(ip2:end)'];%Amplitudska Karakteristika Prop. Opsega
    
    maxbs = max(20*log10(Ha));
    minbp = min(20*log10(Hp));
    if((maxbs > -Ap) || (Aptemp <= 0))    % u nepropusnom iznad -Ap
               %error('Suvise ostri gabariti');ILI SREDITI NA SLEDECI NACIN
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
    if(maxbs > -Aag)
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
    if((NO_OK == 1)&&(PO_OK == 1))
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
% Amplitudska karakteristika digitalnog filtra u linearnoj razmeri
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
            delta_ag = 10^(-0.05*Aag);
            delta_a = 10^(-0.05*Aa);
            delta_p = 1 - 10^(-0.05*Ap);
        hold on
            x_1 = [ 0  fp1];    y_1 = [1-delta_p 1-delta_p];
            x_2 = [fp1 fp1];	y_2 = [1-delta_p delta_a];

            x_3 = [fp2 fs/2];	y_3 = [1-delta_p 1-delta_p];
            x_4 = [fp2 fp2];	y_4 = [1-delta_p delta_a];

            x_5 = [fa1 fa2];	y_5 = [delta_a delta_a];    
            x_6 = [fa1 fa1];	y_6 = [1 delta_a];
            x_7 = [fa2 fa2];	y_7 = [1 delta_a];
            % U slucaju da su previse strogi gabariti dati, 
            % ovo su granice koje je ova funkcija uspela da ispuni
            x_5_g = [fa1 fa2];	y_5_g = [delta_ag delta_ag];
            
            
        % crtanje POOSTRENIH gabarita
            delta_a = 10^(-0.05*Aatemp);
            delta_p = 1 - 10^(-0.05*Aptemp);
            
            x_1v = [  0  fp1]; y_1v = [1-delta_p 1-delta_p];
            x_2v = [fp1 fp1];	y_2v = [1-delta_p delta_a];

            x_3v = [fp2 fs/2]; y_3v = [1-delta_p 1-delta_p];
            x_4v = [fp2 fp2]; y_4v = [1-delta_p delta_a];

            x_5v = [fa1 fa2]; y_5v = [delta_a delta_a];
            x_6v = [fa1 fa1]; y_6v = [1 delta_a];
            x_7v = [fa2 fa2]; y_7v = [1 delta_a];
            
            % crvene - gabariti koji su prosledjeni funkciji
            % zelena - poostreni gabariti
            % magenta(ljubicasta) - slabljenje koje smo promenili da bi
            %               napravili filtar iako su previse ostri gabariti
            plot(x_1,y_1,'r',x_2,y_2,'r',x_3,y_3,'r',x_4,y_4,'r',...
                x_5,y_5,'r',x_6,y_6,'r',x_7,y_7,'r','LineWidth',1.5);
            plot(x_1v,y_1v,':g',x_2v,y_2v,':g',x_3v,y_3v,':g',...
                x_4v,y_4v,':g',x_5v,y_5v,':g',x_6v,y_6v,':g',...
                x_7v,y_7v,':g','LineWidth',1.5);
            plot(x_5_g,y_5_g,'m');
           
        hold off
            
            
% Amplitudska karakteristika digitalnog filtra logaritamskoj razmeri
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
            x_1 = [ 0  fp1];    y_1 = [-Ap -Ap];
            x_2 = [fp1 fp1];    y_2 = [-Ap -Aa];

            x_3 = [fp2 fs/2];   y_3 = [-Ap -Ap];
            x_4 = [fp2 fp2];    y_4 = [-Ap -Aa];

            x_5_g = [fa1 fa2];	y_5_g = [-Aag -Aag];
            x_5 = [fa1 fa2];    y_5 = [-Aa -Aa];
            x_6 = [fa1 fa1];    y_6 = [0 -Aa];
            x_7 = [fa2 fa2];    y_7 = [0 -Aa];
            
        % crtanje POOSTRENIH gabarita
            x_1v = [  0  fp1];     y_1v = [-Aptemp -Aptemp];
            x_2v = [fp1 fp1];     y_2v = [-Aptemp -Aatemp];

            x_3v = [fp2 fs/2];     y_3v = [-Aptemp -Aptemp];
            x_4v = [fp2 fp2];     y_4v = [-Aptemp -Aatemp];

            x_5v = [fa1 fa2];     y_5v = [-Aatemp -Aatemp];
            x_6v = [fa1 fa1];     y_6v = [0 -Aatemp];
            x_7v = [fa2 fa2];     y_7v = [0 -Aatemp];
            
            
            % crvene - gabariti koji su prosledjeni funkciji
            % zelena - poostreni gabariti
            % magenta(ljubicasta) - slabljenje koje smo promenili da bi
            %               napravili filtar iako su previse ostri gabariti
            plot(x_1,y_1,'r',x_2,y_2,'r',x_3,y_3,'r',x_4,y_4,'r',...
                x_5,y_5,'r',x_6,y_6,'r',x_7,y_7,'r','LineWidth',1.5);
            plot(x_1v,y_1v,':g',x_2v,y_2v,':g',x_3v,y_3v,':g',...
                x_4v,y_4v,':g',x_5v,y_5v,':g',x_6v,y_6v,':g',...
                x_7v,y_7v,':g','LineWidth',1.5);
            plot(x_5_g,y_5_g,'m');
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

        %crtanje DATIH gabarita
        hold on
            x_1 = [fp1*1e-1  fp1];    y_1 = [-Ap -Ap];
            x_2 = [fp1 fp1];	y_2 = [-Ap -Aa];

            x_3 = [fp2 fs/2];	y_3 = [-Ap -Ap];
            x_4 = [fp2 fp2];	y_4 = [-Ap -Aa];

            x_5_g = [fa1 fa2];	y_5_g = [-Aag -Aag];
            x_5 = [fa1 fa2];	y_5 = [-Aa -Aa];
            x_6 = [fa1 fa1];	y_6 = [0 -Aa];
            x_7 = [fa2 fa2];	y_7 = [0 -Aa];
            
        % crtanje POOSTRENIH gabarita
            
            x_1v = [fp1*1e-1  fp1];     y_1v = [-Aptemp -Aptemp];
            x_2v = [fp1 fp1];     y_2v = [-Aptemp -Aatemp];

            x_3v = [fp2 fs/2];     y_3v = [-Aptemp -Aptemp];
            x_4v = [fp2 fp2];     y_4v = [-Aptemp -Aatemp];

            x_5v = [fa1 fa2];     y_5v = [-Aatemp -Aatemp];
            x_6v = [fa1 fa1];     y_6v = [0 -Aatemp];
            x_7v = [fa2 fa2];     y_7v = [0 -Aatemp];
            
            
            % crvene - gabariti koji su prosledjeni funkciji
            % zelena - poostreni gabariti
            % magenta(ljubicasta) - slabljenje koje smo promenili da bi
            %               napravili filtar iako su previse ostri gabariti
            plot(x_1,y_1,'r',x_2,y_2,'r',x_3,y_3,'r',x_4,y_4,'r',...
                x_5,y_5,'r',x_6,y_6,'r',x_7,y_7,'r','LineWidth',1.5);
            plot(x_1v,y_1v,':g',x_2v,y_2v,':g',x_3v,y_3v,':g',...
                x_4v,y_4v,':g',x_5v,y_5v,':g',x_6v,y_6v,':g',...
                x_7v,y_7v,':g','LineWidth',1.5);
            plot(x_5_g,y_5_g,'m');
       hold off
end

if (meddled == 1)      % Ispis u slucaju da smo korigovali nepropusni opseg
    str_out = sprintf('Gabariti su korigovani Aa=%d ->> %d',Aa,Aag);
    display(str_out)
end

end
