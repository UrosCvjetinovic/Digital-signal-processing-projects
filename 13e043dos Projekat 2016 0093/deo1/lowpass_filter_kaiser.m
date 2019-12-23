function [h] = lowpass_filter_kaiser(Wc, Bt, Aa, Ap) 
%lowpass_filter_kaiser projektovanje FIR filtra propusnika niskih ucestanosti
%   Funkcija kao argumente prima idealnu grani?nu u?estanost Wc, 
% širinu prelazne zone Bt i odgovaraju?a slabljenja u nepropusnom 
% (Aa) i propusnom (Ap) opsegu. Kao povratnu vrednost, funkcija 
% vra?a impulsni odziv FIR filtra projektovanog odsecanjem idealnog
% impulsnog odziva Kajzerovim prozorom. Ukoliko gabariti digitalnog 
% filtra nisu zadovoljeni, funkcija treba da podešava parametre 
% analognog prototipa dok oni ne postanu zadovoljeni. Ra?unati
% frekvencijske karakteristike u dovoljno velikom broju ta?aka 
% prilikom provere, npr. više od 10000. 2

% provera --->crtaju se amplitudske karakteristike filtra
provera = 0;    
  % indvorsubp ----> crta grafike pojedinacno ili na subplotu
    indvorsubp = 0;

if ((Wc > pi) || (Bt > pi) || (Wc < 0) || (Bt < 0))
    error('Error: Pogresno prosledjeni podaci za FIR filtar');
end
if (Wc - Bt/2 < 0) 
    error('Error: Prevelik propusni opseg u odnosu na granicnu vrednost');
end
if (Ap > Aa)
    error('Error: Slabljenje u propusnom vece nego u nepropusnom')
end

Wp = Wc - Bt/2;
Wa = Wc + Bt/2;

Ap_g = Ap;
Aa_g = Aa;
Aptemp = Ap;
Aatemp = Aa;

delta_a = 10^(-0.05*Aa);
delta_p = (10^(0.05*Ap) - 1) / (10^(0.05*Ap) + 1);

meddled = 0;

while(1)
    delta_atemp = 10^(-0.05*Aatemp);
    delta_ptemp = (10^(0.05*Aptemp) - 1) / (10^(0.05*Aptemp) + 1);

    %---odredjivanje delta
    delta = min(delta_ptemp, delta_atemp);
    if (delta ~= delta_atemp);
        Aatemp = -20*log10(delta);
    end

    beta = 0;
    %---odredjivanje beta
    if((Aatemp >= 21) && (Aatemp <=50))
        beta = 0.5842*(Aatemp - 21) ^ 0.4 + 0.07886*(Aatemp - 21);
    elseif (Aatemp > 50)
        beta = 0.1102*(Aatemp - 8.7);
    end

    %---odredjivanje M
    D = 0.9222;
    if (Aatemp > 21)
        D=(Aatemp-7.95)/14.36; 
    end
    M = ceil(2*pi*D/Bt+1); % filter length / duzina impulsnog odziva
    N = M-1;               % filter order / red filtra
    
    %---generisanje kajzerovog prozora
    window = kaiser(M,beta)';

    %---odredjivanje impulsnog odziva idealnog NF filtra
    n = -(M-1)/2:(M-1)/2;  
    h = sin(n*Wc)./(n*pi);  % ifft od idealnog NF
    if ( mod(M,2) == 1 )
        indeks = (M+1)/2 ;
        h(indeks) = Wc/pi;
    end

    %---mnozenje sa odbircima prozorske funkcije
    h=h.*window;
    
    %---Kreiranje digitalnog
    Nfreqz = 2^16;    % > 10 000
    H = fft(h,Nfreqz);
    Hfa = abs(H(1:floor(Nfreqz/2))); 
%     Hfp = unwrap(angle(H));

    %---------provera gabarita
        NO_OK = 0;  %nepropusni opseg ok (kada je NO_OK = 1)
        PO_OK = 0;  %propusni opseg ok (kada je PO_OK = 1)
        dw = (2*pi)/Nfreqz;
        ia = floor(Wa/dw); 
        ip = ceil(Wp/dw);
        
        Ha = Hfa(ia:floor(Nfreqz/2));%Amplitudska Karakteristika Nepropusnog Opsega
        Hp = Hfa(1:ip);             %Amplitudska Karakteristika Propusnog Opsega
        
        maxbs = max(20*log10(Ha));
        minbp = min(20*log10(Hp));
        
        if((maxbs > -Ap) || (Aptemp <= 0))    
                    % u nepropusnom iznad -Ap
                   %error('Suvise ostri gabariti');ILI SREDITI NA SLEDECI NACIN
                 % MENJAMO ULAZNE GABARITE (u proveri crtamo i njih i pocetne)
                 Aa_g = Aa_g - 0.01*Aa; 
                 Aatemp = Aa_g;
                 Aptemp = Ap;
                 meddled = 1;
        end
        if (Aa_g < 0.5*Aa)
            % ogranicavamo se da slabljenje u neprop. opsegu
            % mozemo promeniti do pola pocetne vrednosti
            error('Previše strogi gabariti');
        end
        if(maxbs > -Aa_g)
            Aatemp = Aatemp + 0.002*Aa;            
        else
            NO_OK=1;
        end
        if(minbp < -Ap) 
            if(Aptemp > 0.01*Ap)
                Aptemp = Aptemp - 0.01*Ap;
            end
           else
                PO_OK=1;
        end
        if((NO_OK==1)&&(PO_OK==1))
            break
        end
end    
    
if (provera == 1)
%-----------------GRAFICKA PROVERA
        w = (0:length(Hfa)-1)*2*pi/Nfreqz;
        Wend = pi;
        delta_a_g = 10^(-0.05*Aa_g);
        delta_p_g = (10^(0.05*Ap_g) - 1) / (10^(0.05*Ap_g) + 1);
%---Amplitudska karakteristika digitalnog filtra u linearnoj razmeri
            figure;
         if(indvorsubp == 1)
         else
            subplot(3,1,1);
         end
            plot(w,Hfa, 'LineWidth', 2);
            grid on
            title('Amplitudska karakteristika u linearnoj razmeri');
            xlabel('Kruzna ucestanost (Hz)');
            ylabel('|H(z)|');

        % crtanje DATIH gabarita\
        hold on
            x_1_g = [Wa Wend];    y_1_g = [delta_a_g delta_a_g];
            x_1 = [Wa Wend];      y_1 = [delta_a delta_a];
            x_2 = [Wa Wa];	y_2 = [delta_a 1-delta_p];
            
            x_3 = [0  Wp];	y_3 = [1-delta_p 1-delta_p];
            x_4 = [Wp Wp];	y_4= [delta_a 1-delta_p];
            
        % crtanje POOSTRENIH gabarita
            delta_a = 10^(-0.05*Aatemp);
            delta_p = 1 - 10^(-0.05*Aptemp);
            
            x_1v = [Wa Wend];      y_1v = [delta_a delta_a];
            x_2v = [Wa Wa];	y_2v = [delta_a 10^(-0.05*Aa)];
            
            x_3v = [0  Wp];	y_3v = [1-delta_p 1-delta_p];
            x_4v = [Wp Wp];	y_4v= [10^(-0.05*Ap) 1-delta_p];
            
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
            plot(w,20*log10(Hfa), 'LineWidth', 2);
            grid on
            title('Amplitudska karakteristika u logaritamskoj razmeri');
            xlabel('Kruzna ucestanost (Hz)');
            ylabel('|H(z)| dB');

        hold on
        % crtanje DATIH gabarita

            x_1_g = [Wa Wend];    y_1_g = [-Aa_g -Aa_g];
            x_1 = [Wa Wend];    y_1 = [-Aa -Aa];
            x_2 = [Wa Wa];      y_2 = [-Aa -Ap];

            x_3 = [ 0 Wp];      y_3 = [-Ap -Ap];
            x_4 = [Wp Wp];      y_4 = [-Aa -Ap];
        % crtanje POOSTRENIH gabarita
            x_1v = [Wa Wend];       y_1v = [-Aatemp -Aatemp];
            x_2v = [Wa Wa];      y_2v = [-Aa -Aatemp];

            x_3v = [ 0 Wp]; 	y_3v = [-Aptemp -Aptemp];
            x_4v = [Wp Wp];      y_4v = [-Ap -Aptemp];
            
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
            semilogx(w,20*log10(Hfa), 'LineWidth', 2);
            grid on
            title('Amplitudska karakteristika semilogx');
            xlabel('Kruzna ucestanost (Hz)');
            ylabel('|H(z)| dB');

            %crtanje DATIH gabarita
        hold on

            x_1_g = [Wa Wend];        y_1_g = [-Aa_g -Aa_g];
            x_1 = [Wa Wend];          y_1 = [-Aa -Aa];
            x_2 = [Wa Wa];            y_2 = [-Aa -Ap];
    
            x_3 = [Wp*1e-1 Wp];        y_3 = [-Ap -Ap];
            x_4 = [Wp Wp];            y_4 = [-Aa -Ap];
        % crtanje POOSTRENIH gabarita
            x_1v = [Wa Wend];          y_1v = [-Aatemp -Aatemp];
            x_2v = [Wa Wa];            y_2v = [-Aatemp -Aa];

            x_3v = [Wp*1e-1 Wp];        y_3v = [-Aptemp -Aptemp];
            x_4v = [Wp Wp];             y_4v = [-Ap -Aptemp];
            
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
    str_out = sprintf('Gabariti su korigovani Aa=%d ->> %d',Aa,Aa_g);
    display(str_out)
end

end

