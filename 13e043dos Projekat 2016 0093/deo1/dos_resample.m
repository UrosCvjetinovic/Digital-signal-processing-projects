function [y] = dos_resample(x, up_down, r)
% dos_resample  menja ucestanost odabiranja ulaznog signala x, za celobrojni
% faktor r. 
% Ucestanost odbiranja se povecava, ako je parametar up_down na vrednosti 
% false, a smanjuje ako je na vrednosti true. Filtar potreban za decimaciju 
% ili interpolaciju ima idealnu grani?nu u?estanost WC = pi/r i prelaznu zonu 
% širine Bt = 0,4pi/r. Filtar se projektuje koriš?enjem funckije iz prethodne
% tacke. Slabljenje u nepropusnom opsegu treba da bude minimalno 70 dB, a
% varijacija amplitude u propusnom opsegu ne sme da pre?e 0,05 dB. Primetiti
% da impulsni odzivi decimacionih/interpolacionih filtara za ovaj slu?aj 
% imaju vrednost 0 na indeksima koji su umnožak faktora r.

if (mod(r,1) ~= 0) 
    error('Error: Funkcija dos_resample prima cele brojeve'); 
end
if ( r == 1)
    y = x;
else
    Wc = pi/r;
    Bt = 0.4 * Wc;
    Aa = 70;
    Ap = 0.05;

    h = lowpass_filter_kaiser(Wc,Bt,Aa,Ap);     % NF filtar
    Nh = length(h);
    Nx = length(x);

    if ( strcmp(up_down, 'false'))   % Interpolacija
        y = zeros(1, (r * Nx));
        for i = 0:Nx - 1
            y(r*i + 1) = x(i+1);
        end
        y = r * conv(y,h); % NF filtriranje / Snaga se smanjuje r puta => 
                           % odabirke treba mnoziti sa sqrt(r). Zbog racuna u 
                           % ekvilizatoru bolje je bez korena
        delay = floor(Nh/2) + 1;
        y = y( delay : r*Nx+delay-1 );

    elseif ( strcmp(up_down, 'true'))  % Decimacija
        z = conv(x,h);              % NF filtriranje
        delay = floor(Nh/2)+1;
        z = z(delay : Nx+delay-1);
        y = zeros( 1, floor(length(z)/r) );
        for i=0:length(y)-1
            y(i+1)=z(r*i+1);
        end
    else 
        error('Error: Invalid up_down');
    end
end
end