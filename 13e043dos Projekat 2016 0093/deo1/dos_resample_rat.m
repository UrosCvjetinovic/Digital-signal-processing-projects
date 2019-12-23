function  [out] = dos_resample_rat(x, p, q) 
%dos_resample_rat menja ucestanost odabiranja ulaznog signala x, 
% za racionalni faktor p/q. 
%  Proceniti granicne ucestanosti filtra potrebnog za decimaciju ili 
% interpolaciju. Filtar se projektuje korišcenjem Parks-MekKlelanovog
% optimizacionog postupka. Slabljenje u nepropusnom opsegu treba da bude
% minimalno 70 dB, a varijacija amplitude u propusnom opsegu ne sme da 
% predje 0,05 dB.

if ((mod(p,1) ~= 0) || (mod(q,1)))
    error('Error: Funkcija dos_resample_rat prima cele brojeve'); 
end

if (( p <= 0) || (q <= 0))
    error('Error: Pogresno pozvan dos_resample_rat');
end

Ap = 0.05; Aa = 70;

pq = p;             % NF od interpolacije stroziji
if( q > p)
    pq = q;         % NF od decimacije stroziji
end

Wc = pi * 1/pq;
Bt = 0.4 * Wc;

h = lowpass_Parks_McClellan( Wc,Bt,Aa,Ap);   % NF filtar
Nh = length(h);
Nx = length(x);

%%---- Dodavanje nula
    yu = zeros(1, (p * Nx));
    for i = 0:Nx - 1
        yu(p*i + 1) = x(i+1);
    end
    
%%--- NF filtriranje
    yi = p * conv(yu,h); % NF filtriranje / Snaga se smanjuje p puta => 
                       % odabirke treba mnoziti sa sqrt(p). Zbog racuna u 
                       % ekvilizatoru bolje je bez korena
    delay = floor(Nh/2) + 1;
    yi = yi( delay : p*Nx+delay-1 );
    Nyi = length(yi);
    
%%--- Odbacivanje odabiraka
    z = conv(yi,h);
    delay = floor(Nh/2)+1; 
    z = z(delay : Nyi+delay-1);
    out = zeros( 1, floor(length(z)/q) );
    for i=0:length(out)-1
        out(i+1)=z(q*i+1);
    end
    
end

