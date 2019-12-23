function [ y ] = IIR_direct_I( b,a,x )
% IIR_direct_I Implementira direktnu I realizaciju IIR
% filtra ?iji su koeficijenti ulazni argumenti b i a i filtrira signal x. 
% Kao povratnu vrednost, funkcija vra?a filtrirani signal y.

    Nx = length(x);
    Nb = length(b);
    Na = length(a);
    y  = zeros(1,Nx);   % Duzina izlazne skvence

    for i = 1:Nx          % i-ti element sekvence y
        x1 = 0;           % Pobuda
        y1 = 0;           % Povratna sprega

        for j = 0:Nb-1   % uticaj poude
            if (i == j) 
                break;
            end
            x1 = x1+b(j+1)*x(i-j);
        end

        for j = 1:Na    % uticaj povratne sprege
             if (i == j) 
                 break;
             end 
             y1 = y1-a(j)*y(i-j);
        end

        y(i)=x1+y1;
    end
end

