function [ y ] = FI_IIR_direct( b, a, x )
% FI_IIR_direct_I Implementira direktnu I realizaciju IIR
% filtra ?iji su koeficijenti ulazni argumenti b i a i filtrira signal x. 
% Kao povratnu vrednost, funkcija vra?a filtrirani signal y.
% Svi koeficijenti i ulazni signal predstavljeni su kao brojevi sa fiksnom
% ta?kom, pa i sva izra?unavanja su brojevi sa fiksnom ta?kom. Tako?e, svi
% me?urezultati imaju preciznost definisanu u opcijama ulaznih argumenata 


    % Sirine izlaznog signala se podesavaju prema sirinama ulaznog signala
    % i prema sirinama koeficijenata (moze i drugacije, ovde je tako usvojeno)
fi_params = struct( 'SIGNAL_BITLENGTH', x.WordLength + b.WordLength, ...
                'SIGNAL_FRAC', x.FractionLength + b.FractionLength);

FixedPointAttributes = fimath( 'RoundingMethod', x.RoundingMethod, 'OverflowAction', 'Saturate', ...
        'ProductMode', 'SpecifyPrecision', 'ProductWordLength', fi_params.SIGNAL_BITLENGTH , 'ProductFractionLength', fi_params.SIGNAL_FRAC, ...
        'SumMode', 'SpecifyPrecision', 'SumWordLength', fi_params.SIGNAL_BITLENGTH, 'SumFractionLength', fi_params.SIGNAL_FRAC ) ;
    
 y = fi (zeros(1,length(x)) , true , fi_params.SIGNAL_BITLENGTH , fi_params.SIGNAL_FRAC, FixedPointAttributes);
 
 FixedPointAttributes.OverflowAction = 'Wrap';
     x = setfimath(x,FixedPointAttributes);
 FixedPointAttributes.RoundingMethod = b.RoundingMethod;
     b = setfimath(b,FixedPointAttributes);
     x1 = fi(zeros(1,1), true, fi_params.SIGNAL_BITLENGTH, fi_params.SIGNAL_FRAC, FixedPointAttributes);
 FixedPointAttributes.RoundingMethod = a.RoundingMethod;
     a = setfimath(a,FixedPointAttributes);
     y1 = fi(zeros(1,1), true, fi_params.SIGNAL_BITLENGTH, fi_params.SIGNAL_FRAC, FixedPointAttributes);
 
    Nx = length(x);
    Nb = length(b);
    Na = length(a);
    y  = zeros(1,Nx);   % Duzina izlazne skvence
    
    for i = 1:Nx   % i-ti element sekvence y
        x1 = 0;           % Pobuda
        y1 = 0;           % Povratna sprega

        for j = 0:Nb-1  
            if (i==j) 
                break;
            end  
            x1 = x1 + b(j+1)*x(i-j);
        end

        for j = 1:Na
             if (i==j) 
                 break;
             end  
             y1 = y1 - a(j)*y(i-j);
        end
         
        y(i) = x1+y1;
    end
    % Komentar sa vežbi:
        % Obratiti paznju da matricno mnozenje radi sumiranje na nama nepoznat
        % nacin. Ako zelimo da imamo bas model hardvera, ovde mora da se
        % raspise petlja u kojoj se sabira element po element onako kako se to
        % radi u hardveru.
    % ZBOG toga je uradjeno da bi imali model hardvera
    
end