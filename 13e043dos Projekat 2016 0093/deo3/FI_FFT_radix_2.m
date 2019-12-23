function [ X ] = FI_FFT_radix_2(x, params)
% FI_FFT_radix_2 implementira brzu Furijeovu
% transofrmaciju iz prvog doma?eg zadatka koriš?enjem kona?ne preciznosti
% FFT_radix_2

    % Promenljiva koja se menja samo pri rekurzivnom pozivanju
    % Predstavlja stepen (stage) fft-a
    persistent stage
    if isempty(stage)
        stage = 0;
    end
    stage = stage + 1;
    
    
    FixedPointAttributes = fimath ( 'RoundingMethod', params.RoundingMethod,...
        'OverflowAction', params.OverflowAction,'ProductMode', 'SpecifyPrecision', 'ProductWordLength',...
        params.SIGNAL_BITLENGTH , 'ProductFractionLength', params.SIGNAL_FRAC ,'SumMode', 'SpecifyPrecision',...
        'SumWordLength', params.SIGNAL_BITLENGTH , 'SumFractionLength',params.SIGNAL_FRAC  ) ;  
    
    N0 = length(x); 
    X = fi( zeros(1,N0) , true , params.SIGNAL_BITLENGTH ,...
            params.SIGNAL_FRAC, FixedPointAttributes);   
    x = fi( x , true , params.SIGNAL_BITLENGTH ,...
            params.SIGNAL_FRAC, FixedPointAttributes);
    x0 = fi( x , true , params.SIGNAL_BITLENGTH ,...            %mozda bi bolje bilo sa wrap nadalje...
            params.SIGNAL_FRAC, FixedPointAttributes);   
    twidle = fi( exp(-2*pi*1j/N0) , true , params.SIGNAL_BITLENGTH ,...
            params.SIGNAL_FRAC, FixedPointAttributes);   
    X1 = fi( zeros(1,N0) , true , params.SIGNAL_BITLENGTH ,...
            params.SIGNAL_FRAC, FixedPointAttributes);   
  
        
%%

    if (length(x0) == 1)                         % Izlazak iz rekurzivnog pozivanja
        X = x0;
    else
        if(stage == 1)
            if( mod(log2(N0),1) ~= 0)         % proveramo da li je N stepen dvojke
                pow = floor(log2(N0)) + 1;          % p pomcna kojom cemo dopuniti niz 
                delta = 2^pow - N0;                % delta nam govori koliko da dopunimo
                x0 = [ x  zeros(1,delta)];          % dopunjen signal
                N0 = length(x0);                    % duzina produzenog signala
            end
        end

        X1 = zeros(1,N0);
        h = 1;                          % h pomaze pri neusaglasenosti 
        for it = 1:N0
            if (mod(it,2) == 1)
                X1(it) = x0(h) + x0(h + N0/2);
            else
                X1(it) = x0(h) - x0(h + N0/2); % * twidle ^ (nesto); to u sledecoj sekciji
                h = h + 1;
            end
        end


        %% Parne clanove mnozimo sa odgovarajucim twidle faktorima 
        p = 0;                                    % inicijalizujemo pomocnu
        twidle = exp(-2*pi*1j/N0);
        for k = 1:N0/2
            X1(2*k) = X1(2*k) * twidle ^ p;       % svaki paran mnozimo sa twidle-om
            if (mod(k,2^(stage - 1)) == 0 )                      
                p = p + 2^(stage - 1);                      
            end
        end
        %%
        if (2^stage ~= N0)
            X = FI_FFT_radix_2(X1,params);
        else
            % Preuredjivanje dobijenog signala
            Xtemp = X1;
            X(bitrevorder(1:N0)) = 1/N0* Xtemp(1:N0);
            stage = 0;
        end

    end
end
