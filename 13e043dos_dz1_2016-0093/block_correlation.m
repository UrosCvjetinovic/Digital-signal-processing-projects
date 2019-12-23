function R = block_correlation(x, pn_seq, block_length)
%block_correlation 

Ls = length(x);
h = pn_seq;
M = length(pn_seq);
L = block_length;

x0 = [ x zeros( 1, ceil(Ls/L) * L - Ls) ];
h0 = fliplr(h);
h0 = [ h0 zeros( 1, L - 1)];

%Umesto pravljenja matrice za sve clanove ci[n] koristimo dve promenljive
%koje pokazuju kako je izgledao c_{i-1}[n] i kako izgleda trenutni ci[n],
%njih koristimo da bi preklapanja sabrali u datim iteracijamatic
for i = 1 : ceil(Ls/L) 
             %blok xi[n]
    B =  [x0((i-1)*L + 1 : i*L) zeros(1, M-1)];
             %blok ci[n]= xi[n]*h[n]
    C = ifft_radix_2( fft_radix_2(B) .* fft_radix_2(h0));
    C = C(1:length(B));
    if( i == 1)
        % pocetni c1[n] pamtimo u promenljivu C_prev koja ce biti prethodni
        % u narednoj iteraciji a njega smo dobili preko korealacije izmedju
        % blokova L duzine i pseudo slucajne sekvence
        C_prev = C;
        c = C_prev(1:L);
    else
        C_cur = C;              % trunutni niz ci[n]
        c = [ c C_prev(L+1:L+M-1) + C_cur(1:M-1) C_cur(M:L)];
        %sadasnji ci[n] nam postaje prethodni za sledecu iteraciju
        C_prev = C_cur;                                                
    end
end
c = [c C_prev(L+1:end)];        % dodavanje poslednjeg clana

R = c;
end

