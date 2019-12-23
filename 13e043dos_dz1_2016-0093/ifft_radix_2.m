function x = ifft_radix_2( X )
%ifft_radix_2 

if(mod(log2(length(X)),1) == 0) 
    x = real( 1 / length(X) *1j * conj(fft_radix_2(1j*conj(X))));
else
    error('Ulazni vektor nije duzine koja je stepen broja 2')
end

end
