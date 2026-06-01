function image = ifft2c(spectrum)
%IFFT2C Centered two-dimensional inverse FFT.

image = fftshift(ifft2(ifftshift(spectrum)));
end
