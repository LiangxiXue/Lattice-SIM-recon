function spectrum = fft2c(image)
%FFT2C Centered two-dimensional FFT.

spectrum = fftshift(fft2(ifftshift(image)));
end
