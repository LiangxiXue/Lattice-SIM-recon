function out = FFT2D(in, inverse)
%FFT2D HiFi-SIM compatible centered two-dimensional FFT helper.

if inverse == true
    out = ifft2(fftshift(in));
else
    out = fftshift(fft2(in));
end
end
