function F=fft2c(f)
F = ifftshift(fft2(fftshift(f)));
