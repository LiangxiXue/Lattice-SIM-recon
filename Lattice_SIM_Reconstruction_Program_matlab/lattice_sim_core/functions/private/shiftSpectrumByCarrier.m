function shifted = shiftSpectrumByCarrier(component, carrierRadPerPixel)
%SHIFTSPECTRUMBYCARRIER Center a carrier component by spatial phase removal.

[h, w] = size(component);
[x, y] = meshgrid(0:w-1, 0:h-1);
shifted = component .* exp(-1i * (carrierRadPerPixel(1) * x + carrierRadPerPixel(2) * y));
end
