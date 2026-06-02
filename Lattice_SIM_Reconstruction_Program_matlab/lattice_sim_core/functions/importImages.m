function out = importImages(images)
%IMPORTIMAGES Apply the HiFi-SIM cosine border fade before deconvolution.

N = size(images, 3);
L = size(images, 1);
for idx = 1:N
    curImg = images(:, :, idx);
    if L > 256
        curImg = fadeBorderCos(curImg, 10);
    else
        curImg = fadeBorderCos(curImg, 0);
    end
    images(:, :, idx) = curImg;
end
out = images;
end

function out = fadeBorderCos(img, px)
[h, w] = size(img);
dat = img;
if px == 0
    out = dat;
    return;
end

fac = 1 / px * pi / 2;
for y = 1:px
    for x = 1:w
        dat(y, x) = dat(y, x) * power(sin((y - 1) * fac), 2);
    end
end
for y = h-px+1:h
    for x = 1:w
        dat(y, x) = dat(y, x) * power(sin((h - y) * fac), 2);
    end
end
for y = 1:h
    for x = 1:px
        dat(y, x) = dat(y, x) * power(sin((x - 1) * fac), 2);
    end
end
for y = 1:h
    for x = w-px+1:w
        dat(y, x) = dat(y, x) * power(sin((w - x) * fac), 2);
    end
end
out = dat;
end
