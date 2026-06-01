function writeTempTiffStack(stack, outputPath)
%WRITETEMPTIFFSTACK Write a numeric stack as a multi-page TIFF.

for idx = 1:size(stack, 3)
    if idx == 1
        imwrite(stack(:, :, idx), outputPath, 'tif');
    else
        imwrite(stack(:, :, idx), outputPath, 'tif', 'WriteMode', 'append');
    end
end
end
