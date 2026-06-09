clc;
clear;

inputPath  = 'C9.tif';
outputPath = 'C9WF.tif';

% 读取后三帧
img3 = imread(inputPath, 3);
img4 = imread(inputPath, 4);
img5 = imread(inputPath, 5);

% 平均
avgImg = (double(img3) + double(img4) + double(img5)) / 3;

% 转回原图类型
avgImg = cast(round(avgImg), class(img3));

% 输出
imwrite(avgImg, outputPath, 'tif');

disp('后三帧平均图像已输出');