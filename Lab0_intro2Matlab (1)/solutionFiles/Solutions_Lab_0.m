clear all;
close all;

pic = imread('cat.tif');
% Convert to double so you can have values extending past 255.
pic = double(pic);

figure('Name','Cat','NumberTitle','off')
imagesc(pic)
colormap('gray')
axis off

[xsize, ysize]=size(pic);
filesize = xsize * ysize;
% Note that you need to convert a number to a string to display it.
disp(['The image has: ', num2str(filesize),' pixels'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End of Example, start of exercises
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pic_transpose=pic';

figure('Name','Transposed Cat','NumberTitle','off')
imagesc(pic_transpose)
colormap('gray')
axis off

figure; tiledlayout(1,2,"TileSpacing","compact",'Padding','compact'); nexttile;
imagesc(pic); colormap("gray");axis off; xticks([]); yticks([]);
title('Original', 'FontSize', 14);
clim([20, 220]); axis image;
nexttile;
imagesc( pic_transpose); colormap("gray");axis off; xticks([]); yticks([]);
title('Transpose', 'FontSize', 14); 
clim([20, 220]); axis image;


%% Number 2
% Max pixel value
C_max=max(pic(:));
C_min=min(pic(:));

%% Number 3:
% Max pixel location
[xmax ymax]=find(pic_transpose==C_max);

%% Number 4: 
X = 25;
img = Lab0CropImage(pic_transpose, X);

figure('Name','Transposed Cat Cropped X = 25','NumberTitle','off')
imagesc(img)
colormap('gray')
axis off

imgC = Lab0CropImage2(pic_transpose, X);
figure('Name','Transposed Cat Cropped X = 25','NumberTitle','off')
imagesc(imgC)
colormap('gray')
axis off

%% Number 5: 
yPlot = pic_transpose(:, ysize/2);

figure; tiledlayout(1,2,"TileSpacing","compact",'Padding','compact'); nexttile;
imagesc(pic_transpose); colormap("gray");
axis off; xticks([]); yticks([]);
xline(ysize/2)
title('Transpose', 'FontSize', 14);
clim([20, 220]); axis image;
nexttile;
plot(yPlot, 'LineWidth', 3);
xlabel('Row'); ylabel('Pixel Intensity')
title('Line Profile', 'FontSize', 14); 
clim([20, 220]); axis image;

ax = gca;    ax.FontSize = 14;

%% Number 6:

Clm = twoImageTiledDisplay(pic_transpose, 25);

temp = pic_transpose;
temp(25, 25) = 600;

Clm = twoImageTiledDisplay(temp, 25);




