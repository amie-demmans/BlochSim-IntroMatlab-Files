clear all;
close all;
addpath(genpath( '/Users/amiedemmans/Documents/GitHub/test'));
pic = imread('cat.tif');
pic = double(pic);
%imwrite(pic, 'cat.tif');
%figure a nicely formatted title.
figure('Name','Cat','NumberTitle','off');
imagesc(pic);
colormap('gray');
axis off;


[xsize, ysize] = size(pic);

filesize = xsize * ysize;

disp(['The image has: ', num2str(filesize),' pixels'])

%% Exercise 1 
transposed_image = pic';
figure('Name','Transposed Cat','NumberTitle','off')
imagesc(transposed_image)
colormap('gray')
axis off

figure;
tiledlayout(1,2,"TileSpacing","compact","Padding",'compact');
nexttile;
imagesc(pic);
colormap('gray');
axis off;
xticks([]);
yticks([]);
title('Original Cat');
clim([20, 220]); axis image;

nexttile;
imagesc(transposed_image);
colormap('gray');
axis off;
xticks([]);
yticks([]);
title('Transposed Cat');
clim([20, 220]); axis image;

%% Exercise 2
min_value = min(pic(:)); 
max_value = max(pic(:)); 
fprintf('Minimum pixel value: %d\n', min_value);
fprintf('Maximum pixel value: %d\n', max_value);

%% Exercise 3 
[max_row, max_col] = find(pic == max_value);
fprintf('Max row: %d\n', max_row);
fprintf('Max col: %d\n', max_col);

%% Exercise 4 
X = 25;
cropimg = CropImage(transposed_image,X);

figure('Name','Transposed Cat Cropped X = 25','NumberTitle','off')
imagesc(cropimg)
colormap('gray')
axis off

cropimg = CropImage2(transposed_image, X);
figure('Name','Transposed Cat Cropped X = 25','NumberTitle','off')
imagesc(cropimg)
colormap('gray')
axis off

%% Exercise 5
yplot = transposed_image(:,ysize/2);

figure;
tiledlayout(1,2,"TileSpacing","compact","Padding",'compact');
nexttile;
imagesc(transposed_image);
colormap('gray');
axis off;
xticks([]);
yticks([]);
xline(ysize/2)
title('Transpose');
clim([20, 220]); axis image;

nexttile;
plot(yplot, 'LineWidth', 3);
xlabel('Row');
ylabel('Pixel Intensity');
title('Line Profile'); 
clim([20,220]);
axis image;

%% Exercise 6
Clm = twoImageTiledDisplay(transposed_image,25);

temp = transposed_image;
temp (25, 25) = 600;

Clm = twoImageTiledDisplay(temp,25);


