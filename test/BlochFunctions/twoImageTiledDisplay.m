function Clm = twoImageTiledDisplay(transposed_image,column)

min_value = min(transposed_image(:)); 
max_value = max(transposed_image(:)); 
Clm = [min_value, max_value];

yplot = transposed_image(:,column);
%cat picture
figure;
tiledlayout(1,2,"TileSpacing","compact","Padding",'compact');
nexttile;
imagesc(transposed_image);
colormap('gray');
axis off;
xticks([]);
yticks([]);
xline(column);
title('Transpose');
clim(Clm); axis image;
colorbar;

%plot 
nexttile;
plot(yplot, 'LineWidth', 3);
xlabel('Row');
ylabel('Pixel Intensity');
title('Line Profile');

axis equal;

end