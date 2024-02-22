function Clm = twoImageTiledDisplay(img, column)

C_max=max(img(:));
C_min=min(img(:));
Clm = [C_min, C_max];

yPlot = img(:, column);

figure; tiledlayout(1,2,"TileSpacing","compact",'Padding','compact'); nexttile;
imagesc(img); colormap("gray");
axis off; xticks([]); yticks([]);
xline(column)
title('Image', 'FontSize', 14);
clim([C_min, C_max]); axis image;
nexttile;
plot(yPlot, 'LineWidth', 3);
xlabel('Row'); ylabel('Pixel Intensity')
title(['Line Profile at Column: ', num2str(column)], 'FontSize', 14); 

ax = gca;    ax.FontSize = 14;













