% The following script will lead you through some basic operations related
% to images, like opening an image and finding its dimensions and data type.
% A two-dimensional image is an n by m matrix of values, consisting of 'n' rows and
% 'm' columns. The values can be of several different data types.

clear all;
close all;
addpath('C:\path\amiedemmans\documents\MATLAB');
pic = imread('cat.tif');
% Convert to double so you can have values extending past 255.
pic = double(pic);
imwrite(pic, 'cat.tif');
% 'imread' will read the image file saved as 'cat.tif' into the variable 'pic'.
% Note the use of quotation marks with the file name. You can use any variable name
% of your choice.
% Note the use of semicolon at the end of the command to suppress the output.
% Without it, MATLAB will display the entire matrix in the command window.
%% Double comment then space creates a section break
% Create a MATLAB figure and display it.
% This line creates a figure graphics object. The parts in the brackets give the

%figure a nicely formatted title.
figure('Name','Cat','NumberTitle','off');
imagesc(pic);
colormap('gray');
axis off;

% 'figure' generates a new figure object each time you run the command. Otherwise
% MATLAB will replace the first open figure with the second and so on.
% 'imagesc' intensity scales the data and displays the image.
% 'colormap('gray')' maps the image values to a range of gray level intensities.

[xsize, ysize] = size(pic);

% xsize and ysize are now variables holding the n (or x) dimension and the m (or y) dimension respectively. You can use these in your script, for example:

filesize = xsize * ysize;

% Note that you need to convert a number to a string to display it.
disp(['The image has: ', num2str(filesize),' pixels'])

%% Exercise 1 
transposed_image = pic';
imshow(transposed_image);
figure('Name','TransposedCat','NumberTitle','off');
imagesc(pic);
colormap('gray');
axis off;
