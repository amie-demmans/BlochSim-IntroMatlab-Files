function imgC = Lab0CropImage(img, X)
% Can do this by setting values to 0, or by created a zeros matrix and
% copying over the values. I will do the form

imgC = img;
imgC(1:X,:) = 0;
imgC(end-X+1: end,:) = 0;
imgC(:, 1:X) = 0;
imgC(:, end-X+1: end) = 0;

return;