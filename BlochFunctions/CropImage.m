function cropimg = CropImage(image,X)

    % Set pixel values X rows from the top and bottom
    image(1:X, :) = 0;
    image(end - X + 1:end, :) = 0;

    % Set pixel values X columns from the left and right
    image(:, 1:X) = 0;
    image(:, end - X + 1:end) = 0;

    cropimg = image;
return
