function cropimg = CropImage2(image,X)

    % Remove X rows from the top and bottom
    image(1:X, :) = [];
    image(end - X + 1:end, :) = [];

    % Remove X columns from the left and right
    image(:, 1:X) = [];
    image(:, end - X + 1:end) = [];

    cropimg = image;
return
