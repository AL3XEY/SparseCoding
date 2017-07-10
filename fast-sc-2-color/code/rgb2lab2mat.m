function [img] = rgb2lab2mat(img)
    img = rgb2lab(img);
    img(:,:,1) = (img(:,:,1)./100) - 0.5;
    img(:,:,2) = (img(:,:,2)+128)./255 - 0.5;
    img(:,:,3) = (img(:,:,3)+128)./255 - 0.5;
end