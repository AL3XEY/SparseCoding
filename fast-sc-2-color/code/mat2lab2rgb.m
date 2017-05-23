function [img] = mat2lab2rgb(img)
    img(:,:,1) = (img(:,:,1)+0.5).*100;
    img(:,:,2) = ((img(:,:,2)+0.5).*255)-128;
    img(:,:,3) = ((img(:,:,3)+0.5).*255)-128;
    img = lab2rgb(img);
end