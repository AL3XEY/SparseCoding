close all;
clear all;
clc;

pkg load image;

img = imread('../res/lena.pgm');
[h w] = size(img);
figure;
imshow(img);

pyramid(img, 8, 'Laplacian', true);

%phaseRandomization(img);

% see rotv(angle) from the linear algebra package
