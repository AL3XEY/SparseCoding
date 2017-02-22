close all;
clear all;
clc;

pkg load image;

img = imread('../res/lena.pgm');
[h w] = size(img);
figure;
imshow(img);

pyramid(img, 8, 'Laplacian', true);
