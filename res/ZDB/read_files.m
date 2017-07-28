close all;
clear all;
clc;

%pkg load image;

scale = 0.25;

img = imread('DSCN1572.JPG');
img2 = imresize(img, scale);
[h w c] = size(img);
[h2 w2 c2] = size(img2);
zdbS = zeros(h2,w2,128);
zdb = zeros(h,w,128);

for i=1:8
  disp(i+1571)
  img = imread(sprintf('DSCN%d.JPG', i+1571));
  img2 = rgb2gray(img);
  img2 = mat2gray(img2, [0 255])-0.5;
  zdb(:,:,i) = img2(:,:,1);
  img2 = imresize(img2, scale);
  zdbS(:,:,i) = img2(:,:,1);
end

for i=9:11
  disp(i+2056)
  img = imread(sprintf('DSCN%d.JPG', i+2056));
  img2 = rgb2gray(img);
  img2 = mat2gray(img2, [0 255])-0.5;
  zdb(:,:,i) = img2(:,:,1);
  img2 = imresize(img2, scale);
  zdbS(:,:,i) = img2(:,:,1);
end

for i=12:18
  disp(i+2058)
  img = imread(sprintf('DSCN%d.JPG', i+2058));
  img2 = rgb2gray(img);
  img2 = mat2gray(img2, [0 255])-0.5;
  zdb(:,:,i) = img2(:,:,1);
  img2 = imresize(img2, scale);
  zdbS(:,:,i) = img2(:,:,1);
end

for i=19:22
  disp(i+2059)
  img = imread(sprintf('DSCN%d.JPG', i+2059));
  img2 = rgb2gray(img);
  img2 = mat2gray(img2, [0 255])-0.5;
  zdb(:,:,i) = img2(:,:,1);
  img2 = imresize(img2, scale);
  zdbS(:,:,i) = img2(:,:,1);
end

for i=23:23
  disp(i+2060)
  img = imread(sprintf('DSCN%d.JPG', i+2060));
  img2 = rgb2gray(img);
  img2 = mat2gray(img2, [0 255])-0.5;
  zdb(:,:,i) = img2(:,:,1);
  img2 = imresize(img2, scale);
  zdbS(:,:,i) = img2(:,:,1);
end

for i=24:40
  disp(i+2061)
  img = imread(sprintf('DSCN%d.JPG', i+2061));
  img2 = rgb2gray(img);
  img2 = mat2gray(img2, [0 255])-0.5;
  zdb(:,:,i) = img2(:,:,1);
  img2 = imresize(img2, scale);
  zdbS(:,:,i) = img2(:,:,1);
end

for i=41:48
  disp(i+2063)
  img = imread(sprintf('DSCN%d.JPG', i+2063));
  img2 = rgb2gray(img);
  img2 = mat2gray(img2, [0 255])-0.5;
  zdb(:,:,i) = img2(:,:,1);
  img2 = imresize(img2, scale);
  zdbS(:,:,i) = img2(:,:,1);
end

for i=49:128
  disp(i+2064)
  img = imread(sprintf('DSCN%d.JPG', i+2064));
  img2 = rgb2gray(img);
  img2 = mat2gray(img2, [0 255])-0.5;
  zdb(:,:,i) = img2(:,:,1);
  img2 = imresize(img2, scale);
  zdbS(:,:,i) = img2(:,:,1);
end

%figure;
%imshow(zdb(:,:,1));
%figure;
%colormap gray;
%imagesc(zdb(:,:,1), [-0.5 0.5]);
save('zdbScaled.mat', '-v7.3', 'zdb');
zdb = zdbS;
save('zdbSmallScaled.mat', '-v7.3', 'zdb');