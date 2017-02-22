function [  ] = phaseRandomization( img )
%[WIP] Perform a phase randomization on given image
%

IMG = fft2(img);
IMG2 = IMG;%fftshift(IMG);
MODULE = abs(IMG2);
MODULE2 = log(MODULE+1);
PHASE = arg(IMG2);
figure;
imshow(mat2gray(MODULE2));
figure
imshow(mat2gray(PHASE));

figure;
PHASE = rot90(PHASE);
%PHASE = randmatrix(h,w) between min(min(PHASE)) and max(max(PHASE))
imshow(mat2gray(abs(ifft2(MODULE.*(cos(PHASE) + i.*sin(PHASE))))));
%imshow(mat2gray(abs(ifft2(IMG))));
