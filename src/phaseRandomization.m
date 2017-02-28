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
[h w] = size(img);
foo = rand(h).*(max(max(PHASE))-min(min(PHASE)))+min(min(PHASE));
foo = foo(:,1:h/2);
PHASE(:,1:h/2) = foo;
PHASE(:,(h/2)+1:h) = rot90(foo,2);
imshow(mat2gray(abs(ifft2(MODULE.*(cos(PHASE) + i.*sin(PHASE)))))); %z = |z|(cos + i sin )
%imshow(mat2gray(abs(ifft2(IMG))));
