clear all;
close all;
clc;

faces = uint8(zeros(220,220,21));
k=1;
for i=1:21
	A = imread(strcat('1', char(i+96), '/1', char(i+96), '025.pgm'));
	faces(:,:,i) = A(:,:,1);
end
