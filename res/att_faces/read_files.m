clear all;
close all;
clc;

images = uint8(zeros(112,92,400));
k=1;
for i=1:40
	for j=1:10
		A = imread(strcat('s', sprintf('%d',i), '/', sprintf('%d',j), '.pgm'));
		%figure(1);
		%imshow(A);
		%A = reshape(A,112,92,1);
		%figure(2);
		%imshow(A);
		images(:,:,(i-1)*10+j) = A(:,:,1);
		%figure(3);
		%imshow(images(:,:,(i-1)*10+j));
		%strcat('s', sprintf('%d',i), '/', sprintf('%d',j), '.pgm')
		%(i-1)*10+j
		%k = k + 1;
	end
end

%figure(4);
%imshow(images(:,:,400));
