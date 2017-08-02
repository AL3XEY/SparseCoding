close all;
clear all;

size = 8;
dct = dctmtx(size);
M = zeros(size^2,size^2);
cpt = 1;
for i = 1:size
    for j = 1:size
        M(:,cpt) = reshape(dct(i,:)' * dct(j,:), size^2, 1);
        %subplot(Size,Size,((i-1)*Size) + j), imshow(M);
        cpt = cpt + 1;
    end
end
dct = M;
if(is_octave)
	save('../results/dct8.mat', 'dct', '-v7');
else
	save('../results/dct8.mat', 'dct');
end
