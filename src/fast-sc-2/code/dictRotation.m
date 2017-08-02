clear all
close all

load ../results/dict/4/sc_L1_b128_beta0.4_20170228T202450.mat
winsize = sqrt(size(B,1));
for i=1:size(B,2)
	foo = reshape(B(:,i), winsize, winsize);
	bar = reshape(rot90(foo), winsize^2,1);
	B(:,i) = bar(:,1);
end
