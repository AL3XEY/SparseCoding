function [sigma lambda gam nth nscales filter_sz pool_sz] = HMAXparameters()
	gam = 0.3;
	nth = 12;
	%nscales = 8; % = size(...)
	%sigma = [2.8,4.5,6.7,8.2,10.2,12.3,14.6,17.0];
	%lambda = [3.5,5.6,7.9,10.3,12.7,15.5,18.2,21.2];
	%filter_sz = [7,11,15,19,23,27,31,35];
	%pool_sz = [8,10,12,14,16,18,22,24];

	nscales = 4; % = size(...)
	sigma = [2.8,4.5,6.7,8.2];
	lambda = [3.5,5.6,7.9,10.3];
	filter_sz = [7,11,15,19];
	pool_sz = [8,10,12,14];
end
