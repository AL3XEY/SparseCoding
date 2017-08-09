function [HMAXparams] = HMAXparameters()    
	HMAXparams.gam = 0.3;
	%HMAXparams.nth = 12;
	HMAXparams.nth = 4;
    
    %DO descriptors transformation matrix
    %Here, the matrix for human beings
    HMAXparams.inChans = 3;
    HMAXparams.outChansHalf = 4;
    W=ones(HMAXparams.inChans,HMAXparams.outChansHalf);
	W(2,1)=-1;
	W(3,1)=0;
	W(1,2)=2;
	W(2,2)=-1;
	W(3,2)=-1;
	W(3,3)=-2;
	W(:,1)=W(:,1)/sqrt(2);
	W(:,2)=W(:,2)/sqrt(6);
	W(:,3)=W(:,3)/sqrt(6);
	W(:,4)=W(:,4)/sqrt(3);
    HMAXparams.W = W;

	%HMAXparams.nscales = 8; % = size(...)
	%HMAXparams.sigma = [2.8,4.5,6.7,8.2,10.2,12.3,14.6,17.0];
	%HMAXparams.lambda = [3.5,5.6,7.9,10.3,12.7,15.5,18.2,21.2];
	%HMAXparams.filter_sz = [7,11,15,19,23,27,31,35];
	%HMAXparams.pool_sz = [8,10,12,14,16,18,22,24];

	HMAXparams.nscales = 4; % = size(...)
	HMAXparams.sigma = [2.8,4.5,6.7,8.2];
	HMAXparams.lambda = [3.5,5.6,7.9,10.3];
	HMAXparams.filter_sz = [7,11,15,19];
	HMAXparams.pool_sz = [8,10,12,14];
end
