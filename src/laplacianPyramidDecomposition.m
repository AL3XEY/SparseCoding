function [Lplac] = laplacianPyramidDecomposition(ima, Nlev)
%ima : image in range [0 1]

hgf = [1,4,6,4,1]/16; % demi filtre binomial à 5 points
[gfx,gfy] = meshgrid(hgf,hgf);
Gausfilt = gfx.*gfy; % filtre binomial

if nargin < 2 || isempty(Nlev)
	Nlev = 4 % nombre de niveaux de la pyramide.
end

% Construction de la pyramide de Gauss-Laplace :
 Gau{1} = ima;
for i =1:Nlev
	imi = Gau{i};
	[nH,nL] = size(imi); % taille de l'image
	imf = conv2(imi,Gausfilt,'same'); % convolution de l'image par le filtre
	imp = imf(2:2:nH,2:2:nL); % sous echantillonage de l'image
	Gau{i+1} = imp;

	% expansion de l'image avec interpolation
	exp0 = zeros(nH,nL);
	exp0(2:2:nH,2:2:nL) = Gau{i+1}; % on recopie l'image réduite sur les noeuds pairs
	expa = conv2(exp0,4*Gausfilt,'same'); % on convolue ce qui précède avec le filtre binomial.
	Lplac{i} = imi - expa;
end
Lplac{Nlev+1} = Gau{Nlev+1};
