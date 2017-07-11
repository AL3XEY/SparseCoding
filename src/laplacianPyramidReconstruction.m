function [imgout] = laplacianPyramidReconstruction(Lplac)

Nlev = size(Lplac,1)-1;
hgf = [1,4,6,4,1]/16; % demi filtre binomial à 5 points
[gfx,gfy] = meshgrid(hgf,hgf);
Gausfilt = gfx.*gfy; % filtre binomial

% Reconstruction de l'image d'origine :
bld{Nlev+1} = Lplac{Nlev+1}; % sommet de la pyramide Laplacienne
for i =Nlev:-1:1
	ini = bld{i+1};
	[nH,nL] = size(ini);  % sa taille
	exp0 = zeros(2*nH,2*nL); % matrice dilatée
	exp0(2:2:2*nH,2:2:2*nL) = ini; % remplissage des noeuds pairs
	expi = conv2(exp0,4*Gausfilt,'same'); % remplissage des noeuds impairs par interpolation
	bld{i} = Lplac{i} + expi; % reconstruction de l'image au rang i. Au nveau 1 on retrouve l'image initiale.
end

imgout = bld{1};

%check = max(abs(bld{1}-ima)(:));
%prec = 1E-10;
%if(check>prec)
%	disp('l image reconstruite est differente de l image initiale !!!')
%end
