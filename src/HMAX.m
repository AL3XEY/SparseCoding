close all
clear all
clc

tic

% Construit les filtres de Gabor utilisé par HMAX L1. On choisit l'échelle du filtre scal = 1 : 8
% et le programme fabrique les nth = 12 ici filtres correspondant aux différentes orientations.
global filter_sz = [7,11,15,19,23,27,31,35];
global sigma = [2.8,4.5,6.7,8.2,10.2,12.3,14.6,17.0];
global lambda = [3.5,5.6,7.9,10.3,12.7,15.5,18.2,21.2];
pool_sz = [8,10,12,14,16,18,22,24];
global gam = 0.3; global nth = 12;

display = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = gabor(x,y,thet,scale) %
global sigma, global lambda; global gam; global nth;
lamb = lambda(scale); % on choisit l'échelle et le filtre utilise les paramètres définis ci dessus
sig = sigma(scale);
[nx,ny] = size(x);
x0 = reshape(permute(x,[3,1,2])' * cos(thet) + permute(y,[3,1,2])' * sin(thet),nx,ny,nth);
y0 = reshape(permute(y,[3,1,2])' * cos(thet) - permute(x,[3,1,2])' * sin(thet),nx,ny,nth);
out = exp( -0.5 * (x0.^2 + gam * y0.^2) / sig^2) .* cos(2 * pi * x0 / lamb);
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imcol = imread('../res/Californie_m.JPG'); % lecture de l'image
ima = rot90(rgb2hsv(imcol)(:,:,3),0); % conversion en niveaux de gris
[dx,dy] = size(ima);
figure
imshow(uint8(255*ima)) % affichage image originale

%% Construction de la série de filtre pour les différentes orientations %%
%scal = 2; % choix de l'échelle
th = [0:nth-1]*pi/nth; % Les orientations
for scal = 1:8
	nxy = filter_sz(scal); % taille du filtre
	xx = [-(nxy-1)/2:(nxy-1)/2];
	yy = xx;
	[x,y] = meshgrid(xx,yy);

	filt = gabor(x,y,th,scal);
	filt -= mean(mean(filt)); % centrage
	filt ./= sqrt(sum(sumsq(filt))); % normalisation (norme L2)

	%affichage des filtres
	if display
		figure
		for i = [1:nth]
		  imas = filt(:,:,i);
		  nor = max(imas(:));
		  subplot(3,4,i)
		  imshow((imas/nor+1)/2)
		endfor
		%}
	end

	for i = [1:nth]
	  filtre = filt(:,:,i);
	  imf{scal}(:,:,i) = abs(filter2(filtre,ima)); % images filtrées
	endfor

	%affichage images filtrées
	if display
		figure
		for i = [1:nth]
		  imaf = imf{scal}(:,:,i);
		  vis = max(imaf(:));
		  subplot(3,4,i)
		  imshow(uint8(255*(imaf/vis + 0.3)))
		endfor
		%}
	end
end

	% taille des pools
for scal=1:7
	ps = pool_sz(scal);
	pxm = floor(dx/ps); pym = floor(dy/ps);

	for j = [1:nth] % pooling pour tous les angles
	  for px = [0:pxm-1]
	    for py = [0:pym-1]
	      impool{scal}(px+1,py+1,j) = max(max([imf{scal}(px*ps+1:(px+1)*ps,py*ps+1:(py+1)*ps,j) imf{scal+1}(px*ps+1:(px+1)*ps,py*ps+1:(py+1)*ps,j)]));
	    endfor
	  endfor
	endfor

	% affichage images filtrées et poolées
	if display
		figure
		for i = [1:nth]
		  imaf = impool{scal}(:,:,i);
		  vis = max(imaf(:));
		  subplot(3,4,i)
		  imshow(uint8(255*(imaf/vis+0.3)))
		endfor
	end

	%%%%%S2 layer - take random prototypes from C1 and compute their response to every C1 patch%%%%%

	%Taking random patches from C1
	%nHL is the number of prototypes to take
	nHL = 2;
	%nHL = 100;
	%nHL = 1000;
	sz = pool_sz(scal);
	HLfilters{scal} = zeros(sz, sz, nth, nHL);
	for cpt=1:nHL
		rdx = 1+round(rand()*(pxm-sz-1)); %TODO unique %TODO generate 1000 random numbers at once instead of one after the other to gain time
		rdy = 1+round(rand()*(pym-sz-1));
		HLfilters{scal}(1:sz,1:sz,:,cpt) = impool{scal}(rdx:rdx+sz-1,rdy:rdy+sz-1,:);

		%calculate the response of each prototype over each patch of the C1 layer (see Mutch & Lowe 2008)
		sigma2 = 1;
		alpha = (sz/4)^2;
		for x=1:pxm-sz
			for y=1:pym-sz
				X = HLfilters{scal}(:,:,:,cpt);
				P = impool{scal}(x:x+sz-1,y:y+sz-1,:);
				R = X - P;
				R = R.^2;
				R = sum(sum(sum(R)));
				R = sqrt(R);
				R = R/(2*sigma2*alpha);
				L3{scal}(x,y,cpt) = R;
			end
		end
	end
end

%Display a small part of the C1 layer
if display
	figure
	imgL3 = L3{1}(:,:,1);
	vis = max(max(max(imgL3(:,:,:))));
	imshow(uint8(255*(imgL3/vis + 0.3)))
	figure
	imgL3 = L3{1}(:,:,2);
	imshow(uint8(255*(imgL3/vis + 0.3)))
	figure
	imgL3 = L3{6}(:,:,1);
	vis = max(max(max(imgL3(:,:,:))));
	imshow(uint8(255*(imgL3/vis + 0.3)))
	figure
	imgL3 = L3{6}(:,:,2);
	imshow(uint8(255*(imgL3/vis + 0.3)))

%	figure
%	colormap gray
%	imgL3 = L3{1}(:,:,1);
%	imagesc(imgL3, [min(min(min(imgL3))) max(max(max(imgL3)))]);
%	figure
%	colormap gray
%	imgL3 = L3{1}(:,:,2);
%	imagesc(imgL3, [min(min(min(imgL3))) max(max(max(imgL3)))]);
%	figure
%	colormap gray
%	imgL3 = L3{6}(:,:,1);
%	imagesc(imgL3, [min(min(min(imgL3))) max(max(max(imgL3)))]);
%	figure
%	colormap gray
%	imgL3 = L3{6}(:,:,2);
%	imagesc(imgL3, [min(min(min(imgL3))) max(max(max(imgL3)))]);
end

%%%%% C2 layer - max response from the S2 layer %%%%%
L4 = zeros(nHL, 1);
for cpt=1:nHL
	%L3b(scal,1:nHL) = max(L3{scal}(:,:,1:nHL));
	for scal=1:7
		L3b(cpt,scal) = max(max(L3{scal}(:,:,cpt)));
	end
	L4(cpt) = max(L3b(cpt,:));
end

toc %print execution time

L4 %display final layer
