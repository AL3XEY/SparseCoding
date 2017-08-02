function h=display_network_nonsquare(A, channels, numcols, figstart)
%  display_network -- displays the state of the network
%  A = basis function matrix

warning off all

if exist('figstart', 'var') && ~isempty(figstart), figure(figstart); end

[L,M]=size(A);
if ~exist('numcols', 'var')
    numcols = ceil(sqrt(L/channels));
    while mod(L/channels, numcols), numcols= numcols+1; end
end
ysz = numcols;
xsz = ceil(L/channels/ysz);

m=floor(sqrt(M*ysz/xsz));
n=ceil(M/m);

if channels==1
	colormap gray
end

buf=1;
array=zeros(buf+m*(xsz+buf),buf+n*(ysz+buf),3);

k=1;
for i=1:m
    for j=1:n
        if k>M continue; end
		for chan=1:channels
		    clim=max(abs(A((chan-1)*L/channels+1:chan*L/channels,k)));
		    array(buf+(i-1)*(xsz+buf)+[1:xsz],buf+(j-1)*(ysz+buf)+[1:ysz],chan)=...
		        reshape(A((chan-1)*L/channels+1:chan*L/channels,k),xsz,ysz)/clim;
		end
        k=k+1;
    end
end

if channels==3
	array = mat2lab2rgb(array./2);
	mn = min(min(min(array)));
	array = array-mn;
	mx = max(max(max(array)));
	array = array./mx;
	array = array*2-1;
end

%min(min(min(array)))
%max(max(max(array)))

if isreal(array)
    if is_octave
      h=imagesc(array,[-1 1],'EraseMode','none');
    else
      h=imagesc(array,'EraseMode','none',[0 1]);
      %h=imshow(array);
    end
else
    if is_octave
      h=imagesc(20*log10(abs(array)),[-1 1],'EraseMode','none');
    else
      h=imagesc(20*log10(abs(array)),'EraseMode','none',[0 1]);
      %h=imshow(20*log10(abs(array)));
    end
end;
axis image off

drawnow

warning on all
