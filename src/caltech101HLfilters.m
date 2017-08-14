HMAXparams = HMAXparameters();
gaborFilters = getGaborFilters(HMAXparams, false);

% for i=1:92
%    i
%    tic
%    load(sprintf('../res/101_ObjectCategories/caltech101gray-%d.mat',i));
%    HLfilters = getHLfilters(images,1,false);
%    save(sprintf('../res/101_ObjectCategories/caltechHLfiltersgray-%d.mat',i),'HLfilters');
%    toc
% end

HLfilters=cell(1,3);
for i=1:91
   tic
   load(sprintf('../res/101_ObjectCategories/caltech101gray-%d.mat',i));
   HLfilterstmp = getHLfilters(images,1,HMAXparams,gaborFilters,false);
   
   for j=1:size(HLfilterstmp,2)
      j
      %(i-1)*100+1
      %i*100
      %size(HLfilterstmp)
      %size(HLfilterstmp{1})
      HLfilters{j}(:,:,:,(i-1)*100+1:i*100) = HLfilterstmp{j}(:,:,:,1:100); %TODO size images instead of 100
   end
   toc
end
load('../res/101_ObjectCategories/caltech101gray-92.mat');
HLfilterstmp = getHLfilters(images,1,HMAXparams,gaborFilters,false);
for j=1:size(HLfilters,2)
  HLfilters{j}(:,:,:,9101:9101+size(HLfilterstmp{j},4)-1) = HLfilterstmp{j}(:,:,:,:);
end
save('../res/101_ObjectCategories/caltechHLfiltersgray.mat','HLfilters');