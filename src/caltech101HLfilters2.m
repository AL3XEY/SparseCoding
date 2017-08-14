HMAXparams = HMAXparameters();
gaborFilters = getGaborFilters(HMAXparams, false);
for i=1:92
   tic
   load(sprintf('../res/101_ObjectCategories/caltech101gray-%d.mat',i));
   %HLfilterstmp = getHLfilters(images,1,false);
   for j=1:size(images,2)
       imagestmp{((i-1)*100)+j} = images{j};
   end
   toc
end
tic
HLfilters = getHLfilters(images,1000,HMAXparams,gaborFilters,false);
toc
save('../res/101_ObjectCategories/caltechHLfiltersgray.mat','HLfilters');
%save('../res/101_ObjectCategories/caltechHLfiltersgray1.mat','HLfilters');