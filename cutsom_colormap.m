function[map] = cutsom_colormap(N)
% N - number of colours in the map
vec = [100:-10:0];
hex = ['#ffe700';'#ffce00';'#ffb600';'#fe9f00';'#e29300';'#c78600';'#ae7900';'#806e00';'#576002';'#354f0d';'#173d13'];
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
map = interp1(vec,raw,linspace(100,0,N),'pchip');
% %% Check maps
% figure;
% rgbplot(map);
% colormap(map);
% colorbar('Ticks',[]);
% figure;
% colormap(map);
% for iCol = 1:size(map,1)
% scatter3(map(iCol,1), map(iCol,2), map(iCol,3), 50, map(iCol,:), 'filled'); hold on;
% end
% xlabel('Red'); ylabel('Green'); zlabel('Blue');
% colorbar('Ticks',[]);
end

