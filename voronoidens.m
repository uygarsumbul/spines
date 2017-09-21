function area = voronoidens(kx,ky);

% function area = voronoidens(kx,ky); 
% input: kx, ky are k-space trajectories 
% output: area of cells for each point 

[row,column] = size(kx);
kxy = [kx(:),ky(:)];

% returns vertices and cells of 
% voronoi diagram 
[V,C] = voronoin(kxy);
% unpack cell array, compute 
% area of each ploygon 
area = [];
for j = 1:length(kxy) 
  x = V(C{j},1);
  y = V(C{j},2);
  lxy = length(x);
  A = abs(sum(0.5*(x([2:lxy 1])-x(:)).* ...
              (y([2:lxy 1])+y(:))));
  area = [area A];
end

area = reshape(area,row,column);
