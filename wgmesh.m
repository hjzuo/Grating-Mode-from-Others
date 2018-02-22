function [x, y, xc, yc, nx, ny, eps] = wgmesh(n,h,rh,rw,side,dx,dy);

%

ih = round(h/dy);
irh = round (rh/dy);
irw = round (2*rw/dx);
iside = round (side/dx);
nlayers = length(h);

nx = irw+2*iside+1;
ny = sum(ih)+1;

 x = dx*(-(irw/2+iside):1:(irw/2+iside))';
 xc = (x(1:nx-1) + x(2:nx))/2;

 y = (0:(ny-1))*dy;
 yc = (1:(ny-1))*dy - dy/2;
 
 for iii = 1:length(h)-2
     
     y =  y - h(iii);
     
     yc = yc - h(iii);
     
 end

eps = zeros(nx-1,ny-1);

iy = 1;

for jj = 1:nlayers,
  for i = 1:ih(jj),
	eps(:,iy) = n(jj)^2*ones(nx-1,1);
	iy = iy+1;
  end
end

iy = sum(ih)-ih(nlayers);
for i = 1:irh,
  eps(1:iside,iy) = n(nlayers)^2*ones(iside,1);
  eps(irw+iside+1:irw+2*iside,iy) = n(nlayers)^2*ones(iside,1);
  iy = iy-1;
end

 nx = length(xc);
 ny = length(yc);
