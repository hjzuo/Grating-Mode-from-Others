function [x, y, xc, yc, nx, ny, eps] = wgmesh_Slot(n,h,rh,rw,side,rw_s_half,n_s,dx,dy);

%

ih = round(h/dy);
irh = round (rh/dy);
irh_s = round(rw_s_half/dy);
irw = round (2*rw/dx);
iside = round (side/dx);
irw_s_h = round(rw_s_half/dx);
irw_s = irw + 2*irw_s_h;
nlayers = length(h);

nx = irw_s+2*iside+1;
ny = sum(ih)+1;

 x = dx*(-(irw_s/2+iside):1:(irw_s/2+iside))';
 xc = (x(1:nx-1) + x(2:nx))/2;

 y = (0:(ny-1))*dy;
 yc = (1:(ny-1))*dy - dy/2;
 
 for iii = 1:length(h)-3
     
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

% figure(1)
% 
% imagesc(eps');

iy = sum(ih)-ih(nlayers);
for i = 1:irh,
  eps(1:iside+irw_s_h,iy) = n_s^2*ones(iside+irw_s_h,1);
  eps(irw+iside+irw_s_h+1:irw+2*irw_s_h+2*iside,iy) = n_s^2*ones(iside+irw_s_h,1);
  iy = iy-1;
end

% figure(1)
% 
% imagesc(eps');

iy = sum(ih)-ih(nlayers)+irh_s;
for i = 1:irh,
  eps(1:iside,iy) = n(nlayers)^2*ones(iside,1);
  eps(irw_s+iside+1:irw_s+2*iside,iy) = n(nlayers)^2*ones(iside,1);
  iy = iy-1;
end

% figure(1)
% 
% imagesc(eps');

 nx = length(xc);
 ny = length(yc);
