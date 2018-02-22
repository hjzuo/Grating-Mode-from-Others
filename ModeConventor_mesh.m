function [x, y, xc, yc, nx, ny, eps] = ModeConventor_mesh(n,h,rh,rw,side,dx,dy,post_offset,post_h,post_w,post_n)

%

ih = round(h/dy);
irh = round (rh/dy);
irw = round (2*rw/dx);
iside = round (side/dx);
nlayers = length(h);
ipost_offset=round(post_offset/dx);
ipost_h=round(post_h/dy);
ipost_w=round(post_w/dx);


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


iy = sum(ih)-ih(nlayers);
for i = 1:ipost_h,
  %eps(1:iside,iy) = post_n^2*ones(iside,1);
%   irw/2+iside+ipost_offset-ipost_w/2
%   irw/2+iside+ipost_offset+ipost_w/2
%   ipost_w
  eps(irw/2+iside+ipost_offset-ipost_w/2+1:irw/2+iside+ipost_offset+ipost_w/2,iy) = post_n^2*ones(ipost_w,1);
  iy = iy+1;
end



 nx = length(xc);
 ny = length(yc);
