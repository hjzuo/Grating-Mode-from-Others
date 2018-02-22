function  [Hxc, Hyc, Hzc,Ezc,Exc,Eyc] = mode_all(lambda, neff, Hx, Hy, dx, dy, eps);

% This function computes the two transverse magnetic field
% components of a dielectric waveguide, using the finite
% difference method.  



[nx,ny] = size(eps);
nx = nx + 1;
ny = ny + 1;

% now we pad eps on all sides by one grid point
eps = [eps(:,1),eps,eps(:,ny-1)];
eps = [eps(1,1:ny+1);eps;eps(nx-1,1:ny+1)];



k = 2*pi/lambda;                    % free-space wavevector

if isscalar(dx)
  dx = dx*ones(nx+1,1);             % uniform grid
else
  dx = dx(:);                       % convert to column vector
  dx = [dx(1);dx;dx(length(dx))];   % pad dx on top and bottom
end

if isscalar(dy)
  dy = dy*ones(1,ny+1);             % uniform grid
else
  dy = dy(:);                       % convert to column vector
  dy = [dy(1);dy;dy(length(dy))]';  % pad dy on top and bottom
end

% distance to neighboring points to north south east and west,
% relative to point under consideration (P), as shown below.

n = ones(1,nx*ny);      n(:) = ones(nx,1)*dy(2:ny+1);
s = ones(1,nx*ny);      s(:) = ones(nx,1)*dy(1:ny);
e = ones(1,nx*ny);      e(:) = dx(2:nx+1)*ones(1,ny);
w = ones(1,nx*ny);      w(:) = dx(1:nx)*ones(1,ny);


%
%                   ------N------
%                 |       |       |
%                 |   1   n   4   |
%                 |       |       |
%                 W---w---P---e---E
%                 |       |       |
%                 |   2   s   3   |
%                 |       |       |
%                   ------S------

eps1 = ones(1,nx*ny);   eps1(:) = eps(1:nx,2:ny+1);
eps2 = ones(1,nx*ny);   eps2(:) = eps(1:nx,1:ny);
eps3 = ones(1,nx*ny);   eps3(:) = eps(2:nx+1,1:ny);
eps4 = ones(1,nx*ny);   eps4(:) = eps(2:nx+1,2:ny+1);



beta = 2*pi/lambda * neff;




bzxn= ((2.*(-1./2.*n.*eps3-1./2.*s.*eps4)./n.*eps2+2.*(1./2.*n.*eps2+1./2.*s.*eps1)./n.*eps3)./(n.*eps3+s.*eps4)./(n.*eps2+s.*eps1)./(e+w).*w.*e+((eps3.*w+eps2.*e).*(eps4-eps1).*s./n./(n+s)-(w.*eps4+e.*eps1).*(eps3-eps2).*s./n./(n+s))./(eps3.*w+eps2.*e)./(w.*eps4+e.*eps1)./(n+s).*n.*s)./beta;

bzxs =((2.*(-1./2.*n.*eps3-1./2.*s.*eps4)./s.*eps1+2.*(1./2.*n.*eps2+1./2.*s.*eps1)./s.*eps4)./(n.*eps3+s.*eps4)./(n.*eps2+s.*eps1)./(e+w).*w.*e+(-(eps3.*w+eps2.*e).*(eps4-eps1).*n./s./(n+s)+(w.*eps4+e.*eps1).*(eps3-eps2).*n./s./(n+s))./(eps3.*w+eps2.*e)./(w.*eps4+e.*eps1)./(n+s).*n.*s)./beta;

bzxe = (w./e./(e+w))./beta;
 
bzxw = (-1./w.*e./(e+w))./beta;

bzxp = (((-n.*eps3-s.*eps4).*(1./2.*n.*eps2.*(-2./w.^2-2./n.^2+k.^2.*eps1)+1./2.*s.*eps1.*(-2./w.^2-2./s.^2+k.^2.*eps2))+(n.*eps2+s.*eps1).*(1./2.*n.*eps3.*(-2./e.^2-2./n.^2+k.^2.*eps4)+1./2.*s.*eps4.*(-2./e.^2-2./s.^2+k.^2.*eps3)))./(n.*eps3+s.*eps4)./(n.*eps2+s.*eps1)./(e+w).*w.*e+((eps3.*w+eps2.*e).*(eps4-eps1).*(n-s)./n./s-(w.*eps4+e.*eps1).*(eps3-eps2).*(n-s)./n./s)./(eps3.*w+eps2.*e)./(w.*eps4+e.*eps1)./(n+s).*n.*s)./beta;


bzyn = (1./n.*s./(n+s))./beta;
  
bzys = (-n./s./(n+s))./beta;
  
bzye = (((-n.*eps3-s.*eps4).*(eps1-eps2).*w./e./(e+w)+(n.*eps2+s.*eps1).*(eps4-eps3).*w./e./(e+w))./(n.*eps3+s.*eps4)./(n.*eps2+s.*eps1)./(e+w).*w.*e+(2.*(1./2.*eps3.*w+1./2.*eps2.*e).*eps1./e-2.*(1./2.*w.*eps4+1./2.*e.*eps1).*eps2./e)./(eps3.*w+eps2.*e)./(w.*eps4+e.*eps1)./(n+s).*n.*s)./beta;
  
bzyw = ((-(-n.*eps3-s.*eps4).*(eps1-eps2).*e./w./(e+w)-(n.*eps2+s.*eps1).*(eps4-eps3).*e./w./(e+w))./(n.*eps3+s.*eps4)./(n.*eps2+s.*eps1)./(e+w).*w.*e+(2.*(1./2.*eps3.*w+1./2.*eps2.*e).*eps4./w-2.*(1./2.*w.*eps4+1./2.*e.*eps1).*eps3./w)./(eps3.*w+eps2.*e)./(w.*eps4+e.*eps1)./(n+s).*n.*s)./beta;
  
bzyp = (((-n.*eps3-s.*eps4).*(eps1-eps2).*(e-w)./e./w+(n.*eps2+s.*eps1).*(eps4-eps3).*(e-w)./e./w)./(n.*eps3+s.*eps4)./(n.*eps2+s.*eps1)./(e+w).*w.*e+((eps3.*w+eps2.*e).*(1./2.*eps4.*(-2./n.^2-2./w.^2+k.^2.*eps1).*w+1./2.*eps1.*(-2./n.^2-2./e.^2+k.^2.*eps4).*e)-(w.*eps4+e.*eps1).*(1./2.*eps3.*(-2./s.^2-2./w.^2+k.^2.*eps2).*w+1./2.*eps2.*(-2./s.^2-2./e.^2+k.^2.*eps3).*e))./(eps3.*w+eps2.*e)./(w.*eps4+e.*eps1)./(n+s).*n.*s)./beta;
 
 

ii = zeros(nx,ny);
ii(:) = (1:nx*ny); 



% Assemble sparse matrix

iall = zeros(1,nx*ny);          iall(:) = ii;
is = zeros(1,nx*(ny-1));        is(:) = ii(1:nx,1:(ny-1));
in = zeros(1,nx*(ny-1));        in(:) = ii(1:nx,2:ny);
ie = zeros(1,(nx-1)*ny);        ie(:) = ii(2:nx,1:ny);
iw = zeros(1,(nx-1)*ny);        iw(:) = ii(1:(nx-1),1:ny);


Bzx = sparse ([iall,iw,ie,is,in], ...
	[iall,ie,iw,in,is], ...
	[bzxp(iall),bzxe(iw),bzxw(ie),bzxn(is),bzxs(in)]);

Bzy = sparse ([iall,iw,ie,is,in], ...
	[iall,ie,iw,in,is], ...
	[bzyp(iall),bzye(iw),bzyw(ie),bzyn(is),bzys(in)]);


B = [Bzx , Bzy];


beta = 2*pi/lambda * neff;

Hz = zeros(size(Hx));
Hz(:) = B*reshape([Hx,Hy],2*nx*ny,1)/j/beta;







nx = nx-1;
ny = ny-1;

eps_z = eps(2:nx+1, 2:ny+1);


edet = eps_z.^2;

h = dx(2:nx+1)*ones(1,ny);
v = ones(nx,1)*dy(2:ny+1);

i1 = ii(1:nx,2:ny+1);
i2 = ii(1:nx,1:ny);
i3 = ii(2:nx+1,1:ny);
i4 = ii(2:nx+1,2:ny+1);

Dx = +neff*(Hy(i1) + Hy(i2) + Hy(i3) + Hy(i4))/4 + ...
     (Hz(i1) + Hz(i4) - Hz(i2) - Hz(i3))./(j*2*k*v);
Dy = -neff*(Hx(i1) + Hx(i2) + Hx(i3) + Hx(i4))/4 - ...
     (Hz(i3) + Hz(i4) - Hz(i1) - Hz(i2))./(j*2*k*h);
Dz = ((Hy(i3) + Hy(i4) - Hy(i1) - Hy(i2))./(2*h) - ...
      (Hx(i1) + Hx(i4) - Hx(i2) - Hx(i3))./(2*v))/(j*k);

Exc = (eps_z.*Dx)./edet;
Eyc = (eps_z.*Dy)./edet;
Ezc = Dz./eps_z;

Hyc = (Hy(i1) + Hy(i2) + Hy(i3) + Hy(i4))/4;

Hxc = (Hx(i1) + Hx(i2) + Hx(i3) + Hx(i4))/4;

Hzc = (Hz(i1) + Hz(i2) + Hz(i3) + Hz(i4))/4;


