function  [Hx, Hy, neff]=solver_mode(lambda, guess, nmodes, dx, dy, eps);

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



ns21 = n.*eps2+s.*eps1;
ns34 = n.*eps3+s.*eps4;
ew14 = e.*eps1+w.*eps4;
ew23 = e.*eps2+w.*eps3;

axxn = ((2*eps4.*e).*(eps3./eps4)./ns34 + ...
        (2*eps1.*w).*(eps2./eps1)./ns21)./(n.*(e+w));

axxs = ((2*eps3.*e).*(eps4./eps3)./ns34 + ...
        (2*eps2.*w).*(eps1./eps2)./ns21)./(s.*(e+w));

ayye = (2.*n.*eps4).*eps1./eps4./e./ew14./(n+s) + ...
       (2.*s.*eps3).*eps2./eps3./e./ew23./(n+s);

ayyw = (2.*eps1.*n).*eps4./eps1./w./ew14./(n+s) + ...
       (2.*eps2.*s).*eps3./eps2./w./ew23./(n+s);

axxe = 2./(e.*(e+w));
       

axxw = 2./(w.*(e+w)) ;

ayyn = 2./(n.*(n+s)) ;

ayys = 2./(s.*(n+s)) ;





axxp = - axxn - axxs - axxe - axxw  ...
       + k^2*(n+s).*(eps4.*eps3.*e./ns34 + eps1.*eps2.*w./ns21)./(e+w);

ayyp = - ayyn - ayys - ayye - ayyw  ...
       + k^2*(e+w).*(eps1.*eps4.*n./ew14 + eps2.*eps3.*s./ew23)./(n+s);


axye = -2*(eps1.*eps2.*(1./eps1-1./eps2).*w.^2./ns21 + ...
           eps3.*eps4.*(1./eps4-1./eps3).*e.*w./ns34)./e./(e+w).^2;

axyw =  -2*(eps4.*eps3.*(1./eps3-1./eps4).*e.^2./ns34 + ...
          eps2.*eps1.*(1./eps2-1./eps1).*w.*e./ns21)./w./(e+w).^2;

ayxn = - 2*(eps3.*eps2.*(1./eps3-1./eps2).*s.^2./ew23 + ...
          eps1.*eps4.*(1./eps4-1./eps1).*n.*s./ew14)./n./(n+s).^2;

ayxs = - 2*(eps4.*eps1.*(1./eps1-1./eps4).*n.^2./ew14 + ...
          eps2.*eps3.*(1./eps2-1./eps3).*s.*n./ew23)./s./(n+s).^2;


axyp = -( axye + axyw ) ;

ayxp = -( ayxn + ayxs ) ;  

ii = zeros(nx,ny);
ii(:) = (1:nx*ny); 



% Assemble sparse matrix

iall = zeros(1,nx*ny);          iall(:) = ii;
is = zeros(1,nx*(ny-1));        is(:) = ii(1:nx,1:(ny-1));
in = zeros(1,nx*(ny-1));        in(:) = ii(1:nx,2:ny);
ie = zeros(1,(nx-1)*ny);        ie(:) = ii(2:nx,1:ny);
iw = zeros(1,(nx-1)*ny);        iw(:) = ii(1:(nx-1),1:ny);


Axx = sparse ([iall,iw,ie,is,in], ...
	[iall,ie,iw,in,is], ...
	[axxp(iall),axxe(iw),axxw(ie),axxn(is),axxs(in)]);

Axy = sparse ([iall,iw,ie], ...
	[iall,ie,iw], ...
	[axyp(iall),axye(iw),axyw(ie)]);

Ayx = sparse ([iall,is,in], ...
	[iall,in,is], ...
	[ayxp(iall),ayxn(is),ayxs(in)]);

Ayy = sparse ([iall,iw,ie,is,in], ...
	[iall,ie,iw,in,is], ...
	[ayyp(iall),ayye(iw),ayyw(ie),ayyn(is),ayys(in)]);

A = [[Axx Axy];[Ayx Ayy]];

clear Axx Axy Ayx Ayy ...
    
shift = (guess*k)^2;
options.tol = 1e-8;
options.disp = 0;						% suppress output

[H,beta] = eigs(A,speye(size(A)),nmodes,shift,options);
neff = lambda*sqrt(diag(beta))/(2*pi);

Hx = zeros(nx,ny,nmodes);
Hy = zeros(nx,ny,nmodes);


% Normalize modes

temp = zeros(nx*ny,2);
for kk = 1:nmodes;
  temp(:) = H(:,kk);
  [mag,ii] = max(sqrt(sum(abs(temp).^2,2)));
  if abs(temp(ii,1)) > abs(temp(ii,2)),
    jj = 1;
  else 
    jj = 2;
  end
  mag = mag*temp(ii,jj)/abs(temp(ii,jj));
  temp = temp/mag;
  Hx(:,:,kk)   = reshape(temp(:,1),nx,ny);
  Hy(:,:,kk)   = reshape(temp(:,2),nx,ny);
end;

return;

