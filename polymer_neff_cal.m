clear all;
% Ridge waveguide
% %         ^        coat (n3)
% %         |
% %         | h3
% %         ^         -------------       ^
% %         |         |           |       |  rh 
% %         h2        |           |       |
% %         |   -------            ---------
% %         |         core film (n2)
% %         -   ----------------------------
% %         ^
% %         |
% %         h1
% %         |          sub (n1)
% %         |
% %         -   ------------------------------
% % 

dx0=25;                 %Grid size nm
dy0=25;                %grid size nm
side =4000;            % space on side of waveguide (nm) matters
h1 =4000;              % buffer thickness(nm)
h3 =4000;              % coating  thickess - cladding (nm)

lambda=850;      % Wavelength range
index_clad=1.5;
index_core=1.52;
n1=index_clad;%lower cladding
n3=index_clad;%upper cladding
n2=index_core;%core

h2=5000;         % film thickness 
w=5000;
rh=h2;
   
post.offset=0;
post.h=0;
post.w=0;
post.n=3.5;
post_offset=post.offset;
post_h=post.h;
post_w=post.w;
post_n=post.n;


nmodes=8;



 
 

 
    
    ii=1;
 
            dx = dy0;             % grid size (x)
            dy = dx0;             % grid size (y) 
    %         fprintf (1,'generating index mesh...\n');
 
    [x,y,xc,yc,nx,ny,eps] = ModeConventor_mesh([n1(ii),n2(ii),n3(ii)],[h1,h2,h3],...
        rh,w/2,side,dx,dy,post_offset,post_h,post_w,post_n);
    %             waveguidemeshfull([n1(ii),n2(ii),n3(ii)],[h1,h2,h3],rh,w/2,side,dx,dy);
    flipud(y);
    imagesc(x/1000,y/1000,sqrt(eps)');
    set(gca,'Ydir','Normal')
    axis equal;
     [x,y,xc,yc,dx,dy] = stretch(x,y,[20,20,40,40],[1.5,1.5,1.5,1.5]);
% % figure; imagesc(x/1000,y/1000,sqrt(eps));
% 
% 
% 
     [Hx,Hy,neff] = wgmodes (lambda(ii), n2(ii), nmodes, dx, dy, eps, '0000');
%     
 neff
 figure;plot(neff,'d');