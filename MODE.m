clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eps structure

% %                 IPG cladding 
% % 
% %         ^           -----rw------         ^
% %         |           |           |         |  rh 
% %         h2          |           |         |
 %          |-side_half--             ---------
% %         |             core film
% %         -----------------------------------
% %         ^
% %         |
% %         h1
% %         |             SiO2 sub
% %         |
% %         -----------------------------------




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n1 = 1.44;           % Lower cladding
n2 = 2.44;           % Core
n3 = 1.00;           % Upper cladding (air)

% Layer heights:
h1 = 0.5;            % Lower cladding
h2 = 0.9;            % Core thickness
h3 = 0.5;            % Upper cladding

rh = 0.4;            % Ridge height

rw = 4.0;            % Ridge width

rw_half = rw/2;      % Ridge half-width
side_half = 1;       % Space on one side
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The setting of grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx = 0.025;            % grid size (horizontal)
dy = 0.0125;           % grid size (vertical)

lambda = 1.55;         % vacuum wavelength

nmodes = 6;            % number of modes to compute

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x, y, xc, yc, nx, ny, eps] = wgmesh([n1,n2,n3],[h1,h2,h3],rh,rw_half,side_half,dx,dy); %generate the eps pattern


[x, y, xc, yc, dx, dy]      = stretch(x,y,[20,30,20,20],[4,4,3,3],'GGGG');              % stretch  the eps pattern to get bigger calculation windows 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Hx, Hy, neff] = solver_mode(lambda, n2, nmodes, dx, dy, eps);                          % Mode calcualte

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

image_all_fields(nmodes, lambda, neff, Hx, Hy, dx, dy, eps,[h1,h2,h3],[n1,n2,n3],rh,rw_half,side_half,xc,yc);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saving
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mkdir('D:\MATLAB\Xin\smode5','Lambda1550_rw4_rh04');
% 
% s = ['save D:\MATLAB\Xin\smode5\Lambda1550_rw4_rh04\data'];
% 
% eval(s);
