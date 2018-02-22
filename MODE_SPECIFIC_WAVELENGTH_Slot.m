
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eps structure

% %                 IPG cladding 
% % 
% %         ^           -----rw------         ^
% %         |           |           |         |  rh 
% %         h2          |           |         |
% %         |-side_half--             ---------
% %         |             core film
% %         -----------------------------------
% %         ^
% %         |
% %         h1
% %         |             SiO2 sub
% %         |
% %         -----------------------------------




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load As2S3_aneal_130

load Ge_11

load SiO2

load IPG

c0            = 299792458;                   % Speed of light in vacuum

c0            = 299792458;                   % Speed of light in vacuum

wavelength    = [1250:0.5:1645]*1e-9;               % um

wavelength_2  = [1300:20:1650]*1e-9;

wave_min      = min(wavelength);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wavelength_TeO2 = wavelength_2*1e6;                   %um

At    = 2.3771962;
Bt    = 1.7305093;
Ct    = 5.21E-02;
Dt    = 2.2568767;
Et    = 225;


n_TeO2=sqrt(At+Bt./(1-Ct./wavelength_TeO2.^2)+Dt./(1-Et./wavelength_TeO2.^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wavelength_GeAsS75S25 = wavelength_2*1e6;                   %um

B1   = 0.52019;
C1   = 0.1844;
B2   = 4.07291;
C2   = 0.04561;

n_GeA2S75Se25 = sqrt(1+B1*wavelength_GeAsS75S25.^2./(wavelength_GeAsS75S25.^2-C1)+B2*wavelength_GeAsS75S25.^2./(wavelength_GeAsS75S25.^2-C2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[s,pos_lower] = min(abs(n_SiO2(:,1)*1e-9-wave_min));

[s,pos_upper] = min(abs(n_IPG(:,1)*1e-9-wave_min));

[s,pos_core]  = min(abs(n_Ge11(:,1)*1e-9-wave_min));

% [s,pos_as2s3]  = min(abs(n_As2S3(:,1)*1e-9-wave_min));



% [s,pos_Ge11]  = min(abs(n_Ge11(:,1)*1e-9-wave_min));

omega_lower = fliplr(2*pi*c0./(n_SiO2( pos_lower : end, 1 )'*1e-9));

omega_upper = fliplr(2*pi*c0./(n_IPG(pos_upper:end,1)'*1e-9));

omega_core = fliplr(2*pi*c0./(n_Ge11(pos_core:end,1)'*1e-9));

% omega_as2s3 = fliplr(2*pi*c0./(n_As2S3(pos_as2s3:end,1)'*1e-9));

% omega_Ge11 = fliplr(2*pi*c0./(n_Ge11(pos_core:end,1)*1e-9));

omega_dis  = fliplr(2*pi*c0./wavelength_2);


P_lower       = polyfit(omega_lower,fliplr(n_SiO2( pos_lower:end,2)'),7);


P_upper       = polyfit(omega_upper,fliplr(n_IPG(pos_upper:end,2)'),7);


P_core        = polyfit(omega_core,fliplr(n_Ge11(pos_core:end,2)'),7);

% P_as2s3        = polyfit(omega_as2s3,fliplr(n_As2S3(pos_as2s3:end,2)'),7);
% P_Ge11        = polyfit(omega_core,fliplr(n_Ge11(pos_Ge11:end,2)),9);

n_lower_omega      = polyval(P_lower,omega_dis);

n_upper_omega      = polyval(P_upper,omega_dis);

n_core_omega       = polyval(P_core,omega_dis);
 
% n_as2s3_omega       = polyval(P_as2s3,omega_dis);


n_lower      = fliplr(n_lower_omega);

n_upper      = fliplr(n_upper_omega);

n_core       = fliplr(n_core_omega);

% n_as2s3       = fliplr(n_as2s3_omega);




% figure(1),plot(n_SiO2( pos_lower : end, 1 ),n_SiO2( pos_lower:end,2),'*r',wavelength_2*1e9,n_lower);
% 
% 
% figure(2),plot(n_IPG( pos_upper : end, 1 ),n_IPG( pos_upper:end,2),'*r',wavelength_2*1e9,n_upper);
% 
% 
% figure(3), plot(n_Ge11(pos_core:end,1),n_Ge11(pos_core:end,2),'*g',wavelength_2*1e9,n_core,'--r');
% 
% % figure(4), plot(n_As2S3(pos_as2s3:end,1),n_As2S3(pos_as2s3:end,2),'*g',wavelength_2*1e9,n_as2s3,'--r');
% 
% figure(5), plot(wavelength_2*1e9,n_TeO2,'--r');
% 
% figure(6), plot(wavelength_2*1e9,n_GeA2S75Se25,'--r');


n_w         = length(wavelength_2);

rw = 1;
 


h1 = 0.7;            % Lower cladding
     
h2 = 0.5;

h3 = 0;

h4 = 0.2;

h5 = 0;

h6 = 1;            % Upper cladding

h7 = 0.0;

h8 = 0.7;            

% rh = 0.3;            % Ridge height////////////////////////
% rh = h2+h3+h4;            % Ridge height

rh = h7+h6+h5+h4+h3+h2;            % Ridge height
% rw = 4.0;            % Ridge width

rw_half    = rw/2;      % Ridge half-width

rw_s_half  = h7;


side_half = 1;       % Space on one side



[min_wave, min_pos]= min(abs(wavelength_2*1e6 - 1.55));

ii = min_pos;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [in_min, in_pos]=min(abs(wavelength_2*1e9-1550));

n1     =  1.45;                 % Lower clad*

n2     =  3.67;                 % Core*

n3     =  1.62;                 % Insert Protect layer upper

n4     =  1.57;                 % Insert layer*

n5     =  1.62;                 % Insert Protect layer lower

n6     =  3.67;                 % Core*

n7     = 1.62;                  % Upper Protect layer Al2O3 

n8     = 1;                     % Upper cladding air*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% n1     =  n_lower(ii);            % Lower cladding
% n2     =  n_core(ii);             % Core
% n3     =  n_lower(ii);
% n4     =  n_upper(ii);            % Upper cladding (air)

lambda =  wavelength_2(ii)*1e6;         % vacuum wavelength

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The setting of grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx = 0.01;            % grid size (horizontal)
dy = 0.01;           % grid size (vertical)



nmodes = 6;            % number of modes to compute

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x, y, xc, yc, nx, ny, eps] = wgmesh_Slot([n1,n2,n3,n4,n5,n6,n7,n8],[h1,h2,h3,h4,h5,h6,h7,h8],rh,rw_half,side_half,rw_s_half,n7,dx,dy); %generate the eps pattern

% figure(11)
imagesc(xc,yc,eps')
% [x, y, xc, yc, nx, ny, eps] = wgmesh([n1,n2,n3],[h1,h2,h3],rh,rw_half,side_half,dx,dy); %generate the eps pattern

[x, y, xc, yc, dx, dy]      = stretch(x,y,[15,15,15,15],[6,6,6,6],'GGGG');             % stretch  the eps pattern to get bigger calculation windows 
% 
% figure(11)
% imagesc(xc,yc,eps')
% 
% save test dx dy eps x y xc yc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Hx, Hy, neff] = solver_mode(lambda, n2, nmodes, dx, dy, eps);                          % Mode calcualte

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% image_all_fields_slot(nmodes, lambda, neff, Hx, Hy, dx, dy, eps,[h1,h2,h3,h4,h5,h6,h7,h8],[n1,n2,n3,n4,n5,n6,n7,n8],rh,rw_half,side_half,xc,yc,x,y);

image_all_fields(nmodes, lambda, neff, Hx, Hy, dx, dy, eps,[h1,h2,h3,h4,h5,h6,h7,h8],[n1,n2,n3,n4,n5,n6,n7,n8],rh,rw_half,side_half,xc,yc,x,y);




















