


% c = [' .r',' ob',' xg',' +c',' *m',' sy',' dr',' vb',' <g','--c','--m','--y','-.r','-.b','-.g','-.c','-.m','-.y'];

c = ['-*r','-*m','-*y','-*b','-*c','-*g','  r','  b','  g','  c','  m','  y','--r','--b','--g','--c','--m','--y'];

c = reshape(c,3,18);

c = c';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TM1 mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c0         = 299792458;  

omega_dis  = fliplr(2*pi*c0./wavelength_2);

wave_t     = [1350:0.5:1600]*1e-9;

omega      = fliplr(2*pi*c0./wave_t);

P_TM0 = polyfit(omega_dis,fliplr(neff_TM(:,1)'),11);

n_TM0_omega = polyval(P_TM0,omega);

n_TM0 =fliplr(n_TM0_omega);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      beta0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 

% 


% 
k        = omega/c0;

% 
beta0_omega_TM0 = k.*n_TM0_omega;

beta0_TM0 = fliplr(beta0_omega_TM0);

figure(101),
plot(wave_t*1e6, beta0_TM0);


ylabel('TM0 Beta0 (1/m)')

xlabel('WAVELENGTH(\mum)');

box on
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        beta1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(P_TM0) - 1;

P_multi_TM0 = [n+1:-1:1];

P_beta1_TM0 = P_TM0.*P_multi_TM0/c0;

beta1_omega_TM0 = polyval(P_beta1_TM0,omega);

beta1_TM0 = fliplr(beta1_omega_TM0);

figure(102),
plot(wave_t*1e6, beta1_TM0);
ylabel('TM0 Beta1 (s/m)')

xlabel('WAVELENGTH(\mum)');

box on
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  beta2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P_multi_TM0 = [n:-1:1];

P_beta2_TM0 = P_beta1_TM0(1:end-1).*P_multi_TM0;

beta2_omega_TM0 = polyval(P_beta2_TM0,omega);

beta2_TM0 = fliplr(beta2_omega_TM0);

figure(103),
plot(wave_t*1e6, beta2_TM0);
ylabel('TM0 Beta2 (s^2/m)')

xlabel('WAVELENGTH(\mum)');

box on
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  beta3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P_multi_TM0 = [n-1:-1:1];

P_beta3_TM0 = P_beta2_TM0(1:end-1).*P_multi_TM0;

beta3_omega_TM0 = polyval(P_beta3_TM0,omega);

beta3_TM0 = fliplr(beta3_omega_TM0);

figure(104),
plot(wave_t*1e6, beta3_TM0);
ylabel('TM0 Beta3 (s^3/m)')

xlabel('WAVELENGTH(\mum)');

box on
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  beta4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P_multi_TM0 = [n-2:-1:1];

P_beta4_TM0 = P_beta3_TM0(1:end-1).*P_multi_TM0;

beta4_omega_TM0 = polyval(P_beta4_TM0,omega);

beta4_TM0 = fliplr(beta4_omega_TM0);

figure(105),
plot(wave_t*1e6, beta4_TM0);
ylabel('TM0 Beta4 (s^4/m)')

xlabel('WAVELENGTH(\mum)');

box on
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  beta5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P_multi_TM0 = [n-3:-1:1];

P_beta5_TM0 = P_beta4_TM0(1:end-1).*P_multi_TM0;

beta5_omega_TM0 = polyval(P_beta5_TM0,omega);

beta5_TM0 = fliplr(beta5_omega_TM0);

figure(106),
plot(wave_t*1e6, beta5_TM0);
ylabel('TM0 Beta5 (s^5/m)')

xlabel('WAVELENGTH(\mum)');

box on
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  beta6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P_multi_TM0 = [n-4:-1:1];

P_beta6_TM0 = P_beta5_TM0(1:end-1).*P_multi_TM0;

beta6_omega_TM0 = polyval(P_beta6_TM0,omega);

beta6_TM0 = fliplr(beta6_omega_TM0);

figure(107),
plot(wave_t*1e6, beta6_TM0);
ylabel('TM0 Beta6 (s^6/m)')

xlabel('WAVELENGTH(\mum)');

box on
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  beta7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P_multi_TM0 = [n-5:-1:1];

P_beta7_TM0 = P_beta6_TM0(1:end-1).*P_multi_TM0;

beta7_omega_TM0 = polyval(P_beta7_TM0,omega);

beta7_TM0 = fliplr(beta7_omega_TM0);

figure(108),
plot(wave_t*1e6, beta7_TM0);
ylabel('TM0 Beta7 (s^7/m)')

xlabel('WAVELENGTH(\mum)');

box on
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  beta8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P_multi_TM0 = [n-6:-1:1];

P_beta8_TM0 = P_beta7_TM0(1:end-1).*P_multi_TM0;

beta8_omega_TM0 = polyval(P_beta8_TM0,omega);

beta8_TM0 = fliplr(beta8_omega_TM0);

figure(109),
plot(wave_t*1e6, beta8_TM0);
ylabel('TM0 Beta8 (s^8/m)')

xlabel('WAVELENGTH(\mum)');

box on
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D_TM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D_TM0 = - 2*pi*c0./wave_t.^2.*beta2_TM0*1e6;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[min_b,pos_b] = min(abs(wave_t*1e6 - 1.550));

beta = [beta0_TM0',beta1_TM0',beta2_TM0',beta3_TM0',beta4_TM0',beta5_TM0',beta6_TM0',beta7_TM0',beta8_TM0'];


for ii = 1 : 9
    
    
    beta_tt = beta(pos_b,ii);
    fprintf(sprintf('beta%d = %e\n',ii-1,beta_tt));

end


D_TM0_tt = D_TM0(pos_b);
fprintf(sprintf('D = %e\n',D_TM0_tt));

figure(1001),
plot(wave_t(1:10:end)*1e6, D_TM0(1:10:end),'-*r');

ylabel('Dispersion parameter TM0 (ps/km/m)')

xlabel('WAVELENGTH(\mum)');
axis([1.350,1.6,-400,200])

xx = xlim;

yy = ylim;

line([xx(1),xx(2)],[0,0],'color','b');

line([1.55,1.55],[yy(1),yy(2)],'color','g')

box on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  nonlinearity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% n2 = 8.65e-18;                       % m^2/w

n2  = 8.65e-18;

Aeff = interp1(wavelength_2,Aeff_TM(:,1)*1e-12,wave_t);  % um^2 - m^2

gamma = n2*fliplr(omega)/c0./Aeff;

[min_g,pos_g] = min(abs(wave_t*1e6 - 1.550));

gamma_tt = gamma(pos_g);

Aeff_tt  = Aeff(pos_g);
fprintf(sprintf('gamma = %f\n',gamma_tt));
fprintf(sprintf('Aeff = %e\n',Aeff_tt));


figure(301),plot(wave_t*1e6,gamma,'-*c');

yy = ylim;


axis([min(wave_t*1e6),max(wave_t*1e6),yy(1),yy(2)]);

xlabel('WAVELENGTH(\mum)');

ylabel('GAMMA (1/W/M)');


box on
