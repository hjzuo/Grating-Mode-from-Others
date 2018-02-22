clear
load Material_half

% load Material_Raman

pyl = [0.07,0.1,0.15];

P_pyl = pyl/0.005/140;

S_i =[];

S_i_fwm = [];

S_Raman_p = [];

FOM_p     = [];

P_Raman_p=[];

% for ii = 1:length(P_pyl);

% P = 0.1429;

P = P_pyl(2);


% loss = 1.5; %dB/cm
% 
% loss_l = loss/4.343*100;

L = 0.005;

% L = [1-exp(-loss_l*L)]/loss_l;

i = sqrt(-1);

f = 0.13;

% h_i = fftshift(h_w_i_Ge11);
% 
% h_r = fftshift(h_w_r_Ge11);
% 
% w_f_f = fftshift(w_f_Ge11);
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The quoted data is on higher frequency side---anti-stokes side
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h_i = h_w_i_Ge11;

h_r = h_w_r_Ge11;

w_f_f = w_f_Ge11;

[min_s,pos_s]= min(abs(w_f_f-119872004869917));

% figure(101)
% 
% plot(w_f_f(pos_s:end-10000),h_i(pos_s:end-10000))

a = [w_f_f(pos_s:end-10000)',h_i(pos_s:end-10000)'];
 
h_i(pos_s:end)=  -0.005*log10(w_f_f(pos_s:end))./log10(exp(1)) + 0.1752;



% figure(102)
% 
% plot(w_f_f(pos_s-500:end-1000),h_i(pos_s-500:end-1000));

P_ss = polyfit(w_f_f(pos_s-300:5:pos_s+200),h_i(pos_s-300:5:pos_s+200),10);
% 
h_i(pos_s-300:pos_s+200)=polyval(P_ss,w_f_f(pos_s-300:pos_s+200));


% figure(103)
% plot(w_f_f(pos_s-300:pos_s+300),h_i(pos_s-300:pos_s+300));
% 
% figure(104)
% 
% plot(w_f_f,h_i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% q = 1-f+f*h_r-f*h_i*i; % Total Raman
q = 1-f+f*h_r;  % Real Raman
% q = 1-f-f*h_i*i;
q_f = 1 ; % pure FWM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Dispersion and nonlinearity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beta_2_s = [-1.903279e-027, -6.573574e-027,-1.927880e-026];
% beta_4_s = [-6.987405e-055,-6.863090e-055,-4.825687e-055];

beta_2_s = [-1.903279e-027, -6.573574e-027,-1.182290e-026,-1.927880e-026];
beta_4_s = [-6.987405e-055,-6.863090e-055,-5.254010e-055,-4.825687e-055];

% beta_2_s = [1.903279e-027, -6.573574e-027,1.927880e-026];
% beta_4_s = [6.987405e-055,-6.863090e-055,4.825687e-055];


% beta_2_s = [-1.903279e-027, -6.573574e-027];
% beta_4_s = [-6.987405e-055,-6.863090e-055];



for ii = 1:length(beta_2_s)

beta2 = beta_2_s(ii);

beta4 = beta_4_s(ii);
    
% beta2 = -1.903279e-027;  %% zero-dispersion2
% % beta2 = -1.903279e-027;  %% zero-dispersion2
% % beta3 = 1.186135e-039;
% beta4 = -6.987405e-055;
% % beta5 = -3.044915e-069;
% % beta6 = 5.818785e-083;
% % beta7 = -5.971047e-097;
% % beta8 = 6.001905e-111;
% D = 1.501917e+000;

% beta2 = -6.573574e-027;         %% zero-dispersion1
% beta4 = -6.863090e-055;
% D = 5.153933e+000;

% beta2 = -1.927880e-026;    %% 4um
% beta4 = -4.825687e-055;
% D = 1.501827e+001;

% beta2 = -3.712364e-026
% beta4 = 3.995749e-054
% D = 3.066883e+001



gamma = 140;

deta_k = beta2*w_f_f.^2+beta4/12*w_f_f.^4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K = -deta_k/2/gamma/P;


R = sqrt(K.*(2*q-K));

R_fwm = sqrt(K.*(2*q_f-K));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conversion effeciency  Pi(z)/Ps(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Conver_i_anti = abs(q./(K-q+i*R./tanh(gamma*P*R*L))).^2; %signal in stokes and idler in anti-stokes

Conver_i_stokes = abs(q./(K-q-i*R./tanh(gamma*P*R*L))).^2; %signal in anti-stokes and idler in stokes

% Conver_i_stokes = abs((i*q./R.*sinh(gamma*P*R*L))).^2;

Conver_i_stokes_fwm = abs(q_f./(K-q_f-i*R_fwm./tanh(gamma*P*R_fwm*L))).^2;

% G = abs((i*q./R.*sinh(gamma*P*R*L))).^2;                % Pi(z)/Ps(0)

% G = abs((K-q-i*R./tanh(gamma*P*R*L)./q)).^2;
% [w_min,w_pos]=min(abs(w_f_f-0));
% G   = abs((cosh(gamma*P*R*L)+i*(K-q)./R.*sinh(gamma*P*R*L))).^2;
% G_s = abs(q(1:w_pos)./(K(1:w_pos)-q(1:w_pos)+i*R(1:w_pos)./tanh(gamma*P*R(1:w_pos)*L))).^2;
% 
% G_a = abs(q(w_pos+1:end)./(K(w_pos+1:end)-q(w_pos+1:end)-i*R(w_pos+1:end)./tanh(gamma*P*R(w_pos+1:end)*L))).^2;
% 
% G = [G_s,G_a];

% omega_1550 = 2*pi*3e8/1.55e-6;
% 
% omega= omega_1550 + w_f_f;
% 
% lambda = 2*pi*3e8./omega;
% 
% lambda=fliplr(lambda); 
% 
% G = fliplr(G);

Conver_i_anti_dB = 10*log10(Conver_i_anti);

Conver_i_stokes_dB = 10*log10(Conver_i_stokes);

figure(1)
plot(w_f_f/2/pi/1e12,Conver_i_stokes_fwm,'g',w_f_f/2/pi/1e12,Conver_i_stokes,'r')


% axis([0,30,-22,-14])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Photon pair generation rate
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


deta_v = 0.12e12;

S_i_anti = Conver_i_anti*deta_v;

S_i_stokes = Conver_i_stokes*deta_v;

S_i_stokes_fwm = Conver_i_stokes_fwm*deta_v;

S_i = [S_i;S_i_stokes];

S_i_fwm = [S_i_fwm;S_i_stokes_fwm];

figure(2)
plot(w_f_f/2/pi/1e12,S_i_anti,'g',w_f_f/2/pi/1e12,S_i_stokes,'r')

figure(21)
plot(w_f_f/2/pi/1e12,S_i_stokes,'g',w_f_f/2/pi/1e12,S_i_stokes_fwm,'r')



% figure(22) 
% plot(w_f_f/2/pi/1e12,S_i(1,:),'g',w_f_f/2/pi/1e12,S_i_fwm(1,:),'r',w_f_f/2/pi/1e12,S_i(2,:),'g',w_f_f/2/pi/1e12,S_i(3,:),'g',...
%       w_f_f/2/pi/1e12,S_i_fwm(2,:),'r',w_f_f/2/pi/1e12,S_i_fwm(3,:),'r')
%   
%   axis([0,50,0,1.9])
%   
%   legend('FWM with Re[hr]','FWM only')
%   
% ylabel('Photon pair generation rate (Photons G/s)')  
% 
% xlabel('Pump-idler frequency detuning (\omegap-\omegas)/2\pi (THz)')
  
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%Photon pair generation of Raman (at stoke side)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = 6.626068e-34;    %m2*kg/s

kb = 1.3806503e-23;  %m2*kg*s-2*K-1

T=300;

v = w_f_f/2/pi;

n_th = 1./(exp(h*v/kb/T)-1);

g_R = 2*gamma*f*h_i;

f_Raman = P*L*g_R.*(n_th+1);

S_Raman = deta_v*f_Raman;

S_Raman_p = [S_Raman_p;S_Raman]

figure(3)
plot(w_f_f/2/pi/1e12,S_Raman,'g',w_f_f/2/pi/1e12,S_i_stokes,'r')

FOM = S_i_stokes./S_Raman;

FOM_p = [FOM_p;FOM];

figure(4)
plot(w_f_f/2/pi/1e12,FOM,'g')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%Photon pair correlation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ng = 1-f+f*h_r+f*h_i*i;

P_Raman = (abs(gamma*real(ng)).^2 + abs(g_R.*(n_th+0.5)).^2)./(abs(gamma*ng).^2*P*L + abs(g_R.*(n_th+1)))./(abs(gamma*ng).^2*P*L+abs(g_R).*n_th);

P_Raman_p = [P_Raman_p;P_Raman];

figure(5)
plot(w_f_f/2/pi/1e12,P_Raman,'g')

end

figure(21)
plot(w_f_f/2/pi/1e12,S_i/1e9,'r')

figure(31)
plot(w_f_f/2/pi/1e12,S_Raman_p(1,:)/1e9,'-.r',w_f_f/2/pi/1e12,S_i(1,:)/1e9,'g',w_f_f/2/pi/1e12,S_Raman_p/1e9,'-.r',w_f_f/2/pi/1e12,S_i/1e9,'g')

legend('Raman','FWM with Re[hr]') 
ylabel('Photon pair generation rate (Photons G/s)')   
xlabel('Pump-idler frequency detuning (\omegap-\omegas)/2\pi (THz)')

% axis([0,35,0,4])

figure(41)
plot(w_f_f/2/pi/1e12,FOM_p,'r')


ylabel('Figure of Merit')   
xlabel('Pump-idler frequency detuning (\omegap-\omegas)/2\pi (THz)')


% axis([0,35,0,35])

figure(51)
plot(w_f_f/2/pi/1e12,P_Raman_p,'g')

% axis([0,35,0,280])

ylabel('Photon pair correlation')   
xlabel('Pump-idler frequency detuning (\omegap-\omegas)/2\pi (THz)')

