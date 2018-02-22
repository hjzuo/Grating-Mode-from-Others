function image_all_fields(nmode, lambda, neff, Hx, Hy, dx, dy, eps,n,xc,yc,x,y)





Hxc_TE     = [];

Hyc_TE     = [];

Hzc_TE     = [];

Ezc_TE     = [];

Exc_TE     = [];

Eyc_TE     = [];

Ic_TE      = [];

Aeff_TE    = [];

P_upper_TE = []; 

P_wg_TE    = [];

P_lower_TE = [];


Hxc_TM     = [];

Hyc_TM     = [];

Hzc_TM     = [];

Ezc_TM     = [];

Exc_TM     = [];

Eyc_TM     = [];

Ic_TM      = [];

Aeff_TM    = [];

P_upper_TM = []; 

P_wg_TM    = [];

P_lower_TM = [];

TM_n = 0;

TE_n = 0;

for ii = 1 : nmode
    
    



% [Hxc_t, Hyc_t, Hzc_t, Ezc_t, Exc_t, Eyc_t]  = mode_all(lambda, neff(ii), Hx(:,:,ii), Hy(:,:,ii), dx, dy, eps);          % generate all the field components from Hx, Hy

[Hxc_t,Hyc_t,Hzc_t,Exc_t,Eyc_t,Ezc_t] = postprocess (lambda, neff(ii), Hx(:,:,ii), Hy(:,:,ii), dx, dy, eps, '0000');


% [Aeff_t, Ic_t,P_upper_t,P_wg_t,P_lower_t] = wg_area_eff(Hxc_t,Hyc_t,Hzc_t,Exc_t,Eyc_t,Ezc_t,n,dx,dy,eps);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   field Components display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  figure(1),con_mode(x,y,Hx,1,3,60,'Hx TM mode');
figure((2*ii - 1))
colormap(jet(256));


hn1 = max(abs(Hyc_t(:)));

hn2 = max(abs(Hxc_t(:)));

if hn1>hn2
    
    hn = hn1;
    
    mode = 'TE';
    
    TE_n = TE_n +1;
    
    tt = TE_n;
    
Hxc_TE     = [Hxc_TE, Hxc_t];

Hyc_TE     = [Hyc_TE, Hyc_t];

Hzc_TE     = [Hzc_TE, Hzc_t];

Ezc_TE     = [Ezc_TE, Ezc_t];

Exc_TE     = [Exc_TE, Exc_t];

Eyc_TE     = [Eyc_TE, Eyc_t];

% Ic_TE      = [Ic_TE, Ic_t];
% 
% Aeff_TE    = [Aeff_TE, Aeff_t];
% 
% P_upper_TE = [P_upper_TE, P_upper_t]; 
% 
% P_wg_TE    = [P_wg_TE, P_wg_t];
% 
% P_lower_TE = [P_lower_TE,P_lower_t];
    
else
    
    hn = hn2;
    
    mode = 'TM';
    
    TM_n = TM_n +1;
    
    tt = TM_n;
Hxc_TM     = [Hxc_TM, Hxc_t];

Hyc_TM     = [Hyc_TM, Hyc_t];

Hzc_TM     = [Hzc_TM, Hzc_t];

Ezc_TM     = [Ezc_TM, Ezc_t];

Exc_TM     = [Exc_TM, Exc_t];

Eyc_TM     = [Eyc_TM, Eyc_t];

% Ic_TM      = [Ic_TM, Ic_t];
% 
% Aeff_TM    = [Aeff_TM, Aeff_t];
% 
% P_upper_TM = [P_upper_TM, P_upper_t]; 
% 
% P_wg_TM    = [P_wg_TM, P_wg_t];
% 
% P_lower_TM = [P_lower_TM,P_lower_t];
end

fprintf(1,'-----------%c%c%d MODE ----------\n\n',mode,tt);
    
fprintf(1,'neff (%c%c%d mode)= %7.5f\n',mode,tt,neff(ii));

en = hn/n(2);


subplot(2,3,1),image_mode(xc,yc,Hxc_t/hn,sprintf('Hx (%c%c%d mode)',mode,tt));
hold on;
con_mode(xc,yc,Hxc_t/hn,1,3,55,sprintf('Hx (%c%c%d mode)',mode,tt));
hold off;

% v = xlim;

% line(v,[0,0],'color','k');
% 
% line([v(1),-rw_half,-rw_half,rw_half,rw_half,v(2)],[h(2)-rh,h(2)-rh,h(2),h(2),h(2)-rh,h(2)-rh],'color','k');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2),con_mode(x,y,Hy,1,3,60,'Hy TM mode');

subplot(2,3,2),image_mode(xc,yc,Hyc_t/hn,sprintf('Hy (%c%c%d mode)',mode,tt));
hold on;
con_mode(xc,yc,Hyc_t/hn,1,3,55,sprintf('Hy (%c%c%d mode)',mode,tt));
hold off;

% v = xlim;
% 
% line(v,[0,0],'color','k');
% 
% line([v(1),-rw_half,-rw_half,rw_half,rw_half,v(2)],[h(2)-rh,h(2)-rh,h(2),h(2),h(2)-rh,h(2)-rh],'color','k');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,3),image_mode(xc,yc,Hzc_t/hn,sprintf('Hz (%c%c%d mode)',mode,tt));
hold on;
con_mode(xc,yc,Hzc_t/hn,1,3,55,sprintf('Hz (%c%c%d mode)',mode,tt));
hold off;

% v = xlim;
% 
% line(v,[0,0],'color','k');
% 
% line([v(1),-rw_half,-rw_half,rw_half,rw_half,v(2)],[h(2)-rh,h(2)-rh,h(2),h(2),h(2)-rh,h(2)-rh],'color','k');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,4),image_mode(xc,yc,Exc_t/en,sprintf('Ex (%c%c%d mode)',mode,tt));
hold on;
con_mode(xc,yc,Exc_t/en,1,3,55,sprintf('Ex (%c%c%d mode)',mode,tt));
hold off;

% v = xlim;
% 
% line(v,[0,0],'color','k');
% 
% line([v(1),-rw_half,-rw_half,rw_half,rw_half,v(2)],[h(2)-rh,h(2)-rh,h(2),h(2),h(2)-rh,h(2)-rh],'color','k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,5),image_mode(xc,yc,Eyc_t/en,sprintf('Ey (%c%c%d mode)',mode,tt));
hold on;
con_mode(xc,yc,Eyc_t/en,1,3,55,sprintf('Ey (%c%c%d mode)',mode,tt));
hold off;

% v = xlim;
% 
% line(v,[0,0],'color','k');
% 
% line([v(1),-rw_half,-rw_half,rw_half,rw_half,v(2)],[h(2)-rh,h(2)-rh,h(2),h(2),h(2)-rh,h(2)-rh],'color','k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,6),image_mode(xc,yc,Ezc_t/en,sprintf('Ez (%c%c%d mode)',mode,tt));
hold on;
con_mode(xc,yc,Ezc_t/en,1,3,55,sprintf('Ez (%c%c%d mode)',mode,tt));
hold off;

% v = xlim;
% 
% line(v,[0,0],'color','k');
% 
% line([v(1),-rw_half,-rw_half,rw_half,rw_half,v(2)],[h(2)-rh,h(2)-rh,h(2),h(2),h(2)-rh,h(2)-rh],'color','k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The InTMnsity : I = 1/4(E x H* + E* x H) = 1/4 (ExHy*-EyHx* +Ex*Hy-Ey*Hx)
%
% CalculaTM the inTMnsity and effective mode area and the distribution of
% the Power
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% fprintf(1,'Aneff (%c%c%d Mode) = %7.5f um^2\n',mode,tt,Aeff_t);
% 
% fprintf(1,'Power in upper clad (%c%c%d Mode) = %7.5f \n',mode, tt,P_upper_t);
% 
% fprintf(1,'Power in waveguide (%c%c%d Mode) = %7.5f \n',mode,tt,P_wg_t);
% 
% fprintf(1,'Power in lower clad (%c%c%d Mode) = %7.5f \n',mode,tt,P_lower_t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% figure((2*ii));
% 
% colormap(jet(256));
% 
% Ih = max(abs(Ic_t(:)));
% 
% image_mode(xc,yc,Ic_t/Ih,sprintf('Intensity (%c%c%d mode)',mode,tt));
% hold on;
% con_mode(xc,yc,Ic_t/Ih,1,3,55,sprintf('Intensity (%c%c%d mode)',mode,tt));
% hold off;

% v = xlim;
% 
% line(v,[0,0],'color','k');
% 
% line([v(1),-rw_half,-rw_half,rw_half,rw_half,v(2)],[h(2)-rh,h(2)-rh,h(2),h(2),h(2)-rh,h(2)-rh],'color','k');
% 
% fprintf('----------------------\n\n\n');

end




