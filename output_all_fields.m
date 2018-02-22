function [Aeff_TE, P_upper_TE, P_wg_TE, P_lower_TE, neff_TE, Aeff_TM, P_upper_TM, P_wg_TM, P_lower_TM, neff_TM] =  output_all_fields(nmode, lambda, neff, Hx, Hy, dx, dy, eps,h,n,rh,rw_half,side_half,xc,yc,jj,neff_last_te,neff_last_tm)




% 
% Hxc_TE     = [];
% 
% Hyc_TE     = [];
% 
% Hzc_TE     = [];
% 
% Ezc_TE     = [];
% 
% Exc_TE     = [];
% 
% Eyc_TE     = [];
% 
% Ic_TE      = [];

Aeff_TE    = [];

P_upper_TE = []; 

P_wg_TE    = [];

P_lower_TE = [];

neff_TE    = [];


% Hxc_TM     = [];
% 
% Hyc_TM     = [];
% 
% Hzc_TM     = [];
% 
% Ezc_TM     = [];
% 
% Exc_TM     = [];
% 
% Eyc_TM     = [];
% 
% Ic_TM      = [];

Aeff_TM    = [];

P_upper_TM = []; 

P_wg_TM    = [];

P_lower_TM = [];

neff_TM    = [];

TM_n = 0;

TE_n = 0;

% if jj == 1

for ii = 1 : nmode
    
    



[Hxc_t, Hyc_t, Hzc_t, Ezc_t, Exc_t, Eyc_t]  = mode_all(lambda, neff(ii), Hx(:,:,ii), Hy(:,:,ii), dx, dy, eps);          % generate all the field components from Hx, Hy

[Aeff_t, Ic_t,P_upper_t,P_wg_t,P_lower_t]   = wg_area_eff(Hxc_t,Hyc_t,Exc_t,Eyc_t,n,dx,dy,eps);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   field Components display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


hn1 = max(abs(Hyc_t(:)));

hn2 = max(abs(Hxc_t(:)));

if hn1>hn2
    
    hn = hn1;
    
    mode = 'TE';
    
    TE_n = TE_n +1;
    
    tt = TE_n;
    
%     Hxc_TE     = [Hxc_TE, Hxc_t];
% 
%     Hyc_TE     = [Hyc_TE, Hyc_t];
% 
%     Hzc_TE     = [Hzc_TE, Hzc_t];
% 
%     Ezc_TE     = [Ezc_TE, Ezc_t];
% 
%     Exc_TE     = [Exc_TE, Exc_t];
% 
%     Eyc_TE     = [Eyc_TE, Eyc_t];
% 
%     Ic_TE      = [Ic_TE, Ic_t];

    Aeff_TE    = [Aeff_TE, Aeff_t];

    P_upper_TE = [P_upper_TE, P_upper_t]; 

    P_wg_TE    = [P_wg_TE, P_wg_t];

    P_lower_TE = [P_lower_TE,P_lower_t];

    neff_TE =[neff_TE, neff(ii)];
    
else
    
    hn = hn2;
    
    mode = 'TM';
    
    TM_n = TM_n +1;
    
    tt = TM_n;
%     Hxc_TM     = [Hxc_TM, Hxc_t];
% 
%     Hyc_TM     = [Hyc_TM, Hyc_t];
% 
%     Hzc_TM     = [Hzc_TM, Hzc_t];
% 
%     Ezc_TM     = [Ezc_TM, Ezc_t];
% 
%     Exc_TM     = [Exc_TM, Exc_t];
% 
%     Eyc_TM     = [Eyc_TM, Eyc_t];
% 
%     Ic_TM      = [Ic_TM, Ic_t];

    Aeff_TM    = [Aeff_TM, Aeff_t];

    P_upper_TM = [P_upper_TM, P_upper_t]; 

    P_wg_TM    = [P_wg_TM, P_wg_t];

    P_lower_TM = [P_lower_TM,P_lower_t];
    
    neff_TM    = [neff_TM, neff(ii)];

end
    
    
    
end

% % 
% % else
% %     
% %     for ii = 1 : nmode
% %     
% %     
% % 
% % 
% % 
% %    [Hxc_t, Hyc_t, Hzc_t, Ezc_t, Exc_t, Eyc_t]  = mode_all(lambda, neff(ii), Hx(:,:,ii), Hy(:,:,ii), dx, dy, eps);          % generate all the field components from Hx, Hy
% % 
% %    [Aeff_t, Ic_t,P_upper_t,P_wg_t,P_lower_t]   = wg_area_eff(Hxc_t,Hyc_t,Exc_t,Eyc_t,n,dx,dy,eps);
% %    
% %    
% %    [min_te, pos_te] = min(abs(neff_TE - neff(ii)));
% %    
% %    [min_tm, pos_tm] = min(abs(neff_TM - neff(ii)));
% %    
% %    if abs(min_te) < abs(min_tm)
% %    
% %       
% %     
% %          Aeff_TE    = [Aeff_TE, Aeff_t];
% % 
% %          P_upper_TE = [P_upper_TE, P_upper_t]; 
% % 
% %          P_wg_TE    = [P_wg_TE, P_wg_t];
% % 
% %          P_lower_TE = [P_lower_TE,P_lower_t];
% % 
% %          neff_TE =[neff_TE, neff(ii)];
% % 





