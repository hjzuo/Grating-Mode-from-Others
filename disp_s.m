D_te1 = [];

D_tm1 = [];

h_s = [0.2 : 0.005:0.6;];

for ii = 1: n_D
    
    D_e_t = interp1(h,te1_D(ii,:),h_s);
    
    
    D_m_t = interp1(h,tm1_D(ii,:),h_s);
    
    
    D_te1 = [D_te1;D_e_t];
    
    
    D_tm1 = [D_tm1;D_m_t];
    
    
end

D_e = [];

D_m = [];


n_h = length(h_s);


wave = [1201:1:1645];

for ii = 1:n_h
    
    D_e_t = interp1(wave_s, D_te1(:,ii),wave)';
    
    
    D_m_t = interp1(wave_s, D_tm1(:,ii),wave)';
    


    D_e = [D_e D_e_t];
    
    D_m = [D_m D_m_t];

end


n_gamma = max(size(te1_gamma))

gamma_te1 = [];

gamma_tm1 = [];

h_s = [0.2 : 0.005:0.6;];

for ii = 1: n_gamma
    
    gamma_e_t = interp1(h,te1_gamma(ii,:),h_s);
    
    
    gamma_m_t = interp1(h,tm1_gamma(ii,:),h_s);
    
    
    gamma_te1 = [gamma_te1;gamma_e_t];
    
    
    gamma_tm1 = [gamma_tm1;gamma_m_t];
    
    
end

gamma_e = [];
gamma_m = [];


n_h = length(h_s);


wave = [1201:1:1645];

for ii = 1:n_h
    
    gamma_e_t = interp1(wavelength, gamma_te1(:,ii),wave)';
    
    
    gamma_m_t = interp1(wavelength, gamma_tm1(:,ii),wave)';
    


    gamma_e = [gamma_e gamma_e_t];
    
    gamma_m = [gamma_m gamma_m_t];

end


