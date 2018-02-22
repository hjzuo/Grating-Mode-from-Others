load D_material
% 
% load D_wg_500nm
% 
% load D_wg_500nm_cor

[min_a, min_pos]=min(abs(wave_t-1.55*1e-6))

D_re = [];

for ii = 1:7
    
    D_t = D_TM1_ds(ii,min_pos)+D_Ge11(min_pos) - D_TM1(ii,min_pos)
    
    D_c = D_TM1_ds(ii,:)+D_Ge11-D_t;

    plot(wave_t,D_c)
    
    hold on
    
    D_re=[D_re D_c'];
    
end

