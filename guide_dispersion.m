clear
load D_material

load D_wg_500nm

load D_wg_500nm_cor

[min_1550, pos_1550]=min(abs(wave_t-1.55e-6));

for ii = 1:4

D_t = D_TM1_ds(ii,:) + D_Ge11;

Ds = D_TM1(ii,pos_1550)-D_t(pos_1550);


figure(3)

plot(wave_t,D_t+Ds,'r',wave_t,D_TM1(ii,:),'g');

hold on
end
hold off


