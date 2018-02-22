plot(rw_n(5:end-20),Aeff_TM_05(5:end-20),'r',rw_n(3:end-20),Aeff_TM_12(3:end-20),'g',rw_n(4:end-20),Aeff_TM_17(4:end-20),'c')
legend('Index contrast 0.5','Index contrast 1.2','Index contrast 1.7')

xlabel('Waveguide Dimention(\mum)')
ylabel('Effective Mode area(\mum^2)')