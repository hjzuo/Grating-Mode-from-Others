clear all;
h=160;
post.offset=0;
post.h=0;
post.w=0;
post.n=3.5;
nmodes=4;
neff_TE_array=[];neff_TM_array=[];
for i=1:length(h)

    
    post.h=h(i);
    
[neff_TE,neff_TM]=ModeConventor_neff(nmodes,post);
neff_TE_array=[neff_TE_array neff_TE(1)];
neff_TM_array=[neff_TM_array neff_TM(1)];

end