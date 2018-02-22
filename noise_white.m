function  u = noise_white(P,T);





n_t = length(T);

u = zeros(size(T));

u(n_t/2+1)=sqrt(P);


