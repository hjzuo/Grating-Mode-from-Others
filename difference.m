function U_dt  = difference( u , dh);

% this function produce the differeciation of input arrays
% U_dt = du/dt

n    = length ( u );

U_dt = zeros ( size (u) );                  % U_dt 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For the middle of the u, 3:n-2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



U_n1 = u(1 : end-2 );                       % f(x-h)

U_0  = u(2 : end-1 );                       % f(x)

U_p1 = u(3 : end );                       % f(x+h)



U_dt = (U_p1 - U_n1)/2/dh;

