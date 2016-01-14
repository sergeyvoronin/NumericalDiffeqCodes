function actual_sol = logistic_eqn_solution(r, u0, t)
% solves the logistic equation using 
% the analytical solution
%
% du/dt = r*u(1-u);  u(0) = u0;
actual_sol = 1/( 1 + (1/u0 - 1)*exp(-r*t)  );


