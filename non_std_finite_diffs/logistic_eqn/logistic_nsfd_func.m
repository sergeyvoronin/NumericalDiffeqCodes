function [t_i, approx_sol] = logistic_nsfd_func(t_interval, r, u0, h)
% solves the logistic equation using 
% a nonstandard finite difference scheme
%
% du/dt = r*u(1-u);  u(0) = u0;

N = (t_interval(length(t_interval)) - t_interval(1))/h;

% calculate approx sol
t_i(1) = t_interval(1);
approx_sol(1) = u0;
for i=2:(N+1)
	t = t_interval(1) + h*(i-1);
	t_i(i) = t;
	uk = approx_sol(i-1);
	approx_sol(i) = ( uk*(2-exp(-r*h)) )/( 1 + uk*(1-exp(-r*h)) );
end

