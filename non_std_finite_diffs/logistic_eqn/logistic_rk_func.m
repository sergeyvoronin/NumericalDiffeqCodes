function [t_i, approx_sol] = logistic_rk_func(t_interval, r, u0, h)
% solves the logistic equation using 
% the runge-kutta scheme
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
	
	k1 = f_u(r, uk);
	k2 = f_u(r, uk + (h/2)*k1);
	k3 = f_u(r, uk + (h/2)*k2);
	k4 = f_u(r, uk + h*k3);
	
	approx_sol(i) = uk + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
end


function ans=f_u(r,u)
	ans = r*u*(1-u);
end

end


