% solves the decay equation using 
% a nonstandard finite difference scheme

clear all;

t_interval = [0:0.1:10];
h = 1;
u0 = 2;
lambda = 0.5;

N = (t_interval(length(t_interval)) - t_interval(1))/h;

% calculate approx sol
t_i(1) = t_interval(1);
approx_sol(1) = u0;
for i=2:(N+1)
	t = t_interval(1) + h*(i-1)
	t_i(i) = t;
	approx_sol(i) = approx_sol(i-1)*exp(-lambda*h);
end

% calculate actual sol
for i=1:length(t_interval)
	t = t_interval(i);
	actual_sol(i) = u0*exp(-lambda*(t-t_interval(1)));
end


figure(1);
hold on;
plot(t_interval, actual_sol, 'r', 'LineWidth', 2);
plot(t_i, approx_sol, 'b*-');
leg_str = ['approx w/ h=', num2str(h)];
title_str = ['NSFD Method Approximation to Decay Equation'];
tit_obj = title(title_str);
set(tit_obj,'FontSize',18);
legend('actual',leg_str);
xlabel('time (t)');
ylabel('solution u(t)');
hold off;


