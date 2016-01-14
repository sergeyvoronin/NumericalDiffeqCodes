% solves the logistic equation using 
% a nonstandard finite difference scheme
%
% du/dt = r*u(1-u);  u(0) = 2;


clear all;

t_interval = [0:0.1:10];
h = .2;
u0 = 5;
r = 0.5;

N = (t_interval(length(t_interval)) - t_interval(1))/h;

% calculate approx sol
t_i(1) = t_interval(1);
approx_sol(1) = u0;
for i=2:(N+1)
	t = t_interval(1) + h*(i-1)
	t_i(i) = t;
	uk = approx_sol(i-1);
	approx_sol(i) = ( uk*(2-exp(-r*h)) )/( 1 + uk*(1-exp(-r*h)) );
end

% calculate actual sol
for i=1:length(t_interval)
	t = t_interval(i);
	actual_sol(i) = 1/( 1 + (1/u0 - 1)*exp(-r*t)  );
end


figure(1);
command_str = ['print -djpeg ', 'logistic_ode1.jpg'];
hold on;
plot(t_interval, actual_sol, 'r', 'LineWidth', 2);
plot(t_i, approx_sol, 'b*-');
leg_str = ['approx w/ h=', num2str(h)];
title_str = ['NSFD Method Approximation to Logistic Equation'];
tit_obj = title(title_str);
set(tit_obj,'FontSize',18);
legend('actual',leg_str);
xlabel('time (t)');
ylabel('solution u(t)');
hold off;
eval(command_str);


