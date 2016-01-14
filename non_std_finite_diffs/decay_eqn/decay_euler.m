% solves the decay equation using 
% forward euler method

clear all;

t_interval = [0:0.1:10];
h1 = 0.3;
h2 = 0.1;
u0 = 2;
lambda = 0.5;
actual_sol = zeros(length(t_interval),1);
approx_sol1 = zeros(length(t_interval),1);
approx_sol2 = zeros(length(t_interval),1);


% calculate actual sol
for i=1:length(t_interval)
	t = t_interval(i);
	actual_sol(i) = u0*exp(-lambda*(t-t_interval(1)));
end

% calculate approx sol for h1
h = h1;
N = (t_interval(length(t_interval)) - t_interval(1))/h;
t_i(1) = t_interval(1);
approx_sol1(1) = u0;
for i=2:(N+1)
	t = t_interval(1) + h*(i-1)
	t_i(i) = t;
	approx_sol1(i) = (1-lambda*h)*approx_sol1(i-1);
end


% calculate approx sol for h2
h = h2;
N = (t_interval(length(t_interval)) - t_interval(1))/h;
t_i(1) = t_interval(1);
approx_sol2(1) = u0;
for i=2:(N+1)
	t = t_interval(1) + h*(i-1)
	t_i(i) = t;
	approx_sol2(i) = (1-lambda*h)*approx_sol2(i-1);
end


figure(1);
hold on;
plot(t_interval, actual_sol, 'r', 'LineWidth', 2);
plot(t_i, approx_sol1, 'b*-');
leg_str = ['approx w/ h=', num2str(h1)];
title_str = ['Euler Method Approximation to Decay Equation'];
tit_obj = title(title_str);
set(tit_obj,'FontSize',18);
legend('actual',leg_str);
xlabel('time (t)');
ylabel('solution u(t)');
hold off;

figure(2);
hold on;
plot(t_interval, actual_sol, 'r', 'LineWidth', 2);
plot(t_i, approx_sol2, 'b*-');
leg_str = ['approx w/ h=', num2str(h2)];
title_str = ['Euler Method Approximation to Decay Equation'];
tit_obj = title(title_str);
set(tit_obj,'FontSize',18);
legend('actual',leg_str);
xlabel('time (t)');
ylabel('solution u(t)');
hold off;

