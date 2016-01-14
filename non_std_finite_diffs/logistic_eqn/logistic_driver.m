% compares srd and non std methods for the solution of the logistic equation
%
% du/dt = r*u(1-u);  u(0) = u0;
%

clear all;

t_interval = [0:0.1:10];
u0 = 2;
h_set = [0.1:0.1:1];
r_set = [0.1:0.1:2];

nsfd_errors = zeros(length(r_set), length(h_set));
rk_errors = zeros(length(r_set), length(h_set));
for r_ind=1:length(r_set)
	for h_ind=1:length(h_set)
		r = r_set(r_ind);
		h = h_set(h_ind);
		fprintf('(r = %f, h = %f)\n', r, h);
		
		[t_i, nsfd_solution]  = logistic_nsfd_func(t_interval,r,u0,h);
		actual_solution = [];
		for n=1:length(t_i)
			actual_solution(n) = logistic_eqn_solution(r,u0,t_i(n));
		end

		nsfd_error = abs(actual_solution - nsfd_solution);
		nsfd_errors(r_ind, h_ind) = norm(nsfd_error, 2);


		[t_i, rk_solution]  = logistic_rk_func(t_interval,r,u0,h);
		rk_error = abs(actual_solution - rk_solution);
		rk_errors(r_ind, h_ind) = norm(rk_error, 2);

	end
end


R = meshgrid(r_set);
H = meshgrid(h_set);


figure(1);
print_command = ['print -djpeg images/logistic_rk.jpg']; 
surf(r_set, h_set, rk_errors');
ylim([0,1.5]);
title_str =['Error Norms for Runge-Kutta Algorithm'];
tit_obj = title(title_str);
set(tit_obj,'FontSize',18);
xlabel('value of r');
ylabel('value of h');
zlabel('error norm');
eval(print_command);

figure(2);
print_command = ['print -djpeg images/logistic_nsfd.jpg']; 
surf(r_set, h_set, nsfd_errors');
ylim([0,1.5]);
title_str =['Error Norms for NSFD Algorithm'];
tit_obj = title(title_str);
set(tit_obj,'FontSize',18);
xlabel('value of r');
ylabel('value of h');
zlabel('error norm');
eval(print_command);

