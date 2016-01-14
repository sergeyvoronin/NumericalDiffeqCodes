[X,T,Errors] = textread('data_files/out_constant_std.txt','%f %f %f');

x_values = unique(X);
t_values = unique(T);

Error_matrix = zeros(length(x_values), length(t_values));

my_ind = 1;
for i=1:length(x_values)
	for j=1:length(t_values)
		Error_matrix(i,j) = Errors(my_ind);
		my_ind = my_ind + 1;
	end
end

print_command = ['print -djpeg images/wave_eqn_constant1.jpg']; 
figure(1);
surf(x_values,t_values,Error_matrix');
tit_obj = title('Absolute Errors - Constant V - Std FD');
set(tit_obj,'FontSize',18);
xlabel('x values');
ylabel('t values');
zlabel('absolute error');
eval(print_command);

[X,T,Errors] = textread('data_files/out_constant_nonstd.txt','%f %f %f');

x_values = unique(X);
t_values = unique(T);

Error_matrix = zeros(length(x_values), length(t_values));

my_ind = 1;
for i=1:length(x_values)
	for j=1:length(t_values)
		Error_matrix(i,j) = Errors(my_ind);
		my_ind = my_ind + 1;
	end
end

print_command = ['print -djpeg images/wave_eqn_constant2.jpg']; 
figure(2);
surf(x_values,t_values,Error_matrix');
tit_obj=title('Absolute Errors - Constant V - Non Std FD');
set(tit_obj,'FontSize',18);
xlabel('x values');
ylabel('t values');
zlabel('absolute error');
eval(print_command);



[X,T,Errors] = textread('data_files/out_variable_std.txt','%f %f %f');

x_values = unique(X);
t_values = unique(T);

Error_matrix = zeros(length(x_values), length(t_values));

my_ind = 1;
for i=1:length(x_values)
	for j=1:length(t_values)
		Error_matrix(i,j) = Errors(my_ind);
		my_ind = my_ind + 1;
	end
end

print_command = ['print -djpeg images/wave_eqn_variable1.jpg']; 
figure(3);
surf(x_values,t_values,Error_matrix');
tit_obj=title('Absolute Errors - Variable V - Std FD');
set(tit_obj,'FontSize',18);
xlabel('x values');
ylabel('t values');
zlabel('absolute error');
eval(print_command);



[X,T,Errors] = textread('data_files/out_variable_nonstd.txt','%f %f %f');

x_values = unique(X);
t_values = unique(T);

Error_matrix = zeros(length(x_values), length(t_values));

my_ind = 1;
for i=1:length(x_values)
	for j=1:length(t_values)
		Error_matrix(i,j) = Errors(my_ind);
		my_ind = my_ind + 1;
	end
end



print_command = ['print -djpeg images/wave_eqn_variable2.jpg']; 
figure(4);
surf(x_values,t_values,Error_matrix');
tit_obj=title('Absolute Errors - Variable V - Non Std FD');
set(tit_obj,'FontSize',18);
xlabel('x values');
ylabel('t values');
zlabel('absolute error');
eval(print_command);



% plot the actual solutions ----->

[X,T,Values] = textread('data_files/out_constant_exact.txt','%f %f %f');

x_values = unique(X);
t_values = unique(T);

Value_matrix = zeros(length(x_values), length(t_values));

my_ind = 1;
for i=1:length(x_values)
	for j=1:length(t_values)
		Value_matrix(i,j) = Values(my_ind);
		my_ind = my_ind + 1;
	end
end



print_command = ['print -djpeg images/wave_eqn_constant_sol.jpg']; 
figure(5);
surf(x_values,t_values,Value_matrix');
tit_obj=title('Solution Plot - Constant V');
set(tit_obj,'FontSize',18);
xlabel('x values');
ylabel('t values');
zlabel('u(x,t) values');
eval(print_command);



[X,T,Values] = textread('data_files/out_variable_exact.txt','%f %f %f');

x_values = unique(X);
t_values = unique(T);

Value_matrix = zeros(length(x_values), length(t_values));

my_ind = 1;
for i=1:length(x_values)
	for j=1:length(t_values)
		Value_matrix(i,j) = Values(my_ind);
		my_ind = my_ind + 1;
	end
end



print_command = ['print -djpeg images/wave_eqn_variable_sol.jpg']; 
figure(6);
surf(x_values,t_values,Value_matrix');
tit_obj=title('Solution Plot - Variable V');
set(tit_obj,'FontSize',18);
xlabel('x values');
ylabel('t values');
zlabel('u(x,t) values');
eval(print_command);

