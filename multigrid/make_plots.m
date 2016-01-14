% load data
[x,y,u_mg,u_act] = textread('data/output.txt');

problem_num = 2

num_unq = length(unique(x));
h = 1/num_unq;
k = 1/num_unq;
x_vals = repmat(0,num_unq,num_unq);
y_vals = repmat(0,num_unq,num_unq);
u_act_vals = repmat(0,num_unq,num_unq);
u_mg_vals = repmat(0,num_unq,num_unq);
u_errors_vals = repmat(0,num_unq,num_unq);

ind = 1;
for i=1:num_unq
    for j=1:num_unq
        x_vals(i,j) = x(ind);
        y_vals(i,j) = y(ind);
        u_act_vals(i,j) = u_act(ind);
        u_mg_vals(i,j) = u_mg(ind);
        u_errors_vals(i,j) = abs(u_act(ind) - u_mg(ind));
        ind = ind + 1;
    end
end
 
    % set up contour vector
    if problem_num == 1
        contour_vec = linspace(0, 1, 100);
    else
        contour_vec = linspace(0, 2, 100);
    end

    % make some plots --->
    figure(1);
    hold on;
    tit = title('Approximate Solution');
    set(tit,'FontSize', 22);
    contourf(x_vals,y_vals,u_mg_vals, contour_vec);
    xl = xlabel('x values');
    yl = ylabel('y values');
    set(xl,'FontSize',20);
    set(yl,'FontSize',20);
    colorbar;
    set(gca,'FontSize',18);
    hold off;
    file_name = ['images/problem', num2str(problem_num), '_approx_solution_mg_h_', strrep(num2str(h),'.','_'), '_k_', strrep(num2str(k),'.','_'), '.jpg'];
    print_cmd = ['print -djpeg ', file_name];
    eval(print_cmd);


    figure(2)
    hold on;
    tit = title('Actual Solution');
    set(tit,'FontSize', 22);
    contourf(x_vals,y_vals,u_act_vals, contour_vec);
    xl = xlabel('x values');
    yl = ylabel('y values');
    set(xl,'FontSize',20);
    set(yl,'FontSize',20);
    colorbar;
    set(gca,'FontSize',18);
    hold off;
    file_name = ['images/problem', num2str(problem_num), '_actual_solution_h_', strrep(num2str(h),'.','_'), '_k_', strrep(num2str(k),'.','_'), '.jpg'];
    print_cmd = ['print -djpeg ', file_name];
    eval(print_cmd);


    contour_vec = linspace(min(min(u_errors_vals)), max(max(u_errors_vals)), 100);

    figure(3)
    hold on;
    tit = title('Absolute Errors');
    set(tit,'FontSize', 22);
    contourf(x_vals,y_vals,u_errors_vals, contour_vec);
    xl = xlabel('x values');
    yl = ylabel('y values');
    set(xl,'FontSize',20);
    set(yl,'FontSize',20);
    colorbar;
    set(gca,'FontSize', 18);
    hold off;
    file_name = ['images/problem', num2str(problem_num), '_abs_errors_h_', strrep(num2str(h),'.','_'), '_k_', strrep(num2str(k),'.','_'), '.jpg'];
    print_cmd = ['print -djpeg ', file_name];
    eval(print_cmd);


% close all;
% post process --->
%system('./post_process.pl');

