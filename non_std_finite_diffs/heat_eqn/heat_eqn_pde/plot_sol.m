% solve u_t = u_xx on x,t \in [0,1]
% u(x,0) = sin(pi*x) + x
% u(0,t) = 0; u(1,t) = 1;
close all; clear all;

%divfactors = linspace(1.2,0.15,15);
divfactors = [0.05];
ldivfactor = 1.5;
perrors_fd = zeros(length(divfactors),1);
perrors_nsfd = zeros(length(divfactors),1);
hs = divfactors;

for ind=1:length(divfactors)
    divfactor = divfactors(ind)
    ind
    Ufd = heat_fd1(divfactor,ldivfactor);
    Unsfd = heat_nsfd1(divfactor,ldivfactor);
    h = 0.01/divfactor;
    l = h^2/ldivfactor;
    hs(ind) = h;
    
    Utrue = Ufd;
    for i=1:size(Ufd,1)
        for j=1:size(Ufd,2)
            x = h*(i-1)
            t = l*(j-1);
            v = 2;
            w = 8;
            Ufd(i,j) = Ufd(i,j) + x;
            Unsfd(i,j) = Unsfd(i,j) + x;
            val = x + exp(-pi^2*t)*sin(pi*x);
            Utrue(i,j) = val;
        end
    end

    perror1 =norm(Ufd - Utrue)/norm(Utrue)*100;
    perror2 =norm(Unsfd - Utrue)/norm(Utrue)*100;
    perrors_fd(ind) = perror1;
    perrors_nsfd(ind) = perror2;
end


%% uncomment if run over multiple h values
%figure(1);
%title('ERRORS vs increasing step size');
%hold on;
%plot(hs,perrors_fd,'r','linewidth',2);
%plot(hs,perrors_nsfd,'g','linewidth',2);
%xlim([hs(1),hs(end)]);
%xlabel('h');
%ylabel('reconstruction % error');
%legend('fd','nsfd');
%hold off;
%
%return;
%

xvals = zeros(size(Ufd,1),1);
for i=1:size(Ufd,1)
    xvals(i) = (i-1)*h;
end

figure(2);
for ind=1:1:size(Ufd,2)
    clf;
    hold on;
    title('u(x,t_f) from fd, nsfd, actual solution');
    plot(xvals,Ufd(:,ind),'r-*','linewidth',2);
    plot(xvals,Unsfd(:,ind),'g','linewidth',2);
    plot(xvals,Utrue(:,ind),'k--','linewidth',2);
    xlabel('x');
    ylabel('u(x,t_f)');
    legend('fd','nsfd','actual');
    xlim([xvals(1),xvals(end)]);
    ylim([0,1.2]);
    hold off;
    getframe;
end

return;

