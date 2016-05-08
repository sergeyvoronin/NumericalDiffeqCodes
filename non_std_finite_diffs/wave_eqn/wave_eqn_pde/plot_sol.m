% solve u_tt - v^2 u_xx = 0
% v = 2
% u(x,0) = sin(w*pi*x);
% u_t(x,0) = (x+0.5)^2;
% u(0,t) = 0 = u(1,t)
close all; clear all;

w = 20;
v = 2;
%divfactors = linspace(2,0.5,3);
divfactors = [1.5];
ldivfactor = 1.5;
perrors_fd = zeros(length(divfactors),1);
perrors_nsfd = zeros(length(divfactors),1);
hs = divfactors;

for ind=1:length(divfactors)
    divfactor = divfactors(ind)
    ind
    Ufd = wave_fd1(divfactor,ldivfactor); 
    Unsfd = wave_nsfd1(divfactor,ldivfactor);
    h = 0.01/divfactor;
    l = h^2/ldivfactor;
    hs(ind) = h;

    Utrue = Ufd;
    for i=1:size(Ufd,1)
        for j=1:size(Ufd,2)
            x = h*(i-1)
            t = l*(j-1);
            val = 0.5*(sin(w*pi*(x - v*t)) + sin(w*pi*(x + v*t))) + (1/(2*v))*(2*t*v*x^2 + 2*t*v*x + 2*t^3*v^3/3 + t*v/2);
            Utrue(i,j) = val;
        end
    end

    perror1 =norm(Ufd - Utrue)/norm(Utrue)*100;
    perror2 =norm(Unsfd - Utrue)/norm(Utrue)*100;
    perrors_fd(ind) = perror1;
    perrors_nsfd(ind) = perror2;
end


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
for ind=1:15:size(Ufd,2)
    clf;
    hold on;
    plot(xvals,Ufd(:,ind),'r-*','linewidth',2);
    plot(xvals,Unsfd(:,ind),'g','linewidth',2);
    plot(xvals,Utrue(:,ind),'k--','linewidth',2);
    legend('fd','nsfd','actual');
    xlim([xvals(1),xvals(end)]);
    ylim([-1.5,1.5]);
    hold off;
    getframe;
end


return;


