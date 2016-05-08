% solve discretized version of u_tt - v^2 u_xx = 0
% v = abs(cos(x)) + 1;
% u(x,0) = sin(w*pi*x);
% u_t(x,0) = (x+0.5)^2;
% u(0,t) = 0 = u(1,t)
close all; clear all;

w = 20;
%divfactors = linspace(2,0.5,3);
divfactors = [0.4];
ldivfactor = 1;
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
end


figure(2);
for ind=1:20:size(Ufd,2)
    clf;
    hold on;
    plot(xvals,Ufd(:,ind),'r-*','linewidth',2);
    plot(xvals,Unsfd(:,ind),'g','linewidth',2);
    legend('fd','nsfd');
    xlim([xvals(1),xvals(end)]);
    ylim([-2,2]);
    xlabel('x value');
    ylabel('u(x,t) value');
    hold off;
    getframe;
end



return;


