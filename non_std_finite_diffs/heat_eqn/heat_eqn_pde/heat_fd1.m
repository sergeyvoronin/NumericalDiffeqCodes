function U = heat_fd1(divfactor,ldivfactor) 

h = 0.01/divfactor;
l = h^2/ldivfactor;


L = 1;
T = 1;
Ni = round(L/h) + 1;
Nj = round(T/l) + 1;
k = 1; % constant k
U = zeros(Ni,Nj);

% set U(i,1) based on initial condition
%for i=2:(Ni-1)
for i=1:(Ni)
    x = h*(i-1);
    val = sin(pi*x) + x;
    vale = 0 + x*(1 - 0)/(L);
    fi = val - vale;
    U(i,1) = fi;
end
U(1,1) = 0; U(Ni, 1) = 0;


for j=1:(Nj)
    U(1,j) = 0; U(Ni, j) = 0;
    for i=2:(Ni-1)
        si = l/h^2; 
        %si = 0.5*((exp(-pi^2*k*l) - 1)/(cos(pi*h)-1)); 

        U(i,j+1) = (1-2*si)*U(i,j) + si*(U(i+1,j) + U(i-1,j));
    end
end

end


