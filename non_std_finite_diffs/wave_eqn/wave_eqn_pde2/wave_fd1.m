function U = wave_fd1(divfactor,ldivfactor) 

h = 0.01/divfactor;
l = h^2/ldivfactor;
L = 1;
T = 3;
Ni = round(L/h)+1;
Nj = round(T/l)+1;
U = zeros(Ni,Nj);
w = 20;

% set U(i,1), U(i,2) based on initial conditions
for i=2:(Ni-1)
    x = h*(i-1);
    fi = sin(w*pi*x);
    U(i,1) = fi;
end

for i=2:(Ni-1)
    x = h*(i-1);
    gi = (x+0.5)^2;

    x = h*(i-1);
    xm1 = h*(i-2);
    xp1 = h*(i);
    v = abs(cos(x)) + 1;
    vip1 = abs(cos(xp1)) + 1;
    vim1 = abs(cos(xm1)) + 1;

    si = v^2 * l^2/h^2; 
    sip1 = si;
    sim1 = si;
    
    %U(i,2) = fi + gi*l + 0.5*(sip1*fip1 - 2*si*fi + sim1*fim1);
    U(i,2) = 0.5*(2*(1-si)*U(i,1) + sip1*U(i+1,1) + sim1*U(i-1,1) + 2*gi*l);
end

% set the other values based on explicit fd scheme
for j=2:Nj
    for i=2:(Ni-1)
        x = h*(i-1);
        xm1 = h*(i-2);
        xp1 = h*(i);
        v = abs(cos(x)) + 1;
        vip1 = abs(cos(xp1)) + 1;
        vim1 = abs(cos(xm1)) + 1;
        
        si = v^2 * l^2/h^2; 
        sip1 = si;
        sim1 = si;
    

        U(i,j+1) = 2*(1-si)*U(i,j) + sip1*U(i+1,j) + sim1*U(i-1,j) - U(i,j-1);
        U(1,j+1) = 0;
        U(Ni,j+1) = 0;
    end
end

end

