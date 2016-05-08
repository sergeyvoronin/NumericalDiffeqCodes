function U = wave_nsfd1(divfactor,ldivfactor) 

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
    xm1 = h*(i-2);
    xp1 = h*i;

    gi = (x+0.5)^2;

    v = abs(cos(x)) + 1;
    vm1 = abs(cos(xm1)) + 1;
    vp1 = abs(cos(xp1)) + 1;

    ki = w/v;
    kim1 = w/vm1;
    kip1 = w/vp1;


    %val = (sin(w*l/2))^2/(sin(ki*h/2))^2;
    pival = (cos(w*l)-1)/(cos(ki*h)-1);
    pivalp1 = (cos(w*l)-1)/(cos(kip1*h)-1);
    pivalm1 = (cos(w*l)-1)/(cos(kim1*h)-1);
    pivalp1 = pival;
    pivalm1 = pival;

    
    %U(i,2) = fi + gi*l + 0.5*(pivalp1*fip1 - 2*pival*fi + pivalm1*fim1);
    U(i,2) = 0.5*(2*(1-pival)*U(i,1) + pivalp1*U(i+1,1) + pivalm1*U(i-1,1) + 2*gi*l);
end

% set the other values based on explicit fd scheme
for j=2:Nj
    for i=2:(Ni-1)
        x = h*(i-1);
        xm1 = h*(i-2);
        xp1 = h*i;
        
        v = abs(cos(x)) + 1;
        vm1 = abs(cos(xm1)) + 1;
        vp1 = abs(cos(xp1)) + 1;

        ki = w/v;
        kim1 = w/vm1;
        kip1 = w/vp1;


        %val = (sin(w*l/2))^2/(sin(ki*h/2))^2;
        pival = (cos(w*l)-1)/(cos(ki*h)-1);
        pivalp1 = (cos(w*l)-1)/(cos(kip1*h)-1);
        pivalm1 = (cos(w*l)-1)/(cos(kim1*h)-1);
        pivalp1 = pival;
        pivalm1 = pival;


        U(i,j+1) = 2*(1-pival)*U(i,j) + pivalp1*U(i+1,j) + pivalm1*U(i-1,j) - U(i,j-1);
        U(1,j+1) = 0;
        U(Ni,j+1) = 0;
    end
end

end

