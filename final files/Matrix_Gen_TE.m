
function [Q,V,W]=Matrix_Gen_TE(e_m, e_d, f1, f2, Period, e1,lambda, theta, Num_ord)
    ordMin=-Num_ord;
    ordMax=Num_ord;
    ordDif=2*Num_ord+1;
    i = [ordMin:ordMax];    %-40 to 40
    k0 = 2*pi/lambda;
    kxi = k0*(sqrt(e1)*sin(theta)+i*lambda/Period); %kxi (1_by_i vector)
    Kx2 = diag(power(kxi/k0,2)); %diagonal matrix with kxi/k0 (#i_by_#i matrix)
                
    %Fourier expansion of the dielectric function in grating
    %region: positive (1 to 80)
    h = [1:ordMax-ordMin];
    epsilonp = j*(exp(-j*h*2*pi*f2)-exp(-j*h*2*pi*f1))./(2*pi*h) *(e_m-e_d);
 
    %0th order
    epsilon0 = (f2-f1)*(e_m-e_d)+e_d;
 
    %Fourier expansion of the dielectric function in grating
    %region: negative (-80 to -1)
    h = [ordMin-ordMax:-1];
    epsilonn = j*(exp(-j*h*2*pi*f2)-exp(-j*h*2*pi*f1))./(2*pi*h) *(e_m-e_d);
 
    %Fourier expansion of the dielectric function in grating region
    epsilonG = [epsilonn epsilon0 epsilonp];    %(-80 to -1, 0, 1 to 80)

    %matrix of the dielectric function of grating based on Fourier
    %components
    for i=1:ordMax-ordMin+1
        E(i,:) = epsilonG(ordMax-ordMin+2-i:ordMax-ordMin+2-i+2*Num_ord);
    end
    %Note E(1,:)=[0, ..., 80]
    %     E(2,:)=[-1, ..., 79]
    %     E(81,:)=[-80, ..., 0]

    I = eye(ordDif); %Unit matrix
    A = Kx2-E;
    [W,Q] = eig(A);   %Q:eigenvalue, W:corresponding vector
    Q = sqrt(Q);
    V = W*Q;

    