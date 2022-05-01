%bible
function [epsilon]=DrudeAgzz(lambda)

c0=2.99792458e8;
w=2*pi*c0./(lambda*1e-6);

e_inf1=3.4; w_p1=1.39e16; r1=2.70e13;


epsilon=e_inf1-w_p1^2./(w.*(w+1i*r1));

end