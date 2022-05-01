%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is the main program of RCWA method to calculate the spectral directional radiative properties for TE and TM incidence 
% Dr. Zhuomin Zhang's group at Georgia Tech
% Last modified by Bo Zhao (September 2014)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; clc;
c0 = 299792458;                             % speed of light in space, [m/s]
theta =0;                                           % angle of incidence, [rad]
wn = linspace(3500,25000,300) ; % wavenumber, [1/cm]
lambda =1e4./wn;                           % wavelength, [um]
w = 1e6*2*pi*c0./lambda;            % angular frequency, [rad/s]
d = [0.02 0.022];                                           % thickness of each layer from front to back, [um]
N = length(d);                                   % # of layers
Period(1:N) = 1.12;                            % Period of gratings for each layer, [um]
width = [0 0.07];                                    % width of metal strips, [um]
psi = width/Period(1);                     % filling ratio of grating layer
f1 = [0 0];                                              % normalized position for left-end of metal strip
f2 = [psi];                                           % normalized position for right-end of metal strip
Num_ord = 19;                                  % number for the highest diffraction order
for ind = 1:length(lambda)
    % Incidence medium
      e(1) = 1.33;                                                 % Usually is air or vacuum 
    % Layered structure
      
   %   e_m(1) = Palik_SiO2(lambda(ind));   % Ridge material (metal)  
    %  e_d(1) = Palik_SiO2(lambda(ind));
     % e_m(2:N) = Palik_Au(lambda(ind));   % Ridge material (metal)  
      % e_d(2:N) = Palik_SiO2(lambda(ind));                       % Groove material (air)
      
      e_m(1:N) = Palik_Au(lambda(ind));   % Ridge material (metal)  
      e_d(1:N) = Palik_SiO2(lambda(ind));
      
    %Substrate
      e(2)=Palik_SiO2(lambda(ind));                                                    %  air or opaque substrate
    %==========================================
    
    [Ref(ind), Tran(ind)] = RCWA_Multi_TE(N, e_m, e_d, f1, f2, Period, d, e, lambda(ind), theta, Num_ord); 

end
figure;
plot(wn,Ref,wn,Tran,wn,1-Ref-Tran);
legend('Ref','Tran','\alpha');
xlim([3500,25000]);
xlabel('Wavenumber, \nu (cm^-^1)');
ylabel('R,T,\alpha');
