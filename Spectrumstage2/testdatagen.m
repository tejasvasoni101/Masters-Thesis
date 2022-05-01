%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is the main program of RCWA method to calculate the spectral directional radiative properties for TE and TM incidence 
% Dr. Zhuomin Zhang's group at Georgia Tech
% Last modified by Bo Zhao (September 2014)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;
c0 = 299792458;                             % speed of light in space, [m/s]
theta =0;                                           % angle of incidence, [rad]
wn = linspace(5000,10000,100) ; % wavenumber, [1/cm]
lambda =1e4./wn;                           % wavelength, [um]
w = 1e6*2*pi*c0./lambda;            % angular frequency, [rad/s]


Num_ord = 19;                                  % number for the highest diffraction order
div = 3;
for ana = 1:div
for heig = 1:div 
    
    d = [0.02 0.021245]; %10 is max-min   0.015 to 0.025                                        % thickness of each layer from front to back, [um]
    N = length(d);                                   % # of layers
    Period(1:N) = 1.12;                            % Period of gratings for each layer, [um]
    width = [1 1.05];                                    % width of metal strips, [um]
    psi = width/Period(1);                     % filling ratio of grating layer
    f1 = [0 0];                                              % normalized position for left-end of metal strip
    f2 = [psi];                                           % normalized position for right-end of metal strip
for ind = 1:length(lambda)
    % Incidence medium
      e(1) = 1.35149 ;            %1.30 to 1.40                                     % Usually is air or vacuum 
    % Layered structure
      e_m(1) = Palik_SiO2(lambda(ind));   % Ridge material (metal)  
      e_d(1) = Palik_SiO2(lambda(ind)); 
      
      e_m(2) = Palik_Au(lambda(ind));   % Ridge material (metal)  
      e_d(2) = Palik_SiO2(lambda(ind));                       % Groove material (air)
    %Substrate
      e(2)= Palik_SiO2(lambda(ind));                                                    %  air or opaque substrate
    %==========================================
    
    [Ref(ind), Tran(ind)] = RCWA_Multi_TM(N, e_m, e_d, f1, f2, Period, d, e, lambda(ind), theta, Num_ord); 

end
       disp(ana);

      if ana == 1 && heig ==1
        temp=Tran;
      else temp = [temp;Tran];
      end
      
      temp((div*(ana-1))+heig,1)=ana;   
      temp((div*(ana-1))+heig,2)=heig;

end

end
 
%writematrix(DAta2parameter,'data2parameter.csv') 
% legend('Tran');
% xlim([5000,10000]);
% xlabel('Wavenumber, \nu (cm^-^1)');
% ylabel('R,T,\alpha');
