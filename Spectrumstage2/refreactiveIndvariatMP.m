%2nd stage search

c0 = 299792458;                             % speed of light in space, [m/s]
theta =0;                                           % angle of incidence, [rad]
wn = linspace(5000,10000,1000) ;              % wavenumber, [1/cm]
lambda =1e4./wn;                           % wavelength, [um]  %10^7 / wn in nm (2000-1000)
w = 1e6*2*pi*c0./lambda;            % angular frequency, [rad/s]
Num_ord = 11;                                  % number for the highest diffraction order
div=40;
k=0.09;  %h
m=0.704;  %w
l=0.88;  %p

d = [0.02 k];                                           % thickness of each layer from front to back, [um]
N = length(d);                                   % # of layers
Period(1:N) = l;                          % Period of gratings for each layer, [um]
width = [1 m];                                    % width of metal strips, [um]
psi = width/Period(1);                     % filling ratio of grating layer
f1 = [0 0];                                              % normalized position for left-end of metal strip
f2 = [psi];                                           % normalized position for right-end of metal strip

for ref = 1:5

for ind = 1:length(lambda)
    % Incidence medium
      %e(1) = 1.50 + ana/(div*10) ;            %1.30 to 1.40                                     % Usually is air or vacuum 
      e(1) = 1.33 + (ref*0.01)/(5/10);
    % Layered structure
      e_m(1) = Palik_SiO2(lambda(ind));   % Ridge material (metal)  COVER
      e_d(1) = Palik_SiO2(lambda(ind)); 
      
      e_m(2) = Palik_Au(lambda(ind));   % Ridge material (metal)     GRATING
      e_d(2) = Palik_SiO2(lambda(ind));                       % Groove material (air)
    %Substrate
      e(2)= Palik_SiO2(lambda(ind));                                                    %  air or opaque substrate
    %==========================================
    
    [Ref(ind), Tran(ind)] = RCWA_Multi_TM(N, e_m, e_d, f1, f2, Period, d, e, lambda(ind), theta, Num_ord); 

end
plot(wn,Tran,LineWidth=1.1);

hold on
end
set(gca, 'fontweight','bold','fontsize',14)
hold off
