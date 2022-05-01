%2nd stage search

c0 = 299792458;                             % speed of light in space, [m/s]
theta =0;                                           % angle of incidence, [rad]
wn = linspace(5000,10000,85) ;              % wavenumber, [1/cm]
lambda =1e4./wn;                           % wavelength, [um]  %10^7 / wn in nm (2000-1000)
w = 1e6*2*pi*c0./lambda;            % angular frequency, [rad/s]
Num_ord = 7;                                  % number for the highest diffraction order
div = 10;
peak1pos=[];
peak2pos=[];
finalfiltered=[];
faultyref=[];
for j = 1:size(selected)-1
d = [0.02 (0.02 + ((hi(j))/div)*0.08)];                                           % thickness of each layer from front to back, [um]
N = length(d);                                   % # of layers
Period(1:N) = 0.4+(pie(j)/div)*1.2;                          % Period of gratings for each layer, [um]
width = [1 0.05 + wi(j)/(div)*Period(1)];                                    % width of metal strips, [um]
psi = width/Period(1);                     % filling ratio of grating layer
f1 = [0 0];                                              % normalized position for left-end of metal strip
f2 = [psi];                                           % normalized position for right-end of metal strip
lastval=1000;
count=0;
maxcontinuousshift=0;
for ref = 1:20

for ind = 1:length(lambda)
    % Incidence medium
      %e(1) = 1.50 + ana/(div*10) ;            %1.30 to 1.40                                     % Usually is air or vacuum 
      e(1) = 1.33 + (ref*0.01)/2;
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
[pks,locs] = findpeaks(Tran,'MinPeakHeight',0.25,'MinPeakProminence',0.05,'MaxPeakWidth',20,'MinPeakDistance',10)
si = size(locs);
if(si(2)<2)

    peak2pos(end+1) = -1;
    peak1pos(end+1) = -1;
    
    faultyref(end+1)=ref;
else
peak1pos(end+1) = locs(1); %flat/steps are due to accuracy of loc(being integer)
peak2pos(end+1) = locs(2);
end;


    if(peak2pos(end)==-1 ) 
        lastval=1000;
        maxcontinuousshift=max(count,maxcontinuousshift);
        count=1;
    elseif(peak2pos(end)<lastval+1 ) 
        count=count+1;
        maxcontinuousshift=max(count,maxcontinuousshift);
        lastval=peak2pos(end);
    else
        lastval=peak2pos(end);
        maxcontinuousshift=max(count,maxcontinuousshift);
        count=1;
    end
end
if(maxcontinuousshift>10) % 2 DENOTES 3 WORKING RANGE LIKE 1.33 TO 1.35
    finalfiltered(end+1)=j;
end
%j from 1 to selected
end