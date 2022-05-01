%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is the main program of RCWA method to calculate the spectral directional radiative properties for TE and TM incidence 
% Dr. Zhuomin Zhang's group at Georgia Tech
% Last modified by Bo Zhao (September 2014)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;
c0 = 299792458;                             % speed of light in space, [m/s]
theta =0;                                           % angle of incidence, [rad]
wn = linspace(5000,10000,400) ; % wavenumber, [1/cm]
lambda =1e4./wn;                           % wavelength, [um]
w = 1e6*2*pi*c0./lambda;            % angular frequency, [rad/s]
div = 30;             %division for which loop runs


wi=zeros(1,div*div*div);
hi=zeros(1,div*div*div);
pie=zeros(1,div*div*div);
    righ=0;
    lef=1000;
numb = zeros(1,div*div*div);
Num_ord = 11;                                  % number for the highest diffraction order

c=1;
for peri = 1:div
for widi = 1:div      %for loop for width 
for heig = 1:div      %fo height
      %fo height    
    d = [0.02 (0.02 + (heig/div)*0.08)];%from 0.03 to 0.03+0.02  %10 is max-min   0.015 to 0.025                                        % thickness of each layer from front to back, [um]
    N = length(d);                                   % # of layers
    Period(1:N) = 0.4+(peri/div)*1.2; %Period(1:N) = 1.12;        % Period of gratings for each layer, [um]
    width = [1 widi/(div)*Period(1)];  %from 1.05 to 1.05+1
    % width of metal strips, [um]
    psi = width/Period(1);                     % filling ratio of grating layer
    f1 = [0 0];                                              % normalized position for left-end of metal strip
    f2 = [psi];                                           % normalized position for right-end of metal strip
for ind = 1:length(lambda)
    % Incidence medium
      %e(1) = 1.50 + ana/(div*10) ;            %1.30 to 1.40                                     % Usually is air or vacuum 
      e(1) = 1.35;
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

%       disp(widi);

      if widi == 1 && heig ==1 && peri==1
        temp=Tran;      %to initialize temp variable
      else temp = [temp;Tran];   %appened tran spectrum row to temp dataset
      end
      
      temp((div*div*(peri-1))+(div*(widi-1))+heig,1)=peri;   
      temp((div*div*(peri-1))+(div*(widi-1))+heig,2)=widi;   %numbering the spectrum in temp
      temp((div*div*(peri-1))+(div*(widi-1))+heig,3)=heig;

      
      pk = findpeaks(Tran,'MinPeakHeight',0.3,'MinPeakProminence',0.08,'MinPeakDistance',10,'MaxPeakWidth',190,'Threshold',0.05,'WidthReference','halfheight')
      sz = size(pk)          % no. of peaks in pk
      if sz(2)==2             % atleast 1 peak found in a given plot for a design
       wi(c) = widi;     % 1st parameter of selected design
       hi(c) = heig;     % 2nd parameter of selected design
       pie(c) = peri;    % 3rd parameter of selected design
       numb(c) =(div*div*(peri-1))+(div*(widi-1))+heig ;
       phi=[hi;wi];
       phi=[phi;pie];
       phi=[phi;numb];       %numb coresp to the row in temp which is selected
       c=c+1;
      end
    
    
end
end
end

selected=temp(1,4:end);   %matrix of selected profiles (avoid the first)
parfor i = 1:c-1
    selec = numb(1,i);
    selec2 = temp(selec,4:end);
    selected = [selected;selec2];
end

%2nd stage search
selected(1,:) = [];

wi2=zeros(1,div*div*div); 
hi2=zeros(1,div*div*div); 
pie2=zeros(1,div*div*div);
heat1=zeros(1,div*div*div); 
heat2=zeros(1,div*div*div); 
sensi=zeros(1,div*div*div);
selectedstage2=[];
peak1pos=zeros(1,div*div*div);
peak2pos=zeros(1,div*div*div);
finalfiltered=zeros(1,div*div*div);
faultyref=zeros(1,div*div*div);
c2=1;
c3=1;
for j = 1:size(selected)-1
d = [0.02 (0.02 + ((hi(j))/div)*0.08)];                                           % thickness of each layer from front to back, [um]
N = length(d);                                   % # of layers
Period(1:N) = 0.4+(pie(j)/div)*1.2;                          % Period of gratings for each layer, [um]
width = [1 wi(j)/(div)*Period(1)];                                    % width of metal strips, [um]
psi = width/Period(1);                     % filling ratio of grating layer
f1 = [0 0];                                              % normalized position for left-end of metal strip
f2 = [psi];                                           % normalized position for right-end of metal strip
lastval=1000;
lastvalrefpeak = 0;
count=0;
maxcontinuousshift=0;
maxheat1=0;
maxheat2=0;
for ref = 1:20 %change divisor at line 115 accordingly to ratio to 10

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
[pks,locs] = findpeaks(Tran,'MinPeakHeight',0.25,'MinPeakProminence',0.05,'MinPeakDistance',10)
si = size(locs);
if(si(2)<2)

    peak2pos(c2) = -1;
    peak1pos(c2) = -1;
    c2=c2+1;
    %faultyref(end+1)=ref;
elseif(si(2)==2)
peak1pos(c2) = locs(1); %flat/steps are due to accuracy of loc(being integer)
peak2pos(c2) = locs(2);
maxheat1=maxheat1+pks(1);
maxheat2=maxheat2+pks(2);
c2=c2+1;

end

    if(peak2pos(c2-1)==-1 ) 
        lastval=1000;
        maxcontinuousshift=max(count,maxcontinuousshift);
        count=1;
    elseif(peak2pos(c2-1)<lastval+1 ) %check sign
        count=count+1;
        maxcontinuousshift=max(count,maxcontinuousshift);
        lastval=peak2pos(c2-1);
        
        lef = min(lef,locs(2)); 
        righ = max(righ,locs(2));
    else
        lastval=peak2pos(c2-1);
        maxcontinuousshift=max(count,maxcontinuousshift);
        count=1;
        
        lef = min(lef,locs(2)); 
        righ = max(righ,locs(2));
    end
end
if(maxcontinuousshift>10)     %2 DENOTES 3 WORKING RANGE LIKE 1.33 To 1.35
    heat1(c3)=maxheat1;
    heat2(c3)=maxheat2;
    sensi(c3)=righ-lef;
    hi2(c3) =hi(j) ;
    wi2(c3) =wi(j) ;
    pie2(c3)=pie(j);
    tempp = selected(j,:);
    selectedstage2 = [selectedstage2;tempp];
    c3=c3+1;
end
    righ=0;
    lef=1000;
end


% writematrix(DAta2parameter,'data2parameter.csv') 
% legend('Tran');
% xlim([5000,10000]);
% xlabel('Wavenumber, \nu (cm^-^1)');
% ylabel('R,T,\alpha');