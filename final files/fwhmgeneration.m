fwhm=[];
a=[];
fom=[];
for j = 1:size(selectedstage2)
    Tran=selectedstage2(j,:);
    [pks,locs,w] = findpeaks(Tran,'MinPeakHeight',0.25,'MinPeakProminence',0.05,'MinPeakDistance',10,'WidthReference','halfheight');
    
    if(w(2)>50) 
        fwhm = [fwhm;0];
        a = [a;j];
        fom = [fom;0];
    
else
    fwhm = [fwhm;w(2)];
    fom = [fom;sensi(j)/fwhm(j)];
end
end

fwhm=transpose(fwhm);
