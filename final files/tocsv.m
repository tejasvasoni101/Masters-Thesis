sensi(c3:end) = [];
heat1(c3:end) = [];
heat2(c3:end) = [];
hi2(c3:end) = [];
wi2(c3:end) = [];
pie2(c3:end) = [];
wn = linspace(5000,10000,297) ;
div =40; 
for i = 1:length(wi2)
    hi2(i)=0.02 + (hi2(i)/div)*0.08;
    pie2(i)=0.4+(pie2(i)/div)*1.2;
    wi2(i)=(wi2(i)/div)*pie2(i);
end

for i = 1:length(wi)
    hi(i)=0.02 + (hi(i)/div)*0.08;
    pie(i)=0.4+(pie(i)/div)*1.2;
    wi(i)= (wi(i)/div)*pie(i);
end

hiwipi=[];
hiwipi=[hiwipi;hi2];
hiwipi=[hiwipi;wi2];
hiwipi=[hiwipi;pie2];
better = heat1/mean(heat1);
better = better + heat2/mean(heat2);

hiwipi=[hiwipi;better];

hiwipi=transpose(hiwipi);
%writematrix(hiwipi,'hiwipi7.csv')