
function pm = make_pm(allsegs)


centers = [(allsegs(:,1)+allsegs(:,3))/2 (allsegs(:,2)+allsegs(:,4))/2];
npatch = size(allsegs,1);

W = sqrt( (allsegs(:,1)-allsegs(:,3)).^2 + (allsegs(:,2)-allsegs(:,4)).^2 );

[D,index1] = min([allsegs(:,2) allsegs(:,4)],[],2); 
D=-D;

angle = 180/pi*atan2(allsegs(:,2)-allsegs(:,4),allsegs(:,1)-allsegs(:,3));
dip = angle;

index = angle<=0;  dip(index) = -angle(index);
index = angle>0;  dip(index) = 180-angle(index);

segsE = [allsegs(:,1) allsegs(:,3)];
E = zeros(npatch,1); 

i = dip<=90; 
E(i)=max(segsE(i,:),[],2);  
E(~i)=min(segsE(~i,:),[],2);

pm =[10^4*ones(npatch,1) W D dip zeros(npatch,1) E zeros(npatch,1)];

end
