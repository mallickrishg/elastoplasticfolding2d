function [rcv,rcvdeep] = makefault_pm(allsegs,G,nu)
addpath ~/Dropbox/scripts/unicycle/matlab/
import unicycle.*

patchfname = 'dummy.seg';
patchdeepname = 'dummydeep.seg';

n = length(allsegs(:,1));
L = 10000e3;

w = 1e3.*sqrt( (allsegs(:,1)-allsegs(:,3)).^2 + (allsegs(:,2)-allsegs(:,4)).^2 );

% index for left pt is up
index = abs(allsegs(:,2))==min(abs(allsegs(:,[2,4])),[],2);

dip = atan2d(abs(allsegs(:,4))-abs(allsegs(:,2)),allsegs(:,3)-allsegs(:,1));
dipindex = dip<0;

dip(dipindex) = 180 - abs(dip(dipindex));

x0 = 1e3.*allsegs(:,1);
x0(dipindex) = 1e3.*allsegs(dipindex,3);
z0 = 1e3.*allsegs(:,2);
z0(dipindex) = 1e3.*allsegs(dipindex,4);


% create fault and layers
fileID = fopen(patchfname,'w');

for i = 1:n
    
    fprintf(fileID,'%d %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %d %.9f %.9f %d %d\n',...
        1, 1,-L/2, x0(i), -z0(i),L,w(i), 0,   dip(i), 90, L, w(i), 1.0, 1.0);
    
end
fclose(fileID);

if dipindex(end)==1
    xdeep = 1e3.*allsegs(end,1);
    zdeep = 1e3.*allsegs(end,2);
else
    xdeep = 1e3.*allsegs(end,3);
    zdeep = 1e3.*allsegs(end,4);
end

wdeep = 10000e3;
% deep driving fault
fileID = fopen(patchdeepname,'w');
    fprintf(fileID,'%s\n',...
        '# n  Vpl      x1              x2   x3   Length       Width   Strike  Dip  Rake      L0       W0        qL  qW');
fprintf(fileID,'%d %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %d %d\n',...
                1, 1,-L/2, xdeep,    -zdeep,L,wdeep, 0,   0, 90, L, wdeep, 1.0, 1.0);
fclose(fileID);



earthModel = greens.okada92(G,nu);
rcv = geometry.receiver(patchfname,earthModel);
rcvdeep = geometry.receiver(patchdeepname,earthModel);


end