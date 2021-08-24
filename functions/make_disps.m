function [Ux,Uz] = make_disps(pm,pts,slip,nu)
addpath ~/Dropbox/scripts/unicycle/matlab/
import unicycle.*

pm = [pm zeros(size(pm,1),1) slip zeros(size(pm,1),1)];


%round-off errors can place point slightly above free surface
index = pts(:,2)>0;
pts(index,2) = -pts(index,2);

%xloc = [pts(:,1)'; 0*pts(:,1)'; pts(:,2)'];

% [U,D,Sds,flag] = disloc3d(pm',xloc,1,nu);

angle = deg2rad(pm(:,4));
L = pm(:,1);
W = pm(:,2);

xdis = pm(:,6) - W.*cos(angle);
depth = pm(:,3) - W.*sin(angle); % needs to be positive value
depth(depth<0) = 0;

type = 'd';
strike = 0;
%% compute displacements
GUx = zeros(length(pts(:,1)),length(pm(:,1)));
GUz = zeros(length(pts(:,1)),length(pm(:,1)));
x = pts(:,1);
z = pts(:,2);
parfor i = 1:length(pm(:,1))
    %disp(['Completed ' num2str(round(100*i/length(pm(:,1)))) '%'])
    [Ui,~]=unicycle.greens.computeOkada92(1,x-xdis(i),L(i)/2,z,1,nu,...
        depth(i),angle(i),L(i),W(i),type,0,strike);
    GUx(:,i) = squeeze(Ui(:,:,:,1));
    GUz(:,i) = squeeze(Ui(:,:,:,3));
end
Ux = GUx*slip;
Uz = GUz*slip;
end

