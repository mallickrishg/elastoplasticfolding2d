% DRIVER script to run fold growth model for different parameters
% Rishav Mallick, EOS, 2020

clear
addpath ~/Dropbox/scripts/unicycle/matlab/
import unicycle.*
addpath functions/

% set initial deviatoric stress state
mdparams.Sxxi = -0.*1e6;
mdparams.Sxzi = 0;
mdparams.Szzi = 0;
mdparams.lambda = 1;% ratio of Sxx to Szz (horizontal loading)

% friction coefficients
mdparams.mu_F = 0.0;% set main fault friction coefficient to 0
mdparams.mu_L = 0.01;% friction coefficient of bedding plane contacts

% Elastic properties
mdparams.G = 10e9;% shear modulus
mdparams.nu = 0.25;% poisson's ratio

mdparams.C = 1e6;% Cohesive strength (in Pa)

mdparams.totalslip = 2;% in km (cumulative shortening)
mdparams.nslips = 40;% number of slip increments


%% run simulations
nlayers = 5;
nruns = length(nlayers);
for i = 1:nruns
    mdparams.geom = ['faults/fault_smooth_100_' num2str(nlayers(i)) '_layers'];
    mdparams.layerprop = [1:2:(2*nlayers(i)-1)];
    EXPLORE_foldmodel_imposedslip(mdparams);    
end

% Gvec = [1,3,5,10,30,1,3,5,10,30,1,3,5,10,30].*1e9;
% mulvec = [0.01.*ones(5,1);0.1.*ones(5,1);0.4.*ones(5,1)];
% Gvec = [0.01,0.05,0.1,0.5,...
%     0.01,0.05,0.1,0.5...
%     0.01,0.05,0.1,0.5].*1e9;
% mulvec = [0.01.*ones(4,1);...
%     0.1.*ones(4,1);...
%     0.6.*ones(4,1)];

% nruns = 1;
% nlayers = 30;
% for i = 1:nruns
%     %mdparams.G = Gvec(i);
%     %mdparams.mu_L = mulvec(i);
%     mdparams.geom = ['fault_smooth_500d_' num2str(nlayers) '_layers'];
%     mdparams.layerprop = [1:2:(2*nlayers-1)];
%     EXPLORE_foldmodel_imposedslip(mdparams);
% end



