function EXPLORE_foldmodel_imposedslip(mdparams)

% addpath ~/Dropbox/scripts/unicycle/matlab/
% import unicycle.*
% addpath ~/Dropbox/scripts/utils/


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%size and number of slip increments
total_slip = mdparams.totalslip; %km
nslips = mdparams.nslips;

dslips = [0;(total_slip/(nslips-1)).*ones(nslips-1,1)];
%name of geometry file %must be in quotes
geom = mdparams.geom;
ofolder = 'results_imposedslip/';


%initial stress (Pa), negative is compression
Sxx_init = mdparams.Sxxi;  %horizontal
Sxz_init = mdparams.Sxzi;
Szz_init = mdparams.Szzi;  %vertical 
tectoniclambda = mdparams.lambda;

%coefficients of friction
%layer contact friction -- one entry for each layer property
%use 'inf' to specify bonded contacts
mu_L(1)= mdparams.mu_L; 
mu_L(2)= 10;  

%elastic shear modulus
mu = mdparams.G;
% cohesion
C = mdparams.C;
%fault friction
mu_F= mdparams.mu_F;  

%Poisson's ratio
nu = mdparams.nu;

% allow fault to be advected (1) or stationary (0)
advection = 0;
%% end of input section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic
dirfil = dir([ofolder 'modelrun*']);
ifil = length(dirfil)+1;
result_folder = [ofolder 'modelrun_' num2str(ifil) '/'];
mkdir(result_folder)
result_file = [result_folder 'inputs.dat']
T = table([mu,nu,mu_L(1),mu_F,C,Sxx_init,Sxz_init,Szz_init,tectoniclambda,total_slip,nslips]);
writetable(T,result_file,'WriteVariableNames',0,'Delimiter','\t')

%% Run solver
%number of time steps
numsteps = length(dslips);

%NOTE:  First run draw_fault.m to set up fault and layer geometry
load(geom)
if isfield(mdparams,'layerprop')
    LayerProp(mdparams.layerprop) = 2;
end


%apply Blaise filter to smooth fault (avoid kinking)
y1=bfilt(FaultPatches(:,2)')';y1 = FaultPatches(:,2);
y2=bfilt(FaultPatches(:,4)')';y2 = FaultPatches(:,4);
FaultPatches = [FaultPatches(:,1) y1 FaultPatches(:,3) [y1(2:end); y2(end)]];
%make Layer patches
[pts,LayerPatches,LayerProperties,is_end,patchnums,map_pts_patches,layerid] = get_patches_id(LayerEndpts,LayerProp);
allsegs=[LayerPatches;FaultPatches];
layerid = [layerid;zeros(length(FaultPatches),1)];

num_patches = size(allsegs,1);      

%initial stress  %add arbitrary confinement
Sxx = Sxx_init*ones(size(allsegs,1),1);
Sxz = Sxz_init*ones(size(allsegs,1),1);
Szz = Szz_init*ones(size(allsegs,1),1);
slipvec = zeros(num_patches,1);


for steps = 1:numsteps
    tic
    dslip = dslips(steps);
    
    % setup geometry
    pm = make_pm(allsegs);
    
    faultindex = [false(length(LayerPatches),1);true(length(FaultPatches),1)];
    
    [rcv,rcvdeep] = makefault_pm(allsegs,mu,nu);
    rcv.mu0 = [mu_L(LayerProperties)'; mu_F*ones(size(FaultPatches,1),1)];
    disp('Computing stress kernels')

    [Kdeep,~] = computetractionkernels(rcvdeep,rcv);
    [Kself,Ksigma] = computetractionkernels(rcv,rcv);
    Klayers = Kself(~faultindex,~faultindex);
    
    disp('Finished Computing')
    
    %confinement
    rho = 2.6e3;  %kg/m^3
    g = 9.8; % m/s^2
    
    Szz = Szz + rho*g*(allsegs(:,2)+allsegs(:,4))/2*1000;
    Sxx = Sxx + rho*g*(allsegs(:,2)+allsegs(:,4))/2*1000*tectoniclambda;
    
    % resolve background stress onto fault planes as initial tractions
    taudinit = (Szz-Sxx).*sind(rcv.dip).*cosd(rcv.dip) + Sxz.*(cosd(rcv.dip).^2 - sind(rcv.dip).^2);
    tauninit = Sxx.*sind(rcv.dip).^2 + Szz.*cosd(rcv.dip).^2 - 2.*Sxz.*sind(rcv.dip).*cosd(rcv.dip);
    
    % increment in tractions due to slip on deep detachment    
    dtaud = Kdeep*(1e3.*dslip) + Kself(:,faultindex)*(1e3*dslip.*(ones(length(FaultPatches),1))); %dslip is in km
    
    tic
    disp('Running Solver')
    
    % elastic stress from previous time step
    if steps == 1
        taud_0 = taudinit;
        taun = abs(tauninit) + C;
    else
        taun = abs(rho*g*(allsegs(:,2)+allsegs(:,4))/2*1000) + C;
    end
    
    tauyield = rcv.mu0.*taun;
    
    scf = [];
    scf.cf = 1e6;
    scf.Ksigma = Ksigma;
    scf.option = 1;%1- only shear stress, 2- shear stress balance with equality constraint on footwall

    % calculate incremental slip/strain to balance incremental loading
    rcvmod = [];
    rcvmod.N = length(taud_0(~faultindex));
    rcvmod.mu0 = rcv.mu0(~faultindex);
    
    sliplayers = stressbalance_solver_imposedslip(taud_0(~faultindex),dtaud(~faultindex),taun(~faultindex),rcvmod,Klayers,scf);    
    
    slip = (dslip*1e3).*ones(rcv.N,1); % slip should be in (m)
    slip(~faultindex) = sliplayers;
    
    disp('Finished Running Solver')
    
    slipvec = slipvec + slip;   
    
    slip_noflt = zeros(rcv.N,1);
    slip_noflt(~faultindex) = sliplayers;
    taud_0 = taud_0 + dtaud + Kself*(slip_noflt); % calculate stress balance (dtaud already had contribution from fault)
    
    %%%%%%%%%%%%%      compute displacements of layers   %%%%%%%%%%%%%
    disp('Computing Displacements')
    slipindex = slip~=0;% only consider patches that have nonzero slip
    [Ux,Uz] = make_disps(pm(slipindex,:),pts,slip(slipindex)./1e3,nu);   
    
    pm_drive2_lower = [10^4 10^3 -lower_detachxy(2) 180 0 lower_detachxy(1) 0]';
    [Ux1,Uz1] = make_disps(pm_drive2_lower',pts,dslip,nu);
    Ux = Ux + Ux1;
    Uz = Uz + Uz1;
                    
    %sometimes get nan
    index = isnan(Ux); Ux(index)=0;
    index = isnan(Uz); Uz(index)=0;
    
    disp('Computing Fault Displacements')
    
    
    [Ux1,Uz1] = make_disps(pm_drive2_lower',[FaultPatches(:,1:2);FaultPatches(end,3:4)],dslip,nu);
    
    if advection==1
        [Uxf,Uzf] = make_disps(pm(slipindex,:),[FaultPatches(:,1:2);FaultPatches(end,3:4)],slip(slipindex)./1e3,nu);
        Uxf = (Uxf + Ux1);
        Uzf = (Uzf + Uz1);
    else
        % no advection of fault (stationary)
        Uxf = (0.*Ux1);
        Uzf = (0.*Uz1);
    end
    
    %%%%%%%%%%%%%    update geometry
    disp('Updating Fault Geometry')
    [LayerPatches_new,FaultPatches_new] = update_geom_simplify(Ux,Uz,Uxf,Uzf,LayerPatches,FaultPatches,map_pts_patches);
    allsegs_new = [LayerPatches_new;FaultPatches_new];
    
    
    FaultPatches_new(end,4) =FaultPatches_new(end,2);
    allsegs_new(end,4) =FaultPatches_new(end,2);
    lower_detachxy = FaultPatches_new(end,3:4);    
    
    FaultPatches = FaultPatches_new;
    LayerPatches = LayerPatches_new;
    
    %adjust points near faults and above free surface   
    [pts,map_pts_patches,is_end,LayerPatches_new,LayerProperties,patchnums_new] = get_points_simplify(LayerPatches,LayerProperties,FaultPatches,is_end,patchnums);
    
    %identify layer patches that were discarded above
    patch_indices = ismember(patchnums,patchnums_new);
    
    %remove stresses for patches that were removed in step above
    Sxx = Sxx([patch_indices; true(size(FaultPatches,1),1)]);
    Sxz = Sxz([patch_indices; true(size(FaultPatches,1),1)]);
    Szz = Szz([patch_indices; true(size(FaultPatches,1),1)]);
    slipvec = slipvec([patch_indices; true(size(FaultPatches,1),1)]);
    taud_0 = taud_0([patch_indices; true(size(FaultPatches,1),1)]);
    layerid = layerid([patch_indices; true(size(FaultPatches,1),1)]);
    patchnums = patchnums_new;
    LayerPatches = LayerPatches_new;
    allsegs=[LayerPatches_new;FaultPatches];
    tauyield = tauyield([patch_indices; true(size(FaultPatches,1),1)]);
    incslip = slip([patch_indices; true(size(FaultPatches,1),1)]);
    % export data to textfiles
    T = table([allsegs,slipvec,sum(dslips(1:steps)).*ones(length(allsegs(:,1)),1),taud_0,tauyield,incslip,layerid]); 
    
    
    result_file = [result_folder 'dataout_' num2str(steps) '.dat']
    writetable(T,result_file,'WriteVariableNames',0,'Delimiter','\t')
    toc
    disp(['End of Step ' num2str(steps-1)])
    
 
    
end %steps

disp('End of Simulation')


end