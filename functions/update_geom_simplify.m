function [LayerPatches,FaultPatches] = update_geom_simplify(Ux,Uz,Uxf,Uzf,LayerPatches,FaultPatches,map_pts_patches)



LayerPatches(:,1) = LayerPatches(:,1) + Ux(map_pts_patches(:,1));
LayerPatches(:,2) = LayerPatches(:,2) + Uz(map_pts_patches(:,1));

LayerPatches(:,3) = LayerPatches(:,3) + Ux(map_pts_patches(:,2));
LayerPatches(:,4) = LayerPatches(:,4) + Uz(map_pts_patches(:,2));
    
   
if ~isempty(FaultPatches)
    n=length(Uxf);
    FaultPatches = FaultPatches + [Uxf(1:n-1) Uzf(1:n-1) Uxf(2:n) Uzf(2:n)];
end
