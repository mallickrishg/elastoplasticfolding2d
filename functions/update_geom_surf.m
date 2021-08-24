function [allsegs,FaultPatches,LayerEndpts,surfsegs] = update_geom_surf(Ux,Uz,Uxf,Uzf,Uxs,Uzs,LayerEndpts,FaultPatches,surfsegs)


allsegs=[];

for loop=1:length(LayerEndpts)
    
    n=size(LayerEndpts{loop},1);
    temp = [LayerEndpts{loop}(1:end-1,:) + [Ux(1:n-1) Uz(1:n-1)] LayerEndpts{loop}(2:end,:) + [Ux(2:n) Uz(2:n)]];
    allsegs = [allsegs; temp];
    

    LayerEndpts{loop} = LayerEndpts{loop}+ [Ux(1:n) Uz(1:n)];
    
    Ux(1:n)=[];
    Uz(1:n)=[];
    
end

if ~isempty(FaultPatches)
    n=length(Uxf);
    FaultPatches = FaultPatches + [Uxf(1:n-1) Uzf(1:n-1) Uxf(2:n) Uzf(2:n)];
end

allsegs = [allsegs; FaultPatches];



n=length(Uxs);
surfsegs = surfsegs + [Uxs(1:n-1) Uzs(1:n-1) Uxs(2:n) Uzs(2:n)];
