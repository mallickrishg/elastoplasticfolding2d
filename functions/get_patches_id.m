function [pts,LayerPatches,LayerProperties,is_end,patchnums,map_pts_patches,id] = get_patches_id(LayerEndpts,LayerProp)

        
    
    
    
%need to move points slightly off the fault to avoid numerical problems
LayerPatches=[];
LayerProperties=[];
is_end = [];
id = [];

for loop=1:length(LayerEndpts)
    
    
    
    LP = [LayerEndpts{loop}(1:end-1,1) LayerEndpts{loop}(1:end-1,2) LayerEndpts{loop}(2:end,1) LayerEndpts{loop}(2:end,2)];
    
    %identify patches at ends of Layer contacts
    ends = zeros(size(LP,1),1);
    ends(1)=1;
    ends(end)=1;
       
    LayerPatches = [LayerPatches; LP]; 
    is_end = [is_end; ends];
    
    id = [id;LP(:,1).*0 + loop];
    
    if ~isempty(LayerProp)
        LayerProperties = [LayerProperties; LayerProp(loop)*ones(size(LP,1),1)];
    end
    
end


%make list of unique patch end points for computing displacements later on
[pts,i,map_pts_patches] = unique([LayerPatches(:,1:2);LayerPatches(:,3:4)],'rows');
map_pts_patches = [map_pts_patches(1:end/2) map_pts_patches(end/2+1:end)];

patchnums = (1:size(LayerPatches,1))';
end