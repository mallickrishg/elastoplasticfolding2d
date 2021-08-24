function [pts,map_pts_patches,is_end,LayerPatches,LayerProperties,patchnums] = get_points_simplify(LayerPatches,LayerProperties,FaultPatches,is_end,patchnums)


        
index = zeros(size(LayerPatches,1),1);            
if ~isempty(FaultPatches)
   
        %remove layer segments too close to fault  
        dist=0.1;  %shortest distance between layer point and fault


        for loop=1:size(is_end,1)

            
                if is_end(loop)
                    
           
                    ind = (abs(LayerPatches(loop,1) - FaultPatches(:,1))<dist &...
                        abs(LayerPatches(loop,2) - FaultPatches(:,2))<dist)  |...
                        (abs(LayerPatches(loop,3) - FaultPatches(:,1))<dist  &...
                        abs(LayerPatches(loop,4) - FaultPatches(:,2))<dist);
                    
                    if sum(ind)>0; 
                        index(loop)=1;
                    end
             
                
                end
           
        end


end
    
    
index = logical(index);    
LayerPatches(index,:) = [];
patchnums(index) = [];
LayerProperties(index) = [];    
is_end(index) = [];
        

%remove patches entirely above free surface
index = LayerPatches(:,2)>0 & LayerPatches(:,4)>0;
LayerPatches(index,:)=[];
patchnums(index)=[];
LayerProperties(index)=[];
is_end(index) = [];



%shorten patches that extend above free surface
index = LayerPatches(:,2)>0;
out = lineSegmentIntersect(LayerPatches(index,:),[-10^3 0 10^3 0]);
LayerPatches(index,1:2)=[out.intMatrixX(out.intAdjacencyMatrix) out.intMatrixY(out.intAdjacencyMatrix)];

index = LayerPatches(:,4)>0;
out = lineSegmentIntersect(LayerPatches(index,:),[-10^3 0 10^3 0]);
LayerPatches(index,3:4)=[out.intMatrixX(out.intAdjacencyMatrix) out.intMatrixY(out.intAdjacencyMatrix)];

        

%remove very small patches near surface
index = sqrt( (LayerPatches(:,1)-LayerPatches(:,3)).^2 + (LayerPatches(:,2)-LayerPatches(:,4)).^2 )<dist/10;%0.2;    
LayerPatches(index,:)=[];
patchnums(index)=[];
LayerProperties(index)=[];
is_end(index) = [];

  

%remove large patches near surface that result from connecting
%eroded segments
index = sqrt( (LayerPatches(:,1)-LayerPatches(:,3)).^2 + (LayerPatches(:,2)-LayerPatches(:,4)).^2 )>5*mean(sqrt( (LayerPatches(:,1)-LayerPatches(:,3)).^2 + (LayerPatches(:,2)-LayerPatches(:,4)).^2 ));    
LayerPatches(index,:)=[];
patchnums(index)=[];
LayerProperties(index)=[];
is_end(index) = [];


%remove very shallow segments
index = (LayerPatches(:,2)+LayerPatches(:,4))/2 > -0.01;
LayerPatches(index,:)=[];
patchnums(index)=[];
LayerProperties(index)=[];
is_end(index) = [];

        



%make list of unique patch end points for computing displacements later on
[pts,i,map_pts_patches] = unique([LayerPatches(:,1:2);LayerPatches(:,3:4)],'rows');
map_pts_patches = [map_pts_patches(1:end/2) map_pts_patches(end/2+1:end)];

end


