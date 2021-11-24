function [enrdomain]=crackDetect(Cr,istep,tipElem,splitElem,vertexElem,cornerElem,enrdomain)

%Purpose
%identify elements that are in close proximity to the crack tip
if istep == 1       %first step in growth simulation
    for icrack=1:size(Cr,2)
        Tip1 = Cr(icrack).coor(1,:) ;
        Tip2 = Cr(icrack).coor(size(Cr(icrack).coor,1),:) ;
        [enrdomain1,radius] = ENRdomainf(Tip1,Tip2) ;
        [enrdomain2,radius] = ENRdomainf(Tip2,Tip1) ;
        enrdomain = [enrdomain; enrdomain1; enrdomain2] ;
    end
else
    for icrack = 1:size(Cr,2)
        Tip1 = Cr(icrack).coor(1,:) ;
        Tip2 = Cr(icrack).coor(2,:) ;
        [enrdomain1,radius] = ENRdomainf(Tip1,Tip2) ;
        
        Tip1 = Cr(icrack).coor(size(Cr(icrack).coor,1),:) ;
        Tip2 = Cr(icrack).coor(size(Cr(icrack).coor,1)-1,:) ;
        [enrdomain2,radius] = ENRdomainf(Tip1,Tip2) ;
        enrdomain = [enrdomain;tipElem;splitElem;vertexElem;cornerElem;enrdomain1;enrdomain2] ;
    end
end

enrdomain = unique(enrdomain) ;
