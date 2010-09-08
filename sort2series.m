function [SVal1,PVal2,DirSortIndex1,InvSortIndex1,S2PIndex2]...
   =sort2series(Val1,Val2)
% sorting Val1 and Val2 according to the value1 Direct_sort_index
[SVal1 DirSortIndex1]=sort(Val1);
PVal2=Val2(DirSortIndex1);
% calculus of the InvSortIndex1 
[Temp InvSortIndex1]=sort(DirSortIndex1);
% calculus of the index allowing the passage between sorted and phased values for value2
[Temp P2SIndex2]=sort(PVal2);
[Temp S2PIndex2]=sort(P2SIndex2);

   
