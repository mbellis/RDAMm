%=======================
% FUNCTION SORT2SERIES 
%=======================
%
% SORT2SERIES sorts one series of values according to the order of another series of values
%
% INPUT PARAMETERS
%  1 Val1: The first series of values which is sorted
%  2 Val2:  The second series of values sorted as the first one (sayed as 'in phase')
%
% OUTPUT PARAMETERS
%  1         SVal1: sorted first series of values
%  2         PVal2: 'in phase' second series of values
%  3 DirSortIndex1: direct sort index of first series of values (Val1 => SVal1)
%  4 InvSortInde1x: inverse sort index of the first series of values (SVal1 => Val1)
%  5     S2PIndex2: index allowing to find in phase values from sorted values 
%                   of the second series of values (SVal2 => PVal2)


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤%
%                          c) Michel Bellis                                                %
%                          michel.bellis@crbm.cnrs.fr                                      %
%            Affiliation:  CNRS (Centre National de la Recherche Scientifique - France)    %                               
%  Bioinformatic Project:  ARRAYMATIC => http://code.google.com/p/arraymatic               %
%        Code Repository:  GITHUB => http://github.com/mbellis                             %
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤%

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
%  THIS CODE IS DISTRIBUTED UNDER THE CeCILL LICENSE, WHICH IS COMPATIBLE WITH       %
%  THE GNU GENERAL PUBLIC LICENCE AND IN ACCORDANCE WITH THE EUROPEAN LEGISLATION.   %
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

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

   
