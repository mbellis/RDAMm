%=======================
% FUNCTION PRODUNIFORM %
%=======================

% PRODUNIFORM calculates the product of several independant uniform random variables.

% INPUT PARAMETERS
% 1 ProdNb: the number of variables in the product

% OUTPUT PARAMETERS
% 1 PpvList: the list of Ppv
% 2  PPvCdf: the probabilities for Ppv

% COPYRIGHT
% Translated from C program with this Copyright
%       Copyrights Orestis Georgiou   25-8-09
%
%       This program is related to the paper:
%       Product of n independent Uniform Random Variables.

%        When Theorem 1. of the above paper is integrated, the incomplete gamma
%        function, Gamma[n,x] is obtained which for integer n can be expressed
%        analytically. Thus, this program calculates the probability that a random
%        variable X, consisting of the product of n independent and identically
%        distributed uniform [a,b] random variables X_{i}, i=1,2..n, takes on a value
%        less than or equal to tau.


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


function [PpvList,PpvCdf]=produniform(ProdNb)
PpvCdf=[];
%limits of the uniforma distributions (inferior limit can't be equal to 0)
a=eps;
b=1;
%construct the list of Ppv values
%PpvList=(0.01:0.01:1)'.^ProdNb;
PpvList=(0:0.01:1)'.^ProdNb;
for PpvL=1:length(PpvList)
    Ppv=PpvList(PpvL);
    PpvProb=0;
    Test=1;
    while Ppv> (a^(ProdNb-Test))*b^Test   
        F1Val=F(a,b,Test,ProdNb,(a^(ProdNb-Test+1))* (b^(Test-1)));
        F2Val=F(a,b,Test,ProdNb,(a^(ProdNb-Test))* (b^Test));
        Diff=F2Val-F1Val;
        PpvProb=PpvProb+Diff;
        Test=Test+1;
    end
    F1Val=F(a,b,Test,ProdNb,(a^(ProdNb-Test+1))* (b^(Test-1)));
    FPpvVal=F(a,b,Test,ProdNb,Ppv);
    PpvProb=PpvProb+FPpvVal-F1Val;
    PpvCdf=[PpvCdf;PpvProb];
end
%add 0 if necessary
if PpvList(1)~=0
    PpvList=[0;PpvList];
    PpvCdf=[0;PpvCdf];
else
    PpvCdf(1)=0;
end
        

function FVal=F(a,b,Test,ProdNb,Val)
FVal=0;
for j=0:ProdNb-Test
    A=0;
    for m=1:ProdNb-1
        A= A+((log( ((a^j)*(b^(ProdNb-j)))/(Val)))^m) / factorial(m);
    end
    FVal=FVal+(Val)*factorial(ProdNb-1)*(1+A)*ProdNb*((-1)^j)/(factorial(j)*factorial(ProdNb-j)*((b-a)^ProdNb));
end



