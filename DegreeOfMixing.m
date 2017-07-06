function [DoM,LocalMixingDegree] = DegreeOfMixing(CL_1,CL_2,k1,k2)
%% Input Description
%%% CL_1: Nx1 cell list as Cluster list, each cell represents a cluster with at least one elements
%%% CL_2: Nx1 cell list as Cluster list, each cell represents a cluster with at least one elements
%%% k1: an integer represent first k number of clusters are interest for CL_1
%%% k2: an integer represent first k number of clusters are interest for CL_2
%% Outputs description:
% DoM: Is a 1x1 double, represents the global Degree of Mixing measured between first k1 clusters of CL_1 with
%      first k2 clusters of CL_2.
% LocalMixingDegree: Is 1xk2 double array, indicates the local degree of mixing measured for first k2 clusters in CL_2.
%% Read Me:
% This code is Published for "Cluster fusion-fission dynamics in the Singapore stock exchange", 
% by Boon Kin Teh and Siew Ann Cheong.
% Please refer to the paper for more details, and cite the paper when you
% are using this code for significance testing analysis,
% Thank you.

%% Lastest updated date:
% 03 July 2017

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Start here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine Elements appeared in both CL_1 and CL_2
NumberOfElement = 0;
for cl_i = 1:k1
    NumberOfElement = NumberOfElement+ length(CL_1{cl_i,1});
end
Element_CL1 = zeros(1,NumberOfElement);
NumberOfElement = 0;
for cl_i = 1:k1
    E = CL_1{cl_i,1};
    Element_CL1(NumberOfElement+1:NumberOfElement+length(E)) = E;
    NumberOfElement = NumberOfElement + length(E);
end
NumberOfElement = 0;
for cl_j = 1:k2
    NumberOfElement = NumberOfElement+ length(CL_2{cl_j,1});
end
Element_CL2 = zeros(1,NumberOfElement);
NumberOfElement = 0;
for cl_j = 1:k2
    E = CL_2{cl_j,1};
    Element_CL2(NumberOfElement+1:NumberOfElement+length(E)) = E;
    NumberOfElement = NumberOfElement + length(E);
end
ElementInterested = intersect(Element_CL1,Element_CL2);
clearvars -except CL_1 CL_2 k1 k2 ElementInterested

%% Determine Degree of Mixing
LocalMixingDegree = zeros(1,k2);
for cl_j = 1:k2
    %%% Determine Local Mixing Degree
    a= 1/k1; A = (a^3-a^2-a+1/2)/(a^2-a); 
    f_x = zeros(1,k1);
    
    Ecl_j = CL_2{cl_j,1};
    Ecl_j = intersect(Ecl_j,ElementInterested);
    for cl_i = 1:k1
        Ecl_i = CL_1{cl_i,1};
        Ecl_i = intersect(Ecl_i,ElementInterested);
        f_x(1,cl_i) = (1+length(intersect(Ecl_i,Ecl_j)))/(1+length(Ecl_i));
    end
    f_x = f_x/sum(f_x);
    f_x = ( f_x.*(f_x-1).*exp(-(A-f_x).^2) ) ./ ( a*(a-1)*exp(-(A-a)^2) );
    LocalMixingDegree(1,cl_j) = prod(f_x)^a;
end
%%% Determine Global Mixing Degree
DoM = mean(LocalMixingDegree);