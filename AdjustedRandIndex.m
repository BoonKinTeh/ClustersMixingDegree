function ARI = AdjustedRandIndex(CL_1,CL_2)
%% Input Description
%%% CL_1: Nx1 cell list as Cluster list, each cell represents a cluster with at least one elements
%%% CL_2: Nx1 cell list as Cluster list, each cell represents a cluster with at least one elements
%% Outputs description:
% ARI: Is a 1x1 double, represents the adjusted Rand Index measured between clusters CL_1 and CL_2,
%      the global mixing degree hence is 1-ARI
%
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
%% Contruct Contingency Table
NumberOfCluster_CL1 = size(CL_1,1);
NumberOfCluster_CL2 = size(CL_2,1);

ContingencyTable = zeros(NumberOfCluster_CL1,NumberOfCluster_CL2);
for cl_i = 1:NumberOfCluster_CL1
    Ecl_i = CL_1{cl_i,1};
    for cl_j = 1:NumberOfCluster_CL2
        Ecl_j = CL_2{cl_j,1};
        ContingencyTable(cl_i,cl_j) = length(intersect(Ecl_i,Ecl_j));
    end
end
%% Measure Index
Index = 0;
for cl_i = 1:NumberOfCluster_CL1
    for cl_j = 1:NumberOfCluster_CL2
        Index = Index+nchoosek(ContingencyTable(cl_i,cl_j),2);
    end
end
%% Measure Max Index
MaxIndex = 0;
for cl_i = 1:NumberOfCluster_CL1
    MaxIndex = MaxIndex + nchoosek(sum(ContingencyTable(cl_i,:)),2);
end
for cl_j = 1:NumberOfCluster_CL2
    MaxIndex = MaxIndex + nchoosek(sum(ContingencyTable(:,cl_j)),2);
end
MaxIndex = MaxIndex/2;
%% Measure Expected Index
ExpectedIndex = 0;
for cl_i = 1:NumberOfCluster_CL1
    Value = 0;
    for cl_j = 1:NumberOfCluster_CL2
        Value = Value+nchoosek(sum(ContingencyTable(:,cl_j)),2);
    end
    ExpectedIndex = ExpectedIndex + Value*nchoosek(sum(ContingencyTable(cl_i,:)),2);
end
NumberOfElements = sum(sum(ContingencyTable));
ExpectedIndex = ExpectedIndex/nchoosek(NumberOfElements,2);
%% Return ARI
ARI = (Index-ExpectedIndex)/(MaxIndex-ExpectedIndex);
