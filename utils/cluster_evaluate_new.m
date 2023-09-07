function [purity, NMI, RI, Fscore, ARI] = cluster_evaluate_new(clusterID, groundtruthID,varargin)
% Updated by Nhan Dam on April 6, 2018

% evaluate RI and Fscore
TP = 0;
TN = 0;
FP = 0;
FN = 0;

if(numel(clusterID) ~= numel(groundtruthID))
    purity = 0;
    NMI = 0;
    ARI = 0;
    Fscore = 0;
    disp('Warning:: The cluster and groundtruth labels have different length');
    return;
end

for i = 1:length(clusterID)
    for j = i+1:length(clusterID)
        if(groundtruthID(i) == groundtruthID(j))
            if(clusterID(i) == clusterID(j))
                TP = TP + 1;
            else
                FN = FN + 1;
            end
        else
            if(clusterID(i) == clusterID(j))
                FP = FP + 1;
            else
                TN = TN + 1;
            end            
        end
    end
end

precision = TP / (TP+FP);
recall = TP / (TP+FN);

Fscore = 2 * precision * recall / (precision + recall);

RI = (TP + TN)/(TP + FP + TN + FN);
[ARI,~,~,~]=AdjustedRandIndex(groundtruthID,clusterID);

% evaluate purity and NMI

clusterList = unique(clusterID);
clusterNum = length(clusterList);

groundtruthList = unique(groundtruthID);
groundtruthNum = length(groundtruthList);

clusterIdx = cell(1, clusterNum);
clusterLen = zeros(1, clusterNum);

groundtruthIdx = cell(1, groundtruthNum);
groundtruthLen = zeros(1,groundtruthNum);

for i = 1:clusterNum
    clusterIdx{i} = find(clusterID == clusterList(i));
    clusterLen(i) = length(clusterIdx{i});
end

for i = 1:groundtruthNum
    groundtruthIdx{i} = find(groundtruthID == groundtruthList(i));
    groundtruthLen(i) = length(groundtruthIdx{i});
end

intersection = zeros(clusterNum, groundtruthNum);

for i = 1:clusterNum
    for j = 1:groundtruthNum
        intersection(i,j) = length(intersect(clusterIdx{i}, groundtruthIdx{j}));
    end
end


purity = sum(max(intersection, [], 2)) / length(clusterID);

MI = 0;
person_num = length(clusterID);
for i = 1:clusterNum
    for j = 1:groundtruthNum
        if(intersection(i,j) ~= 0)
            MI = MI + intersection(i,j) / person_num * log(person_num * intersection(i,j) / (clusterLen(i) * groundtruthLen(j)));
        end
    end
end

entropy_cluster = 0;
entropy_groundtruth = 0;

for i = 1:clusterNum
    entropy_cluster = entropy_cluster - clusterLen(i) / person_num * log(clusterLen(i) / person_num);
end

for i = 1:groundtruthNum
    entropy_groundtruth = entropy_groundtruth - groundtruthLen(i) / person_num * log(groundtruthLen(i) / person_num);
end

NMI = MI / ((entropy_cluster + entropy_groundtruth)/2);

if isempty(varargin)
    fprintf('Purity=%.4f NMI=%.4f RI=%.4f Fscore=%.4f ARI=%.4f\n',purity,NMI,RI,Fscore,ARI);
end

end


function [AR,RI,MI,HI]=AdjustedRandIndex(c1,c2)
%RANDINDEX - calculates Rand Indices to compare two partitions
% ARI=RANDINDEX(c1,c2), where c1,c2 are vectors listing the 
% class membership, returns the "Hubert & Arabie adjusted Rand index".
% [AR,RI,MI,HI]=RANDINDEX(c1,c2) returns the adjusted Rand index, 
% the unadjusted Rand index, "Mirkin's" index and "Hubert's" index.
%
% See L. Hubert and P. Arabie (1985) "Comparing Partitions" Journal of 
% Classification 2:193-218

%(C) David Corney (2000)   		D.Corney@cs.ucl.ac.uk

if nargin < 2 || min(size(c1)) > 1 || min(size(c2)) > 1
   error('RandIndex: Requires two vector arguments')
end

C=Contingency(c1,c2);	%form contingency matrix

n=sum(sum(C));
nis=sum(sum(C,2).^2);		%sum of squares of sums of rows
njs=sum(sum(C,1).^2);		%sum of squares of sums of columns

t1=nchoosek(n,2);		%total number of pairs of entities
t2=sum(sum(C.^2));	%sum over rows & columnns of nij^2
t3=.5*(nis+njs);

%Expected index (for adjustment)
nc=(n*(n^2+1)-(n+1)*nis-(n+1)*njs+2*(nis*njs)/n)/(2*(n-1));

A=t1+t2-t3;		%no. agreements
D=  -t2+t3;		%no. disagreements

if t1==nc
   AR=0;			%avoid division by zero; if k=1, define Rand = 0
else
   AR=(A-nc)/(t1-nc);		%adjusted Rand - Hubert & Arabie 1985
end

RI=A/t1;			%Rand 1971		%Probability of agreement
MI=D/t1;			%Mirkin 1970	%p(disagreement)
HI=(A-D)/t1;	%Hubert 1977	%p(agree)-p(disagree)
end
function Cont=Contingency(Mem1,Mem2)

if nargin < 2 || min(size(Mem1)) > 1 || min(size(Mem2)) > 1
   error('Contingency: Requires two vector arguments')
end

Cont=zeros(max(Mem1),max(Mem2));

for i = 1:length(Mem1);
   Cont(Mem1(i),Mem2(i))=Cont(Mem1(i),Mem2(i))+1;
end
end