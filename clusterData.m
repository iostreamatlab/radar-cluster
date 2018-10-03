function [clusters,clusterInds,clusterBounds] = clusterData(data, sensitivity)


% SENSITIVITY：敏感度
% 构建三维的矩阵，设置敏感度，阈值在0到1之间
% 高敏感度：分类数少
% 低敏感度：分类数多




[m,n] = size(data);
if nargin < 2
    sensitivity = 0.2*ones(1,n);
else
    if numel(sensitivity) == 1
        sensitivity = sensitivity*ones(1,n);
    end
    if numel(sensitivity) ~= n
        error('参数错误');
    end
    if any(sensitivity < 0 | sensitivity > 1)
        error('定义错误');
    end
end

tmpClusters = cell(1,n);
clusterInds = cell(1,n);
clusterBounds = cell(1,n);
for ii = 1:n
    [tmpClusters{ii},clusterInds{ii},clusterBounds{ii}] = clusterVector(data(:,ii),sensitivity(ii));
end


clusters = zeros(m,n,'single');
for ii = 1:m
    for jj = 1:n
        found = 0;
        while ~found
            for kk = 1:size(clusterInds{jj},1)
                found = ismember(data(ii,jj),tmpClusters{jj}{kk});
                if found, break, end
            end
        end
        clusters(ii,jj) = kk;
    end
end

[clusterVals, ~, clusterInds] = unique(clusters,'rows');
clusters = cell(size(clusterVals,1),1);
inds = unique(clusterInds);
for ii = 1:numel(inds)
    clusters{ii} = data(clusterInds==ii,:);
end


function [colClusters,colClusterInds,colClusterBounds] = clusterVector(vector,coarseness)

vector = vector(:);
% 分选聚类
[b,inds]  = sort(vector);
diffs     = diff([-Inf;b]);
maxDiff   = max(diffs(isfinite(diffs)));
minDiff   = coarseness * maxDiff;

startInds = find(diffs > minDiff);
endInds   = startInds - 1;
endInds   = [endInds(2:end);numel(b)];


colClusters  = cell(numel(startInds),1);
colClusterInds = colClusters;
colClusterBounds = NaN(numel(startInds),2);

for ind = 1:numel(startInds)
   tmp = startInds(ind):endInds(ind);
   colClusters{ind} = b(tmp);
   colClusterInds{ind} = inds(tmp);
   colClusterBounds(ind,1) = b(startInds(ind));
   colClusterBounds(ind,2) = b(endInds(ind));
end
