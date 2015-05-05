% Load correspondence results
%[M, shapes]=loadCorrespondence;

% Convert distance to similairty
s = 1-M;

% Preference 
p = median(1-M);

% Cluster
clusters = apcluster(s,p);

% Visualize
centers = unique(clusters);
k = length(centers);
clustering = cell(1,k);

% Collect shape indicies per cluster
for i=1:k
    cluster_center = centers(i);
    ele_indices = find(clusters == cluster_center);
    w = length(ele_indices);
    clustering{i}.elements = cell(1,w);
    for j=1:w
        % Index of shape belonging to cluster 'i'
        idx = ele_indices(j);        
        clustering{i}.elements{j} = idx;
    end
end

% Display clusters
for c=1:k 
    cur_cluster = clustering{c}.elements;
    representative = imresize(shapes{ centers(i) }, [NaN 128]);
    for i=1:length(cur_cluster)
        cluster_element = shapes{ cur_cluster{i} };
        img = imread(cluster_element.thumb);
    end
end

