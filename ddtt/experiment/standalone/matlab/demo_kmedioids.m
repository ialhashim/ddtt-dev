function [M, shapes, clustering, centers] = demo_kmedioids(filename,shapes_filename)
    close all;

    % Load correspondence results
    [M, shapes]=loadCorrespondence(filename,shapes_filename);

    % Convert distance to similairty
    %s = 1-M;

    % Cluster
    num_clusters = 5;
    [inds,centers] = kmedioids(M, num_clusters);
    clusters = inds;
    for i=1:size(inds, 2)
        clusters(i) = centers(inds(i));
    end
    
    % Visualize
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

    % Setup thumbnails for clusters
    cols = cell(k,1);
    for c=1:k 
        cur_cluster = clustering{c}.elements;

        idx_representative = centers(c);
        class_representative = makeThumb(shapes{ idx_representative }.thumb);
        length_cur_cluster = length(cur_cluster);
        
        row = class_representative;
        
        for i=1:length_cur_cluster
            cur_idx = cur_cluster{i};
            cluster_element = shapes{ cur_idx };
            img = makeThumb(cluster_element.thumb);
            if cur_idx == idx_representative
                continue
            end
            row = [row img];
        end
        
        cols{c} = row;
    end

    % 
    for i = 1:k
        subplot(k,1,i);
        imshow(cols{i});
        axis off; 
    end 
    fclose('all');
end

function thumb_img = makeThumb(filename)
    thumb_img = imresize(imread(filename), [NaN 128]);
end
