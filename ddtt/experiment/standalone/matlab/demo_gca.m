function [M, shapes, clustering, centers] = demo_gca(filename,shapes_filename)
    close all;

    % Load correspondence results
    [M, shapes]=loadCorrespondence(filename,shapes_filename);

    % Convert distance to similairty
    %s = 1-M;

    % Cluster
    k = 5;
    clusters = gacCluster(M, k);
    %clusters = gacCluster(M, k, 'zeta');
    %clusters = gacCluster(M, k, 'path');
    
    % Visualize
    clustering = cell(1,k);
    
    for i=1:k
        clustering{ i }.elements = [];
    end
    
    for i=1:size(clusters)
        class = clusters(i);
        clustering{ class }.elements(end+1) = i;
    end

    % Setup thumbnails for clusters
    cols = cell(k,1);
    for c=1:k 
        cur_cluster = clustering{c}.elements;

        length_cur_cluster = length(cur_cluster);
        
        row = [];
        
        for i=1:length_cur_cluster
            cur_idx = cur_cluster(i);
            cluster_element = shapes{ cur_idx };
            img = makeThumb(cluster_element.thumb);
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
