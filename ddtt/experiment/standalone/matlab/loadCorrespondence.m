function [M, shapes]=loadCorrespondence(filename,shapes_filename)

if ~exist('filename', 'var')
    filename = 'dataset.corr';
end    

if ~exist('shapes_filename', 'var')
    shapes_filename = 'shapes.txt';
end    

% Load correspondence file
tic
opt = {};
opt.FastArrayParser = 1;
opt.ShowProgress = 1;
json = loadjson(filename, opt);
toc

% Load shape names
fid = fopen(shapes_filename);
C = textscan(fid, '%d %s');
shapes = cell(1, length(C{2}));
names = C{2};
for i=1:length(length(C{2}))
    shapes{i}.name = names{i};
end

num_shapes = length(shapes);
n = num_shapes;

% Similairty matrix
M = zeros(n,n);

% Fill matrix with correspondence scores
for i=1:length(json)
    pair = json{i};
    idx_i = pair.i + 1;
    idx_j = pair.j + 1;
    cost = pair.cost;

    M(idx_i,idx_j) = cost;
    M(idx_j,idx_i) = cost;
    
    shapes{idx_i}.graph = pair.source;
    shapes{idx_i}.thumb = [pair.source(1:length(pair.source)-3) 'png'];
end

% Normalize
%mmin = min(M(:));
%mmax = max(M(:));
%M = (M-mmin) ./ (mmax-mmin);

end
