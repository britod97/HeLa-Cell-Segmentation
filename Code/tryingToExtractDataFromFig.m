fig = openfig('D:\GitHub Repos\HeLa-Cell-Segmentation\NewDataCellCIL50051_1.fig');

% Find all line objects
h = findobj(fig, 'Type', 'line');

% Extract data from the first line
x = get(h(1), 'XData');
y = get(h(1), 'YData');

%%
S = load('D:\GitHub Repos\HeLa-Cell-Segmentation\NewDataCellCIL50051_1.fig', '-mat');
% Extract additional data from the loaded structure
% data = S.data; % Assuming 'data' is a field in the loaded .fig file

whos('-file', 'D:\GitHub Repos\HeLa-Cell-Segmentation\NewDataCellCIL50051_1.fig');


ax = S.hgS_070000.children(strcmp({S.hgS_070000.children.type}, 'axes'));
p = ax.children(strcmp({ax.children.type}, 'patch'));


% Extract mesh
V = p.properties.Vertices;
F = p.properties.Faces;

% Create alpha shape from the surface points
shp = alphaShape(V(:,1), V(:,2), V(:,3));

% Choose voxel size
voxelSize = 1;

minXYZ = floor(min(V,[],1));
maxXYZ = ceil(max(V,[],1));

[x,y,z] = ndgrid(minXYZ(1):voxelSize:maxXYZ(1), ...
                 minXYZ(2):voxelSize:maxXYZ(2), ...
                 minXYZ(3):voxelSize:maxXYZ(3));

% Test points inside the object
BW = inShape(shp, x, y, z);