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

% Load mesh from the patch object
V = p.properties.Vertices;
F = p.properties.Faces;

% Make mesh structure in the format VOXELISE expects
meshFV.faces = F;
meshFV.vertices = V;

% Check mesh bounds
xmin = min(V(:,1));
xmax = max(V(:,1));

ymin = min(V(:,2));
ymax = max(V(:,2));

zmin = min(V(:,3));
zmax = max(V(:,3));

disp([xmin ymin zmin])
disp([xmax ymax zmax])

% Define voxel resolution
% Increase these numbers for finer reconstruction
nx = 512;

ny = round(nx * (ymax-ymin)/(xmax-xmin));
nz = round(nx * (zmax-zmin)/(xmax-xmin));

% Create voxel coordinate vectors
gridX = linspace(xmin,xmax,nx);
gridY = linspace(ymin,ymax,ny);
gridZ = linspace(zmin,zmax,nz);


% Ensure they are clean row vectors
gridX = unique(gridX(:))';
gridY = unique(gridY(:))';
gridZ = unique(gridZ(:))';

% Voxelise mesh
BW_1 = VOXELISE(gridX,gridY,gridZ,meshFV,'xyz');


% Create alpha shape from the surface points
% shp = alphaShape(V(:,1), V(:,2), V(:,3));

% Choose voxel size
% voxelSize = 1;
% 
% minXYZ = floor(min(V,[],1));
% maxXYZ = ceil(max(V,[],1));
% 
% [x,y,z] = ndgrid(minXYZ(1):voxelSize:maxXYZ(1), ...
%                  minXYZ(2):voxelSize:maxXYZ(2), ...
%                  minXYZ(3):voxelSize:maxXYZ(3));

% Test points inside the object
% BW = inShape(shp, x, y, z);

% volshow(BW);

% Save the entire volume as a .mat
save("D:\GitHub Repos\HeLa_Cell_Data\CIL50051\Matlab\ChlamydiaCell1.mat", "BW_1")

%%

outputDir = "D:\GitHub Repos\HeLa_Cell_Data\CIL50051\Segmented Tiffs";
% Create output directory if it doesn't exist
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
% Save the slices of the binary volume as Tif files
for k = 1:size(BW_1,3)
    filename = fullfile(outputDir, sprintf('slice_%03d.tif', k));
    imwrite(uint8(BW_1(:,:,k))*255, filename);
end