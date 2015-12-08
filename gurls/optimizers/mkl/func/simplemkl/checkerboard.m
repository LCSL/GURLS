function [X,Y] = checkerboard(nTiles,nPoints,width)

%nTiles = 10;
%nPoints = 1000;
DO_PLOTS = 0;
dim = 2;

if ~exist('width','var')
    width = sqrt(0.5) * ones(dim,1);
    widht = [0.03;0.1];
end

if size(nTiles,1) == 1
    nTiles(1:2,1) = nTiles;
end

center_x = ceil(nTiles(1)*rand(nPoints,1));
center_y = ceil(nTiles(2)*rand(nPoints,1));

Y = mod(center_x+center_y,2);
Y(Y==0) = -1;

dev = repmat(width,1,nPoints)' .* rand(nPoints,2);

X = [center_x,center_y];
X = X + dev;

if (DO_PLOTS)
figure(1);clf; hold on
ind = find(Y==1);
scatter(center_x(ind),center_y(ind),'ro');
ind = find(Y==-1);
scatter(center_x(ind),center_y(ind),'bo');

figure(2);clf; hold on
ind = find(Y==1);
scatter(X(ind,1),X(ind,2),'ro');
ind = find(Y==-1);
scatter(X(ind,1),X(ind,2),'bo');

end