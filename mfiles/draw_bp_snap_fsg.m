clear; clc;

addmypath;

it = 500;

% Magnify
x_1 = 17;
x_2 = 21;
y_1 = 6;
y_2 = 9;


parfnm  ='../example/bp_model/fsg_nofilter/test.json';
fnm_snap='../example/bp_model/fsg_nofilter/output/volume_vel.nc';

%==== get grid info from par =======
if ~exist(parfnm, 'file')
    error([parfnm, ' does not exist']);
end
par = loadjson(parfnm);
nx = par.number_of_total_grid_points_x * 2;
nz = par.number_of_total_grid_points_z * 2;
dh = par.grid_generation_method.cartesian.inteval;
dx = dh(1)/2000;  % half grid size (km)
dz = dh(2)/2000;
xvec = [0:nx-1]*dx; zvec = [0:nz-1]*dz;
titlenm = ['FSG'];
if par.is_filter == 1
    titlenm=[titlenm,'-Filter'];
end

vxdata = ncread(fnm_snap, 'Vx', [1 1 it], [Inf Inf 1]);
vzdata = ncread(fnm_snap, 'Vz', [1 1 it], [Inf Inf 1]);
time = ncread(fnm_snap, 'time', [it], [1], [1]);

[vertices, vx, faces] = patch_pre(nx, nz, xvec, zvec, vxdata);
[vertices, vz, faces] = patch_pre(nx, nz, xvec, zvec, vzdata);

cal = [-1 1] * 3e4;
mycolor1 = [linspace(0,1,100)', linspace(0,1,100)', linspace(1,1,100)'];
mycolor2 = [linspace(1,1,100)', linspace(1,0,100)', linspace(1,0,100)'];
mycolor = [mycolor1; mycolor2];

interv = 5;
ztick = [floor(zvec(end)):-interv:0];
ztick_zoom = [floor(zvec(end))-6:-1:floor(zvec(end))-9];

func_figure(14,2);
set( gcf, 'Color', 'White');

nSeries = 2;
% - Compute #rows/cols, dimensions, and positions of lower-left corners.
nCol = 1;  nRow = ceil( nSeries / nCol );
rowH = 0.8 / nRow ;  
rowY = 0.05 + linspace( 1, 0.0, nRow+1 );  rowY = rowY(2:end);
colW = 0.83 / nCol;
colX = 0.11 + linspace( 0.00, 1, nCol+1 );  colX = colX(1:end-1);
zoomwidth = colW/2.3;
zoomheigh = rowH/2.3;


dId = 1;
rowId = ceil( dId / nCol );
colId = dId - (rowId - 1) * nCol;
axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] );

f = patch('Faces' ,faces, 'Vertices', vertices, 'FaceVertexCData', vx);
f.FaceColor='interp'; f.EdgeColor='none';
colormap(mycolor); caxis(cal);
axis equal; 
title({titlenm; ['Snapshot of V_x at ',num2str(round(time,2)),' s']});
xlabel('x (km)'); %, 'Fontsize', 24
ylabel('z (km)');
box on;
set(gca,'layer','top');
set(gca, 'FontSize', 14);
set(gca, 'YDir', 'reverse');
xticks(0:interv:xvec(end)); 
yticks(0:interv:zvec(end));
yticklabels(ztick);
xlim([xvec(1) xvec(end)]);
ylim([zvec(1) zvec(end)]);

%= zoom
axes( 'Position', [colX(colId)+0.07, rowY(rowId)+0.025, zoomwidth, zoomheigh] );
f = patch('Faces' ,faces, 'Vertices', vertices, 'FaceVertexCData', vx);
f.FaceColor='interp'; f.EdgeColor='none';colormap(mycolor); caxis(cal);
axis equal; 
set(gca, 'YDir', 'reverse');
set(gca,'layer','top');
box on;
xticks(17:21);
yticks(6:9);
yticklabels(ztick_zoom);
xlim([x_1 x_2]);
ylim([y_1 y_2]);


dId = 2;
rowId = ceil( dId / nCol );
colId = dId - (rowId - 1) * nCol;
axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] );
f = patch('Faces', faces, 'Vertices', vertices, 'FaceVertexCData', vz);
f.FaceColor='interp'; f.EdgeColor='none';
colormap(mycolor); caxis(cal);
axis equal; 
box on;
title({titlenm; ['Snapshot of V_z at ',num2str(round(time,2)),' s']});
xlabel('x (km)'); %, 'Fontsize', 24
ylabel('z (km)');
set(gca, 'YDir', 'reverse');
set(gca,'layer','top');
set(gca, 'FontSize', 14);
xticks(0:interv:xvec(end)); 
yticks(0:interv:zvec(end));
yticklabels(ztick);
xlim([xvec(1) xvec(end)]);
ylim([zvec(1) zvec(end)]);

%= zoom
axes( 'Position', [colX(colId)+0.07, rowY(rowId)+0.025, zoomwidth, zoomheigh] );
f = patch('Faces' ,faces, 'Vertices', vertices, 'FaceVertexCData', vz);
f.FaceColor='interp'; f.EdgeColor='none';colormap(mycolor); caxis(cal);
axis equal; 
set(gca, 'YDir', 'reverse');
set(gca,'layer','top');
box on;
xticks(17:21);
yticks(6:9);
yticklabels(ztick_zoom);
xlim([x_1 x_2]);
ylim([y_1 y_2]);
