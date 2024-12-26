clear; clc;close all;

addmypath

it = 450;

parfnm  ='../example/layer_model/fsg_nofilter_loc/test.json';
fnm_snap='../example/layer_model/fsg_nofilter_loc/output/volume_vel.nc';

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
elseif par.medium.equivalent_medium_method == 'tti'
    titlenm=[titlenm,' (SM calculus)'];
end

vxdata = ncread(fnm_snap, 'Vx', [1 1 it], [Inf Inf 1]);
vzdata = ncread(fnm_snap, 'Vz', [1 1 it], [Inf Inf 1]);
time = ncread(fnm_snap, 'time', [it], [1], [1]);

[vertices, vx, faces] = patch_pre(nx, nz, xvec, zvec, vxdata);
[vertices, vz, faces] = patch_pre(nx, nz, xvec, zvec, vzdata);

cal = [-1 1] * 2e3;
mycolor1 = [linspace(0,1,100)', linspace(0,1,100)', linspace(1,1,100)'];
mycolor2 = [linspace(1,1,100)', linspace(1,0,100)', linspace(1,0,100)'];
mycolor = [mycolor1; mycolor2];

figure;
subplot(2,1,1);
f = patch('Faces' ,faces, 'Vertices', vertices, 'FaceVertexCData', vx);
f.FaceColor='interp'; f.EdgeColor='none';
colormap(mycolor); caxis(cal);
axis equal; 
title({titlenm; ['Snapshot of V_x at ',num2str(time),' s']});
xlabel('x (km)'); %, 'Fontsize', 24
ylabel('z (km)');
set(gca,'layer','top');
set(gca, 'FontSize', 14);
xlim([xvec(1) xvec(end)]);
ylim([zvec(1) zvec(end)]);

subplot(2,1,2);
f = patch('Faces', faces, 'Vertices', vertices, 'FaceVertexCData', vz);
f.FaceColor='interp'; f.EdgeColor='none';
colormap(mycolor); caxis(cal);
axis equal; 
title({titlenm; ['Snapshot of V_z at ',num2str(time),' s']});
xlabel('x (km)'); %, 'Fontsize', 24
ylabel('z (km)');
set(gca,'layer','top');
set(gca, 'FontSize', 14);
xlim([xvec(1) xvec(end)]);
ylim([zvec(1) zvec(end)]);

