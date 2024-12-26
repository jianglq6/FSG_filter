clc;
clear;
%close all;

flag_save = 1;

mtype = 'layer_';
ftype = 'profile_';

nt = 450; % the timestep to draw
ns = 700; % total timesteps

len = 60;
x0 = 600-30;
z0 = 620-30;

% coord
dx = 50*1e-3;
dz = 50*1e-3;
nx = 1000; nz = 1000;
xx = (0:(nx-1)) * dx;
zz = (0:(nz-1)) * dz;

[X, Z] = meshgrid(xx, zz);

xl1 = x0;  zl1 = z0;
xl2 = x0 + len; zl2 = z0 + len;
xl1 = xl1*dx; zl1 = zl1*dx;
xl2 = xl2*dx; zl2 = zl2*dx;


dd = 50*sqrt(2) * 1e-3;
nx = 1000;
nz = 1000;
interv = 20;

varnmx = 'Vx';
varnmz = 'Vz';

name1 = 'Collocated-grid';
name2 = 'FSG';
name3 = 'FSG (SM calculus)';
name4 = 'FSG-Filter';

cal = [-1 1] * 1e3;

output_dir = './';

file2='../example/layer_model/fsg_nofilter_loc/output/volume_vel.nc';
file3='../example/layer_model/fsg_nofilter_tti/output/volume_vel.nc';
file4='../example/layer_model/fsg_filter_tti/output/volume_vel.nc';

sem_path = ['../../specfem2d/model1_2layer/OUTPUT_FILES/'];
sem_coord=[sem_path,'wavefield_grid_for_dumps.txt'];
sem_wave =[sem_path,'wavefield0009200_01.txt'];
coord = load(sem_coord); coord = coord/1000;
wave  = load(sem_wave);


startloc1 = [1 1 nt];
startloc2 = [1 1 nt];
count = [Inf Inf 1]; % [nx nz 1]

data2_x = ncread(file2, varnmx, startloc2, count);
data3_x = ncread(file3, varnmx, startloc2, count);
data4_x = ncread(file4, varnmx, startloc2, count);

data2_z = ncread(file2, varnmz, startloc2, count);
data3_z = ncread(file3, varnmz, startloc2, count);
data4_z = ncread(file4, varnmz, startloc2, count);


time_2 = ncread(file2, 'time');
time2 = time_2(1:ns);

vx2 = data2_x(:,:,1);
vx3 = data3_x(:,:,1);
vx4 = data4_x(:,:,1);

vz2 = data2_z(:,:,1);
vz3 = data3_z(:,:,1);
vz4 = data4_z(:,:,1);

% interp
vx_sem = scatteredInterpolant(coord(:,1), coord(:,2), wave(:,1));
vz_sem = scatteredInterpolant(coord(:,1), coord(:,2), wave(:,2));

name1 = 'SEM';
name2 = 'FSG';
name2 = 'FSG (SM calculus)';
name4 = 'FSG-Filter';

line2_x = zeros(1, len);
line4_x = zeros(1, len);
line2_z = zeros(1, len);
line4_z = zeros(1, len);
line4sem_x=[];
line4sem_z=[];
i0 = 1;
for i = 1:len
    line2_x(i) = vx2(x0+i-1, z0+i-1); 
    line3_x(i) = vx3(x0+i-1, z0+i-1);  
    line4_x(i) = vx4(x0+i-1, z0+i-1);
    line2_z(i) = vz2(x0+i-1, z0+i-1);
    line3_z(i) = vz3(x0+i-1, z0+i-1);
    line4_z(i) = vz4(x0+i-1, z0+i-1);
    line4sem_x(i) = xx(x0+i-1);
    line4sem_z(i) = zz(z0+i-1);
end

startpos1 = 1; endpos1 = len;
xscale1 = [startpos1:endpos1]';
xscale1 = xscale1 - xscale1(1);

% ===== for sem ===
line_vx_sem = vx_sem(line4sem_x, line4sem_z);
line_vz_sem = vz_sem(line4sem_x, line4sem_z);

% -------------------------------------------------------------------
% -              draw profile of CG-FDM and LG-Filter
% -------------------------------------------------------------------
func_figure(22, 0.5);
set(gcf, 'color', 'white', 'renderer', 'painters');

%======= vx =======
axes( 'Position', [0.06, 0.12+0.37, 0.44, 0.43] );
plot(xscale1*dd, line_vx_sem, 'r-', 'linewidth', 2); hold on;
plot(xscale1*dd, line3_x, 'b--o', 'linewidth', 2); hold on;
plot(xscale1*dd, line4_x, 'y--o', 'linewidth', 2); hold off;
ylabel('Amplitude', 'Fontsize', 12);
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca, 'FontSize', 12);
legend({name1,name3,name4},'Fontsize',12,'Location','northeast','interpreter','none');
title(['V_x component at ',num2str(time2(nt)),' s'], 'Fontsize', 12);
xlim([0 (xscale1(end)*dd)]);
xticklabels({});
ax = gca;
ax.YAxis.Exponent = 4;
ylim([-1.1, 2]*1e4);

%--- diff
axes( 'Position', [0.06, 0.12, 0.44, 0.3] );
plot(xscale1*dd, line_vx_sem-line3_x, 'b--o', 'linewidth', 2); hold on;
plot(xscale1*dd, line_vx_sem-line4_x, 'y--o', 'linewidth', 2); hold off;
xlabel('Distance (km)' ,   'Fontsize', 12);
ylabel('diff', 'Fontsize', 12);
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca, 'FontSize', 12);
legend({[name1,'-[', name3,']'],[name1,'-[',name4,']']},...
    'Fontsize',12,'interpreter','none', 'Location','southwest');
xlim([0 (xscale1(end)*dd)])
ax = gca;
ax.YAxis.Exponent = 3;
ylim([-0.7, 0.3]*1e4);


%====== vz =======
axes( 'Position', [0.55, 0.12+0.37, 0.44, 0.43] );
set(gcf, 'color', 'white', 'renderer', 'painters');
plot(xscale1*dd, line_vz_sem,  'r-', 'linewidth', 2); hold on;
plot(xscale1*dd, line3_z, 'b--o', 'linewidth', 2); hold on;
plot(xscale1*dd, line4_z, 'y--o', 'linewidth', 2); hold off;
set(gca,'LooseInset',get(gca,'TightInset'));
xticklabels({});
set(gca, 'FontSize', 12);
legend({name1,name3,name4},'Fontsize',12,'Location','northeast','interpreter','none');
title(['V_z component at ',num2str(time2(nt)),' s'], 'Fontsize', 12);
xlim([0 (xscale1(end)*dd)])
ylim([-1.1, 2]*1e4);

%--- diff
axes( 'Position', [0.55, 0.12, 0.44, 0.3] );
plot(xscale1*dd, line_vz_sem-line3_z, 'b--o', 'linewidth', 2); hold on;
plot(xscale1*dd, line_vz_sem-line4_z, 'y--o', 'linewidth', 2); hold off;
xlabel('Distance (km)' ,   'Fontsize', 12);
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca, 'FontSize', 12);
legend({[name1,'-[', name3,']'],[name1,'-[',name4,']']},...
    'Fontsize',12,'interpreter','none', 'Location','southwest');
xlim([0 (xscale1(end)*dd)])
ax = gca;
ax.YAxis.Exponent = 3;
ylim([-0.7, 0.3]*1e4);

% -------------------------------------------------------------------
% -              save profile fig
% -------------------------------------------------------------------

if flag_save == 1
   print('./planar_wave_diag_sem','-dpng', '-r300');
end
    



