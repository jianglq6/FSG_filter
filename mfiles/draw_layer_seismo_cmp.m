clc;
clear;
%close all;

nx = 1000;
nz = 1000;

varnmx = 'Vx';
varnmz = 'Vz';

name1 = 'SEM';
name2 = 'FSG';
name3 = 'FSG (SM calculus)';
name4 = 'FSG-Filter';

file2='../example/layer_model/fsg_nofilter_loc/output/volume_vel.nc';
file3='../example/layer_model/fsg_nofilter_tti/output/volume_vel.nc';
file4='../example/layer_model/fsg_filter_tti/output/volume_vel.nc';

fnm = ['../../specfem2d/model1_2layer/OUTPUT_FILES/AA.S0001.BXX.semv'];
sem_v = load(fnm);
sem_vx = sem_v(:,2);

fnm = ['../../specfem2d/model1_2layer/OUTPUT_FILES/AA.S0001.BXZ.semv'];
sem_v = load(fnm);
sem_vz = sem_v(:,2);
sem_t = sem_v(:,1) + 0.5;  % sem need to + time shift t0

% receiver
x = 601; z = 681; rnm = 'R';

startloc = [x z 1];

NT1 = 1; % the start timestep to draw
NT2 = 690; % the end timestep to draw
count = [1 1 Inf]; 

dt = 0.01;
tmin = dt*NT1;
tmax = dt*NT2;

% zoomed-in view
xt1 = 4.7;
xt2 = 5.5;
zt1 = 4.7;
zt2 = 5.5;

data2_x = ncread(file2, varnmx, startloc, count);
data3_x = ncread(file3, varnmx, startloc, count);
data4_x = ncread(file4, varnmx, startloc, count);
data2_z = ncread(file2, varnmz, startloc, count);
data3_z = ncread(file3, varnmz, startloc, count);
data4_z = ncread(file4, varnmz, startloc, count);

time2 = ncread(file2, 'time');
time3 = ncread(file3, 'time');
time4 = ncread(file4, 'time');

locax = (x-1) * 50 * 1e-3;
locaz = (z-1) * 50 * 1e-3;

rev2_x = reshape(data2_x(:,:,NT1:NT2), [1, length(time2(NT1:NT2))]);
rev3_x = reshape(data3_x(:,:,NT1:NT2), [1, length(time2(NT1:NT2))]);
rev4_x = reshape(data4_x(:,:,NT1:NT2), [1, length(time4(NT1:NT2))]);
                                                              
rev2_z = reshape(data2_z(:,:,NT1:NT2), [1, length(time2(NT1:NT2))]);
rev3_z = reshape(data3_z(:,:,NT1:NT2), [1, length(time2(NT1:NT2))]);
rev4_z = reshape(data4_z(:,:,NT1:NT2), [1, length(time4(NT1:NT2))]);

[~, it1_sem] = min(abs(sem_t-time2(NT1)));    % t1(NT1) = 0.005
[~, it2_sem] = min(abs(sem_t-time2(NT1+1)));
intv_sem = it2_sem-it1_sem;
[~, it2_sem] = min(abs(sem_t-time2(NT2)));

% -------------------------------------------------------------------
% -              draw rev of CG-FDM and LG-Filter
% -------------------------------------------------------------------
func_figure(25, 0.45);
set(gcf, 'color', 'white', 'renderer', 'painters');

axes( 'Position', [0.054, 0.17+0.23, 0.45, 0.53] );
plot(sem_t(it1_sem:intv_sem:it2_sem), sem_vx(it1_sem:intv_sem:it2_sem), 'k-', 'linewidth', 2); hold on;
plot(time2(NT1:NT2), rev2_x, 'b--', 'linewidth', 2); hold on;
plot(time2(NT1:NT2), rev3_x, '--','color',[0.4940 0.1840 0.5560], 'linewidth', 2); hold on;
plot(time4(NT1:NT2), rev4_x, 'r--', 'linewidth', 2); hold off;
ylabel('Amplitude', 'Fontsize', 12);
set(gca, 'FontSize', 12);
xticklabels({});
title('V_x component', 'Fontsize', 12);
xlim([tmin tmax]);
ylim([-7 5]*1e4);

%zoom
axes( 'Position', [0.11, 0.25+0.21, 0.12, 0.18] );
set(gcf, 'color', 'white', 'renderer', 'painters');
plot(sem_t(it1_sem:intv_sem:it2_sem), sem_vx(it1_sem:intv_sem:it2_sem), 'k-', 'linewidth', 2); hold on;
plot(time2(NT1:NT2), rev2_x, 'b--', 'linewidth', 2); hold on;
plot(time3(NT1:NT2), rev3_x, '--','color',[0.4940 0.1840 0.5560], 'linewidth', 2); hold on;
plot(time4(NT1:NT2), rev4_x, 'r--', 'linewidth', 2); hold off
set(gca, 'LooseInset', [0.0, 0.0, 0.01, 0.0]);
set(gca, 'FontSize', 12);
xlim([xt1 xt2]);
ylim([-3500, 4000]);


%======= diff
axes( 'Position', [0.054, 0.115, 0.45, 0.2] );
plot(time2(NT1:NT2), sem_vx(it1_sem:intv_sem:it2_sem)'-rev2_x, 'b-', 'linewidth', 2); hold on;
plot(time2(NT1:NT2), sem_vx(it1_sem:intv_sem:it2_sem)'-rev3_x, '-','color',[0.4940 0.1840 0.5560], 'linewidth', 2); hold on;
plot(time4(NT1:NT2), sem_vx(it1_sem:intv_sem:it2_sem)'-rev4_x, 'r-', 'linewidth', 2); hold off;
xlabel('Time (s)' ,   'Fontsize', 12);
ylabel('diff', 'Fontsize', 12);
set(gca, 'FontSize', 12);
xlim([tmin tmax]);
ylim([-1 1]*1e4);



%============ vz ============
axes( 'Position', [0.54, 0.17+0.23, 0.45, 0.53] );
plot(sem_t(it1_sem:intv_sem:it2_sem), sem_vz(it1_sem:intv_sem:it2_sem), 'k-', 'linewidth', 2); hold on;
plot(time2(NT1:NT2), rev2_z, 'b--', 'linewidth', 2); hold on
plot(time3(NT1:NT2), rev3_z, '--','color',[0.4940 0.1840 0.5560], 'linewidth', 2); hold on;
plot(time4(NT1:NT2), rev4_z, 'r--', 'linewidth', 2); hold off
%xlabel('Time (s)' ,   'Fontsize', 12);
set(gca, 'FontSize', 12);
legend({name1,name2,name3,name4},'Fontsize',12,'Location','southeast','interpreter','none');
title('V_z component', 'Fontsize', 12);
xlim([tmin tmax]);
ylim([-7 5]*1e4);
xticklabels({});


axes( 'Position', [0.6, 0.25+0.21, 0.12, 0.18] );
set(gcf, 'color', 'white', 'renderer', 'painters');
plot(sem_t(it1_sem:intv_sem:it2_sem), sem_vz(it1_sem:intv_sem:it2_sem), 'k-', 'linewidth', 2); hold on;
plot(time2(NT1:NT2), rev2_z, 'b--', 'linewidth', 2); hold on
plot(time2(NT1:NT2), rev3_z, '--','color',[0.4940 0.1840 0.5560], 'linewidth', 2); hold on;
plot(time4(NT1:NT2), rev4_z, 'r--', 'linewidth', 2); hold off
set(gca, 'LooseInset', [0.0, 0.0, 0.01, 0.0]);
set(gca, 'FontSize', 12);
xlim([xt1 xt2]);
ylim([-3500, 4000]);


%======= diff
axes( 'Position', [0.54, 0.115, 0.45, 0.2] );
plot(time2(NT1:NT2), sem_vz(it1_sem:intv_sem:it2_sem)'-rev2_z, 'b-', 'linewidth', 2); hold on
plot(time3(NT1:NT2), sem_vz(it1_sem:intv_sem:it2_sem)'-rev3_z, '-','color',[0.4940 0.1840 0.5560], 'linewidth', 2); hold on;
plot(time4(NT1:NT2), sem_vz(it1_sem:intv_sem:it2_sem)'-rev4_z, 'r-', 'linewidth', 2); hold off
xlabel('Time (s)' ,   'Fontsize', 12);
set(gca, 'FontSize', 12);
legend({[name1,'-[',name2,']'],[name1,'-[',name3,']'],[name1,'-[',name4,']']},'Fontsize',12,'Location','southeast','interpreter','none');
set(gca, 'FontSize', 12);
xlim([tmin tmax]);
ylim([-1 1]*1e4);



    




