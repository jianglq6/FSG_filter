clc;
clear;
close all;

nt = 500; % the timestep to draw
ns = 900; % total timesteps

len = 50;
x0 = 620 - 80; 
z0 = 400 - 80;

save_figure = 0;

dh_half = 45*sqrt(2) * 1e-3;
nx = 1000;
nz = 900;

varnmx = 'Vx';
varnmz = 'Vz';

name1 = 'Collocated-grid';
name2 = 'FSG';
name4 = 'FSG-Filter';

output_dir = './';

% reference solution
file1='../example/bp_model/ref/output/volume_vel.nc';
file2='../example/bp_model/fsg_nofilter/output/volume_vel.nc';
file4='../example/bp_model/fsg_filter/output/volume_vel.nc';

startloc1 = [1 1 nt];
startloc2 = [1 1 nt];
count = [Inf Inf 1];  % [nx nz 1]

data1_x = ncread(file1, varnmx, startloc1, count);
data2_x = ncread(file2, varnmx, startloc2, count);
data4_x = ncread(file4, varnmx, startloc2, count);

data1_z = ncread(file1, varnmz, startloc1, count);
data2_z = ncread(file2, varnmz, startloc2, count);
data4_z = ncread(file4, varnmz, startloc2, count);

time2 = ncread(file2, 'time', [nt], [1], [1]);

line1_x = zeros(1, len);
line2_x = zeros(1, len);
line4_x = zeros(1, len);
line1_z = zeros(1, len);
line2_z = zeros(1, len);
line4_z = zeros(1, len);
i0 = 1;
for i = 1:len
    line1_x(i) = data1_x(x0+i-1, z0+i-1);
    line2_x(i) = data2_x(x0+i-1, z0+i-1);
    line4_x(i) = data4_x(x0+i-1, z0+i-1);
    line1_z(i) = data1_z(x0+i-1, z0+i-1);
    line2_z(i) = data2_z(x0+i-1, z0+i-1);
    line4_z(i) = data4_z(x0+i-1, z0+i-1);
end

startpos1 = 1; endpos1 = len;
startpos2 = startpos1; endpos2 = endpos1 * 0.5;

xscale1 = [startpos1:endpos1]';
xscale1 = xscale1 - xscale1(1);

% -------------------------------------------------------------------
% -              draw profile of CG-FDM and LG
% -------------------------------------------------------------------

func_figure(22, 0.3);
set(gcf, 'color', 'white', 'renderer', 'painters');

axes( 'Position', [0.06, 0.19, 0.44, 0.7] );
plot(xscale1*dh_half, line1_x,  'r-', 'linewidth', 2); hold on
plot(xscale1*dh_half, line2_x,  'b--o', 'linewidth', 2); hold on
plot(xscale1*dh_half, line4_x,  'y--o', 'linewidth', 2); hold off
xlabel('Distance (km)' ,   'Fontsize', 12);
ylabel('Amplitude', 'Fontsize', 12);
%yticks([-1 1]*1e4); 
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca, 'FontSize', 12);
legend({name1,name2,name4},'Fontsize',12,'Location','northeast','interpreter','none');
title(['V_x component at ',num2str(time2(nt)),' s'], 'Fontsize', 12);
xlim([0 (xscale1(end)*dh_half)])
%ylim([-4,4]*1e4);

axes( 'Position', [0.55, 0.19, 0.44, 0.7] );
set(gcf, 'color', 'white', 'renderer', 'painters');
plot(xscale1*dh_half, line1_z, 'r-', 'linewidth', 2); hold on
plot(xscale1*dh_half, line2_z, 'b--o', 'linewidth', 2); hold on
plot(xscale1*dh_half, line4_z, 'y--o', 'linewidth', 2); hold off
xlabel('Distance (km)' ,   'Fontsize', 12);
%ylabel('Amplitude', 'Fontsize', 12);
%yticks([-1 0 1]*1e4); 
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca, 'FontSize', 12);
legend({name1,name2,name4},'Fontsize',12,'Location','northeast','interpreter','none');
title(['V_z component at ',num2str(time2(nt)),' s'], 'Fontsize', 12);
xlim([0 (xscale1(end)*dh_half)])
%ylim([-4,4]*1e4);





