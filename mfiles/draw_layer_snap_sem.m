clear all;

%=== grid info ===
% same as fd
x0 = 0.0;   z0 = 0.0;
dx = 100.0; dz = 100.0;
nx = 500;   nz = 500;
xvec = [0:nx-1]*dx + x0;
zvec = [0:nz-1]*dz + z0;
[X, Z] = meshgrid(xvec, zvec);
cal = [-1 1] * 2e3;

mycolor1 = [linspace(0,1,100)', linspace(0,1,100)', linspace(1,1,100)'];
mycolor2 = [linspace(1,1,100)', linspace(1,0,100)', linspace(1,0,100)'];
mycolor = [mycolor1; mycolor2];

sem_path = ['../../specfem2d/model1_2layer/OUTPUT_FILES/'];

coord_path=[sem_path,'wavefield_grid_for_dumps.txt'];
wave_path =[sem_path,'wavefield0009200_01.txt'];

coord = load(coord_path);
wave  = load(wave_path);

vx_sem = scatteredInterpolant(coord(:,1), coord(:,2), wave(:,1));
vz_sem = scatteredInterpolant(coord(:,1), coord(:,2), wave(:,2));

figure;
subplot(2,1,1);
imagesc(xvec, zvec, vx_sem(X,Z));
shading interp;
shading interp;
colormap(mycolor); caxis(cal);
axis equal;
set(gca,'ydir','normal');
xlabel('x (km)'); 
ylabel('z (km)');
set(gca,'layer','top');
set(gca, 'FontSize', 14);
xlim([xvec(1) xvec(end)]);
ylim([zvec(1) zvec(end)]);

subplot(2,1,2);
imagesc(xvec, zvec, vz_sem(X, Z));
shading interp;
colormap(mycolor); caxis(cal);
axis equal;
set(gca,'ydir','normal');
xlabel('x (km)'); 
ylabel('z (km)');
set(gca,'layer','top');
set(gca, 'FontSize', 14);
xlim([xvec(1) xvec(end)]);
ylim([zvec(1) zvec(end)]);
