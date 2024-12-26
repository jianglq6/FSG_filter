function fid = func_figure( wid, hwratio, render )

%=====================================================
% figure setting
%=====================================================
Winch=wid; 

Hinch = Winch * hwratio;
if Hinch > 50
    disp(['Height=',num2str(Hinch),' exceeds the limit, reset to 50cm']);
    Hinch=50; %- largest dep (AGU)
end

Wline=0.5;

fid = figure;
%=== painters for eps file, zbuffer for bitmap file ===
if exist('render','var')
    set(gcf,'renderer',render);
else
    set(gcf,'renderer','painters');
end

%=== paper size ===
%set(gcf,'PaperPositionMode','auto'); % the same size as on screen and centered
set(gcf,'PaperPositionMode','manual'); % honors the PaperPosition
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperSize',[Winch Hinch]);    
%=== printed figure to lower left corner of the page and width, height ===
set(gcf,'PaperPosition',[0 0 Winch Hinch]);
%=== figure size on screen ===
set(gcf,'Units','centimeters');
set(gcf,'Position',[2 2 Winch Hinch]);
%=== retain the background color ===
%set(gcf,'color','blue');  % bg color of figure
%set(gca,'color','red');   % bg color of axes
%set(gcf,'InvertHardCopy','off');

set(gcf,'defaultaxesfontname','Helvetica');
%set(gcf,'defaultaxesfontname','Times New Roman');
set(gcf,'defaultaxesfontsize',14);
set(gcf,'defaulttextfontsize',14);
set(gcf,'defaultaxeslinewidth',0.8);
set(gcf,'defaultlinelinewidth',0.8);


