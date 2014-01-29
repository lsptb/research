% 2-D Projection, streamlines, & origin calculation
% TODO: redo streamline & origin calculations in 3-D without projection.

%% POSITION & DELAY MAP

% PREPARE POSITION
nchan = data.num_sensors;
T     = data.coor_trans.device2head;
pos   = zeros(nchan,3); % (x,y,z) for each channel
for k = 1:nchan
  loc         = data.sensor_info(k).loc;
  loc         = T*loc;
  pos(k,1:3)  = loc(1:3,4);
end
method = 'stereographic'; % gnomic, stereographic, ortographic, inverse, polar
prj    = elproj(pos, method); % * [0 1; -1 0];
          % ELPROJ makes a azimuthal projection of a 3D electrode cloud
          %  on a plane tangent to the sphere fitted through the electrodes
          %  the projection is along the z-axis
X = prj(:,1);   % x-coordinates
Y = prj(:,2);   % y-coordinates
z = [];         % delay vector (one value per channel)

% Scale the data to a circle with x-axis and y-axis: -0.45 to 0.45
y = 0.9*((X-min(X))/(max(X)-min(X))-0.5); % NOTE: x becomes y and y becomes x, in griddata is also reversed which makes x x and y y again
x = 0.9*((Y-min(Y))/(max(Y)-min(Y))-0.5);

% DELAY MAP INTERPOLATION
interplimits = 'head';
  % 'electrodes' to furthest electrode
  % 'head' to edge of head
% Find limits for interpolation:
if strcmp(interplimits,'head')
  xmin = min(-.5,min(x)); xmax = max(0.5,max(x));
  ymin = min(-.5,min(y)); ymax = max(0.5,max(y));
else
  xmin = max(-.5,min(x)); xmax = min(0.5,max(x));
  ymin = max(-.5,min(y)); ymax = min(0.5,max(y));
end

gridscale  = 100;                             % resolution
interp     = 'v4';                            % 'linear','cubic','nearest','v4'
xi         = linspace(xmin,xmax,gridscale);   % x-axis description (row vector)
yi         = linspace(ymin,ymax,gridscale);   % y-axis description (row vector)
% calculate delay map
[Xi,Yi,Zi] = griddata(y,x,z,yi',xi,interp);  % Interpolate data; NOTE: undo the reversal of x & y
  % griddata uses meshgrid to create evenly-spaced XI & YI

%% VELOCITY FIELD & STREAMLINES
% calculate gradient of delay map
[Dx Dy] = gradient(Zi,Xi',Yi);

% calculate components of the velocity field
Vx = gradient(Xi)./Dx;
Vy = gradient(Yi)./Dy;

% calculate streamlines for each channel
figure
wavenum = 1;
for k = 1:nchan
  x0  = y(k); % RECALL: x & y were switched before regular spacing
  y0  = x(k); 
  waves(wavenum).streamlines(k).label     = '';
  waves(wavenum).streamlines(k).vertices  = stream2(XI,YI,Vx,Vy,x0,y0);
  h   = streamline(XI,YI,Vx,Vy,x0,y0);
  set(h,'LineWidth',2,'Color','r')
  title(waves(wavenum).streamlines(k).label)
  pause
end

%% ORIGIN CALCULATION

% for each wave:
% find the longest streamline
% if it occurs within 20ms of the earliest peak, it's IC is the origin


%% TOPOPLOTS
% Take data within head
rmax       = .5;
mask       = (sqrt(Xi.^2+Yi.^2) <= rmax);
Zi(mask==0)= NaN;

% Draw topoplot on head
numcontours = 6;
shading     = 'flat';                         % 'flat' or 'interp'
delta       = xi(2)-xi(1);
contour(Xi,Yi,Zi,numcontours,'k');
h = surface(Xi-delta/2,Yi-delta/2,zeros(size(Zi)),Zi,'EdgeColor','none', 'FaceColor',shading);

% calculate colormap limits
zmin = min(Zi(:));
zmax = max(Zi(:));
caxis([zmin zmax])
    
% Masking data
% use griddata to create mask
% set FaceAlpha to interp and AlphaData to the mask
    
% draw electrodes
    hp2    = plot(y(normal),        x(normal),        cfg.emarker,  'Color', cfg.ecolor,  'markersize', cfg.emarkersize);
    hp2    = plot(y(cfg.highlight), x(cfg.highlight), cfg.hlmarker, 'Color', cfg.hlcolor, 'markersize', cfg.hlmarkersize, ...
% add labels
    text(y(i), x(i), labels(i,:), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                                  'Color', cfg.ecolor, 'FontSize', cfg.efsize);
% draw head
  % Define the outline of the head, ears and nose:
  l     = 0:2*pi/100:2*pi;
  tip   = rmax*1.15; base = rmax-.004;
  EarX  = [.497 .510 .518 .5299 .5419 .54 .547 .532 .510 .489];
  EarY  = [.0555 .0775 .0783 .0746 .0555 -.0055 -.0932 -.1313 -.1384 -.1199];
  % Plot head, ears, and nose:
  plot(cos(l).*rmax, sin(l).*rmax, 'color', cfg.hcolor, 'Linestyle', '-', 'LineWidth', cfg.hlinewidth);
  plot([0.18*rmax;0;-0.18*rmax], [base;tip;base], 'Color', cfg.hcolor, 'LineWidth', cfg.hlinewidth);
  plot( EarX, EarY, 'color', cfg.hcolor, 'LineWidth', cfg.hlinewidth)
  plot(-EarX, EarY, 'color', cfg.hcolor, 'LineWidth', cfg.hlinewidth)

colorbar
hold off
axis off
xlim([-.6 .6]);
ylim([-.6 .6]);

% Add streamlines
