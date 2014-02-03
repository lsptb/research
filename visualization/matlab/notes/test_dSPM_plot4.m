tstart = tic;

if ~any(findstr('/home/jsherfey/svn/dev/onestream',path))
  addpath(genpath('/home/jsherfey/svn/dev/onestream'));
end

% -----------------------------------------------------
c          = 1; % condition index (not event code)
prefix     = 'test_dSPM_WordNPNW_grad_bem';
rootoutdir = '/home/jsherfey/svn/dev/onestream/test_20110125';
% -----------------------------------------------------

% load parms structure from dSPM
load(sprintf('%s/matfiles/%s_parms.mat',rootoutdir,prefix));

% Load avg_data and get time vector
load(sprintf('%s/matfiles/%s_avg_data.mat',rootoutdir,prefix)); % avg_data
t = avg_data.averages(1).time;

cwd = pwd; cd(rootoutdir)

%% Anatomical data (FreeSurfer recon)

% % Load subject recon (created by FreeSurfer)
lh_surf = fs_load_subj(parms.subjname,'lh','white',0,parms.subjdir);
rh_surf = fs_load_subj(parms.subjname,'rh','white',0,parms.subjdir);

% T2: Coordinate transformation for the brain  [rotation+translation]
add       = [0 0 0]';%[0.14 0.13 0.16]';   % translations  (addX, addY, addZ)
rot       = [0 0 0]';%[0 -pi/40 pi/2]; % pi/10 % rotations (rotX, rotY, rotZ)
Rx        = [1 0 0 0; 0 cos(rot(1)) -sin(rot(1)) 0; 0 sin(rot(1)) cos(rot(1)) 0; 0 0 0 1];
Ry        = [cos(rot(2)) 0 sin(rot(2)) 0; 0 1 0 0; -sin(rot(2)) 0 cos(rot(2)) 0; 0 0 0 1];
Rz        = [cos(rot(3)) -sin(rot(3)) 0 0; sin(rot(3)) cos(rot(3)) 0 0; 0 0 1 0; 0 0 0 1];
T2        = Rz*Ry*Rx*eye(4);  % rotation
T2(1:3,4) = T2(1:3,4) + add;  % translation

if ~isfield(lh_surf,'vertices') && isfield(lh_surf,'coords')
  lh_surf.vertices=lh_surf.coords; 
  rh_surf.vertices=rh_surf.coords;
end

% Brain
lh_verts = lh_surf.vertices'; LL = length(lh_surf.vertices);  % [3 x Lnverts]
rh_verts = rh_surf.vertices'; LR = length(rh_surf.vertices);  % [3 x Rnverts]
brain    = [lh_verts rh_verts];                               % [3 x nverts]  where nverts = Lnverts + Rnverts
brain    = brain;
brain    = [brain; ones(1,length(brain))];                    % [4 x nverts]
brain    = T2*brain;  % Apply brain coord transformation

% FreeSurfer left hemi brain notes
  % [lh_surf.coords] = 151179 x 3
  % [lh_surf.nbrs]   = 151179 x 27
  % [lh_surf.nverts] = 151179
% FreeSurfer right hemi brain notes
  % [rh_surf.coords] = 150930 x 3
  % [rh_surf.nbrs]   = 150930 x 27
  % [rh_surf.nverts] = 150930

%% Functional data (dSPM) from MGH files
conds = parms.conditions;
cond  = conds(c);

% MGH files with dSPM results
mghfile_lh_ico7 = sprintf('mghfiles/%s_cond%02i-spsm10-sm10-ico7-sm3-lh.mgh',prefix,cond);
mghfile_lh      = sprintf('mghfiles/%s_cond%02i-spsm10-sm10-lh.mgh',prefix,cond);
mghfile_rh_ico7 = sprintf('mghfiles/%s_cond%02i-spsm10-sm10-ico7-sm3-rh.mgh',prefix,cond);
mghfile_rh      = sprintf('mghfiles/%s_cond%02i-spsm10-sm10-rh.mgh',prefix,cond);

% Load source activity in native space 
[lvol,lM] = fs_load_mgh(mghfile_lh);      % [vol] = 151179 x 1 x 1 x 331
[rvol,lM] = fs_load_mgh(mghfile_rh);      % [vol] = 150930 x 1 x 1 x 331

% xyz coordinates from FreeSurfer brain after transformation
lx = brain(1,1:LL);
ly = brain(2,1:LL);
lz = brain(3,1:LL);
rx = brain(1,LL+1:LL+LR);
ry = brain(2,LL+1:LL+LR);
rz = brain(3,LL+1:LL+LR);

% PARAMETERS ----------------------------------------------------
% Plot parameters
plotfun = @trisurf; % @trisurf or @trimesh
eAlpha  = 0;        % edge transparency (set to 1 for trimesh else <=1)
fAlpha  = 1;        % face transparency (set to 1 for trisurf else <=1)

labels  = {'left' 'top'   'front' 'right' 'bottom'  'back'};
views   = {[90 0] [90 90] [180 0] [-90 0] [-90 -90] [0 0]};
nrow    = 2;                          % number of subplot rows
ncol    = 3;                          % number of subplot columns

% Movie parameters
frames  = 1:length(t);                                                % indices of frames to show in movie
fps     = avg_data.sfreq / 10; %/ 5;%  10;                            % frames per second in the movie
avi     = sprintf('%s/%s_cond%02i_movie2.avi',rootoutdir,prefix,cond); % filename for saving the movie
% ----------------------------------------------------

% set figure properties and check the number of subplots
fpos    = [50 50 500*ncol 400*nrow]; % screen position of the figure
clims   = [min(min(lvol(:)),min(rvol(:))) max(max(lvol(:)),max(rvol(:)))];
xlims   = [min([lx rx]) max([lx rx])];
ylims   = [min([ly ry]) max([ly ry])];
zlims   = [min([lz rz]) max([lz rz])];
nplots  = length(views);             % number of subplots for each frame
nframes = length(frames);            % number of movie frames
if nrow*ncol < length(views)
  ncol  = ceil(nplots/nrow); 
end

tic
figure; set(gcf,'position',fpos);
aviobj  = avifile(avi,'fps',fps);    % create avi object
for k   = 1:nframes
  fprintf('Drawing frame %g of %g\n',k,nframes); clf
  frame = frames(k);
  lC    = lvol(:,1,1,frame)'; % source activity of left hemi for this frame
  rC    = rvol(:,1,1,frame)'; % source activity of right hemi for this frame
  for j = 1:nplots
    subplot(nrow,ncol,j)
    % plot FreeSurfer elements of left brain with colors = source activity
    p1 = plotfun(lh_surf.faces,lx,ly,lz,lC,'EdgeLighting','phong','FaceLighting','phong',...
                 'EdgeAlpha',eAlpha,'FaceAlpha',fAlpha); %,'FaceColor',[238 223 204]./256);
    % hold plot for adding right hemisphere
    hold on
    % plot FreeSurfer elements of right brain with colors = source activity
    p2 = plotfun(rh_surf.faces,rx,ry,rz,rC,'EdgeLighting','phong','FaceLighting','phong',...
                 'EdgeAlpha',eAlpha,'FaceAlpha',fAlpha); %,'FaceColor',[238 223 204]./256);  
    % adjust the lighting
    light('Position',[0 -1 .9]); lighting phong; camlight %headlight
    % adjust axes
    axis equal; box off; grid off; axis off;
    set(gca,'xlim',xlims,'ylim',ylims,'zlim',zlims,'clim',clims);
    view(views{j}); title(labels{j});
  end
  % add annotations to figure
  annotation('textbox','Position',[.3 .9 .4 .1]  ,'String',sprintf('%gms',round(1000*t(frame))),'VerticalAlignment','top','HorizontalAlignment','center','FontSize',14,'Color','k','FitHeightToText','on','LineStyle','none');
  annotation('textbox','Position',[.01 .01 .9 .1],'String',sprintf('subject: %s (condition %02i)\n(%s)',strrep(parms.subjname,'_','\_'),cond,strrep(parms.subjdir,'_','\_')),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10,'Color','k','FitHeightToText','on','LineStyle','none');
  % add frame to movie
%   drawnow
  aviframe = getframe(gcf);
  aviobj   = addframe(aviobj,aviframe);
  toc
end
close
aviobj = close(aviobj);

% movie compression
try
  % check that mencoder is in your linux path
  mencoderpath = '/home/jsherfey/projects/movies/mencoder';
  avi2         = '';
  [s,r]=unix('echo $path');
  if ~any(findstr(mencoderpath,r))
    error('Add mencoder to your linux path: %s',mencoderpath);
    return
  end
  % compress avi using mencoder
  avi2  = sprintf('%s_compressed.avi',avi(1:end-4));
  cmd   = sprintf('mencoder %s -oac lavc -ovc lavc -lavcopts acodec=mp3,1bitrate=128,vcodec=mpeg4,vbitrate=800,vhq,vm4v -o %s',avi,avi2);
  [s,r] = unix(cmd); if s, disp(r); cd(cwd); return; end
%   cmd   = sprintf('rm %s',avi);
%   [s,r] = unix(cmd); if s, disp(r); cd(cwd); return; end
%   avi = [avi(1:end-4) '_lowrate'];
end

% watch uncompressed movie
fps = 10;                 % frames per second during playback
figure('position',fpos);
mov = aviread(avi);
movie(mov,1,fps);

fprintf('New uncompressed avi movie: %s\n',avi);
if avi2, fprintf('New compressed avi movie: %s\n',avi2); end

% -----------------
% ICO 7
% -----------------

% % Load source data in ico7
% [vol7, M7, mr_parms7, volsz7] = fs_load_mgh(mghfile_lh_ico7);
% % [vol] = 163842 x 1 x 1 x 331
% % [M]   = 4x4 transformation matrix

fprintf('Time to generate movie (%g frames): %gmin\n',nframes,toc(tstart)/60);
cd(cwd); return

%% Functional data (dSPM) from MAT files (slow bc requires source interpolation)

% Load avg_data and get time vector
load(sprintf('matfiles/%s_avg_data.mat',prefix)); % avg_data
t = avg_data.averages(1).time;

% MEG source (dipole) locations in MRI space
lh_ix    = find(parms.lh_dec_dips==1); % 3434x1
rh_ix    = find(parms.rh_dec_dips==1); % 3448x1    : sum = 6882
grid_mri = cat(2,parms.lh_dip_info(1:3,lh_ix),parms.rh_dip_info(1:3,rh_ix));

x = grid_mri(1,:); 
y = grid_mri(2,:); 
z = grid_mri(3,:);
X = [x' y' z'];

% load dSPM source data
load(sprintf('matfiles/%s_results_cond%02i.mat',prefix,cond));
% [S_xyz] = 331x20646, [F] = [S] = 331x6882

ds     = 2;         % controls mesh resolution during source interpolation
nplots = length(t);%331; 10       % number of movie frames
toilim = [min(t) max(t)];%[.2 .8];
tplots = toilim(1):diff(toilim)/(nplots-1):toilim(2);

% points in 3D space whose surface value will be fitted
[x0,y0,z0] = meshgrid([min(x):ds:max(x)],[min(y):ds:max(y)],[min(z):ds:max(z)]);
XI  = [x0(:) y0(:) z0(:)];

% source interpolation
tic
YI    = zeros([size(x0) nplots]);
T     = zeros(1,nplots);
for k = 1:nplots
  tplot = tplots(k);
  fprintf('Processing t = %gsec (%g of %g)\n',tplot,k,nplots);
  tix  = nearest(t,tplot);
  T(k) = t(tix);
  Y    = S(tix,:)';
  tmp  = griddatan(X,Y,XI);
  YI(:,:,:,k) = reshape(tmp, size(x0));
  clear tmp; toc
end
% YI(:,:,:,k) = F(XI) @ t = T(k)

save('test_dSPM_plot4_F.mat','T','YI','x0','y0','z0');

% source plot and create slideshow/movie
fps    = avg_data.sfreq / 10; %/ 5;%  10;                               % frames per second in the movie
avi    = sprintf('%s_movie4_F.avi',prefix);  % filename for saving the movie
aviobj = avifile(avi,'fps',fps);          % create avi object

aFaceAlpha = .5; % MRI transparency (anatomical)
fFaceAlpha = .5; % dSPM transparency (functional)
views = {[-80 10],[0 90],[0 0]};
figure
set(gcf,'position',[80 400 1500 425]);
tic
for k = 1:nplots
  fprintf('Plotting figure %g of %g\n',k,nplots);
  clf
  for j = 1:length(views)
    subplot(1,length(views),j)
    % plot FreeSurfer elements of left brain
    p = trimesh(lh_surf.faces,brain(1,1:LL),brain(2,1:LL),brain(3,1:LL),'EdgeColor',[238 223 204]./256,'EdgeAlpha',1,...
        'EdgeLighting','phong','FaceColor',[238 223 204]./256,'FaceAlpha',aFaceAlpha,'FaceLighting','phong'); hold on;
    % plot FreeSurfer elements of right brain
    p = trimesh(rh_surf.faces,brain(1,LL+1:end),brain(2,LL+1:end),brain(3,LL+1:end),'EdgeColor',[238 223 204]./256,'EdgeAlpha',1,...
                                       'EdgeLighting','phong','FaceColor',[238 223 204]./256,'FaceAlpha',aFaceAlpha,'FaceLighting','phong'); hold on;
    % adjust the lighting
    light('Position',[0 -1 .9]); axis equal; %set(gca,'xlim',[-.01 .255],'ylim',[-.01 .255],'zlim',[0 .27]); 
    % overlay dSPM on the brain
    p2 = patch(isosurface(x0,y0,z0,YI(:,:,:,k)));%,0.02));
    isonormals(x0,y0,z0,YI(:,:,:,k),p2);
    set(p2,'FaceAlpha',fFaceAlpha,'EdgeColor','none','FaceColor','blue');
    box off; grid off; axis off; % camlight, lighting phong      
    view(views{j})
    title(sprintf('t = %gsec',T(k)));
    if j~=3, camlight headlight; end
  end
  drawnow
  Fr      = getframe(gcf);
  aviobj = addframe(aviobj,Fr);
  toc
end
aviobj = close(aviobj);

% watch movie
mov = aviread(avi);
fps = 10;             % playback frames per second
movie(mov,1,fps);

%% retired

% %%
% 
% figure;
% aFaceAlpha=.2;
% fFaceAlpha=.8;
% % plot FreeSurfer elements of left brain
% p = trimesh(lh_surf.faces,brain(1,1:LL),brain(2,1:LL),brain(3,1:LL),'EdgeColor',[238 223 204]./256,'EdgeAlpha',1,...
%     'EdgeLighting','phong','FaceColor',[238 223 204]./256,'FaceAlpha',aFaceAlpha,'FaceLighting','phong'); hold on;
% % plot FreeSurfer elements of right brain
% p = trimesh(rh_surf.faces,brain(1,LL+1:end),brain(2,LL+1:end),brain(3,LL+1:end),'EdgeColor',[238 223 204]./256,'EdgeAlpha',1,...
%                                    'EdgeLighting','phong','FaceColor',[238 223 204]./256,'FaceAlpha',aFaceAlpha,'FaceLighting','phong'); hold on;
% % adjust the lighting
% light('Position',[0 -1 .9]); axis equal; %set(gca,'xlim',[-.01 .255],'ylim',[-.01 .255],'zlim',[0 .27]); 
% 
% hold on;%figure
% p2 = patch(isosurface(x0,y0,z0,YI));%,0.02));
% isonormals(x0,y0,z0,YI,p2);
% set(p2,'FaceAlpha',fFaceAlpha,'EdgeColor','blue','FaceColor','blue');
% box off; grid off; axis off; view(3); % camlight, lighting phong      
% 

% 
% 
% % % Get MEG and source locations
% [meg_locs,source_locs] = ts_plot3d_dSPM_sensors_and_sources_jason(prefix);
% % lh_ix    = find(parms.lh_dec_dips==1); % 3434x1
% % rh_ix    = find(parms.rh_dec_dips==1); % 3448x1    : sum = 6882
% % grid_mri = cat(2,parms.lh_dip_info(1:3,lh_ix),parms.rh_dip_info(1:3,rh_ix));
% % 
% % % load dSPM results
% % load(sprintf('matfiles/%s_results_cond%02i.mat',prefix,cond));
% % % [S_xyz] = 331x20646, [F] = [S] = 331x6882
% 
% %%
% 
% % test plot
% % x=grid_mri(1,:); y=grid_mri(2,:); z=grid_mri(3,:);
% % figure; h=plot3(x,y,z,'.'); set(h,'MarkerSize',4,'LineWidth',.4);
% 
% % % load BEM
% % load(sprintf('matfiles/%s_bem.mat',prefix));      % [L_bem] = [U_bem] = 3318x3318
% % 
% % % load forward matrix
% % load(sprintf('matfiles/%s_forward.mat',prefix));  % [G_xyz] = 204x20646
% % 
% % % load inverse matrix
% % load(sprintf('matfiles/%s_inverse.mat',prefix));  % [G_xyz] = 204x20646
% 
% % load dSPM results
% load(sprintf('matfiles/%s_results_cond%02i.mat',prefix,cond));
% % [S_xyz] = 331x20646, [F] = [S] = 331x6882
% 
% % test plot
% tplot = .3; tix = nearest(t,tplot);
% 
% x = grid_mri(1,:); 
% y = grid_mri(2,:); 
% z = grid_mri(3,:);
% C = S(tix,:); % S, F
% 
% 
% % figure; mesh(x,y,z,C);
% % figure; mesh(X,Y,Z,
% 
% figure; h=plot3(x,y,z,'.'); set(h,'MarkerSize',4,'LineWidth',.4);
% 
% 
% % lh_dip_file         = 'bem/lh_white.dip';
% % rh_dip_file         = 'bem/rh_white.dip';
% % lh_dec_file         = 'bem/lh_white_7.dec';
% % rh_dec_file         = 'bem/rh_white_7.dec';
% % bem_surf_files      = {'bem/inner_skull.tri'};
% % 
% % lh_dip_info = ts_read_dip_file(sprintf('%s/%s/%s',parms.subjdir,parms.subjname,lh_dip_file));
% 
% 
% 
% %%
% % tplot = .3; 
% % tix   = nearest(t,tplot);
% % C     = S(tix,:); % S, F
% % 
% % hold on; h=plot3(x,y,z,'.'); set(h,'MarkerSize',5,'LineWidth',.4);
% % 
% % tic
% % X   = [x' y' z'];
% % Y   = C';
% % ds  = 2;
% % [x0,y0,z0] = meshgrid([min(x):ds:max(x)],[min(y):ds:max(y)],[min(z):ds:max(z)]);
% % XI  = [x0(:) y0(:) z0(:)];
% % 
% % YI = griddatan(X,Y,XI);
% % YI = reshape(YI, size(x0));
% % 
% % hold on;%figure
% % p2 = patch(isosurface(x0,y0,z0,YI,0.02));
% % isonormals(x0,y0,z0,YI,p);
% % set(p2,'FaceAlpha',.2,'EdgeColor','none','FaceColor','blue');
% % view(3), axis off, camlight, lighting phong      
% toc