addpath(genpath('/home/jsherfey/svn/dev/onestream'));
% -----------------------------------------------------
% Load dSPM parms:
% -----------------------------------------------------
% load matfiles/test_dSPM_WordNPNW_grad_bem_parms
% prefix     = parms.prefix;
% rootoutdir = parms.rootoutdir;
% -----------------------------------------------------
% or
% -----------------------------------------------------
prefix     = 'test_dSPM_WordNPNW_grad_bem';
rootoutdir = '/home/jsherfey/svn/dev/onestream/test_20110125';
load(sprintf('%s/matfiles/%s_parms.mat',rootoutdir,prefix));
% -----------------------------------------------------
% -----------------------------------------------------

cd(rootoutdir)
conds = parms.conditions;
c     = 1;
cond  = conds(c);

% load avg_data
load(sprintf('matfiles/%s_avg_data.mat',prefix)); % avg_data
t = avg_data.averages(1).time;

% Get meg and source locations
[meg_locs,source_locs] = ts_plot3d_dSPM_sensors_and_sources_jason(prefix);
lh_ix    = find(parms.lh_dec_dips==1); % 3434x1
rh_ix    = find(parms.rh_dec_dips==1); % 3448x1    : sum = 6882
grid_mri = cat(2,parms.lh_dip_info(1:3,lh_ix),parms.rh_dip_info(1:3,rh_ix));

% test plot
x=grid_mri(1,:); y=grid_mri(2,:); z=grid_mri(3,:);
figure; h=plot3(x,y,z,'.'); set(h,'MarkerSize',4,'LineWidth',.4);

% % load BEM
% load(sprintf('matfiles/%s_bem.mat',prefix));      % [L_bem] = [U_bem] = 3318x3318
% 
% % load forward matrix
% load(sprintf('matfiles/%s_forward.mat',prefix));  % [G_xyz] = 204x20646
% 
% % load inverse matrix
% load(sprintf('matfiles/%s_inverse.mat',prefix));  % [G_xyz] = 204x20646

% load dSPM results
load(sprintf('matfiles/%s_results_cond%02i.mat',prefix,cond));
% [S_xyz] = 331x20646, [F] = [S] = 331x6882

% test plot
tplot = .3; tix = nearest(t,tplot);

x = grid_mri(1,:); 
y = grid_mri(2,:); 
z = grid_mri(3,:);
C = S(tix,:); % S, F


% figure; mesh(x,y,z,C);
% figure; mesh(X,Y,Z,

figure; h=plot3(x,y,z,'.'); set(h,'MarkerSize',4,'LineWidth',.4);


% lh_dip_file         = 'bem/lh_white.dip';
% rh_dip_file         = 'bem/rh_white.dip';
% lh_dec_file         = 'bem/lh_white_7.dec';
% rh_dec_file         = 'bem/rh_white_7.dec';
% bem_surf_files      = {'bem/inner_skull.tri'};
% 
% lh_dip_info = ts_read_dip_file(sprintf('%s/%s/%s',parms.subjdir,parms.subjname,lh_dip_file));

%% recon (anatomical)

% % Load subject recon (created by FreeSurfer)
lh_surf = fs_load_subj(parms.subjname,'lh','white',0,parms.subjdir);
rh_surf = fs_load_subj(parms.subjname,'rh','white',0,parms.subjdir);

% T2: Coordinate transformation for the brain  [rotation+translation]
add       = [0 0 0]';%[0.14 0.13 0.16]';   % translations  (addX, addY, addZ)
rot       = [0 0 0]';%[0 -pi/40 pi/2]; % pi/10
Rx        = [1 0 0 0; 0 cos(rot(1)) -sin(rot(1)) 0; 0 sin(rot(1)) cos(rot(1)) 0; 0 0 0 1];
Ry        = [cos(rot(2)) 0 sin(rot(2)) 0; 0 1 0 0; -sin(rot(2)) 0 cos(rot(2)) 0; 0 0 0 1];
Rz        = [cos(rot(3)) -sin(rot(3)) 0 0; sin(rot(3)) cos(rot(3)) 0 0; 0 0 1 0; 0 0 0 1];
T2        = Rz*Ry*Rx*eye(4);  % rotation
T2(1:3,4) = T2(1:3,4) + add;  % translation

if ~isfield(lh_surf,'vertices') && isfield(lh_surf,'coords')
  lh_surf.vertices=lh_surf.coords; 
  rh_surf.vertices=rh_surf.coords;
end

lh_verts = lh_surf.vertices'; LL = length(lh_surf.vertices);  % [3 x Lnverts]
rh_verts = rh_surf.vertices'; LR = length(rh_surf.vertices);  % [3 x Rnverts]
brain    = [lh_verts rh_verts];                               % [3 x nverts]  where nverts = Lnverts + Rnverts
brain    = brain;
brain    = [brain; ones(1,length(brain))];                    % [4 x nverts]
brain    = T2*brain;  % Apply brain coord transformation

% ealph = .0;
% FaceAlpha = .2;
% % front:-90,0; above: 0, 90
% figure;
% % plot FreeSurfer elements of left brain
% p = trimesh(lh_surf.faces,brain(1,1:LL),brain(2,1:LL),brain(3,1:LL),'EdgeColor',[238 223 204]./256,'EdgeAlpha',1,...
%     'EdgeLighting','phong','FaceColor',[238 223 204]./256,'FaceAlpha',FaceAlpha,'FaceLighting','phong'); hold on;
% % plot FreeSurfer elements of right brain
% p = trimesh(rh_surf.faces,brain(1,LL+1:end),brain(2,LL+1:end),brain(3,LL+1:end),'EdgeColor',[238 223 204]./256,'EdgeAlpha',1,...
%                                    'EdgeLighting','phong','FaceColor',[238 223 204]./256,'FaceAlpha',FaceAlpha,'FaceLighting','phong'); hold on;
% % adjust the lighting
% light('Position',[0 -1 .9]); axis equal; view(0,0); %set(gca,'xlim',[-.01 .255],'ylim',[-.01 .255],'zlim',[0 .27]); 
% box off; grid off; axis off;

%% dSPM (functional)

% dipole locations
lh_ix    = find(parms.lh_dec_dips==1); % 3434x1
rh_ix    = find(parms.rh_dec_dips==1); % 3448x1    : sum = 6882
grid_mri = cat(2,parms.lh_dip_info(1:3,lh_ix),parms.rh_dip_info(1:3,rh_ix));

x = grid_mri(1,:); 
y = grid_mri(2,:); 
z = grid_mri(3,:);
X = [x' y' z'];

[x0,y0,z0] = meshgrid([min(x):ds:max(x)],[min(y):ds:max(y)],[min(z):ds:max(z)]);
XI  = [x0(:) y0(:) z0(:)];

% load dSPM
load(sprintf('matfiles/%s_results_cond%02i.mat',prefix,cond));
% [S_xyz] = 331x20646, [F] = [S] = 331x6882

ds     = 2;         % controls mesh resolution during source interpolation
nplots = 10;
toilim = [-.1 .8];
tplots = toilim(1):diff(toilim)/(nplots-1):toilim(2);

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

save('test_dSPM_plot1.mat','T','YI','x0','y0','z0');

% plot and create slideshow/movie
fps    = 1;                       % frames per second in the movie
avi    = 'test_dSPM_movie.avi';   % filename for saving the movie
aviobj = avifile(avi,'fps',fps);  % create avi object

aFaceAlpha = .2; % MRI transparency (anatomical)
fFaceAlpha = .8; % dSPM transparency (functional)

figure
for k = 1:nplots
  
  % plot FreeSurfer elements of left brain
  p = trimesh(lh_surf.faces,brain(1,1:LL),brain(2,1:LL),brain(3,1:LL),'EdgeColor',[238 223 204]./256,'EdgeAlpha',1,...
      'EdgeLighting','phong','FaceColor',[238 223 204]./256,'FaceAlpha',aFaceAlpha,'FaceLighting','phong'); hold on;
  % plot FreeSurfer elements of right brain
  p = trimesh(rh_surf.faces,brain(1,LL+1:end),brain(2,LL+1:end),brain(3,LL+1:end),'EdgeColor',[238 223 204]./256,'EdgeAlpha',1,...
                                     'EdgeLighting','phong','FaceColor',[238 223 204]./256,'FaceAlpha',aFaceAlpha,'FaceLighting','phong'); hold on;
  % adjust the lighting
  light('Position',[0 -1 .9]); axis equal; %set(gca,'xlim',[-.01 .255],'ylim',[-.01 .255],'zlim',[0 .27]); 

  hold on;%figure
  p2 = patch(isosurface(x0,y0,z0,YI));%,0.02));
  isonormals(x0,y0,z0,YI,p2);
  set(p2,'FaceAlpha',fFaceAlpha,'EdgeColor','blue','FaceColor','blue');
  box off; grid off; axis off; view(3); % camlight, lighting phong      

    F      = getframe(gcf);
    aviobj = addframe(aviobj,F);


aviobj = close(aviobj);

% watch movie
mov = aviread(avi);
fps = 10;             % frames per second (dt=100ms ==> fps=1/dt=10)
movie(mov,1,fps);


%%

figure;
aFaceAlpha=.2;
fFaceAlpha=.8;
% plot FreeSurfer elements of left brain
p = trimesh(lh_surf.faces,brain(1,1:LL),brain(2,1:LL),brain(3,1:LL),'EdgeColor',[238 223 204]./256,'EdgeAlpha',1,...
    'EdgeLighting','phong','FaceColor',[238 223 204]./256,'FaceAlpha',aFaceAlpha,'FaceLighting','phong'); hold on;
% plot FreeSurfer elements of right brain
p = trimesh(rh_surf.faces,brain(1,LL+1:end),brain(2,LL+1:end),brain(3,LL+1:end),'EdgeColor',[238 223 204]./256,'EdgeAlpha',1,...
                                   'EdgeLighting','phong','FaceColor',[238 223 204]./256,'FaceAlpha',aFaceAlpha,'FaceLighting','phong'); hold on;
% adjust the lighting
light('Position',[0 -1 .9]); axis equal; %set(gca,'xlim',[-.01 .255],'ylim',[-.01 .255],'zlim',[0 .27]); 

hold on;%figure
p2 = patch(isosurface(x0,y0,z0,YI));%,0.02));
isonormals(x0,y0,z0,YI,p2);
set(p2,'FaceAlpha',fFaceAlpha,'EdgeColor','blue','FaceColor','blue');
box off; grid off; axis off; view(3); % camlight, lighting phong      



%%
% tplot = .3; 
% tix   = nearest(t,tplot);
% C     = S(tix,:); % S, F
% 
% hold on; h=plot3(x,y,z,'.'); set(h,'MarkerSize',5,'LineWidth',.4);
% 
% tic
% X   = [x' y' z'];
% Y   = C';
% ds  = 2;
% [x0,y0,z0] = meshgrid([min(x):ds:max(x)],[min(y):ds:max(y)],[min(z):ds:max(z)]);
% XI  = [x0(:) y0(:) z0(:)];
% 
% YI = griddatan(X,Y,XI);
% YI = reshape(YI, size(x0));
% 
% hold on;%figure
% p2 = patch(isosurface(x0,y0,z0,YI,0.02));
% isonormals(x0,y0,z0,YI,p);
% set(p2,'FaceAlpha',.2,'EdgeColor','none','FaceColor','blue');
% view(3), axis off, camlight, lighting phong      
% toc