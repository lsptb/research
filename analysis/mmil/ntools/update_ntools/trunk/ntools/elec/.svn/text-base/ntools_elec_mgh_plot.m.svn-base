function ntools_elec_mgh_plot(varargin)

% ntools_elec_mgh_plot('loc text','surface mat','mgh overlay')
% surface mat: only one hemi is allowed.

%% colormap
cm =[...
    0.7094    0.7008    0.7000;
    0.7187    0.7617    0.5731;
    0.7281    0.8187    0.4692;
    0.7375    0.8703    0.3842;
    0.7469    0.9149    0.3145;
    0.7562    0.9512    0.2575;
    0.7656    0.9780    0.2108;
    0.7750    0.9945    0.1726;
    0.7844    1.0000    0.1413;
    0.7937    0.9945    0.1157;
    0.8031    0.9780    0.0947;
    0.8125    0.9512    0.0776;
    0.8219    0.9149    0.0635;
    0.8312    0.8703    0.0520;
    0.8406    0.8187    0.0426;
    0.8500    0.7617    0.0349;
    0.8594    0.7008         0;
    0.8687    0.6376         0;
    0.8781    0.5738         0;
    0.8875    0.5106         0;
    0.8969    0.4493         0;
    0.9062    0.3911         0;
    0.9156    0.3366         0;
    0.9250    0.2865         0;
    0.9344    0.2412         0;
    0.9437    0.2008         0;
    0.9531    0.1653         0;
    0.9625    0.1346         0;
    0.9719    0.1084         0;
    0.9812    0.0863         0;
    0.9906    0.0680         0;
    1.0000    0.0529         0;
];

%% Get the elec info
if nargin==0
    [FileName,PathName] = uigetfile('*.txt','Select the electrodes text file','/home/ccarlson/loc/');
    [surfname, surfpath] = uigetfile('*.mat','Select the patient brain surf',PathName);
    [mghname mghpath] = uigetfile('*.mgh','Select mgh overlay file','/home/nyuproj/subjects/');
    surf = strcat(surfpath,surfname);
    mgh = strcat(mghpath,mghname);
elseif nargin==3
    aa = findstr(varargin{1},'/');
    FileName = varargin{1}(aa(end)+1:end);
    PathName = varargin{1}(1:aa(end));
    surf = varargin{2}; 
    mgh = varargin{3};
    bb = findstr(varargin{3},'/');
    mghname = varargin{3}(bb(end)+1:end);

else
    disp('not correct input');
    return;
end

fid = fopen([PathName, FileName]);
elec_all = textscan(fid,'%s %f %f %f %s');
elec_cell = [elec_all{1},num2cell(elec_all{2}),num2cell(elec_all{3}),num2cell(elec_all{4})];

vol = fs_load_mgh(mgh);
% cope = mgh(findstr(mgh,'cope'):findstr(mgh,'.feat')-1);
% zthres = mgh(findstr(mgh,'z'):findstr(mgh,'.gfeat')-1);
% mghname = [zthres '-' cope '-' mghname(1:end-4)];
mghname = mghname(1:end-4);

%% show only positive
vol(find(vol<0)) = 0;

%% Get the filename info
b = findstr(FileName,'_');
Pname = FileName(1:b(1)-1);

if ~isempty(findstr(upper(FileName),'T1'))
    space = '_T1_';
elseif ~isempty(findstr(upper(FileName),'MNI'))
    space = '_MNI_';
else
    space = '_';
end

sph = regexp(surf,'[r,l]h','match');
sph = char(sph{:});

%% get lh(rh) electrodes only

g = strmatch('G',upper(elec_cell(:,1)));
d = strmatch('D',upper(elec_cell(:,1)));
if ~isempty(g) && ~isempty(d)
    elec_grid = elec_cell(g,:);
    elec_depth = elec_cell(d,:);
    elec_cell([g;d],:) = [];    
elseif isempty(d)
    elec_depth = [];
    elec_grid = elec_cell(g,:);
    elec_cell(g,:) = [];
elseif isempty(g)
    elec_grid = [];
    elec_depth = elec_cell(d,:);
    elec_cell(d,:) = [];
end


%% Plot the elecs
colrange = input('Input the overlay threshold range, e.g. [1.8 5] : ');
plt = menu('What part do you want to plot?','Grid only', 'Strip only','Depth Only','Both Grid and Strip');
labelshow = menu('Do you want to show the label?','Yes','No');
genimg = menu('Do you want to save the images?','Yes', 'No');
% colrange = [2 5];
% plt = 4;
% labelshow = 2;
% genimg = 1;

surf_brain = load(surf);


if plt==1 && ~isempty(elec_grid)
    showpart = 'G';
    nyu_plot(surf_brain,sph,cell2mat(elec_grid(:,2:4)),char(elec_grid(:,1)),vol,'r',labelshow);
    colormap(cm)
    set(gca,'CLim',colrange)
    colorbar
    set(colorbar,'YColor','w')
elseif plt==2 && ~isempty(elec_cell)
    showpart = 'S';
    nyu_plot(surf_brain,sph,cell2mat(elec_cell(:,2:4)),char(elec_cell(:,1)),vol,'b',labelshow);
    colormap(cm)
    set(gca,'CLim',colrange)
    colorbar
    set(colorbar,'YColor','w')
elseif plt==3 && ~isempty(elec_depth)
    showpart = 'D';
    nyu_plot(surf_brain,sph,cell2mat(elec_depth(:,2:4)),char(elec_depth(:,1)),vol,'g',labelshow,0.3,6);
    colormap(cm)
    set(gca,'CLim',colrange)
    colorbar
    set(colorbar,'YColor','w')
elseif plt==4 && ~isempty(elec_grid) && ~isempty(elec_cell)
    showpart = 'GS';
    elec = cell2mat(elec_cell(:,2:4));
    elec_name = char(elec_cell(:,1));
    nyu_plot(surf_brain,sph,cell2mat(elec_grid(:,2:4)),char(elec_grid(:,1)),vol,'r',labelshow); hold on;
    colormap(cm)
    set(gca,'CLim',colrange)
    colorbar
    set(colorbar,'YColor','w')
    for i=1:size(elec,1)
        plot3(elec(i,1),elec(i,2),elec(i,3),'bo','MarkerFaceColor','b','MarkerSize',11.3);
        if labelshow==1
            text('Position',[elec(i,1) elec(i,2) elec(i,3)],'String',elec_name(i,:),'Color','w');
        end
    end
    hold off;   
else
    disp('sorry, the electrodes you choose to show are not on the surface you loaded');
    return;
end

%% save images

if genimg==1
    if ~exist([PathName 'images/'],'dir')
        mkdir([PathName 'images/']);
    end
    ext = 'jpg';  
    if strcmp(sph,'lh')
        view(270, 0);
        saveas(gcf,[PathName,'images/',Pname,space,showpart,'_lateral_',mghname,'.',ext]);
        view(90,0);
        saveas(gcf,[PathName,'images/',Pname,space,showpart,'_mesial_',mghname,'.',ext]);
        
    elseif strcmp(sph,'rh')
        view(270, 0);
        saveas(gcf,[PathName,'images/',Pname,space,showpart,'_mesial_',mghname,'.',ext]);
        view(90,0);
        saveas(gcf,[PathName,'images/',Pname,space,showpart,'_lateral_',mghname,'.',ext]);
    end
    view(0,0);
    saveas(gcf,[PathName,'images/',Pname,space,showpart,'_posterior_',mghname,'.',ext]);

    view(180,0);
    saveas(gcf,[PathName,'images/',Pname,space,showpart,'_frontal_',mghname,'.',ext]);

    view(90,90);
    saveas(gcf,[PathName,'images/',Pname,space,showpart,'_dorsal_',mghname,'.',ext]);

    view(90,-90);
    set(light,'Position',[1 0 -1]);
    saveas(gcf,[PathName,'images/',Pname,space,showpart,'_ventral_',mghname,'.',ext]);
else 
    return;
end

end
%% nyu_plot
function nyu_plot(surf_brain,sph,elec,elecname,surfcol,color,label,alpha,marksize)

if ~exist('color','var')
    color = 'w';
end
if ~exist('label','var')
    label = 2;
end
if ~exist('alpha','var')
    alpha = 1;
end
if ~exist('marksize','var')
    marksize = 11.3;
end

figure;

if isfield(surf_brain,'coords')==0
    sub.vert = surf_brain.surf_brain.coords;
    sub.tri = surf_brain.surf_brain.faces;
else
    sub.vert = surf_brain.coords;
    sub.tri = surf_brain.faces;
end

c=zeros([length(sub.vert),1])+surfcol;

% c=(c/max(c));
trisurf(sub.tri, sub.vert(:, 1), sub.vert(:, 2), sub.vert(:, 3),...
    'FaceVertexCData', c,'FaceColor', 'interp','FaceAlpha',alpha);

shading interp;
lighting gouraud;
material dull;
light;
axis off
hold on;
for i=1:length(elec)
    plot3(elec(i,1),elec(i,2),elec(i,3),[color 'o'],'MarkerFaceColor',color,'MarkerSize',marksize);
    if label==1
        text('Position',[elec(i,1) elec(i,2) elec(i,3)],'String',elecname(i,:),'Color','w');
    end
end
set(light,'Position',[-1 0 1]); 
    if strcmp(sph,'lh')
        view(270, 0);      
    elseif strcmp(sph,'rh')
        view(90,0);        
    elseif strcmp(sph,'both')
        view(90,90);
    end
set(gcf, 'color','black','InvertHardCopy', 'off');
axis tight;
axis equal;
end