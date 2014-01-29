function ntools_elec_view

% a stand-alone program that shows ieeg electrode coordinates on the pial
% surface 
% requires:
% elec_text: text file with xyz electrode coordinates
% pial_mat: matlab structure with pial surface
%
% Usage: run nyu_ieeg_elec_plot in command window
% the gui is prompted out for file selection
%
% also see: http://ieeg.pbworks.com/Viewing-Electrode-Locations
%
% written by  Hugh Wang, wangxiuyuan85@gmail.com, May 13, 2009

% modified on May 14, 2009 by Hugh
% make judgement on the input file type and not sensitive to the order of 
% input variable. add the option to show the electrodes' labels or not.

% modified on July 22nd, 2009 by Hugh
% for subjects who has electrodes on both hemisphere, loading the both
% pial.mat file will generate the image with whole brain and save the
% images from 6 views (left, right, top, below, back & front)

% modified on Aug 8th, 2009 by Hugh
% show only lh(rh)'s electrodes if choose lh(rh)_surf.mat

%% Get the elec info

[FileName,PathName] = uigetfile('*.txt','Select the electrodes text file','/home/ccarlson/loc/');
[surfname, surfpath] = uigetfile('*.mat','Select the patient brain surf',PathName,'MultiSelect','on');

[name x y z] = textread([PathName, FileName],'%s %f %f %f');    
elec_cell = [name,num2cell(x),num2cell(y),num2cell(z)];

%% Get the surf info
b = findstr(FileName,'_');
Pname = FileName(1:b(1));
if iscell(surfname)
    sph = 'both';
else
    a = findstr(surfname,'_')+1;
    sph = surfname(a(1):a(1)+1);
end

%% get lh(rh) electrodes only
if strcmp(sph,'lh')
    c = find(cell2mat(elec_cell(:,2))>0);
elseif strcmp(sph,'rh')
    c = find(cell2mat(elec_cell(:,2))<0);
elseif strcmp(sph,'both')
    c = [];
end
elec_cell(c,:) = [];

g = strmatch('G',upper(elec_cell(:,1)));
d = strmatch('D',upper(elec_cell(:,1)));
if isempty(g)
    elec_grid = [];
else
    elec_grid = elec_cell(g,:);
    elec_cell([g;d],:) = [];
end


%% Plot the elecs
plt = menu('What part do you want to plot?','Grid only', 'Strip only',...
    'Both Grid and Strip');
labelshow = menu('Do you want to show the labels?','Yes','No');
genimg = menu('Do you want to save the images?','Yes', 'No');

if strcmp(sph,'both')
    surf_brain.sph1 = load([surfpath surfname{1}]);
    surf_brain.sph2 = load([surfpath surfname{2}]);
else 
    surf_brain = load([surfpath surfname]);
end

if plt==1 && ~isempty(elec_grid)
    nyu_plot(surf_brain,sph,cell2mat(elec_grid(:,2:4)),char(elec_grid(:,1)),'r',labelshow);
elseif plt==2 && ~isempty(elec_cell)
    nyu_plot(surf_brain,sph,cell2mat(elec_cell(:,2:4)),char(elec_cell(:,1)),'b',labelshow);
elseif plt==3 && ~isempty(elec_grid) && ~isempty(elec_cell)
    elec = cell2mat(elec_cell(:,2:4));
    elec_name = char(elec_cell(:,1));
    nyu_plot(surf_brain,sph,cell2mat(elec_grid(:,2:4)),char(elec_grid(:,1)),'r',labelshow); hold on;
    for i=1:length(elec)
        plot3(elec(i,1),elec(i,2),elec(i,3),'bo','MarkerFaceColor','b');
        if labelshow==1
            text('Position',[elec(i,1) elec(i,2) elec(i,3)],'String',elec_name(i,:),'Color','w');
        end
    end
    hold off;   
else
    disp('not correct input data');
    return;
end

%% save images
if genimg==1
    ext = 'png';  
    if strcmp(sph,'lh')
        view(270, 0);
        saveas(gcf,[PathName,Pname,'T1_lateral_',sph,'.',ext]);
        view(90,0);
        saveas(gcf,[PathName,Pname,'T1_mesial_',sph,'.',ext]);
        
    elseif strcmp(sph,'rh')
        view(270, 0);
        saveas(gcf,[PathName,Pname,'T1_mesial_',sph,'.',ext]);
        view(90,0);
        saveas(gcf,[PathName,Pname,'T1_lateral_',sph,'.',ext]);
        
    elseif strcmp(sph,'both')
        view(270, 0);
        saveas(gcf,[PathName,Pname,'T1_left_',sph,'.',ext]);
        view(90,0);
        saveas(gcf,[PathName,Pname,'T1_right_',sph,'.',ext]);
    end
    view(0,0);
    saveas(gcf,[PathName,Pname,'T1_posterior_',sph,'.',ext]);

    view(180,0);
    saveas(gcf,[PathName,Pname,'T1_frontal_',sph,'.',ext]);

    view(90,90);
    saveas(gcf,[PathName,Pname,'T1_dorsal_',sph,'.',ext]);

    view(90,-90);
    set(light,'Position',[1 0 -1]);
    saveas(gcf,[PathName,Pname,'T1_ventral_',sph,'.',ext]);
else 
    return;
end

end
%% nyu_plot
function nyu_plot(surf_brain,sph,elec,elecname,color,label)

col = [.7 .7 .7];

if strcmp(sph,'both')
    sub_sph1.vert = surf_brain.sph1.coords;
    sub_sph1.tri = surf_brain.sph1.faces;

    sub_sph2.vert = surf_brain.sph2.coords;
    sub_sph2.tri = surf_brain.sph2.faces;
    
    col1=repmat(col(:)', [size(sub_sph1.vert, 1) 1]);
    col2=repmat(col(:)', [size(sub_sph2.vert, 1) 1]);
    
    trisurf(sub_sph1.tri, sub_sph1.vert(:, 1), sub_sph1.vert(:, 2),sub_sph1.vert(:, 3),...
        'FaceVertexCData', col1,'FaceColor', 'interp');

    hold on;
    trisurf(sub_sph2.tri, sub_sph2.vert(:, 1), sub_sph2.vert(:, 2), sub_sph2.vert(:, 3),...
        'FaceVertexCData', col2,'FaceColor', 'interp');

else    
    if isfield(surf_brain,'coords')==0
        sub.vert = surf_brain.surf_brain.coords;
        sub.tri = surf_brain.surf_brain.faces;
    else
        sub.vert = surf_brain.coords;
        sub.tri = surf_brain.faces;
    end
    col=repmat(col(:)', [size(sub.vert, 1) 1]);
    trisurf(sub.tri, sub.vert(:, 1), sub.vert(:, 2), sub.vert(:, 3), ...
         'FaceVertexCData', col,'FaceColor', 'interp');

end

shading interp;
lighting gouraud;
material dull;
light;
axis off
hold on;
for i=1:length(elec)
    plot3(elec(i,1),elec(i,2),elec(i,3),[color 'o'],'MarkerFaceColor',color);
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