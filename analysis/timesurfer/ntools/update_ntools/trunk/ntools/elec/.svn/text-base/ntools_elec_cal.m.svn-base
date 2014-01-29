function elec = ntools_elec_cal(name,ini_cell,mode,SubjectPath)
% calculate the elec location from initial cell file. 2 modes: grid and
% depth. for the grid, it will use the sph and Subject_path parameter,
% depth doesn't. 

if ~exist('mode','var')
    mode = 'depth';
    SubjectPath = [];
end

%% depth
if strcmp(mode,'depth')
    for i = 1:length(name)
        elec_temp = cell(size(ini_cell));
%         n = strfind(ini_cell(:,1),name{i});
        n = regexp(ini_cell(:,1),[name{i} '[\d]'],'match');
        k = 1;
        for l = 1:length(n)
            if ~isempty(n{l})
                elec_temp(k,:) = ini_cell(l,:);
                k = k+1;
            end
        end
        elec_temp(all(cellfun(@isempty,elec_temp),2),:) = [];
        elec_num = regexp(elec_temp(:,1),'[^A-Za-z]*[\d*]','match');
        elec_num(all(cellfun(@isempty,elec_num),2),:) = [];
        if length(elec_num)~=2
            error('only 2 initial positions are required');
            return;
        end
        elec_ini_loc = elec_temp(:,2:4);
        E_temp = ntools_elec_position(cell2mat(elec_ini_loc),cell2num(elec_num{1}),...
            cell2num(elec_num{2}));
        elec.(char(name{i})) = E_temp;
        % clear the unnecessary data
        clear elec_temp elec_num E_temp elec_ini_pos elec_ini_loc;
    end
    
%% Grid    
elseif strcmp(mode,'grid')
    %  choose the bem surface
    bem_option = menu('Do you want to create a bem surface or choose one?',...
        'Create a new one','BEM','LGI','I have another bem surface');
    if bem_option==1
        [bem_name,bem_surf_path] = ntools_elec_make_tri_file(SubjectPath);
    elseif bem_option==2
        bem_name = 'inner_skull';
        bem_surf_path = [SubjectPath '/bem/'];
    elseif bem_option==3
        bem_surf_path = [SubjectPath '/surf/'];
    else
        [bem_name bem_surf_path] = uigetfile([SubjectPath '/*'],'Choose a bem surface');
    end
    % get the grid initial positions by name
    for i = 1:length(name)
        elec_temp = cell(size(ini_cell));
%         n = strfind(ini_cell(:,1),name{i});
        n = regexp(ini_cell(:,1),[name{i} '[\d]'],'match');
        k = 1;
        for l = 1:length(n)
            if ~isempty(n{l})
                elec_temp(k,:) = ini_cell(l,:);
                k = k+1;
            end
        end
        elec_temp(all(cellfun(@isempty,elec_temp),2),:) = [];
        elec_num = regexp(elec_temp(:,1),'[^A-Za-z]*[\d*]','match');
        elec_num(all(cellfun(@isempty,elec_num),2),:) = [];
        if length(elec_num)~=2
            error('only 2 initial positions are required');
            break;
        end
        elec_ini_pos = [str2double(cell2mat(elec_num{1})); str2double(cell2mat(elec_num{2}))];
        elec_ini_loc = cell2mat(elec_temp(:,2:4));
        
        % determine the hemisphere that grid locates
        if elec_ini_loc(:,1)>0
            sph = 'rh';
        elseif elec_ini_loc(:,1)<0
            sph = 'lh';
        else
            error(['wrong initial positions for grid ', name{i}])
            break;
        end

        % determin the hemisphere if using LGI
        if bem_option==3
            bem_name = [sph '.pial-outer-smoothed'];
        end

        % input the size of the grid
        a = menu(['What is the size of the grid ',name{i},' ?'], '8*8', '4*8', 'manually input');
        if a==1
            s = [8,8];
        elseif a==2
            s = [4,8];
        else
            s = input(['Please input the size of the grid ',name{i},' [row column]: ']);
        end

        % calcuate and project the grid on the bem
        E_temp = ntools_elec_bemproj(elec_ini_loc,elec_ini_pos,s,sph,bem_name,bem_surf_path);
        elec.(char(name{i})) = E_temp;
        
        % clear the unnecessary data
        clear elec_temp elec_num E_temp elec_ini_pos elec_ini_loc;
    end
else 
    elec = [];
end

