function location = ntools_elec_position(ini,row,col)

% ------------------------------------
% This program is used to calculate the position of 
% every matrix element in a 3D or 2D plane or line
% 
%     A     col     B
%      *   *   *   *
%  row *   *   *   *
%      *   *   *   *
%      *   *   *   *
%     C             D
% 
% Input:
% ini: the initial 3 points of the plane, A B C, row vectors.
%      for example:[0,0,0; 0,4,4; 4,0,0];
% col: the column of the matrix, say, col=5;
% row: the row of the matrix, say, row=5;
% 
% Special Case: if ini is a 2*3 or 2*2 matrix, that means all the elec are in a
% line. for this case, row = 1;
% 
% Output:
% pos: 3D vectors. Each slice contains the positions of the row points. 
% location: 2D vectors, stored by row.
% ------------------------------------

form = size(ini);
D3 = [3,3]; D2 = [3,2]; D1 = [2,3]; D0 = [2,2];

if form==D3 | form==D2
    % Validation: 3 initial points are in the same line, they cann't form a plane
    tan_xy_AB = (ini(2,2)-ini(1,2))/(ini(2,1)-ini(1,1));
    tan_xy_BC = (ini(3,2)-ini(2,2))/(ini(3,1)-ini(2,1));
    if form==D3
        tan_xz_AB = (ini(2,3)-ini(1,3))/(ini(2,1)-ini(1,1));
        tan_xz_BC = (ini(3,3)-ini(2,3))/(ini(3,1)-ini(2,1));
        tan_yz_AB = (ini(2,3)-ini(1,3))/(ini(2,2)-ini(1,2));
        tan_yz_BC = (ini(3,3)-ini(2,3))/(ini(3,2)-ini(2,2));
        if tan_xy_AB==tan_xy_BC & tan_xz_AB==tan_xz_BC & tan_yz_AB==tan_yz_BC
%             error('The initial points are in the same line!');
             location = [];return
        end
    else if tan_xy_AB==tan_xy_BC
%             error('The initial points are in the same line!');
            location = [];return
        end
    end

    % Get the parameter
    pos(1,:,1) = ini(1,:);  % point A
    pos(col,:,1) = ini(2,:);  % point B
    pos(1,:,row) = ini(3,:);  % point C

    % point D
    AB = ini(2,:)-ini(1,:);
    AC = ini(3,:)-ini(1,:);
    AD = AB+AC;
    pos(col,:,row) = AD+ini(1,:);

    % Calculate the position of the first and last elec of each row
    for j = 2:row-1
        l = (j-1)/(row-j);
        pos(1,:,j) = (pos(1,:,1)+l*pos(1,:,row))/(1+l);
        pos(col,:,j) = (pos(col,:,1)+l*pos(col,:,row))/(1+l);
    end

    %% Determine the positions of the rest elec
    for j = 1:row
        for i = 2:col-1
            l = (i-1)/(col-i);
            pos(i,:,j) = (pos(1,:,j)+l*pos(col,:,j))/(1+l);
        end
    end

    %% Restore the pos to location, 2D vector
    k = 1;
    for i=1:row
        for j=1:col
            location(k,:) = pos(j,:,i);
            k = k+1;
        end
    end
%     save location.txt location -ascii -tabs;

elseif form==D1 | form==D0     % calculate the position in line
    if row ~=1
        error('There are only 2 initial positions, please check your row number');
        return;
    end
    location(1,:) = ini(1,:);  % point A
    location(col,:) = ini(2,:);  % point B
    for j=2:col-1
        l=(j-1)/(col-j);
        location(j,:)=(location(1,:)+l*location(col,:))/(1+l);
    end
%     save location.txt location -ascii -tabs;
    
else
    error('ini must be one of a 3*3, 3*2, 2*3, 2*2 matrix');
end

return
