function [f] = invtrad(H, norm);

% INVTRAD creates  a 9 parameter traditional "Talairach-model"
%   transformation from the homogenous spatial transformation matrix
%
% f = traditional(H)
%
% The transformation vector f contains the 
%	x-shift
%	y-shift
%	z-shift
% followed by the
%	pitch (rotation around x-axis)
%	roll  (rotation around y-axis)
%	yaw   (rotation around z-axis)
% followed by the 
%	x-rescaling factor
%	y-rescaling factor
%	z-rescaling factor
%
% See also WARP3D

% Copyright (C) 2005, Arnaud Delorme from the function traditional of 
%                     R. Oostenveld
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% compute the homogenous transformation matrix for the translation
    
    % translation
    % -----------
    f(1) = H(1,4); H(1,4) = 0;
    f(2) = H(2,4); H(2,4) = 0;
    f(3) = H(3,4); H(3,4) = 0;
    R = H(1:3,1:3);
            
    % analitic solution
    % -----------------
    R11 = R(1,1);
    R12 = R(1,2);
    R13 = R(1,3);
    R21 = R(2,1);
    R22 = R(2,2);
    R23 = R(2,3);
    R31 = R(3,1);
    R32 = R(3,2);
    R33 = R(3,3);
    tmpH2 = sum(R.^2);
    if nargin < 2
        f1 = sqrt(tmpH2(1));        
        f3 = 1/6*(6/(f1^2-R21^2)*(-12*R23^4*f1^6*R13^2+12*R23^6*f1^4*R21^2+24*R23^4*f1^6*R33^2+12*R23^6*f1^2*R21^4-24*R23^2*f1^6*R33^4+24*R13^6*R21^4*f1^2-24*R13^6*R21^2*f1^4+8*f1^6*R33^6-8*R13^6*R21^6-8*R21^6*R23^6-12*R23^4*f1^4*R13^2*R21^2-8*R23^6*f1^6+8*f1^6*R13^6+36*R23^2*f1^4*R13^2*R21^2*R33^2- ...
                              24*R23^2*f1^2*R13^2*R33^2*R21^4+60*R23^4*f1^2*R21^4*R33^2-48*R23^2*f1^4*R13^4*R21^2+60*R23^2*f1^2*R13^4*R21^4-48*R13^4*R21^2*f1^4*R33^2+24*R13^4*R21^4*f1^2*R33^2-24*R13^2*R21^2*f1^4*R33^4+12*R23^2*f1^6*R13^4+48*R23^4*f1^2*R13^2*R21^4-60*R23^4*f1^4*R33^2*R21^2-12*R23^2*f1^6*R13^2*R33^2+48*R23^2*f1^4*R33^4*R21^2-24*R13^4*R21^6*R23^2-24*R13^2*R21^6*R23^4-12*R23^2*f1*(-36*R23^4*R33^2*R21^6*R13^2-36*R23^2*R33^2*R21^6*R13^4-3*f1^6*R33^4*R13^4-6*f1^6*R13^6*R23^2-6*f1^6*R33^2*R13^6-3*f1^6*R13^4*R23^4+6*f1^4*R13^8*R21^2-3*f1^2*R23^8*R21^4-3*f1^2*R13^8*R21^4-12*R13^6*R21^6*R33^2+18*f1^4*R13^6*R23^2*R21^2+6*f1^6*R13^4*R23^2*R33^2-3*f1^6*R13^8-12*R23^6*R33^2*R21^6+18*f1^4*R13^4*R23^4*R21^2-24*f1^4*R13^2*R21^2*R23^4*R33^2-12*f1^4*R13^2*R33^6*R21^2-24*f1^4*R13^4*R21^2*R33^4+6*f1^4*R13^2*R23^6*R21^2-6*f1^4*R13^6*R21^2*R33^2-18*f1^2*R13^4*R23^4*R21^4+24*f1^2*R13^4*R21^4*R33^4+6*f1^2*R33^2*R23^6*R21^4+24*f1^2*R13^6*R21^4*R33^2-12*f1^2*R13^6*R23^2*R21^4-12*f1^2*R13^2*R23^6*R21^4-3*f1^2*R23^4*R33^4*R21^4+30*f1^4*R13^2*R23^2*R21^2*R33^4-30*f1^4*R13^4*R23^2*R33^2*R21^2+36*f1^2*R13^2*R23^4*R21^4*R33^2+54*f1^2*R13^4*R23^2*R21^4*R33^2-60*f1^2*R23^2*R33^4*R21^4*R13^2)^(1/2)*R21^2+12*R23^2*f1^3*(-36*R23^4*R33^2*R21^6*R13^2-36*R23^2*R33^2*R21^6*R13^4-3*f1^6*R33^4*R13^4-6*f1^6*R13^6*R23^2-6*f1^6*R33^2*R13^6-3*f1^6*R13^4*R23^4+6*f1^4*R13^8*R21^2-3*f1^2*R23^8*R21^4-3*f1^2*R13^8*R21^4-12*R13^6*R21^6*R33^2+18*f1^4*R13^6*R23^2*R21^2+6*f1^6*R13^4*R23^2*R33^2-3*f1^6*R13^8-12*R23^6*R33^2*R21^6+18*f1^4*R13^4*R23^4*R21^2-24*f1^4*R13^2*R21^2*R23^4*R33^2-12*f1^4*R13^2*R33^6*R21^2-24*f1^4*R13^4*R21^2*R33^4+6*f1^4*R13^2*R23^6*R21^2-6*f1^4*R13^6*R21^2*R33^2-18*f1^2*R13^4*R23^4*R21^4+24*f1^2*R13^4*R21^4*R33^4+6*f1^2*R33^2*R23^6*R21^4+24*f1^2*R13^6*R21^4*R33^2-12*f1^2*R13^6*R23^2*R21^4-12*f1^2*R13^2*R23^6*R21^4-3*f1^2*R23^4*R33^4*R21^4+30*f1^4*R13^2*R23^2*R21^2*R33^4-30*f1^4*R13^4*R23^2*R33^2*R21^2+36*f1^2*R13^2*R23^4*R21^4*R33^2+54*f1^2*R13^4*R23^2*R21^4*R33^2-60*f1^2*R23^2*R33^4*R21^4*R13^2)^(1/2)+24*f1^6*R13^2*R33^4+24*f1^6*R13^4*R33^2)^(1/3)+24*(R23^4*f1^4-R23^4*f1^2*R21^2-2*R23^2*f1^4*R33^2+4*R23^2*f1^2*R33^2*R21^2+R23^2*f1^4*R13^2-3*R23^2*f1^2*R13^2*R21^2-2*R13^4*R21^2*f1^2+2*f1^4*R13^2*R33^2+2*R13^2*R21^4*R23^2+R13^4*R21^4+R21^4*R23^4-2*R13^2*R21^2*f1^2*R33^2+f1^4*R13^4+f1^4*R33^4)/(f1^2-R21^2)/(-12*R23^4*f1^6*R13^2+12*R23^6*f1^4*R21^2+24*R23^4*f1^6*R33^2+12*R23^6*f1^2*R21^4-24*R23^2*f1^6*R33^4+24*R13^6*R21^4*f1^2-24*R13^6*R21^2*f1^4+8*f1^6*R33^6-8*R13^6*R21^6-8*R21^6*R23^6-12*R23^4*f1^4*R13^2*R21^2-8*R23^6*f1^6+8*f1^6*R13^6+36*R23^2*f1^4*R13^2*R21^2*R33^2-24*R23^2*f1^2*R13^2*R33^2*R21^4+60*R23^4*f1^2*R21^4*R33^2-48*R23^2*f1^4*R13^4*R21^2+60*R23^2*f1^2*R13^4*R21^4-48*R13^4*R21^2*f1^4*R33^2+24*R13^4*R21^4*f1^2*R33^2-24*R13^2*R21^2*f1^4*R33^4+12*R23^2*f1^6*R13^4+48*R23^4*f1^2*R13^2*R21^4-60*R23^4*f1^4*R33^2*R21^2-12*R23^2*f1^6*R13^2*R33^2+48*R23^2*f1^4*R33^4*R21^2-24*R13^4*R21^6*R23^2-24*R13^2*R21^6*R23^4-12*R23^2*f1*(-36*R23^4*R33^2*R21^6*R13^2-36*R23^2*R33^2*R21^6*R13^4-3*f1^6*R33^4*R13^4-6*f1^6*R13^6*R23^2-6*f1^6*R33^2*R13^6-3*f1^6*R13^4*R23^4+6*f1^4*R13^8*R21^2-3*f1^2*R23^8*R21^4-3*f1^2*R13^8*R21^4-12*R13^6*R21^6*R33^2+18*f1^4*R13^6*R23^2*R21^2+6*f1^6*R13^4*R23^2*R33^2-3*f1^6*R13^8-12*R23^6*R33^2*R21^6+18*f1^4*R13^4*R23^4*R21^2-24*f1^4*R13^2*R21^2*R23^4*R33^2-12*f1^4*R13^2*R33^6*R21^2-24*f1^4*R13^4*R21^2*R33^4+6*f1^4*R13^2*R23^6*R21^2-6*f1^4*R13^6*R21^2*R33^2-18*f1^2*R13^4*R23^4*R21^4+24*f1^2*R13^4*R21^4*R33^4+6*f1^2*R33^2*R23^6*R21^4+24*f1^2*R13^6*R21^4*R33^2-12*f1^2*R13^6*R23^2*R21^4-12*f1^2*R13^2*R23^6*R21^4-3*f1^2*R23^4*R33^4*R21^4+30*f1^4*R13^2*R23^2*R21^2*R33^4-30*f1^4*R13^4*R23^2*R33^2*R21^2+36*f1^2*R13^2*R23^4*R21^4*R33^2+54*f1^2*R13^4*R23^2*R21^4*R33^2-60*f1^2*R23^2*R33^4*R21^4*R13^2)^(1/2)*R21^2+12*R23^2*f1^3*(-36*R23^4*R33^2*R21^6*R13^2-36*R23^2*R33^2*R21^6*R13^4-3*f1^6*R33^4*R13^4-6*f1^6*R13^6*R23^2-6*f1^6*R33^2*R13^6-3*f1^6*R13^4*R23^4+6*f1^4*R13^8*R21^2-3*f1^2*R23^8*R21^4-3*f1^2*R13^8*R21^4-12*R13^6*R21^6*R33^2+18*f1^4*R13^6*R23^2*R21^2+6*f1^6*R13^4*R23^2*R33^2-3*f1^6*R13^8-12*R23^6*R33^2*R21^6+18*f1^4*R13^4*R23^4*R21^2-24*f1^4*R13^2*R21^2*R23^4*R33^2-12*f1^4*R13^2*R33^6*R21^2-24*f1^4*R13^4*R21^2*R33^4+6*f1^4*R13^2*R23^6*R21^2-6*f1^4*R13^6*R21^2*R33^2-18*f1^2*R13^4*R23^4*R21^4+24*f1^2*R13^4*R21^4*R33^4+6*f1^2*R33^2*R23^6*R21^4+24*f1^2*R13^6*R21^4*R33^2-12*f1^2*R13^6*R23^2*R21^4-12*f1^2*R13^2*R23^6*R21^4-3*f1^2*R23^4*R33^4*R21^4+30*f1^4*R13^2*R23^2*R21^2*R33^4-30*f1^4*R13^4*R23^2*R33^2*R21^2+36*f1^2*R13^2*R23^4*R21^4*R33^2+54*f1^2*R13^4*R23^2*R21^4*R33^2-60*f1^2*R23^2*R33^4*R21^4*R13^2)^(1/2)+24*f1^6*R13^2*R33^4+24*f1^6*R13^4*R33^2)^(1/3)+12*(-R13^2*R21^2+f1^2*R13^2+2*R23^2*f1^2+f1^2*R33^2-R21^2*R23^2)/(f1^2-R21^2))^(1/2);
        f3 = real(f3);
        f2 = -R22*(R13^2+R23^2-f3^2)/(-(R13^2+R23^2-f3^2)*(-R23^2+f3^2))^(1/2)/R33*f3^2/(-R23^2+f3^2)^(1/2); f2 = real(f2);
    else
        f1 = 1;
        f2 = 1;
        f3 = 1;
    end;
    cX = (-R23^2+f3^2)^(1/2)/f3;
    sX = R23/f3;
    cY = -(R13^2+R23^2-f3^2)/(-(R13^2+R23^2-f3^2)*(-R23^2+f3^2))^(1/2);
    sY = R13/(-R23^2+f3^2)^(1/2);
    cZ = -1/(R13^2+R23^2-f3^2)*(-(R13^2+R23^2-f3^2)*(-R23^2+f3^2))^(1/2)*R33/f3;
    sZ = 1/(R13^2+R23^2-f3^2)*((R13^2+R23^2-f3^2)*(R13^2*f3^2-R33^2*R23^2+R33^2*f3^2+f3^2*R23^2-f3^4))^(1/2)/f3;
    
    f(4) = atan2( sX, cX );
    f(5) = atan2( sY, cY );
    f(6) = atan2( sZ, cZ );
    f(7) = f1;
    f(8) = f2;
    f(9) = f3;

    % test for f6 or -f6
    dist1 = sum(sum( (traditional(f) - H).^2 ));
    f(6) = -f(6);
    dist2 = sum(sum( (traditional(f) - H).^2 ));
    if dist1 < dist2, f(6) = -f(6); end;
    return;
 
    if 0
        % solving for rotation
        equ{1} = 'R11 = (cZ*cY + sZ*sX*sY)*f1';
        equ{2} = 'R12 = (sZ*cY + cZ*sX*sY)*f2';
        equ{3} = 'R13 = (           cX*sY)*f3';
        equ{4} = 'R21 = (-sZ*cX)*f1';
        equ{5} = 'R22 = (cZ*cX )*f2';
        equ{6} = 'R23 =       sX*f3';
        equ{7} = 'R31 = ( sZ*sX*cY - cZ*sY)*f1';
        equ{8} = 'R32 = (-cZ*sX*cY - sZ*sY)*f2';
        equ{9} = 'R33 =             (cZ*cY)*f3';
        equ{10} = 'sX^2+cX^2 = 1';
        equ{11} = 'sY^2+cY^2 = 1';
        equ{12} = 'sZ^2+cZ^2 = 1';
        vars = { 'cX' 'sX' 'cY' 'sY' 'cZ' 'sZ' 'f1' 'f2' 'f3' };
        
        S = solve( equ{:}, vars{:}) % unconstrained because of 3x 3-D rotation
        
        clear equ;
        equ{1} = 'R13 =  (cX*sY)*f3';
        equ{2} = 'R21 = (-sZ*cX)*f1';
        equ{3} = 'R22 = (cZ*cX )*f2';
        equ{4} = 'R23 =       sX*f3';
        equ{5} = 'R33 =  (cX*cY)*f3';
        equ{6} = 'sX^2+cX^2 = 1';
        equ{7} = 'sY^2+cY^2 = 1';
        equ{8} = 'sZ^2+cZ^2 = 1';
        %equ{9} = 'R31 = ( sZ*sX*cY - cZ*sY)*f1';
        %equ{10}= 'R32 = (-cZ*sX*cY - sZ*sY)*f2';
        vars = { 'cX' 'sX' 'cY' 'sY' 'cZ' 'sZ' 'f1' 'f2' 'f3' };
        
        S = solve( equ{:}, vars{:}) % OK, f3 has to be fixed
        
        S2 = solve( [ S.f1 ' =f1' ], 'f3'); % solve for f3 (since f3 is returned in solution above)
        
        R11 = R(1,1);
        R12 = R(1,2);
        R13 = R(1,3);
        R21 = R(2,1);
        R22 = R(2,2);
        R23 = R(2,3);
        R31 = R(3,1);
        R32 = R(3,2);
        R33 = R(3,3);
        
        tmpH2 = sum(H.^2);
        f1 = sqrt(tmpH2(1));
        f3 = real(eval(S2(1))); 
        f2 = eval(S.f2(3)); % to be positive
        cX = eval(S.cX(3));
        cY = eval(S.cY(3));
        cZ = eval(S.cZ(3));
        sX = eval(S.sX(3));
        sY = eval(S.sY(3));
        sZ = eval(S.sZ(3));    
        
    end;
    
    return;

    % -----------------------------------------------
    % Other analitic solution (using unitary vectors)
    % -----------------------------------------------
    
    R11 = sym('(cZ*cY + sZ*sX*sY)*f1');
    R12 = sym('(sZ*cY + cZ*sX*sY)*f2');
    R13 = sym('(cX*sY)*f3');
    R21 = sym('(-sZ*cX)*f1');
    R22 = sym('(cZ*cX )*f2');
    R23 = sym('sX*f3');
    R31 = sym('( sZ*sX*cY - cZ*sY)*f1');
    R32 = sym('(-cZ*sX*cY - sZ*sY)*f2');
    R33 = sym('(cZ*cY)*f3');
    R11 = sym('(cZ*cY + sZ*sX*sY)');
    R12 = sym('(sZ*cY + cZ*sX*sY)');
    R13 = sym('(cX*sY)');
    R21 = sym('(-sZ*cX)');
    R22 = sym('(cZ*cX )');
    R23 = sym('sX');
    R31 = sym('( sZ*sX*cY - cZ*sY)');
    R32 = sym('(-cZ*sX*cY - sZ*sY)');
    R33 = sym('(cZ*cY)');
    matrix    = [R11 R12 R13; R21 R22 R23; R31 R32 R33];
    tmparray1 = matrix*[1 0 0]'; realres1 = R*[1 0 0]';
    tmparray2 = matrix*[0 1 0]'; realres2 = R*[0 1 0]';
    tmparray3 = matrix*[0 0 1]'; realres3 = R*[0 0 1]';
    
    eqn{1} = tmparray1(1)-realres1(1);
    eqn{2} = tmparray1(2)-realres1(2);
    eqn{3} = tmparray1(3)-realres1(3);
    eqn{4} = tmparray2(1)-realres2(1);
    eqn{5} = tmparray2(2)-realres2(2);
    eqn{6} = tmparray2(3)-realres2(3);
    eqn{7} = tmparray3(1)-realres3(1);
    eqn{8} = tmparray3(2)-realres3(2);
    eqn{9} = tmparray3(3)-realres3(3);
    eqn{1} = tmparray1(1)-'a1';
    eqn{2} = tmparray1(2)-'a2';
    eqn{3} = tmparray1(3)-'a3';
    eqn{4} = tmparray2(1)-'b1';
    eqn{5} = tmparray2(2)-'b2';
    eqn{6} = tmparray2(3)-'b3';
    eqn{7} = tmparray3(1)-'c1';
    eqn{8} = tmparray3(2)-'c2';
    eqn{9} = tmparray3(3)-'c3';
    vars = { 'cX' 'sX' 'cY' 'sY' 'cZ' 'sZ' 'f1' 'f2' 'f3' };
    dfdsf
    S = solve( eqn{:}, vars{:})
    f1 = 1;
    f2 = 1;
    f3 = 1;
    S = solve( eqn{[2 3 4 5 6 7] }, vars{:});

    allcX = eval(S.cX);
    allsX = eval(S.sX);
    allcY = eval(S.cY);
    allsY = eval(S.sY);
    allcZ = eval(S.cY);
    allsZ = eval(S.sZ);

    % best matrix
    % -----------
    f(1) = 0;
    f(2) = 0;
    f(3) = 0;
    f(7) = 1;
    f(8) = 1;
    f(9) = 1;
    for ind = 1:4
        f(4) = -atan2( allsX(ind), allcX(ind) );
        f(5) = -atan2( allsY(ind), allcY(ind) );
        f(6) = -atan2( allsZ(ind), allcZ(ind) );
        dist(ind) = sum(sum( (traditional(f) - H).^2 ));
    end;
    return
    
    % ---------------------------------------------------
    % axis solutions using unitary vectors (did not work)
    % ---------------------------------------------------
    tmpx  = R*[1 0 0]';
    tmpy  = R*[0 1 0]';
    tmpz  = R*[0 0 1]';
    yaw   = -atan2(tmpx(2), tmpx(1));
    roll  = -atan2(tmpx(3), tmpx(1));
    pitch = -atan2(tmpy(3), tmpy(2));
    f(5) = roll;
 
    %H2 = traditional([ 0 0 0 0 -roll 0 1 1 1])
    %H  = inv(H2)*H; % ok for pitch and roll
    H2 = traditional([ 0 0 0 0 roll 0 1 1 1])
    H  = inv(H2)*H; % ok for roll and yaw; could not find general formula
    R = H(1:3,1:3)
    tmpx  = R*[1 0 0]';
    tmpy  = R*[0 1 0]';
    tmpz  = R*[0 0 1]';
    yaw   = -atan2(tmpx(2), tmpx(1));
    roll  = -atan2(tmpx(3), tmpx(1));
    pitch = -atan2(tmpy(3), tmpy(2));
   
    f(4) = pitch;
    f(6) = yaw;
        
    f(7) = 1;
    f(8) = 1;
    f(9) = 1;
    return;

    % old transformation
    % ------------------
    sX = R(2,3);
    cX = sqrt(1-sX^2); % contrain the sign to be positive
    sZ = -R(2,1)/cX;
    sY = R(1,3)/cX;
    systeq = [ sZ sX*sY; sZ*sX -sY ]; sol = [ R(1,2); R(3,1) ];
    res = systeq\sol;
    cY = res(1);
    cZ = res(2);
    
    f(4) = atan2( sX, cX );
    f(5) = atan2( sY, cY );
    f(6) = atan2( sZ, cZ );
    f(7) = R(1,1)/(cZ*cY + sZ*sX*sY);
    f(8) = R(2,2)/cZ/cX;
    f(9) = R(3,3)/cZ/cY;
