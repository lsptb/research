function [roinums,roilabels] = fs_read_annotation(fname)
%function [roinums,roilabels] = fs_read_annotation(fname)
%
% Purpose: reads annotation file
%
% Input:
%   fname: full or relative path of annotation file
%
% Output:
%   roinums: vector of roi numbers (one for each vertex)
%   roilabels: cell array of roi labels (one for every unique roinum)
%
% Created:  02/01/07 by Don Hagler
%   based on read_annotation.m distributed with freesurfer
% Last Mod: 02/01/07 by Don Hagler
%
%

fp = fopen(fname, 'r', 'b');

if(fp < 0)
  fprintf('%s: ERROR: unable to open file %s\n',mfilename,fname);
  return;
end

A = fread(fp, 1, 'int');

tmp = fread(fp, 2*A, 'int');
vertices = tmp(1:2:end);
codevec = tmp(2:2:end);

bool = fread(fp, 1, 'int');
if(isempty(bool)) %means no colortable
  fprintf('%s: ERROR: No colortable found in %s\n',mfilename,fname);
  colortable = struct([]);
  fclose(fp);
  return; 
end

if(bool)
    %Read colortable
    numEntries = fread(fp, 1, 'int');
    if(numEntries > 0)
        colortable.numEntries = numEntries;
        len = fread(fp, 1, 'int');
        colortable.orig_tab = fread(fp, len, '*char')';
        colortable.orig_tab = colortable.orig_tab(1:end-1);
        colortable.struct_names = cell(numEntries,1);
        colortable.table = zeros(numEntries,5);
        for i = 1:numEntries
            len = fread(fp, 1, 'int');
            colortable.struct_names{i} = fread(fp, len, '*char')';
            colortable.struct_names{i} = colortable.struct_names{i}(1:end-1);
            colortable.table(i,1) = fread(fp, 1, 'int');
            colortable.table(i,2) = fread(fp, 1, 'int');
            colortable.table(i,3) = fread(fp, 1, 'int');
            colortable.table(i,4) = fread(fp, 1, 'int');
            colortable.table(i,5) = colortable.table(i,1) + colortable.table(i,2)*2^8 + colortable.table(i,3)*2^16 + colortable.table(i,4)*2^24;
        end
%        disp(['colortable with ' num2str(colortable.numEntries) ' entries read (originally ' colortable.orig_tab ')']);
    else
        version = -numEntries;
        if(version~=2)    
          fprintf('%s: ERROR: version %d not supported\n',mfilename,version);
        else
          fprintf('%s: reading color table version %d\n',mfilename,version);
        end
        numEntries = fread(fp, 1, 'int');
        colortable.numEntries = numEntries;
        len = fread(fp, 1, 'int');
        colortable.orig_tab = fread(fp, len, '*char')';
        colortable.orig_tab = colortable.orig_tab(1:end-1);
        colortable.struct_names = cell(numEntries,1);
        colortable.table = zeros(numEntries,5);
        numEntriesToRead = fread(fp, 1, 'int');
        for i = 1:numEntriesToRead
            structure = fread(fp, 1, 'int')+1;
            if (structure < 0)
              fprintf('%s: ERROR: bad structure %d\n',mfilename,structure);
            end
            if(~isempty(colortable.struct_names{structure}))
              fprintf('%s: ERROR: duplicated structure %d\n',mfilename,structure);
            end
            len = fread(fp, 1, 'int');
            colortable.struct_names{structure} = fread(fp, len, '*char')';
            colortable.struct_names{structure} = colortable.struct_names{structure}(1:end-1);
            colortable.table(structure,1) = fread(fp, 1, 'int');
            colortable.table(structure,2) = fread(fp, 1, 'int');
            colortable.table(structure,3) = fread(fp, 1, 'int');
            colortable.table(structure,4) = fread(fp, 1, 'int');
            colortable.table(structure,5) = colortable.table(structure,1) + colortable.table(structure,2)*2^8 + colortable.table(structure,3)*2^16 + colortable.table(structure,4)*2^24;       
        end
%        disp(['colortable with ' num2str(colortable.numEntries) ' entries read (originally ' colortable.orig_tab ')']);
    end    
else
  fprintf('%s: ERROR: bad value in annotation file\n',mfilename);
end

fclose(fp);

codes = colortable.table(:,5);
nrois = length(codes);
roilabels = colortable.struct_names;
roinums = zeros(size(codevec));
for i=1:nrois
  roinums(codevec==codes(i))=i;
end;

