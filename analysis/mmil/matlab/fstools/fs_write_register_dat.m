function [status]=fs_write_register_dat(fname,M,subj,inplane,slicethick)
%function [status]=fs_write_register_dat(fname,[M],[subj],[inplane],[slicethick])

if nargin<1
  help(mfilename);
  return;
end;

if ~exist('M','var') || isempty(M), M=eye(4); end;
if ~exist('subj','var') || isempty(subj), subj='unknown'; end;
if ~exist('inplane','var') || isempty(inplane), inplane=1; end;
if ~exist('slicethick','var') || isempty(slicethick), slicethick=1; end;


status = 0;

fid=fopen(fname,'wt');
if fid==-1
  status=1;
  return;
end;
fprintf(fid,'%s\n',subj);
fprintf(fid,'%f\n',inplane);
fprintf(fid,'%f\n',slicethick);
fprintf(fid,'1\n');
fprintf(fid,'%f ',M(1,:));
fprintf(fid,'\n');
fprintf(fid,'%f ',M(2,:));
fprintf(fid,'\n');
fprintf(fid,'%f ',M(3,:));
fprintf(fid,'\n');
fprintf(fid,'%f ',M(4,:));
fprintf(fid,'\n');
fprintf(fid,'round\n');
fclose(fid);

