function ico    = ts_mgh2ico( infile, subject, hemi, varargin )
%
%   Matlab wrapper to facilitate calls to mri_surf2surf
%
% Required Parameters:
%   subject
%   hemi
%
% Optional Parameters:
%
% [I/O args:]
%   prefix 
%   inpath
%   outfile
%   outdir
%
% [Analysis args:]
%   icolevel
%   sphsmoothsteps
%
%
% Created:       05/30/2007 by Ben Cipollini
% Last Modified: 06/27/2007 by Ben Cipollini
%

    
    
    parms           = args2parms( varargin, ...
	                              { 'prefix',     '', [],...
                                    'sphsmoothsteps', 0, [0 Inf],...
                                    'outtype',    'mgh', [],...
                                    'outdir',    pwd, [],...
                                    ...
                                    'subjectsdir',getenv('SUBJECTS_DIR'), [], ...
                                    'outfile',    '', [],...
                                    ...
                                    'icolevel',   7,  [1 10],...
                                    ...
                                    'force',      true, sort([false true]),...
                                    'verbose',    true, sort([false true])}...
                                  );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sample to sphere
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    
    if (isempty(parms.outfile))
        if (parms.sphsmoothsteps ~= 1)
            prefix_out = sprintf('%s-sphere-sm%d',parms.prefix, parms.sphsmoothsteps);
        else
            prefix_out = sprintf('%s-sphere',parms.prefix);
        end;

        parms.outfile = sprintf('%s-%s.%s', prefix_out,   hemi, parms.outtype);
    end;
    
    parms.outfile   = fullfile( parms.outdir, parms.outfile );

    
    %   Run the unix command
    if (parms.force || ~exist(parms.outfile,'file'))

        cmd = sprintf('mri_surf2surf --srcsubject %s --trgsubject ico --hemi %s', subject, hemi);
        cmd = sprintf('%s  --trgicoorder %d --sval %s --tval %s', cmd, parms.icolevel, infile, parms.outfile);
        cmd = sprintf('%s --sfmt %s --nsmooth-out %d', cmd, parms.outtype, parms.sphsmoothsteps);
        cmd = sprintf('%s --noreshape', cmd);

        %if (parms.verbose), display(cmd); end;

        oldsdir     = getenv('SUBJECTS_DIR');
        setenv('SUBJECTS_DIR', parms.subjectsdir);
        
        [status,result]=unix(cmd);

        setenv('SUBJECTS_DIR', oldsdir);

        %   Validate the result
        if (status || ~isempty(findstr(result, 'could not open')))
            error('%s: %s\n\ncmd:\n%s',mfilename, result, cmd);
        end;
    end;
    
    
    %   Return the result if necessary
    if (nargout==1)
        ico     = fs_load_mgh( parms.outfile );
    end;
