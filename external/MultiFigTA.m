function HF = MultiFigTA(FuncHandle,Vars)
%
% Creates subplotted figures using any function handle (e.g. plot, imagesc,
% etc) for all the fields of a structure, where each field is a cell array
% of inputs for the chosen function.  It also works for arrays of cells,
% with functions requiring multiple inputs inputted as arrays of nested
% cells.  It also works for simple 2D array inputs, but only on functions
% requiring one input (e.g. plot(y)).
%
% egs:
% FuncHandle = @plot
% VarsStruct.f1...fn = {[input1] [input2]}
% VarsCell{{1x2}...{1x2}} = {input1} {input2}
% VarsDouble[] = columns to subplot individually
%

switch class(Vars)
    
    case 'struct'
        
        fn = fieldnames(Vars);
        nfiles = length(fn);
        nfr = rem(nfiles,8);
        if nfr==0
            nf = floor(nfiles/8);
        else nf = floor(nfiles/8)+1;
        end
        c1 = 1;
        for k = 1:nf
            fm
            if k==nf, subnum = nfr; else subnum = 8; end
            for k2 = 1:subnum
                if c1<=nfiles
                    name = fn{c1};
                    if iscell(Vars.(name))
                        subplot(2,4,k2), feval(FuncHandle,Vars.(name){:});
                    else subplot(2,4,k2), feval(FuncHandle,Vars.(name));
                    end
                    axis tight, title(name);
                    c1 = c1+1;
                end
            end
        end
        
    case 'cell'
        
        nfiles = length(Vars);
        nfr = rem(nfiles,8);
        if nfr==0
            nf = floor(nfiles/8);
        else nf = floor(nfiles/8)+1;
        end
        c1 = 1;
        for k = 1:nf
            fm
            if k==nf, subnum = nfr; else subnum = 8; end
            for k2 = 1:subnum
                name = num2str(c1);
                if iscell(Vars{c1})
                    subplot(2,4,k2), feval(FuncHandle,Vars{c1}{:});
                elseif isstruct(Vars{c1})
                    if length(Vars{c1})>1
                        subplot(2,4,k2), hold on, structfun(@(x) feval(FuncHandle,x),Vars{c1}); axis tight
                    else fn2 = fieldnames(Vars{c1});
                        val_fn2 = zeros(length(fn2),1);
                        val_fn2(:) = structfun(@(x) x,Vars{c1});
                        subplot(2,4,k2), feval(FuncHandle,val_fn2); axis tight
                    end
                else subplot(2,4,k2), feval(FuncHandle,Vars{c1});
                end
                title(name);
                c1 = c1+1;
            end
        end
        
    case 'double'
        
        nfiles = size(Vars,2);
        nfr = rem(nfiles,8);
        if nfr==0
            nf = floor(nfiles/8);
        else nf = floor(nfiles/8)+1;
        end
        c1 = 1;
        for k = 1:nf
            fm
            if k==nf, subnum = nfr; else subnum = 8; end
            for k2 = 1:subnum
                if c1<=nfiles
                    name = num2str(c1);
                    subplot(2,4,k2), feval(FuncHandle,Vars(:,c1));
                    axis tight, title(name);
                    c1 = c1+1;
                end
            end
        end
        
    otherwise, disp([class ': dunno'])
        
end

HF = findobj('Type','figure');
HF = findobj(HF,'Type','axes');

end

