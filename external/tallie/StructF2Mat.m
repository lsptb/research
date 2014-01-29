function [S,excluded] = StructF2Mat(S,fieldname,fieldpath)

if ~isempty(fieldpath{1})
    for k = 1:length(fieldpath)
        Sin1 = structfun(@(x) ~isfield(x,fieldpath{k}),S);
        F1 = fieldnames(S);
        S = rmfield(S,F1(Sin1));
        S = structfun(@(x) x.(fieldpath{k}),S,'Uni',0);
    end
else F1 = []; Sin1 = [];
end
Sin2 = structfun(@(x) ~isfield(x,fieldname),S);
F2 = fieldnames(S);
S = rmfield(S,F2(Sin2));
if any(structfun(@(x) length(x.(fieldname)),S))
    S = structfun(@(x) {x.(fieldname)},S);
else S = structfun(@(x) x.(fieldname),S);
end
excluded = cat(1,F1(Sin1),F2(Sin2));


end

