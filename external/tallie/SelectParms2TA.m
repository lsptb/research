function S = SelectParms2TA(Parms,fieldname,val)

S = Parms;
if isnumeric(val)
    idc1 = arrayfun(@(x) find(x.(fieldname)==val),Parms,'Uni',0);
    idc2 = cellfun(@isempty,idc1);
    S(idc2) = [];
else
    idc1 = arrayfun(@(x) find(strcmp(val,x.(fieldname))),Parms,'Uni',0);
    idc2 = cellfun(@isempty,idc1);
    S(idc2) = [];
end

end

