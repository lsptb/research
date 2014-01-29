function parms   = args2parms(options, options_schema, strict)
% dummy function while we migrate to mmil_args2parms

  if (~exist('strict', 'var'))
    parms = mmil_args2parms(options, options_schema);
  else
    parms = mmil_args2parms(options, options_schema, strict);
  end;
