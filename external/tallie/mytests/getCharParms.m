function parms = getCharParms(source,type)
% purpose: get parameters for cell characterization functions
% inputs:
% source: 
%   if sim: model spec structure
%   if exp: parmfile
% type: 'sim' or 'exp'
% see also: CharHyperpolStepTA, CharDepolStepTA, CharDepolTonicSpikesTA
% created on 13-Feb-2014 by JSS

parms = [];
switch type
  case 'sim' % extract parameters from model specification
    if isstruct(source)
      parms = parsespec(source);
    else
      error('parm source for simulated data must be the model specification structure.');
    end
  case 'exp' % extract parameters from excel spreadsheet
    if ischar(source) && exist(source,'file')
      parms = parsexls(source);
    else
      error('parm source for experimental data must be an excel filename');
    end
  otherwise
    error('unknown source type; must be ''sim'' or ''exp''.');
end

function parms = parsespec(spec)


function parms = parsexls(file)
% use Tallie's function for loading parameters
parms = GetParmsTA(file);


