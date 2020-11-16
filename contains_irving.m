function tf = contains_irving(s, pattern,arg_input)
%CONTAINS True if pattern is found in text.
%   TF = CONTAINS_IRVING(STR,PATTERN) returns 1 (true) if STR contains PATTERN,
%   and returns 0 (false) otherwise.
%
%   TF = CONTAINS_V2(STR,PATTERN,'IGNORECASE') returns 1 (true) if STR
%   contains PATTERN, and returns 0 (false) otherwise. Case insensitive.
%
%   STR can be a string array, a character vector, or a cell array of
%   character vectors. So can PATTERN. PATTERN and STR need not be the same
%   size. If PATTERN is a string array or cell array, then CONTAINS returns
%   true if it finds any element of PATTERN in STR. If STR is a string
%   array or cell array, then TF is a logical array that is the same size.
% Description: created to replace MATLAB builtin "contains.m" in case older
% version of MATLAB.
% Author: Brita Irving <bkirving@alaska.edu> 
%%
narginchk(2, inf);
if nargin == 3 && strcmpi(arg_input,'ignorecase')
  s = lower(s);
  pattern = lower(pattern);
elseif nargin == 3 && ~strcmpi(arg_input,'ignorecase')
  fprintf('3rd input argument must be "ignorecase" but was %s\n',arg_input)
end
if iscell(s)
  tf = ~cellfun(@isempty,strfind(s,pattern));
else
  tf = ~isempty(strfind(s,pattern));
end
end