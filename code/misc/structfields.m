function f=structfields(s)
%STRUCTFIELDS Generate recursive list of fields in a structure.
%
%   STRUCTFIELDS(S) returns a cell array with the names of all
%   fields in S. Fields that are structures are handled
%   recursively.

f=cell(0,1);

if isempty(s) || ~isstruct(s)
    return;
end

fn=fieldnames(s);
for i=1:length(fn)
    field=s.(fn{i});
    if ~isstruct(field)
        f=cat(1,f,fn{i});
    else
        f2=structfields(field);
        f=cat(1,f,strcat(fn{i},'.',f2));
    end
end
