function s=parseops(s,operations)
%PARSEOPS Parse and execute operations section of a DBAT XML file.
%
%    S=PARSEOPS(S,OPS) parses the operations XML block OPS from a DBAT
%    XML script file and executes the operations with the DBAT
%    structure S. The result is returned in the updated DBAT structure S.
%
%See also: PARSEINPUT, RUNDBATSCRIPT.
    
% Known operations
opsFields={'operation'};

% All fields are optional.
[ok,msg]=checkxmlfields(operations,opsFields,false(size(opsFields)));
if ~ok, error('DBAT XML script operations error: %s',msg); end

% Parse and execute each operations, in sequence.
if ~isfield(operations,'operation')
    ops={};
else
    ops=operations.operation;
    if ~iscell(ops)
        % Single operation
        ops={ops};
    end
end

% For each operation...
for i=1:length(ops)
    op=ops{i};
    if isfield(op,'Text')
        % Operation without structure
        switch op.Text
          case 'check_ray_count'
            warning('%s not implemented yet',op.Text)
          case 'spatial_resection'
            warning('%s not implemented yet',op.Text)
          case 'forward_intersection'
            warning('%s not implemented yet',op.Text)
          case 'bundle_adjustment'
            warning('%s not implemented yet',op.Text)
          otherwise
            error('DBAT XML script operations error: Unknown operation %s',...
                  op.Text)
        end
    else
        % Structured operation
        fn=fieldnames(op);
        if length(fn)>1
            joined=join(fn,', ');
            error(['DBAT XML script operations error: Too many fields ' ...
                   'in one operation: %s'],joined{1});
        end
        if iscell(op.(fn{1}))
            error(['DBAT XML script operations error: Repeated field ' ...
                   'in one operation: %s'],fn{1});
        end
        switch fn{1}
          case 'set_initial_values'
            s=parsesetinitialvalues(s,op.(fn{1}));
          case 'set_bundle_estimate_params'
            s=parsesetbundleest(s,op.(fn{1}));
          otherwise
            error('DBAT XML script operations error: Unknown operation %s',...
                  fn{1})
        end
    end
end
