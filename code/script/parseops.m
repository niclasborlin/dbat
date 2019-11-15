function s=parseops(s,operations)
%PARSEOPS Parse and execute operations section of a DBAT XML file.
%
%    S=PARSEOPS(S,OPS) parses the operations XML block OPS from a DBAT
%    XML script file and executes the operations with the DBAT
%    structure S. The result is returned in the updated DBAT structure S.
%
%See also: PARSEINPUT, RUNDBATSCRIPT.
    
% Known operations
opsFields={'operation','c'};

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
            cpId=s.OP.id(s.prior.OP.isCtrl);
            [s,~,fail]=resect(s,'all',cpId,1,0,cpId);
            if fail
                error('Resection failed.');
            end
          case 'forward_intersection'
            s=forwintersect(s,'all',true);
          case 'bundle_adjustment'
            s0=s;
            [s,ok,iters,sigma0,E]=bundle(s0);
            if ok
                fprintf(['Bundle ok after %d iterations with sigma0=%.2f ' ...
                         '(%.2f pixels)\n'],iters,sigma0,s.post.sigmas(1));
            else
                fprintf(['Bundle failed after %d iterations (code=%d). ' ...
                         'Last sigma0 estimate=%.2f (%.2f pixels)\' ...
                         'n'],iters,E.code,sigma0,sigma0*s0.IP.sigmas(1));
            end
            s.bundle.info=E;
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
