classdef GurlsOptions < dynamicprops

    properties
    end
    
    methods
        function obj = GurlsOptions(oldOpt)
            if nargin == 1
                obj.newprops(oldOpt);
            end
        end
        
        function newprops(obj, oldOpt)
            if isa(oldOpt, 'struct')
                names = fieldnames(oldOpt);
            elseif isa(oldOpt, 'GurlsOptions')
                names = properties(oldOpt);
            else
                error('Unrecognizable GURLS option list');
            end
            for name=names'
                obj.newprop(name{1}, oldOpt.(name{1}));
            end
        end
        
        function newprop(obj, name, value)
            if numel(strfind(name, '.')) > 0
                C = strsplit(name, '.');
                if ~isprop(obj, C{1}) || ~isstruct(obj.(C{1}))
                    obj.newprop(C{1}, struct());
                end
                eval(['obj.',name,' = value;']);
            else
                if ~isprop(obj, name)
                    addprop(obj, name);
                end
                obj.(name) = value;
            end
        end

        function s = toStruct(obj)
	    s = struct();
            names = properties(obj);
            for name=names'
                s.(name{1}) = obj.(name{1});
            end
        end
    end 
end

