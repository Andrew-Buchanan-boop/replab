classdef Function < replab.infra.PackageElement
% Describes a MATLAB function
%
% Is also used to parse methods, thus the `parseAbstractBody` parsing node
    
    properties
        declaration
    end
    
    methods
        
        function self = Function(name, declaration, doc, packageNameParts, fullFilename)
            self.name = name;
            self.declaration = declaration;
            self.doc = doc;
            self.packageNameParts = packageNameParts;
            self.kind = 'function';
            self.fullFilename = fullFilename;
        end
        
    end
    
    methods (Static)
                
        function f = fromParseState(ct, packageNameParts, filename)
        % Parses a function and returns a `replab.infra.Function` instance
            pos = 1;
            [pos fld] = replab.infra.FunctionLikeData.parse(ct, pos, []);
            assert(~isempty(pos));
            doc = replab.infra.Doc.leftTrimmed(fld.docLines);
            f = replab.infra.Function(fld.name, fld.declaration, doc, packageNameParts, filename);
        end
        
    end
    
end
