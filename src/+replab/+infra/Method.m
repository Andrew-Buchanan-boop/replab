classdef Method < replab.infra.ClassElement
% Describes a method in a class
    
    properties
        declaration % charstring: Declaration line of the method
    end
    
    methods
        
        function self = Method(name, attributes, declaration, doc, packageNameParts, className, lineNumber)
            self.name = name;
            self.attributes = attributes;
            self.declaration = declaration;
            self.doc = doc;
            self.packageNameParts = packageNameParts;
            self.className = className;
            self.kind = 'method';
            self.lineNumber = lineNumber;
        end
        
        function str = headerStr(self)
            keywords = {};
            if isfield(self.attributes, 'Abstract') && self.attributes.Abstract
                keywords{1,end+1} = 'abstract';
            end
            if isfield(self.attributes, 'Access') && ismember(self.attributes.Access, {'private', 'protected'})
                keywords{1,end+1} = self.attributes.Access;
            end
            if isfield(self.attributes, 'Static') && self.attributes.Static
                keywords{1,end+1} = 'static';
            end
            keywords{1,end+1} = 'method';
            str = [self.name ' (' strjoin(keywords, ' ') ')'];
        end

    end
    
end
