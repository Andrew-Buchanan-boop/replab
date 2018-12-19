classdef FiniteGroupRep < replab.Rep
    
    properties (SetAccess = protected)
        images;
        isUnitary;
    end
    
    methods
        
        function self = FiniteGroupRep(group, images, isUnitary, T)
            assert(length(images) == group.nGenerators);
            self.group = group;
            self.images = images;
            self.isUnitary = isUnitary;
            self.T = T;
            self.d = T.n;
        end
        
        function rho = image(self, g)
            word = self.group.factorization(g);
            rho = self.T.identity;
            for i = 1:length(word.indices)
                g = self.images{word.indices(i)};
                e = word.exponents(i);
                ge = self.T.composeN(g, e);
                rho = self.T.compose(rho, ge);
            end
        end
        
        function s = str(self)
            if self.isUnitary
                t = 'Unitary representation';
            else
                t = 'Representation';
            end
            s = sprintf('%s of dimension %d with generator images', t, self.d);
            for i = 1:length(self.images)
                gen = char('a' + i - 1);
                s = [s char(10) '- ' gen ':' char(10)];
                img = replab.prependLines(replab.strOf(self.images{i}), '    ');
                s = [s img char(10)];
            end
        end
        
        function disp(self)
            disp(self.str);
        end
        
    end
    
    properties (Access = protected)
        % Cached values, see the corresponding methods for help
        centralizerAlgebra_ = [];
        fibers_ = [];
        isotypic_ = [];
    end
    
    methods

        function I = isotypic(self)
            if isempty(self.isotypic_)
                self.isotypic_ = replab.IsotypicDecomposition.ofRep(self);
            end
            I = replab.IsotypicDecomposition.ofRep(self);
        end
        
        function a = centralizerAlgebra(self)
        % Returns the centralizer algebra of this representation,
        % which corresponds to d x d matrices that are invariant
        % under conjugation by this representation
            if isempty(self.centralizerAlgebra_)
                self.centralizerAlgebra_ = replab.rep.Algebra.forRep(self);
            end
            a = self.centralizerAlgebra_;
        end

        function f = fibers(self)
        % Returns the finest partition of 1..d into blocks corresponding
        % to subrepresentations of this representation.
        %
        % When the representation comes from a permutation group
        % this corresponds to the orbits.
            if isempty(self.fibers_)
                d = self.d;
                mask = false(d, d);
                for i = 1:length(self.images)
                    mask = mask | (self.images{i} ~= 0);
                end
                self.fibers_ = replab.Partition.connectedComponents(mask);
            end
            f = self.fibers_;
        end
        
% $$$                 function I = isoDec(self)
% $$$             if isempty(self.isoDec_)
% $$$                 self.isoDec_ = replab.rep.IsoDec.forAlgebra(self.centralizerAlgebra);
% $$$             end
% $$$             I = self.isoDec_;
% $$$         end
% $$$         
% $$$         function n = nIsotypicComponents(self)
% $$$         % Returns the number of isotypic components in this representation
% $$$             n = self.isoDec.nComponents;
% $$$         end
% $$$         
% $$$         function [subrho U dim mul] = isotypicComponent(self, i)
% $$$         % Returns the i-th isotypic component
% $$$         % Outputs:
% $$$         % subrho is the subrepresentation in this isotypic component
% $$$         % U      is the change of basis matrix such that subrho = U'*rho*U
% $$$         % dim    is the representation dimension
% $$$         % mul    is the representation multiplicity 
% $$$             U = self.isoDec.compBasis(i);
% $$$             subImages = cell(1, length(self.images));
% $$$             for i = 1:length(self.images)
% $$$                 subImages{i} = U'*self.images{i}*U;
% $$$             end
% $$$             md = size(U, 2);
% $$$             subT = replab.GeneralLinearGroup(md, 'R15'); %TODO
% $$$             subrho = replab.FiniteGroupRep(self.group, subImages, subT);
% $$$             dim = self.isoDec.repDims(i);
% $$$             mul = self.isoDec.repMuls(i);
% $$$         end
% $$$         
% $$$         function M = irrepMultiplicities(self)
% $$$         end
% $$$         
% $$$         function D = irrepDimensions(self)
% $$$         end
% $$$         
% $$$         function [subrho U] = irrep(self, i, j)
% $$$         end
% $$$         
        
        
        

        
% $$$         function M = sampleGroupElement(self)
% $$$             M = self.image(self.group.sample);
% $$$         end
% $$$ 
% $$$         function M = sampleGroupAlgebra(self, nSamples, nRounds)
% $$$             if nargin < 3
% $$$                 nRounds = 4;
% $$$              end
% $$$              if nargin < 2
% $$$                  nSamples = 5;
% $$$             end                
% $$$             M = self.T.identity;
% $$$             for i = 1:nRounds
% $$$                 M1 = randn * self.image(self.group.sample);
% $$$                 for j = 2:nSamples
% $$$                     M1 = M1 + randn * self.image(self.group.sample);
% $$$                 end
% $$$                 M = M * M1;
% $$$             end
% $$$             M = M / sqrt(nSamples^nRounds);
% $$$         end
        
    end

end