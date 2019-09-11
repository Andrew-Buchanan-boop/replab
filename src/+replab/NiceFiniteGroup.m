classdef NiceFiniteGroup < replab.FiniteGroup
% A nice finite group is a finite group equipped with an injective homomorphism into a permutation group
%
% The injective homomorphism is called a 'nice monomorphism' to keep up with the GAP System notation.
%
% This nice monomorphism enables the computation of a BSGS chain, and thus to transfer group operations
% to operations on permutation groups.
%
% If this group is its own parent, the methods that are delegated to the parent group
% (including `eqv`/`compose`/`inverse`) needs to be overriden.
    
    properties (SetAccess = protected)
        parent % replab.NiceFiniteGroup: Parent nice finite group
        niceMonomorphism % function_handle: Injective group homomorphism from this group into a permutation group
    end
    
    
    properties (Access = protected)
        chain_ % replab.bsgs.Chain: BSGS chain describing this group
    end

    methods % Abstract
        
        function sub = subgroup(self, generators, order)
        % Constructs a subgroup of the current group from generators
        %
        % Args:
        %   generators (row cell array of elements of this group): List of generators
        %   order (vpi, optional): Subgroup order
        %
        % Returns:
        %   replab.NiceFiniteSubgroup: The subgroup generated by `generators`
            error('Abstract');
        end

    end
    
    methods (Access = protected)

        %% FiniteGroup methods
        
        function order = computeOrder(self)
            order = self.chain.order;
        end
        
        function chain = computeChain(self)
            imgId = self.niceMonomorphism(self.identity);
            n = length(imgId);
            nG = self.nGenerators;
            S = zeros(n, nG);
            for i = 1:nG
                S(:,i) = self.niceMonomorphism(self.generator(i));
            end
            chain = replab.bsgs.Chain.makeWithImages(n, S, self, self.generators);
        end

        function E = computeElements(self)
            E = replab.IndexedFamily.lambda(self.order, ...
                                            @(ind) self.chain.elementFromIndex(ind), ...
                                            @(el) self.chain.indexFromElement(el));
        end

        function dec = computeDecomposition(self)
            dec = replab.FiniteGroupDecomposition(self, self.chain.imagesDecomposition);
        end

    end

    methods

        function o = elementOrder(self, g)
        % Returns the order of a group element
        %
        % Args:
        %   g (element): Group element
        %
        % Returns:
        %   integer: The order of `g`, i.e. the smallest `o` such that ``g^o == identity``
            p = self.niceMonomorphism(g);
            orbits = replab.Partition.permutationsOrbits(p);
            orders = unique(orbits.blockSizes);
            o = 1;
            for i = 1:length(orders)
                o = lcm(o, orders(i));
            end
        end
        
        function c = chain(self)
        % Returns the BSGS chain corresponding to this group
        %
        % Returns:
        %   replab.bsgs.Chain: BSGS chain describing this group
            if isempty(self.chain_)
                self.chain_ = self.computeChain;
            end
            c = self.chain_;
        end

        %% Domain methods
        
        function b = eqv(self, x, y)
            b = self.parent.eqv(x, y);
        end
        
        %% Monoid methods
        
        function z = compose(self, x, y)
            z = self.parent.compose(x, y);
        end
 
        %% Group methods
        
        function xInv = inverse(self, x)
            xInv = self.parent.inverse(x);
        end

        %% CompactGroup methods
        
        function g = sampleUniformly(self)
            g = self.chain.sampleUniformly;
        end

        %% Methods enabled by the BSGS algorithms
        
        function b = contains(self, g)
        % Tests whether this group contains the given parent group element
        %
        % Abstract in `replab.FiniteSubgroup`
        %
        % Args:
        %   g (element of `self.parent`): Element to test membership of
        %
        % Returns:
        %   logical: True if this group contains `g` and false otherwise
            b = self.chain.contains(g);
        end
                               
        %% Representation construction
        
        function rho = rep(self, field, dimension, images)
        % Constructs a finite dimensional representation of this group
        %
        % Args:
        %   field ({'R', 'C'}): Whether the representation is real (R) or complex (C)
        %   dimension (integer): Representation dimension
        %   images (row cell array of matrices): Orthonormal/unitary images of the group generators
        %
        % Returns:
        %   replab.Rep: The constructed group representation
            rho = replab.RepByImages(self, field, dimension, images);
        end

        function rho = permutationRep(self, dimension, permutations)
        % Constructs a permutation representation of this group
        %
        % The returned representation is real. Use ``rep.complexification``
        % to obtain a complex representation.
        %
        % Args:
        %   dimension: Dimension of the representation
        %   permutations (row cell array of permutations): Images of the generators as permutations of size "dimension"
        %
        % Returns:
        %   replab.Rep: The constructed group representation
            S = replab.Permutations(dimension);
            f = @(g) S.toMatrix(g);
            images = cellfun(f, permutations, 'uniform', 0);
            rho = self.rep('R', dimension, images);
        end

        function rho = signedPermutationRep(self, dimension, signedPermutations)
        % Returns a real signed permutation representation of this group
        %
        % The returned representation is real. Use ``rep.complexification``
        % to obtain a complex representation.
        %
        % Args:
        %   dimension: Dimension of the representation
        %   signedPermutations (row cell array of signed permutations): Images of the generators as signed permutations of size "dimension"
        %
        % Returns:
        %   replab.Rep: The constructed group representation
            S = replab.SignedPermutations(dimension);
            f = @(g) S.toMatrix(g);
            images = cellfun(f, signedPermutations, 'uniform', 0);
            rho = self.rep('R', dimension, images);
        end

    end

end
