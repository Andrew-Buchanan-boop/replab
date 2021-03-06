classdef MorphismLaws < replab.Laws

    properties (SetAccess = protected)
        morphism % (`+replab.Morphism`): Morphism tested
        S % `.Group`: Source group
        T % `.Group`: Target group
    end

    methods

        function self = MorphismLaws(morphism)
            self.morphism = morphism;
            self.S = morphism.source;
            self.T = morphism.target;
        end

    end

    methods % LAWS

        function law_inverse_S(self, s)
            t = self.morphism.image(s);
            sI = self.S.inverse(s);
            tI1 = self.T.inverse(t);
            tI2 = self.morphism.image(sI);
            self.T.assertEqv(tI1, tI2);
        end

        function law_composition_SS(self, s1, s2)
            s12 = self.S.compose(s1, s2);
            t1 = self.morphism.image(s1);
            t2 = self.morphism.image(s2);
            t12_1 = self.morphism.image(s12);
            t12_2 = self.T.compose(t1, t2);
            self.T.assertEqv(t12_1, t12_2);
        end

        function law_identity_(self)
            self.T.assertEqv(self.T.identity, self.morphism.image(self.S.identity));
        end

    end

end
