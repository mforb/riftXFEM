classdef XFEMTest < matlab.unittest.TestCase
    methods(Test)
        function alwaysworks(testCase)
            actSolution = [2 1];
            expSolution = [2 1];
            testCase.verifyEqualactSolution,expSolution)
        end
    end
end

