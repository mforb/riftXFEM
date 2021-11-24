classdef basic < matlab.unittest.TestCase
    methods(Test)
        function alwaysworks(testCase)
            actSolution = [2 1];
            expSolution = [2 1];
            testCase.verifyEqual(actSolution,expSolution)
        end
        function realSolution(testCase)
          actSolution = quadraticSolver(1,-3,2);
          expSolution = [ 2, 1];
          testCase.verifyEqual(actSolution,expSolution)
        end
        function nonnumericInput(testCase)
          testCase.verifyError(@()quadraticSolver(1,'-3',2), ...
            'quadraticSolver:InputMustBeNumeric')
        end
    end
end

