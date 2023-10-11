classdef TestExperiment < matlab.unittest.TestCase
    %TESTEXPERIMENT Unit-test-class for EXPERIMENT
    
    methods(TestClassSetup)
        function setupPath(testCase)
            addpath('../src')
        end
    end
    
    methods(TestClassTeardown)
        function cclearFiles(testCase)
            fileConditions = fopen('test-experiment-conditions.txt', 'w');
            fileCriteria = fopen('test-experiment-criteria.txt', 'w');
            filePLVs = fopen('test-experiment-plvs.txt', 'w');
            fileEvents = fopen('test-experiment-events.txt', 'w');

            fprintf(fileConditions, '');
            fprintf(fileCriteria, '');
            fprintf(filePLVs, '');
            fprintf(fileEvents, '');
            
            fclose(fileEvents);
            fclose(fileConditions);
            fclose(fileCriteria);
            fclose(filePLVs);
        end
    end
    
    methods(Test)
        function testCorrectNumberOfConditions(testCase)
            nHigh = 30;
            nLow = 50;
            ex = Experiment(nHigh, nLow, 10);
            testCase.verifyEqual(length(ex.Conditions), nHigh+nLow)
            testCase.verifyEqual(sum(ex.Conditions), nHigh)
            testCase.verifyEqual(sum(~ex.Conditions), nLow)
        end
        
        function testRandomOrderOfConditions(testCase)
            nHigh = 400;
            nLow  = 600;
            ex1 = Experiment(nHigh, nLow, 100);
            ex2 = Experiment(nHigh, nLow, 100);
            testCase.verifyFalse(all(ex1.Conditions == ex2.Conditions), ...
                'This test case has a very small chance of randomly failing. Rerun to be sure')
        end
        
        function testInitiallyTrialIsOne(testCase)
            ex = Experiment(10, 20, 301);
            testCase.verifyEqual(ex.CurrentTrial, 1);
        end
        
        function testNonEmptyIsNotDone(testCase)
            ex = Experiment(1, 0, 10);
            testCase.verifyFalse(ex.isDone);
            
            ex = Experiment(10, 230, 5);
            testCase.verifyFalse(ex.isDone);
        end
        
        function testEmptyIsImmediatelyDone(testCase)
            ex = Experiment(0, 0, 105);
            testCase.verifyTrue(ex.isDone);
        end
        
        function testPLVcriterion(testCase)
            initialPLVs = (2+randn(1, 200)*0.4).^2;
            lowPLV = quantile(initialPLVs, 0.1);
            highPLV = quantile(initialPLVs, 0.9);
            ex = Experiment(1, 1, 300);
            ex.storePLV(initialPLVs);
            if ex.Conditions(1) % i.e. high then low
                testCase.verifyTrue(ex.fire(highPLV));
                testCase.verifyFalse(ex.fire(lowPLV));
                
                [done, nextHigh] = next(ex);
                testCase.verifyFalse(done);
                testCase.verifyFalse(nextHigh);
                
                testCase.verifyFalse(ex.fire(highPLV));
                testCase.verifyTrue(ex.fire(lowPLV));
                
                [done, nextHigh] = next(ex);
                testCase.verifyTrue(done);
                testCase.verifyFalse(nextHigh);
            else % low then high
                testCase.verifyFalse(ex.fire(highPLV));
                testCase.verifyTrue(ex.fire(lowPLV));
                
                [done, nextHigh] = next(ex);
                testCase.verifyFalse(done);
                testCase.verifyTrue(nextHigh);
                
                testCase.verifyTrue(ex.fire(highPLV));
                testCase.verifyFalse(ex.fire(lowPLV));
                
                [done, nextHigh] = next(ex);
                testCase.verifyTrue(done);
                testCase.verifyFalse(nextHigh);
                % there is no next one, thus, it is not a high-FC trial
            end
        end
        
        function testNext(testCase)
            ex = Experiment(10, 12, 100);
            for i = 2:22
                [done, nextHigh] = next(ex);
                testCase.verifyEqual(nextHigh, ex.Conditions(i));
                testCase.verifyFalse(done);
            end
            [done, nextHigh] = next(ex);
            testCase.verifyTrue(done);
            testCase.verifyFalse(nextHigh);
        end
        
        function testLogAllGiven(testCase)
            fileConditions = fopen('test-experiment-conditions.txt', 'w');
            fileCriteria = fopen('test-experiment-criteria.txt', 'w');
            filePLVs = fopen('test-experiment-plvs.txt', 'w');
            
            ex = Experiment(10, 5, 12);
            plvs = 1:12;
            ex.storePLV(plvs);
            ex.log(fileConditions, fileCriteria, filePLVs);
            
            fclose(fileConditions);
            fclose(fileCriteria);
            fclose(filePLVs);
            
            expectedConditions = sprintf('%d,', ex.Conditions);
            expectedConditions = expectedConditions(1:end-1);
            loggedConditions = strtrim(fileread('test-experiment-conditions.txt'));
            testCase.verifyEqual(loggedConditions, expectedConditions);
            
            expectedPLVs = '1,2,3,4,5,6,7,8,9,10,11,12';
            loggedPLVs = strtrim(fileread('test-experiment-plvs.txt'));
            testCase.verifyEqual(loggedPLVs, expectedPLVs);
            
            expectedQuantiles = [quantile(1:12,0.25), quantile(1:12,0.75)];
            loggedQuantiles = cellfun(@(c) str2double(strtrim(c)), strsplit(strtrim(fileread('test-experiment-criteria.txt')), ','));
            testCase.verifyEqual(loggedQuantiles, expectedQuantiles);
            
            % Delete file contents after test.
            fileConditions = fopen('test-experiment-conditions.txt', 'w');
            fileCriteria = fopen('test-experiment-criteria.txt', 'w');
            filePLVs = fopen('test-experiment-plvs.txt', 'w');
            fprintf(fileConditions, '');
            fprintf(fileCriteria, '');
            fprintf(filePLVs, '');
            fclose(fileConditions);
            fclose(fileCriteria);
            fclose(filePLVs);
        end
        
        function testLog(testCase)
            fileCriteria = fopen('test-experiment-criteria.txt', 'w');
            ex = Experiment(10, 5, 12);
            plvs = 1:12;
            ex.storePLV(plvs);
            ex.log([], fileCriteria, []);
            ex.storePLV(13);
            ex.storePLV(14);
            ex.log([], fileCriteria, []);
            ex.storePLV([1 1 1 1]);
            ex.log([], fileCriteria, []);
            fclose(fileCriteria);
            
            expectedQuantiles = [quantile(1:12,0.25), quantile(1:12,0.75); 
                quantile(3:14,0.25), quantile(3:14,0.75);
                quantile([7:14 1 1 1 1],0.25), quantile([7:14 1 1 1 1],0.75)];
            lines = strsplit(strtrim(fileread('test-experiment-criteria.txt')), '\n');
            for i = 1:3
                loggedQuantiles = cellfun(@(c) str2double(strtrim(c)), strsplit(strtrim(lines{i}), ','));
                testCase.verifyEqual(loggedQuantiles, expectedQuantiles(i,:));
            end
        end
        
        function testLogEvent(testCase)
            fileEvents = fopen('test-experiment-events.txt', 'w');
            ex = Experiment(11, 15, 24);
            ex.logEvent(fileEvents, 'fired coil')
            ex.logEvent(fileEvents, 'waiting')
            ex.logEvent(fileEvents, '')
            ex.logEvent(fileEvents, 'fired coil, low condition')
            ex.logEvent(fileEvents, 'waiting')
            ex.logEvent(fileEvents, 35)
            fclose(fileEvents);
            
            lines = split(strtrim(fileread('test-experiment-events.txt')), newline);
            testCase.verifyEqual(lines{1}, 'fired coil');
            testCase.verifyEqual(lines{2}, 'waiting');
            testCase.verifyEqual(lines{3}, '');
            testCase.verifyEqual(lines{4}, 'fired coil, low condition');
            testCase.verifyEqual(lines{5}, 'waiting');
            testCase.verifyEqual(lines{6}, '35');
        end
        
        function testGetAllPLVsEmpty(testCase)
            ex = Experiment(9, 17, 14);
            plvs = ex.getAllPLVs();
            testCase.verifyTrue(isempty(plvs));
        end
        
        function testGetAllPLVs(testCase)
            ex = Experiment(9, 17, 14);
            ex.storePLV(1:10)
            plvs = ex.getAllPLVs();
            testCase.verifyEqual(plvs, (1:10)');
        end
    end
end

