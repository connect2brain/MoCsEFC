classdef TestFifoBuffer < matlab.unittest.TestCase
    %TESTFIFOBUFFER Unit-test-class for FIFOBUFFER
    
    methods(TestClassSetup)
        function setupPath(testCase)
            addpath('../src')
        end
    end
    
    methods(Test)
        function testValidLength(testCase)
            length = 14;
            elementsize = 3;
            buffer = FifoBuffer(length, elementsize);
            testCase.verifyEqual(buffer.Length, length, ...
                'Length is not correctly set')
        end
        
        function testInvalidLength(testCase)
            length = -2;
            elementsize = 7;
            verifyError(testCase, @() FifoBuffer(length, elementsize), ...
                'FifoBuffer:InputMustBeNatural')
        end
        
        function testElementSize1D(testCase)
            length = 26;
            elementsize = 8;
            buffer = FifoBuffer(length, elementsize);
            testCase.verifyEqual(buffer.ElementSize, elementsize, ...
                'One-dimensional elementsize is not correctly set')
        end
        
        function testElementSizeND(testCase)
            length = 26;
            elementsize = [3 12 7];
            buffer = FifoBuffer(length, elementsize);
            testCase.verifyEqual(buffer.ElementSize, elementsize, ...
                'N-dimensional elementsize is not correctly set')
        end
        
        function testEmptyInitially(testCase)
            buffer = FifoBuffer(153, [2 9]);
            testCase.verifyTrue(buffer.isEmpty, 'Buffer should start out empty')
        end
        
        function testEmptyThenAdd(testCase)
            buffer = FifoBuffer(291, 64);
            testCase.verifyTrue(buffer.isEmpty, 'Buffer should start out empty')
            buffer.add(ones(1,64));
            testCase.verifyFalse(buffer.isEmpty, 'Buffer should no longer be empty after adding an element')
        end
        
        function testEmptyAfterClear(testCase)
            buffer = FifoBuffer(291, 1);
            testCase.verifyTrue(buffer.isEmpty, 'Buffer should start out empty')
            buffer.add(3);
            buffer.add(8+0.5i);
            buffer.add(-3.12);
            testCase.verifyFalse(buffer.isEmpty, 'Buffer should no longer be empty after adding elements')
            buffer.clearAll();
            testCase.verifyTrue(buffer.isEmpty, 'Buffer should be empty again after calling clearAll')
        end
        
        function testGetAllOnEmpty(testCase)
            buffer = FifoBuffer(5, [3 2]);
            testCase.verifyEqual(buffer.getAll(), zeros(0, 3, 2), ...
                'GetAll should return an array whose first dimension has length zero, but whose remaining dimensions have elementsize')
        end
        
        function testGetAllNoOverwriting(testCase)
            buffer = FifoBuffer(15, [3 2]);
            entries = rand(4, 3, 2);
            buffer.add(squeeze(entries(1,:,:)));
            buffer.add(squeeze(entries(2,:,:)));
            buffer.add(squeeze(entries(3,:,:)));
            buffer.add(squeeze(entries(4,:,:)));
            
            testCase.verifyEqual(buffer.getAll(), entries, ...
                'Add should add the given data, and GetAll retrieve it and only it')
        end
        
        function testGetAllOverwriting(testCase)
            buffer = FifoBuffer(3, [2 4]);
            entries = rand(5, 2, 4);
            buffer.add(squeeze(entries(1,:,:)));
            buffer.add(squeeze(entries(2,:,:)));
            buffer.add(squeeze(entries(3,:,:)));
            buffer.add(squeeze(entries(4,:,:)));
            buffer.add(squeeze(entries(5,:,:)));
            
            testCase.verifyEqual(buffer.getAll(), entries(3:5,:,:), ...
                'Add should add the given data, and overwrite the oldest elements if the buffer is full, and GetAll retrieve it and only it')
        end
        
        function testGetAllDoubleOverwriting(testCase)
            buffer = FifoBuffer(3, [2 4]);
            entries = rand(8, 2, 4);
            buffer.add(squeeze(entries(1,:,:)));
            buffer.add(squeeze(entries(2,:,:)));
            buffer.add(squeeze(entries(3,:,:)));
            buffer.add(squeeze(entries(4,:,:))); % overwrites 1
            buffer.add(squeeze(entries(5,:,:))); % overwrites 2
            buffer.add(squeeze(entries(6,:,:))); % overwrites 3
            buffer.add(squeeze(entries(7,:,:))); % overwrites 4
            buffer.add(squeeze(entries(8,:,:))); % overwrites 5
            
            testCase.verifyEqual(buffer.getAll(), entries(6:8,:,:), ...
                'Add should add the given data, and overwrite the oldest elements if the buffer is full, and GetAll retrieve it and only it')
        end
        
        function testFull(testCase)
            buffer = FifoBuffer(3, 1);
            testCase.verifyFalse(buffer.isFull, 'Buffer should start out empty (not full)')
            buffer.add(0.5);
            testCase.verifyFalse(buffer.isFull, 'Buffer has length 3, should not be full after adding one element')
            buffer.add(nan);
            testCase.verifyFalse(buffer.isFull, 'Buffer has length 3, should not be full after adding two elements')
            buffer.add(-1);
            testCase.verifyTrue(buffer.isFull, 'Buffer has length 3, 3 elements were added, it should be full')
            buffer.add(3.1416);
            testCase.verifyTrue(buffer.isFull, 'Buffer stays full for adding more elements (overwriting)')
        end
        
        function testNumberOfEntries(testCase)
            buffer = FifoBuffer(3, [1 8 3]);
            testCase.verifyEqual(buffer.numberOfEntries, 0, 'Nothing has been added, so number of elements should be 0');
            buffer.add(rand(1, 8, 3));
            testCase.verifyEqual(buffer.numberOfEntries, 1, 'One element has been added, so number of elements should be 1');
            buffer.add(rand(1, 8, 3));
            testCase.verifyEqual(buffer.numberOfEntries, 2, '2 elements have been added, so number of elements should be 2');
            buffer.add(rand(1, 8, 3));
            testCase.verifyEqual(buffer.numberOfEntries, 3, '3 elements have been added, so number of elements should be 3');
            buffer.add(rand(1, 8, 3));
            testCase.verifyEqual(buffer.numberOfEntries, 3, '4 elements have been added, but the buffer has length 3, so number of elements should be 3');
        end
        
        function testGetAllConsistentWithNumberOfEntries(testCase)
            buffer = FifoBuffer(5, 9);
            entries = buffer.getAll();
            testCase.verifyEqual(size(entries, 1), buffer.numberOfEntries)
            
            buffer.add(rand(1,9));
            entries = buffer.getAll();
            testCase.verifyEqual(size(entries, 1), buffer.numberOfEntries)
            
            buffer.add(rand(1,9));
            buffer.add(rand(1,9));
            buffer.add(rand(1,9));
            buffer.add(rand(1,9));
            buffer.add(rand(1,9));
            buffer.add(rand(1,9));
            entries = buffer.getAll();
            testCase.verifyEqual(size(entries, 1), buffer.numberOfEntries)
        end
        
        function testGetAllConsistentWithElementSize(testCase)
            buffer = FifoBuffer(1, [3 7]);
            entries = buffer.getAll();
            testCase.verifyEqual(size(entries, 2:3), buffer.ElementSize)
            
            buffer.add(rand(3, 7))
            buffer.add(rand(3, 7))
            entries = buffer.getAll();
            testCase.verifyEqual(size(entries, 2:3), buffer.ElementSize)
        end
        
        function testClearAllOnEmpty(testCase)
            buffer = FifoBuffer(31, [5 3]);
            buffer.clearAll();
            testCase.verifyTrue(buffer.isEmpty, 'Buffer should be empty after calling clearAll')
            added = rand(1, 5, 3); % it is optional to have the first dimension explicitly
            buffer.add(added);
            testCase.verifyEqual(buffer.getAll(), added);
            testCase.verifyFalse(buffer.isEmpty, 'Elements should be added as usual after call to clearAll')
        end
        
        function testClearAllOnNonEmpty(testCase)
            buffer = FifoBuffer(4, 2);
            buffer.add(rand(1,2));
            buffer.add(rand(1,2));
            buffer.add(rand(1,2));
            buffer.clearAll();
            testCase.verifyEqual(buffer.getAll, zeros(0,2))
            testCase.verifyTrue(buffer.isEmpty)
            
            newOnes = rand(6,2);
            buffer.add(newOnes(1,:));
            buffer.add(newOnes(2,:));
            buffer.add(newOnes(3,:));
            buffer.add(newOnes(4,:));
            buffer.add(newOnes(5,:)); % overwrites 1
            buffer.add(newOnes(6,:)); % overwrites 2
            testCase.verifyEqual(buffer.getAll(), newOnes(3:6, :))
            testCase.verifyTrue(buffer.isFull)
        end
    end
    
end

