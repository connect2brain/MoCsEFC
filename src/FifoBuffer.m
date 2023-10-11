classdef FifoBuffer < handle
    % FIFOBUFFER is a first-in-first-out buffer, containing at most 
    % <length> numeric elements of size <elementSize>.
    % use ADD to put elements into the buffer, and GETALL to retrieve the
    % contents of the buffer (ordered from oldest to newest)
    
    properties(Access = private)
        data
        isSet
        next 
        entryDimensionIndex
    end
    
    properties(SetAccess = private)
        ElementSize
        Length
    end
    
    methods
        function this = FifoBuffer(bufferLength, elementSize)
            % FIFOBUFFER is a buffer, containing at most <length> numeric 
            % elements of size <elementSize>.
            if bufferLength <= 0 || floor(bufferLength) ~= bufferLength || isinf(bufferLength)
                fprintf('\nErrored\n')
                error('FifoBuffer:InputMustBeNatural', ...
                      'Length must be a positive integer')
            end
            this.Length = bufferLength;
            this.ElementSize = elementSize;
            this.data = nan([this.Length, this.ElementSize]);
            this.isSet = false(1, this.Length);
            
            this.entryDimensionIndex = repmat({':'}, [1, length(elementSize)]);
            
            this.next = 1;
        end
        
        function add(this, entry)
            % Add the given entry (of ElementSize) to the buffer.
            % If the buffer is already full, overwrites the oldest element
            this.data(this.next, this.entryDimensionIndex{:}) = entry;
            this.isSet(this.next) = true;
            
            this.next = rem(this.next, this.Length) + 1;
        end
        
        function clearAll(this)
            % Removes all elements from this buffer. ElementSize and length
            % remain unchanged.

            this.next = 1;
            this.isSet(:) = false;
            % Note: do not need to delete the actual contents of the
            % data-array, just declare them unset
        end
        
        function whether = isEmpty(this)
            % Returns whether this buffer is empty (no elements added)
            whether = all(~this.isSet);
        end
        
        function n = numberOfEntries(this)
            % Returns the number of entries in the buffer. This is between
            % 0 (empty) and <Length> (full)
            n = sum(this.isSet);
        end
        
        function whether = isFull(this)
            % Returns whether this buffer is full. When adding to a full
            % buffer, the oldest element in the buffer is replaced by the
            % new one.
            whether = all(this.isSet);
        end
        
        
        
        function [entries] = getAll(this)
            % Returns all elements that are currently stored in the buffer
            % Result has size <numberOfElements> x <ElementSize>
            entryIndices = [this.next:this.Length 1:(this.next - 1)];
            entryIndices = entryIndices(this.isSet(entryIndices));
            entries = this.data(entryIndices, this.entryDimensionIndex{:});
        end
    end
end

