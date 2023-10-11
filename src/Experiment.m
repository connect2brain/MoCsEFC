classdef Experiment < handle
    %EXPERIMENT takes care of picking an order of high vs low condition
    %trials
    
    properties(SetAccess=private)
        Conditions
        CurrentTrial
    end
    
    properties(Access=private)
        plvBuffer
        q25
        q75
    end
    
    methods
        function this = Experiment(nTrialsHigh, nTrialsLow, plvBufferSize)
            %EXPERIMENT Construct an instance of this class
            %   Represents a MoCsEFC Experiment with <n_trials_high>
            %   high-functional-connectivity trials and <n_trials_low>
            %   low-FC trials in a random order
            %   <plvBufferSize> is the number of most recent PLVs that are
            %   stored for computing the adaptive PLV-criterion
            target_high = [false(1, nTrialsLow) true(1, nTrialsHigh)];
            this.Conditions = target_high(randperm(nTrialsHigh + nTrialsLow));
            this.CurrentTrial = 1;
            this.plvBuffer = FifoBuffer(plvBufferSize, 1);
            this.q25 = -Inf;
            this.q75 = Inf;
        end
        
        function whether = isDone(this)
            % ISDONE returns whether the experiment has concluded
            whether = this.CurrentTrial > length(this.Conditions);
        end
        
        function storePLV(this, plv)
            % STOREPLV adds given plv to the stored ones, and recomputes
            % the q25/q75 criteria. PLV can also be a list, which will be
            % added from start to end.
            for p = plv
                this.plvBuffer.add(p);
            end
            this.q25 = quantile(this.plvBuffer.getAll(), 0.25);
            this.q75 = quantile(this.plvBuffer.getAll(), 0.75);
        end
        
        function whether = fire(this, plv)
            % FIRE checks whether given PLV fits the next condition (high
            % xor low). 
            % If the PLV should be added for adaptive thresholding, call
            % storePLV separately!
            if this.isDone()
                whether = false; % don't fire if done.
            elseif this.Conditions(this.CurrentTrial) % High
                whether = plv >= this.q75;
            else % Low
                whether = plv <= this.q25;
            end
        end
        
        function [done, isHigh] = next(this)
            % NEXT steps the trial counter up (do this only when the coil
            % has actually sent a pulse)
            % Returns whether the experiment is done afterwards (done) and
            % whether the next trial (if existant) is a high-FC-trial
            isHigh = false;
            if this.CurrentTrial < length(this.Conditions)
                this.CurrentTrial = this.CurrentTrial + 1;
                isHigh = this.Conditions(this.CurrentTrial);
            elseif this.CurrentTrial == length(this.Conditions)
                this.CurrentTrial = this.CurrentTrial + 1;
            end
            done = this.isDone();
        end
        
        function [] = log(this, conditionFile, quantileFile, plvFile)
            % LOG writes infos about this instance to specified files.
            % Giving empty (i.e. []) for any of the files is possible, to
            % ignore the respective file
            if ~isempty(conditionFile)
                outstr = sprintf('%d,', this.Conditions);
                fprintf(conditionFile, [outstr(1:end-1) '\n']);
            end
            if ~isempty(quantileFile)
                fprintf(quantileFile, '%d,%d\n', this.q25, this.q75);
            end
            if ~isempty(plvFile)
                outstr = sprintf('%d,', this.plvBuffer.getAll);
                fprintf(plvFile, [outstr(1:end-1) '\n']);
            end
        end
        
        function [] = logEvent(this, eventFile, event)
            % LOGEVENT writes the given event as a line into given file
            % Use this in parallel with log to freely annotate the lines in
            % the other files. Make sure to keep the files consistent.
            fprintf(eventFile, sprintf('%s\n', string(event)));
        end
        
        function [plvs] = getAllPLVs(this)
            % GETALLPLVS returns all PLVs currently used for the criteria
            plvs = this.plvBuffer.getAll();
        end
    end
end

