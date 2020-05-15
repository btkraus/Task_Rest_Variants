function [mintasksamps] = mintasksamps(QCFiles_path, subs, SplitHalf, ConcatenateTasks, MatchData, MatchAcrossSubs)

    % sets up variables for sum of sample points in each task of even- and odd-numbered
    % sessions to be used later
    memptsoddsum = [];
    motorptsoddsum = [];
    mixedptsoddsum = [];
    memptsevensum = [];
    motorptsevensum = [];
    mixedptsevensum = [];
    
    %% Matches data points for all tasks across subjects
    if MatchData ==1 && MatchAcrossSubs ==1
        %determines the number of sample points to match
        for n = 1:numel(subs)
            load([QCFiles_path subs{n} '_QCFile.mat']);
            
            % sets up variables for number of sample points available for
            % current subject
            memptsodd = [];
            motorptsodd = [];
            mixedptsodd = [];
            memptseven = [];
            motorptseven = [];
            mixedptseven = [];
            
            % QC files contain SubStruct with # of sample points for each task
         % and whether the session is an even- or odd-numbered session
        for u = 1:length(SubStruct)
            
            % divides sample points between even and odd sessions
            % 1 = Odd Session, 2 = Even Session
            if SubStruct(u).OddEven == 1
                memptsodd = [memptsodd; SubStruct(u).MemSampPts];
                mixedptsodd = [mixedptsodd; SubStruct(u).MixedSampPts];
                if ~strcmp(subs{n}, 'MSC09') % different for MSC09 because not enough motor sample points
                    motorptsodd = [motorptsodd; SubStruct(u).MotorSampPts];
                end
            elseif SubStruct(u).OddEven == 2
                memptseven = [memptseven; SubStruct(u).MemSampPts];
                mixedptseven = [mixedptseven; SubStruct(u).MixedSampPts];
                if ~strcmp(subs{n}, 'MSC09')
                    motorptseven = [motorptseven; SubStruct(u).MotorSampPts];
                end
      
            end
        end
    
        memptsoddsum = [memptsoddsum; sum(memptsodd)];
        mixedptsoddsum = [mixedptsoddsum; sum(mixedptsodd)];
        memptsevensum = [memptsevensum; sum(memptseven)];
        mixedptsevensum = [mixedptsevensum; sum(mixedptseven)];
        
        % sums all sample points for all even- or odd-numbered sessions for
        % current subject
        if strcmp(subs{n}, 'MSC09')
            motorptsoddsum = [motorptsoddsum; 9999];
            motorptsevensum = [motorptsevensum; 9999];
        else
            motorptsoddsum = [motorptsoddsum; sum(motorptsodd)];
            motorptsevensum = [motorptsevensum; sum(motorptseven)];
        end

        
    end
    
    %finds minimum number of sample points to match the other tasks across
    %subjects
    if SplitHalf == 1
        mintasksamps = min(min([memptsoddsum motorptsoddsum mixedptsoddsum memptsevensum motorptsevensum mixedptsevensum]));
    elseif ConcatenateTasks == 1
        mintasksamps = min(min([(memptsoddsum + memptsevensum)  (motorptsoddsum + motorptsevensum)  (mixedptsoddsum + mixedptsevensum)]));
    end    
end
            
end