%% the goal of this function is to get the number of samples from each session of the MSC

% in Diana's function, she found the minimum number of samples present in
% any MSC subject for any particular task and split half. So its min(num
% samples existest for any task/split half/subject). 

% With this function, our goal is to get that number of of samples and take
% it as equally as possible from the sessions which exist. 
% it does not choose here to the minimum number that exist in any
% particular session for each one (there are probably some very bad ones).
% So, instead, we are going to go through each session, take what we think
% we need, and then if one doesn't have enough, hwe have to sample equally
% to make up for it from the rest. 

function [sampspersession] = GetSampsPerSession(mintasksamps, task, sub, QCFiles_path, SplitHalf, MatchData, ConcatenateTasks)

    if MatchData == 1 && SplitHalf == 1
        load ([QCFiles_path sub '_QCFile.mat']);

        memptsodd = [];
        motorptsodd = [];
        mixedptsodd = [];
        memptseven = [];
        motorptseven = [];
        mixedptseven = [];
      for u = 1:length(SubStruct)
                if SubStruct(u).OddEven == 1
                    memptsodd = [memptsodd; SubStruct(u).MemSampPts];
                    mixedptsodd = [mixedptsodd; SubStruct(u).MixedSampPts];
                    motorptsodd = [motorptsodd SubStruct(u).MotorSampPts];
                elseif SubStruct(u).OddEven == 2
                    memptseven = [memptseven; SubStruct(u).MemSampPts];
                    mixedptseven = [mixedptseven; SubStruct(u).MixedSampPts];
                    motorptseven = [motorptseven; SubStruct(u).MotorSampPts];
                end
      end

     % gets minimum number of sample points and then mean needed per sesh
    if strcmp(sub, 'MSC09')  
        minsampspts = min([sum(memptsodd) sum(mixedptsodd) sum(memptseven) sum(mixedptseven)]);
    else
        minsampspts = min([sum(memptsodd) sum(mixedptsodd) sum(motorptsodd) sum(memptseven) sum(mixedptseven) sum(motorptseven)]);
    end

    meansamppts = floor(mintasksamps/5);
    disp(sprintf('For all subjects, the minimum number of sample points in a split-half is %i, with a mean of %i points per sesssion: %s', mintasksamps, meansamppts, datestr(now)));
    meansampptstempodd = meansamppts;
    meansampptstempeven = meansamppts;

    % Sets data points for subject MSC09
    if strcmp(sub, 'MSC09') && strcmp(task, 'motor') && ConcatenateTasks == 1                      
        meansampptstempodd = floor(sum(motorptsodd)/5);
        meansampptstempeven = floor(sum(motorptseven)/5);                        
    elseif strcmp(sub, 'MSC09') && ConcatenateTasks == 1                        
        meansampptstempodd = meansamppts + round((meansamppts - round(sum(motorptsodd)/5))/2);
        meansampptstempeven = meansamppts + round((meansamppts - round(sum(motorptseven)/5))/2);                        
    end

    % get the data for all the other subjects
    remaindertotalodd = 0;
    remaindertotaleven = 0;           
    notenoughdata = zeros(length(SubStruct),1);

    while true                
        remainderodd = 0;
        remaindereven = 0;
        sampspersession = zeros(length(SubStruct),1);

        notenoughdataeven = sum(notenoughdata(2:2:end));
        notenoughdataodd = sum(notenoughdata(1:2:end));

        meansampptstempodd = (meansampptstempodd + round(remaindertotalodd/(5-notenoughdataodd)));
        meansampptstempeven = (meansampptstempeven + round(remaindertotaleven/(5-notenoughdataeven)));

        for v = 1:length(SubStruct)                        
            if strcmp(task, 'mem')                
                if SubStruct(v).OddEven == 1
                    if SubStruct(v).MemSampPts < meansampptstempodd && notenoughdata(v) == 0                   
                        notenoughdata(v) = 1;
                        remainderodd = (remainderodd + (meansampptstempodd - SubStruct(v).MemSampPts));
                        sampspersession(v) = SubStruct(v).MemSampPts;                    
                        disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                    
                    elseif notenoughdata(v) == 0                   
                        sampspersession(v) = meansampptstempodd;                    
                        disp(sprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                               
                    else                                
                        sampspersession(v) = SubStruct(v).MemSampPts;                                
                        disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                    
                    end                    
                else                    
                    if SubStruct(v).MemSampPts < meansampptstempeven && notenoughdata(v) == 0                    
                        notenoughdata(v) = 1;
                        remaindereven = (remaindereven + (meansampptstempeven - SubStruct(v).MemSampPts));
                        sampspersession(v) = SubStruct(v).MemSampPts;                    
                        disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                    
                    elseif notenoughdata(v) == 0                    
                        sampspersession(v) = meansampptstempeven;                    
                        disp(sprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                                
                    else                                
                        sampspersession(v) = SubStruct(v).MemSampPts;                                
                        disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                    
                    end                    
                end                            
            elseif strcmp(task, 'mixed')                            
                 if SubStruct(v).OddEven == 1
                    if SubStruct(v).MixedSampPts < meansampptstempodd && notenoughdata(v) == 0                    
                        notenoughdata(v) = 1;
                        remainderodd = (remainderodd + (meansampptstempodd - SubStruct(v).MixedSampPts));
                        sampspersession(v) = SubStruct(v).MixedSampPts;                    
                        disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                    
                    elseif notenoughdata(v) == 0                    
                        sampspersession(v) = meansampptstempodd;                    
                        disp(sprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                               
                    else                                
                        sampspersession(v) = SubStruct(v).MixedSampPts;                                
                        disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                    
                    end                    
                 else                    
                    if SubStruct(v).MixedSampPts < meansampptstempeven && notenoughdata(v) == 0                    
                        notenoughdata(v) = 1;
                        remaindereven = (remaindereven + (meansampptstempeven - SubStruct(v).MixedSampPts));
                        sampspersession(v) = SubStruct(v).MixedSampPts;                    
                        disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                    
                    elseif notenoughdata(v) == 0                    
                        sampspersession(v) = meansampptstempeven;                    
                        disp(sprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                                
                    else                                
                        sampspersession(v) = SubStruct(v).MixedSampPts;                                
                        disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                    
                    end
                 end                            
            elseif strcmp(task, 'motor')                            
                if SubStruct(v).OddEven == 1
                    if SubStruct(v).MotorSampPts < meansampptstempodd && notenoughdata(v) == 0                    
                        notenoughdata(v) = 1;
                        remainderodd = (remainderodd + (meansampptstempodd - SubStruct(v).MotorSampPts));
                        sampspersession(v) = SubStruct(v).MotorSampPts;                    
                        disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                    
                    elseif notenoughdata(v) == 0                    
                        sampspersession(v) = meansampptstempodd;                    
                        disp(sprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                                
                    else                                
                        sampspersession(v) = SubStruct(v).MotorSampPts;                                
                        disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                    
                    end                    
                else                    
                    if SubStruct(v).MotorSampPts < meansampptstempeven && notenoughdata(v) == 0                    
                        notenoughdata(v) = 1;
                        remaindereven = (remaindereven + (meansampptstempeven - SubStruct(v).MotorSampPts));
                        sampspersession(v) = SubStruct(v).MotorSampPts;                    
                        disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                    
                    elseif notenoughdata(v) == 0                    
                        sampspersession(v) = meansampptstempeven;                    
                        disp(sprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                                
                    else                                
                        sampspersession(v) = SubStruct(v).MotorSampPts;                                
                        disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                    
                    end                            
                end
            end    
        end
        if remainderodd == 0 && remaindereven == 0                    
            disp(sprintf('Data sampling calculated for subject %s for task %s: %s', sub, task, datestr(now)));                
            break                    
        elseif (sum(notenoughdata) == 10) || (sum(notenoughdataeven) == 5 && remainderodd == 0) || (sum(notenoughdataodd) == 5 && remaindereven == 0)                    
            disp(sprintf('Not enough data is present for this task to have enough data for subject %s for task %s: %s', sub, task, datestr(now)));                    
            break                
        else                
            remaindertotalodd = remainderodd;
            remaindertotaleven = remaindereven;                    
            disp(sprintf('Remainder of %i calculated for odd files for subject %s for task %s and session %s: %s', remainderodd, sub, task, SubStruct(v).SessionFile, datestr(now)));
            disp(sprintf('Remainder of %i calculated for even files for subject %s for task %s and session %s: %s', remaindereven, sub, task, SubStruct(v).SessionFile, datestr(now)));                
        end
    end

    elseif MatchData == 1         
        load ([QCFiles_path sub '_QCFile.mat']);

        mempts = [];
        motorpts = [];
        mixedpts = [];

        % get number of data points from QC files for each task
        for u = 1:length(SubStruct)                
            mempts = [mempts; SubStruct(u).MemSampPts];
            mixedpts = [mixedpts; SubStruct(u).MixedSampPts];
            motorpts = [motorpts; SubStruct(u).MotorSampPts];
        end

        minsampspts = min([sum(mempts) sum(mixedpts) sum(motorpts)]);

         % gets mean number of data points for current subjects
        if MatchAcrossSubs == 0
            meansamppts = floor(minsampspts/5);                
            disp(sprintf('For subject %s the minimum number of sample points in a split-half is %i, with a mean of %i points per sesssion: %s', sub, minsampspts, meansamppts, datestr(now)));            
        elseif SplitHalf == 1                
            meansamppts = floor(mintasksamps/5);        
            disp(sprintf('For all subjects, the minimum number of sample points in a split-half is %i, with a mean of %i points per sesssion: %s', mintasksamps, meansamppts, datestr(now)));            
            meansampptstempodd = meansamppts;
            meansampptstempeven = meansamppts;                
        else                
            meansamppts = floor(mintasksamps/10);
            meansampptstemp = meansamppts;                
            disp(sprintf('For all subjects, the minimum number of sample points in a full file is %i, with a mean of %i points per sesssion: %s', mintasksamps, meansamppts, datestr(now)));            
        end

        disp(sprintf('For subject %s the minimum number of sample points is %i, with a mean of %i points per sesssion per task for concatenation: %s', sub, minsampspts, meansamppts/3, datestr(now)));

        % sets data points for MSC09
        if strcmp(sub, 'MSC09') && strcmp(task, 'motor') && ConcatenateTasks == 1  %% Sets motor task data points for subject MSC09                
            if SplitHalf == 1                   
                meansampptstempodd = floor(sum(motorptsodd)/5);
                meansampptstempeven = floor(sum(motorptseven)/5);                    
            else                    
                meansampptstemp = floor(sum(motorpts)/10);                    
            end                        
        elseif strcmp(sub, 'MSC09') && ConcatenateTasks == 1                
            if SplitHalf == 1                        
                meansampptstempodd = meansamppts + round((meansamppts - round(sum(motorptsodd)/5))/2);
                meansampptstempeven = meansamppts + round((meansamppts - round(sum(motorptseven)/5))/2);                    
            else                    
                meansampptstemp = meansamppts + round((meansamppts - round(sum(motorpts)/10))/2);                    
            end                        
        end

        remaindertotal = 0;           
        notenoughdata = zeros(length(SubStruct),1);

        % if not enough data, will add sample points to match other
        % tasks/subjects
        while true                
            remainder = 0;
            sampspersession = zeros(length(SubStruct),1);                
            notenoughdatatotal = sum(notenoughdata);            
            meansampptstemp = (meansampptstemp + round(remaindertotal/(10-notenoughdatatotal)));

            for v = 1:length(SubStruct)                        
                if strcmp(task, 'mem')                
                    if SubStruct(v).MemSampPts < meansampptstemp && notenoughdata(v) == 0                    
                        notenoughdata(v) = 1;
                        remainder = (remainder + (meansampptstemp - SubStruct(v).MemSampPts));
                        sampspersession(v) = SubStruct(v).MemSampPts;                    
                        disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                    
                    elseif notenoughdata(v) == 0                    
                        sampspersession(v) = meansampptstemp;                    
                        disp(sprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                                
                    else                                
                        sampspersession(v) = SubStruct(v).MemSampPts;                                
                        disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                    
                    end                                           
                elseif strcmp(task, 'mixed')                            
                    if SubStruct(v).MixedSampPts < meansampptstemp && notenoughdata(v) == 0                    
                        notenoughdata(v) = 1;
                        remainder = (remainder + (meansampptstemp - SubStruct(v).MixedSampPts));
                        sampspersession(v) = SubStruct(v).MixedSampPts;                    
                        disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                    
                    elseif notenoughdata(v) == 0                    
                        sampspersession(v) = meansampptstemp;                   
                        disp(sprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                                
                    else                                
                        sampspersession(v) = SubStruct(v).MixedSampPts;                                
                        disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                    
                    end                                          
                elseif strcmp(task, 'motor')                            
                    if SubStruct(v).MotorSampPts < meansampptstemp && notenoughdata(v) == 0                    
                        notenoughdata(v) = 1;
                        remainder = (remainder + (meansampptstemp - SubStruct(v).MotorSampPts));
                        sampspersession(v) = SubStruct(v).MotorSampPts;                    
                        disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                    
                    elseif notenoughdata(v) == 0                    
                        sampspersession(v) = meansampptstemp;                    
                        disp(sprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                                
                    else                                
                        sampspersession(v) = SubStruct(v).MotorSampPts;                                
                        disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));                    
                    end   
                end    
            end

            % if needed, samples data from other (sessions?) to match
            % number of sample points for other tasks/sessions
            if remainder == 0                    
                disp(sprintf('Data sampling calculated for subject %s for task %s: %s', sub, task, datestr(now)));                
                break                    
            elseif (sum(notenoughdata) == 10) || (sum(notenoughdata) == 10 && remainder == 0)                    
                disp(sprintf('Not enough data is present for this task to have enough data for subject %s for task %s: %s', sub, task, datestr(now)));                    
                break                
            else                
                remaindertotal = remainder;                    
                disp(sprintf('Remainder of %i calculated for subject %s for task %s and session %s: %s', remaindertotal, sub, task, SubStruct(v).SessionFile, datestr(now)));                
            end
        end            
    end
end



