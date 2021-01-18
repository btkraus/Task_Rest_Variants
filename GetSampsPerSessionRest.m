
%% Main for-loop: matches number of data points to sample from rest to task, makes Dconns
function [sampspersession, restsampspersession, tasksampspersession] = GetSampsPerSessionRest(mintasksamps, MatchAcrossTasks, MatchAcrossSubs, ConcatenateTaskData, ConcatenateTasksMatch, QC_path, SplitHalf, MatchData, tasks, sub)
%sub = 'MSC06';
restsampspersession= [];
sampspersession= [];
tasksampspersession = [];

if ConcatenateTasksMatch == 1    
    taskbackup = tasks;
end
    
    % establishing number of sample points to match for each session
    if ConcatenateTaskData == 1
        if SplitHalf == 1
            sampspersession = floor(repmat(mintasksamps,[1 10])./5);
        else
            sampspersession = floor(repmat(mintasksamps,[1 10])./10);
        end
        load([QC_path sub '_QCFile.mat']);
    else              
        if MatchAcrossSubs == 1                 
            load([QC_path sub '_QCFile.mat']);                
            sampspersession = zeros(length(SubStruct),1);                
            if MatchData == 1 && SplitHalf == 1                    
                for g = 1:length(SubStruct)                    
                    sampspersession(g) = floor(mintasksamps/5);                    
                end                    
            else                
                for g = 1:length(SubStruct)                    
                    sampspersession(g) = floor(minmatch/10);                    
                end                    
            end
        end
    end

    if ConcatenateTaskData == 1
        tasks = {'allTaskCat'};  %% Set variable 'tasks' to only one task for next task loop/filenaming
    end

    if SplitHalf == 1 && ConcatenateTaskData == 1
        % Calculate the amount of task data per session across
        % all tasks
        sampspersessiontasktotal = round(sampspersession/3);
        sampspersession = zeros(length(SubStruct),1);

        for h = 1:numel(taskbackup)
            %sets up variables
            remaindertotaloddtask = 0;
            remaindertotaleventask = 0;
            notenoughdatatask = zeros(length(SubStruct),1);
            sampspersessiontask = zeros(10,1);
            motorptsodd = [];
            motorptseven = [];

            % gets motor sample points from QC file
            for u = 1:length(SubStruct)
                if SubStruct(u).OddEven == 1
                    motorptsodd = [motorptsodd; SubStruct(u).MotorSampPts];
                elseif SubStruct(u).OddEven == 2
                    motorptseven = [motorptseven; SubStruct(u).MotorSampPts];
                end
            end

            % Sets motor task data points for subject MSC09
            if strcmp(sub, 'MSC09') && strcmp(taskbackup{h}, 'motor') && ConcatenateTaskData == 1  
                if SplitHalf == 1
                    sampspersessiontask(1:2:end) = floor(sum(motorptsodd)/5);
                    sampspersessiontask(2:2:end) = floor(sum(motorptseven)/5);
                else
                    sampspersessiontask = floor(sum(motorpts)/10);
                end
            elseif strcmp(sub, 'MSC09') && ConcatenateTaskData == 1
                if SplitHalf == 1
                    sampspersessiontask(1:2:end) = round(mintasksamps/(5*numel(taskbackup))) + round(((mintasksamps/(5*numel(taskbackup))) - round(sum(motorptsodd)/5))/2);
                    sampspersessiontask(2:2:end) = round(mintasksamps/(5*numel(taskbackup))) + round(((mintasksamps/(5*numel(taskbackup))) - round(sum(motorptseven)/5))/2);
                else
                    sampspersessiontask = (mintasksamps/(10*numel(taskbackup))) + round(((mintasksamps/(10*numel(taskbackup))) - round(sum(motorpts)/10))/2);
                end
            else
                sampspersessiontask = sampspersessiontasktotal';
            end

            tasksampspersession = sampspersessiontask;

            % Goes through each task checking number of data points.
            % If not enough data to match other tasks/subs, will sample
            % from other sessions (?) until number of sample points matches
            while true
                remainderoddtask = 0;
                remaindereventask = 0;

                notenoughdata= sum(notenoughdatatask);
                notenoughdataeventask = sum(notenoughdatatask(2:2:end));
                notenoughdataoddtask = sum(notenoughdatatask(1:2:end));

                if notenoughdataeventask > 0 || notenoughdataoddtask > 0
                    for b = 1:length(SubStruct)    
                        if SubStruct(b).OddEven == 1 && notenoughdatatask(b) == 0
                            tasksampspersession(b) = (tasksampspersession(b) + round(remaindertotaloddtask/(5-notenoughdataoddtask)));
                        elseif SubStruct(b).OddEven == 2 && notenoughdatatask(b) == 0
                            tasksampspersession(b) = (tasksampspersession(b) + round(remaindertotaleventask/(5-notenoughdataeventask)));
                        end
                    end
                end

                for c = 1:length(SubStruct)
                    if SubStruct(c).OddEven == 1
                        if strcmp(taskbackup{h},'mem')
                            if SubStruct(c).MemSampPts <  tasksampspersession(c) && notenoughdatatask(c) == 0
                                notenoughdatatask(c) = 1;
                                remainderoddtask = (remainderoddtask + (tasksampspersession(c) - SubStruct(c).MemSampPts));
                                tasksampspersession(c) = SubStruct(c).MemSampPts;

                                disp(sprintf('For subject %s and session %s there are not enough sample points for the task, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            elseif notenoughdatatask(c) == 0
                                tasksampspersession(c) = tasksampspersession(c);

                                disp(sprintf('For subject %s and session %s there are enough sample points for the task, sampling %i sample points from this session: %s', sub, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            else
                                tasksampspersession(c) = SubStruct(c).MemSampPts;

                                disp(sprintf('For subject %s and session %s there are not enough sample points for the task, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            end
                        elseif strcmp(taskbackup{h},'mixed')
                            if SubStruct(c).MixedSampPts <  tasksampspersession(c) && notenoughdatatask(c) == 0
                                notenoughdatatask(c) = 1;
                                remainderoddtask = (remainderoddtask + (tasksampspersession(c) - SubStruct(c).MixedSampPts));
                                tasksampspersession(c) = SubStruct(c).MixedSampPts;

                                disp(sprintf('For subject %s and session %s there are not enough sample points for the task, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            elseif notenoughdatatask(c) == 0
                                tasksampspersession(c) = tasksampspersession(c);

                                disp(sprintf('For subject %s and session %s there are enough sample points for the task, sampling %i sample points from this session: %s', sub, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            else
                                tasksampspersession(c) = SubStruct(c).MixedSampPts

                                disp(sprintf('For subject %s and session %s there are not enough sample points for the task, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            end
                        elseif strcmp(taskbackup{h},'motor')
                            if SubStruct(c).MotorSampPts <  tasksampspersession(c) && notenoughdatatask(c) == 0
                                notenoughdatatask(c) = 1;
                                remainderoddtask = (remainderoddtask + (tasksampspersession(c) - SubStruct(c).MotorSampPts));
                                tasksampspersession(c) = SubStruct(c).MotorSampPts;

                                disp(sprintf('For subject %s and session %s there are not enough sample points for the task, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            elseif notenoughdatatask(c) == 0
                                tasksampspersession(c) = tasksampspersession(c);

                                disp(sprintf('For subject %s and session %s there are enough sample points for the task, sampling %i sample points from this session: %s', sub, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            else
                                tasksampspersession(c) = SubStruct(c).MotorSampPts;

                                disp(sprintf('For subject %s and session %s there are not enough sample points for the task, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            end
                        end
                    else
                        if strcmp(taskbackup{h},'mem')
                            if SubStruct(c).MemSampPts < tasksampspersession(c) && notenoughdatatask(c) == 0
                                notenoughdatatask(c) = 1;
                                remaindereventask = (remaindereventask + (tasksampspersession(c) - SubStruct(c).MemSampPts));
                                tasksampspersession(c) = SubStruct(c).MemSampPts;

                                disp(sprintf('For subject %s and session %s there are not enough sample points for the task, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            elseif notenoughdatatask(c) == 0
                                tasksampspersession(c) = tasksampspersession(c);

                                disp(sprintf('For subject %s and session %s there are enough sample points for the task, sampling %i sample points from this session: %s', sub, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            else
                                tasksampspersession(c) = SubStruct(c).MemSampPts;

                                disp(sprintf('For subject %s and session %s there are not enough sample points in the task data, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            end
                        elseif strcmp(taskbackup{h},'mixed')
                            if SubStruct(c).MixedSampPts < tasksampspersession(c) && notenoughdatatask(c) == 0
                                notenoughdatatask(c) = 1;
                                remaindereventask = (remaindereventask + (tasksampspersession(c) - SubStruct(c).MemSampPts));
                                tasksampspersession(c) = SubStruct(c).MixedSampPts;

                                disp(sprintf('For subject %s and session %s there are not enough sample points for the task, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            elseif notenoughdatatask(c) == 0
                                tasksampspersession(c) = tasksampspersession(c);

                                disp(sprintf('For subject %s and session %s there are enough sample points for the task, sampling %i sample points from this session: %s', sub, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            else
                                tasksampspersession(c) = SubStruct(c).MixedSampPts;

                                disp(sprintf('For subject %s and session %s there are not enough sample points in the task data, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            end
                        elseif strcmp(taskbackup{h},'motor')
                            if SubStruct(c).MixedSampPts < tasksampspersession(c) && notenoughdatatask(c) == 0
                                notenoughdatatask(c) = 1;
                                remaindereventask = (remaindereventask + (tasksampspersession(c) - SubStruct(c).MemSampPts));
                                tasksampspersession(c) = SubStruct(c).MixedSampPts;

                                disp(sprintf('For subject %s and session %s there are not enough sample points for the task, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            elseif notenoughdatatask(c) == 0
                                tasksampspersession(c) = tasksampspersession(c);

                                disp(sprintf('For subject %s and session %s there are enough sample points for the task, sampling %i sample points from this session: %s', sub, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            else
                                tasksampspersession(c) = SubStruct(c).MixedSampPts;

                                disp(sprintf('For subject %s and session %s there are not enough sample points in the task data, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            end
                        end
                    end
                end

                if remainderoddtask == 0 && remaindereventask == 0
                    disp(sprintf('Data sampling calculated for subject %s for task %s: %s', sub, taskbackup{h}, datestr(now)));
                    break
                elseif (sum(notenoughdatatask) == 10) || (sum(notenoughdataeventask) == 5 && remainderoddtask == 0) || (sum(notenoughdataoddtask) == 5 && remaindereventask == 0)
                    disp(sprintf('Not enough data is present for this task to have enough data for subject %s for task %s: %s', sub, datestr(now)));
                    break
                else
                    remaindertotaloddtask = remainderoddtask;
                    remaindertotaleventask = remaindereventask;
                    disp(sprintf('Remainder of %i calculated for odd files for subject %s: %s', remainderoddtask, sub, datestr(now)));
                    disp(sprintf('Remainder of %i calculated for even files for subject %s: %s', remaindereventask, sub, datestr(now)));
                end
            end

            sampspersession = sampspersession + tasksampspersession;

        end


        % Calculate rest data for each session based on task data
        remaindertotalodd = 0;
        remaindertotaleven = 0;
        notenoughdatarest = zeros(length(SubStruct),1);

        % checks if rest runs have enough data to match task data. If not
        % enough, samples from other rest sessions until number of data
        % points matches
        while true

            remainderodd = 0;
            remaindereven = 0;
            restsampspersession = zeros(length(SubStruct),1);

            notenoughdata = sum(notenoughdatarest);
            notenoughdataeven = sum(notenoughdatarest(2:2:end));
            notenoughdataodd = sum(notenoughdatarest(1:2:end));

            if notenoughdataeven > 0 || notenoughdataodd > 0
                for b = 1:length(SubStruct)
                    if SubStruct(b).OddEven == 1 && notenoughdatarest(b) == 0
                        sampspersession(b) = (sampspersession(b) + round(remaindertotalodd/(5-notenoughdataodd)));
                    elseif SubStruct(b).OddEven == 2 && notenoughdatarest(b) == 0
                        sampspersession(b) = (sampspersession(b) + round(remaindertotaleven/(5-notenoughdataeven)));
                    end
                end
            end

            for c = 1:length(SubStruct)
                if SubStruct(c).OddEven == 1
                    if SubStruct(c).RestSampPts <  sampspersession(c) && notenoughdatarest(c) == 0
                        notenoughdatarest(c) = 1;
                        remainderodd = (remainderodd + (sampspersession(c) - SubStruct(c).RestSampPts));
                        restsampspersession(c) = SubStruct(c).RestSampPts;

                        disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    elseif notenoughdatarest(c) == 0
                        restsampspersession(c) = sampspersession(c);

                        disp(sprintf('For subject %s and session %s there are enough sample points in the resting data to match the task, sampling %i sample points from this session: %s', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    else
                        restsampspersession(c) = SubStruct(c).RestSampPts;

                        disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    end
                else
                    if SubStruct(c).RestSampPts < sampspersession(c) && notenoughdatarest(c) == 0
                        notenoughdatarest(c) = 1;
                        remaindereven = (remaindereven + (sampspersession(c) - SubStruct(c).RestSampPts));
                        restsampspersession(c) = SubStruct(c).RestSampPts;

                        disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    elseif notenoughdatarest(c) == 0
                        restsampspersession(c) = sampspersession(c);

                        disp(sprintf('For subject %s and session %s there are enough sample points in the resting data to match the task, sampling %i sample points from this session: %s', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    else
                        restsampspersession(c) = SubStruct(c).RestSampPts;

                        disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    end
                end
            end

            if remainderodd == 0 && remaindereven == 0
                disp(sprintf('Data sampling calculated for subject %s for rest: %s', sub, datestr(now)));
                break
            elseif (sum(notenoughdata) == 10) || (sum(notenoughdataeven) == 5 && remainderodd == 0) || (sum(notenoughdataodd) == 5 && remaindereven == 0)
                disp(sprintf('Not enough data is present for this task to have enough data for subject %s for task %s: %s', sub, datestr(now)));
                break
            else
                remaindertotalodd = remainderodd;
                remaindertotaleven = remaindereven;
                disp(sprintf('Remainder of %i calculated for odd files for subject %s: %s', remainderodd, sub, datestr(now)));
                disp(sprintf('Remainder of %i calculated for even files for subject %s: %s', remaindereven, sub, datestr(now)));
            end
        end
    elseif SplitHalf == 0 && ConcatenateTaskData == 1
        remaindertotal = 0;
        notenoughdatarest = zeros(length(SubStruct),1);

        % checks if rest runs have enough data to match task data. If not
        % enough, samples from other rest sessions until number of data
        % points matches
        while true
            remainder = 0;
            restsampspersession = zeros(length(SubStruct),1);
            notenoughdata = sum(notenoughdatarest);
            if notenoughdata > 0
                for b = 1:length(SubStruct)
                    sampspersession(b) = (sampspersession(b) + round(remaindertotal/(10-notenoughdata)));
                end
            end

            for c = 1:length(SubStruct)
                if SubStruct(c).RestSampPts <  sampspersession(c) && notenoughdatarest(c) == 0
                    notenoughdatarest(c) = 1;
                    remainder = (remainder + (sampspersession(c) - SubStruct(c).RestSampPts));
                    restsampspersession(c) = SubStruct(c).RestSampPts;

                    disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                elseif notenoughdatarest(c) == 0
                    restsampspersession(c) = sampspersession(c);

                    disp(sprintf('For subject %s and session %s there are enough sample points in the resting data to match the task, sampling %i sample points from this session: %s', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                else
                    restsampspersession(c) = SubStruct(c).RestSampPts;

                    disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                end
            end

            if remainder == 0
                disp(sprintf('Data sampling calculated for subject %s for rest: %s', sub, datestr(now)));
                break
            elseif (sum(notenoughdata) == 10)
                disp(sprintf('Not enough data is present for this task to have enough data for subject %s for task %s: %s', sub, datestr(now)));
                break
            else
                remaindertotal = remainder;
                disp(sprintf('Remainder of %i calculated for full file for subject %s: %s', remainder, sub, datestr(now)));
            end
        end
    end

    if MatchData == 1 && ConcatenateTasksMatch == 1 && SplitHalf == 1 && MatchAcrossSubs == 0        
        load ([QC_path sub '_QCFile.mat']);

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

        if strcmp(sub, 'MSC09')  %% Removes motor task from consideration for MSC09
            minsampspts = min([sum(memptsodd) sum(mixedptsodd) sum(memptseven) sum(mixedptseven)]);
        else
            minsampspts = min([sum(memptsodd) sum(mixedptsodd) sum(motorptsodd) sum(memptseven) sum(mixedptseven) sum(motorptseven)]);
        end

        % Get final count for task data per session
        meansamppts = round(minsampspts/5);

        disp(sprintf('For subject %s the minimum number of sample points in a split-half is %i, with a mean of %i points per sesssion per task for concatenation: %s', sub, minsampspts, round(meansamppts/3), datestr(now)));

        sampspersessionfinal = zeros(length(SubStruct),1);

        for j = 1:numel(tasks)

            if strcmp(sub, 'MSC09') && strcmp(tasks{j}, 'motor')  %% Sets motor task data points for subject MSC09
                meansampptstempodd = floor(min([sum(motorptsodd) sum(motorptseven)])/5);
                meansampptstempeven = floor(min([sum(motorptsodd) sum(motorptseven)])/5);
            elseif strcmp(sub, 'MSC09')
                meansampptstempodd = round(meansamppts/3) + (round(meansamppts/3) - sum(motorptsodd));
                meansampptstempeven = round(meansamppts/3) + (round(meansamppts/3) - sum(motorptseven));
            else
                meansampptstempodd = round(meansamppts/3);
                meansampptstempeven = round(meansamppts/3);
            end

            remaindertotalodd = 0;
            remaindertotaleven = 0;
            notenoughdata = zeros(length(SubStruct),1);

            % checks if task runs have enough data to match other task runs. If not
            % enough, samples from other task sessions until number of data
            % points matches
            while true
            remainderodd = 0;
            remaindereven = 0;
            sampspersession = zeros(length(SubStruct),1);

            notenoughdataeven = sum(notenoughdata(2:2:end));
            notenoughdataodd = sum(notenoughdata(1:2:end));

            meansampptstempodd = (meansampptstempodd + round(remaindertotalodd/(5-notenoughdataodd)));
            meansampptstempeven = (meansampptstempeven + round(remaindertotaleven/(5-notenoughdataeven)));

            for v = 1:length(SubStruct)
                if strcmp(tasks{j}, 'mem')
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

                elseif strcmp(tasks{j}, 'mixed')
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
                elseif strcmp(tasks{j}, 'motor')
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
                disp(sprintf('Data sampling calculated for subject %s for task %s: %s', sub, tasks{j}, datestr(now)));
                break
            elseif (sum(notenoughdata) == 10) || (sum(notenoughdataeven) == 5 && remainderodd == 0) || (sum(notenoughdataodd) == 5 && remaindereven == 0)
                disp(sprintf('Not enough data is present for this task to have enough data for subject %s for task %s: %s', sub, tasks{j}, datestr(now)));
                break
            else
                remaindertotalodd = remainderodd;
                remaindertotaleven = remaindereven;
                disp(sprintf('Remainder of %i calculated for odd files for subject %s for task %s and session %s: %s', remainderodd, sub, tasks{j}, SubStruct(v).SessionFile, datestr(now)));
                disp(sprintf('Remainder of %i calculated for even files for subject %s for task %s and session %s: %s', remaindereven, sub, tasks{j}, SubStruct(v).SessionFile, datestr(now)));
            end
         end

        sampspersessionfinal = sampspersessionfinal + sampspersession;    

    end

        tasks = {'allTaskCat'};  %% Set tasks to only one task for next task loop/filenaming
        sampspersession = sampspersessionfinal;  %% Set correct name of file for the rest of the script

        %% Calculate rest data for each session based on task data
        remaindertotalodd = 0;
        remaindertotaleven = 0;

        notenoughdatarest = zeros(length(SubStruct),1);

        % checks if rest runs have enough data to match task data. If not
        % enough, samples from other rest sessions until number of data
        % points matches
        while true
            remainderodd = 0;
            remaindereven = 0;
            restsampspersession = zeros(length(SubStruct),1);

            notenoughdataeven = sum(notenoughdatarest(2:2:end));
            notenoughdataodd = sum(notenoughdatarest(1:2:end));

            if notenoughdataeven > 0 || notenoughdataodd > 0
                for b = 1:length(SubStruct)
                    if SubStruct(b).OddEven == 1 && notenoughdatarest(b) == 0
                       sampspersession(b) = (sampspersession(b) + round(remaindertotalodd/(5-notenoughdataodd)));
                    elseif SubStruct(b).OddEven == 2 && notenoughdatarest(b) == 0
                       sampspersession(b) = (sampspersession(b) + round(remaindertotaleven/(5-notenoughdataeven)));
                    end
                end
            end

            for c = 1:length(SubStruct)
                if SubStruct(c).OddEven == 1
                    if SubStruct(c).RestSampPts <  sampspersession(c) && notenoughdatarest(c) == 0
                        notenoughdatarest(c) = 1;
                        remainderodd = (remainderodd + (sampspersession(c) - SubStruct(c).RestSampPts));
                        restsampspersession(c) = SubStruct(c).RestSampPts;

                        disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    elseif notenoughdatarest(c) == 0
                        restsampspersession(c) = sampspersession(c);

                        disp(sprintf('For subject %s and session %s there are enough sample points in the resting data to match the task, sampling %i sample points from this session: %s', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    else
                        restsampspersession(c) = SubStruct(c).RestSampPts;

                        disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    end
                else
                    if SubStruct(c).RestSampPts < sampspersession(c) && notenoughdatarest(c) == 0
                        notenoughdatarest(c) = 1;
                        remaindereven = (remaindereven + (sampspersession(c) - SubStruct(c).RestSampPts));
                        restsampspersession(c) = SubStruct(c).RestSampPts;

                        disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    elseif notenoughdatarest(c) == 0
                        restsampspersession(c) = sampspersession(c);

                        disp(sprintf('For subject %s and session %s there are enough sample points in the resting data to match the task, sampling %i sample points from this session: %s', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    else
                        restsampspersession(c) = SubStruct(c).RestSampPts;

                        disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    end
                end
            end

            if remainderodd == 0 && remaindereven == 0
                disp(sprintf('Data sampling calculated for subject %s for rest: %s', sub, datestr(now)));
                break
            elseif (sum(notenoughdata) == 10) || (sum(notenoughdataeven) == 5 && remainderodd == 0) || (sum(notenoughdataodd) == 5 && remaindereven == 0)
                disp(sprintf('Not enough data is present for this task to have enough data for subject %s for task %s: %s', sub, tasks{j}, datestr(now)));
                break
            else

                remaindertotalodd = remainderodd;
                remaindertotaleven = remaindereven;

                disp(sprintf('Remainder of %i calculated for odd files for subject %s for task %s and session %s: %s', remainderodd, sub, tasks{j}, SubStruct(v).SessionFile, datestr(now)));
                disp(sprintf('Remainder of %i calculated for even files for subject %s for task %s and session %s: %s', remaindereven, sub, tasks{j}, SubStruct(v).SessionFile, datestr(now)));

            end
        end      
    end
if MatchData == 1 && ConcatenateTasksMatch == 0 && SplitHalf == 1        
        for j=1:length(tasks)
        
            load ([QC_path sub '_QCFile.mat']);

            memptsodd = [];
            motorptsodd = [];
            mixedptsodd = [];
            memptseven = [];
            motorptseven = [];
            mixedptseven = [];

            % gets sample points for tasks from QC files
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

            if strcmp(sub, 'MSC09')  %% Removes motor task from consideration for MSC09                
                minsampspts = min([sum(memptsodd) sum(mixedptsodd) sum(memptseven) sum(mixedptseven)]);                
            else        
                minsampspts = min([sum(memptsodd) sum(mixedptsodd) sum(motorptsodd) sum(memptseven) sum(mixedptseven) sum(motorptseven)]);                
            end
            % gets mean number of sample points per task per session
            if MatchAcrossSubs == 0                    
                if MatchAcrossTasks == 1
                    meansamppts = floor(minsampspts/(5*numel(tasks)));                
                    disp(sprintf('For subject %s the minimum number of sample points in a split-half is %i, with %i points per task per spit-half, with a mean of %i points per rest sesssion: %s', sub, minsampspts, meansamppts, meansamppts*numel(tasks), datestr(now)));                        
                else
                    meansamppts = floor(minsampspts/5);                
                    disp(sprintf('For subject %s the minimum number of sample points in a split-half is %i, with a mean of %i points per rest sesssion: %s', sub, minsampspts, meansamppts, datestr(now)));                        
                end            
            elseif SplitHalf == 1                    
                if MatchAcrossTasks == 1
                    meansamppts = floor(mintasksamps/(5*numel(tasks)));        
                    disp(sprintf('For all subjects, the minimum number of sample points in a split-half is %i, with %i points per task per spit-half, with a mean of %i points per rest sesssion: %s', mintasksamps, meansamppts, meansamppts*numel(tasks), datestr(now)));                        
                else                        
                    meansamppts = floor(mintasksamps/5);        
                    disp(sprintf('For all subjects, the minimum number of sample points in a split-half is %i, with a mean of %i points per rest sesssion: %s', mintasksamps, meansamppts, datestr(now)));                        
                end
                meansampptstempodd = meansamppts;
                meansampptstempeven = meansamppts;                
            else                    
                if MatchAcrossTasks == 1                
                    meansamppts = floor(mintasksamps/(10*numel(tasks)));                    
                    disp(sprintf('For all subjects, the minimum number of sample points in a full file is %i, with %i points per task, with a mean of %i points per rest sesssion: %s', mintasksamps, meansamppts, meansamppts*numel(tasks), datestr(now)));                    
                else                        
                    meansamppts = floor(mintasksamps/10);                    
                    disp(sprintf('For all subjects, the minimum number of sample points in a full file is %i, with a mean of %i points per rest sesssion: %s', mintasksamps, meansamppts, datestr(now)));                        
                end                    
                meansampptstemp = meansamppts;            
            end

            % Get final count for task data per session
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
                    if strcmp(tasks{j}, 'mem')                    
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
                    elseif strcmp(tasks{j}, 'mixed')                            
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
                    elseif strcmp(tasks{j}, 'motor')
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
                disp(sprintf('Data sampling calculated for subject %s for task %s: %s', sub, tasks{j}, datestr(now)));                
                break                    
            elseif (sum(notenoughdata) == 10) || (sum(notenoughdataeven) == 5 && remainderodd == 0) || (sum(notenoughdataodd) == 5 && remaindereven == 0)                    
                disp(sprintf('Not enough data is present for this task to have enough data for subject %s for task %s: %s', sub, tasks{j}, datestr(now)));                    
                break                
            else                
                remaindertotalodd = remainderodd;
                remaindertotaleven = remaindereven;                    
                disp(sprintf('Remainder of %i calculated for odd files for subject %s for task %s and session %s: %s', remainderodd, sub, tasks{j}, SubStruct(v).SessionFile, datestr(now)));
                disp(sprintf('Remainder of %i calculated for even files for subject %s for task %s and session %s: %s', remaindereven, sub, tasks{j}, SubStruct(v).SessionFile, datestr(now)));                
            end
        end

            % Calculate rest data for each session based on task data        
            remaindertotalodd = 0;
            remaindertotaleven = 0;                
            restsampspersession = zeros(length(SubStruct),1);
            notenoughdatarest = zeros(length(SubStruct),1);

            while true            
                remainderodd = 0;
                remaindereven = 0;
                restsampspersession = zeros(length(SubStruct),1);                    
                notenoughdataeven = sum(notenoughdatarest(2:2:end));
                notenoughdataodd = sum(notenoughdatarest(1:2:end));

                if notenoughdataeven > 0 || notenoughdataodd > 0                        
                    for b = 1:length(SubStruct)                            
                        if SubStruct(b).OddEven == 1 && notenoughdatarest(b) == 0                                
                           sampspersession(b) = (sampspersession(b) + round(remaindertotalodd/(5-notenoughdataodd)));                               
                        elseif SubStruct(b).OddEven == 2 && notenoughdatarest(b) == 0                                
                           sampspersession(b) = (sampspersession(b) + round(remaindertotaleven/(5-notenoughdataeven)));                               
                        end
                    end
                end

                for c = 1:length(SubStruct)                
                    if SubStruct(c).OddEven == 1            
                        if SubStruct(c).RestSampPts <  sampspersession(c) && notenoughdatarest(c) == 0                    
                            notenoughdatarest(c) = 1;
                            remainderodd = (remainderodd + (sampspersession(c) - SubStruct(c).RestSampPts));
                            restsampspersession(c) = SubStruct(c).RestSampPts;                    
                            disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));                    
                        elseif notenoughdatarest(c) == 0                    
                            restsampspersession(c) = sampspersession(c);                    
                            disp(sprintf('For subject %s and session %s there are enough sample points in the resting data to match the task, sampling %i sample points from this session: %s', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));                        
                        else                        
                            restsampspersession(c) = SubStruct(c).RestSampPts;                        
                            disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));                    
                        end                    
                    else                    
                        if SubStruct(c).RestSampPts < sampspersession(c) && notenoughdatarest(c) == 0                    
                            notenoughdatarest(c) = 1;
                            remaindereven = (remaindereven + (sampspersession(c) - SubStruct(c).RestSampPts));
                            restsampspersession(c) = SubStruct(c).RestSampPts;                    
                            disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));                    
                        elseif notenoughdatarest(c) == 0                    
                            restsampspersession(c) = sampspersession(c);                    
                            disp(sprintf('For subject %s and session %s there are enough sample points in the resting data to match the task, sampling %i sample points from this session: %s', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));                        
                        else                        
                            restsampspersession(c) = SubStruct(c).RestSampPts;                        
                            disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));                    
                        end                    
                    end
                end

                if remainderodd == 0 && remaindereven == 0                
                    disp(sprintf('Data sampling calculated for subject %s for and session %s: %s', sub, SubStruct(c).SessionFile, datestr(now)));                
                    break                
                else                
                    remaindertotalodd = remainderodd;
                    remaindertotaleven = remaindereven;                
                    disp(sprintf('Remainder of %i calculated for odd rest files for subject %s for session %s: %s', remainderodd, sub, SubStruct(c).SessionFile, datestr(now)));
                    disp(sprintf('Remainder of %i calculated for even rest files for subject %s for session %s: %s', remaindereven, sub, SubStruct(c).SessionFile, datestr(now)));                
                end                
            end                
        end
    end
end