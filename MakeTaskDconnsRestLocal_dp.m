parpool('local', 28)     %% Name of cluster profile for batch job

%/davta/vnil-bluearc/GMT/Caterina/TaskFC/FCProc_MSC01_motor_passvv2/
%/data/nil-bluearc/GMT/Caterina/TaskFC/FCProc_MSC01_mem_pass2/
%/data/nil-bluearc/GMT/Caterina/TaskFC/FCProc_MSC01_mixed_pass2/
clear all

RestOnly = 0; %% Toggles whether to only use resting data in dconn and not consider task data at all
CortexOnly = 1; %% Toggles cortex only template on/off
SplitHalf = 1;  %% Toggles whether to calculate rest for split-half of sessions
SubSample = 1;  %% Toggles subsampling of tmask on/off (for comparing task to rest)
RandSample = 1; %% Toggles random subsampling of tmask on/off
MatchData = 1; %% Toggles whether to match the amount of data per task as the lowest value within each split-half
MatchAcrossSubs = 1;  %% Toggles whether to match the amount of data across subjects
ConcatenateTasksMatch = 1; %% Toggles whether to calculate matched task data for concatenated tasks
ConcatenateTaskData = 1; %% Toggles whether to concatenate task data within subjects
ConcatenateSessionData = 1; %% Toggles whether to concatenate session data within subjects
ConcatenateSplitHalf = 0;  %% Toggles whether to concatenate split halves into one vector
WriteDconn = 0; %% Toggles whether to write dconn to disk
CreateVariant = 1; %% Toggles whether to create a variant map from dconn
SaveTimeseries = 0;     %% Save concatenated timeseries for subject

outdir = '/projects/b1081/Brian_MSC/Analysis_Scripts_Replication/dconn_task_files';
dataLocStem = '/MSC/TaskFC/';
subs = {'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC09','MSC10'};
%subs = {'MSC04'};
tasks = {'motor','mem','mixed'};
%tasks = {'motor'};
%times = [65 70 75 80 85 90 95 100];
times = [];

cd '/projects/b1081';   %% Change CD to root project directory

disp(sprintf('Job Submitted: %s', datestr(now)));


disp(sprintf('Job Started: %s', datestr(now)));

if numel(times) == 0
    
    times = 1;
    
end

if ConcatenateTasksMatch == 1
    
    taskbackup = tasks;
    
end


if RestOnly == 1
    
    tasks = {'rest'};
    
end

if CortexOnly == 1      %% Select correct number of voxels for template
    
    voxnum = 59412;
    
else
    
    voxnum = 65625;
    
end

if MatchData == 1 && MatchAcrossSubs == 1

   	memptsoddsum = [];
   	motorptsoddsum = [];
   	mixedptsoddsum = [];
   	memptsevensum = [];
   	motorptsevensum = [];
   	mixedptsevensum = [];

    for n=1:numel(subs)
    
        load (['/projects/b1081/Brian_MSC/Analysis_Scripts_Replication/QC_files/' subs{n} '_QCFile.mat']);
        
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
            
                if ~strcmp(subs{n}, 'MSC09')
            
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
        
        if strcmp(subs{n}, 'MSC09')
        
            motorptsoddsum = [motorptsoddsum; 9999];
            
        else
            
            motorptsoddsum = [motorptsoddsum; sum(motorptsodd)];
            
        end
        
        mixedptsoddsum = [mixedptsoddsum; sum(mixedptsodd)];
        memptsevensum = [memptsevensum; sum(memptseven)];
        
        if strcmp(subs{n}, 'MSC09')
            
            %if SplitHalf == 1
        
                motorptsevensum = [motorptsevensum; 9999];
                
            %else
                
                %motorptsevensum = [motorptsevensum; sum(motorptseven)];
                
            %end
            
        else
            
            motorptsevensum = [motorptsevensum; sum(motorptseven)];
            
        end
        
        mixedptsevensum = [mixedptsevensum; sum(mixedptseven)];
        
    end
    
    if SplitHalf == 1
    
        mintasksamps = min(min([memptsoddsum motorptsoddsum mixedptsoddsum memptsevensum motorptsevensum mixedptsevensum]));
        
    elseif ConcatenateTaskData == 1
        
        mintasksamps = min(min([(memptsoddsum + memptsevensum)  (motorptsoddsum + motorptsevensum)  (mixedptsoddsum + mixedptsevensum)]));
        
    end
     
    if ConcatenateTaskData == 1
         
         mintasksamps = mintasksamps * numel(taskbackup);
         
    end
    
end


for i=1:numel(subs)
    
    if ConcatenateTasksMatch == 1
    
        tasks = taskbackup;  %% Resets task loop for each subject
        
    end
    
            if ConcatenateTaskData == 1
                
                 sampspersession = floor(repmat(mintasksamps,[1 10])./10);
                 
                 load(['/projects/b1081/Brian_MSC/Analysis_Scripts_Replication/QC_files/' subs{i} '_QCFile.mat']);
                  
            else
                
                if MatchAcrossSubs == 1
                
                    load(['/projects/b1081/Brian_MSC/Analysis_Scripts_Replication/QC_files/' subs{i} '_QCFile.mat']);
                
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
                
                
                tasks = {'allTaskCat'};  %% Set tasks to only one task for next task loop/filenaming
                
                if SplitHalf == 1 && ConcatenateTaskData == 1
       
                    % Calculate rest data for each session based on task data
        
                    remaindertotalodd = 0;
                    remaindertotaleven = 0;

                    notenoughdatarest = zeros(length(SubStruct),1);

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
                    
                                    disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    
                                elseif notenoughdatarest(c) == 0
                    
                                    restsampspersession(c) = sampspersession(c);
                    
                                    disp(sprintf('For subject %s and session %s there are enough sample points in the resting data to match the task, sampling %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                        
                                else
                        
                                    restsampspersession(c) = SubStruct(c).RestSampPts;
                        
                                    disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    
                                end
                    
                            else
                    
                                if SubStruct(c).RestSampPts < sampspersession(c) && notenoughdatarest(c) == 0
                    
                                    notenoughdatarest(c) = 1;
                                    remaindereven = (remaindereven + (sampspersession(c) - SubStruct(c).RestSampPts));
                                    restsampspersession(c) = SubStruct(c).RestSampPts;
                    
                                    disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    
                                elseif notenoughdatarest(c) == 0
                    
                                    restsampspersession(c) = sampspersession(c);
                    
                                    disp(sprintf('For subject %s and session %s there are enough sample points in the resting data to match the task, sampling %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                        
                                else
                        
                                    restsampspersession(c) = SubStruct(c).RestSampPts;
                        
                                    disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    
                                end
                    
                            end
                        end
            
                        if remainderodd == 0 && remaindereven == 0
                    
                            disp(sprintf('Data sampling calculated for subject %s for rest: %s', subs{i}, datestr(now)));
                
                            break
                    
                        elseif (sum(notenoughdata) == 10) || (sum(notenoughdataeven) == 5 && remainderodd == 0) || (sum(notenoughdataodd) == 5 && remaindereven == 0)
                    
                            disp(sprintf('Not enough data is present for this task to have enough data for subject %s for task %s: %s', subs{i}, datestr(now)));
                    
                            break
                
                        else
                
                            remaindertotalodd = remainderodd;
                            remaindertotaleven = remaindereven;
                    
                            disp(sprintf('Remainder of %i calculated for odd files for subject %s: %s', remainderodd, subs{i}, datestr(now)));
                            disp(sprintf('Remainder of %i calculated for even files for subject %s: %s', remaindereven, subs{i}, datestr(now)));
                    
                        end
                
                    end
                    
                elseif SplitHalf == 0 && ConcatenateTaskData == 1
                    
                    remaindertotal = 0;

                    notenoughdatarest = zeros(length(SubStruct),1);

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
                    
                             	disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    
                            elseif notenoughdatarest(c) == 0
                    
                               	restsampspersession(c) = sampspersession(c);
                    
                               	disp(sprintf('For subject %s and session %s there are enough sample points in the resting data to match the task, sampling %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                        
                            else
                        
                             	restsampspersession(c) = SubStruct(c).RestSampPts;
                        
                              	disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    
                            end
                        end
            
                        if remainder == 0
                    
                            disp(sprintf('Data sampling calculated for subject %s for rest: %s', subs{i}, datestr(now)));
                
                            break
                    
                        elseif (sum(notenoughdata) == 10)
                    
                            disp(sprintf('Not enough data is present for this task to have enough data for subject %s for task %s: %s', subs{i}, datestr(now)));
                    
                            break
                
                        else
                
                            remaindertotal = remainder;
                    
                            disp(sprintf('Remainder of %i calculated for full file for subject %s: %s', remainder, subs{i}, datestr(now)));
                    
                        end
                
                    end
                    
                end
 
    
            if MatchData == 1 && ConcatenateTasksMatch == 1 && SplitHalf == 1 && MatchAcrossSubs == 0
        
                load (['/projects/b1081/Brian_MSC/Analysis_Scripts_Replication/QC_files/' subs{i} '_QCFile.mat']);
        
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
        
                if strcmp(subs{i}, 'MSC09')  %% Removes motor task from consideration for MSC09
                
                    minsampspts = min([sum(memptsodd) sum(mixedptsodd) sum(memptseven) sum(mixedptseven)]);
                
                else
        
                    minsampspts = min([sum(memptsodd) sum(mixedptsodd) sum(motorptsodd) sum(memptseven) sum(mixedptseven) sum(motorptseven)]);

                end
                
                meansamppts = round(minsampspts/5);
        
                % Get final count for task data per session
                
                %meansampptstempodd = round(meansamppts/3);
                %meansampptstempeven = round(meansamppts/3);
                

                
                disp(sprintf('For subject %s the minimum number of sample points in a split-half is %i, with a mean of %i points per sesssion per task for concatenation: %s', subs{i}, minsampspts, round(meansamppts/3), datestr(now)));

                sampspersessionfinal = zeros(length(SubStruct),1);
                
                for j = 1:numel(tasks)
                    
                	if strcmp(subs{i}, 'MSC09') && strcmp(tasks{j}, 'motor')  %% Sets motor task data points for subject MSC09
                    
                        meansampptstempodd = floor(min([sum(motorptsodd) sum(motorptseven)])/5);
                        meansampptstempeven = floor(min([sum(motorptsodd) sum(motorptseven)])/5);
                        
                    elseif strcmp(subs{i}, 'MSC09')
                        
                        meansampptstempodd = round(meansamppts/3) + (round(meansamppts/3) - sum(motorptsodd));
                        meansampptstempeven = round(meansamppts/3) + (round(meansamppts/3) - sum(motorptseven));
                        
                    else
                        
                        meansampptstempodd = round(meansamppts/3);
                        meansampptstempeven = round(meansamppts/3);
                    
                    end
                
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
                    
                                    disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                                elseif notenoughdata(v) == 0
                    
                                    sampspersession(v) = meansampptstempodd;
                    
                                    disp(sprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                                
                                else
                                
                                    sampspersession(v) = SubStruct(v).MemSampPts;
                                
                                    disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                                end
                    
                            else
                    
                                if SubStruct(v).MemSampPts < meansampptstempeven && notenoughdata(v) == 0
                    
                                    notenoughdata(v) = 1;
                                    remaindereven = (remaindereven + (meansampptstempeven - SubStruct(v).MemSampPts));
                                    sampspersession(v) = SubStruct(v).MemSampPts;
                    
                                    disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                                elseif notenoughdata(v) == 0
                    
                                    sampspersession(v) = meansampptstempeven;
                    
                                    disp(sprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                                
                                else
                                
                                    sampspersession(v) = SubStruct(v).MemSampPts;
                                
                                    disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                                end
                    
                            end
                            
                        elseif strcmp(tasks{j}, 'mixed')
                            
                            if SubStruct(v).OddEven == 1

                                if SubStruct(v).MixedSampPts < meansampptstempodd && notenoughdata(v) == 0
                    
                                    notenoughdata(v) = 1;
                                    remainderodd = (remainderodd + (meansampptstempodd - SubStruct(v).MixedSampPts));
                                    sampspersession(v) = SubStruct(v).MixedSampPts;
                    
                                    disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                                elseif notenoughdata(v) == 0
                    
                                    sampspersession(v) = meansampptstempodd;
                    
                                    disp(sprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                                
                                else
                                
                                    sampspersession(v) = SubStruct(v).MixedSampPts;
                                
                                    disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                                end
                    
                            else
                    
                                if SubStruct(v).MixedSampPts < meansampptstempeven && notenoughdata(v) == 0
                    
                                    notenoughdata(v) = 1;
                                    remaindereven = (remaindereven + (meansampptstempeven - SubStruct(v).MixedSampPts));
                                    sampspersession(v) = SubStruct(v).MixedSampPts;
                    
                                    disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                                elseif notenoughdata(v) == 0
                    
                                    sampspersession(v) = meansampptstempeven;
                    
                                    disp(sprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                                
                                else
                                
                                    sampspersession(v) = SubStruct(v).MixedSampPts;
                                
                                    disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                                end
                            end
                            
                        elseif strcmp(tasks{j}, 'motor')
                            
                            if SubStruct(v).OddEven == 1

                                if SubStruct(v).MotorSampPts < meansampptstempodd && notenoughdata(v) == 0
                    
                                    notenoughdata(v) = 1;
                                    remainderodd = (remainderodd + (meansampptstempodd - SubStruct(v).MotorSampPts));
                                    sampspersession(v) = SubStruct(v).MotorSampPts;
                    
                                    disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                                elseif notenoughdata(v) == 0
                    
                                    sampspersession(v) = meansampptstempodd;
                    
                                    disp(sprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                                
                                else
                                
                                    sampspersession(v) = SubStruct(v).MotorSampPts;
                                
                                    disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                                end
                    
                            else
                    
                                if SubStruct(v).MotorSampPts < meansampptstempeven && notenoughdata(v) == 0
                    
                                    notenoughdata(v) = 1;
                                    remaindereven = (remaindereven + (meansampptstempeven - SubStruct(v).MotorSampPts));
                                    sampspersession(v) = SubStruct(v).MotorSampPts;
                    
                                    disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                                elseif notenoughdata(v) == 0
                    
                                    sampspersession(v) = meansampptstempeven;
                    
                                    disp(sprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                                
                                else
                                
                                    sampspersession(v) = SubStruct(v).MotorSampPts;
                                
                                    disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                                end
                            
                            end
                        end
                    end
            
                    if remainderodd == 0 && remaindereven == 0
                    
                        disp(sprintf('Data sampling calculated for subject %s for task %s: %s', subs{i}, tasks{j}, datestr(now)));
                
                        break
                    
                    elseif (sum(notenoughdata) == 10) || (sum(notenoughdataeven) == 5 && remainderodd == 0) || (sum(notenoughdataodd) == 5 && remaindereven == 0)
                    
                        disp(sprintf('Not enough data is present for this task to have enough data for subject %s for task %s: %s', subs{i}, tasks{j}, datestr(now)));
                    
                        break
                
                    else
                
                        remaindertotalodd = remainderodd;
                        remaindertotaleven = remaindereven;
                    
                        disp(sprintf('Remainder of %i calculated for odd files for subject %s for task %s and session %s: %s', remainderodd, subs{i}, tasks{j}, SubStruct(v).SessionFile, datestr(now)));
                        disp(sprintf('Remainder of %i calculated for even files for subject %s for task %s and session %s: %s', remaindereven, subs{i}, tasks{j}, SubStruct(v).SessionFile, datestr(now)));
                
                    end
                 end
                 
                sampspersessionfinal = sampspersessionfinal + sampspersession;    
                    
                
                end
                
                tasks = {'allTaskCat'};  %% Set tasks to only one task for next task loop/filenaming
                sampspersession = sampspersessionfinal;  %% Set correct name of file for the rest of the script
       
                % Calculate rest data for each session based on task data
        
                remaindertotalodd = 0;
                remaindertotaleven = 0;

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
                    
                                disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    
                            elseif notenoughdatarest(c) == 0
                    
                                restsampspersession(c) = sampspersession(c);
                    
                                disp(sprintf('For subject %s and session %s there are enough sample points in the resting data to match the task, sampling %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                        
                            else
                        
                                restsampspersession(c) = SubStruct(c).RestSampPts;
                        
                                disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    
                            end
                    
                        else
                    
                            if SubStruct(c).RestSampPts < sampspersession(c) && notenoughdatarest(c) == 0
                    
                                notenoughdatarest(c) = 1;
                                remaindereven = (remaindereven + (sampspersession(c) - SubStruct(c).RestSampPts));
                                restsampspersession(c) = SubStruct(c).RestSampPts;
                    
                                disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    
                            elseif notenoughdatarest(c) == 0
                    
                                restsampspersession(c) = sampspersession(c);
                    
                                disp(sprintf('For subject %s and session %s there are enough sample points in the resting data to match the task, sampling %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                        
                            else
                        
                                restsampspersession(c) = SubStruct(c).RestSampPts;
                        
                                disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    
                            end
                    
                        end
                    end
            
                    if remainderodd == 0 && remaindereven == 0
                    
                        disp(sprintf('Data sampling calculated for subject %s for rest: %s', subs{i}, datestr(now)));
                
                        break
                    
                    elseif (sum(notenoughdata) == 10) || (sum(notenoughdataeven) == 5 && remainderodd == 0) || (sum(notenoughdataodd) == 5 && remaindereven == 0)
                    
                        disp(sprintf('Not enough data is present for this task to have enough data for subject %s for task %s: %s', subs{i}, tasks{j}, datestr(now)));
                    
                        break
                
                    else
                
                        remaindertotalodd = remainderodd;
                        remaindertotaleven = remaindereven;
                    
                        disp(sprintf('Remainder of %i calculated for odd files for subject %s for task %s and session %s: %s', remainderodd, subs{i}, tasks{j}, SubStruct(v).SessionFile, datestr(now)));
                        disp(sprintf('Remainder of %i calculated for even files for subject %s for task %s and session %s: %s', remaindereven, subs{i}, tasks{j}, SubStruct(v).SessionFile, datestr(now)));
                
                    end
                
                end
      
            end
    

    for l=1:numel(times)
        
        if numel(times) == 1
    
            disp(sprintf('Creating dconn for subject %s: %s', subs{i}, datestr(now)));
        
        else
            
            disp(sprintf('Creating dconn for subject %s for time %i minutes: %s', subs{i}, times(l), datestr(now)));
            
        end
    
    
        % Initialize cat data
        
        if SplitHalf == 1
            
            catData1 = [];
            catData2 = [];
            
        else
            
            catData = [];
            
        end
    
    
        for j=1:length(tasks)
        
            disp(sprintf('Creating resting dconn for subject %s task %s: %s', subs{i}, tasks{j}, datestr(now)));
            
            if MatchData == 1 && ConcatenateTasksMatch == 0 && SplitHalf == 1
        
                load (['/projects/b1081/Brian_MSC/QC_files/' subs{i} '_QCFile.mat']);
        
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
        
                if strcmp(subs{i}, 'MSC09')  %% Removes motor task from consideration for MSC09
                
                    minsampspts = min([sum(memptsodd) sum(mixedptsodd) sum(memptseven) sum(mixedptseven)]);
                
                else
        
                    minsampspts = min([sum(memptsodd) sum(mixedptsodd) sum(motorptsodd) sum(memptseven) sum(mixedptseven) sum(motorptseven)]);
                
                end
                
                
                meansamppts = floor(minsampspts/5);
        
                disp(sprintf('For subject %s the minimum number of sample points in a split-half is %i, with a mean of %i points per sesssion: %s', subs{i}, minsampspts, meansamppts, datestr(now)));
        
                % Get final count for task data per session
        
                meansampptstempodd = meansamppts;
                meansampptstempeven = meansamppts;
                
                
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
                    
                            	disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                            elseif notenoughdata(v) == 0
                    
                            	sampspersession(v) = meansampptstempodd;
                    
                                disp(sprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                                
                            else
                                
                                sampspersession(v) = SubStruct(v).MemSampPts;
                                
                                disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                            end
                    
                        else
                    
                        	if SubStruct(v).MemSampPts < meansampptstempeven && notenoughdata(v) == 0
                    
                            	notenoughdata(v) = 1;
                                remaindereven = (remaindereven + (meansampptstempeven - SubStruct(v).MemSampPts));
                                sampspersession(v) = SubStruct(v).MemSampPts;
                    
                                disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                            elseif notenoughdata(v) == 0
                    
                                sampspersession(v) = meansampptstempeven;
                    
                                disp(sprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                                
                            else
                                
                                sampspersession(v) = SubStruct(v).MemSampPts;
                                
                                disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                            end
                    
                        end
                            
                    elseif strcmp(tasks{j}, 'mixed')
                            
                    	 if SubStruct(v).OddEven == 1

                        	if SubStruct(v).MixedSampPts < meansampptstempodd && notenoughdata(v) == 0
                    
                                notenoughdata(v) = 1;
                            	remainderodd = (remainderodd + (meansampptstempodd - SubStruct(v).MixedSampPts));
                            	sampspersession(v) = SubStruct(v).MixedSampPts;
                    
                            	disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                            elseif notenoughdata(v) == 0
                    
                            	sampspersession(v) = meansampptstempodd;
                    
                                disp(sprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                                
                            else
                                
                                sampspersession(v) = SubStruct(v).MixedSampPts;
                                
                                disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                            end
                    
                        else
                    
                        	if SubStruct(v).MixedSampPts < meansampptstempeven && notenoughdata(v) == 0
                    
                            	notenoughdata(v) = 1;
                                remaindereven = (remaindereven + (meansampptstempeven - SubStruct(v).MixedSampPts));
                                sampspersession(v) = SubStruct(v).MixedSampPts;
                    
                                disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                            elseif notenoughdata(v) == 0
                    
                                sampspersession(v) = meansampptstempeven;
                    
                                disp(sprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                                
                            else
                                
                                sampspersession(v) = SubStruct(v).MixedSampPts;
                                
                                disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                            end
                        end
                            
                  	elseif strcmp(tasks{j}, 'motor')
                            
                     	if SubStruct(v).OddEven == 1

                        	if SubStruct(v).MotorSampPts < meansampptstempodd && notenoughdata(v) == 0
                    
                                notenoughdata(v) = 1;
                            	remainderodd = (remainderodd + (meansampptstempodd - SubStruct(v).MotorSampPts));
                            	sampspersession(v) = SubStruct(v).MotorSampPts;
                    
                            	disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                            elseif notenoughdata(v) == 0
                    
                            	sampspersession(v) = meansampptstempodd;
                    
                                disp(sprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                                
                            else
                                
                                sampspersession(v) = SubStruct(v).MotorSampPts;
                                
                                disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                            end
                    
                        else
                    
                        	if SubStruct(v).MotorSampPts < meansampptstempeven && notenoughdata(v) == 0
                    
                            	notenoughdata(v) = 1;
                                remaindereven = (remaindereven + (meansampptstempeven - SubStruct(v).MotorSampPts));
                                sampspersession(v) = SubStruct(v).MotorSampPts;
                    
                                disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                            elseif notenoughdata(v) == 0
                    
                                sampspersession(v) = meansampptstempeven;
                    
                                disp(sprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                                
                            else
                                
                                sampspersession(v) = SubStruct(v).MotorSampPts;
                                
                                disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                            end
                            
                        end
                    end
                end
            
                if remainderodd == 0 && remaindereven == 0
                    
                    disp(sprintf('Data sampling calculated for subject %s for task %s: %s', subs{i}, tasks{j}, datestr(now)));
                
                    break
                    
                elseif (sum(notenoughdata) == 10) || (sum(notenoughdataeven) == 5 && remainderodd == 0) || (sum(notenoughdataodd) == 5 && remaindereven == 0)
                    
                    disp(sprintf('Not enough data is present for this task to have enough data for subject %s for task %s: %s', subs{i}, tasks{j}, datestr(now)));
                    
                    break
                
                else
                
                	remaindertotalodd = remainderodd;
                    remaindertotaleven = remaindereven;
                    
                    disp(sprintf('Remainder of %i calculated for odd files for subject %s for task %s and session %s: %s', remainderodd, subs{i}, tasks{j}, SubStruct(v).SessionFile, datestr(now)));
                    disp(sprintf('Remainder of %i calculated for even files for subject %s for task %s and session %s: %s', remaindereven, subs{i}, tasks{j}, SubStruct(v).SessionFile, datestr(now)));
                
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
                    
                                disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    
                            elseif notenoughdatarest(c) == 0
                    
                                restsampspersession(c) = sampspersession(c);
                    
                                disp(sprintf('For subject %s and session %s there are enough sample points in the resting data to match the task, sampling %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                        
                            else
                        
                                restsampspersession(c) = SubStruct(c).RestSampPts;
                        
                                disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    
                            end
                    
                        else
                    
                            if SubStruct(c).RestSampPts < sampspersession(c) && notenoughdatarest(c) == 0
                    
                                notenoughdatarest(c) = 1;
                                remaindereven = (remaindereven + (sampspersession(c) - SubStruct(c).RestSampPts));
                                restsampspersession(c) = SubStruct(c).RestSampPts;
                    
                                disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    
                            elseif notenoughdatarest(c) == 0
                    
                                restsampspersession(c) = sampspersession(c);
                    
                                disp(sprintf('For subject %s and session %s there are enough sample points in the resting data to match the task, sampling %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                        
                            else
                        
                                restsampspersession(c) = SubStruct(c).RestSampPts;
                        
                                disp(sprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, restsampspersession(c), datestr(now)));
                    
                            end
                    
                        end
                    end
            
                    if remainderodd == 0 && remaindereven == 0
                
                        disp(sprintf('Data sampling calculated for subject %s for and session %s: %s', subs{i}, SubStruct(c).SessionFile, datestr(now)));
                
                        break
                
                    else
                
                        remaindertotalodd = remainderodd;
                        remaindertotaleven = remaindereven;
                
                        disp(sprintf('Remainder of %i calculated for odd rest files for subject %s for session %s: %s', remainderodd, subs{i}, SubStruct(c).SessionFile, datestr(now)));
                        disp(sprintf('Remainder of %i calculated for even rest files for subject %s for session %s: %s', remaindereven, subs{i}, SubStruct(c).SessionFile, datestr(now)));
                
                    end
                
                end
                
            end
        
            if RestOnly == 0 && MatchData == 0 && ConcatenateTaskData == 0
        
                % Load task vcids
                MSCTaskdir = strcat(cd, dataLocStem);
                %[~,vcids,~,~,~] = textread([MSCTaskdir '/' subs{i} '_' tasks{j} '_DATALIST.txt'],'%s%s%s%s%s');
        
        
                % Load Task tmasks
                
                if strcmp(subs{i}, 'MSC03') || strcmp(subs{i}, 'MSC10')
                    
                    MSCcondidir = strcat(strcat(MSCTaskdir, '/'), ['FCProc_' subs{i} '_' tasks{j} '_pass2_FDfilt/']);
                    
                else
                
                    MSCcondidir = strcat(strcat(MSCTaskdir, '/'), ['FCProc_' subs{i} '_' tasks{j} '_pass2/']);
                    
                end
                
                load ([MSCcondidir 'condindices.mat']);
            
            else
            
                % Create template path for resting data
                MSCciftidir = ['/projects/b1081/MSC/TaskFC/FCProc_' subs{i} '_mem_pass2/cifti_timeseries_normalwall_native_freesurf'];
        
            end
        
            % Load rest vcids
            restdir = dir(['/projects/b1081/MSC/MSCdata_v1/' subs{i} '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_normalwall_native_freesurf']);
            vcidlist = restdir(~cellfun(@isempty,strfind({restdir.name},'dtseries')));
            vcidlist = vcidlist(~cellfun(@isempty,strfind({vcidlist.name},'vc')));
            
            disp(sprintf('%i resting sessions found for subject %s: %s', size(vcidlist,1), subs{i}, datestr(now)));
        
            % Load rest tmasks
            load(['/projects/b1081/MSC/MSCdata_v1/' subs{i} '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/QC.mat'])
        


            % Load and concatentate data
            for k=1:length(vcidlist)
            
            vcid = vcidlist(k).name;
            
                disp(sprintf('Loading timeseries for session %s for subject %s task %s: %s', vcid, subs{i}, tasks{j}, datestr(now)));
            
                if RestOnly == 0 && MatchData == 0 && ConcatenateTaskData == 0
            
                    MSCciftidir = strcat(MSCcondidir, 'cifti_timeseries_normalwall_native_freesurf');
            
                end

            
                %if k==1
                    %cd('cifti_timeseries_normalwall_native_freesurf')
                %end
            
                tmask = 0;  % Create empty tmask variable
            
                try
                    data = ft_read_cifti_mod(['/projects/b1081/MSC/MSCdata_v1/' subs{i} '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_normalwall_native_freesurf/' vcid]);
                    disp(sprintf('Loading data size %i by %i, %s', size(data.data,1), size(data.data,2), datestr(now)));
                    data = data.data;
                
                    if RestOnly == 0 && MatchData == 0 && ConcatenateTaskData == 0
                
                        if strcmp(tasks{j},'motor')
                            tmask = TIndFin(k).AllMotor;
                        elseif strcmp(tasks{j},'mem')
                            tmask = TIndFin(k).AllMem;
                        else
                            tmask = TIndFin(k).AllGlass + TIndFin(k).AllSemantic;
                        end
                        
                        disp(sprintf('tmask for task %s file has %i good sample points, %s', tasks{j}, sum(tmask), datestr(now)));
                
                    end
                
                    resttmask = QC(k).tmask;
                    
                    disp(sprintf('tmask for rest file has %i good sample points, %s', sum(resttmask), datestr(now)));
                
                    if RestOnly == 0 && sum(tmask) == 0 && MatchData == 0 && ConcatenateTaskData == 0    % Skip file if no task task data
                    
                        data = [];
                        disp(sprintf('Vcid %s for subject %s has no usable data for task %s', vcid, subs{i}, tasks{j}));
                    
                    else                    % Else get first xth points in resting file corresponding to number of task data points
                 
                    
                        if RestOnly == 1 && SubSample == 1 && RandSample == 1 && numel(times) > 1
                            
                            subsamp = (10/(2*SplitHalf));  %% Samples half of each value in time from each session
                            
                            if subsamp == Inf  %% Sets divisior to 1 if SplitHalf set to 0
                                    
                            	subsamp = 1;
                                    
                            end
                        
                            tmasknum = round((60*times(l)/2.2)/subsamp);     % Get number of samples for each time
                            
                            samplepts = find(resttmask == 1);
                            
                            if tmasknum > sum(resttmask)
                                
                                disp(sprintf('Vcid %s for subject %s does not have %i minutes of resting data', vcid, subs{i}, times(l)));
                                
                                data = data(1:voxnum,logical(resttmask));
                                
                            else
                        
                                resttmask(samplepts) = 0;
                        
                                sampleselect = datasample(samplepts,tmasknum,'Replace',false);
                        
                                resttmask(sampleselect) = 1;
                                
                                disp(sprintf('%i sample points from rest file %s have been randomly selected for task %s, %s', sum(resttmask), vcid, tasks{j}, datestr(now)));
                            
                                data = data(1:voxnum,logical(resttmask));
                            
                            end
                            
                        elseif RestOnly == 1 && SubSample == 1 && numel(times) > 1
                            
                            subsamp = (10/(2*SplitHalf));  %% Samples half of each value in time from each session
                            
                            if subsamp == Inf  %% Sets divisior to 1 if SplitHalf set to 0
                                    
                            	subsamp = 1;
                                    
                            end
                        
                            tmasknum = round((60*times(l)/2.2)/subsamp);     % Get number of samples for each time
                            
                            if tmasknum > sum(resttmask)
                                
                                disp(sprintf('Vcid %s for subject %s does not have %i minutes of resting data', vcid, subs{i}, times(l)));
                                
                                data = data(1:voxnum,logical(resttmask));
                                
                            else
                
                            data = data(1:voxnum,logical(resttmask));        %% Get good data from rest
                            data = data(:,1:tmasknum);          %% Get correct amount of points corresponding to task
                            
                            disp(sprintf('%i sample points from rest file %s have been sequentially selected for task %s, %s', size(data,2), vcid, tasks{j}, datestr(now)));
                            
                            end
                            
                        elseif RestOnly == 1 && numel(times) == 1
                            
                            data = data(1:voxnum,logical(resttmask));
                            
                        elseif SubSample == 1 && RandSample == 1 && MatchData == 1    %% Subsample (rest) data for the desired number of voxels as set by the lowest split-half for a random subset of voxels
                            
                            tmasknum = restsampspersession(k);
                        
                            samplepts = find(resttmask == 1);
                        
                            resttmask(samplepts) = 0;
                        
                            sampleselect = datasample(samplepts,tmasknum,'Replace',false);
                        
                            resttmask(sampleselect) = 1;
                            
                            disp(sprintf('%i sample points from rest file %s have been randomly selected for task %s, %s', sum(resttmask), vcid, tasks{j}, datestr(now)));
                            
                            data = data(1:voxnum,logical(resttmask));
                            
                        elseif SubSample == 1 && RandSample == 1 && ConcatenateTaskData == 1
                            
                            tmasknum = restsampspersession(k);
                        
                            samplepts = find(resttmask == 1);
                        
                            resttmask(samplepts) = 0;
                        
                            sampleselect = datasample(samplepts,tmasknum,'Replace',false);
                        
                            resttmask(sampleselect) = 1;
                            
                            disp(sprintf('%i sample points from rest file %s have been randomly selected for task %s, %s', sum(resttmask), vcid, tasks{j}, datestr(now)));
                            
                            data = data(1:voxnum,logical(resttmask));
                    
                        elseif  SubSample == 1 && RandSample == 1       %% Subsample (rest) data for the desired number of voxels for a random subset of voxels
                        
                            tmasknum = sum(tmask);
                        
                            samplepts = find(resttmask == 1);
                        
                            resttmask(samplepts) = 0;
                        
                            sampleselect = datasample(samplepts,tmasknum,'Replace',false);
                        
                            resttmask(sampleselect) = 1;
                            
                            disp(sprintf('%i sample points from rest file %s have been randomly selected for task %s, %s', sum(resttmask), vcid, tasks{j}, datestr(now)));
                            
                            data = data(1:voxnum,logical(resttmask));    
                        
                        elseif SubSample == 1     %% Subsample data from the first voxel to the desired number of voxels (in this case task number per session)
                    
                            tmasknum = sum(tmask);
                
                            data = data(1:voxnum,logical(resttmask));        %% Get good data from rest
                            data = data(:,1:tmasknum);          %% Get correct amount of points corresponding to task
                    
                        else    %% Select all usable data points from tmask
                    
                            data = data(1:voxnum,logical(tmask));
                    
                        end
                        
                        if SplitHalf == 1
                            
                            if mod(k,2) == 1
                                
                                disp(sprintf('Adding %i data points from session %s for subject %s and time %i minutes to odd data', size(data,2), vcid, subs{i}, times(l)));
                        
                                catData1 = [catData1 data];
                            
                            else
                                
                                disp(sprintf('Adding %i data points from session %s for subject %s and time %i minutes to even data', size(data,2), vcid, subs{i}, times(l)));
                                
                                catData2 = [catData2 data];
                            
                            end
                         
                        else
                        
                            catData = [catData data];
                        
                        end
                
                    end
                
                catch ME
                    if strcmp(ME.message,'Invalid file identifier.  Use fopen to generate a valid file identifier.')
                        continue
                    end
                end
            end         % End of session loop to create concatenated session data   
        end        % End of task loop to create concatenated task data
        %across all tasks within each session

        
        if ConcatenateSplitHalf == 1
            
            catData = [catData1 catData2];
            
        end
        
        if SaveTimeseries == 1          %% Save concatenated FC processed timeseries
        
            timseriestemplate = ft_read_cifti_mod('/projects/b1081/Brian_MSC/Analysis_Scripts_Replication/dconn_scripts/templates/MSC01_allses_mean_native_freesurf_vs_120sub_corr.dtseries.nii');
            timseriestemplate.data = [];
        
            timseriestemplate.data = catData;
            
            disp(sprintf('Data timseries is size %i by %i, %s', size(catData,1), size(catData,2), datestr(now)));
        
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_FCProcessed_Timeseries_Matched_REST_cortex'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_FCProcessed_Timeseries_Matched_REST_cortex'],timseriestemplate)
        
        end
    
        % Make and save rmat
    
        if CortexOnly == 1
        
            disp('Loading template: MSC01_allses_mean_native_freesurf_vs_120sub_corr.dtseries.nii');
            template = ft_read_cifti_mod('/projects/b1081/Brian_MSC/Analysis_Scripts_Replication/dconn_scripts/templates/MSC01_allses_mean_native_freesurf_vs_120sub_corr.dtseries.nii');    
        
        else
    
            disp(sprintf('Reading %s: %s', [MSCciftidir '/' vcid], datestr(now)));
            template = ft_read_cifti_mod([MSCciftidir '/' vcid]);
    
    
        end
        
    
        template.data = [];
        template.dimord = 'pos_pos';
        template.hdr.dim(7) = voxnum;
        template.hdr.dim(6) = template.hdr.dim(7);
        template.hdr.intent_code = 3001;
        template.hdr.intent_name = 'ConnDense';
        
        if SplitHalf == 1 && ConcatenateSplitHalf == 0
            
            template2 = template;
            
        end
        
        %template.hdr.intent_name = 'ConnDense       ';
        %template.hdr.descrip = pad('',80);
        
        if SplitHalf == 1 && WriteDconn == 1 && ConcatenateSplitHalf == 0
            
            disp(sprintf('Running Correlations: on data size %i by %i, %s', size(catData1,1), size(catData1,2), datestr(now)));
            template.data = paircorr_mod(catData1');
            clear catData1
            disp(sprintf('Running Correlations: on data size %i by %i, %s', size(catData2,1), size(catData2,2), datestr(now)));
            template2.data = paircorr_mod(catData2');
            clear catData2
            
        elseif SplitHalf == 0 || ConcatenateSplitHalf == 1
        
            disp(sprintf('Running Correlations: on data size %i by %i, %s', size(catData,1), size(catData,2), datestr(now)));
            template.data = paircorr_mod(catData');
            
        end
        
        if SplitHalf == 1 && WriteDconn == 1 && ConcatenateSplitHalf == 0
            
            catData1 = [];
            catData2 = [];
            
        elseif SplitHalf == 0 || ConcatenateSplitHalf == 1
        
        	catData = [];   % Clear cat data
            
        end
        %disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_allTaskCat'], datestr(now)));
        %ft_write_cifti_mod([outdir '/' subs{i} '_allTaskCat'],template)
        
        if ConcatenateSessionData == 0 && CortexOnly == 1 && RestOnly == 1
            
            if WriteDconn == 1
                
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' vcid '_' 'REST_AllSessions_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' vcid '_' 'REST_AllSessions_cortex'],template)
                
            end
            
            if CreateVariant == 1
                
                outputfile = [subs{i} '_' vcid '_' 'REST_AllSessions'];
                
            end
            
        elseif ConcatenateSessionData == 0 && RestOnly == 1
            
            if WriteDconn == 1
                
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' vcid '_' 'REST_AllSessions'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' vcid '_' 'REST_AllSessions'],template)
                
            end
            
            if CreateVariant == 1
                
                outputfile = [subs{i} '_' vcid '_' 'REST_AllSessions'];
                
            end
            
        elseif CortexOnly == 1 && SplitHalf ==1 && MatchAcrossSubs == 1 && MatchData == 1
            
            if WriteDconn == 1
            
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' tasks{j} '_matcheddata_Final_REST_OddSessions_Rand_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_matcheddata_Final_REST_OddSessions_Rand_cortex'],template)
            
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' tasks{j} '_matcheddata_Final_REST_EvenSessions_Rand_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_matcheddata_Final_REST_EvenSessions_Rand_cortex'],template2)
            
            end
            
            if CreateVariant == 1
                
                outputfile1 = [subs{i} '_' tasks{j} '_matcheddata_Final_REST_OddSessions_Rand'];
                outputfile2 = [subs{i} '_' tasks{j} '_matcheddata_Final_REST_EvenSessions_Rand'];
                
            end
        
        elseif RestOnly == 1 && CortexOnly == 1 && SplitHalf ==1 && numel(times) > 1
            
            if WriteDconn == 1
            
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' num2str(times(l)) '_minutes_total' '_' 'REST_OddSessions_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' num2str(times(l)) '_minutes_total' '_' 'REST_OddSessions_cortex'],template)
            
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' num2str(times(l)) '_minutes_total' '_' 'REST_EvenSessions_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' num2str(times(l)) '_minutes_total' '_' 'REST_EvenSessions_cortex'],template2)
            
            end
            
            if CreateVariant == 1
                
                outputfile1 = [subs{i} '_' num2str(times(l)) '_minutes_total' '_' 'REST_OddSessions'];
                outputfile2 = [subs{i} '_' num2str(times(l)) '_minutes_total' '_' 'REST_EvenSessions'];
                
            end
            
        elseif RestOnly == 1 && CortexOnly == 1 && SplitHalf ==1  && numel(times) == 1
            
            if WriteDconn == 1
            
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' 'All' '_' 'REST_OddSessions_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' 'All' '_' 'REST_OddSessions_cortex'],template)
            
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' 'All' '_' 'REST_EvenSessions_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' 'All' '_' 'REST_EvenSessions_cortex'],template2)
            
            end
            
            if CreateVariant == 1
                
                outputfile1 = [subs{i} '_' 'All' '_' 'REST_OddSessions'];
                outputfile2 = [subs{i} '_' 'All' '_' 'REST_EvenSessions'];
                
            end
    
        elseif RestOnly == 1 && CortexOnly == 1 && numel(times) > 1         %% Ouput files with correct names
            
        	if WriteDconn == 1  
            
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' num2str(times(l)) '_minutes_total' '_' 'REST_AllSessions_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' num2str(times(l)) '_minutes_total' '_' 'REST_AllSessions_cortex'],template)
            
            end
            
            if CreateVariant == 1
                
                outputfile = [subs{i} '_' num2str(times(l)) '_minutes_total' '_' 'REST_AllSessions'];
                
            end
                
        elseif RestOnly == 1 && numel(times) > 1
            
            if WriteDconn == 1
            
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' num2str(times(l)) '_minutes' '_' 'REST_AllSessions'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' num2str(times(l)) '_minutes' '_' 'REST_AllSessions'],template)
                
            end
            
            if CreateVariant == 1
                
                outputfile = [subs{i} '_' num2str(times(l)) '_minutes' '_' 'REST_AllSessions'];
                
            end
        
        elseif RestOnly == 1 && CortexOnly == 1
            
            if WriteDconn == 1
        
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' 'REST_AllSessions_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' 'REST_AllSessions_cortex'],template)
                
            end
            
            if CreateVariant == 1
                
                outputfile = [subs{i} '_' 'REST_AllSessions'];
                
            end
        
        elseif RestOnly == 1
            
            if WriteDconn == 1
        
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' 'REST_AllSessions'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' 'REST_AllSessions'],template)
                
            end
            
             if CreateVariant == 1
                 
                 outputfile = [subs{i} '_' 'REST_AllSessions'];
                 
             end
             
        elseif SplitHalf == 1 && SubSample == 1 && RandSample == 1 && CortexOnly == 1 && MatchData == 1
            
            if WriteDconn == 1
            
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' tasks{j} '_matcheddata_' 'REST_OddSessions_Rand_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_matcheddata_' 'REST_OddSessions_Rand_cortex'],template)
            
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' tasks{j} '_matcheddata_' 'REST_EvenSessions_Rand_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_matcheddata_' 'REST_EvenSessions_Rand_cortex'],template2)
            
            end
            
            if CreateVariant == 1
                
                outputfile1 = [subs{i} '_' tasks{j} '_matcheddata_' 'REST_OddSessions_Rand'];
                outputfile2 = [subs{i} '_' tasks{j} '_matcheddata_' 'REST_EvenSessions_Rand'];
                
            end
            
        elseif SplitHalf == 1 && SubSample == 1 && RandSample == 1 && CortexOnly == 1 && ConcatenateTaskData == 1
            
            if WriteDconn == 1
            
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' tasks{j} '_' 'REST_OddSessions_Rand_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_' 'REST_OddSessions_Rand_cortex'],template)
            
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' tasks{j} '_' 'REST_EvenSessions_Rand_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_' 'REST_EvenSessions_Rand_cortex'],template2)
            
            end
            
            if CreateVariant == 1
                
                outputfile1 = [subs{i} '_' tasks{j} '_' 'REST_OddSessions_Rand'];
                outputfile2 = [subs{i} '_' tasks{j} '_' 'REST_EvenSessions_Rand'];
                
            end
            
             
        elseif SplitHalf == 1 && SubSample == 1 && RandSample == 1 && CortexOnly == 1
            
            if WriteDconn == 1
            
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' tasks{j} '_' 'REST_OddSessions_Rand_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_' 'REST_OddSessions_Rand_cortex'],template)
            
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' tasks{j} '_' 'REST_EvenSessions_Rand_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_' 'REST_EvenSessions_Rand_cortex'],template2)
            
            end
            
            if CreateVariant == 1
                
                outputfile1 = [subs{i} '_' tasks{j} '_' 'REST_OddSessions_Rand'];
                outputfile2 = [subs{i} '_' tasks{j} '_' 'REST_EvenSessions_Rand'];
                
            end
    
        elseif SubSample == 1 && RandSample == 1 && CortexOnly == 1
            
            if WriteDconn == 1
        
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' tasks{j} '_REST_AllSessions_Rand_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_REST_AllSessions_Rand_cortex'],template)
        
            end
            
            if CreateVariant == 1
                
                outputfile = [subs{i} '_' tasks{j} '_REST_AllSessions_Rand'];
                
            end
                
        elseif SubSample == 1 && RandSample == 1
            
            if WriteDconn == 1
        
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' tasks{j} '_REST_AllSessions_Rand'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_REST_AllSessions_Rand'],template)
        
            end
            
            if CreateVariant == 1
                
                outputfile = [subs{i} '_' tasks{j} '_REST_AllSessions_Rand'];
                
            end
 
        elseif SubSample == 1 && CortexOnly == 1
            
            if WriteDconn == 1
    
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' tasks{j} '_REST_AllSessions_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_REST_AllSessions_cortex'],template)
    
            end
            
            if CreateVariant == 1
                
                outputfile = [subs{i} '_' tasks{j} '_REST_AllSessions'];
                
            end
                
        elseif SubSample == 1
            
            if WriteDconn == 1
        
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' tasks{j} '_REST_AllSessions'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_REST_AllSessions'],template)
        
            end
            
            if CreateVariant == 1
                
                outputfile = [subs{i} '_' tasks{j} '_REST_AllSessions'];
                
            end
                
        elseif CortexOnly == 1
            
            if WriteDconn == 1
        
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' tasks{j} '_AllSessions_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_AllSessions_cortex'],template)
        
            end
            
            if CreateVariant == 1
                
                outputfile = [subs{i} '_' tasks{j} '_AllSessions'];
                
            end
                
        else
            
            if WriteDconn == 1
        
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' tasks{j} '_AllSessions'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_AllSessions'],template)
        
            end
            
            if CreateVariant == 1
                
                outputfile = [subs{i} '_' tasks{j} '_AllSessions'];
                
            end
                
        end
        
        if CreateVariant == 1 && SplitHalf == 1 && WriteDconn == 0
            
            disp(sprintf('Running Correlations: on data size %i by %i, %s', size(catData1,1), size(catData1,2), datestr(now)));
            template.data = paircorr_mod(catData1');
            
            catData1 = [];
            
            createSptlcorr_MSCdconns('/projects/b1081/Atlases', '120_allsubs_corr',1,'/projects/b1081/Brian_MSC/Analysis_Scripts_Replication/output_files/variant_maps',template.data, outputfile1)
            
            template.data = [];
            
            disp(sprintf('Running Correlations: on data size %i by %i, %s', size(catData2,1), size(catData2,2), datestr(now)));
            template2.data = paircorr_mod(catData2');
            
            catData2 = [];
            
            createSptlcorr_MSCdconns('/projects/b1081/Atlases', '120_allsubs_corr',1,'/projects/b1081/Brian_MSC/Analysis_Scripts_Replication/output_files/variant_maps',template2.data, outputfile2)
            
            template2.data = [];
        
        elseif CreateVariant == 1 && SplitHalf == 1 && ConcatenateSplitHalf == 0
            
            createSptlcorr_MSCdconns('/projects/b1081/Atlases', '120_allsubs_corr',1,'/projects/b1081/Brian_MSC/Analysis_Scripts_Replication/output_files/variant_maps',template.data, outputfile1)
            
            createSptlcorr_MSCdconns('/projects/b1081/Atlases', '120_allsubs_corr',1,'/projects/b1081/Brian_MSC/Analysis_Scripts_Replication/output_files/variant_maps',template2.data, outputfile2)
            
        elseif CreateVariant == 1
        
            createSptlcorr_MSCdconns('/projects/b1081/Atlases', '120_allsubs_corr',1,'/projects/b1081/Brian_MSC/Analysis_Scripts_Replication/output_files/variant_maps',template.data, outputfile)
            
        end
            
            clear template
            clear template2
            
           % end % End of session loop for saving each session separately
        %end % End of task data loop for saving each task separately
    end
end

disp(sprintf('Job Completed: %s', datestr(now)));