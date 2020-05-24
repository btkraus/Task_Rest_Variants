%This script makes Dconns (if the MakeDconn variable = 1) and variant files
% (if VariantMap = 1)for rest data of each subject's even and odd sessions. 
% QC files need to be na med [subject]_QCFile.mat, ex. 'MSC01_QCFile.mat'

%parpool('local', 28)     %% Name of cluster profile for batch job (how many workers/cores you need to run this job)

clear all

disp(sprintf('Job Submitted: %s', datestr(now)));

%% Paths
outdir = '/Users/diana/Desktop/';%specify output directory
dataLocStem = '/Users/diana/Box/Quest_Backup/MSC/TaskFC/'; %specify location of data
QC_path = '/Users/diana/Box/Quest_Backup/member_directories/bkraus/Brian_MSC/QC_files/'; %specify location of QC files
groupavg_path = '/Users/diana/Box/Quest_Backup/Atlases'; %specify location of group average map
groupavg_name = '120_allsubs_corr'; %specify name of group average map
rest_path = '/Users/diana/Box/Quest_Backup/MSC/MSCdata_v1/'; %specify location of rest data
output_path = '/Users/diana/Desktop/variant_maps'; %%%%%% WHAT IS THE DIFFERENCE BETWEEN THIS AND outdir?
template_path = '/Users/diana/Box/Quest_Backup/member_directories/bkraus/Brian_MSC/dconn_scripts/templates/MSC01_allses_mean_native_freesurf_vs_120sub_corr.dtseries.nii'; %location of correlation map (???) to be used as template
cd '/Users/diana/Box/Quest_Backup';   %% Change CD to root project directory
%% Options - deleted RestOnly (default is 0), SubSample (default is 1), ConcatenateSplitHalf (default is 0)
SplitHalf = 1;  %% Toggles whether to calculate rest for split-half of sessions
RandSample = 1; %% Toggles random subsampling of tmask on/off
MatchData = 1; %% Toggles whether to match the amount of data per task as the lowest value within each split-half
MatchAcrossTasks = 0; %% Toggles whether to match the amount of data across tasks (for concatenated task data)
MatchAcrossSubs = 1;  %% Toggles whether to match the amount of data across subjects
ConcatenateTasksMatch = 1; %% Toggles whether to calculate matched task data for concatenated tasks
ConcatenateTaskData = 1; %% Toggles whether to concatenate task data within subjects
ConcatenateSessionData = 1; %% Toggles whether to concatenate session data within subjects
WriteDconn = 0; %% Toggles whether to write dconn to disk
CreateVariant = 1; %% Toggles whether to create a variant map from dconn
SaveTimeseries = 1;     %% Save concatenated timeseries for subject
%% Variables
subs = {'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC09','MSC10'};
tasks = {'motor','mem','mixed'};
times = 1;
voxnum = 59412;

if ConcatenateTasksMatch == 1    
    taskbackup = tasks;
end

%% Start Analysis 
disp(sprintf('Job Started: %s', datestr(now)));
%% Matches data points for all tasks across subjects
[mintasksamps] = mintasksamps(QC_path, subs, SplitHalf, ConcatenateTaskData, MatchData, MatchAcrossSubs);
%% Main for-loop: matches number of data points to sample from rest to task, makes Dconns
for i=1:numel(subs)
    
    if ConcatenateTasksMatch == 1    
        tasks = taskbackup;  %% Resets task loop for each subject        
    end
    
    % establishing number of sample points to match for each session
    if ConcatenateTaskData == 1
        if SplitHalf == 1
            sampspersession = floor(repmat(mintasksamps,[1 10])./5);
        else
            sampspersession = floor(repmat(mintasksamps,[1 10])./10);
        end
        load([QC_path subs{i} '_QCFile.mat']);
    else              
        if MatchAcrossSubs == 1                
            load([QC_path subs{i} '_QCFile.mat']);                
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
            if strcmp(subs{i}, 'MSC09') && strcmp(taskbackup{h}, 'motor') && ConcatenateTaskData == 1  
                if SplitHalf == 1
                    sampspersessiontask(1:2:end) = floor(sum(motorptsodd)/5);
                    sampspersessiontask(2:2:end) = floor(sum(motorptseven)/5);
                else
                    sampspersessiontask = floor(sum(motorpts)/10);
                end
            elseif strcmp(subs{i}, 'MSC09') && ConcatenateTaskData == 1
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

                                disp(sprintf('For subject %s and session %s there are not enough sample points for the task, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            elseif notenoughdatatask(c) == 0
                                tasksampspersession(c) = tasksampspersession(c);

                                disp(sprintf('For subject %s and session %s there are enough sample points for the task, sampling %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            else
                                tasksampspersession(c) = SubStruct(c).MemSampPts;
                                
                                disp(sprintf('For subject %s and session %s there are not enough sample points for the task, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            end
                        elseif strcmp(taskbackup{h},'mixed')
                            if SubStruct(c).MixedSampPts <  tasksampspersession(c) && notenoughdatatask(c) == 0
                                notenoughdatatask(c) = 1;
                                remainderoddtask = (remainderoddtask + (tasksampspersession(c) - SubStruct(c).MixedSampPts));
                                tasksampspersession(c) = SubStruct(c).MixedSampPts;

                                disp(sprintf('For subject %s and session %s there are not enough sample points for the task, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            elseif notenoughdatatask(c) == 0
                                tasksampspersession(c) = tasksampspersession(c);

                                disp(sprintf('For subject %s and session %s there are enough sample points for the task, sampling %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            else
                                tasksampspersession(c) = SubStruct(c).MixedSampPts

                                disp(sprintf('For subject %s and session %s there are not enough sample points for the task, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            end
                        elseif strcmp(taskbackup{h},'motor')
                            if SubStruct(c).MotorSampPts <  tasksampspersession(c) && notenoughdatatask(c) == 0
                                notenoughdatatask(c) = 1;
                                remainderoddtask = (remainderoddtask + (tasksampspersession(c) - SubStruct(c).MotorSampPts));
                                tasksampspersession(c) = SubStruct(c).MotorSampPts;

                                disp(sprintf('For subject %s and session %s there are not enough sample points for the task, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            elseif notenoughdatatask(c) == 0
                                tasksampspersession(c) = tasksampspersession(c);

                                disp(sprintf('For subject %s and session %s there are enough sample points for the task, sampling %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            else
                                tasksampspersession(c) = SubStruct(c).MotorSampPts;

                                disp(sprintf('For subject %s and session %s there are not enough sample points for the task, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            end
                        end
                    else
                        if strcmp(taskbackup{h},'mem')
                            if SubStruct(c).MemSampPts < tasksampspersession(c) && notenoughdatatask(c) == 0
                                notenoughdatatask(c) = 1;
                                remaindereventask = (remaindereventask + (tasksampspersession(c) - SubStruct(c).MemSampPts));
                                tasksampspersession(c) = SubStruct(c).MemSampPts;

                                disp(sprintf('For subject %s and session %s there are not enough sample points for the task, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            elseif notenoughdatatask(c) == 0
                                tasksampspersession(c) = tasksampspersession(c);

                                disp(sprintf('For subject %s and session %s there are enough sample points for the task, sampling %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            else
                                tasksampspersession(c) = SubStruct(c).MemSampPts;

                                disp(sprintf('For subject %s and session %s there are not enough sample points in the task data, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            end
                        elseif strcmp(taskbackup{h},'mixed')
                            if SubStruct(c).MixedSampPts < tasksampspersession(c) && notenoughdatatask(c) == 0
                                notenoughdatatask(c) = 1;
                                remaindereventask = (remaindereventask + (tasksampspersession(c) - SubStruct(c).MemSampPts));
                                tasksampspersession(c) = SubStruct(c).MixedSampPts;

                                disp(sprintf('For subject %s and session %s there are not enough sample points for the task, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            elseif notenoughdatatask(c) == 0
                                tasksampspersession(c) = tasksampspersession(c);

                                disp(sprintf('For subject %s and session %s there are enough sample points for the task, sampling %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            else
                                tasksampspersession(c) = SubStruct(c).MixedSampPts;

                                disp(sprintf('For subject %s and session %s there are not enough sample points in the task data, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            end
                        elseif strcmp(taskbackup{h},'motor')
                            if SubStruct(c).MixedSampPts < tasksampspersession(c) && notenoughdatatask(c) == 0
                                notenoughdatatask(c) = 1;
                                remaindereventask = (remaindereventask + (tasksampspersession(c) - SubStruct(c).MemSampPts));
                                tasksampspersession(c) = SubStruct(c).MixedSampPts;

                                disp(sprintf('For subject %s and session %s there are not enough sample points for the task, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            elseif notenoughdatatask(c) == 0
                                tasksampspersession(c) = tasksampspersession(c);

                                disp(sprintf('For subject %s and session %s there are enough sample points for the task, sampling %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            else
                                tasksampspersession(c) = SubStruct(c).MixedSampPts;

                                disp(sprintf('For subject %s and session %s there are not enough sample points in the task data, sampling only %i sample points from this session: %s', subs{i}, SubStruct(c).SessionFile, tasksampspersession(c), datestr(now)));
                            end
                        end
                    end
                end

                if remainderoddtask == 0 && remaindereventask == 0
                    disp(sprintf('Data sampling calculated for subject %s for task %s: %s', subs{i}, taskbackup{h}, datestr(now)));
                    break
                elseif (sum(notenoughdatatask) == 10) || (sum(notenoughdataeventask) == 5 && remainderoddtask == 0) || (sum(notenoughdataoddtask) == 5 && remaindereventask == 0)
                    disp(sprintf('Not enough data is present for this task to have enough data for subject %s for task %s: %s', subs{i}, datestr(now)));
                    break
                else
                    remaindertotaloddtask = remainderoddtask;
                    remaindertotaleventask = remaindereventask;
                    disp(sprintf('Remainder of %i calculated for odd files for subject %s: %s', remainderoddtask, subs{i}, datestr(now)));
                    disp(sprintf('Remainder of %i calculated for even files for subject %s: %s', remaindereventask, subs{i}, datestr(now)));
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
        load ([QC_path subs{i} '_QCFile.mat']);

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

        % Get final count for task data per session
        meansamppts = round(minsampspts/5);

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

    %% Make Dconns
    for l=1:numel(times)        
    
        disp(sprintf('Creating dconn for subject %s: %s', subs{i}, datestr(now)));

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
                load ([QC_path subs{i} '_QCFile.mat']);
        
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
        
                if strcmp(subs{i}, 'MSC09')  %% Removes motor task from consideration for MSC09                
                    minsampspts = min([sum(memptsodd) sum(mixedptsodd) sum(memptseven) sum(mixedptseven)]);                
                else        
                    minsampspts = min([sum(memptsodd) sum(mixedptsodd) sum(motorptsodd) sum(memptseven) sum(mixedptseven) sum(motorptseven)]);                
                end
                 %% IS THIS DOING THE SAME THING THAT IT WAS DOING ABOVE??? DETERMINING NUMBER OF SAMPLE POINTS, MATCHING TASK TO REST ETC?
                % gets mean number of sample points per task per session
                if MatchAcrossSubs == 0                    
                    if MatchAcrossTasks == 1
                        meansamppts = floor(minsampspts/(5*numel(tasks)));                
                        disp(sprintf('For subject %s the minimum number of sample points in a split-half is %i, with %i points per task per spit-half, with a mean of %i points per rest sesssion: %s', subs{i}, minsampspts, meansamppts, meansamppts*numel(tasks), datestr(now)));                        
                    else
                        meansamppts = floor(minsampspts/5);                
                        disp(sprintf('For subject %s the minimum number of sample points in a split-half is %i, with a mean of %i points per rest sesssion: %s', subs{i}, minsampspts, meansamppts, datestr(now)));                        
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
        
            if MatchData == 0 && ConcatenateTaskData == 0
        
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
                MSCciftidir = [dataLocStem subs{i} '_mem_pass2/cifti_timeseries_normalwall_native_freesurf'];        
            end
        
            % Load rest vcids
            restdir = dir([rest_path subs{i} '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_normalwall_native_freesurf']);
            vcidlist = restdir(~cellfun(@isempty,strfind({restdir.name},'dtseries')));
            vcidlist = vcidlist(~cellfun(@isempty,strfind({vcidlist.name},'vc')));
            
            disp(sprintf('%i resting sessions found for subject %s: %s', size(vcidlist,1), subs{i}, datestr(now)));
        
            % Load rest tmasks
            load([rest_path subs{i} '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/QC.mat'])
        
            % Load and concatentate data
            for k=1:length(vcidlist)           
                vcid = vcidlist(k).name;            
                disp(sprintf('Loading timeseries for session %s for subject %s task %s: %s', vcid, subs{i}, tasks{j}, datestr(now)));
            
                if MatchData == 0 && ConcatenateTaskData == 0            
                    MSCciftidir = strcat(MSCcondidir, 'cifti_timeseries_normalwall_native_freesurf');            
                end

                tmask = 0;  % Create empty tmask variable
            
                try
                    data = ft_read_cifti_mod([rest_path subs{i} '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_normalwall_native_freesurf/' vcid]);
                    disp(sprintf('Loading data size %i by %i, %s', size(data.data,1), size(data.data,2), datestr(now)));
                    data = data.data;
                
                    if MatchData == 0 && ConcatenateTaskData == 0                
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
                
                    if sum(tmask) == 0 && MatchData == 0 && ConcatenateTaskData == 0    % Skip file if no task task data                    
                        data = [];
                        disp(sprintf('Vcid %s for subject %s has no usable data for task %s', vcid, subs{i}, tasks{j}));                    
                    else                    % Else get first xth points in resting file corresponding to number of task data points                            
                        if RandSample == 1 && MatchData == 1    %% Subsample (rest) data for the desired number of voxels as set by the lowest split-half for a random subset of voxels                            
                            tmasknum = restsampspersession(k);                        
                            samplepts = find(resttmask == 1);                        
                            resttmask(samplepts) = 0;                        
                            sampleselect = datasample(samplepts,tmasknum,'Replace',false);                        
                            resttmask(sampleselect) = 1;                            
                            disp(sprintf('%i sample points from rest file %s have been randomly selected for task %s, %s', sum(resttmask), vcid, tasks{j}, datestr(now)));                            
                            data = data(1:voxnum,logical(resttmask));                            
                        elseif RandSample == 1 && ConcatenateTaskData == 1                            
                            tmasknum = restsampspersession(k);                        
                            samplepts = find(resttmask == 1);                        
                            resttmask(samplepts) = 0;                        
                            sampleselect = datasample(samplepts,tmasknum,'Replace',false);                        
                            resttmask(sampleselect) = 1;                            
                            disp(sprintf('%i sample points from rest file %s have been randomly selected for task %s, %s', sum(resttmask), vcid, tasks{j}, datestr(now)));                            
                            data = data(1:voxnum,logical(resttmask));                    
                        elseif RandSample == 1       %% Subsample (rest) data for the desired number of voxels for a random subset of voxels                        
                            tmasknum = restsampspersession(k);                        
                            samplepts = find(resttmask == 1);                        
                            resttmask(samplepts) = 0;                        
                            sampleselect = datasample(samplepts,tmasknum,'Replace',false);                        
                            resttmask(sampleselect) = 1;                            
                            disp(sprintf('%i sample points from rest file %s have been randomly selected for task %s, %s', sum(resttmask), vcid, tasks{j}, datestr(now)));                            
                            data = data(1:voxnum,logical(resttmask));                            
                        else     %% Subsample data from the first voxel to the desired number of voxels (in this case task number per session)                    
                            tmasknum = restsampspersession(k);                
                            data = data(1:voxnum,logical(resttmask));        %% Get good data from rest
                            data = data(:,1:tmasknum);          %% Get correct amount of points corresponding to task                    
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
        
        if SaveTimeseries == 1          %% Save concatenated FC processed timeseries            
            if SplitHalf == 1                
                timeseriestemplate = ft_read_cifti_mod(template_path);
                timeseriestemplate.data = [];
                timeseriestemplate2 = timeseriestemplate;                        
                timeseriestemplate.data = catData1;
                timeseriestemplate2.data = catData2;            
                disp(sprintf('Data timseries is size %i by %i, %s', size(catData1,1), size(catData1,2), datestr(now)));        
                disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_FCProcessed_Timeseries_Matched_Odd_REST_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_FCProcessed_Timeseries_Matched_REST_Odd_cortex'],timeseriestemplate)                
                disp(sprintf('Data timseries is size %i by %i, %s', size(catData2,1), size(catData2,2), datestr(now)));        
                disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_FCProcessed_Timeseries_Matched_Even_REST_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_FCProcessed_Timeseries_Matched_Even_REST_cortex'],timeseriestemplate2)                
            else        
                timeseriestemplate = ft_read_cifti_mod(template_path);
                timeseriestemplate.data = [];        
                timeseriestemplate.data = catData;            
                disp(sprintf('Data timseries is size %i by %i, %s', size(catData,1), size(catData,2), datestr(now)));        
                disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_FCProcessed_Timeseries_Matched_REST_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_FCProcessed_Timeseries_Matched_REST_cortex'],timeseriestemplate)        
            end
        end
    
        % Make and save rmat
    
        disp('Loading template: MSC01_allses_mean_native_freesurf_vs_120sub_corr.dtseries.nii');
        template = ft_read_cifti_mod(template_path);    
           
        template.data = [];
        template.dimord = 'pos_pos';
        template.hdr.dim(7) = voxnum;
        template.hdr.dim(6) = template.hdr.dim(7);
        template.hdr.intent_code = 3001;
        template.hdr.intent_name = 'ConnDense';
        
        if SplitHalf == 1            
            template2 = template;            
        end        
        
        %% Not doing concatenate split half on this one either? 
        if SplitHalf == 1 && WriteDconn == 1            
            disp(sprintf('Running Correlations: on data size %i by %i, %s', size(catData1,1), size(catData1,2), datestr(now)));
            template.data = paircorr_mod(catData1');
            clear catData1
            disp(sprintf('Running Correlations: on data size %i by %i, %s', size(catData2,1), size(catData2,2), datestr(now)));
            template2.data = paircorr_mod(catData2');
            clear catData2            
        elseif SplitHalf == 0        
            disp(sprintf('Running Correlations: on data size %i by %i, %s', size(catData,1), size(catData,2), datestr(now)));
            template.data = paircorr_mod(catData');            
        end
        
        if SplitHalf == 1 && WriteDconn == 1            
            catData1 = [];
            catData2 = [];            
        elseif SplitHalf == 0        
        	catData = [];   % Clear cat data            
        end
       
        
        if SplitHalf ==1 && MatchAcrossSubs == 1 && MatchData == 1 && RandSample == 0            
            if WriteDconn == 1            
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' tasks{j} '_matcheddata_' 'REST_OddSessions_Consec_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_matcheddata_' 'REST_OddSessions_Consec_cortex'],template)
            
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' tasks{j} '_matcheddata_' 'REST_EvenSessions_Consec_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_matcheddata_' 'REST_EvenSessions_Consec_cortex'],template2)            
            end
            
            if CreateVariant == 1                
                outputfile1 = [subs{i} '_' tasks{j} '_matcheddata_' 'REST_OddSessions_Consec'];
                outputfile2 = [subs{i} '_' tasks{j} '_matcheddata_' 'REST_EvenSessions_Consec'];                
            end            
        elseif SplitHalf ==1 && MatchAcrossSubs == 1 && MatchData == 1            
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
        elseif SplitHalf == 1 && RandSample == 1 && MatchData == 1            
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
        elseif SplitHalf == 1 && RandSample == 1 && ConcatenateTaskData == 1            
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
        elseif SplitHalf == 1 && RandSample == 1            
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
        elseif SplitHalf == 1 && MatchData == 1            
            if WriteDconn == 1            
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' tasks{j} '_matcheddata_' 'REST_OddSessions_Consec_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_matcheddata_' 'REST_OddSessions_Consec_cortex'],template)
            
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' tasks{j} '_matcheddata_' 'REST_EvenSessions_Consec_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_matcheddata_' 'REST_EvenSessions_Consec_cortex'],template2)            
            end
            
            if CreateVariant == 1                
                outputfile1 = [subs{i} '_' tasks{j} '_matcheddata_' 'REST_OddSessions_Consec'];
                outputfile2 = [subs{i} '_' tasks{j} '_matcheddata_' 'REST_EvenSessions_Consec'];                
            end                    
        elseif RandSample == 1            
            if WriteDconn == 1        
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' tasks{j} '_REST_AllSessions_Rand'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_REST_AllSessions_Rand'],template)        
            end
            
            if CreateVariant == 1                
                outputfile = [subs{i} '_' tasks{j} '_REST_AllSessions_Rand'];                
            end 
        else           
            if WriteDconn == 1    
                disp(sprintf('Writing %s: %s', ['/' outdir '/' subs{i} '_' tasks{j} '_REST_AllSessions_cortex'], datestr(now)));
                ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_REST_AllSessions_cortex'],template)    
            end
            
            if CreateVariant == 1                
                outputfile = [subs{i} '_' tasks{j} '_REST_AllSessions'];                
            end                        
        end
        
        %% Creates Spatial Correlation Maps
        % createSptlcorr_MSCdconns(groupAvgLoc,groupAvgName,cortexOnly,outputdir,dconnData,
        % outputname)
        if CreateVariant == 1 && SplitHalf == 1 && WriteDconn == 0            
            disp(sprintf('Running Correlations: on data size %i by %i, %s', size(catData1,1), size(catData1,2), datestr(now)));
            template.data = paircorr_mod(catData1');            
            catData1 = [];            
            createSptlcorr_MSCdconns(groupavg_path, groupavg_name,1,output_path,template.data, outputfile1)            
            template.data = [];            
            disp(sprintf('Running Correlations: on data size %i by %i, %s', size(catData2,1), size(catData2,2), datestr(now)));
            template2.data = paircorr_mod(catData2');
            
            catData2 = [];
            
            createSptlcorr_MSCdconns(groupavg_path, groupavg_name ,1,output_path,template2.data, outputfile2)
            
            template2.data = [];
        
        elseif CreateVariant == 1 && SplitHalf == 1
            
            createSptlcorr_MSCdconns(groupavg_path, groupavg_name,1,output_path,template.data, outputfile1)
            
            createSptlcorr_MSCdconns(groupavg_path, groupavg_name,1,output_path,template2.data, outputfile2)
            
        elseif CreateVariant == 1
        
            createSptlcorr_MSCdconns(groupavg_path, groupavg_name,1,output_path,template.data, outputfile)
            
        end
            
            clear template
            clear template2
            
           % end % End of session loop for saving each session separately
        %end % End of task data loop for saving each task separately
        end
end

disp(sprintf('Job Completed: %s', datestr(now)));