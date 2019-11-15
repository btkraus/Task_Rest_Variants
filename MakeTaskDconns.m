parpool('local', 28) %% Name of cluster profile for batch job, how many workers are allocated to perform this job

clear all

outdir = '/projects/b1081/Brian_MSC/dconn_task_files';%change this to your output directory
dataLocStem = '/MSC/TaskFC/';%directory where the data are located
subs = {'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC09','MSC10'};
tasks = {'motor','mem','mixed'};

CortexOnly = 1; %% Toggles whether to run correlations on cortex only
SplitHalf = 1;  %% Toggles whether to create a separate file for odd/even sessions
MatchData = 1; %% Toggles whether to match the amount of data per task as the lowest value within each split-half
RandSample = 1; %% Toggles whether to randomly sample data from each session
MatchAcrossSubs = 1;  %% Toggles whether to match the amount of data across subjects
MatchAcrossTasks = 0;  %% Toggles whether to match the amount of data across tasks
ConcatenateTasks = 1;   %% Toggles whether to concatenate data for all tasks
MakeDconn = 0;  %% Toggles whether to write a dconn
MakeVariantMap = 1; %% Toggles whether to write a variant map
ConcatenateSplitHalf = 0;  %% Toggles whether to concatenate split halves into one vector
SaveTimeseries = 1;     %% Save concatenated timeseries for subject


cd '/projects/b1081';   %% Change CD to root project directory

disp(sprintf('Job Submitted: %s', datestr(now)));
disp(sprintf('Job Started: %s', datestr(now)));

if CortexOnly == 1      %% Select correct number of voxels for template
    voxnum = 59412;
else
    voxnum = 65625;
end
%%
%%What are these variables? What is the difference between these and the
%%variables in the for loop below
if MatchData == 1 && MatchAcrossSubs == 1
   	memptsoddsum = [];
   	motorptsoddsum = [];
   	mixedptsoddsum = [];
   	memptsevensum = [];
   	motorptsevensum = [];
   	mixedptsevensum = [];

    for n=1:numel(subs)
        %change path to local directory where QC files are located
        load (['/projects/b1081/Brian_MSC/QC_files/' subs{n} '_QCFile.mat']);
        
        memptsodd = [];
        motorptsodd = [];
        mixedptsodd = [];
        memptseven = [];
        motorptseven = [];
        mixedptseven = [];
        
        %What is this for loop doing? Getting values for the variables
        %above from QC files, but what are these variables?
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
        
        %Why does MSC09 get dealt with differently? Not enough motor data?
        if strcmp(subs{n}, 'MSC09')
            motorptsoddsum = [motorptsoddsum; 9999];
        else
            motorptsoddsum = [motorptsoddsum; sum(motorptsodd)];
        end
        
        mixedptsoddsum = [mixedptsoddsum; sum(mixedptsodd)];
        memptsevensum = [memptsevensum; sum(memptseven)];
        
        if strcmp(subs{n}, 'MSC09')
            motorptsevensum = [motorptsevensum; 9999];        
        else
            motorptsevensum = [motorptsevensum; sum(motorptseven)];
        end
        
        mixedptsevensum = [mixedptsevensum; sum(mixedptseven)];
        
    end
    
    %%
    %what is this doing?
    if SplitHalf == 1
    
        mintasksamps = min(min([memptsoddsum motorptsoddsum mixedptsoddsum memptsevensum motorptsevensum mixedptsevensum]));
        
    elseif ConcatenateTasks == 1
        
        mintasksamps = min(min([(memptsoddsum + memptsevensum)  (motorptsoddsum + motorptsevensum)  (mixedptsoddsum + mixedptsevensum)]));
        
    end
    
end


for i=1:numel(subs)
    
    disp(sprintf('Creating dconn for subject %s: %s', subs{i}, datestr(now)));
  
    % Initialize cat data
    % what is cat data?
    if SplitHalf == 1
    	catData1 = [];
    	catData2 = [];
    else
    	catData = [];    
    end
    
    for j=1:length(tasks)
    	if MatchData == 1 && SplitHalf == 1
            %change path to local directory containing QC files
            load (['/projects/b1081/Brian_MSC/QC_files/' subs{i} '_QCFile.mat']);
            %why are we setting up these variables again?
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

            if strcmp(subs{i}, 'MSC09')  %% Removes motor task from consideration for MSC09, why??

                minsampspts = min([sum(memptsodd) sum(mixedptsodd) sum(memptseven) sum(mixedptseven)]);

            else

                minsampspts = min([sum(memptsodd) sum(mixedptsodd) sum(motorptsodd) sum(memptseven) sum(mixedptseven) sum(motorptseven)]);

            end
           
            meansamppts = floor(mintasksamps/5);
            disp(sprintf('For all subjects, the minimum number of sample points in a split-half is %i, with a mean of %i points per sesssion: %s', mintasksamps, meansamppts, datestr(now)));
            meansampptstempodd = meansamppts;
            meansampptstempeven = meansamppts;
            
            if strcmp(subs{i}, 'MSC09') && strcmp(tasks{j}, 'motor') && ConcatenateTasks == 1  %% Sets motor task data points for subject MSC09
                    
            	meansampptstempodd = floor(sum(motorptsodd)/5);
            	meansampptstempeven = floor(sum(motorptseven)/5);
                        
            elseif strcmp(subs{i}, 'MSC09') && ConcatenateTasks == 1
                        
            	meansampptstempodd = meansamppts + round((meansamppts - round(sum(motorptsodd)/5))/2);
            	meansampptstempeven = meansamppts + round((meansamppts - round(sum(motorptseven)/5))/2);
                        
            end
            
           %what are these?
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
            
        elseif MatchData == 1
            

                load (['/projects/b1081/Brian_MSC/QC_files/' subs{i} '_QCFile.mat']);
        
                mempts = [];
                motorpts = [];
                mixedpts = [];

                for u = 1:length(SubStruct)
                
                    mempts = [mempts; SubStruct(u).MemSampPts];
                    mixedpts = [mixedpts; SubStruct(u).MixedSampPts];
                    motorpts = [motorpts; SubStruct(u).MotorSampPts];

                end

                minsampspts = min([sum(mempts) sum(mixedpts) sum(motorpts)]);
                
            if MatchAcrossSubs == 0

                meansamppts = floor(minsampspts/5);
                
                disp(sprintf('For subject %s the minimum number of sample points in a split-half is %i, with a mean of %i points per sesssion: %s', subs{i}, minsampspts, meansamppts, datestr(now)));
            
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

            %disp(sprintf('For subject %s the minimum number of sample points in a split-half is %i, with a mean of %i points per sesssion per task for concatenation: %s', subs{i}, minsampspts, meansamppts/3, datestr(now)));
            
            if strcmp(subs{i}, 'MSC09') && strcmp(tasks{j}, 'motor') && ConcatenateTasks == 1  %% Sets motor task data points for subject MSC09
                
                if SplitHalf == 1
                    
                    meansampptstempodd = floor(sum(motorptsodd)/5);
                    meansampptstempeven = floor(sum(motorptseven)/5);
                    
                else
                    
                 	meansampptstemp = floor(sum(motorpts)/10);
                    
                end
                        
            elseif strcmp(subs{i}, 'MSC09') && ConcatenateTasks == 1
                
                if SplitHalf == 1
                        
                    meansampptstempodd = meansamppts + round((meansamppts - round(sum(motorptsodd)/5))/2);
                    meansampptstempeven = meansamppts + round((meansamppts - round(sum(motorptseven)/5))/2);
                    
                else
                    
                    meansampptstemp = meansamppts + round((meansamppts - round(sum(motorpts)/10))/2);
                    
                end
                        
            end
            
           remaindertotal = 0;
           
           notenoughdata = zeros(length(SubStruct),1);
           
            while true
                
            	remainder = 0;
                sampspersession = zeros(length(SubStruct),1);
                
                notenoughdatatotal = sum(notenoughdata);
            
                meansampptstemp = (meansampptstemp + round(remaindertotal/(10-notenoughdatatotal)));
            
                for v = 1:length(SubStruct)
                        
                	if strcmp(tasks{j}, 'mem')
                
                        	if SubStruct(v).MemSampPts < meansampptstemp && notenoughdata(v) == 0
                    
                                notenoughdata(v) = 1;
                            	remainder = (remainder + (meansampptstemp - SubStruct(v).MemSampPts));
                            	sampspersession(v) = SubStruct(v).MemSampPts;
                    
                            	disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                            elseif notenoughdata(v) == 0
                    
                            	sampspersession(v) = meansampptstemp;
                    
                                disp(sprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                                
                            else
                                
                                sampspersession(v) = SubStruct(v).MemSampPts;
                                
                                disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                            end
                                           
                    elseif strcmp(tasks{j}, 'mixed')
                            
                        	if SubStruct(v).MixedSampPts < meansampptstemp && notenoughdata(v) == 0
                    
                                notenoughdata(v) = 1;
                            	remainder = (remainder + (meansampptstemp - SubStruct(v).MixedSampPts));
                            	sampspersession(v) = SubStruct(v).MixedSampPts;
                    
                            	disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                            elseif notenoughdata(v) == 0
                    
                            	sampspersession(v) = meansampptstemp;
                    
                                disp(sprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                                
                            else
                                
                                sampspersession(v) = SubStruct(v).MixedSampPts;
                                
                                disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                            end
                                          
                  	elseif strcmp(tasks{j}, 'motor')
                            
                        	if SubStruct(v).MotorSampPts < meansampptstemp && notenoughdata(v) == 0
                    
                                notenoughdata(v) = 1;
                            	remainder = (remainder + (meansampptstemp - SubStruct(v).MotorSampPts));
                            	sampspersession(v) = SubStruct(v).MotorSampPts;
                    
                            	disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                            elseif notenoughdata(v) == 0
                    
                            	sampspersession(v) = meansampptstemp;
                    
                                disp(sprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                                
                            else
                                
                                sampspersession(v) = SubStruct(v).MotorSampPts;
                                
                                disp(sprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s', subs{i}, SubStruct(v).SessionFile, sampspersession(v), datestr(now)));
                    
                            end
   
                    end    
              	end
            
                if remainder == 0
                    
                    disp(sprintf('Data sampling calculated for subject %s for task %s: %s', subs{i}, tasks{j}, datestr(now)));
                
                    break
                    
                elseif (sum(notenoughdata) == 10) || (sum(notenoughdata) == 10 && remainder == 0)
                    
                    disp(sprintf('Not enough data is present for this task to have enough data for subject %s for task %s: %s', subs{i}, tasks{j}, datestr(now)));
                    
                    break
                
                else
                
                	remaindertotal = remainder;
                    
                    disp(sprintf('Remainder of %i calculated for subject %s for task %s and session %s: %s', remaindertotal, subs{i}, tasks{j}, SubStruct(v).SessionFile, datestr(now)));
                
                end
            end
            
        end
        
        disp(sprintf('Creating dconn for subject %s task %s: %s', subs{i}, tasks{j}, datestr(now)));
        
        % Load vcids
        MSCTaskdir = strcat(cd, dataLocStem);
        [~,vcids,~,~,~] = textread([MSCTaskdir '/' subs{i} '_' tasks{j} '_DATALIST.txt'],'%s%s%s%s%s');
        
        disp(sprintf('%i sessions found for subject %s task %s: %s', size(vcids,1), subs{i}, tasks{j}, datestr(now)));
        
        % Load tmasks
        if strcmp(subs{i}, 'MSC03') || strcmp(subs{i}, 'MSC10')
                    
        	MSCcondidir = strcat(strcat(MSCTaskdir, '/'), ['FCProc_' subs{i} '_' tasks{j} '_pass2_FDfilt/']);
                    
        else
                
            MSCcondidir = strcat(strcat(MSCTaskdir, '/'), ['FCProc_' subs{i} '_' tasks{j} '_pass2/']);
                    
        end
                
        load ([MSCcondidir 'condindices.mat']);

        % Load and concatentate data
        for k=1:length(vcids)
            
            disp(sprintf('Loading timeseries for session %s for subject %s task %s: %s', vcids{k}, subs{i}, tasks{j}, datestr(now)));
            
            
            MSCciftidir = strcat(MSCcondidir, 'cifti_timeseries_normalwall_native_freesurf');
            
            %if k==1
                %cd('cifti_timeseries_normalwall_native_freesurf')
            %end
            
            try
                data = ft_read_cifti_mod([MSCciftidir '/' vcids{k} '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii']);
                data = data.data;
                if strcmp(tasks{j},'motor')
                	tmask = TIndFin(k).AllMotor;
                elseif strcmp(tasks{j},'mem')
                	tmask = TIndFin(k).AllMem;
                else
                    tmask = TIndFin(k).AllGlass + TIndFin(k).AllSemantic;
                end
                
                if sum(tmask) == 0      % Skip file if no task task data
                    
                    data = [];
                    disp(sprintf('Vcid %s for subject %s has no usable data for task %s', vcids{k}, subs{i}, tasks{j}));
                    
                elseif MatchData == 0
                
                    data = data(1:voxnum,logical(tmask));
                    
                elseif MatchData == 1 && RandSample == 1
                    
                	tmasknum = sampspersession(k);
                        
                	samplepts = find(tmask == 1);
                        
                	tmask(samplepts) = 0;
                        
                	sampleselect = datasample(samplepts,tmasknum,'Replace',false);
                        
                    tmask(sampleselect) = 1;
                            
                 	disp(sprintf('%i sample points from file %s have been randomly selected for task %s, %s', sum(tmask), vcids{k}, tasks{j}, datestr(now)));
                            
                	data = data(1:voxnum,logical(tmask)); 
                    
                else
                    
                 	tmasknum = sampspersession(k);
                
                  	data = data(1:voxnum,logical(tmask));        %% Get good data from task
                 	data = data(:,1:tmasknum);          %% Get correct amount of points corresponding to task
                    
                end
                
                
                if SplitHalf == 1
                
                    disp(sprintf('Data is %i by %i and mod operation equals %i', size(data,1), size(data,2), mod(k,2)));
                   
                    if mod(k,2) == 1
                                
                    	disp(sprintf('Adding %i data points from session %s for subject %s to odd data', size(data,2), vcids{k}, subs{i}));
                        
                    	catData1 = [catData1 data];
                            
                    else
                                
                    	disp(sprintf('Adding %i data points from session %s for subject %s to even data', size(data,2), vcids{k}, subs{i}));
                                
                    	catData2 = [catData2 data];
                            
                    end
                         
                else
                        
                 	catData = [catData data];
                        
                end
                        
            catch ME
                if strcmp(ME.message,'Invalid file identifier.  Use fopen to generate a valid file identifier.')
                    continue
                end
            end
        end 
    end        %% Ending task loop here creates concatenated task data
    
    if ConcatenateSplitHalf == 1
            
    	catData = [catData1 catData2];
            
    end
    
    if SaveTimeseries == 1          %% Save concatenated FC processed timeseries
        
        if SplitHalf == 1
            
            timeseriestemplate = ft_read_cifti_mod('/projects/b1081/Brian_MSC/dconn_scripts/templates/MSC01_allses_mean_native_freesurf_vs_120sub_corr.dtseries.nii');
            timeseriestemplate.data = [];
            timeseriestemplate2 = timeseriestemplate;
        
            timeseriestemplate.data = catData1;
            timeseriestemplate2.data = catData2;
        
            disp(sprintf('Data timseries is size %i by %i, %s', size(catData1,1), size(catData1,2), datestr(now)));
        
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_FCProcessed_Timeseries_Matched_Odd_allTaskCat_cortex'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_FCProcessed_Timeseries_Matched_Odd_allTaskCat_cortex'],timeseriestemplate)
            
            disp(sprintf('Data timseries is size %i by %i, %s', size(catData2,1), size(catData2,2), datestr(now)));
        
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_FCProcessed_Timeseries_Matched_Even_allTaskCat_cortex'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_FCProcessed_Timeseries_Matched_Even_allTaskCat_cortex'],timeseriestemplate2)
            
        else
        
            timeseriestemplate = ft_read_cifti_mod('/projects/b1081/Brian_MSC/dconn_scripts/templates/MSC01_allses_mean_native_freesurf_vs_120sub_corr.dtseries.nii');
            timeseriestemplate.data = [];
        
            timeseriestemplate.data = catData;
        
            disp(sprintf('Data timseries is size %i by %i, %s', size(catData,1), size(catData,2), datestr(now)));
        
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_FCProcessed_Timeseries_Matched_allTaskCat_cortex'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_FCProcessed_Timeseries_Matched_allTaskCat_cortex'],timeseriestemplate)
        
        end
    end
    
    % Make and save rmat
    if CortexOnly == 1
        
    	disp('Loading template: MSC01_allses_mean_native_freesurf_vs_120sub_corr.dtseries.nii');
    	template = ft_read_cifti_mod('/projects/b1081/Brian_MSC/dconn_scripts/templates/MSC01_allses_mean_native_freesurf_vs_120sub_corr.dtseries.nii');    
        
    else
    
    	disp(sprintf('Reading %s: %s', [MSCciftidir '/' vcid], datestr(now)));
    	template = ft_read_cifti_mod([MSCciftidir '/' vcid]);
    
    end
    
    template.data = [];
    template.dimord = 'pos_pos';
    template.hdr.dim(6) = template.hdr.dim(7);
    template.hdr.intent_code = 3001;
    template.hdr.intent_name = 'ConnDense';
    
    if SplitHalf == 1
            
    	template2 = template;
            
    end
    
    if SplitHalf == 1 && MakeDconn == 1
        
    	disp(sprintf('Running Correlations: on data size %i by %i, %s', size(catData1,1), size(catData1,2), datestr(now)));
    	template.data = paircorr_mod(catData1');
        catData1 = [];
    	disp(sprintf('Running Correlations: on data size %i by %i, %s', size(catData2,1), size(catData2,2), datestr(now)));
     	template2.data = paircorr_mod(catData2');
        catData2 = [];
        
    elseif SplitHalf == 0 || ConcatenateSplitHalf == 1
    
        disp(sprintf('Running Correlations: %s', datestr(now)));
        template.data = paircorr_mod(catData');
        catData = [];
    
    end
    
    if SplitHalf == 1 && ConcatenateTasks == 1 && CortexOnly == 1 && MatchData == 1 && MatchAcrossSubs == 1 && RandSample == 1
        
        if MakeDconn == 1 
            
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_matcheddata_Final_allTaskCat_OddSessions_cortex'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_matcheddata_Final_allTaskCat_OddSessions_cortex'],template)
        
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_matcheddata_Final_allTaskCat_EvenSessions_cortex'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_matcheddata_Final_allTaskCat_EvenSessions_cortex'],template2)
        
        end
        
        if MakeVariantMap == 1
            
            outputfile1 = [subs{i} '_matcheddata_Final_allTaskCat_EvenSessions_cortex'];
            outputfile2 = [subs{i} '_matcheddata_Final_allTaskCat_OddSessions_cortex'];
            
        end
    
    elseif SplitHalf == 1 && ConcatenateTasks == 1 && CortexOnly == 1 && MatchData == 1 && RandSample == 1
        
        if MakeDconn == 1 
            
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_matcheddata_allTaskCat_OddSessions_cortex'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_matcheddata_allTaskCat_OddSessions_cortex'],template)
        
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_matcheddata_allTaskCat_EvenSessions_cortex'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_matcheddata_allTaskCat_EvenSessions_cortex'],template2)
        
        end
        
        if MakeVariantMap == 1
            
            outputfile1 = [subs{i} '_matcheddata_allTaskCat_EvenSessions_cortex'];
            outputfile2 = [subs{i} '_matcheddata_allTaskCat_OddSessions_cortex'];
            
        end
        
    elseif SplitHalf == 1 && ConcatenateTasks == 1 && CortexOnly == 1 && MatchData == 1 && numel(subs) == 8
        
        if MakeDconn == 1 
            
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_matcheddata_NoMSC09_Consec_allTaskCat_OddSessions_cortex'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_matcheddata_NoMSC09_Consec_allTaskCat_OddSessions_cortex'],template)
        
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_matcheddata_NoMSC09_Consec_allTaskCat_EvenSessions_cortex'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_matcheddata_NoMSC09_Consec_allTaskCat_EvenSessions_cortex'],template2)
        
        end
        
        if MakeVariantMap == 1
            
            outputfile1 = [subs{i} '_matcheddata_NoMSC09_Consec_allTaskCat_EvenSessions_cortex'];
            outputfile2 = [subs{i} '_matcheddata_NoMSC09_Consec_allTaskCat_OddSessions_cortex'];
            
        end
        
    elseif SplitHalf == 1 && ConcatenateTasks == 1 && CortexOnly == 1 && MatchData == 1
        
        if MakeDconn == 1 
            
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_matcheddata_Consec_allTaskCat_OddSessions_cortex'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_matcheddata_Consec_allTaskCat_OddSessions_cortex'],template)
        
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_matcheddata_Consec_allTaskCat_EvenSessions_cortex'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_matcheddata_Consec_allTaskCat_EvenSessions_cortex'],template2)
        
        end
        
        if MakeVariantMap == 1
            
            outputfile1 = [subs{i} '_matcheddata_Consec_allTaskCat_EvenSessions_cortex'];
            outputfile2 = [subs{i} '_matcheddata_Consec_allTaskCat_OddSessions_cortex'];
            
        end
        
    elseif SplitHalf == 1 && ConcatenateTasks == 1 && CortexOnly == 1 && MatchData == 1 && ConcatenateSplitHalf == 1 && RandSample == 1
        
        if MakeDconn == 1 
        
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_matcheddata_allTaskCat_cortex_Rand'], datestr(now)));
            ft_write_cifti_mod([outdir '/' subs{i} '_matcheddata_allTaskCat_cortex_Rand'],template)
        
        end
        
        if MakeVariantMap == 1
            
            outputfile = [subs{i} '_matcheddata_allTaskCat_cortex_Rand'];
            
        end
    
    elseif SplitHalf == 1 && ConcatenateTasks == 1 && CortexOnly == 1 && ConcatenateSplitHalf == 0 && RandSample == 1
        
        if MakeDconn == 1 
            
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_allTaskCat_OddSessions_cortex'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_allTaskCat_OddSessions_cortex'],template)
        
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_allTaskCat_EvenSessions_cortex'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_allTaskCat_EvenSessions_cortex'],template2)
        
        end
        
        if MakeVariantMap == 1
            
            outputfile1 = [subs{i} '_allTaskCat_EvenSessions_cortex'];
            outputfile2 = [subs{i} '_allTaskCat_OddSessions_cortex'];
            
        end
    
    elseif ConcatenateTasks == 1 && CortexOnly == 1 && RandSample == 1
    
        if MakeDconn == 1 
        
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_allTaskCat_cortex'], datestr(now)));
            ft_write_cifti_mod([outdir '/' subs{i} '_allTaskCat_cortex'],template)
        
        end
        
        if MakeVariantMap == 1
            
            outputfile = [subs{i} '_allTaskCat_cortex'];
            
        end
        
    elseif ConcatenateTasks == 1 && RandSample == 1
        
        if MakeDconn == 1 
        
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_allTaskCat_subcort'], datestr(now)));
            ft_write_cifti_mod([outdir '/' subs{i} '_allTaskCat_subcort'],template)
        
        end
        
        if MakeVariantMap == 1
            
            outputfile = [subs{i} '_allTaskCat_subcort'];
            
        end
        
    elseif SplitHalf == 1 && CortexOnly == 1 && MatchData == 1 && RandSample == 1
        
        if MakeDconn == 1 
            
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_' tasks{j} '_matcheddata_OddSessions_cortex'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_matcheddata_OddSessions_cortex'],template)
        
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_' tasks{j} '_matcheddata_EvenSessions_cortex'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_matcheddata_EvenSessions_cortex'],template2)
        
        end
        
        if MakeVariantMap == 1
            
            outputfile1 = [subs{i} '_' tasks{j} '_matcheddata_EvenSessions_cortex'];
            outputfile2 = [subs{i} '_' tasks{j} '_matcheddata_OddSessions_cortex'];
            
        end
        
    elseif SplitHalf == 1 && CortexOnly == 1 && MatchData == 1 && numel(subs) == 8
        
        if MakeDconn == 1 
            
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_' tasks{j} '_matcheddata_NoMSC09_Consec_OddSessions_cortex'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_matcheddata_NoMSC09_Consec_OddSessions_cortex'],template)
        
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_' tasks{j} '_matcheddata_NoMSC09_Consec_EvenSessions_cortex'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_matcheddata_NoMSC09_Consec_EvenSessions_cortex'],template2)
        
        end
        
        if MakeVariantMap == 1
            
            outputfile1 = [subs{i} '_' tasks{j} '_matcheddata_NoMSC09_Consec_EvenSessions_cortex'];
            outputfile2 = [subs{i} '_' tasks{j} '_matcheddata_NoMSC09_Consec_OddSessions_cortex'];
            
        end
        
    elseif SplitHalf == 1 && CortexOnly == 1 && MatchData == 1
        
        if MakeDconn == 1 
            
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_' tasks{j} '_matcheddata_Consec_OddSessions_cortex'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_matcheddata_Consec_OddSessions_cortex'],template)
        
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_' tasks{j} '_matcheddata_Consec_EvenSessions_cortex'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_matcheddata_Consec_EvenSessions_cortex'],template2)
        
        end
        
        if MakeVariantMap == 1
            
            outputfile1 = [subs{i} '_' tasks{j} '_matcheddata_Consec_EvenSessions_cortex'];
            outputfile2 = [subs{i} '_' tasks{j} '_matcheddata_Consec_OddSessions_cortex'];
            
        end
        
    elseif SplitHalf == 1 && CortexOnly == 1 && RandSample == 1
        
        if MakeDconn == 1 
            
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_' tasks{j} '_OddSessions_cortex'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_OddSessions_cortex'],template)
        
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_' tasks{j} '_EvenSessions_cortex'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_EvenSessions_cortex'],template2)
        
        end
        
        if MakeVariantMap == 1
            
            outputfile1 = [subs{i} '_' tasks{j} '_EvenSessions_cortex'];
            outputfile2 = [subs{i} '_' tasks{j} '_OddSessions_cortex'];
            
        end
        
    elseif SplitHalf == 1 && RandSample == 1
        
        if MakeDconn == 1 
            
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_' tasks{j} '_OddSessions_subcort'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_OddSessions_subcort'],template)
        
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_' tasks{j} '_EvenSessions_subcort'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_EvenSessions_subcort'],template2)
        
        end
        
        if MakeVariantMap == 1
            
            outputfile1 = [subs{i} '_' tasks{j} '_EvenSessions_subcort'];
            outputfile2 = [subs{i} '_' tasks{j} '_OddSessions_subcort'];
            
        end
        
    elseif CortexOnly == 1 && RandSample == 1
        
        if MakeDconn == 1 
        
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_' tasks{j} '_AllSessions_cortex'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_AllSessions_cortex'],template)
            
        end
        
        if MakeVariantMap == 1
            
            outputfile = [subs{i} '_' tasks{j} '_AllSessions_cortex'];
            
        end
        
    elseif RandSample == 1
        
        if MakeDconn == 1 
        
            disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_' tasks{j} '_AllSessions_subcort'], datestr(now)));
            ft_write_cifti_mod(['/' outdir '/' subs{i} '_' tasks{j} '_AllSessions_subcort'],template)
            
        end
        
        if MakeVariantMap == 1
            
            outputfile = [subs{i} '_' tasks{j} '_AllSessions_subcort'];
            
        end
    end
        
        if MakeVariantMap == 1 && SplitHalf == 1 && MakeDconn == 0
            
            disp(sprintf('Running Correlations: on data size %i by %i, %s', size(catData1,1), size(catData1,2), datestr(now)));
            template.data = paircorr_mod(catData1');
            
            catData1 = [];
            
            createSptlcorr_MSCdconns('/projects/b1081/Atlases', '120_allsubs_corr',1,'/home/btk2142/output_files/variant_maps',template.data, outputfile1)
            
            template.data = [];
            
            disp(sprintf('Running Correlations: on data size %i by %i, %s', size(catData2,1), size(catData2,2), datestr(now)));
            template2.data = paircorr_mod(catData2');
            
            catData2 = [];
            
            createSptlcorr_MSCdconns('/projects/b1081/Atlases', '120_allsubs_corr',1,'/home/btk2142/output_files/variant_maps',template2.data, outputfile2)
            
            template2.data = [];
        
        elseif MakeVariantMap == 1 && SplitHalf == 1
            
            createSptlcorr_MSCdconns('/projects/b1081/Atlases', '120_allsubs_corr',1,'/home/btk2142/output_files/variant_maps',template.data, outputfile1)
            
            createSptlcorr_MSCdconns('/projects/b1081/Atlases', '120_allsubs_corr',1,'/home/btk2142/output_files/variant_maps',template2.data, outputfile2)
            
        elseif MakeVariantMap == 1
        
            createSptlcorr_MSCdconns('/projects/b1081/Atlases', '120_allsubs_corr',1,'/home/btk2142/output_files/variant_maps',template.data, outputfile)
            
        end
            
        clear template
        clear template2
    
    %end        %% Ending task loop here creates file for each task
end

disp(sprintf('Job Completed: %s', datestr(now)));
