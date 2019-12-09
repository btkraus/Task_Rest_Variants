%parpool('local', 20)     %% Name of cluster profile for batch job

%/davta/vnil-bluearc/GMT/Caterina/TaskFC/FCProc_MSC01_motor_passvv2/
%/data/nil-bluearc/GMT/Caterina/TaskFC/FCProc_MSC01_mem_pass2/
%/data/nil-bluearc/GMT/Caterina/TaskFC/FCProc_MSC01_mixed_pass2/
clear all

outdir = '/projects/b1081/Brian_MSC/QC_files';
dataLocStem = '/MSC/TaskFC/';
MSC = 0;  %% Toggles whether to create a struct for MSC data
WashU = 1;  %% Toggles whether to create a struct for Poldrome WashU data
Texas = 0;  %% Toggles whether to create a struct for Poldrome Texas data


if MSC == 1

    subs = {'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC09','MSC10'};
    %subs = {'MSC05'};
    tasks = {'motor','mem','mixed'};
    %tasks = {'mem'};

else
    
    subs = 1;
    
end


if ~exist(outdir , 'dir')
mkdir(outdir)
end

cd '/projects/b1081';   %% Change CD to root project directory


if MSC == 1

    for i=1:numel(subs)
    
    
    
    OddCount = 0;
    EvenCount = 0;
    
    % Create Struct

    SubStruct = struct('SessionFile',{},'MemExists',[],'MemSampPts',[], 'MixedExists',[],'MixedSampPts',[], 'MotorExists',[],'MotorSampPts',[], 'AllTaskSampPts', [], 'RestExists', [], 'RestSampPts', [], 'OddEven', [], 'OddEvenNumber', []);
    
        for j=1:length(tasks)
       
        
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
           
            
                % Create template path for resting data
                %MSCciftidir = ['/projects/b1081/MSC/TaskFC/FCProc_' subs{i} '_mem_pass2/cifti_timeseries_normalwall_native_freesurf'];
        
        
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
            
            SubStruct(k).SessionFile = vcid;
            
                disp(sprintf('Loading tmasks for timeseries for session %s for subject %s task %s: %s', vcid, subs{i}, tasks{j}, datestr(now)));

                MSCciftidir = strcat(MSCcondidir, 'cifti_timeseries_normalwall_native_freesurf');

            
                %if k==1
                    %cd('cifti_timeseries_normalwall_native_freesurf')
                %end
            
                tmask = 0;  % Create empty tmask variable
            
                try
                    %data = ft_read_cifti_mod(['/projects/b1081/MSC/MSCdata_v1/' subs{i} '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_normalwall_native_freesurf/' vcid]);
                    %disp(sprintf('Loading data size %i by %i, %s', size(data.data,1), size(data.data,2), datestr(now)));
                    %data = data.data;
                
                
                    if strcmp(tasks{j},'motor')
                    	tmask = TIndFin(k).AllMotor;
                    elseif strcmp(tasks{j},'mem')
                    	tmask = TIndFin(k).AllMem;
                    else
                    	tmask = TIndFin(k).AllGlass + TIndFin(k).AllSemantic;
                    end
                        
                    disp(sprintf('tmask for task %s file has %i good sample points, %s', tasks{j}, sum(tmask), datestr(now)));
                
                    resttmask = QC(k).tmask;
                    
                    if sum(resttmask) == 0
                        
                    	SubStruct(k).RestExists = 0;
                    	SubStruct(k).RestSampPts = 0;
                        
                    else
                        
                        SubStruct(k).RestExists = 1;
                    	SubStruct(k).RestSampPts = sum(resttmask);
                        
                    end
                    
                    disp(sprintf('tmask for rest file has %i good sample points, %s', sum(resttmask), datestr(now)));
                
                    if sum(tmask) == 0      % Skip file if no task task data
                    
                        if strcmp(tasks{j},'motor')
                            SubStruct(k).MotorExists = 0;
                            SubStruct(k).MotorSampPts = 0;
                        elseif strcmp(tasks{j},'mem')
                            SubStruct(k).MemExists = 0;
                            SubStruct(k).MemSampPts = 0;
                        else
                            SubStruct(k).MixedExists = 0;
                            SubStruct(k).MixedSampPts = 0;
                        end
                        
                        disp(sprintf('Vcid %s for subject %s has no usable data for task %s', vcid, subs{i}, tasks{j}));
                    
                    else                    % Else get first xth points in resting file corresponding to number of task data points
                            
                        if strcmp(tasks{j},'motor')
                            SubStruct(k).MotorExists = 1;
                            SubStruct(k).MotorSampPts = sum(tmask);
                        elseif strcmp(tasks{j},'mem')
                            SubStruct(k).MemExists = 1;
                            SubStruct(k).MemSampPts = sum(tmask);
                        else
                            SubStruct(k).MixedExists = 1;
                            SubStruct(k).MixedSampPts = sum(tmask);
                        end
                            
                            
                            if mod(k,2) == 1
                                
                                SubStruct(k).OddEven = 1;
                        
                                OddCount = OddCount + 1;
                                
                                SubStruct(k).OddEvenNumber = (OddCount - 10);
                            
                            else
                                
                                SubStruct(k).OddEven = 2;
                                
                                EvenCount = EvenCount + 1;
                                
                                SubStruct(k).OddEvenNumber = (EvenCount - 10);
                            
                            end
                    end
                
                catch ME
                    if strcmp(ME.message,'Invalid file identifier.  Use fopen to generate a valid file identifier.')
                        continue
                    end
                end
                
                
                
            end         % End of session loop to create concatenated session data
        %end        % End of task loop to create concatenated task data
        %within each task

            
           % end % End of session loop for saving each session separately
        end % End of task data loop for saving each task separately
        
        for l = 1:size(SubStruct,2)  %% Calculate data points for all tasks
            
            SubStruct(l).AllTaskSampPts = (SubStruct(l).MotorSampPts + SubStruct(l).MemSampPts + SubStruct(l).MixedSampPts);
            
        end
        
        % Save final struct for each subject
        
        save(['/' outdir '/' subs{i} '_QCFile.mat'], 'SubStruct');
        
    end
        
elseif WashU == 1
    
    tmask1 = load('/projects/b1081/Poldrome/tmasks_WashU/vc39556_total_tmask.txt');
    tmask2 = load('/projects/b1081/Poldrome/tmasks_WashU/vc39556_2_total_tmask.txt');
        
    WashUtmaskcat = [tmask2; tmask1];
    
    [~,~,condition1,~] = textread(['/projects/b1081/Poldrome/conditions_WashU/vc39556.studies.txt'],'%s%s%s%s');
    [~,~,condition2,~] = textread(['/projects/b1081/Poldrome/conditions_WashU/vc39556_2.studies.txt'],'%s%s%s%s');
    
    conditionlist1 = (~cellfun(@isempty,strfind(condition1,'RSFC')));
    conditionlist1(7) = 0;
    conditionlist2 = (~cellfun(@isempty,strfind(condition2,'RSFC')));
    
    conditionlistfinal = [condition2(conditionlist2); condition1(conditionlist1)];
    
    SubStruct = struct('SessionFile',{}, 'RestExists', [], 'RestSampPts', [], 'Condition', {}, 'Session', [], 'ConditionHalf', []);
    
    eyesopencount = 0;
    eyesclosedcount = 0;
    
    for l = 1:length(conditionlistfinal)
        
        SubStruct(l).SessionFile = conditionlistfinal(l);
        
        if sum(WashUtmaskcat(((l-1)*240)+1:((l-1)*240)+240)) == 0
            
            SubStruct(l).RestExists = 0;
            
        else
            
            SubStruct(l).RestExists = 1;
            
        end
        
        SubStruct(l).RestSampPts = sum(WashUtmaskcat(((l-1)*240)+1:((l-1)*240)+240));
        
        if any(strfind(conditionlistfinal{l},'EO'))
        
            SubStruct(l).Condition = 'Eyes_Open';
            eyesopencount = eyesopencount + 1;
            
            if mod(eyesopencount,2) == 1
                
                SubStruct(l).ConditionHalf = 1;
                
            else
                
                SubStruct(l).ConditionHalf = 2;
                
            end
            
        else
            
            SubStruct(l).Condition = 'Eyes_Closed';
            eyesclosedcount = eyesclosedcount + 1;
            
            if mod(eyesclosedcount,2) == 1
                
                SubStruct(l).ConditionHalf = 1;
                
            else
                
                SubStruct(l).ConditionHalf = 2;
                
            end
            
        end
        
        if l < 11
            
            SubStruct(l).Session = 2;
            
        else
            
            SubStruct(l).Session = 1;
            
        end
    end
    
    % Save final struct
        
    save(['/' outdir '/WashU_QCFile.mat'], 'SubStruct');

elseif Texas == 1
    
	restdir = dir('/projects/b1081/Poldrome/cifti_data_Texas/cifti_timeseries_normalwall');
	vcidlist = restdir(~cellfun(@isempty,strfind({restdir.name},'dtseries')));
    
    SubStruct = struct('SessionFile',{}, 'RestExists', [], 'RestSampPts', [], 'EquivSampPts', [], 'SessionHalf', []);

    sessioncount = 0;
    
    for l = 1:size(vcidlist,1)
        
        sessioncount = sessioncount + 1;
        
        vcid = vcidlist(l).name;
    
        vcidmatch = strtok(vcid,'_z');
        vcidmatch = strcat(vcidmatch,'_total_tmask.txt');
                    
        tmaskimport = load(['/projects/b1081/Poldrome/tmasks_Texas/' vcidmatch]);
        
        SubStruct(l).SessionFile = vcidmatch;
        
        if sum(tmaskimport) == 0
            
            SubStruct(l).RestExists = 0;
            
        else
            
            SubStruct(l).RestExists = 1;
            
        end
        
        SubStruct(l).RestSampPts = sum(tmaskimport);
        SubStruct(l).EquivSampPts = round(sum(tmaskimport)/(2.5/1.16));
        
        if mod(sessioncount,2) == 1
                
        	SubStruct(l).SessionHalf = 1;
                
        else
                
        	SubStruct(l).SessionHalf = 2;
                
        end
 
    end
        
    % Save final struct
        
    save(['/' outdir '/Texas_QCFile.mat'], 'SubStruct');
        
end
        
        


disp(sprintf('Job Completed: %s', datestr(now)));