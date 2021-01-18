%This script makes Dconns (if the MakeDconn variable = 1) and variant files
% (if VariantMap = 1)for rest data of each subject's even and odd sessions. 
% QC files need to be na med [subject]_QCFile.mat, ex. 'MSC01_QCFile.mat'

%parpool('local', 28)     %% Name of cluster profile for batch job (how many workers/cores you need to run this job)

clear all

disp(sprintf('Job Submitted: %s', datestr(now)));

%% Paths

groupavg_path = 'projects/b1081/Atlases'; %specify location of group average map
groupavg_name = '120_allsubs_corr'; %specify name of group average map
rest_path = '/projects/b1081/MSC/MSCdata_v1/'; %specify location of rest data

outdir = '/projects/b1081/member_directories/zladwig/Analysis_Scripts_Replication'; %specify output directory
dataLocStem = '/MSC/TaskFC'; %specify location of data
QC_path = '/projects/b1081/member_directories/bkraus/Brian_MSC/Analysis_Scripts_Replication/QC_files/';
template_path = '/projects/b1081/member_directories/bkraus/Brian_MSC/dconn_scripts/templates/MSC01_allses_mean_native_freesurf_vs_120sub_corr.dtseries.nii';
atlas_path = '/projects/b1081/Atlases';
output_path = '/projects/b1081/member_directories/zladwig/Analysis_Scripts_Replication/variant_maps';
cd '/projects/b1081';   %% Change CD to root project directory
%addpath(genpath('/projects/b1081/member_directories/bkraus/Brian_MSC/'));
addpath(genpath('/projects/b1081/member_directories/dperez/Analysis_Scripts_Replication/'));
addpath(genpath('/projects/b1081/member_directories/zladwig/Task_Rest_Variants/'));
addpath(genpath('/projects/b1081/member_directories/zladwig/dependencies/cifti-matlab-master/'));


%% Options - deleted RestOnly (default is 0), SubSample (default is 1), ConcatenateSplitHalf (default is 0)
SplitHalf = 1;  %% Toggles whether to calculate rest for split-half of sessions
RandSample = 1; %% Toggles random subsampling of tmask on/off
MatchData = 1; %% Toggles whether to match the amount of data per task as the lowest value within each split-half
MatchAcrossTasks = 0; %% Toggles whether to match the amount of data across tasks (for concatenated task data)
MatchAcrossSubs = 1;  %% Toggles whether to match the amount of data across subjects
ConcatenateTasksMatch = 0; %% Toggles whether to calculate matched task data for concatenated tasks
ConcatenateTaskData = 0; %% Toggles whether to concatenate task data within subjects
ConcatenateSessionData = 0; %% Toggles whether to concatenate session data within subjects
WriteDconn = 0; %% Toggles whether to write dconn to disk
CreateVariant = 0; %% Toggles whether to create a variant map from dconn
SaveTimeseries = 0;     %% Save concatenated timeseries for subject
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
        
    [sampspersession, restsampspersession, tasksampspersession] = GetSampsPerSessionRest(mintasksamps, MatchAcrossTasks, MatchAcrossSubs, ConcatenateTaskData, ConcatenateTasksMatch, QC_path, SplitHalf, MatchData, tasks, subs{i});
%%
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