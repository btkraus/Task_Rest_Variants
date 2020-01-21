% This script makes Dconns and variant maps for each subject 
% QC files need to be named [subject]_QCFile.mat, ex. 'MSC01_QCFile.mat'
% Written by Brian Kraus, edited by Diana Perez, January 2020.

parpool('local', 28)     %% Name of cluster profile for batch job (how many workers/cores you need to run this job)

clear all

disp(sprintf('Job Submitted: %s', datestr(now)));
disp(sprintf('Job Started: %s', datestr(now)));

%% Paths
outdir = '/projects/b1081/Brian_MSC/Analysis_Scripts_Replication/dconn_task_files'; %specify output directory
dataLocStem = '/MSC/TaskFC/'; %specify location of data
QCFiles_path = '/projects/b1081/Brian_MSC/Analysis_Scripts_Replication/QC_files/';
cd '/projects/b1081';   %% Change CD to root project directory
%% Options
MakeDconn = 0;  %% Toggles whether to write a dconn
MakeVariantMap = 1; %% Toggles whether to write a variant map
SaveTimeseries = 0;     %% Save concatenated timeseries for subject
%% Variables
subs = {'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC09','MSC10'};
tasks = {'motor','mem','mixed'};
voxnum = 59412; % number of voxels for template

% sets up variables for number of sample points in even- and odd-numbered
% sessions to be used later
memptsoddsum = [];
motorptsoddsum = [];
mixedptsoddsum = [];
memptsevensum = [];
motorptsevensum = [];
mixedptsevensum = [];
%% Main for-loop: Creates Dconn for each subject
for i=1:numel(subs)
    
    disp(sprintf('Creating dconn for subject %s: %s', subs{i}, datestr(now)));
    
    catData1 = [];
    catData2 = [];
    
    for j=1:length(tasks)
   
        load ([QCFiles_path subs{i} '_QCFile.mat']);

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
                motorptsodd = [motorptsodd; SubStruct(u).MotorSampPts];
            elseif SubStruct(u).OddEven == 2
                memptseven = [memptseven; SubStruct(u).MemSampPts];
                mixedptseven = [mixedptseven; SubStruct(u).MixedSampPts];
                motorptseven = [motorptseven; SubStruct(u).MotorSampPts];
            end
        end

        % determines minimum # of sample points
        if strcmp(subs{i}, 'MSC09')  %% Removes motor task from consideration for MSC09
            minsampspts = min([sum(memptsodd) sum(mixedptsodd) sum(memptseven) sum(mixedptseven)]);
        else
            minsampspts = min([sum(memptsodd) sum(mixedptsodd) sum(motorptsodd) sum(memptseven) sum(mixedptseven) sum(motorptseven)]);
        end

        meansamppts = floor(mintasksamps/5);

        disp(sprintf('For all subjects, the minimum number of sample points in a split-half is %i, with a mean of %i points per sesssion: %s', mintasksamps, meansamppts, datestr(now)));

        meansampptstempodd = meansamppts;
        meansampptstempeven = meansamppts;

       %% Sets motor task data points for subject MSC09
       if strcmp(subs{i}, 'MSC09') && strcmp(tasks{j}, 'motor') 
            meansampptstempodd = floor(sum(motorptsodd)/5);
            meansampptstempeven = floor(sum(motorptseven)/5);
       elseif strcmp(subs{i}, 'MSC09')
            meansampptstempodd = meansamppts + round((meansamppts - round(sum(motorptsodd)/5))/2);
            meansampptstempeven = meansamppts + round((meansamppts - round(sum(motorptseven)/5))/2);      
       end
            
       %% NOT SURE WHAT THESE VARIABLES AND WHILE LOOP ARE FOR
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
            
        %% what is the remainder and what are those values for sum of notenoughdata?
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

        disp(sprintf('Creating dconn for subject %s task %s: %s', subs{i}, tasks{j}, datestr(now)));
        
        %% ARE THE DATALIST files included in the MSC data download or do we have to include something in the read-me about making these txt files?
        % Load vcids (session identifiers for current subject)
        MSCTaskdir = strcat(cd, dataLocStem);
        [~,vcids,~,~,~] = textread([MSCTaskdir '/' subs{i} '_' tasks{j} '_DATALIST.txt'],'%s%s%s%s%s');
        
        disp(sprintf('%i sessions found for subject %s task %s: %s', size(vcids,1), subs{i}, tasks{j}, datestr(now)));
        
        %% Why is it different for MSC03 and MSC10? 
        % Load tmasks
        if strcmp(subs{i}, 'MSC03') || strcmp(subs{i}, 'MSC10')
        	MSCcondidir = strcat(strcat(MSCTaskdir, '/'), ['FCProc_' subs{i} '_' tasks{j} '_pass2_FDfilt/']);
        else
            MSCcondidir = strcat(strcat(MSCTaskdir, '/'), ['FCProc_' subs{i} '_' tasks{j} '_pass2/']);   
        end
         
        %what are condindices.mat files?
        load ([MSCcondidir 'condindices.mat']);

        % Load and concatentate data
        for k=1:length(vcids)            
            disp(sprintf('Loading timeseries for session %s for subject %s task %s: %s', vcids{k}, subs{i}, tasks{j}, datestr(now)));
            MSCciftidir = strcat(MSCcondidir, 'cifti_timeseries_normalwall_native_freesurf');
            try
                % reads data files, ft_read_cifti_mod is a function in the
                % "private" folder within cifti-matlab-master
                data = ft_read_cifti_mod([MSCciftidir '/' vcids{k} '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii']);
                data = data.data;
                %% where does TIndFin come from? the cifti files?
                if strcmp(tasks{j},'motor')
                	tmask = TIndFin(k).AllMotor;
                elseif strcmp(tasks{j},'mem')
                	tmask = TIndFin(k).AllMem;
                else
                    tmask = TIndFin(k).AllGlass + TIndFin(k).AllSemantic;
                end
                
                % Skip file if no task data, randomly selects sample
                % points if task data exists
                if sum(tmask) == 0      
                    data = [];
                    disp(sprintf('Vcid %s for subject %s has no usable data for task %s', vcids{k}, subs{i}, tasks{j}));
                else
                	tmasknum = sampspersession(k);    
                	samplepts = find(tmask == 1);
                	tmask(samplepts) = 0;
                	sampleselect = datasample(samplepts,tmasknum,'Replace',false);
                    tmask(sampleselect) = 1;                            
                 	disp(sprintf('%i sample points from file %s have been randomly selected for task %s, %s', sum(tmask), vcids{k}, tasks{j}, datestr(now)));
                	data = data(1:voxnum,logical(tmask));  
                end

                disp(sprintf('Data is %i by %i and mod operation equals %i', size(data,1), size(data,2), mod(k,2)));

                if mod(k,2) == 1      
                    disp(sprintf('Adding %i data points from session %s for subject %s to odd data', size(data,2), vcids{k}, subs{i}));
                    catData1 = [catData1 data];                            
                else                                
                    disp(sprintf('Adding %i data points from session %s for subject %s to even data', size(data,2), vcids{k}, subs{i}));                                
                    catData2 = [catData2 data];
                end

            catch ME
                if strcmp(ME.message,'Invalid file identifier.  Use fopen to generate a valid file identifier.')
                    continue
                end
            end
        end 
    end        %% Ending task loop here creates concatenated task data
    
    
    
    if SaveTimeseries == 1          %% Save concatenated FC processed timeseries       
        %% WILL THE TEMPLATE ALWAYS BE MSC01? I NEED TO MOVE THIS PATH TO THE TOP BUT FIRST FIGURE OUT HOW TO SET THIS UP       
        timseriestemplate = ft_read_cifti_mod('/projects/b1081/Brian_MSC/Analysis_Scripts_Replication/dconn_scripts/templates/MSC01_allses_mean_native_freesurf_vs_120sub_corr.dtseries.nii');
        timseriestemplate.data = [];        
        timseriestemplate.data = catData;
        disp(sprintf('Data timseries is size %i by %i, %s', size(catData,1), size(catData,2), datestr(now)));
     	disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_FCProcessed_Timeseries_Matched_allTaskCat_cortex'], datestr(now)));
    	ft_write_cifti_mod(['/' outdir '/' subs{i} '_FCProcessed_Timeseries_Matched_allTaskCat_cortex'],timseriestemplate)
    end

    % Make and save rmat    
    disp('Loading template: MSC01_allses_mean_native_freesurf_vs_120sub_corr.dtseries.nii');
   	template = ft_read_cifti_mod('/projects/b1081/Brian_MSC/Analysis_Scripts_Replication/dconn_scripts/templates/MSC01_allses_mean_native_freesurf_vs_120sub_corr.dtseries.nii');    
    template.data = [];
    template.dimord = 'pos_pos';
    template.hdr.dim(6) = template.hdr.dim(7);
    template.hdr.intent_code = 3001;
    template.hdr.intent_name = 'ConnDense';
    template2 = template;

    if MakeDconn == 1        
    	disp(sprintf('Running Correlations: on data size %i by %i, %s', size(catData1,1), size(catData1,2), datestr(now)));
    	template.data = paircorr_mod(catData1');
        catData1 = [];
    	disp(sprintf('Running Correlations: on data size %i by %i, %s', size(catData2,1), size(catData2,2), datestr(now)));
     	template2.data = paircorr_mod(catData2');
        catData2 = [];
        disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_matcheddata_Final_allTaskCat_OddSessions_cortex'], datestr(now)));
        ft_write_cifti_mod(['/' outdir '/' subs{i} '_matcheddata_Final_allTaskCat_OddSessions_cortex'],template)        
        disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_matcheddata_Final_allTaskCat_EvenSessions_cortex'], datestr(now)));
        ft_write_cifti_mod(['/' outdir '/' subs{i} '_matcheddata_Final_allTaskCat_EvenSessions_cortex'],template2)
    end
    
    if MakeVariantMap == 1           
        outputfile1 = [subs{i} '_matcheddata_Final_allTaskCat_EvenSessions_cortex'];
        outputfile2 = [subs{i} '_matcheddata_Final_allTaskCat_OddSessions_cortex'];
        
        if MakeDconn == 0             
            disp(sprintf('Running Correlations: on data size %i by %i, %s', size(catData1,1), size(catData1,2), datestr(now)));
            template.data = paircorr_mod(catData1');            
            catData1 = [];           
            %% NEED TO FIND THIS FUNCTION (CAN'T FIND IN PROJECT DIRECTORY) AND PUT PATHS AT THE TOP
            createSptlcorr_MSCdconns('/projects/b1081/Atlases', '120_allsubs_corr',1,'/projects/b1081/Brian_MSC/Analysis_Scripts_Replication/output_files/variant_maps',template.data, outputfile1)            
            template.data = [];            
            disp(sprintf('Running Correlations: on data size %i by %i, %s', size(catData2,1), size(catData2,2), datestr(now)));
            template2.data = paircorr_mod(catData2');            
            catData2 = [];            
            createSptlcorr_MSCdconns('/projects/b1081/Atlases', '120_allsubs_corr',1,'/projects/b1081/Brian_MSC/Analysis_Scripts_Replication/output_files/variant_maps',template2.data, outputfile2)            
            template2.data = [];        
            %% What is the difference between this and the one above????
        else            
            createSptlcorr_MSCdconns('/projects/b1081/Atlases', '120_allsubs_corr',1,'/projects/b1081/Brian_MSC/Analysis_Scripts_Replication/output_files/variant_maps',template.data, outputfile1)            
            createSptlcorr_MSCdconns('/projects/b1081/Atlases', '120_allsubs_corr',1,'/projects/b1081/Brian_MSC/Analysis_Scripts_Replication/output_files/variant_maps',template2.data, outputfile2)            
        end
    end

        
        
            
    clear template
    clear template2
    
    %end        %% Ending task loop here creates file for each task
end

disp(sprintf('Job Completed: %s', datestr(now)));
