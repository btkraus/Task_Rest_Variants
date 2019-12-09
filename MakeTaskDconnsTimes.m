parpool('local', 24)     %% Name of cluster profile for batch job

%/davta/vnil-bluearc/GMT/Caterina/TaskFC/FCProc_MSC01_motor_passvv2/
%/data/nil-bluearc/GMT/Caterina/TaskFC/FCProc_MSC01_mem_pass2/
%/data/nil-bluearc/GMT/Caterina/TaskFC/FCProc_MSC01_mixed_pass2/
clear all

outdir = '/projects/b1081/Brian_MSC/dconn_task_files';
dataLocStem = '/MSC/TaskFC/';
%subs = {'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC10'};
subs = {'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC09','MSC10'};
tasks = {'motor','mem','mixed'};
%tasks = {'mem'};

CortexOnly = 1; %% Toggles whether to run correlations on cortex only
SplitHalf = 0;  %% Toggles whether to create a separate file for odd/even sessions
MatchData = 0; %% Toggles whether to match the amount of data per task as the lowest value within each split-half
RandSample = 0; %% Toggles whether to randomly sample data from each session
ConsecSample = 1;  %% Toggles whether to consecutively sample data
MatchAcrossTasks = 0;  %% Toggles whether to match the amount of data across tasks
ConcatenateTasks = 1;   %% Toggles whether to concatenate data for all tasks
MakeDconn = 0;  %% Toggles whether to write a dconn
MakeVariantMap = 1; %% Toggles whether to write a variant map
ConcatenateSplitHalf = 0;  %% Toggles whether to concatenate split halves into one vector
SaveTimeseries = 0;     %% Save concatenated timeseries for subject

times = [5:5:100];      %% Times to make spatial correlation maps


cd '/projects/b1081';   %% Change CD to root project directory

disp(sprintf('Job Submitted: %s', datestr(now)));




disp(sprintf('Job Started: %s', datestr(now)));

if CortexOnly == 1      %% Select correct number of voxels for template
    
    voxnum = 59412;
    
else
    
    voxnum = 65625;
    
end


for i=1:numel(subs)
    
    for l = 1:numel(times)
    
        disp(sprintf('Creating dconn for subject %s: %s', subs{i}, datestr(now)));


        % Initialize cat data
        if SplitHalf == 1

            catData1 = [];
            catData2 = [];

        else

            catData = [];

        end
        
        if ConsecSample == 1
            
            load(['/projects/b1081/Brian_MSC/QC_files/' subs{i} '_QCFile.mat']);
            
            totalsamppts = round((times(l)*(60/2.2)));
            tasksamppts = ceil(totalsamppts/3);
            
           	ptscountmem = 0;       %% Running sample point count across sessions
            ptscountmixed = 0;
            ptscountmotor = 0;
           	enoughptsmem = 0;      %% Toggles whether enough sample points have been taken consecutively from each file
            enoughptsmixed = 0;
            enoughptsmotor = 0;
            maxptsmem = 0;         %% Toggles whether the maximum amount of data for each task has been sampled
            maxptsmixed = 0;
            maxptsmotor = 0;
           	sampspersessionmem = zeros(length(SubStruct),1);
            sampspersessionmixed = zeros(length(SubStruct),1);
            sampspersessionmotor = zeros(length(SubStruct),1);
            
            remaindertaskpts = 0;
            
            while true

                for u = 1:length(SubStruct)

                    if SubStruct(u).OddEven == 1 && (enoughptsmem == 0 && maxptsmem == 0)
                            
                        if SubStruct(u).MemSampPts + ptscountmem >= tasksamppts + remaindertaskpts
                            
                            sampspersessionmem(u) = sampspersessionmem(u) + (tasksamppts - ptscountmem);
                            ptscountmem = sampspersessionmem(u) + ptscountmem;
                            enoughptsmem = 1;
                            
                        elseif SubStruct(u).MemSampPts + ptscountmem < tasksamppts + remaindertaskpts
                            
                            sampspersessionmem(u) = sampspersessionmem(u) + SubStruct(u).MemSampPts;
                            ptscountmem = SubStruct(u).MemSampPts + ptscountmem;
                                
                        end
                        
                        if u == length(SubStruct) && enoughptsmem == 0
                            
                            maxptsmem = 1;
                            
                        end
                        
                    end
                    
                    if SubStruct(u).OddEven == 1 && enoughptsmixed == 0
                    
                        if SubStruct(u).MixedSampPts + ptscountmixed >= tasksamppts + remaindertaskpts
                            
                            sampspersessionmixed(u) = sampspersessionmixed(u) + (tasksamppts - ptscountmixed);
                            ptscountmixed = sampspersessionmixed(u) + ptscountmixed;
                            enoughptsmixed = 1;
                            
                        elseif SubStruct(u).MixedSampPts + ptscountmixed < tasksamppts + remaindertaskpts
                            
                            sampspersessionmixed(u) = sampspersessionmixed(u) + SubStruct(u).MixedSampPts;
                            ptscountmixed = SubStruct(u).MixedSampPts + ptscountmixed;
                                
                        end
                        
                        if u == length(SubStruct) && enoughptsmixed == 0
                            
                            maxptsmixed = 1;
                            
                        end
                        
                    end
                    
                    if SubStruct(u).OddEven == 1 && enoughptsmotor == 0
                    
                        if SubStruct(u).MotorSampPts + ptscountmotor >= tasksamppts + remaindertaskpts
                            
                            sampspersessionmotor(u) = sampspersessionmotor(u) + (tasksamppts - ptscountmotor);
                            ptscountmotor = sampspersessionmotor(u) + ptscountmotor;
                            enoughptsmotor = 1;
                            
                        elseif SubStruct(u).MotorSampPts + ptscountmotor < tasksamppts + remaindertaskpts
                            
                            sampspersessionmotor(u) = sampspersessionmotor(u) + SubStruct(u).MotorSampPts;
                            ptscountmotor = SubStruct(u).MotorSampPts + ptscountmotor;
                                
                        end
                        
                        if u == length(SubStruct) && enoughptsmotor == 0
                            
                            maxptsmotor = 1;
                            
                        end
                        
                    end
                    
                end
                
                if (ptscountmem + ptscountmixed + ptscountmotor) >= totalsamppts || (maxptsmem == 1 && maxptsmixed == 1 && maxptsmotor == 1)
                    
                    break
                    
                else
                    
                    remainderpts = totalsamppts - (ptscountmem + ptscountmixed + ptscountmotor);
                    
                    if maxptsmem == 1 && maxptsmixed == 1
                        
                        remaindertaskpts = totalsamppts - (ptscountmem + ptscountmixed);
                        
                        enoughptsmotor = 0;
                        sampspersessionmotor = zeros(length(SubStruct),1);
                        
                    elseif maxptsmem == 1 && maxptsmotor == 1
                        
                        remaindertaskpts = totalsamppts - (ptscountmem + ptscountmotor);
                        
                        enoughptsmixed = 0;
                        sampspersessionmixed = zeros(length(SubStruct),1);
                        
                    elseif maxptsmixed == 1 && maxptsmotor == 1
                        
                        remaindertaskpts = totalsamppts - (ptscountmixed + ptscountmotor);
                        
                        enoughptsmem = 0;
                        sampspersessionmem = zeros(length(SubStruct),1);
                        
                    elseif maxptsmem == 1
                        
                        remaindertaskpts = ceil((totalsamppts - ptscountmem)/2);
                        
                        enoughptsmixed = 0;
                        sampspersessionmixed = zeros(length(SubStruct),1);
                        enoughptsmotor = 0;
                        sampspersessionmotor = zeros(length(SubStruct),1);
                        
                    elseif maxptsmixed == 1
                        
                        remaindertaskpts = ceil((totalsamppts - ptscountmixed)/2);
                        
                        enoughptsmem = 0;
                        sampspersessionmem = zeros(length(SubStruct),1);
                        enoughptsmotor = 0;
                        sampspersessionmotor = zeros(length(SubStruct),1);
                        
                    elseif maxptsmotor == 1
                        
                        remaindertaskpts = ceil((totalsamppts - ptscountmotor)/2);
                        
                        enoughptsmem = 0;
                        sampspersessionmem = zeros(length(SubStruct),1);
                        enoughptsmixed = 0;
                        sampspersessionmixed = zeros(length(SubStruct),1);
                    
                    end
                end
                    
            end
        end
        
        if maxptsmem == 1 && maxptsmixed == 1 && maxptsmotor == 1
        
            disp(sprintf('Subject %s does not have %i minutes of total task data. There are %i data points for the memory task, %i data points for the mixed task, and %i data points for the motor task: %s', subs{i}, times(l), ptscountmem, ptscountmixed, ptscountmotor, datestr(now)));
            
        else
            
            disp(sprintf('Subject %s has enough data for %i minutes of total task data. There are %i data points for the memory task, %i data points for the mixed task, and %i data points for the motor task: %s', subs{i}, times(l), ptscountmem, ptscountmixed, ptscountmotor, datestr(now)));
            
        end
        
        for j=1:length(tasks)
        

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
                        
                        sampspersession = sampspersessionmotor(k);
                        
                    elseif strcmp(tasks{j},'mem')
                        
                        tmask = TIndFin(k).AllMem;
                        
                        sampspersession = sampspersessionmem(k);
                        
                    else
                        
                        tmask = TIndFin(k).AllGlass + TIndFin(k).AllSemantic;
                        
                        sampspersession = sampspersessionmixed(k);
                        
                    end

                    if sum(tmask) == 0      % Skip file if no task task data

                        data = [];
                        disp(sprintf('Vcid %s for subject %s has no usable data for task %s', vcids{k}, subs{i}, tasks{j}));
                        
                    elseif ConsecSample == 1 && numel(times) > 1
                        
                     	tmasknum = sampspersession;
                            
                      	if tmasknum > 0
                                
                          	data = data(1:voxnum,logical(tmask));
                                
                          	if tmasknum < sum(tmask)
                                    
                              	data = data(:,1:tmasknum);
                                    
                            end
                        end
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

            %disp(sprintf('Reading %s: %s', [MSCciftidir '/' vcid], datestr(now)));
            %template = ft_read_cifti_mod([MSCciftidir '/' vcid]);
            disp(sprintf('Reading %s: %s', '/projects/b1081/MSC/TaskFC/FCProc_MSC01_mem_pass2/cifti_timeseries_normalwall_native_freesurf/vc38671_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii', datestr(now)));
            template = ft_read_cifti_mod('/projects/b1081/MSC/TaskFC/FCProc_MSC01_mem_pass2/cifti_timeseries_normalwall_native_freesurf/vc38671_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii');

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

            disp(sprintf('Running Correlations: on data size %i by %i, %s', size(catData,1), size(catData,2), datestr(now)));
            template.data = paircorr_mod(catData');
            catData = [];

        end

        if ConcatenateTasks == 1 && CortexOnly == 1 && ConsecSample == 1

            if MakeDconn == 1 
                
                disp(sprintf('Writing %s: %s', [outdir '/' subs{i} '_' num2str(times(l)) '_minutes_total_consec_Odd_allTaskCat_cortex'], datestr(now)));
                ft_write_cifti_mod([outdir '/' subs{i} '_' num2str(times(l)) '_minutes_total_consec_Odd_allTaskCat_cortex'],template)

            end

            if MakeVariantMap == 1

                outputfile = [subs{i} '_' num2str(times(l)) '_minutes_total_consec_Odd_allTaskCat_cortex'];

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
end

disp(sprintf('Job Completed: %s', datestr(now)));
