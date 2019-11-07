
clear all

pearsoncoefs = [];
dicecoefs = [];
taskoutput = {};
threshoutput = [];
outputdir = '/Users/briankraus/Desktop/Correct_Variant_Maps/TaskVariantPlots';

tasktaskpearsoncoefs = [];
taskdifftaskpearsoncoefs = [];
taskrestpearsoncoefs = [];
restrestpearsoncoefs = [];
taskresttaskrestpearsoncoefs = [];
alltaskrestpearsoncoefs = [];
alltaskalltaskpearsoncoefs = [];
alltaskrestalltaskrestpearsoncoefs = [];

tasktaskdicecoefs = [];
taskdifftaskdicecoefs = [];
taskrestdicecoefs = [];
restrestdicecoefs = [];
taskresttaskrestdicecoefs = [];
alltaskrestdicecoefs = [];
alltaskalltaskdicecoefs = [];
alltaskrestalltaskrestdicecoefs = [];

randomizevals = 0;  %% Toggles whether to calculate a randomized value for each task
plotcorrs = 0;      %% Toggles whether to plot and save the correlations for each subject
plotavgcorrs = 0;   %% Toggles whether to plot averaged correlations for pairwise task comparisons
concatplots = 1;    %% Toggles whether to concatenate all the plots and save them as one file
saveindividplots = 0;  %% Toggles whether to save all individual plots
plottaskdiffs = 1;  %% Toggles whether to plot individual tasks versus all tasks
thresholds = 2.5;
%thresholds = [5 10 15 20 25 30 35];  %% Thresholds to run dice correlations (percent lowest correlations)
%thresholds = [.15 .2 .25 .3 .35 .4 .45]; %% Thresholds to run dice correlations (absolute correlation threshold)
absolute = 0;   %% Toggles whether absolute or percent thresholds should be used
fishertransform = 0;  %% Toggles whether to Fisher transform the values before plotting
SNRMasks = 1;  %% Toggles whether to exclude areas of low SNR via SNR masks
SizeExclude = 1;  %% Toggles whether to exclude variants smaller than criterion
minsize = 15;  %% Minimum variant size for size exclusion
%comparison = {'task-task', 'task-difftask', 'rest-rest', 'task-rest', 'taskrest-taskrest', 'alltask-alltask', 'alltask-rest', 'alltaskrest-alltaskrest'};
comparison = {'matched-samedata-task-task', 'matched-samedata-task-difftask', 'matched-samedata-task-rest', 'matched-samedata-taskrest-taskrest', 'matched-samedata-rest-rest', 'matched-samedata-alltask-alltask', 'matched-samedata-alltask-rest'};
%comparison = {'matched-samedata-rest-rest'};



if randomizevals == 1
        
	memfilestask = [];
	mixedfilestask = [];
	motorfilestask = [];
	memfilesrest = [];
	mixedfilesrest = [];
	motorfilesrest = [];
	alltaskfilestask = [];
	alltaskfilesrest = [];
        
end


for g = 1:numel(comparison)
    
    comp = comparison(g);
    
    if strcmp(comp, 'task-task')

        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_task_even_varmaps.txt','%s%s%s');

        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_task_odd_varmaps.txt','%s%s%s');
        
    elseif strcmp(comp, 'task-difftask')
        
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_task_difftask_even_varmaps.txt','%s%s%s');

        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_difftask_difftask_odd_varmaps.txt','%s%s%s');
        
    elseif strcmp(comp, 'rest-rest')
        
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_even_varmaps.txt','%s%s%s');

        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_odd_varmaps.txt','%s%s%s');
        
    elseif strcmp(comp, 'task-rest')
        
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_task_even_varmaps.txt','%s%s%s');

        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_task_even_varmaps.txt','%s%s%s');
        
    elseif strcmp(comp, 'taskrest-taskrest')
        
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_task_odd_varmaps.txt','%s%s%s');
        
        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_task_even_varmaps.txt','%s%s%s');
        
    elseif strcmp(comp, 'alltask-alltask')
        
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_even_varmaps.txt','%s%s%s');

        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_odd_varmaps.txt','%s%s%s');
        
    elseif  strcmp(comp, 'alltask-rest')
        
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_even_varmaps.txt','%s%s%s');
        
        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_alltask_even_varmaps.txt','%s%s%s');
        
    elseif strcmp(comp, 'alltaskrest-alltaskrest')
        
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_alltask_odd_varmaps.txt','%s%s%s');
        
        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_alltask_even_varmaps.txt','%s%s%s');
        
    elseif strcmp(comp, 'matched-samedata-task-task')
        
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_task_samedata_even_varmaps.txt','%s%s%s');
        
        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_task_samedata_odd_varmaps.txt','%s%s%s');
        
        [task_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_task_samedata_even_varmaps_SizeOnly.txt','%s%s%s');
        
        [rest_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_task_samedata_odd_varmaps_SizeOnly.txt','%s%s%s');
        
    elseif strcmp(comp, 'matched-samedata-task-difftask')
        
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_task_difftask_samedata_even_varmaps.txt','%s%s%s');

        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_difftask_difftask_samedata_odd_varmaps.txt','%s%s%s');
        
        [task_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_task_difftask_samedata_even_varmaps_SizeOnly.txt','%s%s%s');
        
        [rest_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_difftask_difftask_samedata_odd_varmaps_SizeOnly.txt','%s%s%s');
        
    elseif strcmp(comp, 'matched-samedata-task-rest')
        
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_task_samedata_even_varmaps.txt','%s%s%s');
        
        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_task_samedata_even_varmaps.txt','%s%s%s');
        
        [task_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_task_samedata_even_varmaps_SizeOnly.txt','%s%s%s');
        
        [rest_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_task_samedata_even_varmaps_SizeOnly.txt','%s%s%s');
    
    elseif strcmp(comp, 'matched-samedata-taskrest-taskrest')
        
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_task_samedata_odd_varmaps.txt','%s%s%s');
        
        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_task_samedata_even_varmaps.txt','%s%s%s');
        
        [task_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_task_samedata_odd_varmaps_SizeOnly.txt','%s%s%s');
        
        [rest_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_task_samedata_even_varmaps_SizeOnly.txt','%s%s%s');
        
    elseif strcmp(comp, 'matched-samedata-rest-rest')
        
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_alltask_samedata_even_varmaps.txt','%s%s%s');
        
        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_alltask_samedata_odd_varmaps.txt','%s%s%s');
        
        [task_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_alltask_samedata_even_varmaps_SizeOnly.txt','%s%s%s');
        
        [rest_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_alltask_samedata_odd_varmaps_SizeOnly.txt','%s%s%s');
        
    elseif strcmp(comp, 'matched-samedata-alltask-alltask')
        
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_samedata_even_varmaps.txt','%s%s%s');
        
        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_samedata_odd_varmaps.txt','%s%s%s');
        
        [task_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_samedata_even_varmaps_SizeOnly.txt','%s%s%s');
        
        [rest_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_samedata_odd_varmaps_SizeOnly.txt','%s%s%s');
        
    elseif strcmp(comp, 'matched-samedata-alltask-rest')
        
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_samedata_even_varmaps.txt','%s%s%s');
        
        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_alltask_samedata_even_varmaps.txt','%s%s%s');
        
        [task_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_samedata_even_varmaps_SizeOnly.txt','%s%s%s');
        
        [rest_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_alltask_samedata_even_varmaps_SizeOnly.txt','%s%s%s');
        
    end
        
    subs = unique(subjects1);
    
    
    for s = 1:numel(subs)
        
        tasktaskpearsoncoefstemp = [];
        taskdifftaskpearsoncoefstemp = [];
        taskrestpearsoncoefstemp = [];
        restrestpearsoncoefstemp = [];
        taskresttaskrestpearsoncoefstemp = [];
        alltaskrestpearsoncoefstemp = [];
        alltaskalltaskpearsoncoefstemp = [];
        alltaskrestalltaskrestpearsoncoefstemp = [];

        tasktaskdicecoefstemp = [];
        taskdifftaskdicecoefstemp = [];
        taskrestdicecoefstemp = [];
        restrestdicecoefstemp = [];
        taskresttaskrestdicecoefstemp = [];
        alltaskrestdicecoefstemp = [];
        alltaskalltaskdicecoefstemp = [];
        alltaskrestalltaskrestdicecoefstemp = [];
        


    
        plotlabels = {};
        plotdata = [];
    
        task_filenums = find(strcmp(subjects1, subs{s}))';    %% Gets file/task indices for current subject
        rest_filenums = find(strcmp(subjects2, subs{s}))';    %% Gets file/task indices for current subject
    
        for t = task_filenums
        
            dicecoefstemp = [];
            pearsoncoefstemp = [];
        
            taskoutputtemp = {subs{s}};
        
            if numel(thresholds) == 1
            
                threshold = thresholds;
        
                subject = subs{s};
                task_file = task_files{t};
                rest_file = rest_files{t};
                task = tasks1{t};
                task2 = tasks2{t};
                
                if SizeExclude == 1
                    
                   	task_mask = task_masks{t};
                    rest_mask = rest_masks{t};
                    
                end
        
                taskoutputtemp = [taskoutputtemp; task task2];

                cifti_task = ft_read_cifti_mod(task_file);
                cifti_rest = ft_read_cifti_mod(rest_file);
                
                if SNRMasks == 1
    
                    SNRmask = ft_read_cifti_mod('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/bottomBrainMask.dtseries.nii');
            
                    LowSNR = find(SNRmask.data == 1);
            
                    cifti_task.data(LowSNR,:) = NaN;
                    cifti_rest.data(LowSNR,:) = NaN;
            
                end
        
                r = corrcoef(cifti_rest.data, cifti_task.data);
            
                if absolute == 1
                
                    cifti_task_threshold = find(cifti_task.data < threshold);
                    cifti_rest_threshold = find(cifti_rest.data < threshold);
                
                else
        
                    cifti_task_threshold = find(cifti_task.data < prctile(cifti_task.data,threshold));
                    cifti_rest_threshold = find(cifti_rest.data < prctile(cifti_rest.data,threshold));

                end
            
                cifti_rest_thresh_dat = zeros(size(cifti_rest.data));
                cifti_rest_thresh_dat(cifti_rest_threshold,1) = 1;
            
                cifti_task_thresh_dat = zeros(size(cifti_task.data));
                cifti_task_thresh_dat(cifti_task_threshold,1) = 1;
                
                if SizeExclude == 1
                    
                    cifti_task_mask = ft_read_cifti_mod(task_mask);
                    cifti_rest_mask = ft_read_cifti_mod(rest_mask);
                    
                    allvars_rest = unique(cifti_rest_mask.data);
                    allvars_rest(1) = [];
                    allvars_task = unique(cifti_task_mask.data);
                    allvars_task(1) = [];
                
                    removevars_rest = [];
                    removevars_task = [];
                
                    for z = 1:length(allvars_rest)
                    
                        if length(find(cifti_rest_thresh_dat == allvars_rest(z))) < minsize
                        
                            removevars_rest = [removevars_rest allvars_rest(z)];
                        
                        end
                    end
                
                    for h = 1:length(cifti_rest_thresh_dat)
                    
                        if ismember(cifti_rest_thresh_dat(h),removevars_rest)
                        
                            cifti_rest_thresh_dat(h) = 0;
                        
                        end
                    end
                
                    for z = 1:length(allvars_task)
                    
                        if length(find(cifti_task_thresh_dat == allvars_task(z))) < minsize
                        
                            removevars_task = [removevars_task allvars_task(z)];
                        
                        end
                    end
                
                    for h = 1:length(cifti_task_thresh_dat)
                    
                        if ismember(cifti_task_thresh_dat(h),removevars_task)
                        
                            cifti_task_thresh_dat(h) = 0;
                        
                        end
                    end
                end
        
                dcorrdata = [];
        
                for q = 1:length(cifti_rest_thresh_dat)
            
                    if cifti_rest_thresh_dat(q) > 0 && cifti_task_thresh_dat(q) > 0
                
                        dcorrdata = [dcorrdata;1 1];
                
                    elseif cifti_rest_thresh_dat(q) > 0
                
                        dcorrdata = [dcorrdata;1 0];
                
                    elseif cifti_task_thresh_dat(q) > 0
                
                        dcorrdata = [dcorrdata;0 1];
                
                    end
            
                end
        
                dc = dice_coefficient_mod(dcorrdata(:,1),dcorrdata(:,2));
                
                if strcmp(comp, 'task-task') || strcmp(comp, 'matched-samedata-task-task')
                
                    tasktaskpearsoncoefstemp = [tasktaskpearsoncoefstemp r(1,2)];
                    tasktaskdicecoefstemp = [tasktaskdicecoefstemp dc];
                    
                elseif strcmp(comp, 'task-difftask') || strcmp(comp, 'matched-samedata-task-difftask')
                    
                    taskdifftaskpearsoncoefstemp = [taskdifftaskpearsoncoefstemp r(1,2)];
                    taskdifftaskdicecoefstemp = [taskdifftaskdicecoefstemp dc];
                    
                elseif strcmp(comp, 'task-rest') || strcmp(comp, 'matched-samedata-task-rest')
                    
                    taskrestpearsoncoefstemp = [taskrestpearsoncoefstemp r(1,2)];
                    taskrestdicecoefstemp = [taskrestdicecoefstemp dc];
                    
                    if strcmp(task, 'mem') && randomizevals == 1
                        
                       memfilestask = [memfilestask; cifti_task_thresh_dat'];
                       memfilesrest = [memfilesrest; cifti_rest_thresh_dat'];
                       
                    elseif strcmp(task, 'mixed') && randomizevals == 1
                        
                       mixedfilestask = [mixedfilestask; cifti_task_thresh_dat'];
                       mixedfilesrest = [mixedfilesrest; cifti_rest_thresh_dat']; 
                       
                    elseif strcmp(task, 'motor') && randomizevals == 1
                        
                       motorfilestask = [motorfilestask; cifti_task_thresh_dat'];
                       motorfilesrest = [motorfilesrest; cifti_rest_thresh_dat'];
                       
                    end
                    
                elseif strcmp(comp, 'rest-rest') || strcmp(comp, 'matched-samedata-rest-rest')
                    
                    restrestpearsoncoefstemp = [restrestpearsoncoefstemp r(1,2)];
                    restrestdicecoefstemp = [restrestdicecoefstemp dc];
                    
                elseif strcmp(comp, 'taskrest-taskrest') || strcmp(comp, 'matched-samedata-taskrest-taskrest')
                    
                    taskresttaskrestpearsoncoefstemp = [taskresttaskrestpearsoncoefstemp r(1,2)];
                    taskresttaskrestdicecoefstemp = [taskresttaskrestdicecoefstemp dc];
                    
                elseif strcmp(comp, 'alltask-alltask') || strcmp(comp, 'matched-samedata-alltask-alltask')
                    
                    alltaskalltaskpearsoncoefstemp = [alltaskalltaskpearsoncoefstemp r(1,2)];
                    alltaskalltaskdicecoefstemp = [alltaskalltaskdicecoefstemp dc];
                    
                elseif strcmp(comp, 'alltask-rest') || strcmp(comp, 'matched-samedata-alltask-rest')
                    
                    alltaskrestpearsoncoefstemp = [alltaskrestpearsoncoefstemp r(1,2)];
                    alltaskrestdicecoefstemp = [alltaskrestdicecoefstemp dc];
                    
                    if randomizevals == 1
                        
                        alltaskfilestask = [alltaskfilestask; cifti_task_thresh_dat];
                        alltaskfilesrest = [alltaskfilesrest; cifti_rest_thresh_dat];
                        
                    end
                    
                elseif strcmp(comp, 'alltaskrest-alltaskrest')
                    
                    alltaskrestalltaskrestpearsoncoefstemp = [alltaskrestalltaskrestpearsoncoefstemp r(1,2)];
                    alltaskrestalltaskrestdicecoefstemp = [alltaskrestalltaskrestdicecoefstemp dc];
                    
                end
                    
                    
                pearsoncoefstemp = [pearsoncoefstemp r(1,2)];
            
                dicecoefstemp = [dicecoefstemp dc];
            
            else
            
                for a = 1:numel(thresholds)
                    
                    threshoutputtemp = [];
                
                    threshold = thresholds(a);
                
                    subject = subs{s};
                    task_file = task_files{t};
                    rest_file = rest_files{t};
                    task = tasks1{t};
                    task2 = tasks2{t};
        
                    taskoutputtemp = [taskoutputtemp task task2];

                    cifti_task = ft_read_cifti_mod(task_file);
                    cifti_rest = ft_read_cifti_mod(rest_file);
                    
                    if a == 1
                    
                        r = corrcoef(cifti_rest.data, cifti_task.data);
                        
                        pearsoncoefstemp = [pearsoncoefstemp r(1,2)];
                    
                        if strcmp(comp, 'task-task') || strcmp(comp, 'matched-samedata-task-task')
                
                            tasktaskpearsoncoefstemp = [tasktaskpearsoncoefstemp r(1,2)];
                            
                        elseif strcmp(comp, 'task-difftask') || strcmp(comp, 'matched-samedata-task-difftask')
                            
                            taskdifftaskpearsoncoefstemp = [taskdifftaskpearsoncoefstemp r(1,2)];
                    
                        elseif strcmp(comp, 'task-rest') || strcmp(comp, 'matched-samedata-task-rest')
                    
                            taskrestpearsoncoefstemp = [taskrestpearsoncoefstemp r(1,2)];
                    
                        elseif strcmp(comp, 'rest-rest') || strcmp(comp, 'matched-samedata-rest-rest')
                    
                            restrestpearsoncoefstemp = [restrestpearsoncoefstemp r(1,2)];
                            
                        elseif strcmp(comp, 'taskrest-taskrest') || strcmp(comp, 'matched-samedata-taskrest-taskrest')
                    
                            taskresttaskrestpearsoncoefstemp = [taskresttaskrestpearsoncoefstemp r(1,2)];
                    
                        elseif strcmp(comp, 'alltask-alltask') || strcmp(comp, 'matched-samedata-alltask-alltask')
                    
                            alltaskalltaskpearsoncoefstemp = [alltaskalltaskpearsoncoefstemp r(1,2)];
                    
                        elseif strcmp(comp, 'alltask-rest') || strcmp(comp, 'matched-samedata-alltask-rest')
                    
                            alltaskrestpearsoncoefstemp = [alltaskrestpearsoncoefstemp r(1,2)];
                            
                        elseif strcmp(comp, 'alltaskrest-alltaskrest')
                            
                            alltaskrestalltaskrestpearsoncoefstemp = [alltaskrestalltaskrestpearsoncoefstemp r(1,2)];
                    
                        end
                    
                    end
    
                    if absolute == 1
        
                        cifti_task_threshold = find(cifti_task.data < threshold);
                        cifti_rest_threshold = find(cifti_rest.data < threshold);

                    else
        
                        cifti_task_threshold = find(cifti_task.data < prctile(cifti_task.data,threshold));
                        cifti_rest_threshold = find(cifti_rest.data < prctile(cifti_rest.data,threshold));
                
                    end
                
                    cifti_task_thresh_dat = zeros(size(cifti_task.data));
                    cifti_task_thresh_dat(cifti_task_threshold,1) = 1;
                
                    cifti_rest_thresh_dat = zeros(size(cifti_rest.data));
                    cifti_rest_thresh_dat(cifti_rest_threshold,1) = 1;
        
                    dc = dice_coefficient_mod(cifti_task_thresh_dat,cifti_rest_thresh_dat);
                    
                    if strcmp(comp, 'task-task') || strcmp(comp, 'matched-samedata-task-task')
                
                        tasktaskdicecoefstemp = [tasktaskdicecoefstemp dc];
                        
                    elseif strcmp(comp, 'task-difftask') || strcmp(comp, 'matched-samedata-task-difftask')
                            
                        taskdifftaskdicecoefstemp = [taskdifftaskdicecoefstemp dc];
                    
                    elseif strcmp(comp, 'task-rest') || strcmp(comp, 'matched-samedata-task-rest')
                    
                        taskrestdicecoefstemp = [taskrestdicecoefstemp dc];
                        
                        if strcmp(task, 'mem') && randomizevals == 1 && thresholds(a) == 10
                        
                            memfilestask = [memfilestask; cifti_task_thresh_dat'];
                            memfilesrest = [memfilesrest; cifti_rest_thresh_dat'];
                       
                        elseif strcmp(task, 'mixed') && randomizevals == 1 && thresholds(a) == 10
                        
                            mixedfilestask = [mixedfilestask; cifti_task_thresh_dat'];
                            mixedfilesrest = [mixedfilesrest; cifti_rest_thresh_dat']; 
                       
                        elseif strcmp(task, 'motor') && randomizevals == 1 && thresholds(a) == 10
                        
                            motorfilestask = [motorfilestask; cifti_task_thresh_dat'];
                            motorfilesrest = [motorfilesrest; cifti_rest_thresh_dat'];
                       
                        end
                    
                    elseif strcmp(comp, 'rest-rest') || strcmp(comp, 'matched-samedata-rest-rest')
                    
                        restrestdicecoefstemp = [restrestdicecoefstemp dc];
                        
                    elseif strcmp(comp, 'taskrest-taskrest') || strcmp(comp, 'matched-samedata-taskrest-taskrest')
                        
                        taskresttaskrestdicecoefstemp = [taskresttaskrestdicecoefstemp dc];
                    
                    elseif strcmp(comp, 'alltask-alltask') || strcmp(comp, 'matched-samedata-alltask-alltask')
                    
                        alltaskalltaskdicecoefstemp = [alltaskalltaskdicecoefstemp dc];
                    
                    elseif strcmp(comp, 'alltask-rest') || strcmp(comp, 'matched-samedata-alltask-rest')
                    
                        alltaskrestdicecoefstemp = [alltaskrestdicecoefstemp dc];
                        
                        if randomizevals == 1 && thresholds(a) == 10
                        
                            alltaskfilestask = [alltaskfilestask; cifti_task_thresh_dat'];
                            alltaskfilesrest = [alltaskfilesrest; cifti_rest_thresh_dat'];
                        
                        end
                        
                    elseif strcmp(comp, 'alltaskrest-alltaskrest')
                        
                        alltaskrestalltaskrestdicecoefstemp = [alltaskrestalltaskrestdicecoefstemp dc];
                    
                    end
            
                    dicecoefstemp = [dicecoefstemp dc];
                
                    threshoutputtemp = [threshoutputtemp thresholds(a)];
            
                end
           
                threshoutput = [threshoutput; threshoutputtemp];
            
            end
        
            taskoutput = [taskoutput; taskoutputtemp];
            dicecoefs = [dicecoefs; dicecoefstemp];
            pearsoncoefs = [pearsoncoefs; pearsoncoefstemp];
        
            if plotcorrs == 1
                
                if strcmp(task, 'mem')
            
                	label = 'Memory';
            
                elseif strcmp(task, 'mixed')
            
                	label = 'Mixed';
            
                elseif strcmp(task, 'motor')
            
                	label = 'Motor';
                    
                elseif strcmp(task, 'rest')
                    
                    label = 'Rest';
                    
                elseif strcmp(task, 'all')
                    
                    label = 'All Tasks';
            
                end
  
            plotlabels = [plotlabels {label}];
            plotdata = [plotdata; dicecoefstemp];
        
            end
        end
    
        if plotcorrs == 1
            
            if strcmp(comp, 'task-task') || strcmp(comp, 'task-rest') || strcmp(comp, 'matched-samedata-task-task') || strcmp(comp, 'matched-samedata-task-rest')
    
                plot(thresholds, plotdata(1,:), thresholds, plotdata(2,:), thresholds, plotdata(3,:))
                
            else
                
                plot(thresholds, plotdata)
                
            end
            
            if strcmp(comp, 'task-task')
                
                title([subject ' Task-Task Dice Correlations'], 'fontsize',18)
                
            elseif strcmp(comp, 'task-rest')
                
                title([subject ' Task-Rest Dice Correlations'], 'fontsize',18)
                
            elseif strcmp(comp, 'rest-rest')
                
                title([subject ' Rest-Rest Dice Correlations'], 'fontsize',18)
                
            elseif strcmp(comp, 'alltask-alltask')
                
                title([subject ' All Tasks-All Tasks Dice Correlations'], 'fontsize',18)
                
            elseif strcmp(comp, 'alltask-rest')
                
                title([subject ' All Tasks-Rest Dice Correlations'], 'fontsize',18)
                
            elseif strcmp(comp, 'matched-samedata-task-task')
                
                title([subject ' Task-Task Dice Correlations Matched'], 'fontsize',18)
                
            elseif strcmp(comp, 'matched-samedata-task-rest')
                
                title([subject ' Task-Rest Dice Correlations Matched'], 'fontsize',18)
                
            elseif strcmp(comp, 'matched-samedata-rest-rest')
                
                title([subject ' Rest-Rest Dice Correlations Matched'], 'fontsize',18)
                
            elseif strcmp(comp, 'matched-samedata-alltask-alltask')
                
                title([subject ' All Tasks-All Tasks Dice Correlations Matched'], 'fontsize',18)
                
            elseif strcmp(comp, 'matched-samedata-alltask-rest')
                
                title([subject ' All Tasks-Rest Dice Correlations Matched'], 'fontsize',18)
                
            end
        
            if absolute == 1
            
                xlabel('Absolute Thresholds', 'fontsize',18)
            
            else
            
                xlabel('Percent Thresholds', 'fontsize',18)
            
            end
        
         
            ylabel('Correlations', 'fontsize',18)
            legend(plotlabels, 'fontsize',18, 'Location', 'southeast')
            ylim([.3 1]);
        
            if absolute == 1
                
                if saveindividplots == 1 && SNRMasks == 1
                    
                    saveas(gcf,[outputdir '/' subject '_absolutethreshold_' task 'SNRMask_dicecorrelations.jpg'])
                
                elseif saveindividplots == 1
        
                    saveas(gcf,[outputdir '/' subject '_absolutethreshold_' task 'dicecorrelations.jpg'])
                    
                end
                
                if concatplots == 0
                
                    close gcf
                
                end
        
            else
                
                if saveindividplots == 1 && SNRMasks == 1
                    
                    saveas(gcf,[outputdir '/' subject '_percentthreshold_' task 'SNRMask_dicecorrelations.jpg'])
                
                elseif saveindividplots == 1
            
                    saveas(gcf,[outputdir '/' subject '_percentthreshold_' task 'dicecorrelations.jpg'])
                    
                end
                
                if concatplots == 0
                
                    close gcf
                    
                end
            
            end
            
            if concatplots == 1
                
                if absolute == 1 && SNRMasks == 1
                    
                    export_fig([outputdir '/AbsoluteThreshold_AllPlots_SNRMask_dicecorrelations.jpg'], '-pdf', '-append')
                    close gcf
                
                elseif absolute == 1
                    
                    export_fig([outputdir '/AbsoluteThreshold_AllPlots_dicecorrelations.jpg'], '-pdf', '-append')
                    close gcf
                    
                elseif SNRMasks == 1
                    
                    export_fig([outputdir '/PercentThreshold_AllPlots_SNRMask_dicecorrelations.jpg'], '-pdf', '-append')
                    close gcf
                    
                else
                
                    export_fig([outputdir '/PercentThreshold_AllPlots_dicecorrelations.jpg'], '-pdf', '-append')
                    close gcf
                
                end
            end
        
        end
        
        tasktaskpearsoncoefs = [tasktaskpearsoncoefs; tasktaskpearsoncoefstemp];
        taskdifftaskpearsoncoefs = [taskdifftaskpearsoncoefs; taskdifftaskpearsoncoefstemp];
        taskrestpearsoncoefs = [taskrestpearsoncoefs; taskrestpearsoncoefstemp];
        restrestpearsoncoefs = [restrestpearsoncoefs; restrestpearsoncoefstemp];
        taskresttaskrestpearsoncoefs = [taskresttaskrestpearsoncoefs; taskresttaskrestpearsoncoefstemp];
        alltaskrestpearsoncoefs = [alltaskrestpearsoncoefs; alltaskrestpearsoncoefstemp];
        alltaskalltaskpearsoncoefs = [alltaskalltaskpearsoncoefs; alltaskalltaskpearsoncoefstemp];
        alltaskrestalltaskrestpearsoncoefs = [alltaskrestalltaskrestpearsoncoefs; alltaskrestalltaskrestpearsoncoefstemp];

        tasktaskdicecoefs = [tasktaskdicecoefs; tasktaskdicecoefstemp];
        taskdifftaskdicecoefs = [taskdifftaskdicecoefs; taskdifftaskdicecoefstemp];
        taskrestdicecoefs = [taskrestdicecoefs; taskrestdicecoefstemp];
        restrestdicecoefs = [restrestdicecoefs; restrestdicecoefstemp];
        taskresttaskrestdicecoefs = [taskresttaskrestdicecoefs; taskresttaskrestdicecoefstemp];
        alltaskrestdicecoefs = [alltaskrestdicecoefs; alltaskrestdicecoefstemp];
        alltaskalltaskdicecoefs = [alltaskalltaskdicecoefs; alltaskalltaskdicecoefstemp];
        alltaskrestalltaskrestdicecoefs = [alltaskrestalltaskrestdicecoefs; alltaskrestalltaskrestdicecoefstemp];
      
    end
end

if isempty(alltaskrestalltaskrestpearsoncoefs) && numel(thresholds) > 1

    exportdata = [tasktaskpearsoncoefs(:,1) tasktaskpearsoncoefs(:,2) tasktaskpearsoncoefs(:,3) taskdifftaskpearsoncoefs(:,2) taskdifftaskpearsoncoefs(:,3) taskdifftaskpearsoncoefs(:,6) taskrestpearsoncoefs(:,1) taskrestpearsoncoefs(:,2) taskrestpearsoncoefs(:,3) taskresttaskrestpearsoncoefs(:,1) taskresttaskrestpearsoncoefs(:,2) taskresttaskrestpearsoncoefs(:,3) restrestpearsoncoefs alltaskrestpearsoncoefs alltaskalltaskpearsoncoefs tasktaskdicecoefs(:,2) tasktaskdicecoefs(:,9) tasktaskdicecoefs(:,16) taskdifftaskdicecoefs(:,9) taskdifftaskdicecoefs(:,16) taskdifftaskdicecoefs(:,37) taskrestdicecoefs(:,2) taskrestdicecoefs(:,9) taskrestdicecoefs(:,16) taskresttaskrestdicecoefs(:,2) taskresttaskrestdicecoefs(:,9) taskresttaskrestdicecoefs(:,16) restrestdicecoefs(:,2) alltaskrestdicecoefs(:,2) alltaskalltaskdicecoefs(:,2)];
    
    if fishertransform == 1
        
        for p = 1:size(exportdata,2)
            
            if p == 3 || p == 5 || p == 6 || p == 9 || p == 12 || p == 18 || p == 20 || p == 21 || p == 24 || p == 27
                
                exportdata([1:7 9],p) = atanh(exportdata([1:7 9],p));
                
            else
            
                exportdata(:,p) = atanh(exportdata(:,p));
            
            end
        end
    end

else
    
    exportdata = [tasktaskpearsoncoefs(:,1) tasktaskpearsoncoefs(:,2) tasktaskpearsoncoefs(:,3) taskdifftaskpearsoncoefs(:,2) taskdifftaskpearsoncoefs(:,3) taskdifftaskpearsoncoefs(:,6) taskrestpearsoncoefs(:,1) taskrestpearsoncoefs(:,2) taskrestpearsoncoefs(:,3) taskresttaskrestpearsoncoefs(:,1) taskresttaskrestpearsoncoefs(:,2) taskresttaskrestpearsoncoefs(:,3) restrestpearsoncoefs alltaskrestpearsoncoefs alltaskalltaskpearsoncoefs alltaskrestalltaskrestpearsoncoefs tasktaskdicecoefs(:,2) tasktaskdicecoefs(:,9) tasktaskdicecoefs(:,16) taskdifftaskdicecoefs(:,9) taskdifftaskdicecoefs(:,16) taskdifftaskdicecoefs(:,37) taskrestdicecoefs(:,2) taskrestdicecoefs(:,9) taskrestdicecoefs(:,16) taskresttaskrestdicecoefs(:,2) taskresttaskrestdicecoefs(:,9) taskresttaskrestdicecoefs(:,16) restrestdicecoefs(:,2) alltaskrestdicecoefs(:,2) alltaskalltaskdicecoefs(:,2) alltaskrestalltaskrestdicecoefs(:,2)];
    
    if fishertransform == 1
        
        for p = 1:size(exportdata,2)
            
            if p == 3 || p == 5 || p == 6 || p == 9 || p == 12 || p == 19 || p == 21 || p == 22 || p == 25 || p == 28
                
                exportdata([1:7 9],p) = atanh(exportdata([1:7 9],p));
                
            else
            
                exportdata(:,p) = atanh(exportdata(:,p));
            
            end
        end
    end
    
end



if randomizevals == 1
    
    nsim = 10000;
    nsubs = size(memfilestask,1);
    
    motorfilestask(8,:) = [];       %% Eliminate MSC09 from motor task
    motorfilesrest(8,:) = [];
    
    memtasksims = zeros(1,nsim);
    mixedtasksims = zeros(1,nsim);
    motortasksims = zeros(1,nsim);
    alltasksims = zeros(1,nsim);

    
    for x = 1:nsim
        
        
        memtaskdc = [];
        mixedtaskdc = [];
        motortaskdc = [];
        alltaskdc = [];
    
        rng('shuffle');  %% Reset random number generator on each iteration
        
        % Randomly shuffle rows (subjects) for each task and rest matrix
        
        memtaskshuffle = memfilestask(randperm(size(memfilestask,1)),:);
        mixedtaskshuffle = mixedfilestask(randperm(size(mixedfilestask,1)),:);
        motortaskshuffle = motorfilestask(randperm(size(motorfilestask,1)),:);
        memrestshuffle = memfilesrest(randperm(size(memfilesrest,1)),:);
        mixedrestshuffle = mixedfilesrest(randperm(size(mixedfilesrest,1)),:);
        motorrestshuffle = motorfilesrest(randperm(size(motorfilesrest,1)),:);
        alltasktaskshuffle = alltaskfilestask(randperm(size(alltaskfilestask,1)),:);
        alltaskrestshuffle = alltaskfilesrest(randperm(size(alltaskfilesrest,1)),:);
        
    
        for d = 1:size(memtaskshuffle,1)
            
            memtaskdc = [memtaskdc; dice_coefficient_mod(memtaskshuffle(d,:), memrestshuffle(d,:))];
            mixedtaskdc = [mixedtaskdc; dice_coefficient_mod(mixedtaskshuffle(d,:), mixedrestshuffle(d,:))];
            alltaskdc = [alltaskdc; dice_coefficient_mod(alltasktaskshuffle(d,:), alltaskrestshuffle(d,:))];
            
            if d < size(memtaskshuffle,1) %% Only 8 subjects for motor task
                
                motortaskdc = [motortaskdc; dice_coefficient_mod(motortaskshuffle(d,:), motorrestshuffle(d,:))];
                
            end
            
        end
        
        if fishertransform == 1
            
            memtaskdc = atanh(memtaskdc);
            mixedtaskdc = atanh(mixedtaskdc);
            motortaskdc = atanh(motortaskdc);
            alltaskdc = atanh(alltaskdc);
            
        end
        
        memtasksims(x) = mean(memtaskdc);
        mixedtasksims(x) = mean(mixedtaskdc);
        motortasksims(x) = mean(motortaskdc);
        alltasksims(x) = mean(alltaskdc);
    
    end
    
    memtaskpval = prctile(memtasksims, 95);
    mixedtaskpval = prctile(mixedtasksims, 95);
    motortaskpval = prctile(motortasksims, 95);
    alltaskpval = prctile(alltasksims, 95);
    
end








if plotavgcorrs == 1
    
    for x = 1:numel(comparison)
        
        comp = comparison{x};
        numthresholds = numel(thresholds);
        
        if strcmp(comp, 'task-task') || strcmp(comp, 'matched-samedata-task-task')
            
            for t = 1:2 % Plot pearson and dice correlations
                
                for y = 1:3  %% loop through tasks
                
                % Plot task-task comparisons
                
                if t == 1
                    
                    if isempty(taskresttaskrestdicecoefs)
                    
                        plotdata = [tasktaskpearsoncoefs(:,y) taskrestpearsoncoefs(:,y) restrestpearsoncoefs];
                        
                    else
                        
                        plotdata = [tasktaskpearsoncoefs(:,y) taskrestpearsoncoefs(:,y) taskresttaskrestpearsoncoefs(:,y)];
                        
                    end
                    
                else
                    
                    if isempty(taskresttaskrestdicecoefs)
                
                        plotdata = [tasktaskdicecoefs(:,(numthresholds*(y-1))+2) taskrestdicecoefs(:,(numthresholds*(y-1))+2) restrestdicecoefs(:,2)];
                        
                    else
                        
                        plotdata = [tasktaskdicecoefs(:,(numthresholds*(y-1))+2) taskrestdicecoefs(:,(numthresholds*(y-1))+2) taskresttaskrestdicecoefs(:,(numthresholds*(y-1))+2)];
                    
                    end
                end
                
                SEplot = zeros(1,size(plotdata,2));
                
                if y == 3   %% Exclude MSC09 from motor task
                    
                	if fishertransform == 1
                            
                    	for w = 1:size(plotdata,2)
                            
                        	plotdata([1:7 9],w) = atanh(plotdata([1:7 9],w));
                                
                        end
                    end
                    
                    meanplotdata = mean(plotdata([1:7 9],:),1);
                    
                else
                    
                	if fishertransform == 1
                            
                    	for w = 1:size(plotdata,2)
                            
                        	plotdata(:,w) = atanh(plotdata(:,w));
                                
                        end
                    end
                    
                    meanplotdata = mean(plotdata,1);
                    
                end
                
                for z = 1:size(plotdata,2)
                    
                    if y == 3   %% Exclude MSC09 from motor task
                        
                        SEplot(z) = std(plotdata([1:7 9],z))/(sqrt(size(plotdata,1)-1));
                        
                    else
                
                        SEplot(z) = std(plotdata(:,z))/sqrt(size(plotdata,1));
                        
                    end
                end
                
                if t == 1
                    
                    titlelabel = 'Pearson';
                    
                else
                    
                    titlelabel = 'Dice';
                    
                end
                
                if t == 2 && y == 1
                    
                    taskrestpval = memtaskpval;
                    
                elseif t == 2 && y == 2
                    
                    taskrestpval = mixedtaskpval;
                    
                elseif t == 2 && y == 3
                    
                    taskrestpval = motortaskpval;
                    
                end
                
                %barwitherr(SEplot, categorical({'Task-Task', 'Task-Rest', 'Rest-Rest'}), meanplotdata);
                bar(categorical({'Task-Task', 'Task-Rest', 'Rest-Rest'}), meanplotdata, 'FaceColor', 'w');
                hold on
                errorbar(categorical({'Task-Task', 'Task-Rest', 'Rest-Rest'}),meanplotdata,SEplot,'.')
                hold on
                
                if t == 2
                
                    plot(xlim,[taskrestpval taskrestpval], 'LineWidth', 3, 'color', [0, 0, 0])
                    hold on
                    
                end
                
                plot(categorical({'Task-Task', 'Task-Rest', 'Rest-Rest'}), plotdata(9,:), '.', 'MarkerSize', 15, 'color', [1, 0.5, 0]);
                hold on
                
                if y ~= 3
                    
                    plot(categorical({'Task-Task', 'Task-Rest', 'Rest-Rest'}), plotdata(8,:), '.', 'MarkerSize', 15, 'color', [0, 0.6, 0.6]);
                    hold on
                    
                end
                
                plot(categorical({'Task-Task', 'Task-Rest', 'Rest-Rest'}), plotdata(7,:), '.', 'MarkerSize', 15, 'color', [1, 0, 1]);
                hold on
                plot(categorical({'Task-Task', 'Task-Rest', 'Rest-Rest'}), plotdata(6,:), '.', 'MarkerSize', 15, 'color', [0.2, 1, 1]);
                hold on
                plot(categorical({'Task-Task', 'Task-Rest', 'Rest-Rest'}), plotdata(5,:), '.', 'MarkerSize', 15, 'color', [0, 0, 1]);
                hold on
                plot(categorical({'Task-Task', 'Task-Rest', 'Rest-Rest'}), plotdata(4,:), '.', 'MarkerSize', 15, 'color', [1, 0, 0]);
                hold on
                plot(categorical({'Task-Task', 'Task-Rest', 'Rest-Rest'}), plotdata(3,:), '.', 'MarkerSize', 15, 'color', [0, 1, 0]);
                hold on
                plot(categorical({'Task-Task', 'Task-Rest', 'Rest-Rest'}), plotdata(2,:), '.', 'MarkerSize', 15, 'color', [0.9, 0.9, 0]);
                hold on
                plot(categorical({'Task-Task', 'Task-Rest', 'Rest-Rest'}), plotdata(1,:), '.', 'MarkerSize', 15, 'color', [0, 0, 0]);
                hold on

                
                if strcmp(comp, 'task-task') && y == 1
                
                    title(['All Subjects Memory Task-Task ' titlelabel ' Correlations'], 'fontsize',18)
                    
                elseif strcmp(comp, 'task-task') && y == 2
                
                    title(['All Subjects Mixed Task-Task ' titlelabel ' Correlations'], 'fontsize',18)
                    
                elseif strcmp(comp, 'task-task') && y == 3
                
                    title(['All Subjects Motor Task-Task ' titlelabel ' Correlations'], 'fontsize',18)
                    
                elseif strcmp(comp, 'matched-samedata-task-task') && y == 1
                    
                    title(['All Subjects Memory Task-Task ' titlelabel ' Correlations Matched'], 'fontsize',18)
                    
                elseif strcmp(comp, 'matched-samedata-task-task') && y == 2
                    
                    title(['All Subjects Mixed Task-Task ' titlelabel ' Correlations Matched'], 'fontsize',18)
                    
                elseif strcmp(comp, 'matched-samedata-task-task') && y == 3
                    
                    title(['All Subjects Motor Task-Task ' titlelabel ' Correlations Matched'], 'fontsize',18)
                    
                end
                
                if fishertransform == 1
                    
                    ylabel('Fisher Transformed Correlations', 'fontsize',18)
                    ylim([.3 2]);
                    
                else
                    
                    ylabel('Correlations', 'fontsize',18)
                    ylim([.3 1]);
                    
                end
                
                h = findobj(gca,'Type','line');
                
                if fishertransform == 1
                    
                    legendpos = 'NorthEast';
                    
                else
                    
                    legendpos = 'SouthEast';
                    
                end
                
                if y == 3
                    
                    legend(h(1:8), 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC10', 'Location', legendpos)
                    
                else
                
                    legend(h(1:9), 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', legendpos)
                    
                end
                

                
                if saveindividplots == 1
            
                    if strcmp(comp, 'task-task') && y == 1
                        
                        if fishertransform == 1 && SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_memory_' titlelabel 'SNRMask__fishercorrelations.jpg'])
                        
                        elseif fishertransform == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_memory_' titlelabel '_fishercorrelations.jpg'])
                            
                        elseif SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_memory_' titlelabel 'SNRMask_correlations.jpg'])
                            
                        else
                
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_memory_' titlelabel '_correlations.jpg'])
                        
                        end
                    
                    elseif strcmp(comp, 'task-task') && y == 2
                        
                        if fishertransform == 1 && SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_mixed_' titlelabel 'SNRMask_fishercorrelations.jpg'])
                        
                        elseif fishertransform == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_mixed_' titlelabel '_fishercorrelations.jpg'])
                            
                        elseif SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_mixed_' titlelabel 'SNRMask_correlations.jpg'])
                            
                        else
                
                        	saveas(gcf,[outputdir '/AllSubjects_percentthreshold_mixed_' titlelabel '_correlations.jpg'])
                        
                        end
                    
                    elseif strcmp(comp, 'task-task') && y == 3
                        
                        if fishertransform == 1 && SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_motor_' titlelabel 'SNRMask_fishercorrelations.jpg'])
                        
                        elseif fishertransform == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_motor_' titlelabel '_fishercorrelations.jpg'])
                            
                        elseif SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_motor_' titlelabel 'SNRMask_correlations.jpg'])
                            
                        else
                
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_motor_' titlelabel '_correlations.jpg'])
                        
                        end
                    
                    elseif strcmp(comp, 'matched-samedata-task-task') && y == 1
                        
                        if fishertransform == 1 && SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_memory_matcheddata_' titlelabel 'SNRMask_fishercorrelations.jpg'])
                        
                        elseif fishertransform == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_memory_matcheddata_' titlelabel '_fishercorrelations.jpg'])
                            
                        elseif SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_memory_matcheddata_' titlelabel 'SNRMask_correlations.jpg'])
                            
                        else
                    
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_memory_matcheddata_' titlelabel '_correlations.jpg'])
                            
                        end
                    
                    elseif strcmp(comp, 'matched-samedata-task-task') && y == 2
                        
                        if fishertransform == 1 && SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_mixed_matcheddata_' titlelabel 'SNRMask_fishercorrelations.jpg'])
                        
                        elseif fishertransform == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_mixed_matcheddata_' titlelabel '_fishercorrelations.jpg'])
                            
                        elseif SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_mixed_matcheddata_' titlelabel 'SNRMask_correlations.jpg'])
                            
                        else
                    
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_mixed_matcheddata_' titlelabel '_correlations.jpg'])
                            
                        end
                    
                    elseif strcmp(comp, 'matched-samedata-task-task') && y == 3
                        
                        if fishertransform == 1 && SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_motor_matcheddata_' titlelabel 'SNRMask_fishercorrelations.jpg'])
                        
                        elseif fishertransform == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_motor_matcheddata_' titlelabel '_fishercorrelations.jpg'])
                            
                        elseif SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_motor_matcheddata_' titlelabel 'SNRMask_correlations.jpg'])
                                                  
                        else
                    
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_motor_matcheddata_' titlelabel '_correlations.jpg'])
                            
                        end
                    
                    end
                end
                
                if concatplots == 0
                
                    close gcf
                    
                end
                
                if concatplots == 1
                    
                    if fishertransform == 1 && SNRMasks == 1
                        
                        export_fig([outputdir '/PercentThreshold_AllPlots_matcheddata_SNRMask_fishercorrelations.pdf'], '-pdf', '-append')
                    
                    elseif fishertransform == 1
                        
                        export_fig([outputdir '/PercentThreshold_AllPlots_matcheddata_fishercorrelations.pdf'], '-pdf', '-append')
                        
                    elseif SNRMasks == 1
                        
                        export_fig([outputdir '/PercentThreshold_AllPlots_matcheddata_SNRMask_correlations.pdf'], '-pdf', '-append')
                        
                    else
                    
                        export_fig([outputdir '/PercentThreshold_AllPlots_matcheddata_correlations.pdf'], '-pdf', '-append')
                        
                    end
                    
                    close gcf 
                    
                end
                    
                end
            end
            
        elseif strcmp(comp, 'task-difftask') || strcmp(comp, 'matched-samedata-task-difftask')
            
            for t = 1:2 % Plot pearson and dice correlations
                
                for y = 1:3  %% loop through tasks
                
                % Plot task-task comparisons
                
                if t == 1
                    
                    if y == 1
                    
                        plotdata = [taskdifftaskpearsoncoefs(:,1) taskdifftaskpearsoncoefs(:,2) taskdifftaskpearsoncoefs(:,3)];
                        
                    elseif y == 2
                        
                        plotdata = [taskdifftaskpearsoncoefs(:,5) taskdifftaskpearsoncoefs(:,4) taskdifftaskpearsoncoefs(:,6)];
                        
                    elseif y == 3
                        
                        plotdata = [taskdifftaskpearsoncoefs(:,9) taskdifftaskpearsoncoefs(:,7) taskdifftaskpearsoncoefs(:,8)];
                        
                    end
                    
                else
                    
                    if y == 1
                
                        plotdata = [taskdifftaskdicecoefs(:,2) taskdifftaskdicecoefs(:,9) taskdifftaskdicecoefs(:,16)];
                        
                    elseif y == 2
                        
                        plotdata = [taskdifftaskdicecoefs(:,30) taskdifftaskdicecoefs(:,23) taskdifftaskdicecoefs(:,37)];
                        
                    elseif y == 3
                        
                        plotdata = [taskdifftaskdicecoefs(:,58) taskdifftaskdicecoefs(:,44) taskdifftaskdicecoefs(:,51)];
                    
                    end
                end
                
                SEplot = zeros(1,size(plotdata,2));
                
                if y == 3   %% Exclude MSC09 from motor task
                    
                    if fishertransform == 1
                            
                    	for w = 1:size(plotdata,2)
                            
                        	plotdata([1:7 9],w) = atanh(plotdata([1:7 9],w));
                                
                        end
                    end
                    
                    meanplotdata = mean(plotdata([1:7 9],:),1);
                    
                else
                    
                	if fishertransform == 1
                            
                    	for w = 1:size(plotdata,2)
                            
                        	plotdata(:,w) = atanh(plotdata(:,w));
                                
                        end
                    end
                    
                    meanplotdata = mean(plotdata,1);
                    
                end
                
                for z = 1:size(plotdata,2)
                    
                    if y == 3   %% Exclude MSC09 from motor task
                        
                        SEplot(z) = std(plotdata([1:7 9],z))/(sqrt(size(plotdata,1)-1));
                        
                    else
                
                        SEplot(z) = std(plotdata(:,z))/sqrt(size(plotdata,1));
                        
                    end
                end
                
                if t == 1
                    
                    titlelabel = 'Pearson';
                    
                else
                    
                    titlelabel = 'Dice';
                    
                end
                
                %barwitherr(SEplot, categorical({'Task-Task', 'Task-Rest', 'Rest-Rest'}), meanplotdata);
                
                if y == 1
                    
                    c = categorical({'Memory-Memory', 'Memory-Mixed', 'Memory-Motor'});
                    
                elseif y == 2
                    
                    c = categorical({'Mixed-Mixed', 'Mixed-Memory', 'Mixed-Motor'});
                    
                elseif y == 3
                    
                    c = categorical({'Motor-Motor', 'Motor-Memory', 'Motor-Mixed'});
                    
                end
                
                bar(c, meanplotdata, 'FaceColor', 'w');
                hold on
                errorbar(c,meanplotdata,SEplot,'.')
                plot(c, plotdata(9,:), '.', 'MarkerSize', 15, 'color', [1, 0.5, 0]);
                hold on
                
                if y ~= 3
                    
                    plot(c, plotdata(8,:), '.', 'MarkerSize', 15, 'color', [0, 0.6, 0.6]);
                    hold on
                    
                end
                
                plot(c, plotdata(7,:), '.', 'MarkerSize', 15, 'color', [1, 0, 1]);
                hold on
                plot(c, plotdata(6,:), '.', 'MarkerSize', 15, 'color', [0.2, 1, 1]);
                hold on
                plot(c, plotdata(5,:), '.', 'MarkerSize', 15, 'color', [0, 0, 1]);
                hold on
                plot(c, plotdata(4,:), '.', 'MarkerSize', 15, 'color', [1, 0, 0]);
                hold on
                plot(c, plotdata(3,:), '.', 'MarkerSize', 15, 'color', [0, 1, 0]);
                hold on
                plot(c, plotdata(2,:), '.', 'MarkerSize', 15, 'color', [0.9, 0.9, 0]);
                hold on
                plot(c, plotdata(1,:), '.', 'MarkerSize', 15, 'color', [0, 0, 0]);
                hold on
                
                if strcmp(comp, 'task-difftask') && y == 1
                
                    title(['All Subjects Memory Task-Diff Task ' titlelabel ' Correlations'], 'fontsize',18)
                    
                elseif strcmp(comp, 'task-difftask') && y == 2
                
                    title(['All Subjects Mixed Task-Diff Task ' titlelabel ' Correlations'], 'fontsize',18)
                    
                elseif strcmp(comp, 'task-difftask') && y == 3
                
                    title(['All Subjects Motor Task-Diff Task ' titlelabel ' Correlations'], 'fontsize',18)
                    
                elseif strcmp(comp, 'matched-samedata-task-difftask') && y == 1
                    
                    title(['All Subjects Memory Task-Diff Task ' titlelabel ' Correlations Matched'], 'fontsize',18)
                    
                elseif strcmp(comp, 'matched-samedata-task-difftask') && y == 2
                    
                    title(['All Subjects Mixed Task-Diff Task ' titlelabel ' Correlations Matched'], 'fontsize',18)
                    
                elseif strcmp(comp, 'matched-samedata-task-difftask') && y == 3
                    
                    title(['All Subjects Motor Task-Diff Task ' titlelabel ' Correlations Matched'], 'fontsize',18)
                    
                end
                
                if fishertransform == 1
                    
                    ylabel('Fisher Transformed Correlations', 'fontsize',18)
                    ylim([.3 2]);
                    
                else
                    
                    ylabel('Correlations', 'fontsize',18)
                    ylim([.3 1]);
                    
                end
                
                h = findobj(gca,'Type','line');
                
                if fishertransform == 1
                    
                    legendpos = 'NorthEast';
                    
                else
                    
                    legendpos = 'SouthEast';
                    
                end
                
                if y == 3
                    
                    legend(h(1:8), 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC10', 'Location', legendpos)
                    
                else
                
                    legend(h(1:9), 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', legendpos)
                    
                end
                
                

                
                if saveindividplots == 1
            
                    if strcmp(comp, 'task-difftask') && y == 1
                        
                        if fishertransform == 1 && SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_memory_difftask_' titlelabel 'SNRMask_fishercorrelations.jpg'])
                        
                        elseif fishertransform == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_memory_difftask_' titlelabel '_fishercorrelations.jpg'])
                            
                        elseif SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_memory_difftask_' titlelabel 'SNRMask_correlations.jpg'])
                            
                        else
                
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_memory_difftask_' titlelabel '_correlations.jpg'])
                            
                        end
                    
                    elseif strcmp(comp, 'task-difftask') && y == 2
                        
                        if fishertransform == 1 && SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_mixed_difftask_' titlelabel 'SNRMask_fishercorrelations.jpg'])
                        
                        elseif fishertransform == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_mixed_difftask_' titlelabel '_fishercorrelations.jpg'])
                            
                        elseif SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_mixed_difftask_' titlelabel 'SNRMask_correlations.jpg'])
                            
                        else
                
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_mixed_difftask_' titlelabel '_correlations.jpg'])
                            
                        end
                    
                    elseif strcmp(comp, 'task-difftask') && y == 3
                        
                        if fishertransform == 1 && SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_motor_difftask_' titlelabel 'SNRMask_fishercorrelations.jpg'])
                        
                        elseif fishertransform == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_motor_difftask_' titlelabel '_fishercorrelations.jpg'])
                            
                        elseif SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_motor_difftask_' titlelabel 'SNRMask_correlations.jpg'])
                            
                        else
                
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_motor_difftask_' titlelabel '_correlations.jpg'])
                            
                        end
                    
                    elseif strcmp(comp, 'matched-samedata-task-difftask') && y == 1
                        
                        if fishertransform == 1 && SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_memory_difftask_matcheddata_' titlelabel 'SNRMask_fishercorrelations.jpg'])
                        
                        elseif fishertransform == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_memory_difftask_matcheddata_' titlelabel '_fishercorrelations.jpg'])
                            
                        elseif SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_memory_difftask_matcheddata_' titlelabel 'SNRMask_correlations.jpg'])
                            
                        else
                    
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_memory_difftask_matcheddata_' titlelabel '_correlations.jpg'])
                            
                        end
                    
                    elseif strcmp(comp, 'matched-samedata-task-difftask') && y == 2
                        
                        if fishertransform == 1 && SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_mixed_difftask_matcheddata_' titlelabel 'SNRMask_fishercorrelations.jpg'])
                        
                        elseif fishertransform == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_mixed_difftask_matcheddata_' titlelabel '_fishercorrelations.jpg'])
                            
                        elseif SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_mixed_difftask_matcheddata_' titlelabel 'SNRMask_correlations.jpg'])
                            
                        else
                    
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_mixed_difftask_matcheddata_' titlelabel '_correlations.jpg'])
                            
                        end
                    
                    elseif strcmp(comp, 'matched-samedata-task-difftask') && y == 3
                        
                        if fishertransform == 1 && SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_motor_difftask_matcheddata_' titlelabel 'SNRMask_fishercorrelations.jpg'])
                        
                        elseif fishertransform == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_motor_difftask_matcheddata_' titlelabel '_fishercorrelations.jpg'])
                            
                        elseif SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_motor_difftask_matcheddata_' titlelabel 'SNRMask_correlations.jpg'])
                            
                        else
                    
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_motor_difftask_matcheddata_' titlelabel '_correlations.jpg'])
                            
                        end
                    
                    end
                end
                
                if concatplots == 0
                
                    close gcf
                    
                end
                
                if concatplots == 1
                    
                    if fishertransform == 1 && SNRMasks == 1
                        
                        export_fig([outputdir '/PercentThreshold_AllPlots_matcheddata_SNRMask_fishercorrelations.pdf'], '-pdf', '-append')
                    
                    elseif fishertransform == 1
                        
                        export_fig([outputdir '/PercentThreshold_AllPlots_matcheddata_fishercorrelations.pdf'], '-pdf', '-append')
                        
                    elseif SNRMasks == 1
                        
                        export_fig([outputdir '/PercentThreshold_AllPlots_matcheddata_SNRMask_correlations.pdf'], '-pdf', '-append')
                        
                    else
                    
                        export_fig([outputdir '/PercentThreshold_AllPlots_matcheddata_correlations.pdf'], '-pdf', '-append')
                        
                    end
                    
                    close gcf 
                    
                end
                    
                end
            end
            
        elseif strcmp(comp, 'alltask-alltask') || strcmp(comp, 'matched-samedata-alltask-alltask')
            
            for t = 1:2
            
                % Plot task-task comparisons
                
                if t == 1
                    
                    if isempty(alltaskrestalltaskrestdicecoefs)
                    
                        plotdata = [alltaskalltaskpearsoncoefs alltaskrestpearsoncoefs restrestpearsoncoefs];
                        
                    else
                        
                        plotdata = [alltaskalltaskpearsoncoefs alltaskrestpearsoncoefs alltaskrestalltaskrestpearsoncoefs];
                        
                    end
                    
                else
                    
                    if isempty(alltaskrestalltaskrestdicecoefs)
                
                        plotdata = [alltaskalltaskdicecoefs(:,2) alltaskrestdicecoefs(:,2) restrestdicecoefs(:,2)];
                        
                    else
                        
                        plotdata = [alltaskalltaskdicecoefs(:,2) alltaskrestdicecoefs(:,2) alltaskrestalltaskrestdicecoefs(:,2)];
                        
                    end
                
                end
                
                SEplot = zeros(1,size(plotdata,2));
                
                if fishertransform == 1
                            
                	for w = 1:size(plotdata,2)
                            
                    	plotdata(:,w) = atanh(plotdata(:,w));
                                
                    end
                end
                
                meanplotdata = mean(plotdata,1);
                    
                for z = 1:size(plotdata,2)

                    SEplot(z) = std(plotdata(:,z))/sqrt(size(plotdata,1));
                        
                end
                
                if t == 1
                    
                    titlelabel = 'Pearson';
                    
                else
                    
                    titlelabel = 'Dice';
                    
                end
                
                if t == 2 && y == 1
                    
                    taskrestpval = memtaskpval;
                    
                elseif t == 2 && y == 2
                    
                    taskrestpval = mixedtaskpval;
                    
                elseif t == 2 && y == 3
                    
                    taskrestpval = motortaskpval;
                    
                end
            
                %barwitherr(SEplot, categorical({'AllTask-AllTask', 'AllTask-Rest', 'Rest-Rest'}), meanplotdata);
                bar(categorical({'Task-Task', 'Task-Rest', 'Rest-Rest'}), meanplotdata, 'FaceColor', 'w');
                hold on
                errorbar(categorical({'Task-Task', 'Task-Rest', 'Rest-Rest'}),meanplotdata,SEplot,'.');
                
                if t == 2
                
                    plot(xlim,[taskrestpval taskrestpval], 'LineWidth', 3, 'color', [0, 0, 0])
                    hold on
                    
                end
                
                plot(categorical({'Task-Task', 'Task-Rest', 'Rest-Rest'}), plotdata(9,:), '.', 'MarkerSize', 15, 'color', [1, 0.5, 0]);
                hold on
                
                if y ~= 3
                    
                    plot(categorical({'Task-Task', 'Task-Rest', 'Rest-Rest'}), plotdata(8,:), '.', 'MarkerSize', 15, 'color', [0, 0.6, 0.6]);
                    hold on
                    
                end
                
                plot(categorical({'Task-Task', 'Task-Rest', 'Rest-Rest'}), plotdata(7,:), '.', 'MarkerSize', 15, 'color', [1, 0, 1]);
                hold on
                plot(categorical({'Task-Task', 'Task-Rest', 'Rest-Rest'}), plotdata(6,:), '.', 'MarkerSize', 15, 'color', [0.2, 1, 1]);
                hold on
                plot(categorical({'Task-Task', 'Task-Rest', 'Rest-Rest'}), plotdata(5,:), '.', 'MarkerSize', 15, 'color', [0, 0, 1]);
                hold on
                plot(categorical({'Task-Task', 'Task-Rest', 'Rest-Rest'}), plotdata(4,:), '.', 'MarkerSize', 15, 'color', [1, 0, 0]);
                hold on
                plot(categorical({'Task-Task', 'Task-Rest', 'Rest-Rest'}), plotdata(3,:), '.', 'MarkerSize', 15, 'color', [0, 1, 0]);
                hold on
                plot(categorical({'Task-Task', 'Task-Rest', 'Rest-Rest'}), plotdata(2,:), '.', 'MarkerSize', 15, 'color', [0.9, 0.9, 0]);
                hold on
                plot(categorical({'Task-Task', 'Task-Rest', 'Rest-Rest'}), plotdata(1,:), '.', 'MarkerSize', 15, 'color', [0, 0, 0]);
                hold on
                
                if strcmp(comp, 'alltask-alltask')
                
                    title(['All Subjects All Tasks-All Tasks ' titlelabel ' Correlations'], 'fontsize',18)
                    
                elseif strcmp(comp, 'matched-samedata-alltask-alltask')
                    
                    title(['All Subjects All Tasks-All Tasks ' titlelabel ' Correlations Matched'], 'fontsize',18)
                    
                end
                
                if fishertransform == 1
                    
                    ylabel('Fisher Transformed Correlations', 'fontsize',18)
                    ylim([.3 2]);
                    
                else
                    
                    ylabel('Correlations', 'fontsize',18)
                    ylim([.3 1]);
                    
                end
                
                h = findobj(gca,'Type','line');
                
                if fishertransform == 1
                    
                    legendpos = 'NorthEast';
                    
                else
                    
                    legendpos = 'SouthEast';
                    
                end
                
                if y == 3
                    
                    legend(h(1:8), 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC10', 'Location', legendpos)
                    
                else
                
                    legend(h(1:9), 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', legendpos)
                    
                end
                
                
            
            
                if saveindividplots == 1
                    
                    if strcmp(comp, 'alltask-alltask')
                        
                        if fishertransform == 1 && SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_alltasks_' titlelabel 'SNRMask_fishercorrelations.jpg'])
                        
                        elseif fishertransform == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_alltasks_' titlelabel '_fishercorrelations.jpg'])
                            
                        elseif SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_alltasks_' titlelabel 'SNRMask_correlations.jpg'])
                            
                        else
                        
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_alltasks_' titlelabel '_correlations.jpg'])
                            
                        end
                        
                    else
                        
                        if fishertransform == 1 && SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_alltasks_matcheddata_' titlelabel 'SNRMask_fishercorrelations.jpg'])
                        
                        elseif fishertransform == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_alltasks_matcheddata_' titlelabel '_fishercorrelations.jpg'])
                            
                        elseif SNRMasks == 1
                            
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_alltasks_matcheddata_' titlelabel 'SNRMask_correlations.jpg'])
                            
                        else
                    
                            saveas(gcf,[outputdir '/AllSubjects_percentthreshold_alltasks_matcheddata_' titlelabel '_correlations.jpg'])
                            
                        end
                    end
                    
                end
                
                if concatplots == 0
                
                    close gcf
                    
                end
                
                if concatplots == 1
                    
                    if fishertransform == 1 && SNRMasks == 1
                        
                        export_fig([outputdir '/PercentThreshold_AllPlots_matcheddata_SNRMask_fishercorrelations.pdf'], '-pdf', '-append')
                    
                    elseif fishertransform == 1
                        
                        export_fig([outputdir '/PercentThreshold_AllPlots_matcheddata_fishercorrelations.pdf'], '-pdf', '-append')
                        
                    elseif SNRMasks == 1
                        
                        export_fig([outputdir '/PercentThreshold_AllPlots_matcheddata_SNRMask_correlations.pdf'], '-pdf', '-append')
                        
                    else
                    
                        export_fig([outputdir '/PercentThreshold_AllPlots_matcheddata_correlations.pdf'], '-pdf', '-append')
                        
                    end
                    
                    close gcf 
                    
                end
            end
        end
    end
end





if plottaskdiffs == 1

	x = [restrestdicecoefs alltaskalltaskdicecoefs alltaskrestdicecoefs tasktaskdicecoefs(:,1) taskrestdicecoefs(:,1) tasktaskdicecoefs(:,2) taskrestdicecoefs(:,2) tasktaskdicecoefs(:,3) taskrestdicecoefs(:,3)];
	group = [repmat(1,1,length(restrestdicecoefs)) repmat(2,1,length(alltaskalltaskdicecoefs)) repmat(3,1,length(alltaskrestdicecoefs)) repmat(4,1,length(tasktaskdicecoefs(:,1))) repmat(5,1,length(taskrestdicecoefs(:,1))) repmat(6,1,length(tasktaskdicecoefs(:,2))) repmat(7,1,length(taskrestdicecoefs(:,2))) repmat(8,1,length(tasktaskdicecoefs(:,3))) repmat(9,1,length(taskrestdicecoefs(:,3)))];
	positions = [1 2 2.33 3 3.33 4 4.33 5 5.33];
    boxplot(x,group, 'BoxStyle', 'outline', 'positions', positions, 'Symbol', '');
	hold on
  	scatter(positions, [restrestdicecoefs(9,:) alltaskalltaskdicecoefs(9,:) alltaskrestdicecoefs(9,:) tasktaskdicecoefs(9,1) taskrestdicecoefs(9,1) tasktaskdicecoefs(9,2) taskrestdicecoefs(9,2) tasktaskdicecoefs(9,3) taskrestdicecoefs(9,3)], 'filled', 'MarkerFaceColor', [1, 0.5, 0]);
	hold on
  	scatter(positions, [restrestdicecoefs(8,:) alltaskalltaskdicecoefs(8,:) alltaskrestdicecoefs(8,:) tasktaskdicecoefs(8,1) taskrestdicecoefs(8,1) tasktaskdicecoefs(8,2) taskrestdicecoefs(8,2) tasktaskdicecoefs(8,3) taskrestdicecoefs(8,3)], 'filled', 'MarkerFaceColor', [0, 0.6, 0.6]);
   	hold on    
  	scatter(positions, [restrestdicecoefs(7,:) alltaskalltaskdicecoefs(7,:) alltaskrestdicecoefs(7,:) tasktaskdicecoefs(7,1) taskrestdicecoefs(7,1) tasktaskdicecoefs(7,2) taskrestdicecoefs(7,2) tasktaskdicecoefs(7,3) taskrestdicecoefs(7,3)], 'filled', 'MarkerFaceColor', [1, 0, 1]);
 	hold on
   	scatter(positions, [restrestdicecoefs(6,:) alltaskalltaskdicecoefs(6,:) alltaskrestdicecoefs(6,:) tasktaskdicecoefs(6,1) taskrestdicecoefs(6,1) tasktaskdicecoefs(6,2) taskrestdicecoefs(6,2) tasktaskdicecoefs(6,3) taskrestdicecoefs(6,3)], 'filled', 'MarkerFaceColor', [0.2, 1, 1]);
  	hold on
  	scatter(positions, [restrestdicecoefs(5,:) alltaskalltaskdicecoefs(5,:) alltaskrestdicecoefs(5,:) tasktaskdicecoefs(5,1) taskrestdicecoefs(5,1) tasktaskdicecoefs(5,2) taskrestdicecoefs(5,2) tasktaskdicecoefs(5,3) taskrestdicecoefs(5,3)], 'filled', 'MarkerFaceColor', [0, 0, 1]);
   	hold on
  	scatter(positions, [restrestdicecoefs(4,:) alltaskalltaskdicecoefs(4,:) alltaskrestdicecoefs(4,:) tasktaskdicecoefs(4,1) taskrestdicecoefs(4,1) tasktaskdicecoefs(4,2) taskrestdicecoefs(4,2) tasktaskdicecoefs(4,3) taskrestdicecoefs(4,3)], 'filled', 'MarkerFaceColor', [1, 0, 0]);
  	hold on
 	scatter(positions, [restrestdicecoefs(3,:) alltaskalltaskdicecoefs(3,:) alltaskrestdicecoefs(3,:) tasktaskdicecoefs(3,1) taskrestdicecoefs(3,1) tasktaskdicecoefs(3,2) taskrestdicecoefs(3,2) tasktaskdicecoefs(3,3) taskrestdicecoefs(3,3)], 'filled', 'MarkerFaceColor', [0, 1, 0]);
 	hold on
  	scatter(positions, [restrestdicecoefs(2,:) alltaskalltaskdicecoefs(2,:) alltaskrestdicecoefs(2,:) tasktaskdicecoefs(2,1) taskrestdicecoefs(2,1) tasktaskdicecoefs(2,2) taskrestdicecoefs(2,2) tasktaskdicecoefs(2,3) taskrestdicecoefs(2,3)], 'filled', 'MarkerFaceColor', [0.9, 0.9, 0]);
  	hold on
  	scatter(positions, [restrestdicecoefs(1,:) alltaskalltaskdicecoefs(1,:) alltaskrestdicecoefs(1,:) tasktaskdicecoefs(1,1) taskrestdicecoefs(1,1) tasktaskdicecoefs(1,2) taskrestdicecoefs(1,2) tasktaskdicecoefs(1,3) taskrestdicecoefs(1,3)], 'filled', 'MarkerFaceColor', [0, 0, 0]);
  	hold on
    
   	set(gca,'xtick',[positions(1) mean(positions(2:3)) mean(positions(4:5)) mean(positions(6:7)) mean(positions(8:9))])
 	set(gca,'xticklabel',{'Rest', 'All Tasks', 'Memory Task', 'Mixed Task', 'Motor Task'}, 'FontSize',20)

  	color = ['b', 'g', 'b', 'g', 'b', 'g', 'b', 'g', 'g'];
  	h = findobj(gca,'Tag','Box');
        
  	for j=1:length(h)
            
        patch(get(h(j),'XData'),get(h(j),'YData'),color(j), 'FaceColor', 'none', 'EdgeColor', color(j));
            
    end
    
 	c = get(gca, 'Children');
            
	m = findobj(gca,'Type','scatter');
            
	hleg1 = legend([m(1:9); c(7:8)], 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Between State Comparison', 'Within State Comparison', 'Location', 'NorthEast');
	hleg1.FontSize = 18;
	ylim([0 1]);
    ylabel('Dice Correlations')
    
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    
 	filename = '/AllSubjects_Boxplot_AcrossTaskComparisons_SNRExclude.jpg';
      
 	saveas(gcf,[outputdir filename])
    
    close gcf
    
end
    