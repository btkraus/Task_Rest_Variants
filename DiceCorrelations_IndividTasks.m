
clear all

pearsoncoefs = [];
dicecoefs = [];
taskoutput = {};
threshoutput = [];
outputdir = '/Users/briankraus/Desktop/Correct_Variant_Maps/TaskVariantPlots';

% tasktaskpearsoncoefs = [];
% taskdifftaskpearsoncoefs = [];
% taskrestpearsoncoefs = [];
% restrestpearsoncoefs = [];
% taskresttaskrestpearsoncoefs = [];
% alltaskrestpearsoncoefs = [];
% alltaskalltaskpearsoncoefs = [];
% alltaskrestalltaskrestpearsoncoefs = [];

memmemdicecoefs = [];
mixedmixeddicecoefs = [];
motormotordicecoefs = [];
memmemrestdicecoefs = [];
mixedmixedrestdicecoefs = [];
motormotorrestdicecoefs = [];
memrestmemrestdicecoefs = [];
mixedrestmixedrestdicecoefs = [];
motorrestmotorrestdicecoefs = [];
alltaskalltaskdicecoefs = [];
alltaskalltaskrestdicecoefs = [];
alltaskrestalltaskrestdicecoefs = [];

randomizevals = 1;  %% Toggles whether to calculate a randomized value for each task
BoxPlot = 0;  %% Toggles whether to plot individual tasks versus all tasks
Subplots = 1;  %% Toggles whether to plot each task as a subplot in a larger figure
AddLegend = 1;  %% Toggles whether to plot legend in figure
threshold = 2.5;
%thresholds = [5 10 15 20 25 30 35];  %% Thresholds to run dice correlations (percent lowest correlations)
%thresholds = [.15 .2 .25 .3 .35 .4 .45]; %% Thresholds to run dice correlations (absolute correlation threshold)
SNRMasks = 1;  %% Toggles whether to exclude areas of low SNR via SNR masks
SizeExclude = 1;  %% Toggles whether to exclude variants smaller than criterion
minsize = 15;  %% Minimum variant size for size exclusion
ExcludeMSC09 = 1;  %% Toggles where to exclude MSC09
%FinalPlot = 1;  %% Toggles whether to plot the final plot for the paper

comparison = {'matched-samedata-task-task', 'matched-samedata-task-rest', 'matched-samedata-taskrest-taskrest', 'matched-samedata-alltask-alltask', 'matched-samedata-alltask-rest', 'matched-samedata-alltaskrest-alltaskrest'};



% Load all files, apply size exclusion if necessary

memfilestaskeven = [];
memfilestaskodd = [];
mixedfilestaskeven = [];
mixedfilestaskodd = [];
motorfilestaskeven = [];
motorfilestaskodd = [];
memfilesresteven = [];
memfilesrestodd = [];
mixedfilesresteven = [];
mixedfilesrestodd = [];
motorfilesresteven = [];
motorfilesrestodd = [];
alltaskfilestaskeven = [];
alltaskfilestaskodd = [];
alltaskfilesresteven = [];
alltaskfilesrestodd = [];

for b = 1:4  %% Number of iterations for each file type
    
    if ExcludeMSC09 == 1
        
        if b == 1   %% Task Files

            [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_task_samedata_consec_even_varmaps_NoMSC09.txt','%s%s%s');

            [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_task_samedata_consec_odd_varmaps_NoMSC09.txt','%s%s%s');

            [task_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_task_samedata_consec_even_varmaps_NoMSC09_SizeOnly.txt','%s%s%s');

            [rest_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_task_samedata_consec_odd_varmaps_NoMSC09_SizeOnly.txt','%s%s%s');
            
        elseif b == 2   %% Task-matched rest files
            
            [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_task_samedata_consec_even_varmaps_NoMSC09.txt','%s%s%s');

            [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_task_samedata_consec_odd_varmaps_NoMSC09.txt','%s%s%s');

            [task_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_task_samedata_consec_even_varmaps_NoMSC09_SizeOnly.txt','%s%s%s');

            [rest_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_task_samedata_consec_odd_varmaps_NoMSC09_SizeOnly.txt','%s%s%s');
            
      	elseif b == 3   %% All Tasks Concatenated Files
            
            [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_samedata_consec_even_varmaps_NoMSC09.txt','%s%s%s');

            [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_samedata_consec_odd_varmaps_NoMSC09.txt','%s%s%s');

            [task_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_samedata_consec_even_varmaps_NoMSC09_SizeOnly.txt','%s%s%s');

            [rest_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_samedata_consec_odd_varmaps_NoMSC09_SizeOnly.txt','%s%s%s');
            
      	elseif b == 4 %% All Tasks Concatenated-matched rest Files
            
            [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_alltask_samedata_consec_even_varmaps_NoMSC09.txt','%s%s%s');

            [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_alltask_samedata_consec_odd_varmaps_NoMSC09.txt','%s%s%s');

            [task_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_alltask_samedata_consec_even_varmaps_NoMSC09_SizeOnly.txt','%s%s%s');

            [rest_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_alltask_samedata_consec_odd_varmaps_NoMSC09_SizeOnly.txt','%s%s%s');
            
        end
        
    else
    
        if b == 1   %% Task Files

            [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_task_samedata_consec_even_varmaps.txt','%s%s%s');

            [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_task_samedata_consec_odd_varmaps.txt','%s%s%s');

            [task_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_task_samedata_consec_even_varmaps_SizeOnly.txt','%s%s%s');

            [rest_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_task_samedata_consec_odd_varmaps_SizeOnly.txt','%s%s%s');

        elseif b == 2   %% Task-matched rest files

            [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_task_samedata_consec_even_varmaps.txt','%s%s%s');

            [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_task_samedata_consec_odd_varmaps.txt','%s%s%s');

            [task_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_task_samedata_consec_even_varmaps_SizeOnly.txt','%s%s%s');

            [rest_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_task_samedata_consec_odd_varmaps_SizeOnly.txt','%s%s%s');

        elseif b == 3   %% All Tasks Concatenated Files

            [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_samedata_consec_even_varmaps.txt','%s%s%s');

            [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_samedata_consec_odd_varmaps.txt','%s%s%s');

            [task_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_samedata_consec_even_varmaps_SizeOnly.txt','%s%s%s');

            [rest_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_samedata_consec_odd_varmaps_SizeOnly.txt','%s%s%s');

        elseif b == 4 %% All Tasks Concatenated-matched rest Files

            [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_alltask_samedata_consec_even_varmaps.txt','%s%s%s');

            [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_alltask_samedata_consec_odd_varmaps.txt','%s%s%s');

            [task_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_alltask_samedata_consec_even_varmaps_SizeOnly.txt','%s%s%s');

            [rest_masks, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_alltask_samedata_consec_odd_varmaps_SizeOnly.txt','%s%s%s');

        end
    end

    for a = 1:size(task_files)

        sub1 = subjects1{a};
        sub2 = subjects2{a};
        task_file = task_files{a};
        rest_file = rest_files{a};
        task1 = tasks1{a};
        task2 = tasks2{a};

        cifti_task = ft_read_cifti_mod(task_file);
        cifti_rest = ft_read_cifti_mod(rest_file);

        if SNRMasks == 1

            SNRmask = ft_read_cifti_mod(['/Users/briankraus/Desktop/MSC_SNR_Maps/CIFTI_Files/MSC_Template/Individual_Session_Map/' sub1 '/' sub1 '__SNRMap__AllDataConcatenated.dscalar.nii']); 

            SNRmask.data = SNRmask.data(1:59412,:);
            LowSNR = find(SNRmask.data < 750);

            cifti_task.data(LowSNR,:) = NaN;
            cifti_rest.data(LowSNR,:) = NaN;

        end
        
    	cifti_task_threshold = find(cifti_task.data < prctile(cifti_task.data,threshold));
      	cifti_rest_threshold = find(cifti_rest.data < prctile(cifti_rest.data,threshold));
        
        cifti_task.data(:,:) = 0;
        cifti_rest.data(:,:) = 0;
        cifti_task.data(cifti_task_threshold,:) = 1;
        cifti_rest.data(cifti_rest_threshold,:) = 1;

        if SizeExclude == 1

            task_mask = task_masks{a};
            rest_mask = rest_masks{a};

            cifti_task_mask = ft_read_cifti_mod(task_mask);
            cifti_rest_mask = ft_read_cifti_mod(rest_mask);

            allvars_rest = unique(cifti_rest_mask.data);
            allvars_rest(1) = [];
            allvars_task = unique(cifti_task_mask.data);
            allvars_task(1) = [];

            removevars_rest = [];
            removevars_task = [];

            for z = 1:length(allvars_rest)

                if length(find(cifti_rest.data == allvars_rest(z))) < minsize

                    removevars_rest = [removevars_rest allvars_rest(z)];

                end
            end

            for h = 1:length(cifti_rest.data)

                if ismember(cifti_rest.data(h),removevars_rest)

                    cifti_rest.data(h) = 0;

                end
            end

            for z = 1:length(allvars_task)

                if length(find(cifti_task.data == allvars_task(z))) < minsize

                    removevars_task = [removevars_task allvars_task(z)];

                end
            end

            for h = 1:length(cifti_task.data)

                if ismember(cifti_task.data(h),removevars_task)

                    cifti_task.data(h) = 0;

                end
            end
        end
    
        if b == 1   %% Task Files
            
            if strcmp(task1,'mem')
                
                memfilestaskeven = [memfilestaskeven; cifti_task.data'];
                memfilestaskodd = [memfilestaskodd; cifti_rest.data'];
                
            elseif strcmp(task1,'mixed')
                
                mixedfilestaskeven = [mixedfilestaskeven; cifti_task.data'];
                mixedfilestaskodd = [mixedfilestaskodd; cifti_rest.data'];
                
            elseif strcmp(task1,'motor')
                
                motorfilestaskeven = [motorfilestaskeven; cifti_task.data'];
                motorfilestaskodd = [motorfilestaskodd; cifti_rest.data'];
                
            end
            
        elseif b == 2   %% Task-matched rest files
                
            if strcmp(task1,'memrest')
                
                memfilesresteven = [memfilesresteven; cifti_task.data'];
                memfilesrestodd = [memfilesrestodd; cifti_rest.data'];
                
            elseif strcmp(task1,'mixedrest')
                
                mixedfilesresteven = [mixedfilesresteven; cifti_task.data'];
                mixedfilesrestodd = [mixedfilesrestodd; cifti_rest.data'];
                
            elseif strcmp(task1,'motorrest')
                
                motorfilesresteven = [motorfilesresteven; cifti_task.data'];
                motorfilesrestodd = [motorfilesrestodd; cifti_rest.data'];
                
            end
            
        elseif b == 3   %% All Tasks Concatenated Files
            
            alltaskfilestaskeven = [alltaskfilestaskeven; cifti_task.data'];
            alltaskfilestaskodd = [alltaskfilestaskodd; cifti_rest.data'];
            
        elseif b == 4   %% All Tasks Concatenated-matched rest Files
            
            alltaskfilesresteven = [alltaskfilesresteven; cifti_task.data'];
            alltaskfilesrestodd = [alltaskfilesrestodd; cifti_rest.data'];
            
        end
    end
end


for g = 1:numel(comparison)
    
    comp = comparison(g);
    
    if strcmp(comp, 'matched-samedata-task-task')
         
     	for c = 1:size(memfilestaskeven,1)
            
            % Mem Task
            
          	dcorrdatamem = [];
        
          	for q = 1:size(memfilestaskeven,2)
            
             	if memfilestaskeven(c,q) == 1 && memfilestaskodd(c,q) == 1
                
                 	dcorrdatamem = [dcorrdatamem;1 1];
                
             	elseif memfilestaskeven(c,q) == 1
                
                  	dcorrdatamem = [dcorrdatamem;1 0];
                
              	elseif memfilestaskodd(c,q) == 1
                
                  	dcorrdatamem = [dcorrdatamem;0 1];
                
                end
            
            end
        
        	dcmem = dice_coefficient_mod(dcorrdatamem(:,1),dcorrdatamem(:,2));
            
            memmemdicecoefs = [memmemdicecoefs; dcmem];
            
            % Mixed Task
            
          	dcorrdatamixed = [];
        
          	for q = 1:size(mixedfilestaskeven,2)
            
             	if mixedfilestaskeven(c,q) == 1 && mixedfilestaskodd(c,q) == 1
                
                 	dcorrdatamixed = [dcorrdatamixed;1 1];
                
             	elseif mixedfilestaskeven(c,q) == 1
                
                  	dcorrdatamixed = [dcorrdatamixed;1 0];
                
              	elseif mixedfilestaskodd(c,q) == 1
                
                  	dcorrdatamixed = [dcorrdatamixed;0 1];
                
                end
            
            end
        
        	dcmixed = dice_coefficient_mod(dcorrdatamixed(:,1),dcorrdatamixed(:,2));
            
            mixedmixeddicecoefs = [mixedmixeddicecoefs; dcmixed];
            
            % Motor Task
            
          	dcorrdatamotor = [];
        
          	for q = 1:size(motorfilestaskeven,2)
            
             	if motorfilestaskeven(c,q) == 1 && motorfilestaskodd(c,q) == 1
                
                 	dcorrdatamotor = [dcorrdatamotor;1 1];
                
             	elseif motorfilestaskeven(c,q) == 1
                
                  	dcorrdatamotor = [dcorrdatamotor;1 0];
                
              	elseif motorfilestaskodd(c,q) == 1
                
                  	dcorrdatamotor = [dcorrdatamotor;0 1];
                
                end
            
            end
        
        	dcmotor = dice_coefficient_mod(dcorrdatamotor(:,1),dcorrdatamotor(:,2));
            
            motormotordicecoefs = [motormotordicecoefs; dcmotor];
            
        end
        
    elseif strcmp(comp, 'matched-samedata-task-rest')
        
     	for c = 1:size(memfilesresteven,1)
            
            % Mem Task
            
          	dcorrdatamemresteven = [];
        
          	for q = 1:size(memfilesresteven,2)
            
             	if memfilesresteven(c,q) == 1 && memfilestaskodd(c,q) == 1
                
                 	dcorrdatamemresteven = [dcorrdatamemresteven;1 1];
                
             	elseif memfilesresteven(c,q) == 1
                
                  	dcorrdatamemresteven = [dcorrdatamemresteven;1 0];
                
              	elseif memfilestaskodd(c,q) == 1
                
                  	dcorrdatamemresteven = [dcorrdatamemresteven;0 1];
                
                end
            
            end
        
        	dcmemresteven = dice_coefficient_mod(dcorrdatamemresteven(:,1),dcorrdatamemresteven(:,2));
            
            
          	dcorrdatamemrestodd = [];
        
          	for q = 1:size(memfilesrestodd,2)
            
             	if memfilesrestodd(c,q) == 1 && memfilestaskeven(c,q) == 1
                
                 	dcorrdatamemrestodd = [dcorrdatamemrestodd;1 1];
                
             	elseif memfilesrestodd(c,q) == 1
                
                  	dcorrdatamemrestodd = [dcorrdatamemrestodd;1 0];
                
              	elseif memfilestaskeven(c,q) == 1
                
                  	dcorrdatamemrestodd = [dcorrdatamemrestodd;0 1];
                
                end
            
            end
        
        	dcmemrestodd = dice_coefficient_mod(dcorrdatamemrestodd(:,1),dcorrdatamemrestodd(:,2));
            
            memmemrestdicecoefs = [memmemrestdicecoefs; dcmemresteven dcmemrestodd];

            % Mixed Task
            
          	dcorrdatamixedresteven = [];
        
          	for q = 1:size(mixedfilesresteven,2)
            
             	if mixedfilesresteven(c,q) == 1 && mixedfilestaskodd(c,q) == 1
                
                 	dcorrdatamixedresteven = [dcorrdatamixedresteven;1 1];
                
             	elseif mixedfilesresteven(c,q) == 1
                
                  	dcorrdatamixedresteven = [dcorrdatamixedresteven;1 0];
                
              	elseif mixedfilestaskodd(c,q) == 1
                
                  	dcorrdatamixedresteven = [dcorrdatamixedresteven;0 1];
                
                end
            
            end
        
        	dcmixedresteven = dice_coefficient_mod(dcorrdatamixedresteven(:,1),dcorrdatamixedresteven(:,2));
            
            
          	dcorrdatamixedrestodd = [];
        
          	for q = 1:size(mixedfilesrestodd,2)
            
             	if mixedfilesrestodd(c,q) == 1 && mixedfilestaskeven(c,q) == 1
                
                 	dcorrdatamixedrestodd = [dcorrdatamixedrestodd;1 1];
                
             	elseif mixedfilesrestodd(c,q) == 1
                
                  	dcorrdatamixedrestodd = [dcorrdatamixedrestodd;1 0];
                
              	elseif mixedfilestaskeven(c,q) == 1
                
                  	dcorrdatamixedrestodd = [dcorrdatamixedrestodd;0 1];
                
                end
            
            end
        
        	dcmixedrestodd = dice_coefficient_mod(dcorrdatamixedrestodd(:,1),dcorrdatamixedrestodd(:,2));
            
            mixedmixedrestdicecoefs = [mixedmixedrestdicecoefs; dcmixedresteven dcmixedrestodd];
            
            % Motor Task
            
          	dcorrdatamotorresteven = [];
        
          	for q = 1:size(motorfilesresteven,2)
            
             	if motorfilesresteven(c,q) == 1 && motorfilestaskodd(c,q) == 1
                
                 	dcorrdatamotorresteven = [dcorrdatamotorresteven;1 1];
                
             	elseif motorfilesresteven(c,q) == 1
                
                  	dcorrdatamotorresteven = [dcorrdatamotorresteven;1 0];
                
              	elseif motorfilestaskodd(c,q) == 1
                
                  	dcorrdatamotorresteven = [dcorrdatamotorresteven;0 1];
                
                end
            
            end
        
        	dcmotorresteven = dice_coefficient_mod(dcorrdatamotorresteven(:,1),dcorrdatamotorresteven(:,2));
            
            
          	dcorrdatamotorrestodd = [];
        
          	for q = 1:size(motorfilesrestodd,2)
            
             	if motorfilesrestodd(c,q) == 1 && motorfilestaskeven(c,q) == 1
                
                 	dcorrdatamotorrestodd = [dcorrdatamotorrestodd;1 1];
                
             	elseif motorfilesrestodd(c,q) == 1
                
                  	dcorrdatamotorrestodd = [dcorrdatamotorrestodd;1 0];
                
              	elseif motorfilestaskeven(c,q) == 1
                
                  	dcorrdatamotorrestodd = [dcorrdatamotorrestodd;0 1];
                
                end
            
            end
        
        	dcmotorrestodd = dice_coefficient_mod(dcorrdatamotorrestodd(:,1),dcorrdatamotorrestodd(:,2));
            
            motormotorrestdicecoefs = [motormotorrestdicecoefs; dcmotorresteven dcmotorrestodd];
            
        end

    elseif strcmp(comp, 'matched-samedata-taskrest-taskrest')
        
     	for c = 1:size(memfilesresteven,1)
            
            % Mem Task
            
          	dcorrdatamemrest = [];
        
          	for q = 1:size(memfilesresteven,2)
            
             	if memfilesresteven(c,q) == 1 && memfilesrestodd(c,q) == 1
                
                 	dcorrdatamemrest = [dcorrdatamemrest;1 1];
                
             	elseif memfilesresteven(c,q) == 1
                
                  	dcorrdatamemrest = [dcorrdatamemrest;1 0];
                
              	elseif memfilesrestodd(c,q) == 1
                
                  	dcorrdatamemrest = [dcorrdatamemrest;0 1];
                
                end
            
            end
        
        	dcmemrest = dice_coefficient_mod(dcorrdatamemrest(:,1),dcorrdatamemrest(:,2));
            
            memrestmemrestdicecoefs = [memrestmemrestdicecoefs; dcmemrest];
            
            % Mixed Task
            
          	dcorrdatamixedrest = [];
        
          	for q = 1:size(mixedfilesresteven,2)
            
             	if mixedfilesresteven(c,q) == 1 && mixedfilesrestodd(c,q) == 1
                
                 	dcorrdatamixedrest = [dcorrdatamixedrest;1 1];
                
             	elseif mixedfilesresteven(c,q) == 1
                
                  	dcorrdatamixedrest = [dcorrdatamixedrest;1 0];
                
              	elseif mixedfilesrestodd(c,q) == 1
                
                  	dcorrdatamixedrest = [dcorrdatamixedrest;0 1];
                
                end
            
            end
        
        	dcmixedrest = dice_coefficient_mod(dcorrdatamixedrest(:,1),dcorrdatamixedrest(:,2));
            
            mixedrestmixedrestdicecoefs = [mixedrestmixedrestdicecoefs; dcmixedrest];
            
            % Motor Task
            
          	dcorrdatamotorrest = [];
        
          	for q = 1:size(motorfilesresteven,2)
            
             	if motorfilesresteven(c,q) == 1 && motorfilesrestodd(c,q) == 1
                
                 	dcorrdatamotorrest = [dcorrdatamotorrest;1 1];
                
             	elseif motorfilesresteven(c,q) == 1
                
                  	dcorrdatamotorrest = [dcorrdatamotorrest;1 0];
                
              	elseif motorfilesrestodd(c,q) == 1
                
                  	dcorrdatamotorrest = [dcorrdatamotorrest;0 1];
                
                end
            
            end
        
        	dcmotorrest = dice_coefficient_mod(dcorrdatamotorrest(:,1),dcorrdatamotorrest(:,2));
            
            motorrestmotorrestdicecoefs = [motorrestmotorrestdicecoefs; dcmotorrest];
            
        end
        
    elseif strcmp(comp, 'matched-samedata-alltask-alltask')
        
     	for c = 1:size(alltaskfilestaskeven,1)
            
            % All Tasks Concatenated
            
          	dcorrdataalltasks = [];
        
          	for q = 1:size(alltaskfilestaskeven,2)
            
             	if alltaskfilestaskeven(c,q) == 1 && alltaskfilestaskodd(c,q) == 1
                
                 	dcorrdataalltasks = [dcorrdataalltasks;1 1];
                
             	elseif alltaskfilestaskeven(c,q) == 1
                
                  	dcorrdataalltasks = [dcorrdataalltasks;1 0];
                
              	elseif alltaskfilestaskodd(c,q) == 1
                
                  	dcorrdataalltasks = [dcorrdataalltasks;0 1];
                
                end
            
            end
        
        	dcalltasks = dice_coefficient_mod(dcorrdataalltasks(:,1),dcorrdataalltasks(:,2));
            
            alltaskalltaskdicecoefs = [alltaskalltaskdicecoefs; dcalltasks];
            
        end
        
    elseif strcmp(comp, 'matched-samedata-alltask-rest')
        
     	for c = 1:size(alltaskfilesresteven,1)
            
            % All Tasks Concatenated
            
          	dcorrdataalltasksresteven = [];
        
          	for q = 1:size(alltaskfilesresteven,2)
            
             	if alltaskfilesresteven(c,q) == 1 && alltaskfilestaskodd(c,q) == 1
                
                 	dcorrdataalltasksresteven = [dcorrdataalltasksresteven;1 1];
                
             	elseif alltaskfilesresteven(c,q) == 1
                
                  	dcorrdataalltasksresteven = [dcorrdataalltasksresteven;1 0];
                
              	elseif alltaskfilestaskodd(c,q) == 1
                
                  	dcorrdataalltasksresteven = [dcorrdataalltasksresteven;0 1];
                
                end
            
            end
        
        	dcalltasksresteven = dice_coefficient_mod(dcorrdataalltasksresteven(:,1),dcorrdataalltasksresteven(:,2));
            
            
          	dcorrdataalltasksrestodd = [];
        
          	for q = 1:size(alltaskfilesrestodd,2)
            
             	if alltaskfilesrestodd(c,q) == 1 && alltaskfilestaskeven(c,q) == 1
                
                 	dcorrdataalltasksrestodd = [dcorrdataalltasksrestodd;1 1];
                
             	elseif alltaskfilesrestodd(c,q) == 1
                
                  	dcorrdataalltasksrestodd = [dcorrdataalltasksrestodd;1 0];
                
              	elseif alltaskfilestaskeven(c,q) == 1
                
                  	dcorrdataalltasksrestodd = [dcorrdataalltasksrestodd;0 1];
                
                end
            
            end
        
        	dcalltasksrestodd = dice_coefficient_mod(dcorrdataalltasksrestodd(:,1),dcorrdataalltasksrestodd(:,2));
            
            alltaskalltaskrestdicecoefs = [alltaskalltaskrestdicecoefs; dcalltasksresteven dcalltasksrestodd];
            
        end
        
    elseif strcmp(comp, 'matched-samedata-alltaskrest-alltaskrest')
        
     	for c = 1:size(alltaskfilesresteven,1)
            
            % All Tasks Concatenated
            
          	dcorrdataalltasksrest = [];
        
          	for q = 1:size(alltaskfilesresteven,2)
            
             	if alltaskfilesresteven(c,q) == 1 && alltaskfilesrestodd(c,q) == 1
                
                 	dcorrdataalltasksrest = [dcorrdataalltasksrest;1 1];
                
             	elseif alltaskfilesresteven(c,q) == 1
                
                  	dcorrdataalltasksrest = [dcorrdataalltasksrest;1 0];
                
              	elseif alltaskfilesrestodd(c,q) == 1
                
                  	dcorrdataalltasksrest = [dcorrdataalltasksrest;0 1];
                
                end
            
            end
        
        	dcalltasksrest = dice_coefficient_mod(dcorrdataalltasksrest(:,1),dcorrdataalltasksrest(:,2));
            
            alltaskrestalltaskrestdicecoefs = [alltaskrestalltaskrestdicecoefs; dcalltasksrest];
            
        end
        
    end
end



if randomizevals == 1
    
  	nsubs = size(memfilestaskeven,1);
            
  	DiceCorrssimmem = [];
  	DiceCorrSubsmem = zeros(2*(nsubs-1),nsubs);
  	DiceCorrssimmixed = [];
  	DiceCorrSubsmixed = zeros(2*(nsubs-1),nsubs);
   	DiceCorrssimmotor = [];
  	DiceCorrSubsmotor = zeros(2*(nsubs-1),nsubs);
   	DiceCorrssimalltasks = [];
  	DiceCorrSubsalltasks = zeros(2*(nsubs-1),nsubs);

 	for l = 1:2     %% Do both versions of split-half combinations
                
      	for m = 1:nsubs
                    
          	loopvals = [1:nsubs];
         	loopvals(m) = [];
                    
         	Count = 0;
                    
          	for n = loopvals
                        
            	Count = Count+1;

              	if l == 1
                            
                  	memdatrest = memfilesrestodd(m,:);
                 	memdat = memfilestaskeven(n,:);
                    
                  	mixeddatrest = mixedfilesrestodd(m,:);
                 	mixeddat = mixedfilestaskeven(n,:);
                    
                  	motordatrest = motorfilesrestodd(m,:);
                 	motordat = motorfilestaskeven(n,:);
                    
                  	alltasksdatrest = alltaskfilesrestodd(m,:);
                 	alltasksdat = alltaskfilestaskeven(n,:);
                            
                else
                            
                  	memdatrest = memfilesresteven(m,:);
                 	memdat = memfilestaskodd(n,:);
                    
                  	mixeddatrest = mixedfilesresteven(m,:);
                 	mixeddat = mixedfilestaskodd(n,:);
                    
                  	motordatrest = motorfilesresteven(m,:);
                 	motordat = motorfilestaskodd(n,:);
                    
                  	alltasksdatrest = alltaskfilesresteven(m,:);
                 	alltasksdat = alltaskfilestaskodd(n,:);
                            
                end
                            
              	dcorrdatamem = [];
        
             	for q = 1:length(memdatrest)
            
                   	if memdatrest(:,q) > 0 && memdat(:,q) > 0
                
                     	dcorrdatamem = [dcorrdatamem;1 1];
                
                  	elseif memdatrest(:,q) > 0
                
                     	dcorrdatamem = [dcorrdatamem;1 0];
                
                	elseif memdat(:,q) > 0
                
                      	dcorrdatamem = [dcorrdatamem;0 1];
                
                    end
            
                end
                    
             	if isempty(dcorrdatamem)
                        
                 	dcmem = 0;
                        
                else
        
                  	dcmem = dice_coefficient_mod(dcorrdatamem(:,1),dcorrdatamem(:,2));
                        
                end
                
              	dcorrdatamixed = [];
        
             	for q = 1:length(mixeddatrest)
            
                   	if mixeddatrest(:,q) > 0 && mixeddat(:,q) > 0
                
                     	dcorrdatamixed = [dcorrdatamixed;1 1];
                
                  	elseif mixeddatrest(:,q) > 0
                
                     	dcorrdatamixed = [dcorrdatamixed;1 0];
                
                	elseif mixeddat(:,q) > 0
                
                      	dcorrdatamixed = [dcorrdatamixed;0 1];
                
                    end
            
                end
                    
             	if isempty(dcorrdatamixed)
                        
                 	dcmixed = 0;
                        
                else
        
                  	dcmixed = dice_coefficient_mod(dcorrdatamixed(:,1),dcorrdatamixed(:,2));
                        
                end
                
              	dcorrdatamotor = [];
        
             	for q = 1:length(motordatrest)
            
                   	if motordatrest(:,q) > 0 && motordat(:,q) > 0
                
                     	dcorrdatamotor = [dcorrdatamotor;1 1];
                
                  	elseif motordatrest(:,q) > 0
                
                     	dcorrdatamotor = [dcorrdatamotor;1 0];
                
                	elseif motordat(:,q) > 0
                
                      	dcorrdatamotor = [dcorrdatamotor;0 1];
                
                    end
            
                end
                    
             	if isempty(dcorrdatamotor)
                        
                 	dcmotor = 0;
                        
                else
        
                  	dcmotor = dice_coefficient_mod(dcorrdatamotor(:,1),dcorrdatamotor(:,2));
                        
                end
                
             	if isempty(dcorrdatamixed)
                        
                 	dcmixed = 0;
                        
                else
        
                  	dcmixed = dice_coefficient_mod(dcorrdatamixed(:,1),dcorrdatamixed(:,2));
                        
                end
                
              	dcorrdataalltasks = [];
        
             	for q = 1:length(alltasksdatrest)
            
                   	if alltasksdatrest(:,q) > 0 && alltasksdat(:,q) > 0
                
                     	dcorrdataalltasks = [dcorrdataalltasks;1 1];
                
                  	elseif alltasksdatrest(:,q) > 0
                
                     	dcorrdataalltasks = [dcorrdataalltasks;1 0];
                
                	elseif alltasksdat(:,q) > 0
                
                      	dcorrdataalltasks = [dcorrdataalltasks;0 1];
                
                    end
            
                end
                    
             	if isempty(dcorrdataalltasks)
                        
                 	dcalltasks = 0;
                        
                else
        
                  	dcalltasks = dice_coefficient_mod(dcorrdataalltasks(:,1),dcorrdataalltasks(:,2));
                        
                end
        
              	DiceCorrssimmem = [DiceCorrssimmem; dcmem];
                DiceCorrssimmixed = [DiceCorrssimmixed; dcmixed];
                DiceCorrssimmotor = [DiceCorrssimmotor; dcmotor];
                DiceCorrssimalltasks = [DiceCorrssimalltasks; dcalltasks];
                        
              	if l == 1
                            
                  	DiceCorrSubsmem(Count,m) = dcmem;
                    DiceCorrSubsmixed(Count,m) = dcmixed;
                    DiceCorrSubsmotor(Count,m) = dcmotor;
                    DiceCorrSubsalltasks(Count,m) = dcalltasks;
                            
                else
                            
                 	DiceCorrSubsmem(Count+(nsubs-1),m) = dcmem;
                    DiceCorrSubsmixed(Count+(nsubs-1),m) = dcmixed;
                    DiceCorrSubsmotor(Count+(nsubs-1),m) = dcmotor;
                    DiceCorrSubsalltasks(Count+(nsubs-1),m) = dcalltasks;
                            
                end
                        
                        
            end
        end
    end
            
	DiceCorrspvalmem = prctile(DiceCorrssimmem, 95);
    DiceCorrspvalmixed = prctile(DiceCorrssimmixed, 95);
    DiceCorrspvalmotor = prctile(DiceCorrssimmotor, 95);
    DiceCorrspvalalltasks = prctile(DiceCorrssimalltasks, 95);
    
end



if Subplots == 1
    
    positions = [1.5:.5:3];
    scatterpos = [1.5:.5:2.5];
    
    figure;
    subplot(2,2,1)  %% All Tasks Plot
    scatter(scatterpos, [alltaskrestalltaskrestdicecoefs(8,:) alltaskalltaskdicecoefs(8,:) mean(alltaskalltaskrestdicecoefs(8,:))], 'filled', 'MarkerFaceColor', [1, 0.5, 0]);
    hold on
    scatter(scatterpos, [alltaskrestalltaskrestdicecoefs(7,:) alltaskalltaskdicecoefs(7,:) mean(alltaskalltaskrestdicecoefs(7,:))], 'filled', 'MarkerFaceColor', [1, 0, 1]);
    hold on
   	scatter(scatterpos, [alltaskrestalltaskrestdicecoefs(6,:) alltaskalltaskdicecoefs(6,:) mean(alltaskalltaskrestdicecoefs(6,:))], 'filled', 'MarkerFaceColor', [0.2, 1, 1]);
  	hold on
   	scatter(scatterpos, [alltaskrestalltaskrestdicecoefs(5,:) alltaskalltaskdicecoefs(5,:) mean(alltaskalltaskrestdicecoefs(5,:))], 'filled', 'MarkerFaceColor', [0, 0, 1]);
  	hold on
 	scatter(scatterpos, [alltaskrestalltaskrestdicecoefs(4,:) alltaskalltaskdicecoefs(4,:) mean(alltaskalltaskrestdicecoefs(4,:))], 'filled', 'MarkerFaceColor', [1, 0, 0]);
  	hold on
  	scatter(scatterpos, [alltaskrestalltaskrestdicecoefs(3,:) alltaskalltaskdicecoefs(3,:) mean(alltaskalltaskrestdicecoefs(3,:))], 'filled', 'MarkerFaceColor', [0, 1, 0]);
  	hold on
   	scatter(scatterpos, [alltaskrestalltaskrestdicecoefs(2,:) alltaskalltaskdicecoefs(2,:) mean(alltaskalltaskrestdicecoefs(2,:))], 'filled', 'MarkerFaceColor', [0.9, 0.9, 0]);
  	hold on
   	scatter(scatterpos, [alltaskrestalltaskrestdicecoefs(1,:) alltaskalltaskdicecoefs(1,:) mean(alltaskalltaskrestdicecoefs(1,:))], 'filled', 'MarkerFaceColor', [0, 0, 0]);
 	hold on
    scatter(repmat(positions(4),1,length(DiceCorrssimalltasks)), DiceCorrssimalltasks, 'jitter','on', 'jitterAmount', 0.1, 'MarkerEdgeColor', [0, 0, 0]);
   	set(gca,'xtick',[positions(1) positions(2) positions(3) positions(4)])
   	set(gca,'xticklabel',{'Rest to Rest', 'Task to Task', 'Rest to Task', 'Across Subjects'}, 'FontSize',12)
  	ylabel('Dice Correlation', 'FontSize',12)
    ylim([0 1]);
    xlim([1 3.5]);
    title('All Tasks Combined')
    subplot(2,2,2)  %% Memory Task Plot
    scatter(scatterpos, [memrestmemrestdicecoefs(8,:) memmemdicecoefs(8,:) mean(memmemrestdicecoefs(8,:))], 'filled', 'MarkerFaceColor', [1, 0.5, 0]);
    hold on
    scatter(scatterpos, [memrestmemrestdicecoefs(7,:) memmemdicecoefs(7,:) mean(memmemrestdicecoefs(7,:))], 'filled', 'MarkerFaceColor', [1, 0, 1]);
    hold on
   	scatter(scatterpos, [memrestmemrestdicecoefs(6,:) memmemdicecoefs(6,:) mean(memmemrestdicecoefs(6,:))], 'filled', 'MarkerFaceColor', [0.2, 1, 1]);
  	hold on
   	scatter(scatterpos, [memrestmemrestdicecoefs(5,:) memmemdicecoefs(5,:) mean(memmemrestdicecoefs(5,:))], 'filled', 'MarkerFaceColor', [0, 0, 1]);
  	hold on
 	scatter(scatterpos, [memrestmemrestdicecoefs(4,:) memmemdicecoefs(4,:) mean(memmemrestdicecoefs(4,:))], 'filled', 'MarkerFaceColor', [1, 0, 0]);
  	hold on
  	scatter(scatterpos, [memrestmemrestdicecoefs(3,:) memmemdicecoefs(3,:) mean(memmemrestdicecoefs(3,:))], 'filled', 'MarkerFaceColor', [0, 1, 0]);
  	hold on
   	scatter(scatterpos, [memrestmemrestdicecoefs(2,:) memmemdicecoefs(2,:) mean(memmemrestdicecoefs(2,:))], 'filled', 'MarkerFaceColor', [0.9, 0.9, 0]);
  	hold on
   	scatter(scatterpos, [memrestmemrestdicecoefs(1,:) memmemdicecoefs(1,:) mean(memmemrestdicecoefs(1,:))], 'filled', 'MarkerFaceColor', [0, 0, 0]);
 	hold on
    scatter(repmat(positions(4),1,length(DiceCorrssimmem)), DiceCorrssimmem, 'jitter','on', 'jitterAmount', 0.1, 'MarkerEdgeColor', [0, 0, 0]);
   	set(gca,'xtick',[positions(1) positions(2) positions(3) positions(4)])
   	set(gca,'xticklabel',{'Rest to Rest', 'Task to Task', 'Rest to Task', 'Across Subjects'}, 'FontSize',12)
  	ylabel('Dice Correlation', 'FontSize',12)
    ylim([0 1]);
    xlim([1 3.5]);
    title('Memory Task')
    subplot(2,2,3)  %% Mixed Task Plot
    scatter(scatterpos, [mixedrestmixedrestdicecoefs(8,:) mixedmixeddicecoefs(8,:) mean(mixedmixedrestdicecoefs(8,:))], 'filled', 'MarkerFaceColor', [1, 0.5, 0]);
    hold on
    scatter(scatterpos, [mixedrestmixedrestdicecoefs(7,:) mixedmixeddicecoefs(7,:) mean(mixedmixedrestdicecoefs(7,:))], 'filled', 'MarkerFaceColor', [1, 0, 1]);
    hold on
   	scatter(scatterpos, [mixedrestmixedrestdicecoefs(6,:) mixedmixeddicecoefs(6,:) mean(mixedmixedrestdicecoefs(6,:))], 'filled', 'MarkerFaceColor', [0.2, 1, 1]);
  	hold on
   	scatter(scatterpos, [mixedrestmixedrestdicecoefs(5,:) mixedmixeddicecoefs(5,:) mean(mixedmixedrestdicecoefs(5,:))], 'filled', 'MarkerFaceColor', [0, 0, 1]);
  	hold on
 	scatter(scatterpos, [mixedrestmixedrestdicecoefs(4,:) mixedmixeddicecoefs(4,:) mean(mixedmixedrestdicecoefs(4,:))], 'filled', 'MarkerFaceColor', [1, 0, 0]);
  	hold on
  	scatter(scatterpos, [mixedrestmixedrestdicecoefs(3,:) mixedmixeddicecoefs(3,:) mean(mixedmixedrestdicecoefs(3,:))], 'filled', 'MarkerFaceColor', [0, 1, 0]);
  	hold on
   	scatter(scatterpos, [mixedrestmixedrestdicecoefs(2,:) mixedmixeddicecoefs(2,:) mean(mixedmixedrestdicecoefs(2,:))], 'filled', 'MarkerFaceColor', [0.9, 0.9, 0]);
  	hold on
   	scatter(scatterpos, [mixedrestmixedrestdicecoefs(1,:) mixedmixeddicecoefs(1,:) mean(mixedmixedrestdicecoefs(1,:))], 'filled', 'MarkerFaceColor', [0, 0, 0]);
 	hold on
    scatter(repmat(positions(4),1,length(DiceCorrssimmixed)), DiceCorrssimmixed, 'jitter','on', 'jitterAmount', 0.1, 'MarkerEdgeColor', [0, 0, 0]);
   	set(gca,'xtick',[positions(1) positions(2) positions(3) positions(4)])
   	set(gca,'xticklabel',{'Rest to Rest', 'Task to Task', 'Rest to Task', 'Across Subjects'}, 'FontSize',12)
  	ylabel('Dice Correlation', 'FontSize',12)
    ylim([0 1]);
    xlim([1 3.5]);
    title('Mixed Task')
    subplot(2,2,4)  %% Motor Task Plot
    scatter(scatterpos, [motorrestmotorrestdicecoefs(8,:) motormotordicecoefs(8,:) mean(motormotorrestdicecoefs(8,:))], 'filled', 'MarkerFaceColor', [1, 0.5, 0]);
    hold on
    scatter(scatterpos, [motorrestmotorrestdicecoefs(7,:) motormotordicecoefs(7,:) mean(motormotorrestdicecoefs(7,:))], 'filled', 'MarkerFaceColor', [1, 0, 1]);
    hold on
   	scatter(scatterpos, [motorrestmotorrestdicecoefs(6,:) motormotordicecoefs(6,:) mean(motormotorrestdicecoefs(6,:))], 'filled', 'MarkerFaceColor', [0.2, 1, 1]);
  	hold on
   	scatter(scatterpos, [motorrestmotorrestdicecoefs(5,:) motormotordicecoefs(5,:) mean(motormotorrestdicecoefs(5,:))], 'filled', 'MarkerFaceColor', [0, 0, 1]);
  	hold on
 	scatter(scatterpos, [motorrestmotorrestdicecoefs(4,:) motormotordicecoefs(4,:) mean(motormotorrestdicecoefs(4,:))], 'filled', 'MarkerFaceColor', [1, 0, 0]);
  	hold on
  	scatter(scatterpos, [motorrestmotorrestdicecoefs(3,:) motormotordicecoefs(3,:) mean(motormotorrestdicecoefs(3,:))], 'filled', 'MarkerFaceColor', [0, 1, 0]);
  	hold on
   	scatter(scatterpos, [motorrestmotorrestdicecoefs(2,:) motormotordicecoefs(2,:) mean(motormotorrestdicecoefs(2,:))], 'filled', 'MarkerFaceColor', [0.9, 0.9, 0]);
  	hold on
   	scatter(scatterpos, [motorrestmotorrestdicecoefs(1,:) motormotordicecoefs(1,:) mean(motormotorrestdicecoefs(1,:))], 'filled', 'MarkerFaceColor', [0, 0, 0]);
 	hold on
    scatter(repmat(positions(4),1,length(DiceCorrssimmotor)), DiceCorrssimmotor, 'jitter','on', 'jitterAmount', 0.1, 'MarkerEdgeColor', [0, 0, 0]);
   	set(gca,'xtick',[positions(1) positions(2) positions(3) positions(4)])
   	set(gca,'xticklabel',{'Rest to Rest', 'Task to Task', 'Rest to Task', 'Across Subjects'}, 'FontSize',12)
  	ylabel('Dice Correlation', 'FontSize',12)
    ylim([0 1]);
    xlim([1 3.5]);
    title('Motor Task')
    
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [.4, .4, .6, .6]);
    
    if AddLegend == 1
    
        m = findobj(gca,'Type','scatter');

        hleg1 = legend([m(2:9)], 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC10', 'Location', 'NorthEast');
        %hleg1 = legend([m(2:10)], 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', 'NorthEast');
        %hleg1 = legend(m(2:10), 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', 'NorthEast');
        hleg1.FontSize = 14;
        
        filename = ['/AllSubjects_Subplots_AllTasksComp_SNRExclude_NoMSC09_Legend_' num2str(threshold) '_Percent_pvals.jpg'];
        
    else
        
        filename = ['/AllSubjects_Subplots_AllTasksComp_SNRExclude_NoMSC09_' num2str(threshold) '_Percent_pvals.jpg'];
 	
    end

	print(gcf,[outputdir filename],'-dpng','-r300');
        
	close gcf

end
    
if BoxPlot == 1
    
    if ExcludeMSC09 == 1
        
        x = [alltaskrestalltaskrestdicecoefs' alltaskalltaskdicecoefs' reshape(alltaskalltaskrestdicecoefs,1,size(alltaskalltaskrestdicecoefs,1)*size(alltaskalltaskrestdicecoefs,2)) DiceCorrssimalltasks' memrestmemrestdicecoefs' memmemdicecoefs' reshape(memmemrestdicecoefs,1,size(memmemrestdicecoefs,1)*size(memmemrestdicecoefs,2)) DiceCorrssimmem' mixedrestmixedrestdicecoefs' mixedmixeddicecoefs' reshape(mixedmixedrestdicecoefs,1,size(mixedmixedrestdicecoefs,1)*size(mixedmixedrestdicecoefs,2)) DiceCorrssimmixed' motorrestmotorrestdicecoefs' motormotordicecoefs' reshape(motormotorrestdicecoefs,1,size(motormotorrestdicecoefs,1)*size(motormotorrestdicecoefs,2)) DiceCorrssimmotor'];
        group = [repmat(1,1,length(alltaskrestalltaskrestdicecoefs)) repmat(1.25,1,length(alltaskalltaskdicecoefs)) repmat(1.5,1,size(alltaskalltaskrestdicecoefs,1)*size(alltaskalltaskrestdicecoefs,2)) repmat(1.75,1,length(DiceCorrssimalltasks)) repmat(2,1,length(memrestmemrestdicecoefs)) repmat(2.25,1,length(memmemdicecoefs)) repmat(2.5,1,size(memmemrestdicecoefs,1)*size(memmemrestdicecoefs,2)) repmat(2.75,1,length(DiceCorrssimmem)) repmat(3,1,length(mixedrestmixedrestdicecoefs)) repmat(3.25,1,length(mixedmixeddicecoefs)) repmat(3.5,1,size(mixedmixedrestdicecoefs,1)*size(mixedmixedrestdicecoefs,2)) repmat(3.75,1,length(DiceCorrssimmixed)) repmat(4,1,length(motorrestmotorrestdicecoefs)) repmat(4.25,1,length(motormotordicecoefs)) repmat(4.5,1,size(motormotorrestdicecoefs,1)*size(motormotorrestdicecoefs,2)) repmat(4.75,1,length(DiceCorrssimmotor))];
        positions = [1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4 4.25 4.5 4.75];
        scatterpos = [1 1.25 1.5 2 2.25 2.5 3 3.25 3.5 4 4.25 4.5];

        boxplot(x,group, 'BoxStyle', 'outline', 'positions', positions, 'Symbol', '');
        hold on
        scatter(scatterpos, [alltaskrestalltaskrestdicecoefs(8,:) alltaskalltaskdicecoefs(8,:) mean(alltaskalltaskrestdicecoefs(8,:)) memrestmemrestdicecoefs(8,:) memmemdicecoefs(8,:) mean(memmemrestdicecoefs(8,:)) mixedrestmixedrestdicecoefs(8,:) mixedmixeddicecoefs(8,:) mean(mixedmixedrestdicecoefs(8,:)) motorrestmotorrestdicecoefs(8,:) motormotordicecoefs(8,:) mean(motormotorrestdicecoefs(8,:))], 'filled', 'MarkerFaceColor', [1, 0.5, 0]);
        hold on    
        scatter(scatterpos, [alltaskrestalltaskrestdicecoefs(7,:) alltaskalltaskdicecoefs(7,:) mean(alltaskalltaskrestdicecoefs(7,:)) memrestmemrestdicecoefs(7,:) memmemdicecoefs(7,:) mean(memmemrestdicecoefs(7,:)) mixedrestmixedrestdicecoefs(7,:) mixedmixeddicecoefs(7,:) mean(mixedmixedrestdicecoefs(7,:)) motorrestmotorrestdicecoefs(7,:) motormotordicecoefs(7,:) mean(motormotorrestdicecoefs(7,:))], 'filled', 'MarkerFaceColor', [1, 0, 1]);
        hold on
        scatter(scatterpos, [alltaskrestalltaskrestdicecoefs(6,:) alltaskalltaskdicecoefs(6,:) mean(alltaskalltaskrestdicecoefs(6,:)) memrestmemrestdicecoefs(6,:) memmemdicecoefs(6,:) mean(memmemrestdicecoefs(6,:)) mixedrestmixedrestdicecoefs(6,:) mixedmixeddicecoefs(6,:) mean(mixedmixedrestdicecoefs(6,:)) motorrestmotorrestdicecoefs(6,:) motormotordicecoefs(6,:) mean(motormotorrestdicecoefs(6,:))], 'filled', 'MarkerFaceColor', [0.2, 1, 1]);
        hold on
        scatter(scatterpos, [alltaskrestalltaskrestdicecoefs(5,:) alltaskalltaskdicecoefs(5,:) mean(alltaskalltaskrestdicecoefs(5,:)) memrestmemrestdicecoefs(5,:) memmemdicecoefs(5,:) mean(memmemrestdicecoefs(5,:)) mixedrestmixedrestdicecoefs(5,:) mixedmixeddicecoefs(5,:) mean(mixedmixedrestdicecoefs(5,:)) motorrestmotorrestdicecoefs(5,:) motormotordicecoefs(5,:) mean(motormotorrestdicecoefs(5,:))], 'filled', 'MarkerFaceColor', [0, 0, 1]);
        hold on
        scatter(scatterpos, [alltaskrestalltaskrestdicecoefs(4,:) alltaskalltaskdicecoefs(4,:) mean(alltaskalltaskrestdicecoefs(4,:)) memrestmemrestdicecoefs(4,:) memmemdicecoefs(4,:) mean(memmemrestdicecoefs(4,:)) mixedrestmixedrestdicecoefs(4,:) mixedmixeddicecoefs(4,:) mean(mixedmixedrestdicecoefs(4,:)) motorrestmotorrestdicecoefs(4,:) motormotordicecoefs(4,:) mean(motormotorrestdicecoefs(4,:))], 'filled', 'MarkerFaceColor', [1, 0, 0]);
        hold on
        scatter(scatterpos, [alltaskrestalltaskrestdicecoefs(3,:) alltaskalltaskdicecoefs(3,:) mean(alltaskalltaskrestdicecoefs(3,:)) memrestmemrestdicecoefs(3,:) memmemdicecoefs(3,:) mean(memmemrestdicecoefs(3,:)) mixedrestmixedrestdicecoefs(3,:) mixedmixeddicecoefs(3,:) mean(mixedmixedrestdicecoefs(3,:)) motorrestmotorrestdicecoefs(3,:) motormotordicecoefs(3,:) mean(motormotorrestdicecoefs(3,:))], 'filled', 'MarkerFaceColor', [0, 1, 0]);
        hold on
        scatter(scatterpos, [alltaskrestalltaskrestdicecoefs(2,:) alltaskalltaskdicecoefs(2,:) mean(alltaskalltaskrestdicecoefs(2,:)) memrestmemrestdicecoefs(2,:) memmemdicecoefs(2,:) mean(memmemrestdicecoefs(2,:)) mixedrestmixedrestdicecoefs(2,:) mixedmixeddicecoefs(2,:) mean(mixedmixedrestdicecoefs(2,:)) motorrestmotorrestdicecoefs(2,:) motormotordicecoefs(2,:) mean(motormotorrestdicecoefs(2,:))], 'filled', 'MarkerFaceColor', [0.9, 0.9, 0]);
        hold on
        scatter(scatterpos, [alltaskrestalltaskrestdicecoefs(1,:) alltaskalltaskdicecoefs(1,:) mean(alltaskalltaskrestdicecoefs(1,:)) memrestmemrestdicecoefs(1,:) memmemdicecoefs(1,:) mean(memmemrestdicecoefs(1,:)) mixedrestmixedrestdicecoefs(1,:) mixedmixeddicecoefs(1,:) mean(mixedmixedrestdicecoefs(1,:)) motorrestmotorrestdicecoefs(1,:) motormotordicecoefs(1,:) mean(motormotorrestdicecoefs(1,:))], 'filled', 'MarkerFaceColor', [0, 0, 0]);
        hold on
        scatter([repmat(1.75,1,length(DiceCorrssimalltasks)) repmat(2.75,1,length(DiceCorrssimmem)) repmat(3.75,1,length(DiceCorrssimmixed)) repmat(4.75,1,length(DiceCorrssimmotor))], [DiceCorrssimalltasks' DiceCorrssimmem' DiceCorrssimmixed' DiceCorrssimmotor'], 'jitter','on', 'jitterAmount', 0.1, 'MarkerEdgeColor', [0, 0, 0]);
        hold on

        set(gca,'xtick',[mean(positions(1:4)) mean(positions(5:8)) mean(positions(9:12)) mean(positions(13:16))])
        set(gca,'xticklabel',{'All Tasks Combined', 'Memory Task', 'Mixed Task', 'Motor Task'}, 'FontSize',12)
        ylabel('Dice Correlation', 'FontSize',12)

        color = ['k', 'b', 'g', 'r', 'k', 'b', 'g', 'r', 'k', 'b', 'g', 'r', 'k', 'b', 'g', 'r'];
        %color = ['k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k'];
        h = findobj(gca,'Tag','Box');

        for j=1:length(h)

            patch(get(h(j),'XData'),get(h(j),'YData'),color(j), 'FaceColor', 'none', 'EdgeColor', color(j));

        end

        c = get(gca, 'Children');

        m = findobj(gca,'Type','scatter');

        hleg1 = legend([m(2:9); c(1:4)], 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC10', 'Rest Comparison', 'Task Comparison', 'Task-Rest Comparison', 'Between Subject Comparison', 'Location', 'NorthEast');
        %hleg1 = legend([m(2:10)], 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', 'NorthEast');
        %hleg1 = legend(m(2:10), 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', 'NorthEast');
        hleg1.FontSize = 14;
        ylim([0 1]);

        set(gcf, 'Units', 'Normalized', 'OuterPosition', [.3, .6, .6, .9]);

        %hleg1 = legend(c(1:3), 'Between State Comparison', 'Within State Comparison', 'Between Subject Comparison', 'Location', 'NorthWest');

        %legendpos = 'NorthEast';

        %hleg2 = legend(m(2:10), 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', legendpos);
        
     	filename = ['/AllSubjects_Boxplot_AllTasksComp_SNRExclude_NoMSC09_' num2str(threshold) '_Percent_pvals.jpg'];
        
        print(gcf,[outputdir filename],'-dpng','-r300');
        
        %saveas(gcf,[outputdir filename])
        
        close gcf
        
    else

        x = [alltaskrestalltaskrestdicecoefs' alltaskalltaskdicecoefs' reshape(alltaskalltaskrestdicecoefs,1,size(alltaskalltaskrestdicecoefs,1)*size(alltaskalltaskrestdicecoefs,2)) DiceCorrssimalltasks' memrestmemrestdicecoefs' memmemdicecoefs' reshape(memmemrestdicecoefs,1,size(memmemrestdicecoefs,1)*size(memmemrestdicecoefs,2)) DiceCorrssimmem' mixedrestmixedrestdicecoefs' mixedmixeddicecoefs' reshape(mixedmixedrestdicecoefs,1,size(mixedmixedrestdicecoefs,1)*size(mixedmixedrestdicecoefs,2)) DiceCorrssimmixed' motorrestmotorrestdicecoefs' motormotordicecoefs' reshape(motormotorrestdicecoefs,1,size(motormotorrestdicecoefs,1)*size(motormotorrestdicecoefs,2)) DiceCorrssimmotor'];
        group = [repmat(1,1,length(alltaskrestalltaskrestdicecoefs)) repmat(1.25,1,length(alltaskalltaskdicecoefs)) repmat(1.5,1,size(alltaskalltaskrestdicecoefs,1)*size(alltaskalltaskrestdicecoefs,2)) repmat(1.75,1,length(DiceCorrssimalltasks)) repmat(2,1,length(memrestmemrestdicecoefs)) repmat(2.25,1,length(memmemdicecoefs)) repmat(2.5,1,size(memmemrestdicecoefs,1)*size(memmemrestdicecoefs,2)) repmat(2.75,1,length(DiceCorrssimmem)) repmat(3,1,length(mixedrestmixedrestdicecoefs)) repmat(3.25,1,length(mixedmixeddicecoefs)) repmat(3.5,1,size(mixedmixedrestdicecoefs,1)*size(mixedmixedrestdicecoefs,2)) repmat(3.75,1,length(DiceCorrssimmixed)) repmat(4,1,length(motorrestmotorrestdicecoefs)) repmat(4.25,1,length(motormotordicecoefs)) repmat(4.5,1,size(motormotorrestdicecoefs,1)*size(motormotorrestdicecoefs,2)) repmat(4.75,1,length(DiceCorrssimmotor))];
        positions = [1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4 4.25 4.5 4.75];
        scatterpos = [1 1.25 1.5 2 2.25 2.5 3 3.25 3.5 4 4.25 4.5];

        boxplot(x,group, 'BoxStyle', 'outline', 'positions', positions, 'Symbol', '');
        hold on
        scatter(scatterpos, [alltaskrestalltaskrestdicecoefs(9,:) alltaskalltaskdicecoefs(9,:) mean(alltaskalltaskrestdicecoefs(9,:)) memrestmemrestdicecoefs(9,:) memmemdicecoefs(9,:) mean(memmemrestdicecoefs(9,:)) mixedrestmixedrestdicecoefs(9,:) mixedmixeddicecoefs(9,:) mean(mixedmixedrestdicecoefs(9,:)) motorrestmotorrestdicecoefs(9,:) motormotordicecoefs(9,:) mean(motormotorrestdicecoefs(9,:))], 'filled', 'MarkerFaceColor', [1, 0.5, 0]);
        hold on
        scatter(scatterpos, [alltaskrestalltaskrestdicecoefs(8,:) alltaskalltaskdicecoefs(8,:) mean(alltaskalltaskrestdicecoefs(8,:)) memrestmemrestdicecoefs(8,:) memmemdicecoefs(8,:) mean(memmemrestdicecoefs(8,:)) mixedrestmixedrestdicecoefs(8,:) mixedmixeddicecoefs(8,:) mean(mixedmixedrestdicecoefs(8,:)) motorrestmotorrestdicecoefs(8,:) motormotordicecoefs(8,:) mean(motormotorrestdicecoefs(8,:))], 'filled', 'MarkerFaceColor', [0, 0.6, 0.6]);
        hold on    
        scatter(scatterpos, [alltaskrestalltaskrestdicecoefs(7,:) alltaskalltaskdicecoefs(7,:) mean(alltaskalltaskrestdicecoefs(7,:)) memrestmemrestdicecoefs(7,:) memmemdicecoefs(7,:) mean(memmemrestdicecoefs(7,:)) mixedrestmixedrestdicecoefs(7,:) mixedmixeddicecoefs(7,:) mean(mixedmixedrestdicecoefs(7,:)) motorrestmotorrestdicecoefs(7,:) motormotordicecoefs(7,:) mean(motormotorrestdicecoefs(7,:))], 'filled', 'MarkerFaceColor', [1, 0, 1]);
        hold on
        scatter(scatterpos, [alltaskrestalltaskrestdicecoefs(6,:) alltaskalltaskdicecoefs(6,:) mean(alltaskalltaskrestdicecoefs(6,:)) memrestmemrestdicecoefs(6,:) memmemdicecoefs(6,:) mean(memmemrestdicecoefs(6,:)) mixedrestmixedrestdicecoefs(6,:) mixedmixeddicecoefs(6,:) mean(mixedmixedrestdicecoefs(6,:)) motorrestmotorrestdicecoefs(6,:) motormotordicecoefs(6,:) mean(motormotorrestdicecoefs(6,:))], 'filled', 'MarkerFaceColor', [0.2, 1, 1]);
        hold on
        scatter(scatterpos, [alltaskrestalltaskrestdicecoefs(5,:) alltaskalltaskdicecoefs(5,:) mean(alltaskalltaskrestdicecoefs(5,:)) memrestmemrestdicecoefs(5,:) memmemdicecoefs(5,:) mean(memmemrestdicecoefs(5,:)) mixedrestmixedrestdicecoefs(5,:) mixedmixeddicecoefs(5,:) mean(mixedmixedrestdicecoefs(5,:)) motorrestmotorrestdicecoefs(5,:) motormotordicecoefs(5,:) mean(motormotorrestdicecoefs(5,:))], 'filled', 'MarkerFaceColor', [0, 0, 1]);
        hold on
        scatter(scatterpos, [alltaskrestalltaskrestdicecoefs(4,:) alltaskalltaskdicecoefs(4,:) mean(alltaskalltaskrestdicecoefs(4,:)) memrestmemrestdicecoefs(4,:) memmemdicecoefs(4,:) mean(memmemrestdicecoefs(4,:)) mixedrestmixedrestdicecoefs(4,:) mixedmixeddicecoefs(4,:) mean(mixedmixedrestdicecoefs(4,:)) motorrestmotorrestdicecoefs(4,:) motormotordicecoefs(4,:) mean(motormotorrestdicecoefs(4,:))], 'filled', 'MarkerFaceColor', [1, 0, 0]);
        hold on
        scatter(scatterpos, [alltaskrestalltaskrestdicecoefs(3,:) alltaskalltaskdicecoefs(3,:) mean(alltaskalltaskrestdicecoefs(3,:)) memrestmemrestdicecoefs(3,:) memmemdicecoefs(3,:) mean(memmemrestdicecoefs(3,:)) mixedrestmixedrestdicecoefs(3,:) mixedmixeddicecoefs(3,:) mean(mixedmixedrestdicecoefs(3,:)) motorrestmotorrestdicecoefs(3,:) motormotordicecoefs(3,:) mean(motormotorrestdicecoefs(3,:))], 'filled', 'MarkerFaceColor', [0, 1, 0]);
        hold on
        scatter(scatterpos, [alltaskrestalltaskrestdicecoefs(2,:) alltaskalltaskdicecoefs(2,:) mean(alltaskalltaskrestdicecoefs(2,:)) memrestmemrestdicecoefs(2,:) memmemdicecoefs(2,:) mean(memmemrestdicecoefs(2,:)) mixedrestmixedrestdicecoefs(2,:) mixedmixeddicecoefs(2,:) mean(mixedmixedrestdicecoefs(2,:)) motorrestmotorrestdicecoefs(2,:) motormotordicecoefs(2,:) mean(motormotorrestdicecoefs(2,:))], 'filled', 'MarkerFaceColor', [0.9, 0.9, 0]);
        hold on
        scatter(scatterpos, [alltaskrestalltaskrestdicecoefs(1,:) alltaskalltaskdicecoefs(1,:) mean(alltaskalltaskrestdicecoefs(1,:)) memrestmemrestdicecoefs(1,:) memmemdicecoefs(1,:) mean(memmemrestdicecoefs(1,:)) mixedrestmixedrestdicecoefs(1,:) mixedmixeddicecoefs(1,:) mean(mixedmixedrestdicecoefs(1,:)) motorrestmotorrestdicecoefs(1,:) motormotordicecoefs(1,:) mean(motormotorrestdicecoefs(1,:))], 'filled', 'MarkerFaceColor', [0, 0, 0]);
        hold on
        scatter([repmat(1.75,1,length(DiceCorrssimalltasks)) repmat(2.75,1,length(DiceCorrssimmem)) repmat(3.75,1,length(DiceCorrssimmixed)) repmat(4.75,1,length(DiceCorrssimmotor))], [DiceCorrssimalltasks' DiceCorrssimmem' DiceCorrssimmixed' DiceCorrssimmotor'], 'jitter','on', 'jitterAmount', 0.1, 'MarkerEdgeColor', [1, .25, 0]);
        hold on

        set(gca,'xtick',[mean(positions(1:4)) mean(positions(5:8)) mean(positions(9:12)) mean(positions(13:16))])
        set(gca,'xticklabel',{'All Tasks Combined', 'Memory Task', 'Mixed Task', 'Motor Task'}, 'FontSize',12)
        ylabel('Dice Correlation', 'FontSize',12)

        color = ['k', 'b', 'g', 'r', 'k', 'b', 'g', 'r', 'k', 'b', 'g', 'r', 'k', 'b', 'g', 'r'];
        %color = ['k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k'];
        h = findobj(gca,'Tag','Box');

        for j=1:length(h)

            patch(get(h(j),'XData'),get(h(j),'YData'),color(j), 'FaceColor', 'none', 'EdgeColor', color(j));

        end

        c = get(gca, 'Children');

        m = findobj(gca,'Type','scatter');

        hleg1 = legend([m(2:10); c(1:4)], 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Rest Comparison', 'Task Comparison', 'Task-Rest Comparison', 'Between Subject Comparison', 'Location', 'NorthEast');
        %hleg1 = legend([m(2:10)], 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', 'NorthEast');
        %hleg1 = legend(m(2:10), 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', 'NorthEast');
        hleg1.FontSize = 14;
        ylim([0 1]);

        set(gcf, 'Units', 'Normalized', 'OuterPosition', [.3, .6, .6, .9]);

        %hleg1 = legend(c(1:3), 'Between State Comparison', 'Within State Comparison', 'Between Subject Comparison', 'Location', 'NorthWest');

        %legendpos = 'NorthEast';

        %hleg2 = legend(m(2:10), 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', legendpos);
        
     	filename = '/AllSubjects_Boxplot_AllTasksComp_SNRExclude_pvals.jpg';
        
        saveas(gcf,[outputdir filename])
        
        close gcf

    end        
end