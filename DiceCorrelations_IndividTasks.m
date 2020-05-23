%% Compute and plot dice correlations across individual tasks

% This function is for calculating the dice correlations for spatial
% overlap of network variants between 3 different tasks and all of the
% tasks combined. It also calculates the spatial overlap between states
% (task vs. rest) across subjects. All values are plotted at the end of the
% script. The script requires separate variant maps for each individual
% task, combined tasks, and corresponding resting state maps. Parcel maps
% are required for excluding variants by size if desired. Data are
% referenced via text files (see below for format).


clear all

%% Set variable values

memmemdicecoefs = [];  %% Set empty variables to hold correlation values
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

memfilestaskeven = [];  %% Set empty variables to hold vectors for dice correlation
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

permutatevals = 1;  %% Toggles whether to calculate permutations across subjects
Subplots = 1;  %% Toggles whether to plot each task as a subplot in a larger figure
AddLegend = 1;  %% Toggles whether to plot legend in figure
threshold = 5;  %% Toggles threshold to create variants
SNRMasks = 1;  %% Toggles whether to exclude areas of low SNR via SNR masks
SizeExclude = 1;  %% Toggles whether to exclude variants smaller than criterion
minsize = 50;  %% Minimum variant size for size exclusion
ExcludeMSC09 = 1;  %% Toggles where to exclude MSC09

outputdir = '/Users/briankraus/Desktop/Correct_Variant_Maps/TaskVariantPlots';  %% Output directory for plots
SNRMaskPath = '/Users/briankraus/Desktop/MSC_SNR_Maps/CIFTI_Files/MSC_Template/Individual_Session_Map/';  %% Path to SNR masks

comparisons = {'matched-samedata-task-task', 'matched-samedata-task-rest', 'matched-samedata-taskrest-taskrest', 'matched-samedata-alltask-alltask', 'matched-samedata-alltask-rest', 'matched-samedata-alltaskrest-alltaskrest'};  %% Toggles through comparisons
% 1. individual task to individual task
% 2. individual task to rest
% 3. rest to rest (number of samples per session matched to each individual task)
% 4. all tasks combined to all tasks combined
% 5. all tasks combined to rest
% 6. rest to rest (number of samples per session matched to all tasks combined)

%% Run Script

% Load variant maps, exclude by size and SNR, store variant vertices to
% calculate spatial overlap

% Text files should have the format 'pathtofile subjectid taskid' with
% spaces as delimiters. Text lists for odd/even_files should reference
% variants maps and odd/even_masks should reference parcel files which have
% a unique identifier for each vertex within each variant.

for filetype = 1:4  %% Number of iterations for each file type (see labels below)

    if filetype == 1   %% Individual task files
        
        [even_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_task_samedata_consec_even_varmaps_NoMSC09.txt','%s%s%s');
        
        [odd_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_task_samedata_consec_odd_varmaps_NoMSC09.txt','%s%s%s');
        
        [even_masks, ~, ~] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_task_samedata_consec_even_varmaps_NoMSC09_' num2str(threshold) '_SizeOnly.txt'],'%s%s%s');
        
        [odd_masks, ~, ~] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_task_samedata_consec_odd_varmaps_NoMSC09_' num2str(threshold) '_SizeOnly.txt'],'%s%s%s');
        
    elseif filetype == 2   %% Task-matched rest files (individual tasks)
        
        [even_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_task_samedata_consec_even_varmaps_NoMSC09.txt','%s%s%s');
        
        [odd_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_task_samedata_consec_odd_varmaps_NoMSC09.txt','%s%s%s');
        
        [even_masks, ~, ~] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_task_samedata_consec_even_varmaps_NoMSC09_' num2str(threshold) '_SizeOnly.txt'],'%s%s%s');
        
        [odd_masks, ~, ~] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_task_samedata_consec_odd_varmaps_NoMSC09_' num2str(threshold) '_SizeOnly.txt'],'%s%s%s');
        
    elseif filetype == 3   %% All Tasks Combined files
        
        [even_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_samedata_consec_even_varmaps_NoMSC09.txt','%s%s%s');
        
        [odd_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_samedata_consec_odd_varmaps_NoMSC09.txt','%s%s%s');
        
        [even_masks, ~, ~] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_samedata_consec_even_varmaps_NoMSC09_' num2str(threshold) '_SizeOnly.txt'],'%s%s%s');
        
        [odd_masks, ~, ~] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_samedata_consec_odd_varmaps_NoMSC09_' num2str(threshold) '_SizeOnly.txt'],'%s%s%s');
        
    elseif filetype == 4 %% All Tasks Combined-matched rest files (all tasks combined)
        
        [even_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_alltask_samedata_consec_even_varmaps_NoMSC09.txt','%s%s%s');
        
        [odd_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_alltask_samedata_consec_odd_varmaps_NoMSC09.txt','%s%s%s');
        
        [even_masks, ~, ~] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_alltask_samedata_consec_even_varmaps_NoMSC09_' num2str(threshold) '_SizeOnly.txt'],'%s%s%s');
        
        [odd_masks, ~, ~] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_alltask_samedata_consec_odd_varmaps_NoMSC09_' num2str(threshold) '_SizeOnly.txt'],'%s%s%s');
        
    end

    for a = 1:length(even_files)  %% Loop through each file in text list

        sub1 = subjects1{a};
        sub2 = subjects2{a};
        even_file = even_files{a};
        odd_file = odd_files{a};
        task1 = tasks1{a};
        task2 = tasks2{a};

        cifti_even = ft_read_cifti_mod(even_file);      %% Load each file from textlist
        cifti_odd = ft_read_cifti_mod(odd_file);

        if SNRMasks == 1        %% Remove low SNR areas from consideration

            SNRmask = ft_read_cifti_mod([SNRMaskPath sub1 '/' sub1 '__SNRMap__AllDataConcatenated.dscalar.nii']); 

            SNRmask.data = SNRmask.data(1:59412,:);
            LowSNR = find(SNRmask.data < 750);

            cifti_even.data(LowSNR,:) = NaN;
            cifti_odd.data(LowSNR,:) = NaN;

        end
        
    	cifti_even_threshold = find(cifti_even.data < prctile(cifti_even.data,threshold));  %% Select data below threshold chosen
      	cifti_odd_threshold = find(cifti_odd.data < prctile(cifti_odd.data,threshold));
        
        cifti_even.data(:,:) = 0;       %% Binarize data
        cifti_odd.data(:,:) = 0;
        cifti_even.data(cifti_even_threshold,:) = 1;
        cifti_odd.data(cifti_odd_threshold,:) = 1;

        if SizeExclude == 1     %% Exclude variants below size threshold

            even_mask = even_masks{a};      %% Load masks
            odd_mask = odd_masks{a};

            cifti_even_mask = ft_read_cifti_mod(even_mask);
            cifti_odd_mask = ft_read_cifti_mod(odd_mask);

            allvars_odd = unique(cifti_odd_mask.data);  %% Get list of variants
            allvars_odd(1) = [];
            allvars_even = unique(cifti_even_mask.data);
            allvars_even(1) = [];

            removevars_odd = [];
            removevars_even = [];

            for z = 1:length(allvars_odd)   %% Check length of each variant

                if length(find(cifti_odd.data == allvars_odd(z))) < minsize

                    removevars_odd = [removevars_odd allvars_odd(z)];

                end
            end

            for h = 1:length(cifti_odd.data)    %% Remove variants below size threshold

                if ismember(cifti_odd.data(h),removevars_odd)

                    cifti_odd.data(h) = 0;

                end
            end

            for z = 1:length(allvars_even)   %% Check length of each variant

                if length(find(cifti_even.data == allvars_even(z))) < minsize

                    removevars_even = [removevars_even allvars_even(z)];

                end
            end

            for h = 1:length(cifti_even.data)    %% Remove variants below size threshold

                if ismember(cifti_even.data(h),removevars_even)

                    cifti_even.data(h) = 0;

                end
            end
        end
        
        % Store data for dice correlations after creating final variants
    
        if filetype == 1   %% Individual task files
            
            if strcmp(task1,'mem')
                
                memfilestaskeven = [memfilestaskeven; cifti_even.data'];
                memfilestaskodd = [memfilestaskodd; cifti_odd.data'];
                
            elseif strcmp(task1,'mixed')
                
                mixedfilestaskeven = [mixedfilestaskeven; cifti_even.data'];
                mixedfilestaskodd = [mixedfilestaskodd; cifti_odd.data'];
                
            elseif strcmp(task1,'motor')
                
                motorfilestaskeven = [motorfilestaskeven; cifti_even.data'];
                motorfilestaskodd = [motorfilestaskodd; cifti_odd.data'];
                
            end
            
        elseif filetype == 2   %% Individual task-matched rest files
                
            if strcmp(task1,'memrest')
                
                memfilesresteven = [memfilesresteven; cifti_even.data'];
                memfilesrestodd = [memfilesrestodd; cifti_odd.data'];
                
            elseif strcmp(task1,'mixedrest')
                
                mixedfilesresteven = [mixedfilesresteven; cifti_even.data'];
                mixedfilesrestodd = [mixedfilesrestodd; cifti_odd.data'];
                
            elseif strcmp(task1,'motorrest')
                
                motorfilesresteven = [motorfilesresteven; cifti_even.data'];
                motorfilesrestodd = [motorfilesrestodd; cifti_odd.data'];
                
            end
            
        elseif filetype == 3   %% All tasks combined files
            
            alltaskfilestaskeven = [alltaskfilestaskeven; cifti_even.data'];
            alltaskfilestaskodd = [alltaskfilestaskodd; cifti_odd.data'];
            
        elseif filetype == 4   %% All tasks combined-matched rest files
            
            alltaskfilesresteven = [alltaskfilesresteven; cifti_even.data'];
            alltaskfilesrestodd = [alltaskfilesrestodd; cifti_odd.data'];
            
        end
    end
end


% Calculate dice overlap for appropriate comparisons

for g = 1:numel(comparisons)
    
    comp = comparisons(g);
    
    for dataset = 1:length(subjects1)    %% Loop through each dataset

        if strcmp(comp, 'matched-samedata-task-task')   %% Individual task to individual task

            dcorrdatamem = binarize_vectors(memfilestaskeven,memfilestaskodd,dataset);          %% Binarize vectors

            dcmem = dice_coefficient_mod(dcorrdatamem(:,1),dcorrdatamem(:,2));      %% Calculate dice correlation on binary data

            memmemdicecoefs = [memmemdicecoefs; dcmem];        %% Store dice values

            % Mixed Task

            dcorrdatamixed = binarize_vectors(mixedfilestaskeven,mixedfilestaskodd,dataset);

            dcmixed = dice_coefficient_mod(dcorrdatamixed(:,1),dcorrdatamixed(:,2));

            mixedmixeddicecoefs = [mixedmixeddicecoefs; dcmixed];

            % Motor Task

            dcorrdatamotor = binarize_vectors(motorfilestaskeven,motorfilestaskodd,dataset);

            dcmotor = dice_coefficient_mod(dcorrdatamotor(:,1),dcorrdatamotor(:,2));

            motormotordicecoefs = [motormotordicecoefs; dcmotor];

        elseif strcmp(comp, 'matched-samedata-task-rest')       %% Individual task to rest (matched to individual tasks)

            % Mem Task

            dcorrdatamemresteven = binarize_vectors(memfilesresteven,memfilestaskodd,dataset);          %% Binarize vectors
            dcorrdatamemrestodd = binarize_vectors(memfilestaskeven,memfilesrestodd,dataset);

            dcmemresteven = dice_coefficient_mod(dcorrdatamemresteven(:,1),dcorrdatamemresteven(:,2));      %% Calculate dice correlation on binary data
            dcmemrestodd = dice_coefficient_mod(dcorrdatamemrestodd(:,1),dcorrdatamemrestodd(:,2));

            memmemrestdicecoefs = [memmemrestdicecoefs; dcmemresteven dcmemrestodd];        %% Store dice values

            % Mixed Task

            dcorrdatamixedresteven = binarize_vectors(mixedfilesresteven,mixedfilestaskodd,dataset);
            dcorrdatamixedrestodd = binarize_vectors(mixedfilestaskeven,mixedfilesrestodd,dataset);

            dcmixedresteven = dice_coefficient_mod(dcorrdatamixedresteven(:,1),dcorrdatamixedresteven(:,2));
            dcmixedrestodd = dice_coefficient_mod(dcorrdatamixedrestodd(:,1),dcorrdatamixedrestodd(:,2));

            mixedmixedrestdicecoefs = [mixedmixedrestdicecoefs; dcmixedresteven dcmixedrestodd];

            % Motor Task

            dcorrdatamotorresteven = binarize_vectors(motorfilesresteven,motorfilestaskodd,dataset);
            dcorrdatamotorrestodd = binarize_vectors(motorfilestaskeven,motorfilesrestodd,dataset);

            dcmotorresteven = dice_coefficient_mod(dcorrdatamotorresteven(:,1),dcorrdatamotorresteven(:,2));
            dcmotorrestodd = dice_coefficient_mod(dcorrdatamotorrestodd(:,1),dcorrdatamotorrestodd(:,2));

            motormotorrestdicecoefs = [motormotorrestdicecoefs; dcmotorresteven dcmotorrestodd];

        elseif strcmp(comp, 'matched-samedata-taskrest-taskrest')       %% Rest to rest (matched to individual tasks)

            % Mem Task

            dcorrdatamemrest = binarize_vectors(memfilesresteven,memfilesrestodd,dataset);

            dcmemrest = dice_coefficient_mod(dcorrdatamemrest(:,1),dcorrdatamemrest(:,2));

            memrestmemrestdicecoefs = [memrestmemrestdicecoefs; dcmemrest];

            % Mixed Task

            dcorrdatamixedrest = binarize_vectors(mixedfilesresteven,mixedfilesrestodd,dataset);

            dcmixedrest = dice_coefficient_mod(dcorrdatamixedrest(:,1),dcorrdatamixedrest(:,2));

            mixedrestmixedrestdicecoefs = [mixedrestmixedrestdicecoefs; dcmixedrest];

            % Motor Task

            dcorrdatamotorrest = binarize_vectors(motorfilesresteven,motorfilesrestodd,dataset);

            dcmotorrest = dice_coefficient_mod(dcorrdatamotorrest(:,1),dcorrdatamotorrest(:,2));

            motorrestmotorrestdicecoefs = [motorrestmotorrestdicecoefs; dcmotorrest];

        elseif strcmp(comp, 'matched-samedata-alltask-alltask')         %% All tasks combined to all tasks combined

            % All Tasks Combined

            dcorrdataalltasks = binarize_vectors(alltaskfilestaskeven,alltaskfilestaskodd,dataset);

            dcalltasks = dice_coefficient_mod(dcorrdataalltasks(:,1),dcorrdataalltasks(:,2));

            alltaskalltaskdicecoefs = [alltaskalltaskdicecoefs; dcalltasks];

        elseif strcmp(comp, 'matched-samedata-alltask-rest')         %% All tasks combined to rest (matched to all tasks combined)

            % All Tasks Combined

            dcorrdataalltasksresteven = binarize_vectors(alltaskfilesresteven,alltaskfilestaskodd,dataset);
            dcorrdataalltasksrestodd = binarize_vectors(alltaskfilestaskeven,alltaskfilesrestodd,dataset);

            dcalltasksresteven = dice_coefficient_mod(dcorrdataalltasksresteven(:,1),dcorrdataalltasksresteven(:,2));
            dcalltasksrestodd = dice_coefficient_mod(dcorrdataalltasksrestodd(:,1),dcorrdataalltasksrestodd(:,2));

            alltaskalltaskrestdicecoefs = [alltaskalltaskrestdicecoefs; dcalltasksresteven dcalltasksrestodd];

        elseif strcmp(comp, 'matched-samedata-alltaskrest-alltaskrest')         %% Rest to rest  (matched to all tasks combined)

            % All Tasks Combined Rest

            dcorrdataalltasksrest = binarize_vectors(alltaskfilesresteven,alltaskfilesrestodd,dataset);

            dcalltasksrest = dice_coefficient_mod(dcorrdataalltasksrest(:,1),dcorrdataalltasksrest(:,2));

            alltaskrestalltaskrestdicecoefs = [alltaskrestalltaskrestdicecoefs; dcalltasksrest];

        end
    end
end


% Get all spatial permutations across subjects and states

if permutatevals == 1
    
  	nsubs = size(memfilestaskeven,1);
            
  	DiceCorrssimmem = [];       %% Create empty variables to store results
  	DiceCorrSubsmem = zeros(2*(nsubs-1),nsubs);
  	DiceCorrssimmixed = [];
  	DiceCorrSubsmixed = zeros(2*(nsubs-1),nsubs);
   	DiceCorrssimmotor = [];
  	DiceCorrSubsmotor = zeros(2*(nsubs-1),nsubs);
   	DiceCorrssimalltasks = [];
  	DiceCorrSubsalltasks = zeros(2*(nsubs-1),nsubs);

 	for splithalf = 1:2     %% Do both versions of split-half combinations
                
      	for subiter = 1:nsubs   %% Iterate through subjects
                    
          	loopsubs = [1:nsubs];       %% Loop through each subject
         	loopsubs(subiter) = [];     %% Remove each subject from across-subject permutation
                    
         	Count = 0;      %% Count number of loops to store data in correct place
                    
          	for subcompare = loopsubs    %% Get dice correlation for each participant to every other participant
                        
            	Count = Count+1;
                
                % Get binarized map for each participant and loop through
                % other participants to compare spatial overlap across
                % subjects

              	if splithalf == 1       %% First set of split-half comparisons
                            
                  	memdatrest = memfilesrestodd(subiter,:);
                 	memdat = memfilestaskeven(subcompare,:);
                    
                  	mixeddatrest = mixedfilesrestodd(subiter,:);
                 	mixeddat = mixedfilestaskeven(subcompare,:);
                    
                  	motordatrest = motorfilesrestodd(subiter,:);
                 	motordat = motorfilestaskeven(subcompare,:);
                    
                  	alltasksdatrest = alltaskfilesrestodd(subiter,:);
                 	alltasksdat = alltaskfilestaskeven(subcompare,:);
                            
                else        %% Second set of split-half comparisons
                            
                  	memdatrest = memfilesresteven(subiter,:);
                 	memdat = memfilestaskodd(subcompare,:);
                    
                  	mixeddatrest = mixedfilesresteven(subiter,:);
                 	mixeddat = mixedfilestaskodd(subcompare,:);
                    
                  	motordatrest = motorfilesresteven(subiter,:);
                 	motordat = motorfilestaskodd(subcompare,:);
                    
                  	alltasksdatrest = alltaskfilesresteven(subiter,:);
                 	alltasksdat = alltaskfilestaskodd(subcompare,:);
                            
                end
                
                % Across-subject permutations for mem task
                
                dcorrdatamem = binarize_vectors_perm(memdat,memdatrest);     
                    
             	if isempty(dcorrdatamem)    %% Catch if no vertices overlap
                        
                 	dcmem = 0;
                        
                else
        
                  	dcmem = dice_coefficient_mod(dcorrdatamem(:,1),dcorrdatamem(:,2));
                        
                end
                
                % Across-subject permutations for mixed task
                
                dcorrdatamixed = binarize_vectors_perm(mixeddat,mixeddatrest);

             	if isempty(dcorrdatamixed)
                        
                 	dcmixed = 0;
                        
                else
        
                  	dcmixed = dice_coefficient_mod(dcorrdatamixed(:,1),dcorrdatamixed(:,2));
                        
                end
                
                % Across-subject permutations for motor task
                
                dcorrdatamotor = binarize_vectors_perm(motordat,motordatrest);

             	if isempty(dcorrdatamotor)
                        
                 	dcmotor = 0;
                        
                else
        
                  	dcmotor = dice_coefficient_mod(dcorrdatamotor(:,1),dcorrdatamotor(:,2));
                        
                end
                
                % Across-subject permutations for combined tasks
                
                dcorrdataalltasks = binarize_vectors_perm(alltasksdat,alltasksdatrest);

             	if isempty(dcorrdataalltasks)
                        
                 	dcalltasks = 0;
                        
                else
        
                  	dcalltasks = dice_coefficient_mod(dcorrdataalltasks(:,1),dcorrdataalltasks(:,2));
                        
                end
        
              	DiceCorrssimmem = [DiceCorrssimmem; dcmem];
                DiceCorrssimmixed = [DiceCorrssimmixed; dcmixed];
                DiceCorrssimmotor = [DiceCorrssimmotor; dcmotor];
                DiceCorrssimalltasks = [DiceCorrssimalltasks; dcalltasks];
                        
              	if splithalf == 1       %% Store values from dice correlations
                            
                  	DiceCorrSubsmem(Count,subiter) = dcmem;
                    DiceCorrSubsmixed(Count,subiter) = dcmixed;
                    DiceCorrSubsmotor(Count,subiter) = dcmotor;
                    DiceCorrSubsalltasks(Count,subiter) = dcalltasks;
                            
                else
                            
                 	DiceCorrSubsmem(Count+(nsubs-1),subiter) = dcmem;
                    DiceCorrSubsmixed(Count+(nsubs-1),subiter) = dcmixed;
                    DiceCorrSubsmotor(Count+(nsubs-1),subiter) = dcmotor;
                    DiceCorrSubsalltasks(Count+(nsubs-1),subiter) = dcalltasks;
                            
                end
                        
                        
            end
        end
    end
end



if Subplots == 1    %% Plot data using subplots
    
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
    
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [.2, .2, .6, .6]);
    
    if AddLegend == 1       %% Add legend if desired
    
        subiter = findobj(gca,'Type','scatter');

        hleg1 = legend([subiter(2:9)], 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC10', 'Location', 'NorthEast');
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



%% Functions for binarizing variant maps

function outmat = binarize_vectors(evenvector,oddvector,dataset)

outmat = [];

for vertices = 1:size(oddvector,2)
    
    if oddvector(dataset,vertices) == 1 && evenvector(dataset,vertices) == 1
        
        outmat = [outmat;1 1];
        
    elseif oddvector(dataset,vertices) == 1
        
        outmat = [outmat;1 0];
        
    elseif evenvector(dataset,vertices) == 1
        
        outmat = [outmat;0 1];
        
    end
end


end


function outmatperm = binarize_vectors_perm(taskvector,restvector)

outmatperm = [];

for vertices = 1:size(taskvector,2)
    
    if restvector(:,vertices) > 0 && taskvector(:,vertices) > 0
        
        outmatperm = [outmatperm;1 1];
        
    elseif restvector(:,vertices) > 0
        
        outmatperm = [outmatperm;1 0];
        
    elseif taskvector(:,vertices) > 0
        
        outmatperm = [outmatperm;0 1];
        
    end
end

end
