
%%% Get center of mass for all variants, then calculate Euclidean distance
%%% between each variant and sort by closest variants

clear all

SizeExclude = 1;    %% Toggles whether to exclude small variants
usemasks = 0;        %% Toggles whether to use masks for excluding variants
minsize = 15;       %% Toggles minimum size for variants if not excluding with maps
SplitHalf = 0;      %% Toggles whether to use split-halves
plothists = 0;      %% Toggles whether to plot histograms of minimum distances
plotmats = 0;       %% Toggles whether to plot matrices of minimum distances
threshold = 2.5;    %% Threshold of variants to use for maps
savemaps = 0;       %% Toggles whether to save out the variant maps after size exclusion
savesurfaces = 0;   %% Toggles whether to save the map of all matched/unmatched variants across subjects
separatetaskandrest = 1;    %% Toggles whether to save out separate maps for task/rest variants
excludeunmatched = 1;   %% Toggles whether to exclude variants that don't match
matchbydist = 0;        %% Toggles whether to match by centroid distance
distexclusion = 10;     %% Toggles threshold for excluding umatched variants by distance (mm)
matchbyoverlap = 1;     %% Toggles whether to match by percent overlap
percentoverlap = 33;    %% Toggles amount of overlap to match on
savestruct = 1;         %% Toggles whether to save a structure of mathed variants

DiceNetworks = 1;       %% Toggles whether to match variants to networks
separatematchedvars = 1;    %% Toggles whether to separate matched/unmatched variants for network assignment
subjectdiceplots = 0;   %% Toggles whether to plot variants for each subject
getmatchbystate = 1;    %% Toggles whether to get matched variants across states
getmaxmatch = 0;    %% Toggles whether to get the max network template match for each variant
excludemaxmatch = 1;    %% Toggles whether to exclude variants with low matches to template

outputdir = '/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/COMHistogramPlots';

cifti_coords = ft_read_cifti_mod('/Users/briankraus/Desktop/Quest_Files/Atlases/32k_ConteAtlas_v2_distribute/CIFTI_MetricSurface_File.dtseries.nii');


if savestruct == 1
    
    Variant_Match_Struct = struct('Subject',{},'Matched',{},'State',{},'Network',{},'Rest_Var_Num',{},'Task_Var_Num',{});
    
end


if DiceNetworks == 1

    network_names = {'DMN'	'Vis'	'FP'	'Unassigned'	'DAN'	'Unassigned2'	'VAN'	'Sal'	'CO'	'SMd'	'SMl'	'Aud'	'Tpole'	'MTL'	'PMN'	'PON'};
    wb_colors = [1 2 3 5 7 8 9 10 11 12 13 14 15 16];
    rgb_colors = [1 0 0;
                  0 0 .6;
                  1 1 0;
                  .67 .67 .67;
                  0 .8 0;
                  .67 .67 .67;
                  0 .6 .6;
                  0 0 0;
                  .3 0 .6;
                  .2 1 1;
                  1 .5 0;
                  .6 .2 1;
                  .2 1 .2;
                  0 .2 .4;
                  0 0 1;
                  .8 .8 .6];

    if getmaxmatch == 1
        
        maxmatchvals = [];

    end
    
    if excludemaxmatch == 1
        
        minmatchval = .1345;    %% Lowest 5 percent of matching values
        
    else
        
        minmatchval = 0;
        
    end
    
    %load templates- matrix of 16 networks x XX,XXX vertices
    nettemplates = ft_read_cifti_mod('/Users/briankraus/Desktop/Quest_Files/120_templates.dtseries.nii');
    nettemplates = nettemplates.data(1:59412,:)';

    %threshold templates at top 5%
    temp = sort(nettemplates(:),'descend');
    templateThresh = temp(round(0.05 * numel(temp)));
    cifti_template_mat_full = nettemplates > templateThresh;
    clear temp
    
    %%% group-average network assignments (surface only), 59412 x 1
    consensus = ft_read_cifti_mod('/Users/briankraus/Desktop/Quest_Files/120_colorassn_minsize400_manualconsensus.dtseries.nii');
    
    if SplitHalf == 1
        
        [task_timeseries_even,~,~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_timeseries_matched_rest.txt','%s%s%s');

        [rest_timeseries_even,~,~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_timeseries_matched_task.txt','%s%s%s');
        
        [task_timeseries_odd,~,~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_timeseries_matched_rest.txt','%s%s%s');

        [rest_timeseries_odd,~,~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_timeseries_matched_task.txt','%s%s%s');
        
    else
    
        [task_timeseries,~,~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_timeseries_matched_rest.txt','%s%s%s');

        [rest_timeseries,~,~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_timeseries_matched_task.txt','%s%s%s');
        
    end
    
    if separatematchedvars == 1
        
        restmatchednetnums = zeros(1,16);
     	restunmatchednetnums = zeros(1,16);
     	taskmatchednetnums = zeros(1,16);
     	taskunmatchednetnums = zeros(1,16);
        
        restunmatchednetnumsAllSubs = zeros(size(rest_timeseries,1),16);
        restmatchednetnumsAllSubs = zeros(size(rest_timeseries,1),16);
        taskunmatchednetnumsAllSubs = zeros(size(rest_timeseries,1),16);
        taskmatchednetnumsAllSubs = zeros(size(rest_timeseries,1),16);
        
        if getmatchbystate == 1
        
            matchedvarsnetsrest = {};
            matchedvarsnetstask = {};
            
        end
        
    else
        
        restnetnums = zeros(1,16);
        tasknetnums = zeros(1,16);
        
        restnetnumsAllSubs = zeros(size(rest_timeseries,1),16);
        tasknetnumsAllSubs = zeros(size(rest_timeseries,1),16);
        
    end
    
end

if SizeExclude == 1 && SplitHalf == 1
    
	[task_files_even, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_splithalf_even_matched_variants_SizeOnly.txt','%s%s%s');

	[rest_files_even, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_splithalf_even_matched_variants_SizeOnly.txt','%s%s%s');
    
	[task_files_odd, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_splithalf_odd_matched_variants_SizeOnly.txt','%s%s%s');

	[rest_files_odd, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_splithalf_odd_matched_variants_SizeOnly.txt','%s%s%s');
    
	[task_masks_even, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_splithalf_even_matched_variants_SizeOnly_sizeexclude.txt','%s%s%s');

	[rest_masks_even, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_splithalf_even_matched_variants_SizeOnly_sizeexclude.txt','%s%s%s');
    
	[task_masks_odd, sub1, t1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_splithalf_odd_matched_variants_SizeOnly_sizeexclude.txt','%s%s%s');

	[rest_masks_odd, sub2, t2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_splithalf_odd_matched_variants_SizeOnly_sizeexclude.txt','%s%s%s');
    
elseif SizeExclude == 1
    
	[task_files, subjects1, ~] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_' num2str(threshold) '_matched_variants_SNRexclude.txt'],'%s%s%s');

	[rest_files, subjects2, ~] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_' num2str(threshold) '_matched_variants_SNRexclude.txt'],'%s%s%s');
    
	[task_masks, sub1, t1] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_' num2str(threshold) '_matched_variants_SNRexclude_sizeexclude.txt'],'%s%s%s');

	[rest_masks, sub2, t2] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_' num2str(threshold) '_matched_variants_SNRexclude_sizeexclude.txt'],'%s%s%s');

elseif SplitHalf == 1
    
 	[task_files_even, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_splithalf_even_matched_variants_SNRexclude.txt','%s%s%s');

   	[rest_files_even, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_splithalf_even_matched_variants_SNRexclude.txt','%s%s%s');
    
  	[task_files_odd, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_splithalf_odd_matched_variants_SNRexclude.txt','%s%s%s');

  	[rest_files_odd, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_splithalf_odd_matched_variants_SNRexclude.txt','%s%s%s');
    
else
    
	[task_files, subjects1, ~] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_' num2str(threshold) '_variants_SNRexclude_matched.txt'],'%s%s%s');

	[rest_files, subjects2, ~] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_' num2str(threshold) '_variants_SNRexclude_matched.txt'],'%s%s%s');
    
end


if SplitHalf == 1
        
    COMOverlapsTaskRest = [];
  	COMOverlapsTaskTask = [];
  	COMOverlapsRestRest = [];
        
else
    
 	COMOverlaps = [];
        
end

if SplitHalf == 1
    
	nfiles = length(rest_files_even);
    
else
    
 	nfiles = length(rest_files);
    
end

COMDistanceMatrixAll = {};

for x = 1:nfiles
    
    if separatematchedvars == 1
        
        restmatchednetnumsSub = zeros(1,16);
     	restunmatchednetnumsSub = zeros(1,16);
     	taskmatchednetnumsSub = zeros(1,16);
     	taskunmatchednetnumsSub = zeros(1,16);
        
    else
        
        restnetnumsSub = zeros(1,16);
        tasknetnumsSub = zeros(1,16);
        
    end

	subject = subjects2{x};
    
    disp(['Processing data for subject ' subject])
    
  	if SplitHalf == 1
        
    	cifti_rest_even = ft_read_cifti_mod(rest_files_even{x});
      	cifti_task_even = ft_read_cifti_mod(task_files_even{x});
       	cifti_rest_odd = ft_read_cifti_mod(rest_files_odd{x});
       	cifti_task_odd = ft_read_cifti_mod(task_files_odd{x});
        
    else
        
      	cifti_rest = ft_read_cifti_mod(rest_files{x});
      	cifti_task = ft_read_cifti_mod(task_files{x});
        
    end
    
    if separatetaskandrest == 1 && separatematchedvars == 1 && x == 1 %&& savesurfaces == 1 
        
        unmatchedvarsrest = zeros(size(cifti_rest.data));
        matchedvarsrest = zeros(size(cifti_rest.data));
        unmatchedvarstask = zeros(size(cifti_rest.data));
        matchedvarstask = zeros(size(cifti_rest.data));
    
    elseif separatematchedvars == 1 && x == 1 %&& savesurfaces == 1
    
        unmatchedvars = zeros(size(cifti_rest.data));
        matchedvars = zeros(size(cifti_rest.data));
        
    else
        
        alltaskvars = zeros(size(cifti_rest.data));
        allrestvars = zeros(size(cifti_rest.data));
    
    end
    
   	SNRmask = ft_read_cifti_mod(['/Users/briankraus/Desktop/MSC_SNR_Maps/CIFTI_Files/MSC_Template/Individual_Session_Map/' subject '/' subject '__SNRMap__AllDataConcatenated.dscalar.nii']); 
    
    SNRmask.data = SNRmask.data(1:59412,:);
	LowSNR = find(SNRmask.data < 750);
    
    if SplitHalf == 1
            
        cifti_rest_even.data(LowSNR,:) = 0;
        cifti_task_even.data(LowSNR,:) = 0;
        cifti_rest_odd.data(LowSNR,:) = 0;
        cifti_task_odd.data(LowSNR,:) = 0;
        
    else

        cifti_rest.data(LowSNR,:) = 0;
        cifti_task.data(LowSNR,:) = 0;
        
    end
        
        
  	if SizeExclude == 1         %% Apply exclusion masks for maps that exclude size
        
        if usemasks == 1
        
            if SplitHalf == 1
            
                cifti_rest_mask_even = ft_read_cifti_mod(rest_masks_even{x});
                cifti_task_mask_even = ft_read_cifti_mod(task_masks_even{x});
                cifti_rest_mask_odd = ft_read_cifti_mod(rest_masks_odd{x});
                cifti_task_mask_odd = ft_read_cifti_mod(task_masks_odd{x});
        
                for d = 1:length(cifti_rest_mask_even.data)
            
                    if cifti_rest_mask_even.data(d) == 0
                
                        cifti_rest_even.data(d) = 0;
                
                    end
                
                    if cifti_rest_mask_odd.data(d) == 0
                    
                        cifti_rest_odd.data(d) = 0;
                    
                    end  
                end
        
                for e = 1:length(cifti_task_mask_even.data)
            
                    if cifti_task_mask_even.data(e) == 0
                
                        cifti_task_even.data(e) = 0;
                
                    end
                
                    if cifti_task_mask_odd.data(e) == 0
                    
                        cifti_task_odd.data(e) = 0;
                    
                    end
                end
            
            else
        
                cifti_rest_mask = ft_read_cifti_mod(rest_masks{x});
                cifti_task_mask = ft_read_cifti_mod(task_masks{x});
        
                for d = 1:length(cifti_rest_mask.data)
            
                    if cifti_rest_mask.data(d) == 0
                
                        cifti_rest.data(d) = 0;
                
                    end
                end
        
                for e = 1:length(cifti_task_mask.data)
            
                    if cifti_task_mask.data(e) == 0
                
                        cifti_task.data(e) = 0;
                
                    end
                end
            end
            
        else
            
            if SplitHalf == 1
                
                allvars_rest_even = unique(cifti_rest_even.data);
                allvars_rest_even(1) = [];
                allvars_rest_odd = unique(cifti_rest_odd.data);
                allvars_rest_odd(1) = [];
                allvars_task_even = unique(cifti_task_even.data);
                allvars_task_even(1) = [];
                allvars_task_odd = unique(cifti_task_odd.data);
                allvars_task_odd(1) = [];
                
                removevars_rest_even = [];
                removevars_rest_odd = [];
                removevars_task_even = [];
                removevars_task_odd = [];
                
                for g = 1:length(allvars_rest_even)
                    
                    if length(find(cifti_rest_even.data == allvars_rest_even(g))) < minsize
                        
                        removevars_rest_even = [removevars_rest_even allvars_rest_even(g)];
                        
                    end
                end
                
                for h = 1:length(cifti_rest_even.data)
                    
                    if ismember(cifti_rest_even.data(h),removevars_rest_even)
                        
                        cifti_rest_even.data(h) = 0;
                        
                    end
                end
                
                for g = 1:length(allvars_rest_odd)
                    
                    if length(find(cifti_rest_odd.data == allvars_rest_odd(g))) < minsize
                        
                        removevars_rest_odd = [removevars_rest_odd allvars_rest_odd(g)];
                        
                    end
                end
                
                for h = 1:length(cifti_rest_odd.data)
                    
                    if ismember(cifti_rest_odd.data(h),removevars_rest_odd)
                        
                        cifti_rest_odd.data(h) = 0;
                        
                    end
                end
                
                for g = 1:length(allvars_task_even)
                    
                    if length(find(cifti_task_even.data == allvars_task_even(g))) < minsize
                        
                        removevars_task_even = [removevars_task_even allvars_task_even(g)];
                        
                    end
                end
                
                for h = 1:length(cifti_task_even.data)
                    
                    if ismember(cifti_task_even.data(h),removevars_task_even)
                        
                        cifti_task_even.data(h) = 0;
                        
                    end
                end
                
                for g = 1:length(allvars_task_odd)
                    
                    if length(find(cifti_task_odd.data == allvars_task_odd(g))) < minsize
                        
                        removevars_task_odd = [removevars_task_odd allvars_task_odd(g)];
                        
                    end
                end
                
                for h = 1:length(cifti_task_odd.data)
                    
                    if ismember(cifti_task_odd.data(h),removevars_task_odd)
                        
                        cifti_task_odd.data(h) = 0;
                        
                    end
                end
                
            else
                
            	allvars_rest = unique(cifti_rest.data);
                allvars_rest(1) = [];
                allvars_task = unique(cifti_task.data);
                allvars_task(1) = [];
                
                removevars_rest = [];
                removevars_task = [];
                
                for g = 1:length(allvars_rest)
                    
                    if length(find(cifti_rest.data == allvars_rest(g))) < minsize
                        
                        removevars_rest = [removevars_rest allvars_rest(g)];
                        
                    end
                end
                
                for h = 1:length(cifti_rest.data)
                    
                    if ismember(cifti_rest.data(h),removevars_rest)
                        
                        cifti_rest.data(h) = 0;
                        
                    end
                end
                
                for g = 1:length(allvars_task)
                    
                    if length(find(cifti_task.data == allvars_task(g))) < minsize
                        
                        removevars_task = [removevars_task allvars_task(g)];
                        
                    end
                end
                
                for h = 1:length(cifti_task.data)
                    
                    if ismember(cifti_task.data(h),removevars_task)
                        
                        cifti_task.data(h) = 0;
                        
                    end
                end
                
            end
        end
    end
    
    if savemaps == 1
        
     	template = ft_read_cifti_mod('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Rest_Data/Full_Rest/Full_Data/MSC01/MSC01_REST_AllSessions_vs_120_allsubs_corr_cortex_corr.dtseries.nii');
       	template.data = [];
        template2 = template;
        
        template.data = cifti_rest.data;
        template2.data = cifti_task.data;
        
      	ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Rest_Data/TaskCat_Matched_Rest/' subject '/' subject '_Rest_Thresholded_BothSNRAndSizeExcluded_Map_' num2str(threshold) '_pct'],template)
     	clear template

      	ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Task_Data/TaskCat_Matched_Data/' subject '/' subject '_Task_Thresholded_BothSNRAndSizeExcluded_Map_' num2str(threshold) '_pct'],template2)
     	clear template2
        
    end
    
	if SplitHalf == 1       %% Get number of variants
        
     	vars_rest_even = unique(cifti_rest_even.data);
    	vars_rest_even(1) = [];
      	vars_rest_odd = unique(cifti_rest_odd.data);
        vars_rest_odd(1) = [];
    	vars_task_even = unique(cifti_task_even.data);
    	vars_task_even(1) = [];
     	vars_task_odd = unique(cifti_task_odd.data);
     	vars_task_odd(1) = [];
        
    else
        
      	vars_rest = unique(cifti_rest.data);
       	vars_rest(1) = [];
      	vars_task = unique(cifti_task.data);
       	vars_task(1) = [];
        
    end
    
   	if SplitHalf == 1       %% Generate COMs for each file
            
     	for z = 1:3     %% Loops over 3 comparisons (rest-task, rest-rest, task-task)
            
         	%meanvals = [];
          	COMOverlapsresttemp = [];
          	COMOverlapstasktemp = [];
                
         	if z == 1  %% Task-Rest
        
              	taskverts = find(cifti_task_odd.data > 0);
              	restverts = find(cifti_rest_even.data > 0);
                    
          	elseif z == 2  %% Task-Task
                    
              	taskverts = find(cifti_task_even.data > 0);
               	restverts = find(cifti_task_odd.data > 0);
                    
           	elseif z == 3  %% Rest-Rest
                    
               	taskverts = find(cifti_rest_even.data > 0);
               	restverts = find(cifti_rest_odd.data > 0);
                    
            end
        
         	for a = 1:2   %% loop over task and rest
            
             	if a == 1   %% Rest
                        
                	if z == 1   %% Task-Rest
                
                    	nvars = length(vars_rest_even);
                        
                  	elseif z == 2   %% Task-Task
                            
                      	nvars = length(vars_task_odd);
                            
                 	elseif z == 3   %% Rest-Rest
                            
                      	nvars = length(vars_rest_odd);
                            
                 	end
                
               	else       %% Task
                
                   	if z == 1   %% Task-Rest
                
                     	nvars = length(vars_task_odd);
                        
                   	elseif z == 2   %% Task-Task
                            
                       	nvars = length(vars_task_even);
                            
                   	elseif z == 3   %% Rest-Rest
                            
                      	nvars = length(vars_rest_even);
                            
                    end
                
                end
        
                for s = 1:nvars
                
                 	if a == 1   %% Rest
                            
                        if z == 1   %% Task-Rest
                
                          	currentvariant = find(cifti_rest_even.data == vars_rest_even(s));
                        
                       	elseif z == 2   %% Task-Task

                          	currentvariant = find(cifti_task_odd.data == vars_task_odd(s));
                            
                     	elseif z == 3   %% Rest-Rest

                           	currentvariant = find(cifti_rest_odd.data == vars_rest_odd(s));
                            
                        end
                
                    else   %% Task
                    
                       	if z == 1   %% Task-Rest
                
                          	currentvariant = find(cifti_task_odd.data == vars_task_odd(s));
                        
                       	elseif z == 2   %% Task-Task

                           	currentvariant = find(cifti_task_even.data == vars_task_even(s));
                            
                     	elseif z == 3   %% Rest-Rest

                           	currentvariant = find(cifti_rest_even.data == vars_rest_even(s));
                            
                        end
                    
                    end
                
                    COMtemp = cifti_coords.data(currentvariant,:);
                
                  	if length(currentvariant) > 1
                    
                     	centroid = [];
                    
                      	for t = 1:length(currentvariant)
                        
                         	jacknifevals = 1:length(currentvariant);        %% Leaves current point out of calculation
                           	jacknifevals(t) = [];
                    
                           	% Euclidean distance
                           	dist = mean(sqrt(sum((COMtemp(jacknifevals,:) - COMtemp(t,:)).^2, 2)));
                        
                          	%meanvals = [meanvals; dist];
                        
                         	if isempty(centroid) || dist < centroid   %% If new distance less than current least distance, use as center
                            
                              	centroid = currentvariant(t);
                            
                            end
                        
                        end
                    
                       	centervariant = centroid;
                    
                    else
                    
                     	centervariant = currentvariant;
                    
                    end

                	if a == 1
                
                     	COMOverlapsresttemp = [COMOverlapsresttemp; centervariant];
                    
                    else
                    
                     	COMOverlapstasktemp = [COMOverlapstasktemp; centervariant];
                    
                    end 
                end
            end

        
         	%Taskoverlaps = length(unique(intersect(COMOverlapstasktemp, restverts)));
         	%Restoverlaps = length(unique(intersect(COMOverlapsresttemp, taskverts)));
                
        	if z == 1   %% Task-Rest
        
                %COMOverlapsTaskRest = [COMOverlapsTaskRest; (Taskoverlaps + Restoverlaps)/(length(vars_rest_even) + length(vars_task_odd))];
                    
          	elseif z == 2  %% Task-Task
                    
              	%COMOverlapsTaskTask = [COMOverlapsTaskTask; (Taskoverlaps + Restoverlaps)/(length(vars_task_even) + length(vars_task_odd))];
                    
         	elseif z == 3  %% Rest-Rest
                    
              	%COMOverlapsRestRest = [COMOverlapsRestRest; (Taskoverlaps + Restoverlaps)/(length(vars_rest_even) + length(vars_rest_odd))];
                    
            end
                
        end
            
    else
    
       	%meanvals = [];
       	COMOverlapsresttemp = [];
       	COMOverlapstasktemp = [];

       	taskverts = find(cifti_task.data > 0);
      	restverts = find(cifti_rest.data > 0);
        
       	for a = 1:2   %% loop over task and rest
            
         	if a == 1
                
             	nvars = length(vars_rest);
                
            else
                
             	nvars = length(vars_task);
                
            end
        
          	for s = 1:nvars
                
              	if a == 1
            
                	currentvariant = find(cifti_rest.data == vars_rest(s));
                
                else
                    
                   	currentvariant = find(cifti_task.data == vars_task(s));
                    
                end
                
                 	COMtemp = cifti_coords.data(currentvariant,:);
                
             	if length(currentvariant) > 1
                    
                   	centroid = [];
                    
                  	for t = 1:length(currentvariant)
                        
                     	jacknifevals = 1:length(currentvariant);        %% Leaves current point out of calculation
                      	jacknifevals(t) = [];
                    
                       	% Euclidean distance
                       	dist = mean(sqrt(sum((COMtemp(jacknifevals,:) - COMtemp(t,:)).^2, 2)));
                        
                     	%meanvals = [meanvals; dist];
                        
                      	if isempty(centroid) || dist < centroid   %% If new distance less than current least distance, use as center
                            
                          	centroid = currentvariant(t);
                            
                        end
                        
                    end
                    
                    centervariant = centroid;
                    
                else
                    
                  	centervariant = currentvariant;
                    
                end
                
            	if a == 1
                
                  	COMOverlapsresttemp = [COMOverlapsresttemp; centervariant];
                    
                else
                    
                   	COMOverlapstasktemp = [COMOverlapstasktemp; centervariant];
                    
                end 
            end
        end

        
    	%Taskoverlaps = length(unique(intersect(COMOverlapstasktemp, restverts)));
     	%Restoverlaps = length(unique(intersect(COMOverlapsresttemp, taskverts)));
        
     	%COMOverlaps = [COMOverlaps; (Taskoverlaps + Restoverlaps)/(length(vars_rest) + length(vars_task))];
        
    end
    
    
   	COMCoordsresttemp = [];     %% Get coordinates for each centroid
   	COMCoordstasktemp = [];
        
   	for a = 1:length(COMOverlapsresttemp)
            
     	COMCoordsresttemp = [COMCoordsresttemp; cifti_coords.data(COMOverlapsresttemp(a),:)];
        
    end
    
   	for b = 1:length(COMOverlapstasktemp)
            
     	COMCoordstasktemp = [COMCoordstasktemp; cifti_coords.data(COMOverlapstasktemp(b),:)];
        
    end
    
    COMDistanceMatrix = zeros(size(COMCoordsresttemp,1),size(COMCoordstasktemp,1));     %% Calculate distance between each variant in task and rest
    
    for f = 1:size(COMCoordsresttemp,1)
        
        for g = 1:size(COMCoordstasktemp,1)
            
            COMDistanceMatrix(f,g) = sqrt(sum((COMCoordsresttemp(f,:) - COMCoordstasktemp(g,:)).^2));
            
        end
    end
    

    if matchbydist == 1
    
        %[M,I] = min(COMDistanceMatrix,[],2);
    

        COMCoordsresttemp1 = COMCoordsresttemp;
        COMCoordstasktemp1 = COMCoordstasktemp;
    
     
        %%% First set of values must be smaller than the second set in call
        %%% to command
   
        
        if size(COMCoordsresttemp,1) < size(COMCoordstasktemp,1)
            
            [TaskVariantMatches, COMKMatrix] = knnsearch(COMCoordsresttemp1,COMCoordstasktemp1, 'K', size(COMCoordsresttemp,1));      %% Use nearest neighbors algorithm to find all closest distances
            
        else

            [TaskVariantMatches, COMKMatrix] = knnsearch(COMCoordstasktemp1,COMCoordsresttemp1, 'K', size(COMCoordstasktemp1,1));      %% Use nearest neighbors algorithm to find all closest distances
            
        end


            COMDistanceMatrixAll = [COMDistanceMatrixAll; COMKMatrix(:,1)];
        
            COMTaskVarsMatched = [[1:size(TaskVariantMatches,1)]' TaskVariantMatches(:,1) COMKMatrix(:,1)];
        
            COMTaskVarsMatched = sortrows(COMTaskVarsMatched, [2 3]);
        
            if excludeunmatched == 1
    
                %mindistmm = prctile(COMKMatrix(:,1),distexclusion);
                mindistmm = distexclusion;
            
                removedvars = find(COMTaskVarsMatched(:,3) > mindistmm);
            
                COMTaskVarsMatchedOrig = COMTaskVarsMatched;
            
                COMTaskVarsMatched(removedvars,:) = [];
            
                finalremovedvars = setxor(COMTaskVarsMatchedOrig(:,2),  COMTaskVarsMatched(:,2));
            
                %COMDistanceMatrix2 = COMDistanceMatrix;
                %COMDistanceMatrix2(find(COMTaskVarsMatched(:,3) > mindistmm),:) = []; 
            
                if size(COMCoordsresttemp,1) < size(COMCoordstasktemp,1)        %% Redo with leftover variants to ensure that all possible matches are found
                
                    COMCoordsresttemp2 = COMCoordsresttemp1(finalremovedvars,:);
                    %COMCoordsresttemp2(finalremovedvars,:) = []; 
                
                    [TaskVariantMatches2, COMKMatrix2] = knnsearch(COMCoordstasktemp1,COMCoordsresttemp2, 'K', size(COMCoordsresttemp,1));      %% Use nearest neighbors algorithm to find all closest distances
            
                else
                
                    COMCoordstasktemp2 = COMCoordstasktemp1(finalremovedvars,:);
                    %COMCoordstasktemp2(COMTaskVarsMatched(:,2),:) = []; 
                
                    [TaskVariantMatches2, COMKMatrix2] = knnsearch(COMCoordsresttemp1,COMCoordstasktemp2, 'K', size(COMCoordstasktemp1,1));      %% Use nearest neighbors algorithm to find all closest distances
            
                end
            
                if length(find(COMKMatrix2(:,1) < mindistmm)) > 0
                
                    addvarnums = find(COMKMatrix2(:,1) < mindistmm);
                
                    for r = 1:length(addvarnums)
                
                        COMTaskVarsMatched = [COMTaskVarsMatched; [TaskVariantMatches2(addvarnums(r),1) finalremovedvars(addvarnums(r)) COMKMatrix2(addvarnums(r),1)]];

                        %COMDistanceMatrix2 = [COMDistanceMatrix2; COMDistanceMatrix(TaskVariantMatches2(addvarnums(r),1),:)];
                    
                    end
                
                    COMTaskVarsMatched = sortrows(COMTaskVarsMatched, [2 3]);
                
                end
            
                if size(COMCoordsresttemp,1) < size(COMCoordstasktemp,1)
                
                    FinalMatchedVars = [vars_rest(TaskVariantMatches(COMTaskVarsMatched(:,1),1),:) vars_task(COMTaskVarsMatched(:,1),:)];
                
                else
                
                    FinalMatchedVars = [vars_rest(COMTaskVarsMatched(:,1),:) vars_task(TaskVariantMatches(COMTaskVarsMatched(:,1),1),:)];
                
                end
    
            end
            
    else
        
        COMTaskVarsMatched = [];
        FinalMatchedVars = [];
        
    end
        
        if matchbyoverlap == 1
            
            overlapmatchesfinal = [];
            overlapmatches = [];
            
            for q = 1:length(vars_rest)
                
                for w = 1:length(vars_task)
                    
                    restpts = find(cifti_rest.data == vars_rest(q));
                    taskpts = find(cifti_task.data == vars_task(w));
                    
                    if length(intersect(restpts,taskpts)) >= (length(restpts)*(percentoverlap/100)) || length(intersect(restpts,taskpts)) >= (length(taskpts)*(percentoverlap/100))
                    
                        overlapmatchesfinal = [overlapmatchesfinal; vars_rest(q) vars_task(w)];
                        
                        if size(COMCoordsresttemp,1) < size(COMCoordstasktemp,1)
                            
                            overlapmatches = [overlapmatches; w q];
                            
                        else
                        
                            overlapmatches = [overlapmatches; q w];
                            
                        end
                    end
                end
            end
            
            if matchbydist == 1
            
                addmatches = setxor(COMTaskVarsMatched(:,1:2),overlapmatches,'rows');
                
                COMTaskVarsMatched = [COMTaskVarsMatched; [addmatches zeros(size(addmatches,1),1)]];
                
                addmatchesfinal = setxor(FinalMatchedVars,overlapmatchesfinal,'rows');
                
            else
                
                addmatchesfinal = overlapmatchesfinal;
                
            end
            
          	FinalMatchedVars = [FinalMatchedVars; addmatchesfinal];
            
        end
                
        %if savesurfaces == 1        %% Save subject surface maps out by border shifts/new variants
            
            if DiceNetworks == 1
                
                
                cifti_rest_timeseries = ft_read_cifti_mod(rest_timeseries{x});
                cifti_task_timeseries = ft_read_cifti_mod(task_timeseries{x});
                
                    if separatematchedvars == 1
                
                      	varIDsreassign_restmatched = ft_read_cifti_mod('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Rest_Data/Full_Rest/Full_Data/MSC01/MSC01_REST_AllSessions_vs_120_allsubs_corr_cortex_corr.dtseries.nii');
                        varIDsreassign_restmatched.data = zeros(size(cifti_rest.data));
                      	varIDsreassign_restunmatched = ft_read_cifti_mod('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Rest_Data/Full_Rest/Full_Data/MSC01/MSC01_REST_AllSessions_vs_120_allsubs_corr_cortex_corr.dtseries.nii');
                        varIDsreassign_restunmatched.data = zeros(size(cifti_rest.data));
                        
                        %%% loop through variants, get avg seedmap, match to network templates %%%
                        for var = 1:length(vars_rest)
                        
                            inds = find(cifti_rest.data == vars_rest(var));        
                            varRmat = paircorr_mod(cifti_rest_timeseries.data(inds,:)',cifti_rest_timeseries.data'); %correlate each vertex within variant to all vertices in brain
                            varmean = mean(varRmat,1); %average seedmap of variant
                            %fprintf([num2str(var) '  '])
        
                            temp = sort(varmean,'descend');
                            varthresh = temp(round(0.05 * numel(temp)));
                            varmean_thresh = varmean > varthresh;

                            for net = 1:14
                            
                                %%% dice calculation b/w variant & each network %%%            
                                netTemplate = cifti_template_mat_full(net,:);
                                var_dice_to_template(net,:) = (sum((netTemplate & varmean_thresh),2)*2) ./ sum((netTemplate | varmean_thresh),2);
            
                            end
            
                            %%% pick highest eta value and assign variant to that network %%%
                            [maxval, winner] = max(var_dice_to_template,[],1);
                            
                        	if getmaxmatch == 1
                                    
                             	maxmatchvals = [maxmatchvals; maxval];
                                    
                         	end
                            
                            if maxval > 0 && (excludemaxmatch == 1 && maxval > minmatchval)
                            
                                winner = wb_colors(winner);
                                
                            else
                                
                                winner = 4;
                                
                            end
                            
                            if ismember(vars_rest(var),FinalMatchedVars(:,1))
                            
                                varIDsreassign_restmatched.data(inds) = winner;
                                
                                %diceOfWinner.data(inds)=x;

                                %%% remove variant if over half its vertices overlap w consensus %%%
%                                 if length(find(ismember(inds,find(consensus.data==winner)))) > (length(inds)*0.5)
%                             
%                                     varIDsreassign_restmatched.data(inds) = 0;
%                                     %diceOfWinner.data(inds) = 0;
% 
%                                 else
                                    
                                    restmatchednetnums(:,winner) = restmatchednetnums(:,winner) + 1;
                                    restmatchednetnumsSub(:,winner) = restmatchednetnumsSub(:,winner) + 1;
                            
%                                 end

                              	if savestruct == 1
                                    
                                    StructAdd = size(Variant_Match_Struct,2)+1;
                                    
                                    Variant_Match_Struct(StructAdd).Subject = subject;
                                    Variant_Match_Struct(StructAdd).Matched = 1;
                                    Variant_Match_Struct(StructAdd).State = 'Rest';
                                    Variant_Match_Struct(StructAdd).Network = network_names(winner);
                                    Variant_Match_Struct(StructAdd).Rest_Var_Num = vars_rest(var);
                                    Variant_Match_Struct(StructAdd).Task_Var_Num = overlapmatchesfinal(find(vars_rest(var) == overlapmatchesfinal(:,1)),2);
                                    
                                end
                                
                            else
                                
                                varIDsreassign_restunmatched.data(inds) = winner;
                                
                                %diceOfWinner.data(inds)=x;

                                %%% remove variant if over half its vertices overlap w consensus %%%
%                                 if length(find(ismember(inds,find(consensus.data==winner)))) > (length(inds)*0.5)
%                             
%                                     varIDsreassign_restunmatched.data(inds) = 0;
%                                     %diceOfWinner.data(inds) = 0;
%                                     
%                                 else
                                    
                                    restunmatchednetnums(:,winner) = restunmatchednetnums(:,winner) + 1;
                                    restunmatchednetnumsSub(:,winner) = restunmatchednetnumsSub(:,winner) + 1;
                            
%                                 end
                              	if savestruct == 1
                                    
                                    StructAdd = size(Variant_Match_Struct,2)+1;
                                    
                                    Variant_Match_Struct(StructAdd).Subject = subject;
                                    Variant_Match_Struct(StructAdd).Matched = 0;
                                    Variant_Match_Struct(StructAdd).State = 'Rest';
                                    Variant_Match_Struct(StructAdd).Network = network_names(winner);
                                    Variant_Match_Struct(StructAdd).Rest_Var_Num = vars_rest(var);
                                    Variant_Match_Struct(StructAdd).Task_Var_Num = 0;
                                    
                                end

                            end
        
                        end
                        
                        if getmatchbystate == 1
                            
                            for v = 1:size(FinalMatchedVars,1)
                                
                                inds = find(cifti_rest.data == FinalMatchedVars(v,1));        
                                varRmat = paircorr_mod(cifti_rest_timeseries.data(inds,:)',cifti_rest_timeseries.data'); %correlate each vertex within variant to all vertices in brain
                                varmean = mean(varRmat,1); %average seedmap of variant
                                %fprintf([num2str(var) '  '])
        
                                temp = sort(varmean,'descend');
                                varthresh = temp(round(0.05 * numel(temp)));
                                varmean_thresh = varmean > varthresh;

                                for net = 1:14
                            
                                    %%% dice calculation b/w variant & each network %%%            
                                    netTemplate = cifti_template_mat_full(net,:);
                                    var_dice_to_template(net,:) = (sum((netTemplate & varmean_thresh),2)*2) ./ sum((netTemplate | varmean_thresh),2);
            
                                end
            
                                %%% pick highest eta value and assign variant to that network %%%
                                [maxval, winner] = max(var_dice_to_template,[],1);
                                
                                if getmaxmatch == 1
                                    
                                    maxmatchvals = [maxmatchvals; maxval];
                                    
                                end
                            
                                if maxval > 0 && (excludemaxmatch == 1 && maxval > minmatchval)
                            
                                    winner = wb_colors(winner);
                                
                                else
                                
                                    winner = 4;
                                
                                end
                            
                                matchedvarsnetsrest = [matchedvarsnetsrest; network_names(winner)];
                                
                                %diceOfWinner.data(inds)=x;

                                %%% remove variant if over half its vertices overlap w consensus %%%
%                                 if length(find(ismember(inds,find(consensus.data==winner)))) > (length(inds)*0.5)
%                             
%                                     matchedvarsnetsrest{end} = ' ';
%                                     
%                                 end

                            end
                        end
                        
                        if SizeExclude == 1
                            
                            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Variant_Network_Matching/' subject '_MatchedVariants_Reassigned_Rest_' num2str(threshold) '_pct_SizeExclude'],varIDsreassign_restmatched)
                            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Variant_Network_Matching/' subject '_UnmatchedVariants_Reassigned_Rest_' num2str(threshold) '_pct_SizeExclude'],varIDsreassign_restunmatched)
                            
                        else
                        
                            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Variant_Network_Matching/' subject '_MatchedVariants_Reassigned_Rest_' num2str(threshold) '_pct'],varIDsreassign_restmatched)
                            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Variant_Network_Matching/' subject '_UnmatchedVariants_Reassigned_Rest_' num2str(threshold) '_pct'],varIDsreassign_restunmatched)
                            
                        end
                        
                        clear varIDsreassign_restmatched
                        clear varIDsreassign_restunmatched
                        
                      	varIDsreassign_taskmatched = ft_read_cifti_mod('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Rest_Data/Full_Rest/Full_Data/MSC01/MSC01_REST_AllSessions_vs_120_allsubs_corr_cortex_corr.dtseries.nii');
                        varIDsreassign_taskmatched.data = zeros(size(cifti_task.data));
                      	varIDsreassign_taskunmatched = ft_read_cifti_mod('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Rest_Data/Full_Rest/Full_Data/MSC01/MSC01_REST_AllSessions_vs_120_allsubs_corr_cortex_corr.dtseries.nii');
                        varIDsreassign_taskunmatched.data = zeros(size(cifti_task.data));
                        
                        %%% loop through variants, get avg seedmap, match to network templates %%%
                        for var = 1:length(vars_task)
                        
                            inds = find(cifti_task.data == vars_task(var));        
                            varRmat = paircorr_mod(cifti_task_timeseries.data(inds,:)',cifti_task_timeseries.data'); %correlate each vertex within variant to all vertices in brain
                            varmean = mean(varRmat,1); %average seedmap of variant
                            %fprintf([num2str(var) '  '])
        
                            temp = sort(varmean,'descend');
                            varthresh = temp(round(0.05 * numel(temp)));
                            varmean_thresh = varmean > varthresh;

                            for net = 1:14
                            
                                %%% dice calculation b/w variant & each network %%%            
                                netTemplate = cifti_template_mat_full(net,:);
                                var_dice_to_template(net,:) = (sum((netTemplate & varmean_thresh),2)*2) ./ sum((netTemplate | varmean_thresh),2);
            
                            end
            
                            %%% pick highest eta value and assign variant to that network %%%
                            [maxval, winner] = max(var_dice_to_template,[],1);
                            
                        	if getmaxmatch == 1
                                    
                             	maxmatchvals = [maxmatchvals; maxval];
                                    
                         	end
                            
                            if maxval > 0 && (excludemaxmatch == 1 && maxval > minmatchval)
                            
                                winner = wb_colors(winner);
       
                            else
                                
                                winner = 4;
                                
                            end
                            
                            if ismember(vars_task(var),FinalMatchedVars(:,2))
                            
                                varIDsreassign_taskmatched.data(inds) = winner;
                                %diceOfWinner.data(inds)=x;

                                %%% remove variant if over half its vertices overlap w consensus %%%
%                                 if length(find(ismember(inds,find(consensus.data==winner)))) > (length(inds)*0.5)
%                             
%                                     varIDsreassign_taskmatched.data(inds) = 0;
%                                     %diceOfWinner.data(inds) = 0;
% 
%                                 else
                                    
                                    taskmatchednetnums(:,winner) = taskmatchednetnums(:,winner) + 1;
                                    taskmatchednetnumsSub(:,winner) = taskmatchednetnumsSub(:,winner) + 1;
                            
%                                 end

                                if savestruct == 1

                                  	StructAdd = size(Variant_Match_Struct,2)+1;

                                	Variant_Match_Struct(StructAdd).Subject = subject;
                                 	Variant_Match_Struct(StructAdd).Matched = 1;
                                    Variant_Match_Struct(StructAdd).State = 'Task';
                                  	Variant_Match_Struct(StructAdd).Network = network_names(winner);
                                 	Variant_Match_Struct(StructAdd).Rest_Var_Num = overlapmatchesfinal(find(vars_task(var) == overlapmatchesfinal(:,2)),1);
                                  	Variant_Match_Struct(StructAdd).Task_Var_Num = vars_task(var);
                                        
                                end
                                
                            else
                                
                                varIDsreassign_taskunmatched.data(inds) = winner;
                                
                                %diceOfWinner.data(inds)=x;

                                %%% remove variant if over half its vertices overlap w consensus %%%
%                                 if length(find(ismember(inds,find(consensus.data==winner)))) > (length(inds)*0.5)
%                             
%                                     varIDsreassign_taskunmatched.data(inds) = 0;
%                                     %diceOfWinner.data(inds) = 0;
%                                     
%                                 else
                                    
                                    taskunmatchednetnums(:,winner) = taskunmatchednetnums(:,winner) + 1;
                                    taskunmatchednetnumsSub(:,winner) = taskunmatchednetnumsSub(:,winner) + 1;
                            
%                                 end
                                if savestruct == 1

                                    StructAdd = size(Variant_Match_Struct,2)+1;

                                    Variant_Match_Struct(StructAdd).Subject = subject;
                                    Variant_Match_Struct(StructAdd).Matched = 0;
                                    Variant_Match_Struct(StructAdd).State = 'Task';
                                    Variant_Match_Struct(StructAdd).Network = network_names(winner);
                                    Variant_Match_Struct(StructAdd).Rest_Var_Num = 0;
                                    Variant_Match_Struct(StructAdd).Task_Var_Num = vars_task(var);

                                end
                            end
        
                        end
                        
                        if getmatchbystate == 1
                            
                            for v = 1:size(FinalMatchedVars,1)
                                
                                inds = find(cifti_task.data == FinalMatchedVars(v,2));        
                                varRmat = paircorr_mod(cifti_task_timeseries.data(inds,:)',cifti_task_timeseries.data'); %correlate each vertex within variant to all vertices in brain
                                varmean = mean(varRmat,1); %average seedmap of variant
                                %fprintf([num2str(var) '  '])
        
                                temp = sort(varmean,'descend');
                                varthresh = temp(round(0.05 * numel(temp)));
                                varmean_thresh = varmean > varthresh;

                                for net = 1:14
                            
                                    %%% dice calculation b/w variant & each network %%%            
                                    netTemplate = cifti_template_mat_full(net,:);
                                    var_dice_to_template(net,:) = (sum((netTemplate & varmean_thresh),2)*2) ./ sum((netTemplate | varmean_thresh),2);
            
                                end
            
                                %%% pick highest eta value and assign variant to that network %%%
                                [maxval, winner] = max(var_dice_to_template,[],1);
                                
                                if getmaxmatch == 1
                                    
                                    maxmatchvals = [maxmatchvals; maxval];
                                    
                                end
                            
                                if maxval > 0 && (excludemaxmatch == 1 && maxval > minmatchval)
                            
                                    winner = wb_colors(winner);
                                
                                else
                                
                                    winner = 4;
                                
                                end
                            
                                matchedvarsnetstask = [matchedvarsnetstask; network_names(winner)];
                                
                                %diceOfWinner.data(inds)=x;

                                %%% remove variant if over half its vertices overlap w consensus %%%
%                                 if length(find(ismember(inds,find(consensus.data==winner)))) > (length(inds)*0.5)
%                             
%                                     matchedvarsnetstask{end} = ' ';
%                                     
%                                 end
                            end
                        end
                        
                        if SizeExclude == 1
                            
                            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Variant_Network_Matching/' subject '_MatchedVariants_Reassigned_Task_' num2str(threshold) '_pct_SizeExclude'],varIDsreassign_taskmatched)
                            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Variant_Network_Matching/' subject '_UnmatchedVariants_Reassigned_Task_' num2str(threshold) '_pct_SizeExclude'],varIDsreassign_taskunmatched)
                            
                        else
                        
                            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Variant_Network_Matching/' subject '_MatchedVariants_Reassigned_Task_' num2str(threshold) '_pct'],varIDsreassign_taskmatched)
                            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Variant_Network_Matching/' subject '_UnmatchedVariants_Reassigned_Task_' num2str(threshold) '_pct'],varIDsreassign_taskunmatched)
                        
                        end
                            
                        clear varIDsreassign_taskmatched
                        clear varIDsreassign_taskunmatched
                        
                        if subjectdiceplots == 1
                            
                            outputdir2 = '/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Variant_Network_Matching/';

                            if SizeExclude == 1
            
                                filename = ['Variant_Plot_By_State_And_Match_' num2str(threshold) 'pct_SizeExclude_' subject '.jpg'];
            
                            else
        
                                filename = ['Variant_Plot_By_State_And_Match_' num2str(threshold) 'pct_' subject '.jpg'];
            
                            end
        
                            plot(1:16, restunmatchednetnumsSub, 1:16, restmatchednetnumsSub, 1:16, taskunmatchednetnumsSub, 1:16, taskmatchednetnumsSub)
                            set(gca,'xtick',[1:16],'xticklabel',network_names)
                            legend('Rest Unmatched Variants', 'Rest Matched Variants', 'Task Unmatched Variants', 'Task Matched Variants')
        
                            saveas(gcf,[outputdir2 filename])
                            
                            close gcf
                            
                        end
                        
                      	restunmatchednetnumsAllSubs(x,:) = restunmatchednetnumsSub;
                       	restmatchednetnumsAllSubs(x,:) = restmatchednetnumsSub;
                      	taskunmatchednetnumsAllSubs(x,:) = taskunmatchednetnumsSub;
                       	taskmatchednetnumsAllSubs(x,:) = taskmatchednetnumsSub;
                
%                     else
%             
%                         %%% Unclear how to correlate this timeseries (i.e.
%                         %%% variants matched across states have two
%                         %%% timeseries)
%                         
%                         unmatchedvarsnets = zeros(size(cifti_rest.data));
%                         matchedvarsnets = zeros(size(cifti_rest.data));
%                 
%                     end
                    
                else
                       
                    varIDsreassign = ft_read_cifti_mod('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Rest_Data/Full_Rest/Full_Data/MSC01/MSC01_REST_AllSessions_vs_120_allsubs_corr_cortex_corr.dtseries.nii');
                    varIDsreassign.data = zeros(size(cifti_rest.data));
                    
                    %%% loop through variants, get avg seedmap, match to network templates %%%
                    for var = 1:length(vars_rest)
                        
                        inds = find(cifti_rest.data == vars_rest(var));        
                        varRmat = paircorr_mod(cifti_rest_timeseries.data(inds,:)',cifti_rest_timeseries.data'); %correlate each vertex within variant to all vertices in brain
                        varmean = mean(varRmat,1); %average seedmap of variant
                        %fprintf([num2str(var) '  '])
        
                        temp = sort(varmean,'descend');
                        varthresh = temp(round(0.05 * numel(temp)));
                        varmean_thresh = varmean > varthresh;

                        for net = 1:14
                            
                            %%% dice calculation b/w variant & each network %%%            
                            netTemplate = cifti_template_mat_full(net,:);
                            var_dice_to_template(net,:) = (sum((netTemplate & varmean_thresh),2)*2) ./ sum((netTemplate | varmean_thresh),2);
            
                        end
            
                        %%% pick highest eta value and assign variant to that network %%%
                     	[maxval, winner] = max(var_dice_to_template,[],1);
                        
                      	if getmaxmatch == 1
                                    
                         	maxmatchvals = [maxmatchvals; maxval];
                                    
                      	end
                            
                      	if maxval > 0 && (excludemaxmatch == 1 && maxval > minmatchval)
                            
                        	winner = wb_colors(winner);
                                
                        else
                                
                          	winner = 4;
                                
                       	end
                            
                        varIDsreassign.data(inds) = winner;
        
                        %diceOfWinner.data(inds)=x;

                        %%% remove variant if over half its vertices overlap w consensus %%%
%                         if length(find(ismember(inds,find(consensus.data==winner)))) > (length(inds)*0.5)
%                             
%                             varIDsreassign.data(inds) = 0;
%                             %diceOfWinner.data(inds) = 0;
%                             
%                         else
                            
                            restnetnums(:,winner) = restnetnums(:,winner) + 1;
                            restnetnumsSub(:,winner) = restnetnumsSub(:,winner) + 1;
                            
%                         end
        
                    end
                    
                    if SizeExclude == 1
                        
                        ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Variant_Network_Matching/' subject '_AllVariants_Reassigned_Rest_' num2str(threshold) '_pct_SizeExclude'],varIDsreassign)
                        
                    else
                    
                        ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Variant_Network_Matching/' subject '_AllVariants_Reassigned_Rest_' num2str(threshold) '_pct'],varIDsreassign)
                    
                    end
                        
                    varIDsreassign.data = zeros(size(cifti_rest.data));
                    
                    %%% loop through variants, get avg seedmap, match to network templates %%%
                    for var = 1:length(vars_task)
                        
                        inds = find(cifti_task.data == vars_task(var));        
                        varRmat = paircorr_mod(cifti_task_timeseries.data(inds,:)',cifti_task_timeseries.data'); %correlate each vertex within variant to all vertices in brain
                        varmean = mean(varRmat,1); %average seedmap of variant
                        %fprintf([num2str(var) '  '])
        
                        temp = sort(varmean,'descend');
                        varthresh = temp(round(0.05 * numel(temp)));
                        varmean_thresh = varmean > varthresh;

                        for net = 1:14
                            
                            %%% dice calculation b/w variant & each network %%%            
                            netTemplate = cifti_template_mat_full(net,:);
                            var_dice_to_template(net,:) = (sum((netTemplate & varmean_thresh),2)*2) ./ sum((netTemplate | varmean_thresh),2);
            
                        end
            
                        %%% pick highest eta value and assign variant to that network %%%
                      	[maxval, winner] = max(var_dice_to_template,[],1);
                        
                        if getmaxmatch == 1
                                    
                          	maxmatchvals = [maxmatchvals; maxval];
                                    
                      	end
                            
                      	if maxval > 0 && (excludemaxmatch == 1 && maxval > minmatchval)
                            
                        	winner = wb_colors(winner);
                                
                        else
                                
                         	winner = 4;
                                
                       	end
                            
                        varIDsreassign.data(inds) = winner;
        
                        %diceOfWinner.data(inds)=x;

                        %%% remove variant if over half its vertices overlap w consensus %%%
%                         if length(find(ismember(inds,find(consensus.data==winner)))) > (length(inds)*0.5)
%                             
%                             varIDsreassign.data(inds) = 0;
%                             %diceOfWinner.data(inds) = 0;
%                             
%                         else
                            
                            tasknetnums(:,winner) = tasknetnums(:,winner) + 1;
                            tasknetnumsSub(:,winner) = tasknetnumsSub(:,winner) + 1;
                            
%                         end

                    end
                    
                    if SizeExclude == 1
                        
                        ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Variant_Network_Matching/' subject '_AllVariants_Reassigned_Task_' num2str(threshold) '_pct_SizeExclude'],varIDsreassign)
                        
                    else
                    
                        ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Variant_Network_Matching/' subject '_AllVariants_Reassigned_Task_' num2str(threshold) '_pct'],varIDsreassign)
                    
                    end
                    
                    if subjectdiceplots == 1
                        
                        outputdir2 = '/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Variant_Network_Matching/';
                    
                        if SizeExclude == 1
            
                            filename = ['Variant_Plot_By_State_' num2str(threshold) 'pct_SizeExclude_' subject '.jpg'];

                        else
        
                            filename = ['Variant_Plot_By_State_' num2str(threshold) 'pct_' subject '.jpg'];
            
                        end
    
                        plot(1:16, restnetnumsSub, 1:16, tasknetnumsSub)
                        set(gca,'xtick',[1:16],'xticklabel',network_names)
                        legend('Rest Variants','Task Variants')
        
                        saveas(gcf,[outputdir2 filename])
                        
                        close gcf
                        
                        restnetnumsAllSubs(x,:) = restnetnumsSub;
                        tasknetnumsAllSubs(x,:) = tasknetnumsSub;
                    
                    end
                end
                
                
            end
                      
                if separatetaskandrest == 1 && separatematchedvars == 1
                
                    unmatchedvarsresttemp = zeros(size(cifti_rest.data));
                    matchedvarsresttemp = zeros(size(cifti_rest.data));
                    unmatchedvarstasktemp = zeros(size(cifti_rest.data));
                    matchedvarstasktemp = zeros(size(cifti_rest.data));
                    
                elseif separatematchedvars == 1
            
                    unmatchedvarstemp = zeros(size(cifti_rest.data));
                    matchedvarstemp = zeros(size(cifti_rest.data));
                    
                else
                    
                   	alltaskvarstemp = zeros(size(cifti_rest.data));
                    allrestvarstemp = zeros(size(cifti_rest.data));
                
                end
            
                for p = 1:length(vars_rest)
                
                    if ismember(vars_rest(p),FinalMatchedVars(:,1))
                    
                        overlapadd = find(vars_rest(p) == cifti_rest.data);
                    
                        if separatetaskandrest == 1 && separatematchedvars == 1
                        
                            matchedvarsresttemp(overlapadd,:) = 1;
                        
                        elseif separatematchedvars == 1
                    
                            matchedvarstemp(overlapadd,:) = 1;
                            
                        else
                            
                            allrestvarstemp(overlapadd,:) = 1;
                        
                        end
                    
                    else
                    
                        nonoverlapadd = find(vars_rest(p) == cifti_rest.data);
                    
                        if separatetaskandrest == 1 && separatematchedvars == 1
                        
                            unmatchedvarsresttemp(nonoverlapadd,:) = 1;
                        
                        elseif separatematchedvars == 1
                    
                            unmatchedvarstemp(nonoverlapadd,:) = 1;
                            
                        else 
                            
                            allrestvarstemp(nonoverlapadd,:) = 1;
                        
                        end
                    
                    end
                end
                    
                for r = 1:length(vars_task)
                
                    if ismember(vars_task(r),FinalMatchedVars(:,2))
                    
                        overlapadd = find(vars_task(r) == cifti_task.data);
                    
                        if separatetaskandrest == 1 && separatematchedvars == 1
                        
                            matchedvarstasktemp(overlapadd,:) = 1;
                        
                        elseif separatematchedvars == 1
                    
                            matchedvarstemp(overlapadd,:) = 1;
                            
                        else
                            
                            alltaskvarstemp(overlapadd,:) = 1;
                            
                        end
                    
                    else
                    
                        nonoverlapadd = find(vars_task(r) == cifti_task.data);
                    
                        if separatetaskandrest == 1 && separatematchedvars == 1
                    
                            unmatchedvarstasktemp(nonoverlapadd,:) = 1;
                        
                        elseif separatematchedvars == 1
                        
                            unmatchedvarstemp(nonoverlapadd,:) = 1;
                            
                        else
                            
                            alltaskvarstemp(nonoverlapadd,:) = 1;
                        
                        end
                    
                    end
                end
            
                if separatetaskandrest == 1 && separatematchedvars == 1
                
                    unmatchedvarsrest = unmatchedvarsrest + unmatchedvarsresttemp;
                    matchedvarsrest = matchedvarsrest + matchedvarsresttemp;
                    unmatchedvarstask = unmatchedvarstask + unmatchedvarstasktemp;
                    matchedvarstask = matchedvarstask + matchedvarstasktemp;
                
                elseif separatematchedvars == 1
            
                    unmatchedvars = unmatchedvars + unmatchedvarstemp;
                    matchedvars = matchedvars + matchedvarstemp;
                    
                else
                    
                    alltaskvars = alltaskvars + alltaskvarstemp;
                    allrestvars = allrestvars + allrestvarstemp;
                
                end   
    
    % Fill in matrix with new values
    
        if matchbydist == 1

            if size(COMCoordsresttemp,1) < size(COMCoordstasktemp,1)

                if excludeunmatched == 1

                    COMDistanceMatrixTemp = zeros(size(FinalMatchedVars,1),size(COMCoordstasktemp,1));     %% Match closest variants between task and rest
                    %COMTaskVarsMatched = zeros(size(COMCoordstasktemp,1),2);

                else

                    COMDistanceMatrixTemp = zeros(size(FinalMatchedVars,1),size(COMCoordstasktemp,1));     %% Match closest variants between task and rest
                    %COMTaskVarsMatched = zeros(size(COMCoordstasktemp,1),2);

                end

            else

                if excludeunmatched == 1

                    COMDistanceMatrixTemp = zeros(size(FinalMatchedVars,1),size(COMCoordstasktemp,1));     %% Match closest variants between task and rest
                    %COMTaskVarsMatched = zeros(size(COMCoordsresttemp,1),2);

                else

                    COMDistanceMatrixTemp = zeros(size(FinalMatchedVars,1),size(COMCoordsresttemp,1));     %% Match closest variants between task and rest
                    %COMTaskVarsMatched = zeros(size(COMCoordsresttemp,1),2);

                end

            end

            %COMDistanceMatrixTemp = COMDistanceMatrixTemp';

            if size(COMDistanceMatrixTemp,2) ~= size(COMDistanceMatrix,2)

                COMDistanceMatrixTemp = COMDistanceMatrixTemp';

            end

            for e = 1:size(COMDistanceMatrixTemp,2)

                    addvars = find(COMTaskVarsMatched(:,2) == e);

                if length(addvars) > 1

                    for c = 1:length(addvars)

                        if size(COMCoordsresttemp,1) < size(COMCoordstasktemp,1)

                            COMDistanceMatrixTemp(addvars(c),:) = COMDistanceMatrix(COMTaskVarsMatched(addvars(c),2),:);

                        else

                            COMDistanceMatrixTemp(addvars(c),:) = COMDistanceMatrix(COMTaskVarsMatched(addvars(c),1),:);

                        end

                    end

                else

                    if size(COMCoordsresttemp,1) < size(COMCoordstasktemp,1)

                        COMDistanceMatrixTemp(addvars,:) = COMDistanceMatrix(COMTaskVarsMatched(addvars,2),:);

                    else

                        COMDistanceMatrixTemp(addvars,:) = COMDistanceMatrix(COMTaskVarsMatched(addvars,1),:);

                    end

                end



            end
        end
    
    %clf
    %heatmap(COMDistanceMatrixTemp, 'Colormap', jet);
    
    %h = gca;
    
    if plotmats == 1
    
        outputdir = '/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/COMHistogramPlots';
    
        if SizeExclude == 1 && excludeunmatched == 1
          
            filename = ['/' subject '_VariantDistances_ExcludeUnmatched_' num2str(threshold) '_SNRSizeExclude_Matrix.jpg'];
        
        elseif SizeExclude == 1
        
            filename = ['/' subject '_VariantDistances_' num2str(threshold) '_SNRSizeExclude_Matrix.jpg'];
        
        elseif excludeunmatched == 1
        
            filename = ['/' subject '_VariantDistances_ExcludeUnmatched_' num2str(threshold) '_SNRExclude_Matrix.jpg'];
        
        else
    
            filename = ['/' subject '_VariantDistances_' num2str(threshold) '_SNRExclude_Matrix.jpg'];
        
        end
    
        if size(COMCoordsresttemp,1) < size(COMCoordstasktemp,1)
        
            imagesc(COMDistanceMatrixTemp')
            colorbar
            h = colorbar;
            ylabel(h, 'Distance between COMs (mm)', 'FontSize',20)
        
            xlabel('Task Variants','FontSize',24);
            ylabel('Rest Variants','FontSize',24);

        
            set(gca,'YTick',mean(1:length(COMTaskVarsMatched)), 'YTickLabel',{'Matched Variants'}, 'FontSize',20)
            ytickangle(90)
        
        
        else
    
            imagesc(COMDistanceMatrixTemp)
            colorbar
            h = colorbar;
            ylabel(h, 'Distance between COMs (mm)', 'FontSize',14)
        
            xlabel('Task Variants');
            ylabel('Rest Variants');

        
            set(gca,'XTick',mean(1:length(COMTaskVarsMatched)), 'XTickLabel',{'Matched Variants'}, 'FontSize',20)
            %xtickangle(90)
            
        end
    
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [.2, 0.1, .8, 0.9]);
    
        saveas(gcf,[outputdir filename])
    
        close gcf
    
    end
        %end
        
        disp(['Finished data for subject ' subject])
end


if savesurfaces == 1
    
	template = ft_read_cifti_mod('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Rest_Data/Full_Rest/Full_Data/MSC01/MSC01_REST_AllSessions_vs_120_allsubs_corr_cortex_corr.dtseries.nii');
	template.data = [];
            
	template2 = template;
    
    if separatetaskandrest == 1
        
        template3 = template;
        template4 = template;
        
        template.data = matchedvarsrest;
        template2.data = unmatchedvarsrest;
        template3.data = matchedvarstask;
        template4.data = unmatchedvarstask;
        
        if matchbydist == 0 && excludemaxmatch == 1 && SizeExclude == 1
            
            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Task_vs_Rest_Data/AllSubjects_RestOnly_MatchedVariants_BothSNRAndSizeExcluded_OverlapOnly_ExcludeMatch_Map_' num2str(threshold) '_pct'],template)
            clear template

            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Task_vs_Rest_Data/AllSubjects_RestOnly_UnmatchedVariants_BothSNRAndSizeExcluded_OverlapOnly_ExcludeMatch_Map_' num2str(threshold) '_pct'],template2)
            clear template2
            
            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Task_vs_Rest_Data/AllSubjects_TaskOnly_MatchedVariants_BothSNRAndSizeExcluded_OverlapOnly_ExcludeMatch_Map_' num2str(threshold) '_pct'],template3)
            clear template3

            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Task_vs_Rest_Data/AllSubjects_TaskOnly_UnmatchedVariants_BothSNRAndSizeExcluded_OverlapOnly_ExcludeMatch_Map_' num2str(threshold) '_pct'],template4)
            clear template4
            
        elseif matchbydist == 0 && SizeExclude == 1
            
            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Task_vs_Rest_Data/AllSubjects_RestOnly_MatchedVariants_BothSNRAndSizeExcluded_OverlapOnly_Map_' num2str(threshold) '_pct'],template)
            clear template

            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Task_vs_Rest_Data/AllSubjects_RestOnly_UnmatchedVariants_BothSNRAndSizeExcluded_OverlapOnly_Map_' num2str(threshold) '_pct'],template2)
            clear template2
            
            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Task_vs_Rest_Data/AllSubjects_TaskOnly_MatchedVariants_BothSNRAndSizeExcluded_OverlapOnly_Map_' num2str(threshold) '_pct'],template3)
            clear template3

            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Task_vs_Rest_Data/AllSubjects_TaskOnly_UnmatchedVariants_BothSNRAndSizeExcluded_OverlapOnly_Map_' num2str(threshold) '_pct'],template4)
            clear template4
        
        elseif SizeExclude == 1
    
            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Task_vs_Rest_Data/AllSubjects_RestOnly_MatchedVariants_BothSNRAndSizeExcluded_Map_' num2str(threshold) '_pct'],template)
            clear template

            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Task_vs_Rest_Data/AllSubjects_RestOnly_UnmatchedVariants_BothSNRAndSizeExcluded_Map_' num2str(threshold) '_pct'],template2)
            clear template2
            
            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Task_vs_Rest_Data/AllSubjects_TaskOnly_MatchedVariants_BothSNRAndSizeExcluded_Map_' num2str(threshold) '_pct'],template3)
            clear template3

            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Task_vs_Rest_Data/AllSubjects_TaskOnly_UnmatchedVariants_BothSNRAndSizeExcluded_Map_' num2str(threshold) '_pct'],template4)
            clear template4
    
        else
        
            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Task_vs_Rest_Data/AllSubjects_RestOnly_MatchedVariants_SNRExcluded_Map_' num2str(threshold) '_pct'],template)
            clear template

            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Task_vs_Rest_Data/AllSubjects_RestOnly_UnmatchedVariants_SNRExcluded_Map_' num2str(threshold) '_pct'],template2)
            clear template2
            
            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Task_vs_Rest_Data/AllSubjects_TaskOnly_MatchedVariants_SNRExcluded_Map_' num2str(threshold) '_pct'],template3)
            clear template3

            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Task_vs_Rest_Data/AllSubjects_TaskOnly_UnmatchedVariants_SNRExcluded_Map_' num2str(threshold) '_pct'],template4)
            clear template4
        
        end
        
    elseif separatematchedvars == 1
    
        template.data = matchedvars;
        template2.data = unmatchedvars;
    
        if SizeExclude == 1
    
            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Task_vs_Rest_Data/AllSubjects_TaskvsRest_MatchedVariants_BothSNRAndSizeExcluded_Map_' num2str(threshold) '_pct'],template)
            clear template

            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Task_vs_Rest_Data/AllSubjects_TaskvsRest_UnmatchedVariants_BothSNRAndSizeExcluded_Map_' num2str(threshold) '_pct'],template2)
            clear template2
    
        else
        
            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Task_vs_Rest_Data/AllSubjects_TaskvsRest_MatchedVariants_SNRExcluded_Map_' num2str(threshold) '_pct'],template)
            clear template

            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Task_vs_Rest_Data/AllSubjects_TaskvsRest_UnmatchedVariants_SNRExcluded_Map_' num2str(threshold) '_pct'],template2)
            clear template2
        
        end
        
    else
        
        template.data = matchedvars;
        template2.data = unmatchedvars;
    
        if SizeExclude == 1
    
            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Task_vs_Rest_Data/AllSubjects_TaskvsRest_MatchedVariants_BothSNRAndSizeExcluded_Map_' num2str(threshold) '_pct'],template)
            clear template

            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Task_vs_Rest_Data/AllSubjects_TaskvsRest_UnmatchedVariants_BothSNRAndSizeExcluded_Map_' num2str(threshold) '_pct'],template2)
            clear template2
    
        else
        
            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Task_vs_Rest_Data/AllSubjects_TaskvsRest_MatchedVariants_SNRExcluded_Map_' num2str(threshold) '_pct'],template)
            clear template

            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Task_vs_Rest_Data/AllSubjects_TaskvsRest_UnmatchedVariants_SNRExcluded_Map_' num2str(threshold) '_pct'],template2)
            clear template2
        
        end
        
        
    end
end

if plothists == 1

    for p = 1:nfiles
    
        subject = subjects2{p};

        if SizeExclude == 1
                
            filename = ['/' subject '_VariantDistances_SNRSizeExclude_Histogram.jpg'];

        else 
            
            filename = ['/' subject '_VariantDistances_SNRExclude_Histogram.jpg'];
                
        end
    
        nhist(COMDistanceMatrixAll{p,:}, 'binfactor', 2, 'samebins');
        title(['Variant Minimum Distances ' subject], 'fontsize',18)
        xlabel('Distance (mm)')
        ax = gca;
        ax.FontSize = 14;
        ylim([0 50]);
    
        saveas(gcf,[outputdir filename])
    
        close gcf
        
    end



    if SizeExclude == 1

        filename = ['/AllSubjects_VariantDistances_SNRSizeExclude_Histogram.jpg'];

    else
    
        filename = ['/AllSubjects_VariantDistances_SNRExclude_Histogram.jpg'];
    
    end

    nhist(COMDistanceMatrixAll, 'binfactor', 2, 'samebins');
    title('Variant Distances All Subjects', 'fontsize',18)
    xlabel('Distance (mm)')
    ax = gca;
    ax.FontSize = 14;
    
    saveas(gcf,[outputdir filename])
    
    close gcf
    
end

if DiceNetworks == 1
    
   % Plot line graphs
    
    outputdir2 = '/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Variant_Network_Matching/';
    
    if separatematchedvars == 1
        
        if SizeExclude == 1
            
            filename = ['Variant_Plot_By_State_And_Match_' num2str(threshold) 'pct_SizeExclude_AllSubjects.jpg'];
            
        else
        
            filename = ['Variant_Plot_By_State_And_Match_' num2str(threshold) 'pct_AllSubjects.jpg'];
            
        end
        
        plot(1:16, restunmatchednetnums, 1:16, restmatchednetnums, 1:16, taskunmatchednetnums, 1:16, taskmatchednetnums)
        
      	set(gca,'xtick',[1:16],'xticklabel',network_names)
        legend('Rest Unmatched Variants', 'Rest Matched Variants', 'Task Unmatched Variants', 'Task Matched Variants')
        
      	saveas(gcf,[outputdir2 filename])
        
        close gcf
        
        if SizeExclude == 1
            
            filename = ['Variant_Plot_By_State_And_Match_MeanErrorBars' num2str(threshold) 'pct_SizeExclude_AllSubjects.jpg'];
            
        else
        
            filename = ['Variant_Plot_By_State_And_Match_MeanErrorBars' num2str(threshold) 'pct_AllSubjects.jpg'];
            
        end
        
    	SEplotrestunmatched = std(restunmatchednetnumsAllSubs,0,1)./sqrt(size(restunmatchednetnumsAllSubs,1));
     	SEplotrestmatched = std(restmatchednetnumsAllSubs,0,1)./sqrt(size(restmatchednetnumsAllSubs,1));
    	SEplottaskunmatched = std(taskunmatchednetnumsAllSubs,0,1)./sqrt(size(taskunmatchednetnumsAllSubs,1));
     	SEplottaskmatched = std(taskmatchednetnumsAllSubs,0,1)./sqrt(size(taskmatchednetnumsAllSubs,1));
        
        errorbar(1:16, mean(restunmatchednetnumsAllSubs,1), SEplotrestunmatched)
        hold on
        errorbar(1:16, mean(restmatchednetnumsAllSubs,1), SEplotrestmatched)
        hold on
        errorbar(1:16, mean(taskunmatchednetnumsAllSubs,1), SEplottaskunmatched)
        hold on
        errorbar(1:16, mean(taskmatchednetnumsAllSubs,1), SEplottaskmatched)
        hold on
        set(gca,'xtick',[1:16],'xticklabel',network_names)
        ylim([0 10])

        legend({'Rest Unmatched Variants', 'Rest Matched Variants', 'Task Unmatched Variants', 'Task Matched Variants'})
        
        saveas(gcf,[outputdir2 filename])
        
        close gcf
        
        % Plot correlations between lines
        
        if SizeExclude == 1
            
            filename = ['Variant_Correlations_Across_States_' num2str(threshold) 'pct_SizeExclude_AllSubjects.jpg'];
            
        else
        
            filename = ['Variant_Correlations_Across_States_' num2str(threshold) 'pct_AllSubjects.jpg'];
            
        end
        
        coefmatmatched = [];
        coefmatunmatched = [];
        coefrestmatmatched = [];
        coeftaskmatmatched = [];
        
        for z = 1:size(restunmatchednetnumsAllSubs,1)
            
            rmatch = corrcoef(restmatchednetnumsAllSubs(z,:),taskmatchednetnumsAllSubs(z,:));
            runmatch = corrcoef(restunmatchednetnumsAllSubs(z,:),taskunmatchednetnumsAllSubs(z,:));
            restmatch = corrcoef(restunmatchednetnumsAllSubs(z,:),restmatchednetnumsAllSubs(z,:));
            taskmatch = corrcoef(taskunmatchednetnumsAllSubs(z,:),taskmatchednetnumsAllSubs(z,:));
            
            coefrestmatmatched = [coefrestmatmatched;restmatch(1,2)];
            coeftaskmatmatched = [coeftaskmatmatched;taskmatch(1,2)];
         	coefmatmatched = [coefmatmatched;rmatch(1,2)];
            coefmatunmatched = [coefmatunmatched;runmatch(1,2)];
            
        end
        
    	SEplotcoefrestmatched = std(coefrestmatmatched)./sqrt(size(coefrestmatmatched,1));
     	SEplotcoeftaskmatched = std(coeftaskmatmatched)./sqrt(size(coeftaskmatmatched,1));
    	SEplotcoefunmatched = std(coefmatunmatched)./sqrt(size(coefmatunmatched,1));
     	SEplotcoefmatched = std(coefmatmatched)./sqrt(size(coefmatmatched,1));
        
        bar(categorical({'Matched Variants (Task vs. Rest)', 'Unmatched Variants (Task vs. Rest)', 'Rest Variants (Matched vs. Unmatched)', 'Task Variants (Matched vs. Unmatched)'}), [mean(coefmatmatched) mean(coefmatunmatched) mean(coefrestmatmatched) mean(coeftaskmatmatched)], 'FaceColor', 'w', 'EdgeColor', 'k');
        hold on
        errorbar(categorical({'Matched Variants (Task vs. Rest)', 'Unmatched Variants (Task vs. Rest)', 'Rest Variants (Matched vs. Unmatched)', 'Task Variants (Matched vs. Unmatched)'}),[mean(coefmatmatched) mean(coefmatunmatched) mean(coefrestmatmatched) mean(coeftaskmatmatched)],[SEplotcoefmatched SEplotcoefunmatched SEplotcoefrestmatched SEplotcoeftaskmatched],'.');
        
     	%saveas(gcf,[outputdir2 filename])
        
        print(gcf,[outputdir2 filename],'-dpng','-r300');
        
        close gcf
        
        % Create pie charts
        
        % Create legend first
        
        networksexist = sum([restmatchednetnums;restunmatchednetnums;taskmatchednetnums;taskunmatchednetnums],1);
        netsplot = find(networksexist > 0);
        network_legend = network_names(netsplot);
        color_legend = rgb_colors(netsplot,:);

        filename = ['Variant_PieChart_Legend_' num2str(threshold) 'pct_SizeExclude_AllSubjects_OverlapOnly_ExcludeMatch.jpg'];
        
        h = pie(networksexist);
        delete(findobj(h,'Type','text'))
        
        for k = 1:length(h)/2
            
            set(h(k*2-1), 'FaceColor', color_legend(k,:));
            
        end
        
        legend(network_legend,'Orientation','vertical','Location','eastoutside', 'FontSize',20)
        
        %saveas(gcf,[outputdir2 filename])
        
        print(gcf,[outputdir2 filename],'-dpng','-r300');
        
        close gcf
        
        
        if SizeExclude == 1
            
            if matchbydist == 0 && excludemaxmatch == 1
                
                filename = ['Variant_PieChart_Unmatched_Rest_' num2str(threshold) 'pct_SizeExclude_AllSubjects_OverlapOnly_ExcludeMatch.jpg'];
            
            elseif matchbydist == 0
                
                filename = ['Variant_PieChart_Unmatched_Rest_' num2str(threshold) 'pct_SizeExclude_AllSubjects_OverlapOnly.jpg'];
                
            else
            
                filename = ['Variant_PieChart_Unmatched_Rest_' num2str(threshold) 'pct_SizeExclude_AllSubjects.jpg'];
            
            end
                
        else
        
            filename = ['Variant_PieChart_Unmatched_Rest_' num2str(threshold) 'pct_AllSubjects.jpg'];
            
        end
        
        %restunmatchednetnumsplot = [restunmatchednetnums(1) restunmatchednetnums(16) restunmatchednetnums(2) restunmatchednetnums(3) restunmatchednetnums(4) restunmatchednetnums(15) restunmatchednetnums(5) restunmatchednetnums(6) restunmatchednetnums(7) restunmatchednetnums(8) restunmatchednetnums(9) restunmatchednetnums(10) restunmatchednetnums(11) restunmatchednetnums(12) restunmatchednetnums(13) restunmatchednetnums(14)];
        %network_namesplot = {'DMN'	'PON'  'Vis'	'FP'	''	'PMN'   'DAN'	''	'VAN'	'Sal'	'CO'	'SMd'	'SMl'	'Aud'	'Tpole'	'MTL'};
        
        % create plot to get colors for final chart
        
        h = pie(restunmatchednetnums, network_names);
        %title('', 'fontsize',14)
        set(findobj(h,'type','text'),'fontsize',14)
        
        pieColorMap = [];
        
        for k = 1:length(h)/2
            
            for l = 1:length(network_names)
                
                if strcmp(h(k*2).String,char(network_names(l)))
            
                    % Create a color for this sector of the pie
                    pieColorMap = [pieColorMap; rgb_colors(l,:)];  % Color for this segment.
                    
                end
            
                % Apply the colors we just generated to the pie chart.
%                 set(hPieComponentHandles(k*2-1), 'FaceColor', pieColorMap);
%                 set(hPieComponentHandles(k*2), 'String', num2str(X(k)), 'FontSize', fontSize );
             
            end
        end
        
        close gcf
        
        % Create new pie chart and apply colors
        
        h = pie(restunmatchednetnums);
        delete(findobj(h,'Type','text'))
        %set(findobj(h,'type','text'),'fontsize',11)
        
        for k = 1:length(h)/2
            
            set(h(k*2-1), 'FaceColor', pieColorMap(k,:));
            
        end
        
        
        %saveas(gcf,[outputdir2 filename])
        
        print(gcf,[outputdir2 filename],'-dpng','-r300');
        
        close gcf
        
        if SizeExclude == 1
            
            if matchbydist == 0 && excludemaxmatch == 1
                
                filename = ['Variant_PieChart_Matched_Rest_' num2str(threshold) 'pct_SizeExclude_AllSubjects_OverlapOnly_ExcludeMatch.jpg'];
            
            elseif matchbydist == 0
                
                filename = ['Variant_PieChart_Matched_Rest_' num2str(threshold) 'pct_SizeExclude_AllSubjects_OverlapOnly.jpg'];
                
            else
            
                filename = ['Variant_PieChart_Matched_Rest_' num2str(threshold) 'pct_SizeExclude_AllSubjects.jpg'];
                
            end
            
        else
        
            filename = ['Variant_PieChart_Matched_Rest_' num2str(threshold) 'pct_AllSubjects.jpg'];
            
        end
        
        %restmatchednetnumsplot = [restmatchednetnums(1) restmatchednetnums(16) restmatchednetnums(2) restmatchednetnums(3) restmatchednetnums(4) restmatchednetnums(15) restmatchednetnums(5) restmatchednetnums(6) restmatchednetnums(7) restmatchednetnums(8) restmatchednetnums(9) restmatchednetnums(10) restmatchednetnums(11) restmatchednetnums(12) restmatchednetnums(13) restmatchednetnums(14)];
        %network_namesplot = {'DMN'	'PON'  'Vis'	'FP'	''	'PMN'   'DAN'	''	'VAN'	'Sal'	'CO'	'SMd'	'SMl'	'Aud'	'Tpole'	'MTL'};
        
        h = pie(restmatchednetnums,network_names);
        %title('', 'fontsize',14)
        set(findobj(h,'type','text'),'fontsize',11)
        
        pieColorMap = [];
        
        for k = 1:length(h)/2
            
            for l = 1:length(network_names)
                
                if strcmp(h(k*2).String,char(network_names(l)))
            
                    % Create a color for this sector of the pie
                    pieColorMap = [pieColorMap; rgb_colors(l,:)];  % Color for this segment.
                    
                end
            
                % Apply the colors we just generated to the pie chart.
%                 set(hPieComponentHandles(k*2-1), 'FaceColor', pieColorMap);
%                 set(hPieComponentHandles(k*2), 'String', num2str(X(k)), 'FontSize', fontSize );
             
            end
        end
        
        close gcf
        
        % Create new pie chart and apply colors
        
        h = pie(restmatchednetnums);
        delete(findobj(h,'Type','text'))
        %set(findobj(h,'type','text'),'fontsize',11)
        
        for k = 1:length(h)/2
            
            set(h(k*2-1), 'FaceColor', pieColorMap(k,:));
            
        end
        
        %saveas(gcf,[outputdir2 filename])
        
        print(gcf,[outputdir2 filename],'-dpng','-r300');
        
        close gcf
        
        if SizeExclude == 1
            
            if matchbydist == 0 && excludemaxmatch == 1
                
                filename = ['Variant_PieChart_Unmatched_Task_' num2str(threshold) 'pct_SizeExclude_AllSubjects_OverlapOnly_ExcludeMatch.jpg'];
            
            elseif matchbydist == 0
                
                filename = ['Variant_PieChart_Unmatched_Task_' num2str(threshold) 'pct_SizeExclude_AllSubjects_OverlapOnly.jpg'];
                
            else
            
                filename = ['Variant_PieChart_Unmatched_Task_' num2str(threshold) 'pct_SizeExclude_AllSubjects.jpg'];
            
            end
            
        else
        
            filename = ['Variant_PieChart_Unmatched_Task_' num2str(threshold) 'pct_AllSubjects.jpg'];
            
        end
        
        %taskunmatchednetnumsplot = [taskunmatchednetnums(1) taskunmatchednetnums(16) taskunmatchednetnums(2) taskunmatchednetnums(3) taskunmatchednetnums(4) taskunmatchednetnums(15) taskunmatchednetnums(5) taskunmatchednetnums(6) taskunmatchednetnums(7) taskunmatchednetnums(8) taskunmatchednetnums(9) taskunmatchednetnums(10) taskunmatchednetnums(11) taskunmatchednetnums(12) taskunmatchednetnums(13) taskunmatchednetnums(14)];
        %network_namesplot = {'DMN'	'PON'  'Vis'	'FP'	''	'PMN'   'DAN'	''	'VAN'	'Sal'	'CO'	'SMd'	'SMl'	'Aud'	'Tpole'	'MTL'};
        
        h = pie(taskunmatchednetnums,network_names);
        %title('', 'fontsize',14)
        set(findobj(h,'type','text'),'fontsize',11)
        
        pieColorMap = [];
        
        for k = 1:length(h)/2
            
            for l = 1:length(network_names)
                
                if strcmp(h(k*2).String,char(network_names(l)))
            
                    % Create a color for this sector of the pie
                    pieColorMap = [pieColorMap; rgb_colors(l,:)];  % Color for this segment.
                    
                end
            
                % Apply the colors we just generated to the pie chart.
%                 set(hPieComponentHandles(k*2-1), 'FaceColor', pieColorMap);
%                 set(hPieComponentHandles(k*2), 'String', num2str(X(k)), 'FontSize', fontSize );
             
            end
        end
        
        close gcf
        
        % Create new pie chart and apply colors
        
        h = pie(taskunmatchednetnums);
        delete(findobj(h,'Type','text'))
        %set(findobj(h,'type','text'),'fontsize',11)
        
        for k = 1:length(h)/2
            
            set(h(k*2-1), 'FaceColor', pieColorMap(k,:));
            
        end
        
        %saveas(gcf,[outputdir2 filename])
        
        print(gcf,[outputdir2 filename],'-dpng','-r300');
        
        close gcf
        
        if SizeExclude == 1
            
            if matchbydist == 0 && excludemaxmatch == 1
                
                filename = ['Variant_PieChart_Matched_Task_' num2str(threshold) 'pct_SizeExclude_AllSubjects_OverlapOnly_ExcludeMatch.jpg'];
            
            elseif matchbydist == 0
                
                filename = ['Variant_PieChart_Matched_Task_' num2str(threshold) 'pct_SizeExclude_AllSubjects_OverlapOnly.jpg'];
                
            else
            
                filename = ['Variant_PieChart_Matched_Task_' num2str(threshold) 'pct_SizeExclude_AllSubjects.jpg'];
                
            end
            
        else
        
            filename = ['Variant_PieChart_Matched_Task_' num2str(threshold) 'pct_AllSubjects.jpg'];
            
        end
        
        %taskmatchednetnumsplot = [taskmatchednetnums(1) taskmatchednetnums(16) taskmatchednetnums(2) taskmatchednetnums(3) taskmatchednetnums(4) taskmatchednetnums(15) taskmatchednetnums(5) taskmatchednetnums(6) taskmatchednetnums(7) taskmatchednetnums(8) taskmatchednetnums(9) taskmatchednetnums(10) taskmatchednetnums(11) taskmatchednetnums(12) taskmatchednetnums(13) taskmatchednetnums(14)];
        %network_namesplot = {'DMN'	'PON'  'Vis'	'FP'	''	'PMN'   'DAN'	''	'VAN'	'Sal'	'CO'	'SMd'	'SMl'	'Aud'	'Tpole'	'MTL'};
        
        h = pie(taskmatchednetnums,network_names);
        %title('', 'fontsize',14)
        set(findobj(h,'type','text'),'fontsize',11)
        
        pieColorMap = [];
        
        for k = 1:length(h)/2
            
            for l = 1:length(network_names)
                
                if strcmp(h(k*2).String,char(network_names(l)))
            
                    % Create a color for this sector of the pie
                    pieColorMap = [pieColorMap; rgb_colors(l,:)];  % Color for this segment.
                    
                end
            
                % Apply the colors we just generated to the pie chart.
%                 set(hPieComponentHandles(k*2-1), 'FaceColor', pieColorMap);
%                 set(hPieComponentHandles(k*2), 'String', num2str(X(k)), 'FontSize', fontSize );
             
            end
        end
        
        close gcf
        
        % Create new pie chart and apply colors
        
        h = pie(taskmatchednetnums);
        delete(findobj(h,'Type','text'))
        %set(findobj(h,'type','text'),'fontsize',11)
        
        for k = 1:length(h)/2
            
            set(h(k*2-1), 'FaceColor', pieColorMap(k,:));
            
        end
        
        %saveas(gcf,[outputdir2 filename])
        
        print(gcf,[outputdir2 filename],'-dpng','-r300');
        
        close gcf
        
    else
        
        if SizeExclude == 1
            
            if matchbydist == 0 && excludemaxmatch == 1
                
                filename = ['Variant_Plot_By_State_' num2str(threshold) 'pct_SizeExclude_AllSubjects_OverlapOnly_ExcludeMatch.jpg'];
            
            elseif matchbydist == 0
                
                filename = ['Variant_Plot_By_State_' num2str(threshold) 'pct_SizeExclude_AllSubjects_OverlapOnly.jpg'];
                
            else
            
                filename = ['Variant_Plot_By_State_' num2str(threshold) 'pct_SizeExclude_AllSubjects.jpg'];
                
            end

        else
        
            filename = ['Variant_Plot_By_State_' num2str(threshold) 'pct_AllSubjects.jpg'];
            
        end
    
        plot(1:16, restnetnums, 1:16, tasknetnums)
        set(gca,'xtick',[1:16],'xticklabel',network_names)
        legend('Rest Variants','Task Variants')
        
        %saveas(gcf,[outputdir2 filename])
        
        print(gcf,[outputdir2 filename],'-dpng','-r300');
        
        close gcf
        
        % Create pie charts
        
        if SizeExclude == 1
            
            if matchbydist == 0 && excludemaxmatch == 1
                
                filename = ['Variant_PieChart_All_Rest_' num2str(threshold) 'pct_SizeExclude_AllSubjects_OverlapOnly_ExcludeMatch.jpg'];
            
            elseif matchbydist == 0
                
                filename = ['Variant_PieChart_All_Rest_' num2str(threshold) 'pct_SizeExclude_AllSubjects_OverlapOnly.jpg'];
                
            else
            
                filename = ['Variant_PieChart_All_Rest_' num2str(threshold) 'pct_SizeExclude_AllSubjects.jpg'];
                
            end
            
        else
        
            filename = ['Variant_PieChart_All_Rest_' num2str(threshold) 'pct_AllSubjects.jpg'];
            
        end
        
        pie(restnetnums,network_names)
        title('All Rest Variants', 'fontsize',18)
        
        %saveas(gcf,[outputdir2 filename])
        
        print(gcf,[outputdir2 filename],'-dpng','-r300');
        
        close gcf
        
        if SizeExclude == 1
            
            if matchbydist == 0 && excludemaxmatch == 1
                
                filename = ['Variant_PieChart_All_Task_' num2str(threshold) 'pct_SizeExclude_AllSubjects_OverlapOnly_ExcludeMatch.jpg'];
            
            elseif matchbydist == 0
                
                filename = ['Variant_PieChart_All_Task_' num2str(threshold) 'pct_SizeExclude_AllSubjects_OverlapOnly.jpg'];
                
            else
            
                filename = ['Variant_PieChart_All_Task_' num2str(threshold) 'pct_SizeExclude_AllSubjects.jpg'];
                
            end
            
        else
        
            filename = ['Variant_PieChart_All_Task_' num2str(threshold) 'pct_AllSubjects.jpg'];
            
        end
        
        pie(tasknetnums,network_names)
        title('All Task Variants', 'fontsize',18)
        
        %saveas(gcf,[outputdir2 filename])
        
        print(gcf,[outputdir2 filename],'-dpng','-r300');
        
        close gcf
        
    end
    
     if getmatchbystate == 1
         
         netmatch = {'DMN','Vis','FP','DAN','VAN','Sal','CO','SMd','SMl','Aud','Tpole','MTL','PMN','PON','Unmatched'};
         
         netplots = {'DMN','FP','DAN','VAN','Sal','CO'};
         
         for nets = 1:length(netplots)
             
             plotnets = {};
             
             for y = 1:length(matchedvarsnetsrest)
                 
                 if strcmp(matchedvarsnetsrest(y),netplots(nets))
                     
                     if ismember(matchedvarsnetstask(y),netmatch)
                         
                         matchednet = netmatch(ismember(netmatch,matchedvarsnetstask(y)));
                         
                     else
                         
                         matchednet = 'Unmatched';
                         
                     end
                     
                     plotnets = [plotnets; matchednet];
                     
                 end
             end
             
             outputdir2 = '/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Variant_Network_Matching/';
             
             cellsum = cellfun(@(x) sum(ismember(plotnets,x)),netmatch,'un',0);
             plot(1:length(cellsum), cell2mat(cellsum))
             set(gca,'xtick',[1:16],'xticklabel',netmatch)
             title(['Task Variants - Rest assigned as ' netplots(nets)]);

             if SizeExclude == 1
                 
                 if matchbydist == 0 && excludemaxmatch == 1
                     
                     filename = ['AllSubjects_RestToTask_Variants_' char(netplots(nets)) '_Network_' num2str(threshold) 'pct_SizeExclude_OverlapOnly_ExcludeMatch.jpg'];
                 
                 elseif matchbydist == 0
                     
                     filename = ['AllSubjects_RestToTask_Variants_' char(netplots(nets)) '_Network_' num2str(threshold) 'pct_SizeExclude_OverlapOnly.jpg'];
                    
                 else
                     
                   	filename = ['AllSubjects_RestToTask_Variants_' char(netplots(nets)) '_Network_' num2str(threshold) 'pct_SizeExclude.jpg'];
             
                 end
                 
             else
                 
              	filename = ['AllSubjects_RestToTask_Variants_' char(netplots(nets)) '_Network_' num2str(threshold) 'pct.jpg'];
                
             end
             
           	 %saveas(gcf,[outputdir2 filename])
             
             print(gcf,[outputdir2 filename],'-dpng','-r300');
                     
          	 close gcf
   
         end
         
         for nets = 1:length(netplots)
             
             plotnets = {};
             
             for y = 1:length(matchedvarsnetstask)
                 
                 if strcmp(matchedvarsnetstask(y),netplots(nets))
                     
                     if ismember(matchedvarsnetsrest(y),netmatch)
                         
                         matchednet = netmatch(ismember(netmatch,matchedvarsnetsrest(y)));
                         
                     else
                         
                         matchednet = 'Unmatched';
                         
                     end
                     
                     plotnets = [plotnets; matchednet];
                     
                 end
             end
             
             outputdir2 = '/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Variant_Network_Matching/';
             
             cellsum = cellfun(@(x) sum(ismember(plotnets,x)),netmatch,'un',0);
             plot(1:length(cellsum), cell2mat(cellsum))
             set(gca,'xtick',[1:16],'xticklabel',netmatch)
             title(['Rest Variants - Task assigned as ' netplots(nets)]);
             
             if SizeExclude == 1
                 
                if matchbydist == 0 && excludemaxmatch == 1
                    
                    filename = ['AllSubjects_TaskToRest_Variants_' char(netplots(nets)) '_Network_' num2str(threshold) 'pct_SizeExclude_OverlapOnly_ExcludeMatch.jpg'];
                 
                elseif matchbydist == 0
                    
                    filename = ['AllSubjects_TaskToRest_Variants_' char(netplots(nets)) '_Network_' num2str(threshold) 'pct_SizeExclude_OverlapOnly.jpg'];

                else
                    
                    filename = ['AllSubjects_TaskToRest_Variants_' char(netplots(nets)) '_Network_' num2str(threshold) 'pct_SizeExclude.jpg'];
                    
                end
             
             else
                 
              	filename = ['AllSubjects_TaskToRest_Variants_' char(netplots(nets)) '_Network_' num2str(threshold) 'pct.jpg'];
                
             end
             
             print(gcf,[outputdir2 filename],'-dpng','-r300');
             
           	 %saveas(gcf,[outputdir2 filename])
                     
          	 close gcf    
             
         end
     end
    
end
    

if savestruct == 1
    
    outputstructpath = '/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Variant_Network_Matching/';
    filename2 = 'Matched_Variant_Struct';
    
    save([outputstructpath filename2], 'Variant_Match_Struct');
    
end



    
    