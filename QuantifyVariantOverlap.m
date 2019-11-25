
clear all

DiceCorr = 1;  %% Toggles whether to calculate dice correlations for all subjects
AnyOverlap = 0;  %% Toggles whether to calculate if any portion of a variant overlaps for all subjects
COMOverlap = 0;  %% Toggles whether to calculate if the center of mass of each variant shows overlap
plotresults = 1;  %% Toggles whether to plot results for each comparison
FullMaps = 0;  %% Toggles whether to use all task and rest data for each subject
SplitHalf = 1;  %% Toggles whether to calculate a split-half 
MatchedMaps = 1;  %% Toggles whether to use a matched amount of task and rest data for each subject
SNRExclude = 1;  %% Toggles whether to use data excluded for low SNR
SizeExclude = 1;  %% Toggles whether to exclude small variants
randomizevals = 1;  %% Toggles whether to calculate a null distribution across subjects
permvals = 1; %% Toggles whether to calculate a permutation of all possible combinations instead of bootstrapping
BarGraph = 0;  %% Toggles whether to plot a bar graph for the group level results
BoxPlot = 1;  %% Toggles whether to plot a box plot for the group level results
BoxPlotBySubject = 1;  %% Toggles whether to display boxplot by subjects instead of by comparison
FinalFigure = 1;  %% Toggles whether to plot the final figure (Figure 1) for the paper
SubjectPlots = 0;  %% Toggles whether to plot subject level results


AbsoluteThresholds = 0;  %% Toggles whether thresholds are absolute values or percents
%%
%%Do we want to leave 10 as a threshold?
thresholds = [2.5 10];  %% Sets the thresholds of files to load
%%
% Specify output directory    
outputdir = '/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Variant_Overlap_Plots/Thresholded_Plots';
%outputdirsub = '/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Variant_Overlap_Plots/Thresholded_Plots/Subject_Plots';
   
% Create variables for number of variants/average size
    
    NumVariantsTaskOdd = [];
    NumVariantsTaskEven = [];
    NumVariantsRestOdd = [];
    NumVariantsRestEven = [];
    MeanSizeVariantsTaskOdd = [];
    MeanSizeVariantsTaskEven = [];
    MeanSizeVariantsRestOdd = [];
    MeanSizeVariantsRestEven = [];
    MedianSizeVariantsTaskOdd = [];
    MedianSizeVariantsTaskEven = [];
    MedianSizeVariantsRestOdd = [];
    MedianSizeVariantsRestEven = [];

% this for loop sets up variables to run dice correlations between task and
% rest data for each threshold
for v = 1:numel(thresholds)
    
%% Create temp variables for number of variants/average size
    
        NumVariantsTaskOddTemp = [];
        NumVariantsTaskEvenTemp = [];
        NumVariantsRestOddTemp = [];
        NumVariantsRestEvenTemp = [];
        MeanSizeVariantsTaskOddTemp = [];
        MeanSizeVariantsTaskEvenTemp = [];
        MeanSizeVariantsRestOddTemp = [];
        MeanSizeVariantsRestEvenTemp = [];
        MedianSizeVariantsTaskOddTemp = [];
        MedianSizeVariantsTaskEvenTemp = [];
        MedianSizeVariantsRestOddTemp = [];
        MedianSizeVariantsRestEvenTemp = [];

        %what are these?
        alltaskfilestaskeven = [];
        alltaskfilestaskodd = [];
        alltaskfilesresteven = [];
        alltaskfilesrestodd = [];
        
       

        % specify .txt files 
        [task_files_even, ~, ~] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_' num2str(thresholds(v)) '_splithalf_even_matched_variants_SNRexclude.txt'],'%s%s%s');

        [rest_files_even, ~, ~] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_' num2str(thresholds(v)) '_splithalf_even_matched_variants_SNRexclude.txt'],'%s%s%s');
    
        [task_files_odd, subjects1, tasks1] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_' num2str(thresholds(v)) '_splithalf_odd_matched_variants_SNRexclude.txt'],'%s%s%s');

        [rest_files_odd, subjects2, tasks2] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_' num2str(thresholds(v)) '_splithalf_odd_matched_variants_SNRexclude.txt'],'%s%s%s');
    
        [task_masks_even, ~, ~] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_' num2str(thresholds(v)) '_splithalf_even_matched_variants_SNRexclude_sizeexclude.txt'],'%s%s%s');

        [rest_masks_even, ~, ~] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_' num2str(thresholds(v)) '_splithalf_even_matched_variants_SNRexclude_sizeexclude.txt'],'%s%s%s');
    
        [task_masks_odd, sub1, t1] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_' num2str(thresholds(v)) '_splithalf_odd_matched_variants_SNRexclude_sizeexclude.txt'],'%s%s%s');

        [rest_masks_odd, sub2, t2] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_' num2str(thresholds(v)) '_splithalf_odd_matched_variants_SNRexclude_sizeexclude.txt'],'%s%s%s');
        
 

        % create temp variables for dice correlations
        DiceCorrsTaskRest = [];
        DiceCorrsTaskTask = [];
        DiceCorrsRestRest = [];
        
        
        
       
    end

  

% sets up number of files??
nfiles = length(rest_files_even);
    
  

    for x = 1:nfiles

        subject = subjects2{x};
   
        % reads cifti files from .txt file paths
        cifti_rest_even = ft_read_cifti_mod(rest_files_even{x});
        cifti_task_even = ft_read_cifti_mod(task_files_even{x});
        cifti_rest_odd = ft_read_cifti_mod(rest_files_odd{x});
        cifti_task_odd = ft_read_cifti_mod(task_files_odd{x});
      
         %% Apply exclusion masks for maps that exclude size
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
            
    
  %% Get number of variants and mean/median size
        
            vars_rest_even = unique(cifti_rest_even.data);
            vars_rest_even(1) = [];
            vars_rest_odd = unique(cifti_rest_odd.data);
            vars_rest_odd(1) = [];
            vars_task_even = unique(cifti_task_even.data);
            vars_task_even(1) = [];
            vars_task_odd = unique(cifti_task_odd.data);
            vars_task_odd(1) = [];
            
            NumVariantsTaskOddTemp = [NumVariantsTaskOddTemp; length(vars_task_odd)];
            NumVariantsTaskEvenTemp = [NumVariantsTaskEvenTemp; length(vars_task_even)];
            NumVariantsRestOddTemp = [NumVariantsRestOddTemp; length(vars_rest_odd)];
            NumVariantsRestEvenTemp = [NumVariantsRestEvenTemp; length(vars_rest_even)];
            MeanSizeVariantsTaskOddTemp = [MeanSizeVariantsTaskOddTemp; length(find(cifti_task_odd.data > 0))/length(vars_task_odd)];
            MeanSizeVariantsTaskEvenTemp = [MeanSizeVariantsTaskEvenTemp; length(find(cifti_task_even.data > 0))/length(vars_task_even)];
            MeanSizeVariantsRestOddTemp = [MeanSizeVariantsRestOddTemp; length(find(cifti_rest_odd.data > 0))/length(vars_rest_odd)];
            MeanSizeVariantsRestEvenTemp = [MeanSizeVariantsRestEvenTemp; length(find(cifti_rest_even.data > 0))/length(vars_rest_even)];
            
            vars_rest_size_even = zeros(length(vars_rest_even),1);
            vars_rest_size_odd = zeros(length(vars_rest_odd),1);
            vars_task_size_even = zeros(length(vars_task_even),1);
            vars_task_size_odd = zeros(length(vars_task_odd),1);
            
            for k = 1:length(vars_rest_even)
            
                vars_rest_size_even(k) = length(find(cifti_rest_even.data == vars_rest_even(k)));
                
            end
            
            for k = 1:length(vars_rest_odd)
            
                vars_rest_size_odd(k) = length(find(cifti_rest_odd.data == vars_rest_odd(k)));
                
            end
            
            for k = 1:length(vars_task_even)
            
                vars_task_size_even(k) = length(find(cifti_task_even.data == vars_task_even(k)));
                
            end
            
            for k = 1:length(vars_task_odd)
            
                vars_task_size_odd(k) = length(find(cifti_task_odd.data == vars_task_odd(k)));
                
            end
            
            MedianSizeVariantsTaskOddTemp = [MedianSizeVariantsTaskOddTemp; median(vars_task_size_odd)];
            MedianSizeVariantsTaskEvenTemp = [MedianSizeVariantsTaskEvenTemp; median(vars_task_size_even)];
            MedianSizeVariantsRestOddTemp = [MedianSizeVariantsRestOddTemp; median(vars_rest_size_odd)];
            MedianSizeVariantsRestEvenTemp = [MedianSizeVariantsRestEvenTemp; median(vars_rest_size_even)];
        
       
            %% dice correlations
                dcorrdatataskrest1 = [];
                dcorrdatataskrest2 = [];
                dcorrdatarestrest = [];
                dcorrdatatasktask = [];

                for q = 1:length(cifti_rest_even.data)      %% Task-Rest Comparison 1
            
                    if cifti_rest_even.data(q) > 0 && cifti_task_odd.data(q) > 0
                
                        dcorrdatataskrest1 = [dcorrdatataskrest1;1 1];
                
                    elseif cifti_rest_even.data(q) > 0
                
                        dcorrdatataskrest1 = [dcorrdatataskrest1;1 0];
                
                    elseif cifti_task_odd.data(q) > 0
                
                        dcorrdatataskrest1 = [dcorrdatataskrest1;0 1];
                
                    end
            
                end
                
                for q = 1:length(cifti_rest_odd.data)      %% Task-Rest Comparison 2
            
                    if cifti_rest_odd.data(q) > 0 && cifti_task_even.data(q) > 0
                
                        dcorrdatataskrest2 = [dcorrdatataskrest2;1 1];
                
                    elseif cifti_rest_even.data(q) > 0
                
                        dcorrdatataskrest2 = [dcorrdatataskrest2;1 0];
                
                    elseif cifti_task_odd.data(q) > 0
                
                        dcorrdatataskrest2 = [dcorrdatataskrest2;0 1];
                
                    end
            
                end
            
                for q = 1:length(cifti_rest_odd.data)       %% Rest-Rest Comparison
            
                    if cifti_rest_odd.data(q) > 0 && cifti_rest_even.data(q) > 0
                
                        dcorrdatarestrest = [dcorrdatarestrest;1 1];
                
                    elseif cifti_rest_odd.data(q) > 0
                
                        dcorrdatarestrest = [dcorrdatarestrest;1 0];
                
                    elseif cifti_rest_even.data(q) > 0
                
                        dcorrdatarestrest = [dcorrdatarestrest;0 1];
                
                    end
            
                end
            
                for q = 1:length(cifti_task_odd.data)       %% Task-Task Comparison            
                    if cifti_task_odd.data(q) > 0 && cifti_task_even.data(q) > 0
                        dcorrdatatasktask = [dcorrdatatasktask;1 1];
                    elseif cifti_task_odd.data(q) > 0
                        dcorrdatatasktask = [dcorrdatatasktask;1 0];
                    elseif cifti_task_even.data(q) > 0
                        dcorrdatatasktask = [dcorrdatatasktask;0 1];
                    end
            
                end
                
                if isempty(dcorrdatataskrest1)
                    dctaskrest1 = 0;
                else                   
                    dctaskrest1 = dice_coefficient_mod(dcorrdatataskrest1(:,1),dcorrdatataskrest1(:,2));                    
                end
                
                if isempty(dcorrdatataskrest2)                    
                    dctaskrest2 = 0;
                else
                    dctaskrest2 = dice_coefficient_mod(dcorrdatataskrest2(:,1),dcorrdatataskrest2(:,2));
                end
                
                if isempty(dcorrdatarestrest)
                    dcrestrest = 0;
                else
                    dcrestrest = dice_coefficient_mod(dcorrdatarestrest(:,1),dcorrdatarestrest(:,2));
                end
                
                if isempty(dcorrdatatasktask)
                    dctasktask = 0;
                else
                    dctasktask = dice_coefficient_mod(dcorrdatatasktask(:,1),dcorrdatatasktask(:,2));
                end
        
                DiceCorrsTaskRest = [DiceCorrsTaskRest; [dctaskrest1 dctaskrest2]];
                DiceCorrsTaskTask = [DiceCorrsTaskTask; dctasktask];
                DiceCorrsRestRest = [DiceCorrsRestRest; dcrestrest];
            

%         if AnyOverlap == 1
%         
%             if SplitHalf == 1
%             
%                 overlapvarstaskrest = [];
%                 overlapvarstasktask = [];
%                 overlapvarsrestrest = [];
%         
%                 for r = 1:length(vars_rest_even)     %% Task-Rest Comparison
%             
%                     for s = 1:length(cifti_rest_even.data)
%             
%                         if cifti_rest_even.data(s) == vars_rest_even(r) && cifti_task_odd.data(s) > 0
%                 
%                             overlapvarstaskrest = [overlapvarstaskrest;vars_rest_even(r) cifti_task_odd.data(s)];
%                 
%                         end
%                 
%                     end
%                 end
%             
%                 for r = 1:length(vars_task_odd)        %% Task-Task Comparison
%             
%                     for s = 1:length(cifti_task_odd.data)
%             
%                         if cifti_task_odd.data(s) == vars_task_odd(r) && cifti_task_even.data(s) > 0
%                 
%                             overlapvarstasktask = [overlapvarstasktask;vars_task_odd(r) cifti_task_even.data(s)];
%                 
%                         end
%                 
%                     end
%                 end
%             
%                 for r = 1:length(vars_rest_odd)        %% Rest-Rest Comparison
%             
%                     for s = 1:length(cifti_rest_odd.data)
%             
%                         if cifti_rest_odd.data(s) == vars_rest_odd(r) && cifti_rest_even.data(s) > 0
%                 
%                             overlapvarsrestrest = [overlapvarsrestrest;vars_rest_odd(r) cifti_rest_even.data(s)];
%                 
%                         end
%                 
%                     end
%                 end
%                 
%                 if isempty(overlapvarstasktask) && isempty(overlapvarsrestrest) && isempty(overlapvarstaskrest)
%                         
%                     NumOverlapTaskRest = [NumOverlapTaskRest; 0];
%                     NumOverlapTaskTask = [NumOverlapTaskTask; 0];
%                     NumOverlapRestRest = [NumOverlapRestRest; 0];
%                         
%                 elseif isempty(overlapvarsrestrest) && isempty(overlapvarstaskrest)
%                         
%                     NumOverlapTaskRest = [NumOverlapTaskRest; 0];
%                     NumOverlapRestRest = [NumOverlapRestRest; 0];
%                     NumOverlapTaskTask = [NumOverlapTaskTask; (length(unique(overlapvarstasktask(:,1))) + length(unique(overlapvarstasktask(:,2))))/(length(vars_task_odd) + length(vars_task_even))];
%                         
%                 elseif isempty(overlapvarstasktask) && isempty(overlapvarstaskrest)
%                         
%                  	NumOverlapTaskRest = [NumOverlapTaskRest; 0];
%                     NumOverlapTaskTask = [NumOverlapTaskTask; 0];
%                     NumOverlapRestRest = [NumOverlapRestRest; (length(unique(overlapvarsrestrest(:,1))) + length(unique(overlapvarsrestrest(:,2))))/(length(vars_rest_odd) + length(vars_rest_even))];
%                   
%                 elseif isempty(overlapvarsrestrest) && isempty(overlapvarstasktask)
%                     
%                     NumOverlapTaskTask = [NumOverlapTaskTask; 0];
%                     NumOverlapRestRest = [NumOverlapRestRest; 0];
%                     NumOverlapTaskRest = [NumOverlapTaskRest; (length(unique(overlapvarstaskrest(:,1))) + length(unique(overlapvarstaskrest(:,2))))/(length(vars_rest_odd) + length(vars_task_even))];
%                     
%                 elseif isempty(overlapvarsrestrest)
%                     
%                     NumOverlapRestRest = [NumOverlapRestRest; 0];
%                     NumOverlapTaskTask = [NumOverlapTaskTask; (length(unique(overlapvarstasktask(:,1))) + length(unique(overlapvarstasktask(:,2))))/(length(vars_task_odd) + length(vars_task_even))];
%                     NumOverlapTaskRest = [NumOverlapTaskRest; (length(unique(overlapvarstaskrest(:,1))) + length(unique(overlapvarstaskrest(:,2))))/(length(vars_rest_odd) + length(vars_task_even))];
%                     
%                 elseif isempty(overlapvarstasktask)
%                     
%                     NumOverlapTaskTask = [NumOverlapTaskTask; 0];
%                     NumOverlapTaskRest = [NumOverlapTaskRest; (length(unique(overlapvarstaskrest(:,1))) + length(unique(overlapvarstaskrest(:,2))))/(length(vars_rest_odd) + length(vars_task_even))];
%                     NumOverlapRestRest = [NumOverlapRestRest; (length(unique(overlapvarsrestrest(:,1))) + length(unique(overlapvarsrestrest(:,2))))/(length(vars_rest_odd) + length(vars_rest_even))];
%                     
%                 elseif isempty(overlapvarstaskrest)
%                     
%                     NumOverlapTaskRest = [NumOverlapTaskRest; 0];
%                     NumOverlapTaskTask = [NumOverlapTaskTask; (length(unique(overlapvarstasktask(:,1))) + length(unique(overlapvarstasktask(:,2))))/(length(vars_task_odd) + length(vars_task_even))];
%                     NumOverlapRestRest = [NumOverlapRestRest; (length(unique(overlapvarsrestrest(:,1))) + length(unique(overlapvarsrestrest(:,2))))/(length(vars_rest_odd) + length(vars_rest_even))];
%                     
%                 else
%         
%                     NumOverlapTaskRest = [NumOverlapTaskRest; (length(unique(overlapvarstaskrest(:,1))) + length(unique(overlapvarstaskrest(:,2))))/(length(vars_rest_odd) + length(vars_task_even))];
%                     NumOverlapTaskTask = [NumOverlapTaskTask; (length(unique(overlapvarstasktask(:,1))) + length(unique(overlapvarstasktask(:,2))))/(length(vars_task_odd) + length(vars_task_even))];
%                     NumOverlapRestRest = [NumOverlapRestRest; (length(unique(overlapvarsrestrest(:,1))) + length(unique(overlapvarsrestrest(:,2))))/(length(vars_rest_odd) + length(vars_rest_even))];
%                         
%                 end
%             
%             else
   
                overlapvars = [];
        
                for r = 1:length(vars_rest)
            
                    for s = 1:length(cifti_rest.data)
            
                        if cifti_rest.data(s) == vars_rest(r) && cifti_task.data(s) > 0
                
                            overlapvars = [overlapvars;vars_rest(r) cifti_task.data(s)];
                
                        end
                
                    end
                end
            
                NumOverlap = [NumOverlap; (length(unique(overlapvars(:,1))) + length(unique(overlapvars(:,2))))/(length(vars_rest) + length(vars_task))];
            
%             end
        end

%         if COMOverlap == 1
%         
%             if SplitHalf == 1
%             
%                 for z = 1:3     %% Loops over 3 comparisons (rest-task, rest-rest, task-task)
%             
%                     %meanvals = [];
%                     COMOverlapsresttemp = [];
%                     COMOverlapstasktemp = [];
%                 
%                     if z == 1  %% Task-Rest
%         
%                         taskverts = find(cifti_task_odd.data > 0);
%                         restverts = find(cifti_rest_even.data > 0);
%                     
%                     elseif z == 2  %% Task-Task
%                     
%                         taskverts = find(cifti_task_even.data > 0);
%                         restverts = find(cifti_task_odd.data > 0);
%                     
%                     elseif z == 3  %% Rest-Rest
%                     
%                         taskverts = find(cifti_rest_even.data > 0);
%                         restverts = find(cifti_rest_odd.data > 0);
%                     
%                     end
%         
%                     for a = 1:2   %% loop over task and rest
%             
%                         if a == 1   %% Rest
%                         
%                             if z == 1   %% Task-Rest
%                 
%                                 nvars = length(vars_rest_even);
%                         
%                             elseif z == 2   %% Task-Task
%                             
%                                 nvars = length(vars_task_odd);
%                             
%                             elseif z == 3   %% Rest-Rest
%                             
%                                 nvars = length(vars_rest_odd);
%                             
%                             end
%                 
%                         else       %% Task
%                 
%                             if z == 1   %% Task-Rest
%                 
%                                 nvars = length(vars_task_odd);
%                         
%                             elseif z == 2   %% Task-Task
%                             
%                                 nvars = length(vars_task_even);
%                             
%                             elseif z == 3   %% Rest-Rest
%                             
%                                 nvars = length(vars_rest_even);
%                             
%                             end
%                 
%                         end
%         
%                         for s = 1:nvars
%                 
%                             if a == 1   %% Rest
%                             
%                                 if z == 1   %% Task-Rest
%                 
%                                     currentvariant = find(cifti_rest_even.data == vars_rest_even(s));
%                         
%                                 elseif z == 2   %% Task-Task
% 
%                                     currentvariant = find(cifti_task_odd.data == vars_task_odd(s));
%                             
%                                 elseif z == 3   %% Rest-Rest
% 
%                                     currentvariant = find(cifti_rest_odd.data == vars_rest_odd(s));
%                             
%                                 end
%                 
%                             else   %% Task
%                     
%                                 if z == 1   %% Task-Rest
%                 
%                                     currentvariant = find(cifti_task_odd.data == vars_task_odd(s));
%                         
%                                 elseif z == 2   %% Task-Task
% 
%                                     currentvariant = find(cifti_task_even.data == vars_task_even(s));
%                             
%                                 elseif z == 3   %% Rest-Rest
% 
%                                     currentvariant = find(cifti_rest_even.data == vars_rest_even(s));
%                             
%                                 end
%                     
%                             end
%                 
%                             COMtemp = cifti_coords.data(currentvariant,:);
%                 
%                             if length(currentvariant) > 1
%                     
%                                 centroid = [];
%                     
%                                 for t = 1:length(currentvariant)
%                         
%                                     jacknifevals = 1:length(currentvariant);        %% Leaves current point out of calculation
%                                     jacknifevals(t) = [];
%                     
%                                     % Euclidean distance
%                                     dist = mean(sqrt(sum((COMtemp(jacknifevals,:) - COMtemp(t,:)).^2, 2)));
%                         
%                                     %meanvals = [meanvals; dist];
%                         
%                                     if isempty(centroid) || dist < centroid   %% If new distance less than current least distance, use as center
%                             
%                                         centroid = currentvariant(t);
%                             
%                                     end
%                         
%                                 end
%                     
%                                 centervariant = centroid;
%                     
%                             else
%                     
%                                 centervariant = currentvariant;
%                     
%                             end
%                 
%                             if a == 1
%                 
%                                 COMOverlapsresttemp = [COMOverlapsresttemp; centervariant];
%                     
%                             else
%                     
%                                 COMOverlapstasktemp = [COMOverlapstasktemp; centervariant];
%                     
%                             end 
%                         end
%                     end
% 
%         
%                     Taskoverlaps = length(unique(intersect(COMOverlapstasktemp, restverts)));
%                     Restoverlaps = length(unique(intersect(COMOverlapsresttemp, taskverts)));
%                 
%                     if z == 1   %% Task-Rest
%         
%                         COMOverlapsTaskRest = [COMOverlapsTaskRest; (Taskoverlaps + Restoverlaps)/(length(vars_rest_even) + length(vars_task_odd))];
%                     
%                     elseif z == 2  %% Task-Task
%                     
%                         COMOverlapsTaskTask = [COMOverlapsTaskTask; (Taskoverlaps + Restoverlaps)/(length(vars_task_even) + length(vars_task_odd))];
%                     
%                     elseif z == 3  %% Rest-Rest
%                     
%                         COMOverlapsRestRest = [COMOverlapsRestRest; (Taskoverlaps + Restoverlaps)/(length(vars_rest_even) + length(vars_rest_odd))];
%                     
%                     end
%                 
%                 end
%             
%             else
%     
%                 %meanvals = [];
%                 COMOverlapsresttemp = [];
%                 COMOverlapstasktemp = [];
% 
%                 taskverts = find(cifti_task.data > 0);
%                 restverts = find(cifti_rest.data > 0);
%         
%                 for a = 1:2   %% loop over task and rest
%             
%                     if a == 1
%                 
%                         nvars = length(vars_rest);
%                 
%                     else
%                 
%                         nvars = length(vars_task);
%                 
%                     end
%         
%                     for s = 1:nvars
%                 
%                         if a == 1
%             
%                             currentvariant = find(cifti_rest.data == vars_rest(s));
%                 
%                         else
%                     
%                             currentvariant = find(cifti_task.data == vars_task(s));
%                     
%                         end
%                 
%                         COMtemp = cifti_coords.data(currentvariant,:);
%                 
%                         if length(currentvariant) > 1
%                     
%                             centroid = [];
%                     
%                             for t = 1:length(currentvariant)
%                         
%                                 jacknifevals = 1:length(currentvariant);        %% Leaves current point out of calculation
%                                 jacknifevals(t) = [];
%                     
%                                 % Euclidean distance
%                                 dist = mean(sqrt(sum((COMtemp(jacknifevals,:) - COMtemp(t,:)).^2, 2)));
%                         
%                                 %meanvals = [meanvals; dist];
%                         
%                                 if isempty(centroid) || dist < centroid   %% If new distance less than current least distance, use as center
%                             
%                                     centroid = currentvariant(t);
%                             
%                                 end
%                         
%                             end
%                     
%                             centervariant = centroid;
%                     
%                         else
%                     
%                             centervariant = currentvariant;
%                     
%                         end
%                 
%                         if a == 1
%                 
%                             COMOverlapsresttemp = [COMOverlapsresttemp; centervariant];
%                     
%                         else
%                     
%                             COMOverlapstasktemp = [COMOverlapstasktemp; centervariant];
%                     
%                         end 
%                     end
%                 end
% 
%         
%                 Taskoverlaps = length(unique(intersect(COMOverlapstasktemp, restverts)));
%                 Restoverlaps = length(unique(intersect(COMOverlapsresttemp, taskverts)));
%         
%                 COMOverlaps = [COMOverlaps; (Taskoverlaps + Restoverlaps)/(length(vars_rest) + length(vars_task))];
%         
%             end
        
        
        % Find Euclidean distance between task and rest centroids
        
%         minval = [];
%         
%         for y = 1:length(COMOverlapsresttemp)
%             
%             % Euclidean distance
%           	dist = sqrt(sum((cifti_coords.data(COMOverlapstasktemp,:) - cifti_coords.data(COMOverlapsresttemp(y),:)).^2, 2));
%            
%             minval = [minval dist];
%             
%         end
        
        % Exclude centroids from task/rest variants that are overlapping
%         
%         taskverts = find(cifti_task.data > 0);
%         
%         taskdata = cifti_task.data;
%         restdata = cifti_rest.data;
%         
%         restvertsexclude = [];
%         taskvertsexclude = [];
% 
%         for z = 1:length(COMOverlapsresttemp)
%             
%             if ismember(COMOverlapsresttemp(z), taskverts)
%                 
%                 taskvarianttemp = taskdata(COMOverlapsresttemp(z));
%                 restvarianttemp = restdata(COMOverlapsresttemp(z));
%                 
%                 taskvarcount = length(find(taskdata == taskvarianttemp));
%                 restvarcount = length(find(restdata == restvarianttemp));
%                 
%                 if taskvarcount < restvarcount
%                     
%                     taskvertsexclude = [taskvertsexclude; find(taskdata == taskvarianttemp)];
%                     
%                 elseif taskvarcount > restvarcount
%                     
%                     restvertsexclude = [restvertsexclude; find(restdata == restvarianttemp)];
%                 end 
%             end
%         end
%         
%         for b = 1:length(COMOverlapsresttemp)
%             
%             ExcludeCOMrest = [];
%             
%             if ismember(COMOverlapsresttemp(b), restvertsexclude)
%                 
%                 ExcludeCOMrest = [ExcludeCOMrest; COMOverlapsresttemp(b)];
%                 
%             end
%         end        
%                 
%         COMOverlapsresttemp(ExcludeCOM,:) = [];
%         taskdata(taskvertsexclude,:) = [];
%         restdata(restvertsexclude,:) = [];
%         
%         finaltaskvarcount = length(unique(find(taskdata > 0)));
%         finalrestvarcount = length(unique(find(restdata > 0)));
    
%         end
        

%what are these variables?
      	if randomizevals == 1 || SubjectPlots == 1
            
         	if SplitHalf == 1
                
              	alltaskfilestaskeven = [alltaskfilestaskeven; cifti_task_even.data'];
              	alltaskfilestaskodd = [alltaskfilestaskodd; cifti_task_odd.data'];
              	alltaskfilesresteven = [alltaskfilesresteven; cifti_rest_even.data'];
             	alltaskfilesrestodd = [alltaskfilesrestodd; cifti_rest_odd.data'];
                
            else
                        
             	alltaskfilestask = [alltaskfilestask; cifti_task.data'];
              	alltaskfilesrest = [alltaskfilesrest; cifti_rest.data']; 
            
            end
                        
      	end

%     end
    
    
     %% Add temp variables to final variables
    
        NumVariantsTaskOdd = [NumVariantsTaskOdd NumVariantsTaskOddTemp];
        NumVariantsTaskEven = [NumVariantsTaskEven NumVariantsTaskEvenTemp];
        NumVariantsRestOdd = [NumVariantsRestOdd NumVariantsRestOddTemp];
        NumVariantsRestEven = [NumVariantsRestEven NumVariantsRestEvenTemp];
        MeanSizeVariantsTaskOdd = [MeanSizeVariantsTaskOdd MeanSizeVariantsTaskOddTemp];
        MeanSizeVariantsTaskEven = [MeanSizeVariantsTaskEven MeanSizeVariantsTaskEvenTemp];
        MeanSizeVariantsRestOdd = [MeanSizeVariantsRestOdd MeanSizeVariantsRestOddTemp];
        MeanSizeVariantsRestEven = [MeanSizeVariantsRestEven MeanSizeVariantsRestEvenTemp];
        MedianSizeVariantsTaskOdd = [MedianSizeVariantsTaskOdd MedianSizeVariantsTaskOddTemp];
        MedianSizeVariantsTaskEven = [MedianSizeVariantsTaskEven MedianSizeVariantsTaskEvenTemp];
        MedianSizeVariantsRestOdd = [MedianSizeVariantsRestOdd MedianSizeVariantsRestOddTemp];
        MedianSizeVariantsRestEven = [MedianSizeVariantsRestEven MedianSizeVariantsRestEvenTemp];
    

%     if randomizevals == 1
%         
%         if permvals == 1
%             
            nsubs = nfiles;
            
            DiceCorrssim = [];
            DiceCorrSubs = zeros(2*(nsubs-1),nsubs);

            for l = 1:2     %% Do both versions of split-half combinations
                
                for m = 1:nsubs
                    
                    loopvals = [1:nsubs];
                    loopvals(m) = [];
                    
                    Count = 0;
                    
                    for n = loopvals
                        
                        Count = Count+1;

                        if l == 1
                            
                            restdat = alltaskfilesrestodd(m,:);
                            taskdat = alltaskfilestaskeven(n,:);
                            
                        else
                            
                            restdat = alltaskfilesresteven(m,:);
                            taskdat = alltaskfilestaskodd(n,:);
                            
                        end
                            
                    	dcorrdata = [];
        
                     	for q = 1:length(restdat)
            
                         	if restdat(:,q) > 0 && taskdat(:,q) > 0
                
                             	dcorrdata = [dcorrdata;1 1];
                
                            elseif restdat(:,q) > 0
                
                                dcorrdata = [dcorrdata;1 0];
                
                            elseif taskdat(:,q) > 0
                
                                dcorrdata = [dcorrdata;0 1];
                
                            end
            
                        end
                    
                        if isempty(dcorrdata)
                        
                            dc = 0;
                        
                        else
        
                            dc = dice_coefficient_mod(dcorrdata(:,1),dcorrdata(:,2));
                        
                        end
        
                        DiceCorrssim = [DiceCorrssim; dc];
                        
                        if l == 1
                            
                            DiceCorrSubs(Count,m) = dc;
                            
                        else
                            
                            DiceCorrSubs(Count+(nsubs-1),m) = dc;
                            
                        end
                        
                        
                    end
                end
            end
            
            DiceCorrspval = prctile(DiceCorrssim, 95);
            
%         else
%     
%             nsim = 1000;
%             nsubs = nfiles;
%     
%             DiceCorrssim = zeros(1,nsim);
%         
%             if AnyOverlap == 1
%         
%                 NumOverlapsim = zeros(1,nsim);
%                 
%             end
%         
%             if COMOverlap == 1
%             
%                 COMOverlapssim = zeros(1,nsim);
%             
%             end
%     
%             DiceCorrssimAll = zeros(nsubs,nsim);
%     
%             tasktaskshuffleall = [];
%             taskrestshuffleall = [];
%     
%     
%             for x = 1:nsim
%         
% 
%                 alltaskdc = [];
%         
%                 while true      %% Make sure no subjects match each other
%             
%                     reshuffle = 0;
%     
%                     rng('shuffle');  %% Reset random number generator on each iteration
%         
%                     % Randomly shuffle rows (subjects) for each task and rest matrix
%         
%                     tasktaskshuffle = 1:nsubs;
%                     taskrestshuffle = 1:nsubs;
%         
%                     tasktaskshuffle = tasktaskshuffle(randperm(length(1:nsubs)));
%                     taskrestshuffle = taskrestshuffle(randperm(length(1:nsubs)));
%             
%                     for h = 1:length(tasktaskshuffle)
%                 
%                         if reshuffle == 0 && tasktaskshuffle(h) == taskrestshuffle(h)
%                     
%                             reshuffle = 1;
%                     
%                         end
%                     end
%             
%                     if reshuffle == 0  %% Break loop if all subjects are shuffled
%                 
%                         tasktaskshuffleall = [tasktaskshuffleall tasktaskshuffle'];
%                         taskrestshuffleall = [taskrestshuffleall taskrestshuffle'];
%                 
%                         break
%                 
%                     end
%                 end
%         
%                 if SplitHalf == 1
%             
%                     alltasktaskshuffle = alltaskfilestaskodd(tasktaskshuffle,:);
%                     alltaskrestshuffle = alltaskfilesresteven(taskrestshuffle,:);
%                     %alltasktaskshuffle = alltaskfilesresteven(randperm(size(alltaskfilestaskodd,1)),:);
%                     %alltaskrestshuffle = alltaskfilestaskodd(randperm(size(alltaskfilesresteven,1)),:);
%             
%                 else
%         
%                     alltasktaskshuffle = alltaskfilestask(tasktaskshuffle,:);
%                     alltaskrestshuffle = alltaskfilesrest(taskrestshuffle,:);
%                     %alltasktaskshuffle = alltaskfilestask(randperm(size(alltaskfilestask,1)),:);
%                     %alltaskrestshuffle = alltaskfilesrest(randperm(size(alltaskfilesrest,1)),:);
%             
%                 end
%         
%                 DiceCorrssimtemp = zeros(1,size(alltasktaskshuffle,1));
%                 NumOverlapsimtemp = zeros(1,size(alltasktaskshuffle,1));
%                 COMOverlapssimtemp = zeros(1,size(alltasktaskshuffle,1));
%     
%         
%                 for d = 1:size(alltasktaskshuffle,1)
%             
%                     vars_rest = unique(alltaskrestshuffle(d,:));
%                     vars_rest(1) = [];
%                     vars_task = unique(alltasktaskshuffle(d,:));
%                     vars_task(1) = [];
%     
%                     if DiceCorr == 1
%         
%                         dcorrdata = [];
%         
%                         for q = 1:length(alltaskrestshuffle(d,:))
%             
%                             if alltaskrestshuffle(d,q) > 0 && alltasktaskshuffle(d,q) > 0
%                 
%                                 dcorrdata = [dcorrdata;1 1];
%                 
%                             elseif alltaskrestshuffle(d,q) > 0
%                 
%                                 dcorrdata = [dcorrdata;1 0];
%                 
%                             elseif alltasktaskshuffle(d,q) > 0
%                 
%                                 dcorrdata = [dcorrdata;0 1];
%                 
%                             end
%             
%                         end
%                     
%                         if isempty(dcorrdata)
%                         
%                             dc = 0;
%                         
%                         else
%         
%                             dc = dice_coefficient_mod(dcorrdata(:,1),dcorrdata(:,2));
%                         
%                         end
%         
%                         DiceCorrssimtemp(d) = dc;
%     
%                     end
% 
%                     if AnyOverlap == 1
%    
%                         overlapvars = [];
%         
%                         for r = 1:length(vars_rest)
%             
%                             for s = 1:length(alltaskrestshuffle(d,:))
%             
%                                 if alltaskrestshuffle(d,s) == vars_rest(r) && alltasktaskshuffle(d,s) > 0
%                 
%                                     overlapvars = [overlapvars;vars_rest(r) alltasktaskshuffle(d,s)];
%                 
%                                 end
%                 
%                             end
%             
%                         end
%                     
%                         if isempty(overlapvars)
%                         
%                             NumOverlapsimtemp(d) = 0;
%                         
%                         else
%         
%                             NumOverlapsimtemp(d) = (length(unique(overlapvars(:,1))) + length(unique(overlapvars(:,2))))/(length(vars_rest) + length(vars_task));
%                         
%                         end
%     
%                     end
% 
%                     if COMOverlap == 1
%     
%                         %meanvals = [];
%                         COMOverlapsresttemp = [];
%                         COMOverlapstasktemp = [];
%         
%         
%                         taskverts = find(alltasktaskshuffle(d,:) > 0);
%                         restverts = find(alltaskrestshuffle(d,:) > 0);
%         
%                         for a = 1:2   %% loop over task and rest
%             
%                             if a == 1
%                 
%                                 nvars = length(vars_rest);
%                 
%                             else
%                 
%                                 nvars = length(vars_task);
%                 
%                             end
%         
%                             for s = 1:nvars
%                 
%                                 if a == 1
%             
%                                     currentvariant = find(alltaskrestshuffle(d,:) == vars_rest(s));
%                 
%                                 else
%                     
%                                     currentvariant = find(alltasktaskshuffle(d,:) == vars_task(s));
%                     
%                                 end
%                 
%                                 COMtemp = cifti_coords.data(currentvariant,:);
%                 
%                                 if length(currentvariant) > 1
%                     
%                                     centroid = [];
%                     
%                                     for t = 1:length(currentvariant)
%                         
%                                         jacknifevals = 1:length(currentvariant);        %% Leaves current point out of calculation
%                                         jacknifevals(t) = [];
%                     
%                                         % Euclidean distance
%                                         dist = mean(sqrt(sum((COMtemp(jacknifevals,:) - COMtemp(t,:)).^2, 2)));
%                         
%                                         %meanvals = [meanvals; dist];
%                         
%                                         if isempty(centroid) || dist < centroid   %% If new distance less than current least distance, use as center
%                             
%                                             centroid = currentvariant(t);
%                             
%                                         end
%                         
%                                     end
%                     
%                                     centervariant = centroid;
%                     
%                                 else
%                     
%                                     centervariant = currentvariant;
%                     
%                                 end
%                 
%                                 if a == 1
%                 
%                                     COMOverlapsresttemp = [COMOverlapsresttemp; centervariant];
%                     
%                                 else
%                     
%                                     COMOverlapstasktemp = [COMOverlapstasktemp; centervariant];
%                     
%                                 end 
%                             end
%                         end
%         
%                         Taskoverlaps = length(unique(intersect(COMOverlapstasktemp, restverts)));
%                         Restoverlaps = length(unique(intersect(COMOverlapsresttemp, taskverts)));
%         
%                         COMOverlapssimtemp(d) = (Taskoverlaps + Restoverlaps)/(length(vars_rest) + length(vars_task));
%             
%             
%                     end
%                 end
%     
%                 DiceCorrssim(x) = mean(DiceCorrssimtemp);
%             
%                 if AnyOverlap == 1
%                 
%                     NumOverlapsim(x) = mean(NumOverlapsimtemp);
%             
%                 end
%             
%                 if COMOverlap == 1
%             
%                     COMOverlapssim(x) = mean(COMOverlapssimtemp);
%                 
%                 end
%         
%                 DiceCorrssimAll(:,x) = DiceCorrssimtemp';
%         
%         
%                 if mod(x,100) == 0
%             
%                     disp([num2str(x) ' iterations completed'])
%             
%                 end
%     
%             end
%     
%             DiceCorrspval = prctile(DiceCorrssim, 95);
%         
%             if AnyOverlap == 1
%             
%                 NumOverlappval = prctile(NumOverlapsim, 95);
%         
%             end
%         
%             if COMOverlap == 1
%             
%                 COMOverlapspval = prctile(COMOverlapssim, 95);
%         
%             end
%     
%             DiceCorrssub = zeros(nsubs,1);      %% Calculate averaged values per subject for (pseudo) effect size
%     
%             for m = 1:size(DiceCorrssimAll,2)
%         
%                 for n = 1:size(DiceCorrssimAll,1)
%             
%                     for o = 1:nsubs
%             
%                         if tasktaskshuffleall(n,m) == o || taskrestshuffleall(n,m) == o
%                 
%                             DiceCorrssub(o) = DiceCorrssub(o) + DiceCorrssimAll(n,m);
%                 
%                         end
%                     end
%                 end
%             end
%     
%             DiceCorrssub = DiceCorrssub./(size(DiceCorrssimAll,2)*2);
%         
% 
%         end
%     end


%     if SubjectPlots == 1
%     
%         SubjectDiceCorrs = [];
%         
%         if AnyOverlap == 1
%         
%             SubjectNumOverlap = [];
%             
%         end
%         
%         if COMOverlap == 1
%             
%             SubjectCOMOverlap = [];
%             
%         end
%     
%         for j = 1:nfiles
%         
%             if SplitHalf == 1
%             
%                 alltasktasksubject = alltaskfilestaskodd(j,:);
%                 alltaskrestsubject = alltaskfilesresteven(j,:);
%             
%                 alltasktaskcomp = alltaskfilestaskodd;
%                 alltasktaskcomp(j,:) = [];
%                 alltaskrestcomp = alltaskfilesresteven;
%                 alltaskrestcomp(j,:) = [];
%             
%             else
%         
%                 alltasktasksubject = alltaskfilestask(j,:);
%                 alltaskrestsubject = alltaskfilesrest(j,:);
%             
%                 alltasktaskcomp = alltaskfilestask;
%                 alltasktaskcomp(j,:) = [];
%                 alltaskrestcomp = alltaskfilesrest;
%                 alltaskrestcomp(j,:) = [];
%             
%             end
%         
%             vars_rest = unique(alltaskrestsubject);
%             vars_rest(1) = [];
%             vars_task = unique(alltasktasksubject);
%             vars_task(1) = [];
% 
%             for d = 1:size(alltasktaskcomp,1)
%             
%                 vars_rest_comp = unique(alltaskrestcomp(d,:));
%                 vars_rest_comp(1) = [];
%                 vars_task_comp = unique(alltasktaskcomp(d,:));
%                 vars_task_comp(1) = [];
%             
%                 SubjectDiceCorrstemp = [];
%                 SubjectNumOverlaptemp = [];
%                 SubjectCOMOverlaptemp = [];
%             
%                 if DiceCorr == 1
%         
%                     dcorrdatataskcomp = [];
%                     dcorrdatarestcomp = [];
%         
%                     for q = 1:length(alltasktasksubject)
%             
%                         if alltaskrestsubject(:,q) > 0 && alltasktaskcomp(d,q) > 0
%                 
%                             dcorrdatataskcomp = [dcorrdatataskcomp;1 1];
%                 
%                         elseif alltaskrestsubject(:,q) > 0
%                 
%                             dcorrdatataskcomp = [dcorrdatataskcomp;1 0];
%                 
%                         elseif alltasktaskcomp(d,q) > 0
%                 
%                             dcorrdatataskcomp = [dcorrdatataskcomp;0 1];
%                 
%                         end
%             
%                     end
%                 
%                     for r = 1:length(alltaskrestsubject)
%             
%                         if alltasktasksubject(:,r) > 0 && alltaskrestcomp(d,r) > 0
%                 
%                             dcorrdatarestcomp = [dcorrdatarestcomp;1 1];
%                 
%                         elseif alltasktasksubject(:,r) > 0
%                 
%                             dcorrdatarestcomp = [dcorrdatarestcomp;1 0];
%                 
%                         elseif alltaskrestcomp(d,r) > 0
%                 
%                             dcorrdatarestcomp = [dcorrdatarestcomp;0 1];
%                 
%                         end
%             
%                     end
%                     
%                     if isempty(dcorrdatataskcomp) && isempty(dcorrdatarestcomp)
%                         
%                         dctaskcomp = 0;
%                         dcrestcomp = 0;
%                         
%                     elseif isempty(dcorrdatataskcomp)
%                         
%                         dctaskcomp = 0;
%                         dcrestcomp = dice_coefficient_mod(dcorrdatarestcomp(:,1),dcorrdatarestcomp(:,2));
%                         
%                     elseif isempty(dcorrdatarestcomp)
%                         
%                         dcrestcomp = 0;
%                         dctaskcomp = dice_coefficient_mod(dcorrdatataskcomp(:,1),dcorrdatataskcomp(:,2));
%                         
%                     else
%                         
%                         dctaskcomp = dice_coefficient_mod(dcorrdatataskcomp(:,1),dcorrdatataskcomp(:,2));
%                         dcrestcomp = dice_coefficient_mod(dcorrdatarestcomp(:,1),dcorrdatarestcomp(:,2));
%                         
%                     end
%         
%                     SubjectDiceCorrstemp = [SubjectDiceCorrstemp dctaskcomp dcrestcomp];
%     
%                 end
%             
%             
%                 if AnyOverlap == 1
%    
%                     overlapvarstaskcomp = [];
%                     overlapvarsrestcomp = [];
%         
%                     for r = 1:length(vars_rest)
%             
%                         for s = 1:length(alltaskrestsubject)
%             
%                             if alltaskrestsubject(:,s) == vars_rest(r) && alltasktaskcomp(d,s) > 0
%                 
%                                 overlapvarstaskcomp = [overlapvarstaskcomp;vars_rest(r) alltasktaskcomp(d,s)];
%                 
%                             end
%                 
%                         end
%             
%                     end
%                 
%                     for p = 1:length(vars_task)
%             
%                         for q = 1:length(alltasktasksubject)
%             
%                             if alltasktasksubject(:,q) == vars_task(p) && alltaskrestcomp(d,q) > 0
%                 
%                                 overlapvarsrestcomp = [overlapvarsrestcomp;vars_task(p) alltaskrestcomp(d,q)];
%                 
%                             end
%                 
%                         end
%             
%                     end
%                     
%                     if isempty(overlapvarsrestcomp) && isempty(overlapvarstaskcomp)
%                         
%                         SubjectNumOverlaptemp = [SubjectNumOverlaptemp 0 0];
%                         
%                     elseif isempty(overlapvarsrestcomp)
%                         
%                         SubjectNumOverlaptemp = [SubjectNumOverlaptemp (length(unique(overlapvarstaskcomp(:,1))) + length(unique(overlapvarstaskcomp(:,2))))/(length(vars_rest) + length(vars_task_comp)) 0];
%                         
%                     elseif isempty(overlapvarstaskcomp)
%                         
%                         SubjectNumOverlaptemp = [SubjectNumOverlaptemp 0 (length(unique(overlapvarsrestcomp(:,1))) + length(unique(overlapvarsrestcomp(:,2))))/(length(vars_task) + length(vars_rest_comp))];
%                         
%                     else
%         
%                         SubjectNumOverlaptemp = [SubjectNumOverlaptemp (length(unique(overlapvarstaskcomp(:,1))) + length(unique(overlapvarstaskcomp(:,2))))/(length(vars_rest) + length(vars_task_comp)) (length(unique(overlapvarsrestcomp(:,1))) + length(unique(overlapvarsrestcomp(:,2))))/(length(vars_task) + length(vars_rest_comp))];
%                         
%                     end
%     
%                 end
%             
%                 if COMOverlap == 1
%                 
%                     for z = 1:2     %% Loops over subject task - all rest and subject rest - all task comparisons
%                     
%                         %meanvals = [];
%                         COMOverlapsresttemp = [];
%                         COMOverlapstasktemp = [];
%         
%                         if z == 1     %% subject task - all rest comparison
%                         
%                             taskverts = find(alltasktasksubject > 0);
%                             restverts = find(alltaskrestcomp(d,:) > 0);
%                         
%                         elseif z == 2	  %% subject rest - all task comparison
%                         
%                             taskverts = find(alltasktaskcomp(d,:) > 0);
%                             restverts = find(alltaskrestsubject > 0);
%                         
%                         end
%         
%                         for a = 1:2   %% loop over task and rest
%             
%                             if a == 1    %% Rest loop
%                             
%                                 if z == 1     %% subject task - all rest comparison
%                 
%                                     nvars = length(vars_rest_comp);
%                                 
%                                 elseif z == 2	  %% subject rest - all task comparison
%                                 
%                                     nvars = length(vars_rest);
%                                 
%                                 end
%                 
%                             elseif a == 2   %% Task loop
%                             
%                                 if z == 1     %% subject task - all rest comparison
%                 
%                                     nvars = length(vars_task);
%                                 
%                                 elseif z == 2	  %% subject rest - all task comparison
%                                                                 
%                                     nvars = length(vars_task_comp);
%                                 
%                                 end
%                 
%                             end
%         
%                             for s = 1:nvars
%                 
%                                 if a == 1   %% Rest loop
%                                 
%                                     if z == 1     %% subject task - all rest comparison
%             
%                                         currentvariant = find(alltaskrestcomp(d,:) == vars_rest_comp(s));
%                                     
%                                     elseif z == 2	  %% subject rest - all task comparison
%                                     
%                                         currentvariant = find(alltaskrestsubject == vars_rest(s));
%                                     
%                                     end
%                 
%                                 elseif a == 2   %% Task loop
%                                 
%                                     if z == 1     %% subject task - all rest comparison
%                     
%                                         currentvariant = find(alltasktasksubject == vars_task(s));
%                                     
%                                     elseif z == 2	  %% subject rest - all task comparison
%                                     
%                                         currentvariant = find(alltasktaskcomp(d,:) == vars_task_comp(s));
%                                     
%                                     end
%                     
%                                 end
%                 
%                                 COMtemp = cifti_coords.data(currentvariant,:);
%                 
%                                 if length(currentvariant) > 1
%                     
%                                     centroid = [];
%                     
%                                     for t = 1:length(currentvariant)
%                         
%                                         jacknifevals = 1:length(currentvariant);        %% Leaves current point out of calculation
%                                         jacknifevals(t) = [];
%                     
%                                         % Euclidean distance
%                                         dist = mean(sqrt(sum((COMtemp(jacknifevals,:) - COMtemp(t,:)).^2, 2)));
%                         
%                                         %meanvals = [meanvals; dist];
%                         
%                                         if isempty(centroid) || dist < centroid   %% If new distance less than current least distance, use as center
%                             
%                                             centroid = currentvariant(t);
%                             
%                                         end
%                         
%                                     end
%                     
%                                     centervariant = centroid;
%                     
%                                 else
%                     
%                                     centervariant = currentvariant;
%                     
%                                 end
%                 
%                                 if a == 1
%                 
%                                     COMOverlapsresttemp = [COMOverlapsresttemp; centervariant];
%                     
%                                 else
%                     
%                                     COMOverlapstasktemp = [COMOverlapstasktemp; centervariant];
%                     
%                                 end 
%                             end
%                         end
%         
%                         Taskoverlaps = length(unique(intersect(COMOverlapstasktemp, restverts)));
%                         Restoverlaps = length(unique(intersect(COMOverlapsresttemp, taskverts)));
%                     
%                         if z == 1     %% subject task - all rest comparison
%                         
%                             restdivisor = length(vars_rest_comp);
%                             taskdivisor = length(vars_task);
%                         
%                         elseif z == 2	  %% subject rest - all task comparison
%                         
%                             restdivisor = length(vars_rest);
%                             taskdivisor = length(vars_task_comp);
%                         
%                         end   
%         
%                         SubjectCOMOverlaptemp = [SubjectCOMOverlaptemp (Taskoverlaps + Restoverlaps)/(restdivisor + taskdivisor)];
%             
%                     end
%             
%                 end
%             
%             end
%         
%             SubjectDiceCorrs = [SubjectDiceCorrs; SubjectDiceCorrstemp];
%             
%             if AnyOverlap == 1
%             
%                 SubjectNumOverlap = [SubjectNumOverlap; SubjectNumOverlaptemp];
%                 
%             end
%             
%             if COMOverlap == 1
%                 
%                 SubjectCOMOverlap = [SubjectCOMOverlap; SubjectCOMOverlaptemp];
%                 
%             end
%             
%         end
%         
%     end


%% plot results
%haven't started working on this part
    if plotresults == 1
    
       % if SubjectPlots == 1
        

%        		if SplitHalf == 1
%                 
%                 if AnyOverlap == 1 && COMOverlap == 1
%         
%                     plotdatawithinstate = [mean([DiceCorrsRestRest DiceCorrsTaskTask],2) mean([NumOverlapRestRest NumOverlapTaskTask],2) mean([COMOverlapsRestRest COMOverlapsTaskTask],2)];       
%                     plotdatabetweenstate = [DiceCorrsTaskRest NumOverlapTaskRest COMOverlapsTaskRest];
%                     plotdatanull = [mean(SubjectDiceCorrs,2) mean(SubjectNumOverlap,2) mean(SubjectCOMOverlap,2)];
%                     
%                 elseif AnyOverlap == 1
%                     
%                     plotdatawithinstate = [mean([DiceCorrsRestRest DiceCorrsTaskTask],2) mean([NumOverlapRestRest NumOverlapTaskTask],2)];       
%                     plotdatabetweenstate = [DiceCorrsTaskRest NumOverlapTaskRest];
%                     plotdatanull = [mean(SubjectDiceCorrs,2) mean(SubjectNumOverlap,2)];
%                     
%                 elseif COMOverlap == 1
%                     
%                     plotdatawithinstate = [mean([DiceCorrsRestRest DiceCorrsTaskTask],2) mean([COMOverlapsRestRest COMOverlapsTaskTask],2)];       
%                     plotdatabetweenstate = [DiceCorrsTaskRest COMOverlapsTaskRest];
%                     plotdatanull = [mean(SubjectDiceCorrs,2) mean(SubjectCOMOverlap,2)];
%                     
%                 else
%                     
%                    	plotdatawithinstate = [mean([DiceCorrsRestRest DiceCorrsTaskTask],2)];       
%                     plotdatabetweenstate = DiceCorrsTaskRest;
%                     plotdatanull = mean(SubjectDiceCorrs,2);
%                     
%                 end
%             
%             else
                
%                 if AnyOverlap == 1 && COMOverlap == 1
%             
%                     plotdatabetweenstate = [DiceCorrs NumOverlap COMOverlaps];
%                     plotdatanull = [mean(SubjectDiceCorrs,2) mean(SubjectNumOverlap,2) mean(SubjectCOMOverlap,2)];
%                     
%                 elseif AnyOverlap == 1
%                     
%                     plotdatabetweenstate = [DiceCorrs NumOverlap];
%                     plotdatanull = [mean(SubjectDiceCorrs,2) mean(SubjectNumOverlap,2)];
%                     
%                 elseif COMOverlap == 1
%                     
%                     plotdatabetweenstate = [DiceCorrs COMOverlaps];
%                     plotdatanull = [mean(SubjectDiceCorrs,2) mean(SubjectCOMOverlap,2)];
%                     
%                 else
                    
                    plotdatabetweenstate = DiceCorrs;
                    plotdatanull = mean(SubjectDiceCorrs,2);
                    
%                 end
            
           %end
        
            for z = 1:size(plotdatanull,1)
            
%                 if SplitHalf == 1
%                     
%                     if AnyOverlap == 1 && COMOverlap == 1
%             
%                         SEplotDicewithin(z) = std([DiceCorrsRestRest(z,:); DiceCorrsTaskTask(z,:)])/sqrt(size([DiceCorrsRestRest(z,:); DiceCorrsTaskTask(z,:)],2));
%                         SEplotAnywithin(z) = std([NumOverlapRestRest(z,:); NumOverlapTaskTask(z,:)])/sqrt(size([NumOverlapRestRest(z,:); NumOverlapTaskTask(z,:)],2));
%                         SEplotCOMwithin(z) = std([COMOverlapsRestRest(z,:); COMOverlapsTaskTask(z,:)])/sqrt(size([COMOverlapsRestRest(z,:); COMOverlapsTaskTask(z,:)],2));
%                         
%                     elseif AnyOverlap == 1
%                         
%                      	SEplotDicewithin(z) = std([DiceCorrsRestRest(z,:); DiceCorrsTaskTask(z,:)])/sqrt(size([DiceCorrsRestRest(z,:); DiceCorrsTaskTask(z,:)],2));
%                         SEplotAnywithin(z) = std([NumOverlapRestRest(z,:); NumOverlapTaskTask(z,:)])/sqrt(size([NumOverlapRestRest(z,:); NumOverlapTaskTask(z,:)],2));
%                         
%                     elseif COMOverlap == 1
%                         
%                         SEplotDicewithin(z) = std([DiceCorrsRestRest(z,:); DiceCorrsTaskTask(z,:)])/sqrt(size([DiceCorrsRestRest(z,:); DiceCorrsTaskTask(z,:)],2));
%                         SEplotCOMwithin(z) = std([COMOverlapsRestRest(z,:); COMOverlapsTaskTask(z,:)])/sqrt(size([COMOverlapsRestRest(z,:); COMOverlapsTaskTask(z,:)],2));
%                         
%                     else
                        
                        SEplotDicewithin(z) = std([DiceCorrsRestRest(z,:); DiceCorrsTaskTask(z,:)])/sqrt(size([DiceCorrsRestRest(z,:); DiceCorrsTaskTask(z,:)],2));
                        
%                     end
%             
%                 end
                
%                 if AnyOverlap == 1 && COMOverlap == 1
%                 
%                     SEplotDicenull(z) = std(SubjectDiceCorrs(z,:))/sqrt(size(SubjectDiceCorrs(z,:),2));
%                     SEplotAnynull(z) = std(SubjectNumOverlap(z,:))/sqrt(size(SubjectNumOverlap(z,:),2));
%                     SEplotCOMnull(z) = std(SubjectCOMOverlap(z,:))/sqrt(size(SubjectCOMOverlap(z,:),2));
%                     
%                 elseif AnyOverlap == 1
%             
%                     SEplotDicenull(z) = std(SubjectDiceCorrs(z,:))/sqrt(size(SubjectDiceCorrs(z,:),2));
%                     SEplotAnynull(z) = std(SubjectNumOverlap(z,:))/sqrt(size(SubjectNumOverlap(z,:),2));
%                     
%                 elseif COMOverlap == 1
%                     
%                     SEplotDicenull(z) = std(SubjectDiceCorrs(z,:))/sqrt(size(SubjectDiceCorrs(z,:),2));
%                     SEplotCOMnull(z) = std(SubjectCOMOverlap(z,:))/sqrt(size(SubjectCOMOverlap(z,:),2));
%                     
%                 else
                    
                    SEplotDicenull(z) = std(SubjectDiceCorrs(z,:))/sqrt(size(SubjectDiceCorrs(z,:),2));
%                     
%                 end
%             end
            
            for f = 1:3     %% Loops over Dice, Any overlap, and COM overlap
            
                if f == 1       %% Plot Dice
                
                    SEPlotwithin = SEplotDicewithin;
                    SEPlotnull = SEplotDicenull;
                
                elseif f == 2 && AnyOverlap == 1       %% Plot Any overlap
                
                    SEPlotwithin = SEplotAnywithin;
                    SEPlotnull = SEplotAnynull;
                
                elseif f == 3 && COMOverlap == 1      %% Plot COM overlap
                
                    SEPlotwithin = SEplotCOMwithin;
                    SEPlotnull = SEplotCOMnull;
                
                end
        
        
                %%% Plot Dice correlations by subject
                
                if f == 1 || (f == 2 && AnyOverlap == 1) || (f == 3 && COMOverlap == 1)
        
                    if SplitHalf == 1
            
                        bar(categorical({'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10'}), plotdatawithinstate(:,f), 'FaceColor', 'w', 'EdgeColor', 'g');
                        hold on
                        errorbar(categorical({'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10'}),plotdatawithinstate(:,f),SEPlotwithin,'.', 'Color', 'b');
                        hold on
            
                    end
        
                    bar(categorical({'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10'}), plotdatabetweenstate(:,f), 'FaceColor', 'w', 'FaceAlpha', .5);
                    hold on
        
                    bar(categorical({'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10'}), plotdatanull(:,f), 'FaceColor', 'k');
                    hold on
                    errorbar(categorical({'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10'}),plotdatanull(:,f),SEPlotnull,'.', 'Color', 'r');
                    hold on
        
                    set(gca, 'FontSize',16)
                    c = get(gca, 'Children');
                    hleg1 = legend(c([3 5 2]), 'Between State Comparison', 'Within State Comparison', 'Between Subject Comparison', 'Location', 'NorthEast');
                    hleg1.FontSize = 14;
                    ylim([0 1]);
            
                    if f == 1       %% Plot Dice
                
                        title(['All Subjects Dice Correlations Matched Data'], 'fontsize',18);
                        fstring = 'DiceCorrelations';
                
                    elseif f == 2    %% Plot Any overlap
                
                        title(['All Subjects Any Overlap Matched Data'], 'fontsize',18);
                        fstring = 'AnyOverlap';
                
                    elseif f == 3    %% Plot COM overlap
                
                        title(['All Subjects COM Overlap Matched Data'], 'fontsize',18);
                        fstring = 'COMOverlap';
                
                    end
            
            
                    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
                
                    if MatchedMaps == 1 && SNRExclude == 1 && SizeExclude == 1 && SplitHalf == 1 && numel(thresholds) > 1 && AbsoluteThresholds == 1
                    
                        filename = ['/BySubjectPlots_' fstring '_MatchedData_SplitHalf_Overlap_SNRAndSizeExclude_' num2str(thresholds(v)) '_Abs.jpg'];
                
                    elseif MatchedMaps == 1 && SNRExclude == 1 && SizeExclude == 1 && SplitHalf == 1 && numel(thresholds) > 1
        
                        filename = ['/BySubjectPlots_' fstring '_MatchedData_SplitHalf_Overlap_SNRAndSizeExclude_' num2str(thresholds(v)) '_Percent.jpg'];
        
                    elseif MatchedMaps == 1 && SNRExclude == 1 && SizeExclude == 1 && SplitHalf == 1
        
                        filename = ['/BySubjectPlots_' fstring '_MatchedData_SplitHalf_Overlap_SNRAndSizeExclude.jpg'];
                        
                    elseif MatchedMaps == 1 && SNRExclude == 1 && SizeExclude == 1
        
                        filename = ['/BySubjectPlots_' fstring '_MatchedData_Overlap_SNRAndSizeExclude.jpg'];
                    
                    elseif MatchedMaps == 1 && SNRExclude == 1 && SplitHalf == 1 && numel(thresholds) > 1 && AbsoluteThresholds == 1
                    
                        filename = ['/BySubjectPlots_' fstring '_MatchedData_SplitHalf_Overlap_SNRExclude_' num2str(thresholds(v)) '_Abs.jpg'];
                    
                    elseif MatchedMaps == 1 && SNRExclude == 1 && SplitHalf == 1 && numel(thresholds) > 1
        
                        filename = ['/BySubjectPlots_' fstring '_MatchedData_SplitHalf_Overlap_SNRExclude_' num2str(thresholds(v)) '_Percent.jpg'];
        
                    elseif MatchedMaps == 1 && SNRExclude == 1 && SplitHalf == 1
        
                        filename = ['/BySubjectPlots_' fstring '_MatchedData_SplitHalf_Overlap_SNRExclude.jpg'];
                            
                    elseif MatchedMaps == 1 && SNRExclude == 1
        
                        filename = ['/BySubjectPlots_' fstring '_MatchedData_Overlap_SNRExclude.jpg'];
        
                    elseif MatchedMaps == 1 && SizeExclude == 1 && SplitHalf == 1
        
                        filename = ['/BySubjectPlots_' fstring '_MatchedData_SplitHalf_Overlap_SizeExclude.jpg'];
                            
                    elseif MatchedMaps == 1 && SizeExclude == 1
        
                        filename = ['/BySubjectPlots_' fstring '_MatchedData_Overlap_SizeExclude.jpg'];
                
                    elseif MatchedMaps == 1 && SplitHalf == 1
        
                        filename = ['/BySubjectPlots_' fstring '_MatchedData_SplitHalf_Overlap_NoExclude.jpg'];    
        
                    else
        
                        filename = ['/BySubjectPlots_' fstring '_MatchedData_Overlap_NoExclude.jpg'];
                            
                    end
                    
                    saveas(gcf,[outputdir filename])
    
                    close gcf
        
                end
            end
        end
    
        if BarGraph == 1
    
            if randomizevals == 1
                
                if AnyOverlap == 1 && COMOverlap == 1
        
                    simplotdata = [DiceCorrspval NumOverlappval COMOverlapspval];
                    
                elseif AnyOverlap == 1
                    
                    simplotdata = [DiceCorrspval NumOverlappval];
                    
                elseif COMOverlap == 1
                    
                    simplotdata = [DiceCorrspval COMOverlapspval];
                    
                else
                    
                    simplotdata = DiceCorrspval;
                    
                end

            end
    
            if SplitHalf == 1
                
                if AnyOverlap == 1 && COMOverlap == 1
        
                    plotdata = [DiceCorrsTaskRest NumOverlapTaskRest COMOverlapsTaskRest];
                    plotdata2 = [[DiceCorrsTaskTask;DiceCorrsRestRest] [NumOverlapTaskTask;NumOverlapRestRest] [COMOverlapsTaskTask;COMOverlapsRestRest]];
                    
                elseif AnyOverlap == 1
                    
                    plotdata = [DiceCorrsTaskRest NumOverlapTaskRest];
                    plotdata2 = [[DiceCorrsTaskTask;DiceCorrsRestRest] [NumOverlapTaskTask;NumOverlapRestRest]];
                    
                elseif COMOverlap == 1
                    
                    plotdata = [DiceCorrsTaskRest COMOverlapsTaskRest];
                    plotdata2 = [[DiceCorrsTaskTask;DiceCorrsRestRest] [COMOverlapsTaskTask;COMOverlapsRestRest]];
                    
                else
                    
                    plotdata = DiceCorrsTaskRest;
                    plotdata2 = [DiceCorrsTaskTask;DiceCorrsRestRest];
                    
                end
        
            else
                
                if AnyOverlap == 1 && COMOverlap == 1
    
                    plotdata = [DiceCorrs NumOverlap COMOverlaps];
                    
                elseif AnyOverlap == 1
                    
                    plotdata = [DiceCorrs NumOverlap];
                    
                elseif COMOverlap == 1
                    
                    plotdata = [DiceCorrs COMOverlaps];
                    
                else
                    
                    plotdata = DiceCorrs;
                    
                end
    
            end
    
            for z = 1:size(plotdata,2)
        
                if SplitHalf == 1
            
                    SEplot(z) = std(plotdata(:,z))/sqrt(size(plotdata,1));
                    SEplot2(z) = std(plotdata2(:,z))/sqrt(size(plotdata2,1));
            
                else

                    SEplot(z) = std(plotdata(:,z))/sqrt(size(plotdata,1));
                        
                end
            end
            
            if AnyOverlap == 1 && COMOverlap == 1
                
                xlabels = {'Dice Correlations', 'Any Overlap', 'COM Overlap'};
                
            elseif AnyOverlap == 1
                
                xlabels = {'Dice Correlations', 'Any Overlap'};
                
            elseif COMOverlap == 1
                
                xlabels = {'Dice Correlations', 'COM Overlap'};
                
            else
                
                xlabels = {'Dice Correlations'};
                
            end
    
            if SplitHalf == 1
                
                bar(categorical(xlabels), mean(plotdata2,1), 'FaceColor', 'w', 'EdgeColor', 'g');
               	hold on
               	errorbar(categorical(xlabels),mean(plotdata2,1),SEplot2,'.');
             	hold on
                    
            end

          	bar(categorical(xlabels), mean(plotdata,1), 'FaceColor', 'w');
         	hold on
          	errorbar(categorical(xlabels),mean(plotdata,1),SEplot,'.');
          	hold on

            if randomizevals == 1

             	bar(categorical(xlabels), simplotdata, 'FaceColor', 'k');
               	hold on
                    
            end
      
            plot(categorical(xlabels), plotdata(9,:), '.', 'MarkerSize', 15, 'color', [1, 0.5, 0]);
            hold on
            plot(categorical(xlabels), plotdata(8,:), '.', 'MarkerSize', 15, 'color', [0, 0.6, 0.6]);
            hold on    
            plot(categorical(xlabels), plotdata(7,:), '.', 'MarkerSize', 15, 'color', [1, 0, 1]);
            hold on
            plot(categorical(xlabels), plotdata(6,:), '.', 'MarkerSize', 15, 'color', [0.2, 1, 1]);
            hold on
            plot(categorical(xlabels), plotdata(5,:), '.', 'MarkerSize', 15, 'color', [0, 0, 1]);
            hold on
            plot(categorical(xlabels), plotdata(4,:), '.', 'MarkerSize', 15, 'color', [1, 0, 0]);
            hold on
            plot(categorical(xlabels), plotdata(3,:), '.', 'MarkerSize', 15, 'color', [0, 1, 0]);
            hold on
            plot(categorical(xlabels), plotdata(2,:), '.', 'MarkerSize', 15, 'color', [0.9, 0.9, 0]);
            hold on
            plot(categorical(xlabels), plotdata(1,:), '.', 'MarkerSize', 15, 'color', [0, 0, 0]);
            hold on
    
            ylim([0 1]);
                
            if FullMaps == 1
                
                title(['All Subjects Variant Overlap Unmatched Data'], 'fontsize',18)
                    
            elseif MatchedMaps == 1
                    
                title(['All Subjects Variant Overlap Matched Data'], 'fontsize',18)
                    
            end
                
            h = findobj(gca,'Type','line');

            legendpos = 'SouthEast';
                
            legend(h(1:9), 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', legendpos)
    
            if FullMaps == 1
        
                filename = '/AllSubjects_AllData_Overlap_SNRExclude.jpg';
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end
                
            elseif MatchedMaps == 1 && SNRExclude == 1 && SizeExclude == 1 && SplitHalf == 1 && numel(thresholds) > 1 && AbsoluteThresholds == 1
                
                filename = ['/AllSubjects_MatchedData_SplitHalf_Overlap_SNRAndSizeExclude_' num2str(thresholds(v)) '_Abs.jpg'];
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end
                
            elseif MatchedMaps == 1 && SNRExclude == 1 && SizeExclude == 1 && SplitHalf == 1 && numel(thresholds) > 1
        
                filename = ['/AllSubjects_MatchedData_SplitHalf_Overlap_SNRAndSizeExclude_' num2str(thresholds(v)) '_Percent.jpg'];
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end
                
            elseif MatchedMaps == 1 && SNRExclude == 1 && SizeExclude == 1 && SplitHalf == 1
        
                filename = '/AllSubjects_MatchedData_SplitHalf_Overlap_SNRAndSizeExclude.jpg';
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end
                            
                saveas(gcf,[outputdir filename])
                        
            elseif MatchedMaps == 1 && SNRExclude == 1 && SizeExclude == 1
        
                filename = '/AllSubjects_MatchedData_Overlap_SNRAndSizeExclude.jpg';
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end
                            
                saveas(gcf,[outputdir filename])
                
            elseif MatchedMaps == 1 && SNRExclude == 1 && SplitHalf == 1 && numel(thresholds) > 1 && AbsoluteThresholds == 1
        
                filename = ['/AllSubjects_MatchedData_SplitHalf_Overlap_SNRExclude_' num2str(thresholds(v)) '_Abs.jpg'];
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end
                
            elseif MatchedMaps == 1 && SNRExclude == 1 && SplitHalf == 1 && numel(thresholds) > 1
        
                filename = ['/AllSubjects_MatchedData_SplitHalf_Overlap_SNRExclude_' num2str(thresholds(v)) '_Percent.jpg'];
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end
        
            elseif MatchedMaps == 1 && SNRExclude == 1 && SplitHalf == 1
        
                filename = '/AllSubjects_MatchedData_SplitHalf_Overlap_SNRExclude.jpg';
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end
                            
            elseif MatchedMaps == 1 && SNRExclude == 1
        
                filename = '/AllSubjects_MatchedData_Overlap_SNRExclude.jpg';
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end
        
            elseif MatchedMaps == 1 && SizeExclude == 1 && SplitHalf == 1
        
                filename = '/AllSubjects_MatchedData_SplitHalf_Overlap_SizeExclude.jpg';
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end
                            
            elseif MatchedMaps == 1 && SizeExclude == 1
        
                filename = '/AllSubjects_MatchedData_Overlap_SizeExclude.jpg';
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end
        
            else
        
                filename = '/AllSubjects_MatchedData_Overlap_NoExclude.jpg';
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end
                            
            end
            
            if AnyOverlap == 1 && COMOverlap == 1
                
                saveas(gcf,[outputdir filename])
                
            else
                
                filename = strrep(filename, '.jpg', '_DiceOnly.jpg');

                saveas(gcf,[outputdir filename])
                
            end
    
            close gcf
        
        end
    
        if BoxPlot == 1
            
      	if SplitHalf == 1
                
         	if AnyOverlap == 1 && COMOverlap == 1
        
             	plotdata = [DiceCorrsTaskRest NumOverlapTaskRest COMOverlapsTaskRest];
             	plotdata2 = [[DiceCorrsTaskTask;DiceCorrsRestRest] [NumOverlapTaskTask;NumOverlapRestRest] [COMOverlapsTaskTask;COMOverlapsRestRest]];
                    
          	elseif AnyOverlap == 1
                    
              	plotdata = [DiceCorrsTaskRest NumOverlapTaskRest];
            	plotdata2 = [[DiceCorrsTaskTask;DiceCorrsRestRest] [NumOverlapTaskTask;NumOverlapRestRest]];
                    
            elseif COMOverlap == 1
                    
             	plotdata = [DiceCorrsTaskRest COMOverlapsTaskRest];
              	plotdata2 = [[DiceCorrsTaskTask;DiceCorrsRestRest] [COMOverlapsTaskTask;COMOverlapsRestRest]];
                    
            else
                    
             	plotdata = DiceCorrsTaskRest;
             	plotdata2 = [DiceCorrsTaskTask;DiceCorrsRestRest];
                    
            end
        
        else
                
         	if AnyOverlap == 1 && COMOverlap == 1
    
             	plotdata = [DiceCorrs NumOverlap COMOverlaps];
                    
        	elseif AnyOverlap == 1
                    
             	plotdata = [DiceCorrs NumOverlap];
                    
          	elseif COMOverlap == 1
                    
                plotdata = [DiceCorrs COMOverlaps];
                    
            else
                    
             	plotdata = DiceCorrs;
                    
            end
        end
        
            if randomizevals == 1
                
                if FinalFigure == 1
                    
                        x = [];
                        group = [];
                        positions = [1 2 3 4 5 6 7 8 9 10 10.25 10.5];
                        positionssubs = [1 2 3 4 5 6 7 8 9 10];
                        scatterposbetween = [sort(repmat([1:1:9],1,2))];
                        scatterposwithin = [sort(repmat([1:1:9],1,2))];
                        scatterposacross = [sort(repmat([1:1:9],1,16))];

                        for g = 1:nsubs+1
                            
                            if g < nsubs+1

                                %x = [x [DiceCorrsTaskRest(g,:) DiceCorrsTaskTask(g,:) DiceCorrsRestRest(g,:) DiceCorrSubs(:,1)']];
                                
                                x = [x 0];
                                
                            else
                                
                                x = [x [reshape(DiceCorrsTaskRest,1,size(DiceCorrsTaskRest,1)*size(DiceCorrsTaskRest,2)) DiceCorrsTaskTask' DiceCorrsRestRest' DiceCorrssim']];
                                
                            end
                        end

                        for i = 1:length(positionssubs)
                            
                            if i < length(positionssubs)

                                %group = [group [repmat(positionssubs(i),1,size(DiceCorrsTaskRest,2)) repmat(positionssubs(i)+.25,1,size(DiceCorrsTaskTask,2)+size(DiceCorrsRestRest,2)) repmat(positionssubs(i)+.5,1,size(DiceCorrSubs,1),1)]];
                                
                                group = [group positionssubs(i)];
                                
                            else
                            
                                group = [group [repmat(positionssubs(i),1,size(DiceCorrsTaskRest,1)*size(DiceCorrsTaskRest,2)) repmat(positionssubs(i)+.25,1,size(DiceCorrsTaskTask,1)+size(DiceCorrsRestRest,1)) repmat(positionssubs(i)+.5,1,size(DiceCorrssim,1),1)]];
                                
                            end
                        end

                        %group = [repmat(1,1,size(DiceCorrsTaskRest,2)+size(DiceCorrsTaskTask,2)+size(DiceCorrsRestRest,2)+size(DiceCorrSubs,1)) repmat(2,1,length(DiceCorrsTaskTask)+length(DiceCorrsRestRest)) repmat(3,1,length(DiceCorrssim))];
                        %positions = [1 2 3 4 5 6 7 8 9];
                        %scatterpos = [1 2];

                        boxplot(x,group, 'BoxStyle', 'outline', 'positions', positions, 'Symbol', '');
                        hold on
                        scatter(scatterposbetween, [DiceCorrsTaskRest(1,:) DiceCorrsTaskRest(2,:) DiceCorrsTaskRest(3,:) DiceCorrsTaskRest(4,:) DiceCorrsTaskRest(5,:) DiceCorrsTaskRest(6,:) DiceCorrsTaskRest(7,:) DiceCorrsTaskRest(8,:) DiceCorrsTaskRest(9,:)], 100, 'filled', 'MarkerFaceColor', 'b');
                        hold on
                        scatter(scatterposwithin, [DiceCorrsRestRest(1,:) DiceCorrsTaskTask(1,:) DiceCorrsRestRest(2,:) DiceCorrsTaskTask(2,:) DiceCorrsRestRest(3,:) DiceCorrsTaskTask(3,:) DiceCorrsRestRest(4,:) DiceCorrsTaskTask(4,:) DiceCorrsRestRest(5,:) DiceCorrsTaskTask(5,:) DiceCorrsRestRest(6,:) DiceCorrsTaskTask(6,:) DiceCorrsRestRest(7,:) DiceCorrsTaskTask(7,:) DiceCorrsRestRest(8,:) DiceCorrsTaskTask(8,:) DiceCorrsRestRest(9,:) DiceCorrsTaskTask(9,:)], 100, 'filled', 'MarkerFaceColor', 'g');
                        hold on
                        scatter(scatterposacross, [DiceCorrSubs(:,1)' DiceCorrSubs(:,2)' DiceCorrSubs(:,3)' DiceCorrSubs(:,4)' DiceCorrSubs(:,5)' DiceCorrSubs(:,6)' DiceCorrSubs(:,7)' DiceCorrSubs(:,8)' DiceCorrSubs(:,9)'], 100, 'jitter','on', 'jitterAmount', 0.1, 'MarkerEdgeColor', [0, 0, 0]);
                        hold on
                        set(gca,'xtick',[positions(1) positions(2) positions(3) positions(4) positions(5) positions(6) positions(7) positions(8) positions(9) mean(positions(10:12))])
                        set(gca,'xticklabel',{'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Group Average'}, 'FontSize',16)
                        ylabel('Dice Correlation', 'FontSize',16)

                        color = ['k', 'g', 'b'];

                        h = findobj(gca,'Tag','Box');

                        for j=1:3

                            patch(get(h(j),'XData'),get(h(j),'YData'),color(j), 'FaceColor', 'none', 'EdgeColor', color(j));

                        end

                        c = get(gca, 'Children');

                        hleg1 = legend(c(1:3), 'Between State Comparison', 'Within State Comparison', 'Between Subject Comparison', 'Location', 'NorthEast');
                        hleg1.FontSize = 16;
                        ylim([0 1]);
                        
                        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5, 0.7, 0.5, 0.7]);
                
                elseif AnyOverlap == 1 && COMOverlap == 1
            
                    x = [DiceCorrsTaskRest' DiceCorrsTaskTask' DiceCorrsRestRest' DiceCorrssim NumOverlapTaskRest' NumOverlapTaskTask' NumOverlapRestRest' NumOverlapsim COMOverlapsTaskRest' COMOverlapsTaskTask' COMOverlapsRestRest' COMOverlapssim];
                    group = [repmat(1,1,length(DiceCorrsTaskRest)) repmat(2,1,length(DiceCorrsTaskTask)+length(DiceCorrsRestRest)) repmat(3,1,length(DiceCorrssim)) repmat(4,1,length(NumOverlapTaskRest)) repmat(5,1,length(NumOverlapTaskTask)+length(NumOverlapRestRest)) repmat(6,1,length(NumOverlapsim)) repmat(7,1,length(COMOverlapsTaskRest)) repmat(8,1,length(COMOverlapsTaskTask)+length(COMOverlapsRestRest)) repmat(9,1,length(COMOverlapssim))];
                    positions = [1 1.25 1.5 2 2.25 2.5 3 3.25 3.5];
                    scatterpos = [1 2 3 1.25 2.25 3.25];
                    
                    boxplot(x,group, 'BoxStyle', 'outline', 'positions', positions, 'Symbol', '');
                    hold on
                    scatter(scatterpos, [plotdata(9,:) mean([plotdata2(9,1) plotdata2(18,1)]) mean([plotdata2(9,2) plotdata2(18,2)]) mean([plotdata2(9,3) plotdata2(18,3)])], 'filled', 'MarkerFaceColor', [1, 0.5, 0]);
                    hold on
                    scatter(scatterpos, [plotdata(8,:) mean([plotdata2(8,1) plotdata2(17,1)]) mean([plotdata2(8,2) plotdata2(17,2)]) mean([plotdata2(8,3) plotdata2(17,3)])], 'filled', 'MarkerFaceColor', [0, 0.6, 0.6]);
                    hold on    
                    scatter(scatterpos, [plotdata(7,:) mean([plotdata2(7,1) plotdata2(16,1)]) mean([plotdata2(7,2) plotdata2(16,2)]) mean([plotdata2(7,3) plotdata2(16,3)])], 'filled', 'MarkerFaceColor', [1, 0, 1]);
                    hold on
                    scatter(scatterpos, [plotdata(6,:) mean([plotdata2(6,1) plotdata2(15,1)]) mean([plotdata2(6,2) plotdata2(15,2)]) mean([plotdata2(6,3) plotdata2(15,3)])], 'filled', 'MarkerFaceColor', [0.2, 1, 1]);
                    hold on
                    scatter(scatterpos, [plotdata(5,:) mean([plotdata2(5,1) plotdata2(14,1)]) mean([plotdata2(5,2) plotdata2(14,2)]) mean([plotdata2(5,3) plotdata2(14,3)])], 'filled', 'MarkerFaceColor', [0, 0, 1]);
                    hold on
                    scatter(scatterpos, [plotdata(4,:) mean([plotdata2(4,1) plotdata2(13,1)]) mean([plotdata2(4,2) plotdata2(13,2)]) mean([plotdata2(4,3) plotdata2(13,3)])], 'filled', 'MarkerFaceColor', [1, 0, 0]);
                    hold on
                    scatter(scatterpos, [plotdata(3,:) mean([plotdata2(3,1) plotdata2(12,1)]) mean([plotdata2(3,2) plotdata2(12,2)]) mean([plotdata2(3,3) plotdata2(12,3)])], 'filled', 'MarkerFaceColor', [0, 1, 0]);
                    hold on
                    scatter(scatterpos, [plotdata(2,:) mean([plotdata2(2,1) plotdata2(11,1)]) mean([plotdata2(2,2) plotdata2(11,2)]) mean([plotdata2(2,3) plotdata2(11,3)])], 'filled', 'MarkerFaceColor', [0.9, 0.9, 0]);
                    hold on
                    scatter(scatterpos, [plotdata(1,:) mean([plotdata2(1,1) plotdata2(10,1)]) mean([plotdata2(1,2) plotdata2(10,2)]) mean([plotdata2(1,3) plotdata2(10,3)])], 'filled', 'MarkerFaceColor', [0, 0, 0]);
                    hold on
                    scatter([repmat(1.5,1,length(DiceCorrssim)) repmat(2.5,1,length(NumOverlapsim)) repmat(3.5,1,length(COMOverlapssim))], [DiceCorrssim NumOverlapsim COMOverlapssim], 'jitter','on', 'jitterAmount', 0.1, 'MarkerEdgeColor', [1, .25, 0]);
                    hold on

                    set(gca,'xtick',[mean(positions(1:3)) mean(positions(4:6)) mean(positions(7:9)) ])
                    set(gca,'xticklabel',{'Dice Correlations', 'Any Overlap', 'COM Overlap'}, 'FontSize',12)

                    color = ['k', 'g', 'b', 'k', 'g', 'b', 'k', 'g', 'b'];
                    %color = ['k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k'];
                    h = findobj(gca,'Tag','Box');
        
                    for j=1:length(h)
            
                        patch(get(h(j),'XData'),get(h(j),'YData'),color(j), 'FaceColor', 'none', 'EdgeColor', color(j));
            
                    end

                    c = get(gca, 'Children');
            
                    m = findobj(gca,'Type','scatter');
            
                    hleg1 = legend([m(2:10); c(1:3)], 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Between State Comparison', 'Within State Comparison', 'Between Subject Comparison', 'Location', 'NorthEast');
                    %hleg1 = legend(m(2:10), 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', 'NorthEast');
                    hleg1.FontSize = 14;
                    ylim([0 1]);
            
                    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
            
                    %hleg1 = legend(c(1:3), 'Between State Comparison', 'Within State Comparison', 'Between Subject Comparison', 'Location', 'NorthWest');

                    %legendpos = 'NorthEast';
                
                    %hleg2 = legend(m(2:10), 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', legendpos);
            
     
                else
                    
                    if ~(BoxPlotBySubject == 1)
                    
                        if permvals == 1
                        
                            DiceCorrssim = DiceCorrssim';
                        
                        end
                    
                        x = [DiceCorrsTaskRest' DiceCorrsTaskTask' DiceCorrsRestRest' DiceCorrssim];
                        group = [repmat(1,1,length(DiceCorrsTaskRest)) repmat(2,1,length(DiceCorrsTaskTask)+length(DiceCorrsRestRest)) repmat(3,1,length(DiceCorrssim))];
                        positions = [1 2 3];
                        scatterpos = [1 2];

                        boxplot(x,group, 'BoxStyle', 'outline', 'positions', positions, 'Symbol', '');
                        hold on
                        scatter(scatterpos, [plotdata(9,:) mean([plotdata2(9,1) plotdata2(18,1)])], 'filled', 'MarkerFaceColor', [1, 0.5, 0]);
                        hold on
                        scatter(scatterpos, [plotdata(8,:) mean([plotdata2(8,1) plotdata2(17,1)])], 'filled', 'MarkerFaceColor', [0, 0.6, 0.6]);
                        hold on    
                        scatter(scatterpos, [plotdata(7,:) mean([plotdata2(7,1) plotdata2(16,1)])], 'filled', 'MarkerFaceColor', [1, 0, 1]);
                        hold on
                        scatter(scatterpos, [plotdata(6,:) mean([plotdata2(6,1) plotdata2(15,1)])], 'filled', 'MarkerFaceColor', [0.2, 1, 1]);
                        hold on
                        scatter(scatterpos, [plotdata(5,:) mean([plotdata2(5,1) plotdata2(14,1)])], 'filled', 'MarkerFaceColor', [0, 0, 1]);
                        hold on
                        scatter(scatterpos, [plotdata(4,:) mean([plotdata2(4,1) plotdata2(13,1)])], 'filled', 'MarkerFaceColor', [1, 0, 0]);
                        hold on
                        scatter(scatterpos, [plotdata(3,:) mean([plotdata2(3,1) plotdata2(12,1)])], 'filled', 'MarkerFaceColor', [0, 1, 0]);
                        hold on
                        scatter(scatterpos, [plotdata(2,:) mean([plotdata2(2,1) plotdata2(11,1)])], 'filled', 'MarkerFaceColor', [0.9, 0.9, 0]);
                        hold on
                        scatter(scatterpos, [plotdata(1,:) mean([plotdata2(1,1) plotdata2(10,1)])], 'filled', 'MarkerFaceColor', [0, 0, 0]);
                        hold on
                        scatter(repmat(3,1,length(DiceCorrssim)), DiceCorrssim, 'jitter','on', 'jitterAmount', 0.1, 'MarkerEdgeColor', [1, .25, 0]);
                        hold on

                        set(gca,'xtick',[positions(1) positions(2) positions(3)])
                        set(gca,'xticklabel',{'Between State Comparison', 'Within State Comparison', 'Across Subject Comparison'}, 'FontSize',12)
                        ylabel('Dice Correlation', 'FontSize',12)

                        color = ['k', 'g', 'b'];
                        %color = ['k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k'];
                        h = findobj(gca,'Tag','Box');

                        for j=1:length(h)

                            patch(get(h(j),'XData'),get(h(j),'YData'),color(j), 'FaceColor', 'none', 'EdgeColor', color(j));

                        end

                        c = get(gca, 'Children');

                        m = findobj(gca,'Type','scatter');

                        %hleg1 = legend([m(2:10); c(1:3)], 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Between State Comparison', 'Within State Comparison', 'Between Subject Comparison', 'Location', 'NorthEast');
                        hleg1 = legend([m(2:10)], 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', 'NorthEast');
                        %hleg1 = legend(m(2:10), 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', 'NorthEast');
                        hleg1.FontSize = 14;
                        ylim([0 1]);

                        %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

                        %hleg1 = legend(c(1:3), 'Between State Comparison', 'Within State Comparison', 'Between Subject Comparison', 'Location', 'NorthWest');

                        %legendpos = 'NorthEast';

                        %hleg2 = legend(m(2:10), 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', legendpos);

                    else
                        
                        x = [];
                        group = [];
                        positions = [1 1.25 1.5 2 2.25 2.5 3 3.25 3.5 4 4.25 4.5 5 5.25 5.5 6 6.25 6.5 7 7.25 7.5 8 8.25 8.5 9 9.25 9.5 10 10.25 10.5];
                        positionssubs = [1 2 3 4 5 6 7 8 9 10];
                        scatterposbetween = [sort(repmat([1:1:9],1,2)) repmat(10,1,size(DiceCorrsTaskRest,1)*size(DiceCorrsTaskRest,2))];
                        scatterposwithin = [sort(repmat([1.25:1:9.25],1,2)) repmat(10.25,1,size(DiceCorrsTaskTask,1)+size(DiceCorrsRestRest,1))];
                        scatterposacross = [sort(repmat([1.5:1:9.5],1,16)) repmat(10.5,1,size(DiceCorrssim,1))];

                        for g = 1:nsubs+1
                            
                            if g < nsubs+1

                                x = [x [DiceCorrsTaskRest(g,:) DiceCorrsTaskTask(g,:) DiceCorrsRestRest(g,:) DiceCorrSubs(:,1)']];
                                
                            else
                                
                                x = [x [reshape(DiceCorrsTaskRest,1,size(DiceCorrsTaskRest,1)*size(DiceCorrsTaskRest,2)) DiceCorrsTaskTask' DiceCorrsRestRest' DiceCorrssim']];
                                
                            end
                        end

                        for i = 1:length(positionssubs)
                            
                            if i < length(positionssubs)

                                group = [group [repmat(positionssubs(i),1,size(DiceCorrsTaskRest,2)) repmat(positionssubs(i)+.25,1,size(DiceCorrsTaskTask,2)+size(DiceCorrsRestRest,2)) repmat(positionssubs(i)+.5,1,size(DiceCorrSubs,1),1)]];
                                
                            else
                            
                                group = [group [repmat(positionssubs(i),1,size(DiceCorrsTaskRest,1)*size(DiceCorrsTaskRest,2)) repmat(positionssubs(i)+.25,1,size(DiceCorrsTaskTask,1)+size(DiceCorrsRestRest,1)) repmat(positionssubs(i)+.5,1,size(DiceCorrssim,1),1)]];
                                
                            end
                        end

                        %group = [repmat(1,1,size(DiceCorrsTaskRest,2)+size(DiceCorrsTaskTask,2)+size(DiceCorrsRestRest,2)+size(DiceCorrSubs,1)) repmat(2,1,length(DiceCorrsTaskTask)+length(DiceCorrsRestRest)) repmat(3,1,length(DiceCorrssim))];
                        %positions = [1 2 3 4 5 6 7 8 9];
                        %scatterpos = [1 2];

                        boxplot(x,group, 'BoxStyle', 'outline', 'positions', positions, 'Symbol', '');
                        hold on
                        scatter(scatterposbetween, [DiceCorrsTaskRest(1,:) DiceCorrsTaskRest(2,:) DiceCorrsTaskRest(3,:) DiceCorrsTaskRest(4,:) DiceCorrsTaskRest(5,:) DiceCorrsTaskRest(6,:) DiceCorrsTaskRest(7,:) DiceCorrsTaskRest(8,:) DiceCorrsTaskRest(9,:) reshape(DiceCorrsTaskRest,1,size(DiceCorrsTaskRest,1)*size(DiceCorrsTaskRest,2))], 'filled', 'MarkerFaceColor', 'b');
                        hold on
                        scatter(scatterposwithin, [DiceCorrsRestRest(1,:) DiceCorrsTaskTask(1,:) DiceCorrsRestRest(2,:) DiceCorrsTaskTask(2,:) DiceCorrsRestRest(3,:) DiceCorrsTaskTask(3,:) DiceCorrsRestRest(4,:) DiceCorrsTaskTask(4,:) DiceCorrsRestRest(5,:) DiceCorrsTaskTask(5,:) DiceCorrsRestRest(6,:) DiceCorrsTaskTask(6,:) DiceCorrsRestRest(7,:) DiceCorrsTaskTask(7,:) DiceCorrsRestRest(8,:) DiceCorrsTaskTask(8,:) DiceCorrsRestRest(9,:) DiceCorrsTaskTask(9,:) DiceCorrsRestRest' DiceCorrsTaskTask'], 'filled', 'MarkerFaceColor', 'g');
                        hold on
                        scatter(scatterposacross, [DiceCorrSubs(:,1)' DiceCorrSubs(:,2)' DiceCorrSubs(:,3)' DiceCorrSubs(:,4)' DiceCorrSubs(:,5)' DiceCorrSubs(:,6)' DiceCorrSubs(:,7)' DiceCorrSubs(:,8)' DiceCorrSubs(:,9)' DiceCorrssim'], 'jitter','on', 'jitterAmount', 0.1, 'MarkerEdgeColor', [1, .25, 0]);
                        hold on
                        set(gca,'xtick',[mean(positions(1:3)) mean(positions(4:6)) mean(positions(7:9)) mean(positions(10:12)) mean(positions(13:15)) mean(positions(16:18)) mean(positions(19:21)) mean(positions(22:24)) mean(positions(25:27)) mean(positions(28:30))])
                        set(gca,'xticklabel',{'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Group Average'}, 'FontSize',12)
                        ylabel('Dice Correlation', 'FontSize',12)

                        color = ['k', 'g', 'b', 'k', 'g', 'b', 'k', 'g', 'b', 'k', 'g', 'b', 'k', 'g', 'b', 'k', 'g', 'b', 'k', 'g', 'b', 'k', 'g', 'b', 'k', 'g', 'b', 'k', 'g', 'b', 'k', 'g', 'b'];

                        h = findobj(gca,'Tag','Box');

                        for j=1:length(h)

                            patch(get(h(j),'XData'),get(h(j),'YData'),color(j), 'FaceColor', 'none', 'EdgeColor', color(j));

                        end

                        c = get(gca, 'Children');

                        hleg1 = legend(c(1:3), 'Between Subject Comparison', 'Within State Comparison', 'Between State Comparison', 'Location', 'NorthEast');
                        hleg1.FontSize = 14;
                        ylim([0 1]);
                        
                        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.4, 0.65, 0.4, 0.65]);
                        
                    end

               	end
            
            else
                
                if AnyOverlap == 1 && COMOverlap == 1
        
                    x = [DiceCorrsTaskRest DiceCorrsTaskTask DiceCorrsRestRest NumOverlapTaskRest NumOverlapTaskTask NumOverlapRestRest COMOverlapsTaskRest COMOverlapsTaskTask COMOverlapsRestRest];
                    group = [repmat(1,1,length(DiceCorrsTaskRest)) repmat(2,1,length(DiceCorrsTaskTask)+length(DiceCorrsRestRest)) repmat(3,1,length(NumOverlapTaskRest)) repmat(4,1,length(NumOverlapTaskTask)+length(NumOverlapRestRest)) repmat(5,1,length(COMOverlapsTaskRest)) repmat(6,1,length(COMOverlapsTaskTask)+length(COMOverlapsRestRest))];
                    positions = [1 1.25 2 2.25 3 3.25];
                    scatterpos = [1 2 3 1.25 2.25 3.25];
                    boxplot(x,group, 'BoxStyle', 'outline', 'positions', positions, 'Symbol', '');
                    hold on
                    scatter(scatterpos, [plotdata(9,:) mean([plotdata2(9,1) plotdata2(18,1)]) mean([plotdata2(9,2) plotdata2(18,2)]) mean([plotdata2(9,3) plotdata2(18,3)])], 'filled', 'MarkerFaceColor', [1, 0.5, 0]);
                    hold on
                    scatter(scatterpos, [plotdata(8,:) mean([plotdata2(8,1) plotdata2(17,1)]) mean([plotdata2(8,2) plotdata2(17,2)]) mean([plotdata2(8,3) plotdata2(17,3)])], 'filled', 'MarkerFaceColor', [0, 0.6, 0.6]);
                    hold on    
                    scatter(scatterpos, [plotdata(7,:) mean([plotdata2(7,1) plotdata2(16,1)]) mean([plotdata2(7,2) plotdata2(16,2)]) mean([plotdata2(7,3) plotdata2(16,3)])], 'filled', 'MarkerFaceColor', [1, 0, 1]);
                    hold on
                    scatter(scatterpos, [plotdata(6,:) mean([plotdata2(6,1) plotdata2(15,1)]) mean([plotdata2(6,2) plotdata2(15,2)]) mean([plotdata2(6,3) plotdata2(15,3)])], 'filled', 'MarkerFaceColor', [0.2, 1, 1]);
                    hold on
                    scatter(scatterpos, [plotdata(5,:) mean([plotdata2(5,1) plotdata2(14,1)]) mean([plotdata2(5,2) plotdata2(14,2)]) mean([plotdata2(5,3) plotdata2(14,3)])], 'filled', 'MarkerFaceColor', [0, 0, 1]);
                    hold on
                    scatter(scatterpos, [plotdata(4,:) mean([plotdata2(4,1) plotdata2(13,1)]) mean([plotdata2(4,2) plotdata2(13,2)]) mean([plotdata2(4,3) plotdata2(13,3)])], 'filled', 'MarkerFaceColor', [1, 0, 0]);
                    hold on
                    scatter(scatterpos, [plotdata(3,:) mean([plotdata2(3,1) plotdata2(12,1)]) mean([plotdata2(3,2) plotdata2(12,2)]) mean([plotdata2(3,3) plotdata2(12,3)])], 'filled', 'MarkerFaceColor', [0, 1, 0]);
                    hold on
                    scatter(scatterpos, [plotdata(2,:) mean([plotdata2(2,1) plotdata2(11,1)]) mean([plotdata2(2,2) plotdata2(11,2)]) mean([plotdata2(2,3) plotdata2(11,3)])], 'filled', 'MarkerFaceColor', [0.9, 0.9, 0]);
                    hold on
                    scatter(scatterpos, [plotdata(1,:) mean([plotdata2(1,1) plotdata2(10,1)]) mean([plotdata2(1,2) plotdata2(10,2)]) mean([plotdata2(1,3) plotdata2(10,3)])], 'filled', 'MarkerFaceColor', [0, 0, 0]);
                    hold on

                    set(gca,'xtick',[mean(positions(1:2)) mean(positions(3:4)) mean(positions(5:6)) ])
                    set(gca,'xticklabel',{'Dice Correlations', 'Any Overlap', 'COM Overlap'}, 'FontSize',12)

                    color = ['b', 'g', 'b', 'g', 'b', 'g'];
                    h = findobj(gca,'Tag','Box');
        
                    for j=1:length(h)
            
                        patch(get(h(j),'XData'),get(h(j),'YData'),color(j), 'FaceColor', 'none', 'EdgeColor', color(j));
            
                    end

                    c = get(gca, 'Children');
            
                    m = findobj(gca,'Type','scatter');
            
                    hleg1 = legend([m(1:9); c(1:2)], 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Between State Comparison', 'Within State Comparison', 'Location', 'NorthEast');
                    hleg1.FontSize = 14;
                    ylim([0 1]);
            
                    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
            
                    %hleg1 = legend(c(1:2), 'Between State Comparison', 'Within State Comparison' );
                    
                else
                    
                    x = [DiceCorrsTaskRest DiceCorrsTaskTask DiceCorrsRestRest];
                    group = [repmat(1,1,length(DiceCorrsTaskRest)) repmat(2,1,length(DiceCorrsTaskTask)+length(DiceCorrsRestRest))];
                    positions = [1 2];
                    scatterpos = [1 2];
                    boxplot(x,group, 'BoxStyle', 'outline', 'positions', positions, 'Symbol', '');
                    hold on
                    scatter(scatterpos, [plotdata(9,:) mean([plotdata2(9,1) plotdata2(18,1)])], 'filled', 'MarkerFaceColor', [1, 0.5, 0]);
                    hold on
                    scatter(scatterpos, [plotdata(8,:) mean([plotdata2(8,1) plotdata2(17,1)])], 'filled', 'MarkerFaceColor', [0, 0.6, 0.6]);
                    hold on    
                    scatter(scatterpos, [plotdata(7,:) mean([plotdata2(7,1) plotdata2(16,1)])], 'filled', 'MarkerFaceColor', [1, 0, 1]);
                    hold on
                    scatter(scatterpos, [plotdata(6,:) mean([plotdata2(6,1) plotdata2(15,1)])], 'filled', 'MarkerFaceColor', [0.2, 1, 1]);
                    hold on
                    scatter(scatterpos, [plotdata(5,:) mean([plotdata2(5,1) plotdata2(14,1)])], 'filled', 'MarkerFaceColor', [0, 0, 1]);
                    hold on
                    scatter(scatterpos, [plotdata(4,:) mean([plotdata2(4,1) plotdata2(13,1)])], 'filled', 'MarkerFaceColor', [1, 0, 0]);
                    hold on
                    scatter(scatterpos, [plotdata(3,:) mean([plotdata2(3,1) plotdata2(12,1)])], 'filled', 'MarkerFaceColor', [0, 1, 0]);
                    hold on
                    scatter(scatterpos, [plotdata(2,:) mean([plotdata2(2,1) plotdata2(11,1)])], 'filled', 'MarkerFaceColor', [0.9, 0.9, 0]);
                    hold on
                    scatter(scatterpos, [plotdata(1,:) mean([plotdata2(1,1) plotdata2(10,1)])], 'filled', 'MarkerFaceColor', [0, 0, 0]);
                    hold on

                    set(gca,'xtick',[positions(1) positions(2)])
                    set(gca,'xticklabel',{'Between Subject Comparison', 'Within Subject Comparison'}, 'FontSize',12)

                    color = ['b', 'g'];
                    h = findobj(gca,'Tag','Box');
        
                    for j=1:length(h)
            
                        patch(get(h(j),'XData'),get(h(j),'YData'),color(j), 'FaceColor', 'none', 'EdgeColor', color(j));
            
                    end

                    c = get(gca, 'Children');
            
                    m = findobj(gca,'Type','scatter');
            
                    hleg1 = legend([m(1:9); c(1:2)], 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Between State Comparison', 'Within State Comparison', 'Location', 'NorthEast');
                    hleg1.FontSize = 14;
                    ylim([0 1]);
            
                    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
            
                    %hleg1 = legend(c(1:2), 'Between State Comparison', 'Within State Comparison' );
                    
                end
            end
        
            if FullMaps == 1
        
                filename = '/AllSubjects_Boxplot_AllData_Overlap_SNRExclude.jpg';
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end

            elseif MatchedMaps == 1 && SNRExclude == 1 && SizeExclude == 1 && SplitHalf == 1 && numel(thresholds) > 1 && AbsoluteThresholds == 1
                
                filename = ['/AllSubjects_Boxplot_MatchedData_SplitHalf_Overlap_SNRAndSizeExclude_' num2str(thresholds(v)) '_Abs.jpg'];
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end
                
            elseif MatchedMaps == 1 && SNRExclude == 1 && SizeExclude == 1 && SplitHalf == 1 && FinalFigure == 1
                
             	filename = ['/AllSubjects_Boxplot_FinalPlotForPaper_' num2str(thresholds(v)) '_Percent.jpg'];
                
            elseif MatchedMaps == 1 && SNRExclude == 1 && SizeExclude == 1 && SplitHalf == 1 && numel(thresholds) > 1
                
                filename = ['/AllSubjects_Boxplot_MatchedData_SplitHalf_Overlap_SNRAndSizeExclude_' num2str(thresholds(v)) '_Percent.jpg'];
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end
                
            elseif MatchedMaps == 1 && SNRExclude == 1 && SizeExclude == 1 && SplitHalf == 1 && permvals == 1
        
                filename = ['/AllSubjects_Boxplot_MatchedData_Overlap_Perm_SNRAndSizeExclude_' num2str(thresholds(v)) '.jpg'];
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end

            elseif MatchedMaps == 1 && SNRExclude == 1 && SizeExclude == 1 && SplitHalf == 1
        
                filename = '/AllSubjects_Boxplot_MatchedData_SplitHalf_Overlap_SNRAndSizeExclude.jpg';
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end
      
            elseif MatchedMaps == 1 && SNRExclude == 1 && SizeExclude == 1
        
                filename = '/AllSubjects_Boxplot_MatchedData_Overlap_SNRAndSizeExclude.jpg';
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end

            elseif MatchedMaps == 1 && SNRExclude == 1 && SplitHalf == 1 && numel(thresholds) > 1 && AbsoluteThresholds == 1
                
                filename = ['/AllSubjects_Boxplot_MatchedData_SplitHalf_Overlap_SNRExclude_' num2str(thresholds(v)) '_Abs.jpg'];
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end

            elseif MatchedMaps == 1 && SNRExclude == 1 && SplitHalf == 1 && numel(thresholds) > 1
                
                filename = ['/AllSubjects_Boxplot_MatchedData_SplitHalf_Overlap_SNRExclude_' num2str(thresholds(v)) '_Percent.jpg'];
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end
                
            elseif MatchedMaps == 1 && SNRExclude == 1 && SplitHalf == 1 && permvals == 1
        
                filename = ['/AllSubjects_Boxplot_MatchedData_SplitHalf_Overlap_Perm_SNRExclude_' num2str(thresholds) '.jpg'];
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end

            elseif MatchedMaps == 1 && SNRExclude == 1 && SplitHalf == 1
        
                filename = '/AllSubjects_Boxplot_MatchedData_SplitHalf_Overlap_SNRExclude.jpg';
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end

            elseif MatchedMaps == 1 && SNRExclude == 1
        
                filename = '/AllSubjects_Boxplot_MatchedData_Overlap_SNRExclude.jpg';
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end

            elseif MatchedMaps == 1 && SizeExclude == 1 && SplitHalf == 1
        
                filename = '/AllSubjects_Boxplot_MatchedData_SplitHalf_Overlap_SizeExclude.jpg';
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end

            elseif MatchedMaps == 1 && SizeExclude == 1
        
                filename = '/AllSubjects_Boxplot_MatchedData_Overlap_SizeExclude.jpg';
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end
        
            else
        
                filename = '/AllSubjects_Boxplot_MatchedData_Overlap_NoExclude.jpg';
        
                if randomizevals == 1
            
                    filename = strrep(filename, '.jpg', '_pvals.jpg');
            
                end
                            
            end
            
            if BoxPlotBySubject == 1
                
                filename = strrep(filename, 'Boxplot', 'BoxplotBySubject');
                
            end
            
            if AnyOverlap == 1 && COMOverlap == 1
                
                print(gcf,[outputdir filename],'-dpng','-r300');
                
                %saveas(gcf,[outputdir filename])
                
            else
                
                filename = strrep(filename, '.jpg', '_DiceOnly.jpg');
                
                print(gcf,[outputdir filename],'-dpng','-r300');
                
                %saveas(gcf,[outputdir filename])
                
            end
    
            close gcf
        
        
        end
                    
    end
    
end
    