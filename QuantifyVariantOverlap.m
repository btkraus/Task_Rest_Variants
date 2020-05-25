%% QuantifyVariantOverlap.m 
%This script compares variant locations between rest-rest, task-task, rest-task
%It makes plots that quantify how much they overlap between conditions
%Written by Brian Kraus. Edited by Diana Perez.

clear all

%% Paths
% Specify output directories    
outputdir = '/Users/dianaperez/Box/Latest_Analysis_Replication/Variant_Overlap_Plots/Thresholded_Plots/';
dirpath = '/Users/dianaperez/Documents/GitHub/Task_Rest_Variants/'; %location of txt files
% names of txt files
spCorrRestEvenTxt = 'MSC_rest_spCorrMaps_Even.txt';
spCorrRestOddTxt = 'MSC_rest_spCorrMaps_Odd.txt';
spCorrTaskEvenTxt = 'MSC_task_spCorrMaps_Even.txt';
spCorrTaskOddTxt = 'MSC_task_spCorrMaps_Odd.txt';
varMapsRestEvenTxt = 'MSC_rest_varMaps_Even.txt';
varMapsRestOddTxt = 'MSC_rest_varMaps_Odd.txt';
varMapsTaskEvenTxt = 'MSC_task_varMaps_Even.txt';
varMapsTaskOddTxt = 'MSC_task_varMaps_Odd.txt';
if ~isfolder(outputdir)
    mkdir(outputdir) % creates output directory if it doesn't already exist
end
%% Specify these conditions
DiceCorr = 1;  %% Toggles whether to calculate dice correlations for all subjects
plotresults = 1;  %% Toggles whether to plot results for each comparison
FullMaps = 0;  %% Toggles whether to use all task and rest data for each subject
SplitHalf = 1;  %% Toggles whether to calculate a split-half 
MatchedMaps = 1;  %% Toggles whether to use a matched amount of task and rest data for each subject
SNRExclude = 1;  %% Toggles whether to use data excluded for low SNR
randomizevals = 1;  %% Toggles whether to calculate a null distribution across subjects
permvals = 1; %% Toggles whether to calculate a permutation of all possible combinations instead of bootstrapping
BoxPlot = 1;  %% Toggles whether to plot a box plot for the group level results
BoxPlotBySubject = 1;  %% Toggles whether to display boxplot by subjects instead of by comparison
FinalFigure = 1;  %% Toggles whether to plot the final figure (Figure 1) for the paper
AbsoluteThresholds = 0;  %% Toggles whether thresholds are absolute values or percents
thresholds = [2.5];  %% Sets the thresholds of files to load
%% Create variables for number of variants/average size    
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

% this for-loop sets up variables to run dice correlations between task and
% rest data for each threshold
for v = 1:numel(thresholds)    
    % Create temp variables
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
    alltaskfilestaskeven = [];
    alltaskfilestaskodd = [];
    alltaskfilesresteven = [];
    alltaskfilesrestodd = [];
    DiceCorrsTaskRest = [];
    DiceCorrsTaskTask = [];
    DiceCorrsRestRest = [];

    % read .txt files 
    [task_files_even, ~, ~] = textread([dirpath spCorrTaskEvenTxt],'%s%s%s');
    [rest_files_even, ~, ~] = textread([dirpath spCorrRestEvenTxt],'%s%s%s');    
    [task_files_odd, subjects1, tasks1] = textread([dirpath spCorrTaskOddTxt],'%s%s%s');
    [rest_files_odd, subjects2, tasks2] = textread([dirpath spCorrRestOddTxt],'%s%s%s');
    [task_masks_even, ~, ~] = textread([dirpath varMapsTaskEvenTxt],'%s%s%s');
    [rest_masks_even, ~, ~] = textread([dirpath varMapsRestEvenTxt],'%s%s%s');
    [task_masks_odd, sub1, t1] = textread([dirpath varMapsTaskOddTxt],'%s%s%s');
    [rest_masks_odd, sub2, t2] = textread([dirpath varMapsRestOddTxt],'%s%s%s');

% sets up number of files for for-loop
nfiles = length(rest_files_even);   

    for x = 1:nfiles

        subject = subjects2{x};

        %% Reads cifti files from .txt file paths
        % spatial correlation maps
        cifti_rest_even = ft_read_cifti_mod(rest_files_even{x});
        cifti_task_even = ft_read_cifti_mod(task_files_even{x});
        cifti_rest_odd = ft_read_cifti_mod(rest_files_odd{x});
        cifti_task_odd = ft_read_cifti_mod(task_files_odd{x});
        % variant maps  
        cifti_rest_mask_even = ft_read_cifti_mod(rest_masks_even{x});
        cifti_task_mask_even = ft_read_cifti_mod(task_masks_even{x});
        cifti_rest_mask_odd = ft_read_cifti_mod(rest_masks_odd{x});
        cifti_task_mask_odd = ft_read_cifti_mod(task_masks_odd{x});

        %% Leaves only correlations for vertices belonging to variants
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


        %% Get number of variants and mean/median size, (put it at end of variant making script)        
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
        %% end of descriptive stats part

        %% Dice correlations
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

        % if empty (due to low threshold or high size exclusion), then set dice corr to 0
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
    end

    %% Saves SNR excluded data for each subject if you want to do permutations
    if randomizevals == 1 
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
    %% Add temp variables to final variables, put these with descriptive stats

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

    %% CHECK STARTING HERE!!!! COMPARE TO OLD SCRIPT
    if randomizevals == 1
        if permvals == 1
            
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
        else
            nsim = 1000;
            nsubs = nfiles;
            DiceCorrssim = zeros(1,nsim);
            DiceCorrssimAll = zeros(nsubs,nsim);
            tasktaskshuffleall = [];
            taskrestshuffleall = [];

            for x = 1:nsim
                alltaskdc = [];
        
                while true      %% Make sure no subjects match each other            
                    reshuffle = 0;    
                    rng('shuffle');  %% Reset random number generator on each iteration
        
                    % Randomly shuffle rows (subjects) for each task and rest matrix        
                    tasktaskshuffle = 1:nsubs;
                    taskrestshuffle = 1:nsubs;        
                    tasktaskshuffle = tasktaskshuffle(randperm(length(1:nsubs)));
                    taskrestshuffle = taskrestshuffle(randperm(length(1:nsubs)));
            
                    for h = 1:length(tasktaskshuffle)                
                        if reshuffle == 0 && tasktaskshuffle(h) == taskrestshuffle(h)                    
                            reshuffle = 1;                    
                        end
                    end
            
                    if reshuffle == 0  %% Break loop if all subjects are shuffled                
                        tasktaskshuffleall = [tasktaskshuffleall tasktaskshuffle'];
                        taskrestshuffleall = [taskrestshuffleall taskrestshuffle'];                
                        break                
                    end
                end
        
                if SplitHalf == 1            
                    alltasktaskshuffle = alltaskfilestaskodd(tasktaskshuffle,:);
                    alltaskrestshuffle = alltaskfilesresteven(taskrestshuffle,:);                    
                else        
                    alltasktaskshuffle = alltaskfilestask(tasktaskshuffle,:);
                    alltaskrestshuffle = alltaskfilesrest(taskrestshuffle,:);
                end
        
                DiceCorrssimtemp = zeros(1,size(alltasktaskshuffle,1));
                NumOverlapsimtemp = zeros(1,size(alltasktaskshuffle,1));
                COMOverlapssimtemp = zeros(1,size(alltasktaskshuffle,1));
    
                for d = 1:size(alltasktaskshuffle,1)
            
                    vars_rest = unique(alltaskrestshuffle(d,:));
                    vars_rest(1) = [];
                    vars_task = unique(alltasktaskshuffle(d,:));
                    vars_task(1) = [];
    
                    if DiceCorr == 1        
                        dcorrdata = [];        
                        for q = 1:length(alltaskrestshuffle(d,:))            
                            if alltaskrestshuffle(d,q) > 0 && alltasktaskshuffle(d,q) > 0                
                                dcorrdata = [dcorrdata;1 1];                
                            elseif alltaskrestshuffle(d,q) > 0                
                                dcorrdata = [dcorrdata;1 0];                
                            elseif alltasktaskshuffle(d,q) > 0                
                                dcorrdata = [dcorrdata;0 1];                
                            end            
                        end
                    
                        if isempty(dcorrdata                        
                            dc = 0;                        
                        else        
                            dc = dice_coefficient_mod(dcorrdata(:,1),dcorrdata(:,2));                        
                        end        
                        DiceCorrssimtemp(d) = dc;    
                    end
                end
    
                DiceCorrssim(x) = mean(DiceCorrssimtemp);       
                DiceCorrssimAll(:,x) = DiceCorrssimtemp';        
        
                if mod(x,100) == 0           
                    disp([num2str(x) ' iterations completed'])           
                end    
            end
    
            DiceCorrspval = prctile(DiceCorrssim, 95);
            DiceCorrssub = zeros(nsubs,1);      %% Calculate averaged values per subject for (pseudo) effect size
    
            for m = 1:size(DiceCorrssimAll,2)        
                for n = 1:size(DiceCorrssimAll,1)            
                    for o = 1:nsubs            
                        if tasktaskshuffleall(n,m) == o || taskrestshuffleall(n,m) == o                
                            DiceCorrssub(o) = DiceCorrssub(o) + DiceCorrssimAll(n,m);               
                        end
                    end
                end
            end
            
            DiceCorrssub = DiceCorrssub./(size(DiceCorrssimAll,2)*2);
            
        end

        DiceCorrspval = prctile(DiceCorrssim, 95);
    
        %% Plotting Results

        if plotresults == 1
            if BoxPlot == 1
                if SplitHalf == 1   
                    plotdata = DiceCorrsTaskRest;
                    plotdata2 = [DiceCorrsTaskTask;DiceCorrsRestRest];
                else
                    plotdata = DiceCorrs;
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
                                x = [x 0];
                            else
                                x = [x [reshape(DiceCorrsTaskRest,1,size(DiceCorrsTaskRest,1)*size(DiceCorrsTaskRest,2)) DiceCorrsTaskTask' DiceCorrsRestRest' DiceCorrssim']];
                            end
                        end

                        for i = 1:length(positionssubs)
                            if i < length(positionssubs)
                                group = [group positionssubs(i)];
                            else
                                group = [group [repmat(positionssubs(i),1,size(DiceCorrsTaskRest,1)*size(DiceCorrsTaskRest,2)) repmat(positionssubs(i)+.25,1,size(DiceCorrsTaskTask,1)+size(DiceCorrsRestRest,1)) repmat(positionssubs(i)+.5,1,size(DiceCorrssim,1),1)]];
                            end
                        end

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
                            h = findobj(gca,'Tag','Box');

                            for j=1:length(h)
                                patch(get(h(j),'XData'),get(h(j),'YData'),color(j), 'FaceColor', 'none', 'EdgeColor', color(j));
                            end

                            c = get(gca, 'Children');
                            m = findobj(gca,'Type','scatter');
                            hleg1 = legend([m(2:10)], 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', 'NorthEast');
                            hleg1.FontSize = 14;
                            ylim([0 1]);

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
                end

                if FullMaps == 1
                    filename = '/AllSubjects_Boxplot_AllData_Overlap_SNRExclude.jpg';
                    if randomizevals == 1
                        filename = strrep(filename, '.jpg', '_pvals.jpg');
                    end
                elseif MatchedMaps == 1 && SNRExclude == 1 && SplitHalf == 1 && numel(thresholds) > 1 && AbsoluteThresholds == 1
                    filename = ['/AllSubjects_Boxplot_MatchedData_SplitHalf_Overlap_SNRAndSizeExclude_' num2str(thresholds(v)) '_Abs.jpg'];
                    if randomizevals == 1
                        filename = strrep(filename, '.jpg', '_pvals.jpg');
                    end
                elseif MatchedMaps == 1 && SNRExclude == 1 && SplitHalf == 1 && FinalFigure == 1
                    filename = ['/AllSubjects_Boxplot_FinalPlotForPaper_' num2str(thresholds(v)) '_Percent.jpg'];
                elseif MatchedMaps == 1 && SNRExclude == 1 && SplitHalf == 1 && numel(thresholds) > 1
                    filename = ['/AllSubjects_Boxplot_MatchedData_SplitHalf_Overlap_SNRAndSizeExclude_' num2str(thresholds(v)) '_Percent.jpg'];
                    if randomizevals == 1
                        filename = strrep(filename, '.jpg', '_pvals.jpg');
                    end
                elseif MatchedMaps == 1 && SNRExclude == 1 && SplitHalf == 1 && permvals == 1
                    filename = ['/AllSubjects_Boxplot_MatchedData_Overlap_Perm_SNRAndSizeExclude_' num2str(thresholds(v)) '.jpg'];
                    if randomizevals == 1
                        filename = strrep(filename, '.jpg', '_pvals.jpg');
                    end
                elseif MatchedMaps == 1 && SNRExclude == 1 && SplitHalf == 1
                    filename = '/AllSubjects_Boxplot_MatchedData_SplitHalf_Overlap_SNRAndSizeExclude.jpg';
                    if randomizevals == 1
                        filename = strrep(filename, '.jpg', '_pvals.jpg');
                    end
                elseif MatchedMaps == 1 && SNRExclude == 1 
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
                elseif MatchedMaps == 1 && SplitHalf == 1
                    filename = '/AllSubjects_Boxplot_MatchedData_SplitHalf_Overlap_SizeExclude.jpg';
                    if randomizevals == 1
                        filename = strrep(filename, '.jpg', '_pvals.jpg');
                    end
                elseif MatchedMaps == 1 
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

                filename = strrep(filename, '.jpg', '_DiceOnly.jpg');
                print(gcf,[outputdir filename],'-dpng','-r300');

                close gcf
            end
        end
    end
end