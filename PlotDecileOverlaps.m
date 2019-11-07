

clear all

SNRExclude = 1;        %% Toggles whether to exclude low SNR areas
SplitHalf = 1;        %% Toggles whether to use split-half data
AbsoluteThresholds = 0;        %% Toggles whether threshold is set by percentile or absolute value
ResidualDeciles = 0;        %% Matches data left over from each decile in subsequent deciles
PercentOverlap = 0;     %% Toggles whether to calculate percent overlap in each decile
% deciles = [0 10;      %% Percentage Deciles
%            10 20;
%            20 30;
%            30 40;
%            40 50;
%            50 60;
%            60 70;
%            70 80;
%            80 90;
%            90 100];

deciles = [0 2.5;      %% Percentage Deciles
           2.5 5;
           5 7.5;
           7.5 10;
           10 12.5;
           12.5 15;
           15 17.5;
           17.5 20;
           20 22.5;
           22.5 25
           25 27.5
           27.5 30
           30 32.5
           32.5 35
           35 37.5
           37.5 40
           40 42.5
           42.5 45
           45 47.5
           47.5 50
           50 52.5
           52.5 55
           55 57.5
           57.5 60
           60 62.5
           62.5 65
           65 67.5
           67.5 70
           70 72.5
           72.5 75
           75 77.5
           77.5 80
           80 82.5
           82.5 85
           85 87.5
           87.5 90
           90 92.5
           92.5 95
           95 97.5
           97.5 100];

% deciles = [-1 .1;       %% Absolute Deciles
%            .1 .2;
%            .2 .3;
%            .3 .4;
%            .4 .5;
%            .5 .6;
%            .6 .7;
%            .7 1];

% deciles = [-1 .1;       %% Absolute Deciles
%            .1 .125;
%            .125 .15;
%            .15 .175;
%            .175 .2;
%            .2 .225;
%            .225 .25;
%            .25 .275;
%            .275 .3;
%            .3 .325;
%            .325 .35;
%            .35 .375;
%            .375 .4;
%            .4 .425;
%            .425 .45;
%            .45 .475;
%            .475 .5;
%            .5 .525;
%            .525 .55;
%            .55 .575;
%            .575 .6;
%            .6 .625;
%            .625 .65;
%            .65 .675;
%            .675 .7;
%            .7 1];

% deciles = [-1 .1;       %% Absolute Deciles
%            .1 .15;
%            .15 .2;
%            .2 .25;
%            .25 .3;
%            .3 .35;
%            .35 .4;
%            .4 .45;
%            .45 .5;
%            .5 .55;
%            .55 .6;
%            .6 .65;
%            .65 .7;
%            .7 1];
       
%deciles = [-1 .45];       %% Absolute Deciles
%            .3 .6;
%            .6 .9];
           
outputdir = '/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Decile_Variant_Plots';
       
       
if SplitHalf == 1
    
    if PercentOverlap == 1
       
        pctoverlaptaskrest = [];
        pctoverlaprestrest = [];
        pctoverlaptasktask = [];
        
    end
    
	DiceCorrsTaskRest = [];
	DiceCorrsTaskTask = [];
	DiceCorrsRestRest = [];
        
else
    
	DiceCorrs = [];
        
end

if SplitHalf == 1
        
	alltaskfilestaskeven = {};
	alltaskfilestaskodd = {};
	alltaskfilesresteven = {};
	alltaskfilesrestodd = {};
        
else
                        
	alltaskfilestask = {};
	alltaskfilesrest = {};   
        
end

if AbsoluteThresholds ~= 1
    
    cifti_rest_even_vals = [];
    cifti_task_even_vals = [];
    cifti_rest_odd_vals = [];
    cifti_task_odd_vals = [];
    
end

if ResidualDeciles == 1
    
    ResidVerticesRestRest = cell(9,2);
    ResidVerticesTaskTask = cell(9,2);
    ResidVerticesRestTask = cell(9,2);
    
    MatchedVerticesRestRest = cell(9,2);
    MatchedVerticesTaskTask = cell(9,2);
    MatchedVerticesRestTask = cell(9,2);
    
end

for g = 1:size(deciles,1)

    if SplitHalf == 1 && SNRExclude == 1

        [task_files_even, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_splithalf_even_matched.txt','%s%s%s');

        [rest_files_even, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_splithalf_even_matched.txt','%s%s%s');
    
        [task_files_odd, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_splithalf_odd_matched.txt','%s%s%s');

        [rest_files_odd, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_splithalf_odd_matched.txt','%s%s%s');
        

    end


    subs = unique(subjects1);

    if SplitHalf == 1
        
        if PercentOverlap == 1
        
            pctoverlapcounttaskrestTemp = [];
            pctoverlapcountrestrestTemp = [];
            pctoverlapcounttasktaskTemp = [];
            
        end
        
        DiceCorrsTaskRestTemp = [];
        DiceCorrsTaskTaskTemp = [];
        DiceCorrsRestRestTemp = [];

    else
    
        DiceCorrsTemp = [];
        
    end
    
    if SplitHalf == 1
        
        alltaskfilestaskeventemp = [];
        alltaskfilestaskoddtemp = [];
        alltaskfilesresteventemp = [];
        alltaskfilesrestoddtemp = [];
        
    else
                        
        alltaskfilestasktemp = [];
        alltaskfilesresttemp = [];   
        
    end
    
	if AbsoluteThresholds ~= 1
    
        cifti_rest_even_valsSub = [];
        cifti_task_even_valsSub = [];
        cifti_rest_odd_valsSub = [];
        cifti_task_odd_valsSub = [];
    
	end

    for s = 1:numel(subs)
    
        subject = subs{s};
        
        if ResidualDeciles == 1 && (~isempty(ResidVerticesRestRest{s,1}) || ~isempty(ResidVerticesTaskTask{s,1}) || ~isempty(ResidVerticesRestTask{s,1}))
            
            ResidVerticesRestRest{s,1} = unique(ResidVerticesRestRest{s,1});
            ResidVerticesTaskTask{s,1} = unique(ResidVerticesTaskTask{s,1});
            ResidVerticesRestTask{s,1} = unique(ResidVerticesRestTask{s,1});
            
            ResidVerticesRestRest{s,2} = unique(ResidVerticesRestRest{s,2});
            ResidVerticesTaskTask{s,2} = unique(ResidVerticesTaskTask{s,2});
            ResidVerticesRestTask{s,2} = unique(ResidVerticesRestTask{s,2});
            
            [~,intersectrestrestodd,~] = intersect(ResidVerticesRestRest{s,1}, MatchedVerticesRestRest{s,1});
            [~,intersecttasktaskodd,~] = intersect(ResidVerticesTaskTask{s,1}, MatchedVerticesTaskTask{s,1});
            [~,intersectresttaskodd,~] = intersect(ResidVerticesRestTask{s,1}, MatchedVerticesRestTask{s,1});
            
            [~,intersectrestresteven,~] = intersect(ResidVerticesRestRest{s,2}, MatchedVerticesRestRest{s,2});
            [~,intersecttasktaskeven,~] = intersect(ResidVerticesTaskTask{s,2}, MatchedVerticesTaskTask{s,2});
            [~,intersectresttaskeven,~] = intersect(ResidVerticesRestTask{s,2}, MatchedVerticesRestTask{s,2});
            
            ResidVerticesRestRest{s,1}(:,intersectrestrestodd) = [];
            ResidVerticesTaskTask{s,1}(:,intersecttasktaskodd) = [];
            ResidVerticesRestTask{s,1}(:,intersectresttaskodd) = [];
            
            ResidVerticesRestRest{s,2}(:,intersectrestresteven) = [];
            ResidVerticesTaskTask{s,2}(:,intersecttasktaskeven) = [];
            ResidVerticesRestTask{s,2}(:,intersectresttaskeven) = [];
            
        end
    
        if SplitHalf == 1
        
            cifti_rest_even = ft_read_cifti_mod(rest_files_even{s});
            cifti_task_even = ft_read_cifti_mod(task_files_even{s});
            cifti_rest_odd = ft_read_cifti_mod(rest_files_odd{s});
            cifti_task_odd = ft_read_cifti_mod(task_files_odd{s});

        else
        
            cifti_rest = ft_read_cifti_mod(rest_files{s});
            cifti_task = ft_read_cifti_mod(task_files{s});
        
        end
    
    
        if SNRExclude == 1
            
            SNRmask = ft_read_cifti_mod(['/Users/briankraus/Desktop/MSC_SNR_Maps/CIFTI_Files/MSC_Template/Individual_Session_Map/' subject '/' subject '__SNRMap__AllDataConcatenated.dscalar.nii']); 
            
            SNRmask.data = SNRmask.data(1:59412,:);
            
            LowSNR = find(SNRmask.data < 750);
        
            if SplitHalf == 1
            
                cifti_rest_even.data(LowSNR,:) = NaN;
                cifti_task_even.data(LowSNR,:) = NaN;
                cifti_rest_odd.data(LowSNR,:) = NaN;
                cifti_task_odd.data(LowSNR,:) = NaN;
            
            else
            
                cifti_rest(LowSNR,:) = NaN;
                cifti_task(LowSNR,:) = NaN;
            
            end
        
        end
        
        if SplitHalf == 1
            
            if AbsoluteThresholds == 1
                
                cifti_rest_even_thresh = deciles(g,:);
                cifti_task_even_thresh = deciles(g,:);
                cifti_rest_odd_thresh = deciles(g,:);
                cifti_task_odd_thresh = deciles(g,:);
                
            else
            
                cifti_rest_even_thresh = prctile(cifti_rest_even.data,deciles(g,:));
                cifti_task_even_thresh = prctile(cifti_task_even.data,deciles(g,:));
                cifti_rest_odd_thresh = prctile(cifti_rest_odd.data,deciles(g,:));
                cifti_task_odd_thresh = prctile(cifti_task_odd.data,deciles(g,:));
                
                cifti_rest_even_valsSub = [cifti_rest_even_valsSub; cifti_rest_even_thresh(2)];
                cifti_task_even_valsSub = [cifti_task_even_valsSub; cifti_task_even_thresh(2)];
                cifti_rest_odd_valsSub = [cifti_rest_odd_valsSub; cifti_rest_odd_thresh(2)];
                cifti_task_odd_valsSub = [cifti_task_odd_valsSub; cifti_task_odd_thresh(2)];
                
            end
            
%             if g == 1
%                 
%                 cifti_rest_even_remove = find(cifti_rest_even.data >= cifti_rest_even_thresh(2));
%                 cifti_task_even_remove = find(cifti_task_even.data >= cifti_task_even_thresh(2));
%                 cifti_rest_odd_remove = find(cifti_rest_odd.data >= cifti_rest_odd_thresh(2));
%                 cifti_task_odd_remove = find(cifti_task_odd.data >= cifti_task_odd_thresh(2));
                 
            %elseif g > 1 && g < size(deciles,1)
            
            %if g < size(deciles,1)
            
                cifti_rest_even_remove = find(cifti_rest_even.data < cifti_rest_even_thresh(1) | cifti_rest_even.data >= cifti_rest_even_thresh(2));
                cifti_task_even_remove = find(cifti_task_even.data < cifti_task_even_thresh(1) | cifti_task_even.data >= cifti_task_even_thresh(2));
                cifti_rest_odd_remove = find(cifti_rest_odd.data < cifti_rest_odd_thresh(1) | cifti_rest_odd.data >= cifti_rest_odd_thresh(2));
                cifti_task_odd_remove = find(cifti_task_odd.data < cifti_task_odd_thresh(1) | cifti_task_odd.data >= cifti_task_odd_thresh(2));
            
%             elseif g == size(deciles,1)
%                 
%                 cifti_rest_even_remove = find(cifti_rest_even.data < cifti_rest_even_thresh(1));
%                 cifti_task_even_remove = find(cifti_task_even.data < cifti_task_even_thresh(1));
%                 cifti_rest_odd_remove = find(cifti_rest_odd.data < cifti_rest_odd_thresh(1));
%                 cifti_task_odd_remove = find(cifti_task_odd.data < cifti_task_odd_thresh(1));
%                 
            %end
            
            cifti_rest_even.data(cifti_rest_even_remove,:) = 0;
          	cifti_task_even.data(cifti_task_even_remove,:) = 0;
           	cifti_rest_odd.data(cifti_rest_odd_remove,:) = 0;
          	cifti_task_odd.data(cifti_task_odd_remove,:) = 0;
            
        else
            
           	cifti_rest_thresh = prctile(cifti_rest.data,deciles(g,:));
          	cifti_task_thresh = prctile(cifti_task.data,deciles(g,:));
            
            cifti_rest_remove = find(cifti_rest.data < cifti_rest_thresh(1) | cifti_rest.data >= cifti_rest_thresh(2));
            cifti_task_remove = find(cifti_task.data < cifti_task_thresh(1) | cifti_task.data >= cifti_task_thresh(2));
            
            cifti_rest.data(cifti_rest_remove,:) = 0;
            cifti_task.data(cifti_task_remove,:) = 0;
            
        end
        
        if deciles(g,2) == 10
            
            rest_plot_dat = cifti_rest_even.data;
            task_plot_dat = cifti_task_even.data;
            
            rest_plot_dat(cifti_rest_even_remove,:) = NaN;
          	task_plot_dat(cifti_task_even_remove,:) = NaN;
            
            %plotdat = {[rest_plot_dat; task_plot_dat]};
            
            plotdat = {rest_plot_dat};
            
            filename = ['/Subject_' subject '_RestHistogram_' num2str(deciles(g,2)) '_Decile.jpg'];
            
            plotitle = [subject ' Rest Histogram ' num2str(deciles(g,2)) ' Decile'];
            
            nhist(plotdat, 'binfactor', 2, 'samebins', 'separateplots');
            title(plotitle, 'fontsize',18)
            xlabel('Correlation Values (r)')
            ax = gca;
            ax.FontSize = 14;
            ylim([0 1000]);
    
            saveas(gcf,[outputdir filename])
    
            close gcf
            
        end

        if SplitHalf == 1
            
            if PercentOverlap == 1
            
                pctoverlapcounttaskrest = 0;
                pctoverlapcountrestrest = 0;
                pctoverlapcounttasktask = 0;
                
            end
            
            dcorrdatataskrest = [];
            dcorrdatarestrest = [];
            dcorrdatatasktask = [];
            
            if ResidualDeciles == 1     %% Create temporary data from left over vertices
                
                tempciftitaskptsodd = ResidVerticesRestTask{s,1};
                tempciftirestptseven = ResidVerticesRestTask{s,2};
                
                cifti_rest_even.data(tempciftirestptseven,:) = 1;
                cifti_task_odd.data(tempciftitaskptsodd,:) = 1;
                
            end
        
            for q = 1:length(cifti_rest_even.data)      %% Task-Rest Comparison
            
                if (isnan(cifti_rest_even.data(q)) == 0 && cifti_rest_even.data(q) ~= 0) && (isnan(cifti_task_odd.data(q)) == 0 && cifti_task_odd.data(q) ~= 0)
                
                    dcorrdatataskrest = [dcorrdatataskrest;1 1];
                    
                    if PercentOverlap == 1
                    
                        pctoverlapcounttaskrest = pctoverlapcounttaskrest + 1;
                        
                    end
                    
                    if ResidualDeciles == 1
                    
                        MatchedVerticesRestTask{s,1} = [MatchedVerticesRestTask{s,1} q];
                        MatchedVerticesRestTask{s,2} = [MatchedVerticesRestTask{s,2} q];
                        
                    end
                
                elseif (isnan(cifti_rest_even.data(q)) == 0 && cifti_rest_even.data(q) ~= 0)
                
                    dcorrdatataskrest = [dcorrdatataskrest;1 0];
                    
                    if ResidualDeciles == 1 && isnan(cifti_task_odd.data(q)) == 0
                    
                        ResidVerticesRestTask{s,2} = [ResidVerticesRestTask{s,2} q];
                        
                    end
                
                elseif (isnan(cifti_task_odd.data(q)) == 0 && cifti_task_odd.data(q) ~= 0)
                
                    dcorrdatataskrest = [dcorrdatataskrest;0 1];
                    
                    if ResidualDeciles == 1 && isnan(cifti_rest_even.data(q)) == 0
                    
                        ResidVerticesRestTask{s,1} = [ResidVerticesRestTask{s,1} q];
                        
                    end
                
                end
            
            end
            
            if ResidualDeciles == 1     %% Reset data and create new temporary data
                
                cifti_rest_even.data(tempciftirestptseven,:) = 0;
                cifti_task_odd.data(tempciftitaskptsodd,:) = 0;
                
                tempciftirestptsodd = ResidVerticesRestRest{s,1};
                tempciftirestptseven = ResidVerticesRestRest{s,2};
                
                cifti_rest_even.data(tempciftirestptseven,:) = 1;
                cifti_rest_odd.data(tempciftirestptsodd,:) = 1;
                
            end
            
            
            for q = 1:length(cifti_rest_odd.data)       %% Rest-Rest Comparison
            
                if (isnan(cifti_rest_odd.data(q)) == 0 && cifti_rest_odd.data(q) ~= 0) && (isnan(cifti_rest_even.data(q)) == 0 && cifti_rest_even.data(q) ~= 0)
                
                    dcorrdatarestrest = [dcorrdatarestrest;1 1];
                    
                    if PercentOverlap == 1
                    
                        pctoverlapcountrestrest = pctoverlapcountrestrest + 1;
                        
                    end
                    
                    if ResidualDeciles == 1
                    
                        MatchedVerticesRestRest{s,1} = [MatchedVerticesRestRest{s,1} q];
                        MatchedVerticesRestRest{s,2} = [MatchedVerticesRestRest{s,2} q];
                        
                    end
                
                elseif (isnan(cifti_rest_odd.data(q)) == 0 && cifti_rest_odd.data(q) ~= 0)
                
                    dcorrdatarestrest = [dcorrdatarestrest;1 0];
                    
                    if ResidualDeciles == 1 && isnan(cifti_rest_even.data(q)) == 0
                    
                        ResidVerticesRestRest{s,1} = [ResidVerticesRestRest{s,1} q];
                        
                    end
                
                elseif (isnan(cifti_rest_even.data(q)) == 0 && cifti_rest_even.data(q) ~= 0)
                
                    dcorrdatarestrest = [dcorrdatarestrest;0 1];
                    
                    if ResidualDeciles == 1 && isnan(cifti_rest_odd.data(q)) == 0
                    
                        ResidVerticesRestRest{s,2} = [ResidVerticesRestRest{s,2} q];
                        
                    end
                
                end
            
            end
            
            if ResidualDeciles == 1     %% Reset data and create new temporary data
                
                cifti_rest_even.data(tempciftirestptseven,:) = 0;
                cifti_rest_odd.data(tempciftirestptsodd,:) = 0;
                
                tempciftitaskptsodd = ResidVerticesTaskTask{s,1};
                tempciftitaskptseven = ResidVerticesTaskTask{s,2};
                
                cifti_task_even.data(tempciftitaskptseven,:) = 1;
                cifti_task_odd.data(tempciftitaskptsodd,:) = 1;
                
            end
            
            
            for q = 1:length(cifti_task_odd.data)       %% Task-Task Comparison
            
                if (isnan(cifti_task_odd.data(q)) == 0 && cifti_task_odd.data(q) ~= 0) && (isnan(cifti_task_even.data(q)) == 0 && cifti_task_even.data(q) ~= 0)
                
                    dcorrdatatasktask = [dcorrdatatasktask;1 1];
                    
                    if PercentOverlap == 1
                    
                        pctoverlapcounttasktask = pctoverlapcounttasktask + 1;
                        
                    end
                    
                    if ResidualDeciles == 1
                        
                        MatchedVerticesTaskTask{s,1} = [MatchedVerticesTaskTask{s,1} q];
                        MatchedVerticesTaskTask{s,2} = [MatchedVerticesTaskTask{s,2} q];
                        
                    end
                
                elseif (isnan(cifti_task_odd.data(q)) == 0 && cifti_task_odd.data(q) ~= 0)
                
                    dcorrdatatasktask = [dcorrdatatasktask;1 0];
                    
                    if ResidualDeciles == 1 && isnan(cifti_task_even.data(q)) == 0
                        
                        ResidVerticesTaskTask{s,1} = [ResidVerticesTaskTask{s,1} q];
                        
                    end
                
                elseif (isnan(cifti_task_even.data(q)) == 0 && cifti_task_even.data(q) ~= 0)
                
                    dcorrdatatasktask = [dcorrdatatasktask;0 1];
                    
                    if ResidualDeciles == 1 && isnan(cifti_task_odd.data(q)) == 0
                        
                        ResidVerticesTaskTask{s,2} = [ResidVerticesTaskTask{s,2} q];
                        
                    end
                
                end
            
            end
            
            if ResidualDeciles == 1     %% Reset data
                
                cifti_task_even.data(tempciftitaskptseven,:) = 0;
                cifti_task_odd.data(tempciftitaskptsodd,:) = 0;
                
            end 
                
            if isempty(dcorrdatataskrest) && isempty(dcorrdatarestrest) && isempty(dcorrdatatasktask)
                    
                dctaskrest = 0;
                dcrestrest = 0;
                dctasktask = 0;
                    
            elseif isempty(dcorrdatataskrest) && isempty(dcorrdatarestrest)
                    
                dctaskrest = 0;
                dcrestrest = 0;
                dctasktask = dice_coefficient_mod(dcorrdatatasktask(:,1),dcorrdatatasktask(:,2));
                    
            elseif isempty(dcorrdatataskrest) && isempty(dcorrdatatasktask)
                    
                dctaskrest = 0;
                dctasktask = 0;
                dcrestrest = dice_coefficient_mod(dcorrdatarestrest(:,1),dcorrdatarestrest(:,2));
                    
            elseif isempty(dcorrdatarestrest) && isempty(dcorrdatatasktask)
                    
                dcrestrest = 0;
                dctasktask = 0;
                dctaskrest = dice_coefficient_mod(dcorrdatataskrest(:,1),dcorrdatataskrest(:,2));
                    
            elseif isempty(dcorrdatarestrest)
                    
                dcrestrest = 0;
                dctaskrest = dice_coefficient_mod(dcorrdatataskrest(:,1),dcorrdatataskrest(:,2));
                dctasktask = dice_coefficient_mod(dcorrdatatasktask(:,1),dcorrdatatasktask(:,2));
                    
            elseif isempty(dcorrdatatasktask)
                    
                dctasktask = 0;
                dctaskrest = dice_coefficient_mod(dcorrdatataskrest(:,1),dcorrdatataskrest(:,2));
                dcrestrest = dice_coefficient_mod(dcorrdatarestrest(:,1),dcorrdatarestrest(:,2));
                    
            elseif isempty(dcorrdatataskrest)
                    
                dctaskrest = 0;
                dcrestrest = dice_coefficient_mod(dcorrdatarestrest(:,1),dcorrdatarestrest(:,2));
                dctasktask = dice_coefficient_mod(dcorrdatatasktask(:,1),dcorrdatatasktask(:,2));
                    
            else

                dctaskrest = dice_coefficient_mod(dcorrdatataskrest(:,1),dcorrdatataskrest(:,2));
                dcrestrest = dice_coefficient_mod(dcorrdatarestrest(:,1),dcorrdatarestrest(:,2));
                dctasktask = dice_coefficient_mod(dcorrdatatasktask(:,1),dcorrdatatasktask(:,2));
                    
            end
            
            if PercentOverlap == 1
            
                pctoverlapcounttaskrestTemp = [pctoverlapcounttaskrestTemp; pctoverlapcounttaskrest/size(dcorrdatataskrest,1)];
                pctoverlapcountrestrestTemp = [pctoverlapcountrestrestTemp; pctoverlapcountrestrest/size(dcorrdatarestrest,1)];
                pctoverlapcounttasktaskTemp = [pctoverlapcounttasktaskTemp; pctoverlapcounttasktask/size(dcorrdatatasktask,1)];
            
            end
        
            DiceCorrsTaskRestTemp = [DiceCorrsTaskRestTemp; dctaskrest];
            DiceCorrsTaskTaskTemp = [DiceCorrsTaskTaskTemp; dctasktask];
            DiceCorrsRestRestTemp = [DiceCorrsRestRestTemp; dcrestrest];
            
        else
        
            dcorrdata = [];
        
            for q = 1:length(cifti_rest.data)
            
                if cifti_rest.data(q) > 0 && cifti_task.data(q) > 0
                
                    dcorrdata = [dcorrdata;1 1];
                
                elseif cifti_rest.data(q) > 0
                
                    dcorrdata = [dcorrdata;1 0];
                
                elseif cifti_task.data(q) > 0
                
                    dcorrdata = [dcorrdata;0 1];
                
                end
            
            end
        
            dc = dice_coefficient_mod(dcorrdata(:,1),dcorrdata(:,2));
        
            DiceCorrsTemp = [DiceCorrsTemp;dc];
            
        end
        
        if SplitHalf == 1
                
        	alltaskfilestaskeventemp = [alltaskfilestaskeventemp; cifti_task_even.data'];
          	alltaskfilestaskoddtemp = [alltaskfilestaskoddtemp; cifti_task_odd.data'];
        	alltaskfilesresteventemp = [alltaskfilesresteventemp; cifti_rest_even.data'];
          	alltaskfilesrestoddtemp = [alltaskfilesrestoddtemp; cifti_rest_odd.data'];
                
        else
                        
          	alltaskfilestasktemp = [alltaskfilestasktemp; cifti_task.data'];
          	alltaskfilesresttemp = [alltaskfilesresttemp; cifti_rest.data']; 
            
        end
        
    end
    
    if SplitHalf == 1
        
        if PercentOverlap == 1
        
            pctoverlaptaskrest = [pctoverlaptaskrest pctoverlapcounttaskrestTemp];
            pctoverlaprestrest = [pctoverlaprestrest pctoverlapcountrestrestTemp];
            pctoverlaptasktask = [pctoverlaptasktask pctoverlapcounttasktaskTemp];
            
        end
        
        DiceCorrsTaskRest = [DiceCorrsTaskRest DiceCorrsTaskRestTemp];
        DiceCorrsTaskTask = [DiceCorrsTaskTask DiceCorrsTaskTaskTemp];
        DiceCorrsRestRest = [DiceCorrsRestRest DiceCorrsRestRestTemp];
        
    else
    
        DiceCorrs = [DiceCorrs DiceCorrsTemp];
        
    end
    
	if SplitHalf == 1
                
     	alltaskfilestaskeven = [alltaskfilestaskeven alltaskfilestaskeventemp];
     	alltaskfilestaskodd = [alltaskfilestaskodd alltaskfilestaskoddtemp];
     	alltaskfilesresteven = [alltaskfilesresteven alltaskfilesresteventemp];
      	alltaskfilesrestodd = [alltaskfilesrestodd alltaskfilesrestoddtemp];
                
    else
                        
      	alltaskfilestask = [alltaskfilestask alltaskfilestasktemp];
     	alltaskfilesrest = [alltaskfilesrest alltaskfilesresttemp]; 
            
    end
    
    
    if AbsoluteThresholds ~= 1
    
        cifti_rest_even_vals = [cifti_rest_even_vals cifti_rest_even_valsSub];
        cifti_task_even_vals = [cifti_task_even_vals cifti_task_even_valsSub];
        cifti_rest_odd_vals = [cifti_rest_odd_vals cifti_rest_odd_valsSub];
        cifti_task_odd_vals = [cifti_task_odd_vals cifti_task_odd_valsSub];
    
    end
    
end


if ResidualDeciles == 1
    
    SubjectResidVerticesRest = cell(length(subs)-1,length(subs));
    SubjectResidVerticesTask = cell(length(subs)-1,length(subs));
    AllResidVerticesRest = cell(length(subs)-1,length(subs));
    AllResidVerticesTask = cell(length(subs)-1,length(subs));
    
    SubjectMatchedVerticesRest = cell(length(subs)-1,length(subs));
    SubjectMatchedVerticesTask = cell(length(subs)-1,length(subs));
    AllMatchedVerticesRest = cell(length(subs)-1,length(subs));
    AllMatchedVerticesTask = cell(length(subs)-1,length(subs));
    
end

SubjectOverlapFinal = [];
SubjectDiceCorrsFinal = [];

for a = 1:size(deciles,1)
    
    if PercentOverlap == 1

        SubjectOverlap = [];
        
    end
    
    SubjectDiceCorrs = [];
    
    for j = 1:length(subs)
        
        subject = subs{j};
        
        if SplitHalf == 1
            
            alltasktasksubject = alltaskfilestaskodd{:,a}(j,:);
            alltaskrestsubject = alltaskfilesresteven{:,a}(j,:);
            
            alltasktaskcomp = alltaskfilestaskodd{:,a};
            alltasktaskcomp(j,:) = [];
            alltaskrestcomp = alltaskfilesresteven{:,a};
            alltaskrestcomp(j,:) = [];
            
        else
        
            alltasktasksubject = alltaskfilestask{:,a}(j,:);
            alltaskrestsubject = alltaskfilesrest{:,a}(j,:);
            
            alltasktaskcomp = alltaskfilestask{:,a};
            alltasktaskcomp(j,:) = [];
            alltaskrestcomp = alltaskfilesrest{:,a};
            alltaskrestcomp(j,:) = [];
            
        end
        
        for d = 1:size(alltasktaskcomp,1)
            
            if ResidualDeciles == 1 && (~isempty(SubjectResidVerticesRest{d,j}) || ~isempty(SubjectResidVerticesTask{d,j}))
            
            	SubjectResidVerticesRest{d,j} = unique(SubjectResidVerticesRest{d,j});
                SubjectResidVerticesTask{d,j} = unique(SubjectResidVerticesTask{d,j});
            	AllResidVerticesRest{d,j} = unique(AllResidVerticesRest{d,j});
                AllResidVerticesTask{d,j} = unique(AllResidVerticesTask{d,j});
                
                [~,intersectrestsub,~] = intersect(SubjectResidVerticesRest{d,j}, SubjectMatchedVerticesRest{d,j});
                [~,intersecttasksub,~] = intersect(SubjectResidVerticesTask{d,j}, SubjectMatchedVerticesTask{d,j});
                [~,intersectrestall,~] = intersect(AllResidVerticesRest{d,j}, AllMatchedVerticesRest{d,j});
                [~,intersecttaskall,~] = intersect(AllResidVerticesTask{d,j}, AllMatchedVerticesTask{d,j});

            
            	SubjectResidVerticesRest{d,j}(:,intersectrestsub) = [];
                SubjectResidVerticesTask{d,j}(:,intersecttasksub) = [];
            	AllResidVerticesRest{d,j}(:,intersectrestall) = [];
                AllResidVerticesTask{d,j}(:,intersecttaskall) = [];

            end
            
            if PercentOverlap == 1

                SubjectOverlapTemp = [];
                
            end
            
            SubjectDiceCorrstemp = [];
            
            if PercentOverlap == 1
            
                overlapcounttask = 0;
                overlapcountrest = 0;
                
            end
            
          	dcorrdatataskcomp = [];
          	dcorrdatarestcomp = [];
            
            if ResidualDeciles == 1     %% Create temporary data from left over vertices
                
                subtemptaskpts = SubjectResidVerticesTask{d,j};
                alltemptaskpts = AllResidVerticesTask{d,j};
                
                alltaskrestsubject(:,subtemptaskpts) = 1;
                alltasktaskcomp(:,alltemptaskpts) = 1;
                
            end
        
           	for q = 1:length(alltasktasksubject)
            
             	if (isnan(alltaskrestsubject(:,q)) == 0 && alltaskrestsubject(:,q) ~= 0) && (isnan(alltasktaskcomp(d,q)) == 0 && alltasktaskcomp(d,q) ~= 0)
                
                  	dcorrdatataskcomp = [dcorrdatataskcomp;1 1];
                    
                    if PercentOverlap == 1
                    
                        overlapcounttask = overlapcounttask + 1;
                        
                    end
                    
                    if ResidualDeciles == 1
                    
                        SubjectMatchedVerticesTask{d,j} = [SubjectMatchedVerticesTask{d,j} q];
                        AllMatchedVerticesTask{d,j} = [AllMatchedVerticesTask{d,j} q];
                        
                    end
                
              	elseif (isnan(alltaskrestsubject(:,q)) == 0 && alltaskrestsubject(:,q) ~= 0)
                
                  	dcorrdatataskcomp = [dcorrdatataskcomp;1 0];
                    
                    if ResidualDeciles == 1 && isnan(alltaskrestsubject(:,q)) == 0
                    
                        SubjectResidVerticesTask{d,j} = [SubjectResidVerticesTask{d,j} q];
                        
                    end
                
               	elseif (isnan(alltasktaskcomp(d,q)) == 0 && alltasktaskcomp(d,q) ~= 0)
                
                   	dcorrdatataskcomp = [dcorrdatataskcomp;0 1];
                    
                   if ResidualDeciles == 1 && isnan(alltasktaskcomp(d,q)) == 0
                    
                        AllResidVerticesTask{d,j} = [AllResidVerticesTask{d,j} q];
                        
                   end
                
                end
            
            end
            
            if ResidualDeciles == 1     %% Reset data and create new temporary data
                
                alltaskrestsubject(:,subtemptaskpts) = 0;
                alltasktaskcomp(:,alltemptaskpts) = 0;
                
                subtemprestpts = SubjectResidVerticesRest{d,j};
                alltemprestpts = AllResidVerticesRest{d,j};
                
                alltasktasksubject(:,subtemprestpts) = 1;
                alltaskrestcomp(:,alltemprestpts) = 1;
                
            end
                
         	for r = 1:length(alltaskrestsubject)
            
            	if (isnan(alltasktasksubject(:,r)) == 0 && alltasktasksubject(:,r) ~= 0) && (isnan(alltaskrestcomp(d,r)) == 0 && alltaskrestcomp(d,r) ~= 0)
                
                 	dcorrdatarestcomp = [dcorrdatarestcomp;1 1];
                    
                    if PercentOverlap == 1
                    
                        overlapcountrest = overlapcountrest + 1;
                        
                    end
                    
                    if ResidualDeciles == 1
                        
                        SubjectMatchedVerticesRest{d,j} = [SubjectMatchedVerticesRest{d,j} r];
                        AllMatchedVerticesRest{d,j} = [AllMatchedVerticesRest{d,j} r];
                        
                    end
                
               	elseif (isnan(alltasktasksubject(:,r)) == 0 && alltasktasksubject(:,r) ~= 0)
                
                   	dcorrdatarestcomp = [dcorrdatarestcomp;1 0];
                    
                    if ResidualDeciles == 1 && isnan(alltasktasksubject(:,r)) == 0
                        
                        SubjectResidVerticesRest{d,j} = [SubjectResidVerticesRest{d,j} r];
                        
                    end
                
              	elseif (isnan(alltaskrestcomp(d,r)) == 0 && alltaskrestcomp(d,r) ~= 0)
                
                 	dcorrdatarestcomp = [dcorrdatarestcomp;0 1];
                    
                    if ResidualDeciles == 1 && isnan(alltaskrestcomp(d,r)) == 0
                        
                        AllResidVerticesRest{d,j} = [AllResidVerticesRest{d,j} r];
                        
                    end
                
                end
            
            end
            
            if ResidualDeciles == 1     %% Reset data
                
                alltasktasksubject(:,subtemprestpts) = 0;
                alltaskrestcomp(:,alltemprestpts) = 0;
                
            end 
                    
           	if isempty(dcorrdatataskcomp) && isempty(dcorrdatarestcomp)
                        
               	dctaskcomp = 0;
               	dcrestcomp = 0;
                        
           	elseif isempty(dcorrdatataskcomp)
                        
               	dctaskcomp = 0;
               	dcrestcomp = dice_coefficient_mod(dcorrdatarestcomp(:,1),dcorrdatarestcomp(:,2));
                        
           	elseif isempty(dcorrdatarestcomp)
                        
               	dcrestcomp = 0;
              	dctaskcomp = dice_coefficient_mod(dcorrdatataskcomp(:,1),dcorrdatataskcomp(:,2));
                        
            else
                        
               	dctaskcomp = dice_coefficient_mod(dcorrdatataskcomp(:,1),dcorrdatataskcomp(:,2));
              	dcrestcomp = dice_coefficient_mod(dcorrdatarestcomp(:,1),dcorrdatarestcomp(:,2));
                        
            end
        
          	SubjectDiceCorrstemp = [SubjectDiceCorrstemp dctaskcomp dcrestcomp];
            
            if PercentOverlap == 1
            
                SubjectOverlapTemp = [SubjectOverlapTemp overlapcounttask/size(dcorrdatataskcomp,1) overlapcountrest/size(dcorrdatarestcomp,1)];
            
            end
        end
        
        SubjectDiceCorrs = [SubjectDiceCorrs; SubjectDiceCorrstemp];
        
        if PercentOverlap == 1
            
            SubjectOverlap = [SubjectOverlap; SubjectOverlapTemp];
            
        end
    end
    
    SubjectDiceCorrsFinal = [SubjectDiceCorrsFinal mean(SubjectDiceCorrs,2)];
    
    if PercentOverlap == 1
    
        SubjectOverlapFinal = [SubjectOverlapFinal mean(SubjectOverlap,2)];
    
    end
end





if AbsoluteThresholds == 1
    
    %plotthresholds = [.1:.1:.8];
    
    %plotthresholds = [.1:.05:.75];
    plotthresholds = [.1:.025:.725];
    %plotthresholds = [.3:.3:.9];
    
else

    %plotthresholds = [10:10:100];
    plotthresholds = [2.5:2.5:100];
    
    if numel(plotthresholds) == 10
    
        plotticks = {[num2str(plotthresholds(1)) '%' ' ' num2str(mean(cifti_rest_even_vals(:,1),1)) '(r)']
                    [num2str(plotthresholds(2)) '%' ' ' num2str(mean(cifti_rest_even_vals(:,2),1)) '(r)']
                    [num2str(plotthresholds(3)) '%' ' ' num2str(mean(cifti_rest_even_vals(:,3),1)) '(r)']
                    [num2str(plotthresholds(4)) '%' ' ' num2str(mean(cifti_rest_even_vals(:,4),1)) '(r)']
                    [num2str(plotthresholds(5)) '%' ' ' num2str(mean(cifti_rest_even_vals(:,5),1)) '(r)']
                    [num2str(plotthresholds(6)) '%' ' ' num2str(mean(cifti_rest_even_vals(:,6),1)) '(r)']
                    [num2str(plotthresholds(7)) '%' ' ' num2str(mean(cifti_rest_even_vals(:,7),1)) '(r)']
                    [num2str(plotthresholds(8)) '%' ' ' num2str(mean(cifti_rest_even_vals(:,8),1)) '(r)']
                    [num2str(plotthresholds(9)) '%' ' ' num2str(mean(cifti_rest_even_vals(:,9),1)) '(r)']
                    [num2str(plotthresholds(10)) '%' ' ' num2str(mean(cifti_rest_even_vals(:,10),1)) '(r)']};
                
    end
            

end

if PercentOverlap == 1

    temp = [DiceCorrsRestRest; DiceCorrsTaskTask];
    tempa = [pctoverlaprestrest; pctoverlaptasktask];
    temp2 = DiceCorrsTaskRest;
    temp2a = pctoverlaptaskrest;
    temp3 = SubjectDiceCorrsFinal;
    temp3a = SubjectOverlapFinal;

    shadedErrorBar(plotthresholds,mean(temp,1),std(temp,1)/sqrt(size(temp,1)),'lineprops','g');
    shadedErrorBar(plotthresholds,mean(tempa,1),std(tempa,1)/sqrt(size(tempa,1)),'lineprops','g--');
    shadedErrorBar(plotthresholds,mean(temp2,1),std(temp2,1)/sqrt(size(temp2,1)),'lineprops','b');
    shadedErrorBar(plotthresholds,mean(temp2a,1),std(temp2a,1)/sqrt(size(temp2a,1)),'lineprops','b--');
    shadedErrorBar(plotthresholds,mean(temp3,1),std(temp3,1)/sqrt(size(temp3,1)),'lineprops','k');
    shadedErrorBar(plotthresholds,mean(temp3a,1),std(temp3a,1)/sqrt(size(temp3a,1)),'lineprops','k--');
    c = get(gca, 'Children');
    hleg1 = legend(c(1:6), 'Across Subjects Percent Overlap', 'Across Subjects Dice', 'Across State Percent Overlap', 'Across State Dice', 'Within State Percent Overlap', 'Within State Dice', 'Location', 'North');
    hleg1.FontSize = 14;
    %title('Decile Overlap Plot', 'fontsize',18)
    ylabel('Dice Correlations')
    xlabel('Similarity to Group Functional Connectivity')
    ax = gca;
    ax.FontSize = 14;
    ylim([0 .75])

    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.3, 0.3, 0.3, 0.3]);
    
else
    
    temp = [DiceCorrsRestRest; DiceCorrsTaskTask];
    temp2 = DiceCorrsTaskRest;
    temp3 = SubjectDiceCorrsFinal;

    shadedErrorBar(plotthresholds,mean(temp,1),std(temp,1)/sqrt(size(temp,1)),'lineprops','g');
    shadedErrorBar(plotthresholds,mean(temp2,1),std(temp2,1)/sqrt(size(temp2,1)),'lineprops','b');
    shadedErrorBar(plotthresholds,mean(temp3,1),std(temp3,1)/sqrt(size(temp3,1)),'lineprops','k');
    c = get(gca, 'Children');
    hleg1 = legend(c(1:3), 'Across Subjects Overlap', 'Across State Overlap', 'Within State Overlap', 'Location', 'North');
    hleg1.FontSize = 14;
    %title('Decile Overlap Plot', 'fontsize',18)
    ylabel('Dice Correlations')
    xlabel('Similarity to Group Functional Connectivity')
    ax = gca;
    ax.FontSize = 14;
    ylim([0 .75])

    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.3, 0.3, 0.3, 0.3]);
    
end

if AbsoluteThresholds == 1
    
elseif ResidualDeciles == 0 && numel(plotthresholds) == 10

    set(gca,'xticklabel',plotticks')

end

if AbsoluteThresholds == 1 && ResidualDeciles == 1
    
    print(gcf,[outputdir '/AllSubjects_Decile_Overlap_Absolute_Resids_Plot.jpg'],'-dpng','-r300');
    
	%saveas(gcf,[outputdir '/AllSubjects_Decile_Overlap_Absolute_Resids_Plot.jpg'])
    
elseif AbsoluteThresholds == 1
    
    print(gcf,[outputdir '/AllSubjects_Decile_Overlap_Absolute_Plot.jpg'],'-dpng','-r300');
    
   %saveas(gcf,[outputdir '/AllSubjects_Decile_Overlap_Absolute_Plot.jpg'])
    
elseif ResidualDeciles == 1
    
    print(gcf,[outputdir '/AllSubjects_Decile_Overlap_Resids_Plot.jpg'],'-dpng','-r300');
    
	%saveas(gcf,[outputdir '/AllSubjects_Decile_Overlap_Resids_Plot.jpg'])
    
elseif PercentOverlap == 1
    
    print(gcf,[outputdir '/AllSubjects_Decile_Overlap_Plot.jpg'],'-dpng','-r300');

    %saveas(gcf,[outputdir '/AllSubjects_Decile_Overlap_Plot.jpg'])
    
else
    
    print(gcf,[outputdir '/AllSubjects_Decile_Overlap_Dice_Only_Plot.jpg'],'-dpng','-r300');

end

close gcf



