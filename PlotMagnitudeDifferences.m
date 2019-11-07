
%%% Plots absolute differences in magnitude between task and rest variants
%%% for overlapping/non-overlapping vertices

clear all

outputdir = '/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Magnitude_Histograms';
SNRExclude = 1;
SizeExclude = 1;
threshold = 2.5;



if SNRExclude == 1 && SizeExclude == 1 && threshold == 10
    
	[magnitude_files, subjects1, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_magnitude_varmaps_matched.txt','%s%s%s');

	[overlap_files, subjects2, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_overlap_varmaps_matched.txt','%s%s%s');
    
elseif SNRExclude == 1 && SizeExclude == 1 && threshold == 2.5
    
	[magnitude_files, subjects1, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_magnitude_varmaps_matched_2.5.txt','%s%s%s');

	[overlap_files, subjects2, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_overlap_varmaps_matched_2.5.txt','%s%s%s');
    
end

subs = unique(subjects1);

OverlapMagnitudeAllVerts = [];
OverlapMagnitudeAll = [];
RestMagnitudeAll = [];
TaskMagnitudeAll = [];

for s = 1:numel(subs)
    
	subject = subs{s};
	magnitude_file = magnitude_files{s};
	overlap_file = overlap_files{s};
        
	cifti_magnitude = ft_read_cifti_mod(magnitude_file);
	cifti_overlap = ft_read_cifti_mod(overlap_file);
    
    for z = 1:4     %% For task/rest/overlap portions of map
        
        if z == 1   %% Overlap
            
            magnitudepts = find(cifti_overlap.data == 1);
            
            filename = ['/' subject '_OverlapMagnitudeDiff_SNRSizeExclude_Histogram_' num2str(threshold) '_pct.jpg'];
            plotitle = ['Overlapping Variants Magnitude Differences ' subject];
            
        elseif z == 2   %% Rest
            
            magnitudepts = find(round(cifti_overlap.data,1) == .7);
              
            filename = ['/' subject '_RestMagnitudeDiff_SNRSizeExclude_Histogram_' num2str(threshold) '_pct.jpg'];
            plotitle = ['Rest Only Variants Magnitude Differences ' subject];
            
        elseif z == 3   %% Task
            
            magnitudepts = find(round(cifti_overlap.data,1) == .5);
            
            filename = ['/' subject '_TaskMagnitudeDiff_SNRSizeExclude_Histogram_' num2str(threshold) '_pct.jpg'];
            plotitle = ['Task Only Variants Magnitude Differences ' subject];
            
        elseif z == 4   %% All
            
            magnitudepts = find(round(cifti_overlap.data,1) > 0);
            
            filename = ['/' subject '_AllMagnitudeDiff_SNRSizeExclude_Histogram_' num2str(threshold) '_pct.jpg'];
            plotitle = ['All Variants Magnitude Differences ' subject];
            
        end
        
        magnitudedata = cifti_magnitude.data(magnitudepts,:);
        
        if z == 1   %% Overlap Only
            
            OverlapMagnitudeAll = [OverlapMagnitudeAll; magnitudedata];
            
        elseif z == 2   %% Rest
            
            RestMagnitudeAll = [RestMagnitudeAll; magnitudedata];
            
        elseif z == 3   %% Task
            
            TaskMagnitudeAll = [TaskMagnitudeAll; magnitudedata];
            
        elseif z == 4   %% All vertices
            
            OverlapMagnitudeAllVerts = [OverlapMagnitudeAllVerts; magnitudedata];
            
        end
        
        nhist(magnitudedata, 'binfactor', 2, 'samebins');
        title(plotitle, 'fontsize',18)
        xlabel('Absolute Differences (r)')
        ax = gca;
        ax.FontSize = 14;
        ylim([0 350]);
    
        saveas(gcf,[outputdir filename])
    
        close gcf
        
    end
end

for z = 1:4     %% For task/rest/overlap portions of map
    
    if z == 1   %% Overlap Only

        filename = ['/AllSubjects_OverlapMagnitudeDiff_SNRSizeExclude_Histogram_' num2str(threshold) '_pct.jpg'];
        plotitle = 'Overlapping Variants Magnitude Differences All Subjects';
        magnitudedata = OverlapMagnitudeAll;
        
    elseif z == 2   %% Rest
        
        filename = ['/AllSubjects_RestMagnitudeDiff_SNRSizeExclude_Histogram_' num2str(threshold) '_pct.jpg'];
        plotitle = 'Rest Only Variants Magnitude Differences All Subjects';
        magnitudedata = RestMagnitudeAll;
        
    elseif z == 3   %% Task
        
        filename = ['/AllSubjects_TaskMagnitudeDiff_SNRSizeExclude_Histogram_' num2str(threshold) '_pct.jpg'];
        plotitle = 'Task Only Variants Magnitude Differences All Subjects';
        magnitudedata = TaskMagnitudeAll;
        
    elseif z == 4
        
        filename = ['/AllSubjects_AllMagnitudeDiff_SNRSizeExclude_Histogram_' num2str(threshold) '_pct.jpg'];
        plotitle = '';
        magnitudedata = OverlapMagnitudeAllVerts;
        
    end

nhist(magnitudedata, 'binfactor', 2, 'samebins', 'noerror', 'stdtimes', 7);
title(plotitle, 'fontsize',18)
xlabel('Absolute Differences (r)')
ylabel('Number of Vertices')
ax = gca;
ax.FontSize = 14;
ylim([0 600]);
    
%saveas(gcf,[outputdir filename])

print(gcf,[outputdir filename],'-dpng','-r300');
    
close gcf


end
        
        