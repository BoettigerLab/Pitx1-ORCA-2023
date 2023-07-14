% Code for "" Fig.1 
% Please run "Load_data_00.m" first.

% enter your directories here:
data_folder = 'Z:/TzuChiao/Manuscript/Codes/data/';
save_folder = 'U:\Manuscripts\TzuChiao Paper\NatureGenetics\Revision1\Revised_images/'; 

%% b: ChIPseq tracks
% other tracks were uploaded to UCSC as big-wig bigurl links, then
% downloaded as an EPS file using the image export, after matching the
% chromosomal range. 

track_names = {{'mm10GSM2251432_CTCF-FL-E105-Wt-Mm-merge.bw.wig','mm10GSM2251438_CTCF-HL-E105-Wt-Mm-merge.bw.wig'},...
    {'mm10GSM2251474_H3K27Ac-FL-E105-Wt-Mm-Merge.bw.wig','mm10GSM2251480_H3K27Ac-HL-E105-Wt-Mm-Merge.bw.wig'},...
    {'mm10GSM2251488_H3K27me3-FL-E105-Wt-Mm-Merge.bw.wig','mm10GSM2251494_H3K27me3-HL-E105-Wt-Mm-Merge.bw.wig'},...
    {'mm10GSM2251502_RAD21-FL-E105-Wt-Mm-merge-qfrags.bw.wig','mm10GSM2251508_RAD21-HL-E105-Wt-Mm-merge-qfrags.bw.wig'}};

x= 55615000:200:(56365000);
% Plotting ChIPseq tracks
ylims = [30, 10, 10, 80];    % This is determined ad-hoc for better visualization.

off_limit_marker_size = 1.075;
delta_size = 1/40;
fL_color = [0,86,128]/255;
fL_color_off_limit = [0,200,255]/255;
hL_color = [170,33,41]/255;
hL_color_off_limit = [255,0,255]/255;
line_w = 2;
for s = 1:length(track_names)
    ft = readtable([data_folder, track_names{1,s}{1}],'FileType', 'text','delimiter','\t','HeaderLines',0);
    ht = readtable([data_folder, track_names{1,s}{2}],'FileType', 'text','delimiter','\t','HeaderLines',0);
    
    fCt = interval_maximum(ft{:,2}, ft{:,3}, ft{:,4}, x);
    hCt = interval_maximum(ht{:,2}, ht{:,3}, ht{:,4}, x);
    fCt2 = interp1(ft{:,2}, ft{:,4}, x);
    hCt2 = interp1(ht{:,2}, ht{:,4}, x);
 
    f1 = figure(s); clf; 
    max_y = ylims(s);
    delta = max_y * delta_size; 
    % delta creates a gap at the middle. Otherwise hindlimb
    % tracks, plotted latter, would be laying on top of the forelimb
    % tracks, creating a biased impression that hindlimb track has more
    % signal.
    fy = fCt > max_y;
    fCt_cut = fCt;
    fCt_cut(fy) = max_y;
    hy = hCt > max_y;
    hCt_cut = hCt;
    hCt_cut(hy) = max_y;
    
    area(x,fCt_cut+delta,'FaceColor',fL_color,'EdgeColor',fL_color, 'LineWidth', line_w); hold on; 
    area(x,-hCt_cut,'FaceColor', hL_color,'EdgeColor', hL_color, 'LineWidth', line_w); hold on;
    xlim([x(1),x(end)]);
    set(gcf, 'Position', [10 500 1100 150], 'Color', 'white');
    set(gca, 'ylim', off_limit_marker_size*[-max_y, max_y+delta]);
    set(gca, 'ytick', [-max_y max_y+delta]);
    set(gca, 'yticklabel', {string(max_y),string(max_y)});
    pause(1);
    ax = gca;
    ax.YRuler.Axle.LineStyle='none';
    ax.XRuler.Axle.LineStyle='none';
    set(gca,'XColor','none');
    
    % Create white bands that can be further adjusted in Illustrator to
    % block off confusing graphics resulted from the thickened outline of
    % the tracks.
    rectangle('Position',[x(1), 0, x(end) - x(1), delta], 'FaceColor', 'white', 'EdgeColor','white', 'LineWidth', line_w);  hold on;
    rectangle('Position',[x(1), max_y, x(end) - x(1), delta], 'FaceColor', 'white', 'EdgeColor','white', 'LineWidth', line_w);    hold on;
    rectangle('Position',[x(1), -max_y-delta, x(end) - x(1), delta], 'FaceColor', 'white', 'EdgeColor','white', 'LineWidth', line_w);hold on;
    
    % Plot the out-of-limit signals.
    scatter(find(fy) * (200) + x(1), max_y * off_limit_marker_size + delta, 1, '|', "MarkerEdgeColor", fL_color_off_limit, 'LineWidth', line_w);   hold on;
    scatter(find(hy) * (200) + x(1), -max_y * off_limit_marker_size, 1, '|', "MarkerEdgeColor", hL_color_off_limit,  'LineWidth', line_w');

    %Extract ChIP target from file name;
    target = split(track_names{1,s}{1}, '_');
    target = split(target{2}, '-');
    target = target{1};
    
    exportgraphics(f1,strcat(save_folder,"20230703_chip",target,".eps"),...
    'BackgroundColor','white','ContentType','vector');

end

%% a: ORCA maps
[hLcontact,n1] = ContactFrac(hLmaps,'threshold',th);
hLcontact = InterpMapNans(hLcontact,'badHybes',badHybes,'badPixels',[12,24]);
[fLcontact,n2] = ContactFrac(fLmaps,'threshold',th);
fLcontact = InterpMapNans(fLcontact,'badHybes',badHybes,'badPixels',[12,24]);
hLcontact(hLcontact<0) = 0; fLcontact(fLcontact<0)=0;
rm = GetColorMap('redToWhiteK');

%
figure(1); clf; 
imagesc(fLcontact); 
set(gcf,'color','w'); 
colormap(flipud(rm));
axis square;  
clim([.15 .3]);
saveas(gcf, strcat(save_folder, "fig1_FL_ORCA.epsc"));

%
figure(2); clf; 
imagesc(hLcontact); 
set(gcf,'color','w'); 
colormap(flipud(rm));
axis square;  
clim([.15 .3]);
saveas(gcf, strcat(save_folder, "fig1_HL_ORCA.epsc"));

%%
f2 = figure(3); clf; 
imagesc(log2(hLcontact./fLcontact)); colorbar;
title('log_2( hindlimb / forelimb ) contact freq.');
set(gcf,'color','w'); 
axis image;  
% play with some colormaps
    
load cutoffBluewhitered2.mat B;
colormap(B);
clim([-.3 .3]);
saveas(gcf, strcat(save_folder, "fig1_HF_diff.epsc"));
%% d: correlation of Hi-C and ORCA.  
hL_chic = load([data_folder, 'table_mm10HL-E115-Wt-Mm-cHiC-merged.hicup.MAPQ30.KR_5kb.WashU_5kb.mat']);
hL_chic = hL_chic.mat;
fL_chic = load([data_folder, 'table_mm10FL-E115-Wt-Mm-cHiC-merged.hicup.MAPQ30.KR_5kb.WashU_5kb.mat']);
fL_chic = fL_chic.mat;
chic_FL = imresize(fL_chic(1:150,1:150), [75,75], "nearest");
chic_HL = imresize(hL_chic(1:150,1:150), [75,75], "nearest");
corrS = zeros(0,4); %chicF, chicH, contF, contH
for d = 1:75
    dVec = horzcat(diag(chic_FL, d), diag(chic_HL, d), diag(fLcontact, d), diag(hLcontact, d));
    corrS = vertcat(corrS, dVec);
end
f1 = figure(4); clf; 
subplot(1,2,1); PlotCorr(corrS(:,1),corrS(:,3),'MarkerSize',1, 'color', [0,86,128]/255); axis square; xlabel('Forelimb cHi-C (norm. reads)'); ylabel('Forelimb ORCA (frac of cells)'); 
subplot(1,2,2); PlotCorr(corrS(:,2),corrS(:,4),'MarkerSize',1, 'color', [170,33,41]/255); axis square; xlabel('Hindlimb cHi-C (norm. reads)'); ylabel('Hindlimb ORCA (frac of cells)'); 
saveas(gcf, strcat(save_folder, "fig1_chicORCA_corr.epsc"));
%% e: single chromatins
% Will need to rotate the polymer to find the angle in the paper.
shift_H = -11607;
shift_F = -10087;
hlexp = 20949 + shift_H;
flexp = 31332 + shift_F;
plotPolymerWrapper(hLpol(:,:,hlexp), 'neurog', true);
plotPolymerWrapper(fLpol(:,:,flexp), 'neurog', true);
% Rotate before saving
%% saving in vector format
exportgraphics(gcf,strcat(save_folder,"1_05_FLpol.eps"),...
   'BackgroundColor','white','ContentType','vector');

%% f: HL vs FL, contact frequency of Pitx1 and enhancers, no bootstrap.
ens = [pel, ra3, ra4, pen, neurog];
cFreq = zeros(length(ens),4);
for e=1:length(ens)
    pp_dist = squeeze(hLmaps(pitx,ens(e),:));
    n1 = sum(pp_dist<th);
    n2 = sum(pp_dist<inf);
    cFreq(e,1) = n1;
    cFreq(e,2) = n2;

    pp_dist = squeeze(fLmaps(pitx,ens(e),:));
    n1 = sum(pp_dist<th);
    n2 = sum(pp_dist<inf);
    cFreq(e,3) = n1;
    cFreq(e,4) = n2;
end

% save for ggplot
writematrix(cFreq, [save_folder, 'Fig1f.txt']); % output with Rscript

%%
function output = interval_maximum(starts, ends, values, intervals)
% given a bed-like file for starts, ends as coords and values,
% Put the maximum value of each interval in intervals.
% For segments crossing multiple intervals, the segment is divided into two
% and considered separately for the intervals.
    output = zeros(1, length(intervals));
    ind = 1;
    if starts(1) > intervals(1)
        current_max = values(1);
    else
        current_max = 0;
    end

    for i = 1:size(starts)
        % when end exceeds interval border
        if ends(i) > intervals(ind + 1)
            % register current max
            output(ind) = current_max;
            % carry on current value to the next interval
            current_max = values(i);
            % stop if over all intervals
            if ends(i) > intervals(end)
                output(ind+1 : end) = current_max;
                break;
            end
            ind = ind + 1;
            % keep advancing if end still exceeds border
            while ends(i) > intervals(ind + 1)
                if starts(i) < intervals(ind + 1)
                    output(ind) = current_max;
                end
                ind = ind + 1;
            end
        end
        
        if values(i) > current_max
            current_max = values(i);
        end
    end
end