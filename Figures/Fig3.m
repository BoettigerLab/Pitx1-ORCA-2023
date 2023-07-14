%% 20230712
%% Examine centrality for all
save_folder = 'U:\Manuscripts\TzuChiao Paper\NatureGenetics\Revision1\Revised_images/'; 

%% One needs to run Fig2.m to get the traces categorized into Merge, Stack or Other.
%% b: Distance plot
% Determine distance to center
[hDist, hrG] = getCDistRg(hLpol);
v_reconstructed = zeros(75,1); 
v_reconstructed(badHybes) = 1;

% get ns
sum(hasData & ~isStack)
sum(isStack)
% Boot strapping for median estimate. This part is slow
stack_dists = hDist(isStack, :);
non_dists = hDist(hasData & ~isStack, :);
all_dists = hDist(hasData, :);

boot_time = 1000;
distMeds = zeros(boot_time, 75, 3);

% 95% interval
% 1000 times. %resampling N = N
rng(26); 
s_inds = randi(size(stack_dists, 1), boot_time, size(stack_dists, 1));
n_inds = randi(size(non_dists, 1), boot_time, size(non_dists, 1));
a_inds = randi(size(all_dists, 1), boot_time, size(all_dists, 1));

for i = 1:boot_time
    distMeds(i,:,1) = median(stack_dists(s_inds(i,:), :), 1, 'omitnan');
    distMeds(i,:,2) = median(non_dists(n_inds(i,:), :), 1, 'omitnan');
    distMeds(i,:,3) = median(all_dists(a_inds(i,:), :), 1, 'omitnan');
end
% median plot
distMedSum = zeros(75, 3, 3); %[L, M, U] * [stack, non, all]
Prts = [5, 50, 95];
for i = 1:3
    distMedSum(:,i,:) = prctile(distMeds, Prts(i));
    for j = 1:3
        distMedSum(:,i,j) = interp1(setdiff(1:75, badHybes), distMedSum(find(~v_reconstructed), i, j), 1:75);
    end
end
% writematrix(distMedSum, [save_folder, 'Fig3_th100_CenterDist.txt']); %
% output with Rscript
%% Toy plot if one doesn't want to plot in R.
figure(); clf; 
stack_plot = distMedSum(:,2,1)'; 
all_plot = distMedSum(:,2,3)'; 
non_plot = distMedSum(:,2,2)'; 
stack_interp = interp1(setdiff(1:75, badHybes), stack_plot(find(~v_reconstructed)), 1:75);
all_interp = interp1(setdiff(1:75, badHybes), all_plot(find(~v_reconstructed)), 1:75);
non_interp = interp1(setdiff(1:75, badHybes), non_plot(find(~v_reconstructed)), 1:75);

plot(stack_interp, 'r'); hold on;
plot(all_interp, 'k');hold on;
plot(non_interp, 'b');
ylabel('median distance from center (nm)');
set(gcf,'color','w');
title('');  
xlim([1,75]);
%ylim([0.4, 0.85]);



%% c: Produce graph for the scatter plot
th = 200; % change to 100 to produce Extended Data Fig. 6; Remember to change back!
stripe = 1:75;
hybRate = 1-squeeze(sum(squeeze(sum(isnan(hLmaps(stripe,stripe,:)))) == (length(stripe)-1)))/(length(stripe)-1);
goodInd = hybRate > 0.8;
goodHL = hLmaps(:,:,goodInd);
goodCen = reshape(hDist(goodInd,~v_reconstructed),1,[]);
ContactRate = reshape(squeeze(sum(goodHL(~v_reconstructed, ~v_reconstructed,:) < th))',1,[]);

% Scatter plot
figure();
ids = ~isnan(goodCen);
dscatter(goodCen(ids)', ContactRate(ids)', 'MARKER', 's', ...
    'smoothing', 'None', 'MSIZE', 10, 'FILLED', true, ...
    'BINS', [200, length(unique(ContactRate(ids)))]);
sum(ids) %get n
set(gcf,'position',[200,200,500,300]);
colorbar;

% Obtain log-log correlation.
[R,P, RL, RU] = corrcoef(log(goodCen(ids)),log(ContactRate(ids)'));
R % r value
P % p value
%%
exportgraphics(gcf,strcat(save_folder,"Fig3_scatter.eps"),...
   'BackgroundColor','white','ContentType','vector');

%% c: examine the quartiles:
bracket_n = 4;
lr_range = 20; % 20 barcodes * 10kb per barcode
hDist_clean = hDist(goodInd,:);
hDist_clean(:, badHybes) = nan;

qs = [1:(bracket_n - 1)] /bracket_n;
Qs = quantile(hDist_clean, qs, "all");
Qs = [0; Qs; max(hDist_clean, [], "all")];
% prepare for plotting
cbounds = [0.1, 0.4];
rm = flipud(GetColorMap('redToWhiteK'));

% pad hLmaps
pad_size = size(hLmaps,1) + lr_range*2;
hLmaps_pad = nan(pad_size, pad_size, size(hLmaps,3));
hLmaps_pad((lr_range+1:lr_range+size(hLmaps,1)), (lr_range+1:lr_range+size(hLmaps,1)), :) = hLmaps;
hLmaps_pad(badHybes + lr_range, badHybes + lr_range, :) = nan;
%only good Ind
hLmaps_pad = hLmaps_pad(:,:,goodInd);
% find indices
for b = 1:bracket_n
    b_inds = (hDist_clean > Qs(b)) & (hDist_clean < Qs(b+1));
    n_b = sum(b_inds, 'all');
    %
    [row_inds, col_inds] = find(b_inds);
    col_inds = col_inds + lr_range;
    
    % Calculate the size of the submatrix
    subSize = lr_range*2 + 1;
    output = nan(subSize, subSize, n_b);
    tic;
    for ind = 1:n_b
        barcodes = (col_inds(ind) - lr_range):(col_inds(ind) + lr_range);
        output(:,:,ind) = hLmaps_pad(barcodes, barcodes, row_inds(ind));
    end
    toc;
    %
    [output1, ~] = ContactFrac(output,'threshold',th);
    f1 = figure(b);
    imagesc(output1);
    colormap(rm);
    axis square;
    clim(cbounds);

   %  exportgraphics(f1,strcat(save_folder,"Fig3_th100_PR", num2str(b),"_4.eps"),...
   % 'BackgroundColor','white','ContentType','vector');
end

%% Fig. 3d: Plotting a subpopulation
stack_map = ContactFrac(hLmaps(:,:,isStack), 'threshold', th);
stack_map = InterpMapNans(stack_map,'badHybes',badHybes,'badPixels',[12,24]);
nonStack_map = ContactFrac(hLmaps(:,:,hasData & ~isStack), 'threshold', th);
nonStack_map = InterpMapNans(nonStack_map,'badHybes',badHybes,'badPixels',[12,24]);

cmap = GetColorMap('redToWhiteK');
cbounds = [0.1, 0.4];


figure(1);
imagesc(stack_map);
colormap(flipud(cmap));
clim(cbounds);
colorbar();
axis square;
% saveas(gcf, strcat(save_folder, "fig3_stackTraces.epsc"));

figure(2);
imagesc(nonStack_map); colorbar;
set(gcf,'color','w'); 
axis square; 
colormap(flipud(cmap));
clim(cbounds);
% saveas(gcf, strcat(save_folder, "fig3_nonStack.epsc"));

%% f: RA4-Pitx1: graph
hasData = hLmaps(ra4, pitx, :) < inf;
PB2 = hLmaps(ra4, pitx, :) < th;
sum(hasData) %15057
sum(PB2) %2981
PB2map = ContactFrac(hLmaps(:,:,PB2),'threshold',th);
PB2map = InterpMapNans(PB2map,'badHybes',badHybes,'badPixels',[12,24]);
pop_map = ContactFrac(hLmaps(:,:,hasData),'threshold',th);
pop_map = InterpMapNans(pop_map,'badHybes',badHybes,'badPixels',[12,24]);
% writematrix([PB2map(:,ra4), pop_map(:,ra4)], [save_folder, 'Fig3_th100_PB2.txt']); % output with Rscript

%% g: individual contact maps
cmap = [255 0 0;255 255 255;200 200 200]./255;
IDs = [54,77,147];

N = 0;
for showID = IDs
    sc_map = hLmaps(:,:,showID);
    nanInd = find(sum(isnan(sc_map)) >= 74);
    
    for i = 1:size(sc_map,1)
        if ismember(i,nanInd)
            sc_map(i,i) = nan;
        end
    end
    
    show_map = sc_map;
    show_map(sc_map < th) = 2;
    show_map(sc_map > th) = 1;
    show_map(isnan(sc_map)) = 0;
    
    figure(showID); imagesc(show_map); colormap(flipud(cmap));
    axis square; %clim(axis_dim);
    hold on
    hm = mesh([45.4:1.2:46.6], [0.4:46.2:46.6], zeros([2, 2]));
    hm.FaceColor = 'none';
    hm.EdgeColor = 'k';
    
end