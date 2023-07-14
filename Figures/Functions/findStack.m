function isStack = findStack(Lmaps, th)
    pitx = 22;
    ra3 = 35;
    ra4 = 46;
    pen = 55;

    pb1_inds = Lmaps(pitx, ra3, :) < th;
    pb2_inds = Lmaps(pitx, ra4, :) < th;
    eb1_inds = Lmaps(pen, ra3, :) < th;
    eb2_inds = Lmaps(pen, ra4, :) < th;
    
    isStack = squeeze(pb1_inds | pb2_inds | eb1_inds | eb2_inds);
%     if lenient
%         hasP1 = squeeze(nanMaps(pitx + 1, :) < nanCrit);
%         hasP2 = squeeze(nanMaps(pitx - 1, :) < nanCrit);
%         hasE1 = squeeze(nanMaps(pen + 1, :) < nanCrit);
%         hasE2 = squeeze(nanMaps(pen - 1, :) < nanCrit);
%         hasDataStack = (hasP | hasP1 | hasP2) & (hasE | hasE1 | hasE2) & hasB1B2;
%     else
%         hasDataStack = hasP & hasE & hasB1B2;
%     end

    % substitute pitx value if not present
%     if lenient
%         % find those without P
%         % if has both, average, if only one, use that
%         PFixInds = ~hasP & (hasP1 | hasP2);
%         Lmaps(pitx,:,PFixInds) = mean(cat(1,Lmaps(pitx-1,:,PFixInds), Lmaps(pitx+1,:,PFixInds)), 1, 'omitnan');
%         Lmaps(:,pitx,PFixInds) = mean(cat(2,Lmaps(:,pitx-1,PFixInds), Lmaps(:,pitx+1,PFixInds)), 2, 'omitnan');
%         EFixInds = ~hasE & (hasE1 | hasE2);
%         Lmaps(pen,:,EFixInds) = mean(cat(1,Lmaps(pen-1,:,EFixInds), Lmaps(pen+1,:,EFixInds)), 1, 'omitnan');
%         Lmaps(:,pen,EFixInds) = mean(cat(2,Lmaps(:,pen-1,EFixInds), Lmaps(:,pen+1,EFixInds)), 2, 'omitnan');
%     end
end