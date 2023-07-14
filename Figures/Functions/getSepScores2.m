function sep = getSepScores2(a, border, sz, triuInter)
    % a is pairwise distancemap normalized by pairwise distances.
    % border is the barcode where sep score is to be obtained.
    % sz is the size on each side to consider for intra/inter interaction
    % freq.
    % we compare border-sz+1:border and border+1:border+sz
    % triuInter, if true, takes only the lower triangle of the interacting
    % box.
    b1 = (border - sz + 1):border;
    b2 = (border + 1):(border + sz);
    L  = a(b1,b1,:);
    L(L==0) = nan;
    R = a(b2,b2,:);
    R(R==0) = nan;
    X = a(b1,b2,:);

    if triuInter
        for i = 1:min(size(X, 1), size(X, 2))
            X(i, i:size(X, 2), :) = nan;
        end
    end

    X(X==0) = nan;

    inter = squeeze(median(X, [1,2], 'omitnan'));
    intra = squeeze(median([reshape(L, [], size(L,3)); reshape(R, [], size(R,3))], 1, 'omitnan'))';
    sep = inter./intra;
end