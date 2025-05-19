function [pc_score_mat,explained_M] = pca2d_nanMat(mat_stack,n_components)

if nargin<2
    n_components = 3;
end

M_per = permute(mat_stack,[3 1 2]);
M_vec = M_per(:,:)';
nan_positions = isnan(M_vec);

% Step 2: Replace NaNs with column mean (ignoring NaNs)
col_means = mean(M_vec, 2,'omitnan'); % Compute mean for each feature
col_means(isnan(col_means))=0;
M_filled = M_vec;
for i = 1:size(M_filled, 2)
    nan_idx = isnan(M_filled(:, i));
    M_filled(nan_idx, i) = col_means(nan_idx); % Replace NaNs with mean
end

% PCA
[~, score_M,~ , ~, explained_M] = pca(M_filled);

score_PCs = score_M(:, 1:n_components); % (1600 x n_components)
score_2d = reshape(score_PCs, size(mat_stack,1), size(mat_stack,2), n_components);
pc_score_mat = score_2d;

pc_score_mat(repmat(nan_positions(:,1),[1 n_components])) = NaN;
 
end
