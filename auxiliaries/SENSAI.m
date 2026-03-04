function [SIGNAL_subspace_similarity, NOISE_subspace_similarity, SENSAI_score] = SENSAI(artifact_threshold, refCOV, Eval, Evec, noise_multiplier, cov_total, refCOV_triu, signal_type)

                       %   Evaluates GEDAI cleaning quality for a given threshold.
%%   Creative Commons License
%
%   Credits:  Tomas Ros & Abele Michela 
%             NeuroTuning Lab [ https://github.com/neurotuning ]
%             Center for Biomedical Imaging
%             University of Geneva
%             Switzerland
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% 3. Neither the name of the copyright holder nor the names of its CONTRIBUTORS
% may be used to endorse or promote products derived from this software without
% specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

[cov_signal_epoched, cov_noise_epoched] = clean_SENSAI(artifact_threshold, refCOV, Eval, Evec, cov_total, signal_type);

%% Estimate Signal Quality
%% Estimate Signal Quality
num_chans = size(refCOV, 1);

% Top eigenvectors of reference covariance (Calculated outside and passed in)

num_epochs = size(cov_signal_epoched, 3);

SIGNAL_subspace_similarity_distribution = zeros(1, num_epochs);
NOISE_subspace_similarity_distribution = zeros(1, num_epochs);

if strcmpi(signal_type, 'meg')
    SSI_top_PCs = 4;
else
    SSI_top_PCs = 3;
end

[Vref, Dref] = eig(refCOV);
[~, idxRef] = sort(diag(Dref), 'descend');
basis_ref = Vref(:, idxRef(1:SSI_top_PCs));

clean_angles = zeros(num_epochs, SSI_top_PCs);
artifact_angles = zeros(num_epochs, SSI_top_PCs);

for epoch = 1:num_epochs
    % SIGNAL SUBSPACE
    cov_signal = cov_signal_epoched(:,:,epoch);
    [Vsig, Dsig] = eig(cov_signal);
    [~, idxSig] = sort(diag(Dsig), 'descend');
    basis_sig = Vsig(:, idxSig(1:SSI_top_PCs));
    cos_theta_sig = subspace_angles(basis_sig, basis_ref);
    clean_angles(epoch, :) = cos_theta_sig(:)';
    
    % NOISE SUBSPACE
    cov_noise = cov_noise_epoched(:,:,epoch);
    [Vnoise, Dnoise] = eig(cov_noise);
    [~, idxNoise] = sort(diag(Dnoise), 'descend');
    basis_noise = Vnoise(:, idxNoise(1:SSI_top_PCs));
    cos_theta_noise = subspace_angles(basis_noise, basis_ref);
    artifact_angles(epoch, :) = cos_theta_noise(:)';
end

%% Compute SENSAI Score (Classic SSI - Principal Angle Product)
SIGNAL_subspace_similarity = mean(prod(clean_angles, 2));
NOISE_subspace_similarity  = mean(prod(artifact_angles, 2));
SENSAI_score = SIGNAL_subspace_similarity - noise_multiplier * NOISE_subspace_similarity;

% %% Compute SENSAI Score (LDA AUC)
% 
% % Combine clean and artifact angles for classification
% X_lda = [clean_angles; artifact_angles];
% 
% % Labels: 1 = Clean, 0 = Artifact
% Y_lda = [ones(num_epochs, 1); zeros(num_epochs, 1)];
% 
% % Train a cross-validated LDA model using adaptive K based on minority class size.
% % (If too few artifact epochs exist, hard-coded K=5 causes stratification failures.)
% n_clean = num_epochs;
% n_artifact = num_epochs;
% n_total = n_clean + n_artifact;
% 
% % Cap K so that:
% % 1. K <= 5 (max folds)
% % 2. K <= min class size (stratification requirement)
% % 3. K <= floor(n_total / 3): ensures each training fold has at least 3 samples (> 2 classes)
% 
% K = min([5, min(n_clean, n_artifact), floor(n_total / 3)]);
% 
% if K < 2
%     % Not enough samples in one class for any cross-validation
%     lda_auc = 0.5;
% else
%     try
%         lda_model = fitcdiscr(X_lda, Y_lda, 'CrossVal', 'on', 'KFold', K, 'DiscrimType', 'pseudoLinear');
% 
%         % Predict out-of-fold probabilities
%         [~, cv_posteriors] = kfoldPredict(lda_model);
% 
%         % cv_posteriors(:,2) contains probabilities for class 1 (Clean)
%         scores = cv_posteriors(:, 2);
% 
%         % Calculate cross-validated Area Under the ROC Curve
%         [~, ~, ~, lda_auc] = perfcurve(Y_lda, scores, 1);
%     catch
%         % Fallback if CV still fails (e.g. all predictions are identical)
%         lda_auc = 0.5;
%     end
% end
% 
% % Pass dummy outputs for external logging
% SIGNAL_subspace_similarity = lda_auc; 
% NOISE_subspace_similarity = lda_auc;
% 
% % Calculate the Euclidean distances to the (1,1,1) leadfield reference point
% ref_point = ones(1, SSI_top_PCs);
% mu_artifact = mean(artifact_angles, 1);
% mu_clean = mean(clean_angles, 1);
% 
% distance_artifact_to_leadfield = sqrt(sum((mu_artifact - ref_point).^2));
% distance_clean_to_leadfield = sqrt(sum((mu_clean - ref_point).^2));
% 
% % The wrapper "SENSAIObjective" in SENSAI_fminbnd.m minimizes the NEGATIVE of this output.
% % Therefore, we MUST return a mathematically positive score to be maximized.
% % 1. MAXIMIZE the cross-validated LDA AUC (Separability)
% % 2. MAXIMIZE the distance of artifacts from the leadfield (Artifacts are far from signal)
% % 3. MINIMIZE the distance of clean epochs from the leadfield (Clean data targets the signal)
% 
% SENSAI_score = lda_auc + (noise_multiplier * distance_artifact_to_leadfield) - distance_clean_to_leadfield;
end