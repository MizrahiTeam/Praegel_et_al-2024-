

function [fitted_curve_all,mean_fitted_curve, stderr_fitted_curve] = fit_norm_sigmoid (freqs,lick_all)
for m=1:height(lick_all)
            
            [param,stat]=sigm_fit(1:length(freqs),(lick_all(m,:)),[],[],0);
            
            fitted_curve = stat.ypred;
            fitted_curve = mat2gray(smoothdata(fitted_curve));
            fitted_curve_all (m,:)= fitted_curve;

            all_param (m,:) = param ;

end

    [mean_fitted_curve, stderr_fitted_curve] = mean_stderr (fitted_curve_all);

        mean_fitted_curve = mat2gray(mean_fitted_curve);
      %  stderr_fitted_curve = mat2gray (stderr_fitted_curve)
end   