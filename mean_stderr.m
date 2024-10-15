function [means, stderr] = mean_stderr (data)

                stderr = std(data)./sqrt(size(data,1));
                means = mean (data);
end 