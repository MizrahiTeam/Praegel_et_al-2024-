        function [zHit, zFA, dPrime, c_bias] = dprime(binsize, go_licks, ngo_licks, mice);

            for tt=1:length(go_licks)
                if go_licks(tt)==1;
                    go_licks(tt)=1-1/binsize;
                end
                if ngo_licks(tt)==0;
                    ngo_licks(tt)=1/binsize;
                end
                if go_licks(tt)==0;
                    go_licks(tt)=1/binsize;
                end
                if ngo_licks(tt)==1;
                    ngo_licks(tt)=1-1/binsize;
                end
            end

            %normalize accroding to signal detection theory
            zHit = norminv(go_licks) ;
            zFA = norminv(ngo_licks) ;
            dPrime = zHit - zFA ; %dprime as substractionb of norm hit and fa
            c_bias = -((zHit + zFA)./2); % gain bias as opposite of dprime

            
            
          
            
        end
        
        
    