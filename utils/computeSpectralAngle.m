function [ang_value] = computeSpectralAngle(x,y)
% =========================================================================
% 
% 
% 
% 
% =========================================================================


% ang_value = acos( (squeeze(S_SCLSU_lex2(:,i,n))' * squeeze(Mn_true(:,i,n))) / ...
%                          (norm(squeeze(S_SCLSU_lex2(:,i,n))) * norm(squeeze(Mn_true(:,i,n)))) ) ;
%                      
%                      
                     
ang_value = acos( (squeeze(x)' * squeeze(y)) / ...
             (norm(squeeze(x)) * norm(squeeze(y))) ) ;

