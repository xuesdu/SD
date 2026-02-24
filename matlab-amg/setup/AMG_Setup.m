function [ amgData ] = AMG_Setup( Af, amgParam )
% Setup phase for AMG method
%
% @ Xiaozhe Hu, Tufts University

%----------------
% local variable
%----------------
amg_type = amgParam.amg_type;

switch amg_type
   
    case 'UA'
        [ amgData ] = AMG_Setup_UA( Af, amgParam );
        
    case 'SA'
        [ amgData ] = AMG_Setup_SA( Af, amgParam );
        
    case 'C'
        [ amgData ] = AMG_Setup_Classical( Af, amgParam );
        
    otherwise 
        disp('Wrong amg type!');
        disp('Use default type: unsmoothed aggregation');
        [ amgData ] = AMG_Setup_UA( Af, amgParam );
    
end

end

