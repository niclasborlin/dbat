function  angles=decons_EulerRotMat(R, seq, isFixed)
% DECONS_EULERROTMAT deconstructs rotation matrix R to
% returns three Euler angles as a 3-by-1 vector. 
	
	% deconstruct R
 
 if nargin<3, isFixed=true; end  
 
switch seq
   case 'xyz' 
   if isFixed
        omega = atan2(R(3,2), R(3,3));
        phi = atan2(-R(3,1), sqrt(R(3,2)*R(3,2) + R(3,3)*R(3,3)));
        kappa = atan2(R(2,1), R(1,1));
        
        angles=[omega;phi;kappa];
   else
       omega =- atan2(R(2,3), R(3,3));
        phi = -atan2(-R(1,3), sqrt(R(2,3)^2 + R(3,3)^2));
        kappa =- atan2(R(1,2), R(1,1));
        
        angles=[omega;phi;kappa];
   end
   
    case 'zxz'
     if isFixed
        alpha = atan2(R(3,1), R(3,2));
        beta = acos(R(3,3));
        gamma = atan2(R(1,3), -R(2,3));
        
        angles=[alpha;beta;gamma];  
    else
       alpha = atan2(R(1,3), -R(2,3));
        beta = acos(R(3,3));
        gamma = atan2(R(3,1), R(3,2));
        
        angles=[alpha;beta;gamma];  
    end
    
    case 'ats'
   
       azimuth = atan2(-R(3,1), -R(3,2));
       % tilt= atan2(sqrt(R(3,1)^2 + R(3,2)^2),R(3,3));
       tilt=acos(R(3,3));
        swing = atan2(-R(1,3),- R(2,3));
        
        angles=[azimuth;tilt;swing];    
end
