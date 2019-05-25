	#include "udf.h"

#define thetadegr 19.

#define R 1e-03
#define R2 1e-06
#define pi 3.14159265

#define L 918600.
#define us 0.146
#define D 11.81e-06

DEFINE_PROFILE(flux_profile, t, i)
{ 
	real x[ND_ND]; 
	real y, la, jotheta, theta;
	face_t f;
	
	theta = (thetadegr*pi)/180.;
	la = 0.5 - theta/pi; 
	jotheta = ((D*us)/R)*(0.27*theta*theta + 1.3)*(0.6381-0.2239*(theta-pi/4.)*(theta-pi/4.));
		
	begin_f_loop(f, t)
	{
		F_CENTROID(x, f, t);
		y = x[1];
		F_PROFILE(f, t, i) = -L*jotheta*pow(1.-y*y/R2, -la);	
	}
	end_f_loop(f, t)
} 

