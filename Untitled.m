% Goal: Calculate forces moments and motions for the STAFLO semi-submersible listed in Hooft 1972 Thesis. 
% After programming, the results will be compared with figures IV - 4 on page 51. 
% The numbers on page 51 will be read by me, the idea is to make a rough comparison, because I do not have the original data. 
% Info about the dimensions of the semi will be limited as well, but, Prof said I must take into consideration of the main elements. For things not clear from the graph included on page
% 50, do not include. 

% Below is a summary of the information I get from the graph on page 50. I assume the graph is to scale.
% It is helpful to give the elements some names in programming, descriptions and discussions.
% Naming section; Picture stored at	 C:\Users\Fulin\Desktop\Thesis\September Progress Meeting\STAFLO	; this picture also features the coordinate system.
% The semi has left and right symmetry (underwater body); the semi almost have front aft symmetry for the underbody, the part that made it not fore-aft symmetric is a few pipes
% whose dimensions are unknown, however, I will make smething up to break the fore-aft symmetry.
% Total displacement volume is 12,700 square meters. The elements (in my program) are pontoons, columns and braces, and a deck box with calculated dimensions (to balance the semi). 
% four pontoons in the shape of cylinders; on graph each cylinder measures 4 mm in diameter, so I assume the pontoons have a uniform diameter of 6.03 meters, which is taken from the graph
% Four corners are four big cylinders, each with a diameter of 7.05 meters;
% most other cylinders have diameter of 4.05 meters. c8 and c13 have a smaller diameter of 2.7 m (on the graph, it measures 2mm and c15 measures 3 mm which equals to 4.05 m in diameter.)
% Draft: 18.00 meters
% 
%

% Coordiante System
%	Origin: Centr point of the keel plane, lies in the middle of the pontoons
%	Axis: 	x-from origin to bow section;
%			y-from origin to the right/starboard side;
%			z-from origin downwards. (following a right-handed system)



% Basic Numbers [SI]
g=9.81;
rho=1025;	%density of sea water 



% Some [assumptions] made
%	The elements in real life case intersect each other, and the total volume is smaller than the summation of all of the cylinder volumes; the difference will be disregarded.

% This time the dimension are not going to be loaded from an external file, instead they shall be defined below:
%	
%	[1] dimensions of and about [pontoons]
Lp1 = 64.01;
Lp2 = Lp1  ;
Lp3 = Lp1  ;
Lp4 = Lp1  ;

Rp1 = 6.03/2;
Rp2 = Rp1	;
Rp3 = Rp1	;
Rp4 = Rp1	;

Dp1 = 53.34;
Dp2 = 25.91;

%	[2] draft
d = 18;
%	[3] dimension of and about the [columns] (only the underwater part)
Lc1  = d-6.03;
Lc5  = d-6.03;
Lc16 = d-6.03;
Lc20 = d-6.03;

Rc1  = 7.05/2;
Rc5  = 7.05/2;
Rc16 = 7.05/2;
Rc20 = 7.05/2;

Lc2   = d-6.03;
Lc3   = d-6.03;
Lc4   = d-6.03;
Lc6   = d-6.03;
Lc7   = d-6.03;
Lc9   = d-6.03;
Lc10  = d-6.03;
Lc11  = d-6.03;
Lc12  = d-6.03;
Lc14  = d-6.03;
Lc15  = d-6.03;
Lc17  = d-6.03;
Lc18  = d-6.03;
Lc19  = d-6.03;

Lc8	  = d-6.03;
Lc13  = d-6.03;

Rc2   = 4.05/2;
Rc3   = 4.05/2;
Rc4   = 4.05/2;
Rc6   = 4.05/2;
Rc7   = 4.05/2;
Rc9   = 4.05/2;
Rc10  = 4.05/2;
Rc11  = 4.05/2;
Rc12  = 4.05/2;
Rc14  = 4.05/2;
Rc15  = 4.05/2;
Rc17  = 4.05/2;
Rc18  = 4.05/2;
Rc19  = 4.05/2;

Rc8	  = 2.7/2;
Rc13  = 2.7/2;

Dc1  = 64.01; 
Dc2  = 32.00;

%	[4] dimensions of [braces]
Lb1 = Dp1;
Lb2 = Lb1 		 ;

Rb1  =1.68/2;
Rb2  =1.68/2;

Lb3 = (Dp1^2 + Dc1^2)^0.5;
Lb4 = Lb3 		 ;

Rb3  = 1.68/2;
Rb4  = 1.68/2;

% Matrix L will be defined below for ease of use later, the same will be done for radius numbers
% L(1,1)  till L(1,20) holds column length; L(1,21) till L(1,24) holds pontoon length; L(1,25) till L(1,28) holds brace length
% R(1,1)  till R(1,20) holds column radius; R(1,21) till R(1,24) holds pontoon radius; R(1,25) till R(1,28) holds brace radius
L(1,1)		=	Lc1		;     R(1,1)		=	Rc1		;
L(1,2)		=	Lc2		;     R(1,2)		=	Rc2		;
L(1,3)		=	Lc3		;     R(1,3)		=	Rc3		;
L(1,4)		=	Lc4		;     R(1,4)		=	Rc4		;
L(1,5)		=	Lc5		;     R(1,5)		=	Rc5		;
L(1,6)		=	Lc6		;     R(1,6)		=	Rc6		;
L(1,7)		=	Lc7		;     R(1,7)		=	Rc7		;
L(1,8)		=	Lc8		;     R(1,8)		=	Rc8		;
L(1,9)		=	Lc9		;     R(1,9)		=	Rc9		;
L(1,10)		=	Lc10	;     R(1,10)		=	Rc10	;
L(1,11)		=	Lc11	;     R(1,11)		=	Rc11	;
L(1,12)		=	Lc12	;     R(1,12)		=	Rc12	;
L(1,13)		=	Lc13	;     R(1,13)		=	Rc13	;
L(1,14)		=	Lc14	;     R(1,14)		=	Rc14	;
L(1,15)		=	Lc15	;     R(1,15)		=	Rc15	;
L(1,16)		=	Lc16	;     R(1,16)		=	Rc16	;
L(1,17)		=	Lc17	;     R(1,17)		=	Rc17	;
L(1,18)		=	Lc18	;     R(1,18)		=	Rc18	;
L(1,19)		=	Lc19	;     R(1,19)		=	Rc19	;
L(1,20)		=	Lc20	;     R(1,20)		=	Rc20	;
L(1,21)		=	Lp1  	;     R(1,21)		=	Rp1  	;
L(1,22)		=	Lp2  	;     R(1,22)		=	Rp2  	;
L(1,23)		=	Lp3  	;     R(1,23)		=	Rp3  	;
L(1,24)		=	Lp4  	;     R(1,24)		=	Rp4  	;
L(1,25)		=	Lb1  	;     R(1,25)		=	Rb1  	;
L(1,26)		=	Lb2  	;     R(1,26)		=	Rb2  	;
L(1,27)		=	Lb3  	;     R(1,27)		=	Rb3  	;
L(1,28)		=	Lb4  	;     R(1,28)		=	Rb4  	;
% With all of these dimensions defined, it is time to come up with the [end points] of the elements, which are all cylinders.

%	Pontoons are always in the x direction, and the end ppints are defined with either a (aft) or b (bow). The things in the () are meant for understanding.
%	Columns aer always in the z direction in the semi, end points are defined with l (low) or h(high).
% 	braces are always in the y direction here, end points are defined with p (port side) or s (starboard side).

% 	End points for [pontoons] [4 pontoons in total]
%	The name of the points are in the form of "p1a(1,2)", each point is a (1 * 3)(one row, three column)matrix. 
np = 4;
p1a(1,1) = -Dc1/2	;	p1a(1,2) = -Dp1/2	;	p1a(1,3) = -Rp1		;
p2a(1,1) = -Dc1/2	;	p2a(1,2) = -Dp2/2	;	p2a(1,3) = -Rp2		;
p3a(1,1) = -Dc1/2	;	p3a(1,2) = Dp2/2	;	p3a(1,3) = -Rp3		;
p4a(1,1) = -Dc1/2	;	p4a(1,2) = Dp1/2	;	p4a(1,3) = -Rp4		;

p1b(1,1) =  Dc1/2	;	p1b(1,2) = -Dp1/2	;	p1b(1,3) = -Rp1		;
p2b(1,1) =  Dc1/2	;	p2b(1,2) = -Dp2/2	;	p2b(1,3) = -Rp2		;
p3b(1,1) =  Dc1/2	;	p3b(1,2) = Dp2/2	;	p3b(1,3) = -Rp3		;
p4b(1,1) =  Dc1/2	;	p4b(1,2) = Dp1/2	;	p4b(1,3) = -Rp4		;

% 	End points for [columns] [20 columns in total]
%	The name of the points are in the form of "c1h(1,2)", each point is a (1 * 3)(one row, three column)matrix.
%	The 'h' in the denotation means high, closer to the deck. 
nc=20;
c1h(1,1)	= -Dc1/2	;	 c1h(1,2)	= Dp1/2		;   c1h(1,3)	= -d	;
c2h(1,1)	= -Dc2/2	;	 c2h(1,2)	= Dp1/2		;   c2h(1,3)	= -d	;
c3h(1,1)	= 0			;	 c3h(1,2)	= Dp1/2		;   c3h(1,3)	= -d	;
c4h(1,1)	= Dc2/2	    ;	 c4h(1,2)	= Dp1/2		;   c4h(1,3)	= -d	;
c5h(1,1)	= Dc1/2	    ;	 c5h(1,2)	= Dp1/2		;   c5h(1,3)	= -d	;
c6h(1,1)	= -Dc1/2	;	 c6h(1,2)	= Dp2/2		;   c6h(1,3)	= -d	;
c7h(1,1)	= -Dc2/2	;	 c7h(1,2)	= Dp2/2		;   c7h(1,3)	= -d	;
c8h(1,1)	= 0			;	 c8h(1,2)	= Dp2/2		;   c8h(1,3)	= -d	;
c9h(1,1)	= Dc2/2		;	 c9h(1,2)	= Dp2/2		;   c9h(1,3)	= -d	;
c10h(1,1)	= Dc1/2		;	 c10h(1,2)	= Dp2/2		;   c10h(1,3)	= -d	;
c11h(1,1)	= -Dc1/2	;	 c11h(1,2)	= -Dp2/2	;   c11h(1,3)	= -d	;
c12h(1,1)	= -Dc2/2	;	 c12h(1,2)	= -Dp2/2	;   c12h(1,3)	= -d	;
c13h(1,1)	= 0			;	 c13h(1,2)	= -Dp2/2	;   c13h(1,3)	= -d	;
c14h(1,1)	= Dc2/2	 	;	 c14h(1,2)	= -Dp2/2	;   c14h(1,3)	= -d	;
c15h(1,1)	= Dc1/2	 	;	 c15h(1,2)	= -Dp2/2	;   c15h(1,3)	= -d	;
c16h(1,1)	= -Dc1/2	;	 c16h(1,2)	= -Dp1/2	;   c16h(1,3)	= -d	;
c17h(1,1)	= -Dc2/2	;	 c17h(1,2)	= -Dp1/2	;   c17h(1,3)	= -d	;
c18h(1,1)	= 0			;	 c18h(1,2)	= -Dp1/2	;   c18h(1,3)	= -d	;
c19h(1,1)	= Dc2/2		;	 c19h(1,2)	= -Dp1/2	;   c19h(1,3)	= -d	;
c20h(1,1)	= Dc1/2		;	 c20h(1,2)	= -Dp1/2	;   c20h(1,3)	= -d	;
                                                      
c1l(1,1)	= -Dc1/2	;	 c1l(1,2)	= Dp1/2 ;   	c1l(1,3)	= -6.03;
c2l(1,1)	= -Dc2/2	;	 c2l(1,2)	= Dp1/2 ;   	c2l(1,3)	= -6.03;
c3l(1,1)	= 0			;	 c3l(1,2)	= Dp1/2 ;   	c3l(1,3)	= -6.03;
c4l(1,1)	= Dc2/2	    ;	 c4l(1,2)	= Dp1/2 ;   	c4l(1,3)	= -6.03;
c5l(1,1)	= Dc1/2	    ;	 c5l(1,2)	= Dp1/2 ;   	c5l(1,3)	= -6.03;
c6l(1,1)	= -Dc1/2	;	 c6l(1,2)	= Dp2/2 ;  		c6l(1,3)	= -6.03;
c7l(1,1)	= -Dc2/2	;	 c7l(1,2)	= Dp2/2 ;  		c7l(1,3)	= -6.03;
c8l(1,1)	= 0			;	 c8l(1,2)	= Dp2/2 ;  		c8l(1,3)	= -6.03;
c9l(1,1)	= Dc2/2		;	 c9l(1,2)	= Dp2/2 ;  		c9l(1,3)	= -6.03;
c10l(1,1)	= Dc1/2		;	 c10l(1,2)	= Dp2/2 ;  		c10l(1,3)	= -6.03;
c11l(1,1)	= -Dc1/2	;	 c11l(1,2)	= -Dp2/2	;   c11l(1,3)	= -6.03;
c12l(1,1)	= -Dc2/2	;	 c12l(1,2)	= -Dp2/2	;   c12l(1,3)	= -6.03;
c13l(1,1)	= 0			;	 c13l(1,2)	= -Dp2/2	;   c13l(1,3)	= -6.03;
c14l(1,1)	= Dc2/2	 	;	 c14l(1,2)	= -Dp2/2	;   c14l(1,3)	= -6.03;
c15l(1,1)	= Dc1/2	 	;	 c15l(1,2)	= -Dp2/2	;   c15l(1,3)	= -6.03;
c16l(1,1)	= -Dc1/2	;	 c16l(1,2)	= -Dp1/2	;   c16l(1,3)	= -6.03;
c17l(1,1)	= -Dc2/2	;	 c17l(1,2)	= -Dp1/2	;   c17l(1,3)	= -6.03;
c18l(1,1)	= 0			;	 c18l(1,2)	= -Dp1/2	;   c18l(1,3)	= -6.03;
c19l(1,1)	= Dc2/2		;	 c19l(1,2)	= -Dp1/2	;   c19l(1,3)	= -6.03;
c20l(1,1)	= Dc1/2		;	 c20l(1,2)	= -Dp1/2	;   c20l(1,3)	= -6.03;

% 	End points for [Braces] [4 braces in total]
%	The name of the points are in the form of "b1p(1,2)", each point is a (1 * 3)(one row, three column)matrix.
%	No direct reading of z position of the braces, again, assuming the graph is to scale, I took a measure and decided the position is 9 meters above keel plane.
nb=4;
b1p(1,1) = -Dc2/2	;	b1p(1,2) = -Dp1/2	;	b1p(1,3) = -9		;
b2p(1,1) = Dc2/2	;	b2p(1,2) = -Dp1/2	;	b2p(1,3) = -9		;
b3p(1,1) = Dc1/2	;	b3p(1,2) = -Dp1/2	;	b3p(1,3) = -9		;
b4p(1,1) = -Dc1/2	;	b4p(1,2) = -Dp1/2	;	b4p(1,3) = -9		;
                                                            
b1s(1,1) =  -Dc2/2	;	b1s(1,2) = Dp1/2	;	b1s(1,3) = -9		;
b2s(1,1) =  Dc2/2	;	b2s(1,2) = Dp1/2	;	b2s(1,3) = -9		;
b3s(1,1) =  Dc1/2	;	b3s(1,2) = Dp1/2	;	b3s(1,3) = -9		;
b4s(1,1) =  -Dc1/2	;	b4s(1,2) = Dp1/2	;	b4s(1,3) = -9		;



% The volume of a cylinder is calculated: r^2 * L; the summation of all elements in this model will result in the displacement of STAFLO.
%	This variable is named "dis"; dis will be the summation of all elements underwater, which will be the summation of sum_disp, sum_disb and sum_disc.
%	sum_disp, sum_disb and sum_disc will not be calcualted with a for loop.

sum_disp = (Lp1*Rp1*Rp1+Lp2*Rp2*Rp2+Lp3*Rp3*Rp3+Lp4*Rp4*Rp4)*pi	;

sum_disb = (Lb1*Rb1*Rb1	+Lb2*Rb2*Rb2	+Lb3*Rb3*Rb3	+Lb4*Rb4*Rb4)*pi	;

sum_disc = (Lc1 *Rc1 *Rc1 	+Lc5 *Rc5 *Rc5 	+Lc16*Rc16*Rc16	+Lc20*Rc20*Rc20	+Lc2 *Rc2 *Rc2 	+Lc3 *Rc3 *Rc3 	+Lc4 *Rc4 *Rc4 	+Lc6 *Rc6 *Rc6 + Lc7 *Rc7 *Rc7 	+Lc9 *Rc9 *Rc9 	+Lc11*Rc11*Rc11+Lc15*Rc15*Rc15+Lc12*Rc12*Rc12	+Lc10*Rc10*Rc10	+Lc14*Rc14*Rc14	+Lc17*Rc17*Rc17	+Lc18*Rc18*Rc18	+Lc19*Rc19*Rc19	+Lc8*Rc8*	Rc8	+Lc13*Rc13*Rc13)*pi	;

dis = sum_disp+sum_disb+sum_disc	;
% A brief discussion: first run of the code give a good result which is very close to the displayed number on page 50 Hooft thesis. So, dis will be equalized to 12700.
dis = 12700;

% The next step, itis handy to get the [center of buoyancy]
%	The vessel features, starboard and port side symmetry and with the main elements, it feaures fore and aft symmetry. So, the x and y coordinates are zeros.
%	The task remains to obtain the z coordinate of the semi (a negative number in our coordinate system).
%	General Idea: calculate the first moment of mass(volume is also good in this case) with respect to the keel plane (origin as well) of each element, make a summation and devide by the total mass, the outcome of 
%	this will be the z coordinate of the COB. If the correct coordinates of each element is plugged in, we will arrive at a number that is around 7 meters above the keel plane. 

%cob(1,3)= (-sum_disp*Rp1 -pi*(Lc1 *Rc1 *Rc1 	+Lc5 *Rc5 *Rc5 	+Lc16*Rc16*Rc16	+Lc20*Rc20*Rc20	)*d/2-pi*(Lc8*Rc8*	Rc8	+Lc13*Rc13*Rc13)*(d-0.5*Rp1)-pi*(Lc2 *Rc2 *Rc2 	+Lc3 *Rc3 *Rc3 	+Lc4 *Rc4 *Rc4 	+Lc7 *Rc7 *Rc7 	+Lc9 *Rc9 *Rc9 	+Lc12*Rc12*Rc12	+Lc14*Rc14*Rc14	+Lc17*Rc17*Rc17	+Lc18*Rc18*Rc18)*(d-0.5*Rp1)-9*sum_disc)	;

cob(1,3)=(	sum_disp * p1a(1,3) + ...
			sum_disb * b1p(1,3) + ...
			pi*(Lc1*Rc1*Rc1+Lc5*Rc5*Rc5+Lc16*Rc16*Rc16+Lc20*Rc20*Rc20)*0.5*(c1l(1,3)+c1h(1,3)) + ...
			pi*(Lc2*Rc2*Rc2+Lc3*Rc3*Rc3+Lc4*Rc4*Rc4+Lc7*Rc7*Rc7+Lc9*Rc9*Rc9+Lc12*Rc12*Rc12+Lc14*Rc14*Rc14+Lc17*Rc17*Rc17+...
			Lc18*Rc18*Rc18+Lc19*Rc19*Rc19+Lc8*Rc8*Rc8+Lc13*Rc13*Rc13)* 0.5*(c2l(1,3)+c2h(1,3))) / dis	;
%Brief Record: first run give out number of around -5.80 which is reasonalbe(the other cob(1,3) is not correctly programmed and was commented out.).


% [Metacentric Height] (transverse and longitudinal) will be calcualted below.
%	All of the elements contributing to this number are braces with round profile. Off diagnal elements are not taken into consideration because 
%	lack of knowledge about their dimensions (this will not be mentioned further into the code).
%	Calculation of this number also requires knowledge of area moment of area formula for circles. They could be found on the internet with ease. 
%	Given a circle with a radius of r and an x axis going through the center, the formulas are:
%		I_x = I_y = pi/4*r^4		I_z = pi/2*r^4 
%	The variables are defined with the name I4xx, the '4' here indicates the dimension of value is to the 4th power of Length. 

I4xx=	pi/4*Rc1^4  + pi* Rc1^2	 * c1h(1,2).^2  	+...
        pi/4*Rc2^4  + pi* Rc2^2	 * c2h(1,2).^2  	+...
        pi/4*Rc3^4  + pi* Rc3^2	 * c3h(1,2).^2  	+...
        pi/4*Rc4^4  + pi* Rc4^2	 * c4h(1,2).^2  	+...
        pi/4*Rc5^4  + pi* Rc5^2	 * c5h(1,2).^2  	+...
        pi/4*Rc6^4  + pi* Rc6^2	 * c6h(1,2).^2  	+...
        pi/4*Rc7^4  + pi* Rc7^2	 * c7h(1,2).^2  	+...
        pi/4*Rc8^4  + pi* Rc8^2	 * c8h(1,2).^2  	+...
        pi/4*Rc9^4  + pi* Rc9^2	 * c9h(1,2).^2  	+...
        pi/4*Rc10^4	+ pi* Rc10^2 * c10h(1,2).^2 	+...
        pi/4*Rc11^4	+ pi* Rc11^2 * c11h(1,2).^2 	+...
        pi/4*Rc12^4	+ pi* Rc12^2 * c12h(1,2).^2 	+...
        pi/4*Rc13^4	+ pi* Rc13^2 * c13h(1,2).^2 	+...
        pi/4*Rc14^4	+ pi* Rc14^2 * c14h(1,2).^2 	+...
        pi/4*Rc15^4	+ pi* Rc15^2 * c15h(1,2).^2 	+...
        pi/4*Rc16^4	+ pi* Rc16^2 * c16h(1,2).^2 	+...
        pi/4*Rc17^4	+ pi* Rc17^2 * c17h(1,2).^2 	+...
        pi/4*Rc18^4	+ pi* Rc18^2 * c18h(1,2).^2 	+...
        pi/4*Rc19^4	+ pi* Rc19^2 * c19h(1,2).^2 	+...
        pi/4*Rc20^4	+ pi* Rc20^2 * c20h(1,2).^2 	;
		
I4yy=	pi/4*Rc1^4  + pi* Rc1^2	 * c1h(1,1).^2  	+...
        pi/4*Rc2^4  + pi* Rc2^2	 * c2h(1,1).^2  	+...
        pi/4*Rc3^4  + pi* Rc3^2	 * c3h(1,1).^2  	+...
        pi/4*Rc4^4  + pi* Rc4^2	 * c4h(1,1).^2  	+...
        pi/4*Rc5^4  + pi* Rc5^2	 * c5h(1,1).^2  	+...
        pi/4*Rc6^4  + pi* Rc6^2	 * c6h(1,1).^2  	+...
        pi/4*Rc7^4  + pi* Rc7^2	 * c7h(1,1).^2  	+...
        pi/4*Rc8^4  + pi* Rc8^2	 * c8h(1,1).^2  	+...
        pi/4*Rc9^4  + pi* Rc9^2	 * c9h(1,1).^2  	+...
        pi/4*Rc10^4	+ pi* Rc10^2 * c10h(1,1).^2 	+...
        pi/4*Rc11^4	+ pi* Rc11^2 * c11h(1,1).^2 	+...
        pi/4*Rc12^4	+ pi* Rc12^2 * c12h(1,1).^2 	+...
        pi/4*Rc13^4	+ pi* Rc13^2 * c13h(1,1).^2 	+...
        pi/4*Rc14^4	+ pi* Rc14^2 * c14h(1,1).^2 	+...
        pi/4*Rc15^4	+ pi* Rc15^2 * c15h(1,1).^2 	+...
        pi/4*Rc16^4	+ pi* Rc16^2 * c16h(1,1).^2 	+...
        pi/4*Rc17^4	+ pi* Rc17^2 * c17h(1,1).^2 	+...
        pi/4*Rc18^4	+ pi* Rc18^2 * c18h(1,1).^2 	+...
        pi/4*Rc19^4	+ pi* Rc19^2 * c19h(1,1).^2 	+...
        pi/4*Rc20^4	+ pi* Rc20^2 * c20h(1,1).^2 	;
		

BMt = I4xx / dis;
BMl = I4yy / dis;
GMt = 6.58;
KB  = cob(1,3);
KG = abs(KB) + BMt - GMt	;
% The Hooft thesis did not give any information on the vertical position of COG, but this number could evaluated from GMt values as did above.



% Next, the mass moment of inertias should be calculated. However, the article did not give any information on the estimated density of the different elements, 
% hence, I did not attempt to calculate the mass moment of inertias (same as the radius of gyrations). The codes above already allows me to calcualte the 
% added mass coefficients, which is the whole idea behind this part of the codes. So, the inertias will be skipped at least for now.
% 	Added mass coefficients are calculated according to page 107 of Hooft thesis.
%	!!!	Due to lack of info on density of the semi, all terms in the code that has md in it are subsituted with ad.
%	For the pontoons, they are all horizontal. So, angle bata and angle gamma equal to 0. and most of the added mass coefficients are zero as one would expect. 
%	For the sake of throughness, I will do the whole thing for a single pontoon, and all the others pontoon will be the same situation. 
%	In Hooft thesis, the sequencese of the points inserted is sensitive. The bow section point will be inserted first, as in Xd2 on page 107, 
%	and the the aft point will be in the position of Xd1.
%	(x0,y0,z0) is the position of COG, in which x0 = y0 = 0, and z0 = -KG
z0=-KG;
x0=0;
y0=0;
a = zeros(6,6);
%	All Angles are Listed Below.
alpha_c1  = asin( (c1l(1,1) -  c1h(1,1) ) / Lc1 )	;
alpha_c2  = asin( (c2l(1,1) -  c2h(1,1) ) / Lc2 )	;
alpha_c3  = asin( (c3l(1,1) -  c3h(1,1) ) / Lc3 )	;
alpha_c4  = asin( (c4l(1,1) -  c4h(1,1) ) / Lc4 )	;
alpha_c5  = asin( (c5l(1,1) -  c5h(1,1) ) / Lc5 )	;
alpha_c6  = asin( (c6l(1,1) -  c6h(1,1) ) / Lc6 )	;
alpha_c7  = asin( (c7l(1,1) -  c7h(1,1) ) / Lc7 )	;
alpha_c8  = asin( (c8l(1,1) -  c8h(1,1) ) / Lc8 )	;
alpha_c9  = asin( (c9l(1,1) -  c9h(1,1) ) / Lc9 )	;
alpha_c10 = asin( (c10l(1,1) - c10h(1,1) ) / Lc10)	;
alpha_c11 = asin( (c11l(1,1) - c11h(1,1) ) / Lc11)	;
alpha_c12 = asin( (c12l(1,1) - c12h(1,1) ) / Lc12)	;
alpha_c13 = asin( (c13l(1,1) - c13h(1,1) ) / Lc13)	;
alpha_c14 = asin( (c14l(1,1) - c14h(1,1) ) / Lc14)	;
alpha_c15 = asin( (c15l(1,1) - c15h(1,1) ) / Lc15)	;
alpha_c16 = asin( (c16l(1,1) - c16h(1,1) ) / Lc16)	;
alpha_c17 = asin( (c17l(1,1) - c17h(1,1) ) / Lc17)	;
alpha_c18 = asin( (c18l(1,1) - c18h(1,1) ) / Lc18)	;
alpha_c19 = asin( (c19l(1,1) - c19h(1,1) ) / Lc19)	;
alpha_c20 = asin( (c20l(1,1) - c20h(1,1) ) / Lc20)	;
alpha_p1  = asin( (p1b(1,1) -  p1a(1,1) ) / Lp1 )	;
alpha_p2  = asin( (p2b(1,1) -  p2a(1,1) ) / Lp2 )	;
alpha_p3  = asin( (p3b(1,1) -  p3a(1,1) ) / Lp3 )	;
alpha_p4  = asin( (p4b(1,1) -  p4a(1,1) ) / Lp4 )	;
alpha_b1  = asin( (b1p(1,1) -  b1s(1,1) ) / Lb1 )	;
alpha_b2  = asin( (b2p(1,1) -  b2s(1,1) ) / Lb2 )	;
alpha_b3  = asin( (b3p(1,1) -  b3s(1,1) ) / Lb3 )	;
alpha_b4  = asin( (b4p(1,1) -  b4s(1,1) ) / Lb4 )	;

beta_c1   = asin( (c1l(1,2) -  c1h(1,2) ) / Lc1 )	;
beta_c2   = asin( (c2l(1,2) -  c2h(1,2) ) / Lc2 )	;
beta_c3   = asin( (c3l(1,2) -  c3h(1,2) ) / Lc3 )	;
beta_c4   = asin( (c4l(1,2) -  c4h(1,2) ) / Lc4 )	;
beta_c5   = asin( (c5l(1,2) -  c5h(1,2) ) / Lc5 )	;
beta_c6   = asin( (c6l(1,2) -  c6h(1,2) ) / Lc6 )	;
beta_c7   = asin( (c7l(1,2) -  c7h(1,2) ) / Lc7 )	;
beta_c8   = asin( (c8l(1,2) -  c8h(1,2) ) / Lc8 )	;
beta_c9   = asin( (c9l(1,2) -  c9h(1,2) ) / Lc9 )	;
beta_c10  = asin( (c10l(1,2) - c10h(1,2) ) / Lc10)	;
beta_c11  = asin( (c11l(1,2) - c11h(1,2) ) / Lc11)	;
beta_c12  = asin( (c12l(1,2) - c12h(1,2) ) / Lc12)	;
beta_c13  = asin( (c13l(1,2) - c13h(1,2) ) / Lc13)	;
beta_c14  = asin( (c14l(1,2) - c14h(1,2) ) / Lc14)	;
beta_c15  = asin( (c15l(1,2) - c15h(1,2) ) / Lc15)	;
beta_c16  = asin( (c16l(1,2) - c16h(1,2) ) / Lc16)	;
beta_c17  = asin( (c17l(1,2) - c17h(1,2) ) / Lc17)	;
beta_c18  = asin( (c18l(1,2) - c18h(1,2) ) / Lc18)	;
beta_c19  = asin( (c19l(1,2) - c19h(1,2) ) / Lc19)	;
beta_c20  = asin( (c20l(1,2) - c20h(1,2) ) / Lc20)	;
beta_p1   = asin( (p1b(1,2) -  p1a(1,2) ) / Lp1 )	;
beta_p2   = asin( (p2b(1,2) -  p2a(1,2) ) / Lp2 )	;
beta_p3   = asin( (p3b(1,2) -  p3a(1,2) ) / Lp3 )	;
beta_p4   = asin( (p4b(1,2) -  p4a(1,2) ) / Lp4 )	;
beta_b1   = asin( (b1p(1,2) -  b1s(1,2) ) / Lb1 )	;
beta_b2   = asin( (b2p(1,2) -  b2s(1,2) ) / Lb2 )	;
beta_b3   = asin( (b3p(1,2) -  b3s(1,2) ) / Lb3 )	;
beta_b4   = asin( (b4p(1,2) -  b4s(1,2) ) / Lb4 )	;
          
gamma_c1  = asin( (c1l(1,3) - c1h(1,3) ) / Lc1 )	;
gamma_c2  = asin( (c2l(1,3) - c2h(1,3) ) / Lc2 )	;
gamma_c3  = asin( (c3l(1,3) - c3h(1,3) ) / Lc3 )	;
gamma_c4  = asin( (c4l(1,3) - c4h(1,3) ) / Lc4 )	;
gamma_c5  = asin( (c5l(1,3) - c5h(1,3) ) / Lc5 )	;
gamma_c6  = asin( (c6l(1,3) - c6h(1,3) ) / Lc6 )	;
gamma_c7  = asin( (c7l(1,3) - c7h(1,3) ) / Lc7 )	;
gamma_c8  = asin( (c8l(1,3) - c8h(1,3) ) / Lc8 )	;
gamma_c9  = asin( (c9l(1,3) - c9h(1,3) ) / Lc9 )	;
gamma_c10 = asin( (c10l(1,3) - c10h(1,3) ) / Lc10)	;
gamma_c11 = asin( (c11l(1,3) - c11h(1,3) ) / Lc11)	;
gamma_c12 = asin( (c12l(1,3) - c12h(1,3) ) / Lc12)	;
gamma_c13 = asin( (c13l(1,3) - c13h(1,3) ) / Lc13)	;
gamma_c14 = asin( (c14l(1,3) - c14h(1,3) ) / Lc14)	;
gamma_c15 = asin( (c15l(1,3) - c15h(1,3) ) / Lc15)	;
gamma_c16 = asin( (c16l(1,3) - c16h(1,3) ) / Lc16)	;
gamma_c17 = asin( (c17l(1,3) - c17h(1,3) ) / Lc17)	;
gamma_c18 = asin( (c18l(1,3) - c18h(1,3) ) / Lc18)	;
gamma_c19 = asin( (c19l(1,3) - c19h(1,3) ) / Lc19)	;
gamma_c20 = asin( (c20l(1,3) - c20h(1,3) ) / Lc20)	;
gamma_p1  = asin( (p1b(1,3) - p1a(1,3) ) / Lp1 )	;
gamma_p2  = asin( (p2b(1,3) - p2a(1,3) ) / Lp2 )	;
gamma_p3  = asin( (p3b(1,3) - p3a(1,3) ) / Lp3 )	;
gamma_p4  = asin( (p4b(1,3) - p4a(1,3) ) / Lp4 )	;
gamma_b1  = asin( (b1p(1,3) - b1s(1,3) ) / Lb1 )	;
gamma_b2  = asin( (b2p(1,3) - b2s(1,3) ) / Lb2 )	;
gamma_b3  = asin( (b3p(1,3) - b3s(1,3) ) / Lb3 )	;
gamma_b4  = asin( (b4p(1,3) - b4s(1,3) ) / Lb4 )	;

%	Below the coeffifients are defined.
ca_c	=1;
ca_p	=1;
ca_b	=1;

%	Below [added mass] of (28) elements perpendicular to their main axis is defined. 
a_c1  	=ca_c * rho * (pi * Rc1^2) * Lc1 		;
a_c2  	=ca_c * rho * (pi * Rc2^2) * Lc2 		;
a_c3  	=ca_c * rho * (pi * Rc3^2) * Lc3 		;
a_c4  	=ca_c * rho * (pi * Rc4^2) * Lc4 		;
a_c5  	=ca_c * rho * (pi * Rc5^2) * Lc5 		;
a_c6  	=ca_c * rho * (pi * Rc6^2) * Lc6 		;
a_c7  	=ca_c * rho * (pi * Rc7^2) * Lc7 		;
a_c8  	=ca_c * rho * (pi * Rc8^2) * Lc8 		;
a_c9  	=ca_c * rho * (pi * Rc9^2) * Lc9 		;
a_c10 	=ca_c * rho * (pi * Rc10^2) * Lc10		;
a_c11 	=ca_c * rho * (pi * Rc11^2) * Lc11		;
a_c12 	=ca_c * rho * (pi * Rc12^2) * Lc12		;
a_c13 	=ca_c * rho * (pi * Rc13^2) * Lc13		;
a_c14 	=ca_c * rho * (pi * Rc14^2) * Lc14		;
a_c15 	=ca_c * rho * (pi * Rc15^2) * Lc15		;
a_c16 	=ca_c * rho * (pi * Rc16^2) * Lc16		;
a_c17 	=ca_c * rho * (pi * Rc17^2) * Lc17		;
a_c18 	=ca_c * rho * (pi * Rc18^2) * Lc18		;
a_c19 	=ca_c * rho * (pi * Rc19^2) * Lc19		;
a_c20 	=ca_c * rho * (pi * Rc20^2) * Lc20		;
a_p1	=ca_p * rho * (pi * Rp1^2) * Lp1		;
a_p2	=ca_p * rho * (pi * Rp2^2) * Lp2		;
a_p3	=ca_p * rho * (pi * Rp3^2) * Lp3		;
a_p4	=ca_p * rho * (pi * Rp4^2) * Lp4		;
a_b1	=ca_b * rho * (pi * Rb1^2) * Lb1		;
a_b2	=ca_b * rho * (pi * Rb2^2) * Lb2		;
a_b3	=ca_b * rho * (pi * Rb3^2) * Lb3		;
a_b4	=ca_b * rho * (pi * Rb4^2) * Lb4		;
 
%	Below the mass of the elements are defined.
%	density will be defined first and appriximated as 500 kg/m^3.
density = 500;
m_c1  	=  density * rho * (pi * Rc1^2) * Lc1 		;
m_c2  	=  density * rho * (pi * Rc2^2) * Lc2 		;
m_c3  	=  density * rho * (pi * Rc3^2) * Lc3 		;
m_c4  	=  density * rho * (pi * Rc4^2) * Lc4 		;
m_c5  	=  density * rho * (pi * Rc5^2) * Lc5 		;
m_c6  	=  density * rho * (pi * Rc6^2) * Lc6 		;
m_c7  	=  density * rho * (pi * Rc7^2) * Lc7 		;
m_c8  	=  density * rho * (pi * Rc8^2) * Lc8 		;
m_c9  	=  density * rho * (pi * Rc9^2) * Lc9 		;
m_c10 	=  density * rho * (pi * Rc10^2) * Lc10		;
m_c11 	=  density * rho * (pi * Rc11^2) * Lc11		;
m_c12 	=  density * rho * (pi * Rc12^2) * Lc12		;
m_c13 	=  density * rho * (pi * Rc13^2) * Lc13		;
m_c14 	=  density * rho * (pi * Rc14^2) * Lc14		;
m_c15 	=  density * rho * (pi * Rc15^2) * Lc15		;
m_c16 	=  density * rho * (pi * Rc16^2) * Lc16		;
m_c17 	=  density * rho * (pi * Rc17^2) * Lc17		;
m_c18 	=  density * rho * (pi * Rc18^2) * Lc18		;
m_c19 	=  density * rho * (pi * Rc19^2) * Lc19		;
m_c20 	=  density * rho * (pi * Rc20^2) * Lc20		;
m_p1	=  density * rho * (pi * Rp1^2) * Lp1		;
m_p2	=  density * rho * (pi * Rp2^2) * Lp2		;
m_p3	=  density * rho * (pi * Rp3^2) * Lp3		;
m_p4	=  density * rho * (pi * Rp4^2) * Lp4		;
m_b1	=  density * rho * (pi * Rb1^2) * Lb1		;
m_b2	=  density * rho * (pi * Rb2^2) * Lb2		;
m_b3	=  density * rho * (pi * Rb3^2) * Lb3		;
m_b4	=  density * rho * (pi * Rb4^2) * Lb4		;

%	Below, [first moment of added mass] is calculated.
s_c1  	= 0.5 * ca_c * rho * (pi * Rc1^2)  * Lc1^2		;
s_c2  	= 0.5 * ca_c * rho * (pi * Rc2^2)  * Lc2^2		;
s_c3  	= 0.5 * ca_c * rho * (pi * Rc3^2)  * Lc3^2		;
s_c4  	= 0.5 * ca_c * rho * (pi * Rc4^2)  * Lc4^2		;
s_c5  	= 0.5 * ca_c * rho * (pi * Rc5^2)  * Lc5^2		;
s_c6  	= 0.5 * ca_c * rho * (pi * Rc6^2)  * Lc6^2		;
s_c7  	= 0.5 * ca_c * rho * (pi * Rc7^2)  * Lc7^2		;
s_c8  	= 0.5 * ca_c * rho * (pi * Rc8^2)  * Lc8^2		;
s_c9  	= 0.5 * ca_c * rho * (pi * Rc9^2)  * Lc9^2		;
s_c10 	= 0.5 * ca_c * rho * (pi * Rc10^2) * Lc10^2		;
s_c11 	= 0.5 * ca_c * rho * (pi * Rc11^2) * Lc11^2		;
s_c12 	= 0.5 * ca_c * rho * (pi * Rc12^2) * Lc12^2		;
s_c13 	= 0.5 * ca_c * rho * (pi * Rc13^2) * Lc13^2		;
s_c14 	= 0.5 * ca_c * rho * (pi * Rc14^2) * Lc14^2		;
s_c15 	= 0.5 * ca_c * rho * (pi * Rc15^2) * Lc15^2		;
s_c16 	= 0.5 * ca_c * rho * (pi * Rc16^2) * Lc16^2		;
s_c17 	= 0.5 * ca_c * rho * (pi * Rc17^2) * Lc17^2		;
s_c18 	= 0.5 * ca_c * rho * (pi * Rc18^2) * Lc18^2		;
s_c19 	= 0.5 * ca_c * rho * (pi * Rc19^2) * Lc19^2		;
s_c20 	= 0.5 * ca_c * rho * (pi * Rc20^2) * Lc20^2		;
s_p1	= 0.5 * ca_p * rho * (pi * Rp1^2)  * Lp1^2		;
s_p2	= 0.5 * ca_p * rho * (pi * Rp2^2)  * Lp2^2		;
s_p3	= 0.5 * ca_p * rho * (pi * Rp3^2)  * Lp3^2		;
s_p4	= 0.5 * ca_p * rho * (pi * Rp4^2)  * Lp4^2		;
s_b1	= 0.5 * ca_b * rho * (pi * Rb1^2)  * Lb1^2		;
s_b2	= 0.5 * ca_b * rho * (pi * Rb2^2)  * Lb2^2		;
s_b3	= 0.5 * ca_b * rho * (pi * Rb3^2)  * Lb3^2		;
s_b4	= 0.5 * ca_b * rho * (pi * Rb4^2)  * Lb4^2		;

%	Below, [second moment of added mass] is calculated.
I_c1  	= 1/3 * ca_c * rho * (pi * Rc1^2)  * Lc1^3		;
I_c2  	= 1/3 * ca_c * rho * (pi * Rc2^2)  * Lc2^3		;
I_c3  	= 1/3 * ca_c * rho * (pi * Rc3^2)  * Lc3^3		;
I_c4  	= 1/3 * ca_c * rho * (pi * Rc4^2)  * Lc4^3		;
I_c5  	= 1/3 * ca_c * rho * (pi * Rc5^2)  * Lc5^3		;
I_c6  	= 1/3 * ca_c * rho * (pi * Rc6^2)  * Lc6^3		;
I_c7  	= 1/3 * ca_c * rho * (pi * Rc7^2)  * Lc7^3		;
I_c8  	= 1/3 * ca_c * rho * (pi * Rc8^2)  * Lc8^3		;
I_c9  	= 1/3 * ca_c * rho * (pi * Rc9^2)  * Lc9^3		;
I_c10 	= 1/3 * ca_c * rho * (pi * Rc10^2) * Lc10^3		;
I_c11 	= 1/3 * ca_c * rho * (pi * Rc11^2) * Lc11^3		;
I_c12 	= 1/3 * ca_c * rho * (pi * Rc12^2) * Lc12^3		;
I_c13 	= 1/3 * ca_c * rho * (pi * Rc13^2) * Lc13^3		;
I_c14 	= 1/3 * ca_c * rho * (pi * Rc14^2) * Lc14^3		;
I_c15 	= 1/3 * ca_c * rho * (pi * Rc15^2) * Lc15^3		;
I_c16 	= 1/3 * ca_c * rho * (pi * Rc16^2) * Lc16^3		;
I_c17 	= 1/3 * ca_c * rho * (pi * Rc17^2) * Lc17^3		;
I_c18 	= 1/3 * ca_c * rho * (pi * Rc18^2) * Lc18^3		;
I_c19 	= 1/3 * ca_c * rho * (pi * Rc19^2) * Lc19^3		;
I_c20 	= 1/3 * ca_c * rho * (pi * Rc20^2) * Lc20^3		;
I_p1	= 1/3 * ca_p * rho * (pi * Rp1^2)  * Lp1^3		;
I_p2	= 1/3 * ca_p * rho * (pi * Rp2^2)  * Lp2^3		;
I_p3	= 1/3 * ca_p * rho * (pi * Rp3^2)  * Lp3^3		;
I_p4	= 1/3 * ca_p * rho * (pi * Rp4^2)  * Lp4^3		;
I_b1	= 1/3 * ca_b * rho * (pi * Rb1^2)  * Lb1^3		;
I_b2	= 1/3 * ca_b * rho * (pi * Rb2^2)  * Lb2^3		;
I_b3	= 1/3 * ca_b * rho * (pi * Rb3^2)  * Lb3^3		;
I_b4	= 1/3 * ca_b * rho * (pi * Rb4^2)  * Lb4^3		;

%	Below a11 added mass of elements are defined in a single array (so that later it is easier to add them all up.)
a11_c1 		= a_c1  * ( cos( alpha_c1  ) )^2	;
a11_c2 		= a_c2  * ( cos( alpha_c2  ) )^2	;
a11_c3 		= a_c3  * ( cos( alpha_c3  ) )^2	;
a11_c4 		= a_c4  * ( cos( alpha_c4  ) )^2	;
a11_c5 		= a_c5  * ( cos( alpha_c5  ) )^2	;
a11_c6 		= a_c6  * ( cos( alpha_c6  ) )^2	;
a11_c7 		= a_c7  * ( cos( alpha_c7  ) )^2	;
a11_c8 		= a_c8  * ( cos( alpha_c8  ) )^2	;
a11_c9 		= a_c9  * ( cos( alpha_c9  ) )^2	;
a11_c10		= a_c10 * ( cos( alpha_c10 ) )^2	;
a11_c11		= a_c11 * ( cos( alpha_c11 ) )^2	;
a11_c12		= a_c12 * ( cos( alpha_c12 ) )^2	;
a11_c13		= a_c13 * ( cos( alpha_c13 ) )^2	;
a11_c14		= a_c14 * ( cos( alpha_c14 ) )^2	;
a11_c15		= a_c15 * ( cos( alpha_c15 ) )^2	;
a11_c16		= a_c16 * ( cos( alpha_c16 ) )^2	;
a11_c17		= a_c17 * ( cos( alpha_c17 ) )^2	;
a11_c18		= a_c18 * ( cos( alpha_c18 ) )^2	;
a11_c19		= a_c19 * ( cos( alpha_c19 ) )^2	;
a11_c20		= a_c20 * ( cos( alpha_c20 ) )^2	;
a11_p1 		= a_p1	* ( cos( alpha_p1  ) )^2	;
a11_p2 		= a_p2	* ( cos( alpha_p2  ) )^2	;
a11_p3 		= a_p3	* ( cos( alpha_p3  ) )^2	;
a11_p4 		= a_p4	* ( cos( alpha_p4  ) )^2	;
a11_b1 		= a_b1	* ( cos( alpha_b1  ) )^2	;
a11_b2 		= a_b2	* ( cos( alpha_b2  ) )^2	;
a11_b3 		= a_b3	* ( cos( alpha_b3  ) )^2	;
a11_b4 		= a_b4	* ( cos( alpha_b4  ) )^2	;

a11(1,1) = a11_c1 		;
a11(1,2) = a11_c2 		;
a11(1,3) = a11_c3 		;
a11(1,4) = a11_c4 		;
a11(1,5) = a11_c5 		;
a11(1,6) = a11_c6 		;
a11(1,7) = a11_c7 		;
a11(1,8) = a11_c8 		;
a11(1,9) = a11_c9 		;
a11(1,10)= a11_c10		;
a11(1,11)= a11_c11		;
a11(1,12)= a11_c12		;
a11(1,13)= a11_c13		;
a11(1,14)= a11_c14		;
a11(1,15)= a11_c15		;
a11(1,16)= a11_c16		;
a11(1,17)= a11_c17		;
a11(1,18)= a11_c18		;
a11(1,19)= a11_c19		;
a11(1,20)= a11_c20		;
a11(1,21)= a11_p1 		;
a11(1,22)= a11_p2 		;
a11(1,23)= a11_p3 		;
a11(1,24)= a11_p4 		;
a11(1,25)= a11_b1 		;
a11(1,26)= a11_b2 		;
a11(1,27)= a11_b3 		;
a11(1,28)= a11_b4 		;

a12_c1 		= -a_c1  * sin( beta_c1  )	* sin( alpha_c1  )	;
a12_c2 		= -a_c2  * sin( beta_c2  )	* sin( alpha_c2  )	;
a12_c3 		= -a_c3  * sin( beta_c3  )	* sin( alpha_c3  )	;
a12_c4 		= -a_c4  * sin( beta_c4  )	* sin( alpha_c4  )	;
a12_c5 		= -a_c5  * sin( beta_c5  )	* sin( alpha_c5  )	;
a12_c6 		= -a_c6  * sin( beta_c6  )	* sin( alpha_c6  )	;
a12_c7 		= -a_c7  * sin( beta_c7  )	* sin( alpha_c7  )	;
a12_c8 		= -a_c8  * sin( beta_c8  )	* sin( alpha_c8  )	;
a12_c9 		= -a_c9  * sin( beta_c9  )	* sin( alpha_c9  )	;
a12_c10		= -a_c10 * sin( beta_c10 )	* sin( alpha_c10 )	;
a12_c11		= -a_c11 * sin( beta_c11 )	* sin( alpha_c11 )	;
a12_c12		= -a_c12 * sin( beta_c12 )	* sin( alpha_c12 )	;
a12_c13		= -a_c13 * sin( beta_c13 )	* sin( alpha_c13 )	;
a12_c14		= -a_c14 * sin( beta_c14 )	* sin( alpha_c14 )	;
a12_c15		= -a_c15 * sin( beta_c15 )	* sin( alpha_c15 )	;
a12_c16		= -a_c16 * sin( beta_c16 )	* sin( alpha_c16 )	;
a12_c17		= -a_c17 * sin( beta_c17 )	* sin( alpha_c17 )	;
a12_c18		= -a_c18 * sin( beta_c18 )	* sin( alpha_c18 )	;
a12_c19		= -a_c19 * sin( beta_c19 )	* sin( alpha_c19 )	;
a12_c20		= -a_c20 * sin( beta_c20 )	* sin( alpha_c20 )	;
a12_p1 		= -a_p1	 * sin( beta_p1  )	* sin( alpha_p1  )	;
a12_p2 		= -a_p2	 * sin( beta_p2  )	* sin( alpha_p2  )	;
a12_p3 		= -a_p3	 * sin( beta_p3  )	* sin( alpha_p3  )	;
a12_p4 		= -a_p4	 * sin( beta_p4  )	* sin( alpha_p4  )	;
a12_b1 		= -a_b1	 * sin( beta_b1  )	* sin( alpha_b1  )	;
a12_b2 		= -a_b2	 * sin( beta_b2  )	* sin( alpha_b2  )	;
a12_b3 		= -a_b3	 * sin( beta_b3  )	* sin( alpha_b3  )	;
a12_b4 		= -a_b4	 * sin( beta_b4  )	* sin( alpha_b4  )	;
  
a12(1,1) = a12_c1 		;
a12(1,2) = a12_c2 		;
a12(1,3) = a12_c3 		;
a12(1,4) = a12_c4 		;
a12(1,5) = a12_c5 		;
a12(1,6) = a12_c6 		;
a12(1,7) = a12_c7 		;
a12(1,8) = a12_c8 		;
a12(1,9) = a12_c9 		;
a12(1,10)= a12_c10		;
a12(1,11)= a12_c11		;
a12(1,12)= a12_c12		;
a12(1,13)= a12_c13		;
a12(1,14)= a12_c14		;
a12(1,15)= a12_c15		;
a12(1,16)= a12_c16		;
a12(1,17)= a12_c17		;
a12(1,18)= a12_c18		;
a12(1,19)= a12_c19		;
a12(1,20)= a12_c20		;
a12(1,21)= a12_p1 		;
a12(1,22)= a12_p2 		;
a12(1,23)= a12_p3 		;
a12(1,24)= a12_p4 		;
a12(1,25)= a12_b1 		;
a12(1,26)= a12_b2 		;
a12(1,27)= a12_b3 		;
a12(1,28)= a12_b4 		;

a21_c1 		= a12_c1 		;
a21_c2 		= a12_c2 		;
a21_c3 		= a12_c3 		;
a21_c4 		= a12_c4 		;
a21_c5 		= a12_c5 		;
a21_c6 		= a12_c6 		;
a21_c7 		= a12_c7 		;
a21_c8 		= a12_c8 		;
a21_c9 		= a12_c9 		;
a21_c10		= a12_c10		;
a21_c11		= a12_c11		;
a21_c12		= a12_c12		;
a21_c13		= a12_c13		;
a21_c14		= a12_c14		;
a21_c15		= a12_c15		;
a21_c16		= a12_c16		;
a21_c17		= a12_c17		;
a21_c18		= a12_c18		;
a21_c19		= a12_c19		;
a21_c20		= a12_c20		;
a21_p1 		= a12_p1 		;
a21_p2 		= a12_p2 		;
a21_p3 		= a12_p3 		;
a21_p4 		= a12_p4 		;
a21_b1 		= a12_b1 		;
a21_b2 		= a12_b2 		;
a21_b3 		= a12_b3 		;
a21_b4 		= a12_b4 		;
  
a21(1,1) = a21_c1 		;
a21(1,2) = a21_c2 		;
a21(1,3) = a21_c3 		;
a21(1,4) = a21_c4 		;
a21(1,5) = a21_c5 		;
a21(1,6) = a21_c6 		;
a21(1,7) = a21_c7 		;
a21(1,8) = a21_c8 		;
a21(1,9) = a21_c9 		;
a21(1,10)= a21_c10		;
a21(1,11)= a21_c11		;
a21(1,12)= a21_c12		;
a21(1,13)= a21_c13		;
a21(1,14)= a21_c14		;
a21(1,15)= a21_c15		;
a21(1,16)= a21_c16		;
a21(1,17)= a21_c17		;
a21(1,18)= a21_c18		;
a21(1,19)= a21_c19		;
a21(1,20)= a21_c20		;
a21(1,21)= a21_p1 		;
a21(1,22)= a21_p2 		;
a21(1,23)= a21_p3 		;
a21(1,24)= a21_p4 		;
a21(1,25)= a21_b1 		;
a21(1,26)= a21_b2 		;
a21(1,27)= a21_b3 		;
a21(1,28)= a21_b4 		;

a13_c1 		= -a_c1  * sin( gamma_c1  )	* sin( alpha_c1  )	;
a13_c2 		= -a_c2  * sin( gamma_c2  )	* sin( alpha_c2  )	;
a13_c3 		= -a_c3  * sin( gamma_c3  )	* sin( alpha_c3  )	;
a13_c4 		= -a_c4  * sin( gamma_c4  )	* sin( alpha_c4  )	;
a13_c5 		= -a_c5  * sin( gamma_c5  )	* sin( alpha_c5  )	;
a13_c6 		= -a_c6  * sin( gamma_c6  )	* sin( alpha_c6  )	;
a13_c7 		= -a_c7  * sin( gamma_c7  )	* sin( alpha_c7  )	;
a13_c8 		= -a_c8  * sin( gamma_c8  )	* sin( alpha_c8  )	;
a13_c9 		= -a_c9  * sin( gamma_c9  )	* sin( alpha_c9  )	;
a13_c10		= -a_c10 * sin( gamma_c10 )	* sin( alpha_c10 )	;
a13_c11		= -a_c11 * sin( gamma_c11 )	* sin( alpha_c11 )	;
a13_c12		= -a_c12 * sin( gamma_c12 )	* sin( alpha_c12 )	;
a13_c13		= -a_c13 * sin( gamma_c13 )	* sin( alpha_c13 )	;
a13_c14		= -a_c14 * sin( gamma_c14 )	* sin( alpha_c14 )	;
a13_c15		= -a_c15 * sin( gamma_c15 )	* sin( alpha_c15 )	;
a13_c16		= -a_c16 * sin( gamma_c16 )	* sin( alpha_c16 )	;
a13_c17		= -a_c17 * sin( gamma_c17 )	* sin( alpha_c17 )	;
a13_c18		= -a_c18 * sin( gamma_c18 )	* sin( alpha_c18 )	;
a13_c19		= -a_c19 * sin( gamma_c19 )	* sin( alpha_c19 )	;
a13_c20		= -a_c20 * sin( gamma_c20 )	* sin( alpha_c20 )	;
a13_p1 		= -a_p1	 * sin( gamma_p1  )	* sin( alpha_p1  )	;
a13_p2 		= -a_p2	 * sin( gamma_p2  )	* sin( alpha_p2  )	;
a13_p3 		= -a_p3	 * sin( gamma_p3  )	* sin( alpha_p3  )	;
a13_p4 		= -a_p4	 * sin( gamma_p4  )	* sin( alpha_p4  )	;
a13_b1 		= -a_b1	 * sin( gamma_b1  )	* sin( alpha_b1  )	;
a13_b2 		= -a_b2	 * sin( gamma_b2  )	* sin( alpha_b2  )	;
a13_b3 		= -a_b3	 * sin( gamma_b3  )	* sin( alpha_b3  )	;
a13_b4 		= -a_b4	 * sin( gamma_b4  )	* sin( alpha_b4  )	;
  
a13(1,1) = a13_c1 		;
a13(1,2) = a13_c2 		;
a13(1,3) = a13_c3 		;
a13(1,4) = a13_c4 		;
a13(1,5) = a13_c5 		;
a13(1,6) = a13_c6 		;
a13(1,7) = a13_c7 		;
a13(1,8) = a13_c8 		;
a13(1,9) = a13_c9 		;
a13(1,10)= a13_c10		;
a13(1,11)= a13_c11		;
a13(1,12)= a13_c12		;
a13(1,13)= a13_c13		;
a13(1,14)= a13_c14		;
a13(1,15)= a13_c15		;
a13(1,16)= a13_c16		;
a13(1,17)= a13_c17		;
a13(1,18)= a13_c18		;
a13(1,19)= a13_c19		;
a13(1,20)= a13_c20		;
a13(1,21)= a13_p1 		;
a13(1,22)= a13_p2 		;
a13(1,23)= a13_p3 		;
a13(1,24)= a13_p4 		;
a13(1,25)= a13_b1 		;
a13(1,26)= a13_b2 		;
a13(1,27)= a13_b3 		;
a13(1,28)= a13_b4 		;

a31_c1 		= a13_c1 		;
a31_c2 		= a13_c2 		;
a31_c3 		= a13_c3 		;
a31_c4 		= a13_c4 		;
a31_c5 		= a13_c5 		;
a31_c6 		= a13_c6 		;
a31_c7 		= a13_c7 		;
a31_c8 		= a13_c8 		;
a31_c9 		= a13_c9 		;
a31_c10		= a13_c10		;
a31_c11		= a13_c11		;
a31_c12		= a13_c12		;
a31_c13		= a13_c13		;
a31_c14		= a13_c14		;
a31_c15		= a13_c15		;
a31_c16		= a13_c16		;
a31_c17		= a13_c17		;
a31_c18		= a13_c18		;
a31_c19		= a13_c19		;
a31_c20		= a13_c20		;
a31_p1 		= a13_p1 		;
a31_p2 		= a13_p2 		;
a31_p3 		= a13_p3 		;
a31_p4 		= a13_p4 		;
a31_b1 		= a13_b1 		;
a31_b2 		= a13_b2 		;
a31_b3 		= a13_b3 		;
a31_b4 		= a13_b4 		;

a31(1,1) = a31_c1 		;
a31(1,2) = a31_c2 		;
a31(1,3) = a31_c3 		;
a31(1,4) = a31_c4 		;
a31(1,5) = a31_c5 		;
a31(1,6) = a31_c6 		;
a31(1,7) = a31_c7 		;
a31(1,8) = a31_c8 		;
a31(1,9) = a31_c9 		;
a31(1,10)= a31_c10		;
a31(1,11)= a31_c11		;
a31(1,12)= a31_c12		;
a31(1,13)= a31_c13		;
a31(1,14)= a31_c14		;
a31(1,15)= a31_c15		;
a31(1,16)= a31_c16		;
a31(1,17)= a31_c17		;
a31(1,18)= a31_c18		;
a31(1,19)= a31_c19		;
a31(1,20)= a31_c20		;
a31(1,21)= a31_p1 		;
a31(1,22)= a31_p2 		;
a31(1,23)= a31_p3 		;
a31(1,24)= a31_p4 		;
a31(1,25)= a31_b1 		;
a31(1,26)= a31_b2 		;
a31(1,27)= a31_b3 		;
a31(1,28)= a31_b4 		;

a22_c1 		= a_c1  * ( cos( beta_c1  ) )^2	;
a22_c2 		= a_c2  * ( cos( beta_c2  ) )^2	;
a22_c3 		= a_c3  * ( cos( beta_c3  ) )^2	;
a22_c4 		= a_c4  * ( cos( beta_c4  ) )^2	;
a22_c5 		= a_c5  * ( cos( beta_c5  ) )^2	;
a22_c6 		= a_c6  * ( cos( beta_c6  ) )^2	;
a22_c7 		= a_c7  * ( cos( beta_c7  ) )^2	;
a22_c8 		= a_c8  * ( cos( beta_c8  ) )^2	;
a22_c9 		= a_c9  * ( cos( beta_c9  ) )^2	;
a22_c10		= a_c10 * ( cos( beta_c10 ) )^2	;
a22_c11		= a_c11 * ( cos( beta_c11 ) )^2	;
a22_c12		= a_c12 * ( cos( beta_c12 ) )^2	;
a22_c13		= a_c13 * ( cos( beta_c13 ) )^2	;
a22_c14		= a_c14 * ( cos( beta_c14 ) )^2	;
a22_c15		= a_c15 * ( cos( beta_c15 ) )^2	;
a22_c16		= a_c16 * ( cos( beta_c16 ) )^2	;
a22_c17		= a_c17 * ( cos( beta_c17 ) )^2	;
a22_c18		= a_c18 * ( cos( beta_c18 ) )^2	;
a22_c19		= a_c19 * ( cos( beta_c19 ) )^2	;
a22_c20		= a_c20 * ( cos( beta_c20 ) )^2	;
a22_p1 		= a_p1	* ( cos( beta_p1  ) )^2	;
a22_p2 		= a_p2	* ( cos( beta_p2  ) )^2	;
a22_p3 		= a_p3	* ( cos( beta_p3  ) )^2	;
a22_p4 		= a_p4	* ( cos( beta_p4  ) )^2	;
a22_b1 		= a_b1	* ( cos( beta_b1  ) )^2	;
a22_b2 		= a_b2	* ( cos( beta_b2  ) )^2	;
a22_b3 		= a_b3	* ( cos( beta_b3  ) )^2	;
a22_b4 		= a_b4	* ( cos( beta_b4  ) )^2	;
 
a22(1,1) = a22_c1 		;
a22(1,2) = a22_c2 		;
a22(1,3) = a22_c3 		;
a22(1,4) = a22_c4 		;
a22(1,5) = a22_c5 		;
a22(1,6) = a22_c6 		;
a22(1,7) = a22_c7 		;
a22(1,8) = a22_c8 		;
a22(1,9) = a22_c9 		;
a22(1,10)= a22_c10		;
a22(1,11)= a22_c11		;
a22(1,12)= a22_c12		;
a22(1,13)= a22_c13		;
a22(1,14)= a22_c14		;
a22(1,15)= a22_c15		;
a22(1,16)= a22_c16		;
a22(1,17)= a22_c17		;
a22(1,18)= a22_c18		;
a22(1,19)= a22_c19		;
a22(1,20)= a22_c20		;
a22(1,21)= a22_p1 		;
a22(1,22)= a22_p2 		;
a22(1,23)= a22_p3 		;
a22(1,24)= a22_p4 		;
a22(1,25)= a22_b1 		;
a22(1,26)= a22_b2 		;
a22(1,27)= a22_b3 		;
a22(1,28)= a22_b4 		;

a23_c1 		= -a_c1  * sin( gamma_c1  )	* sin( beta_c1  )	;
a23_c2 		= -a_c2  * sin( gamma_c2  )	* sin( beta_c2  )	;
a23_c3 		= -a_c3  * sin( gamma_c3  )	* sin( beta_c3  )	;
a23_c4 		= -a_c4  * sin( gamma_c4  )	* sin( beta_c4  )	;
a23_c5 		= -a_c5  * sin( gamma_c5  )	* sin( beta_c5  )	;
a23_c6 		= -a_c6  * sin( gamma_c6  )	* sin( beta_c6  )	;
a23_c7 		= -a_c7  * sin( gamma_c7  )	* sin( beta_c7  )	;
a23_c8 		= -a_c8  * sin( gamma_c8  )	* sin( beta_c8  )	;
a23_c9 		= -a_c9  * sin( gamma_c9  )	* sin( beta_c9  )	;
a23_c10		= -a_c10 * sin( gamma_c10 )	* sin( beta_c10 )	;
a23_c11		= -a_c11 * sin( gamma_c11 )	* sin( beta_c11 )	;
a23_c12		= -a_c12 * sin( gamma_c12 )	* sin( beta_c12 )	;
a23_c13		= -a_c13 * sin( gamma_c13 )	* sin( beta_c13 )	;
a23_c14		= -a_c14 * sin( gamma_c14 )	* sin( beta_c14 )	;
a23_c15		= -a_c15 * sin( gamma_c15 )	* sin( beta_c15 )	;
a23_c16		= -a_c16 * sin( gamma_c16 )	* sin( beta_c16 )	;
a23_c17		= -a_c17 * sin( gamma_c17 )	* sin( beta_c17 )	;
a23_c18		= -a_c18 * sin( gamma_c18 )	* sin( beta_c18 )	;
a23_c19		= -a_c19 * sin( gamma_c19 )	* sin( beta_c19 )	;
a23_c20		= -a_c20 * sin( gamma_c20 )	* sin( beta_c20 )	;
a23_p1 		= -a_p1	 * sin( gamma_p1  )	* sin( beta_p1  )	;
a23_p2 		= -a_p2	 * sin( gamma_p2  )	* sin( beta_p2  )	;
a23_p3 		= -a_p3	 * sin( gamma_p3  )	* sin( beta_p3  )	;
a23_p4 		= -a_p4	 * sin( gamma_p4  )	* sin( beta_p4  )	;
a23_b1 		= -a_b1	 * sin( gamma_b1  )	* sin( beta_b1  )	;
a23_b2 		= -a_b2	 * sin( gamma_b2  )	* sin( beta_b2  )	;
a23_b3 		= -a_b3	 * sin( gamma_b3  )	* sin( beta_b3  )	;
a23_b4 		= -a_b4	 * sin( gamma_b4  )	* sin( beta_b4  )	;

a23(1,1) = a23_c1 		;
a23(1,2) = a23_c2 		;
a23(1,3) = a23_c3 		;
a23(1,4) = a23_c4 		;
a23(1,5) = a23_c5 		;
a23(1,6) = a23_c6 		;
a23(1,7) = a23_c7 		;
a23(1,8) = a23_c8 		;
a23(1,9) = a23_c9 		;
a23(1,10)= a23_c10		;
a23(1,11)= a23_c11		;
a23(1,12)= a23_c12		;
a23(1,13)= a23_c13		;
a23(1,14)= a23_c14		;
a23(1,15)= a23_c15		;
a23(1,16)= a23_c16		;
a23(1,17)= a23_c17		;
a23(1,18)= a23_c18		;
a23(1,19)= a23_c19		;
a23(1,20)= a23_c20		;
a23(1,21)= a23_p1 		;
a23(1,22)= a23_p2 		;
a23(1,23)= a23_p3 		;
a23(1,24)= a23_p4 		;
a23(1,25)= a23_b1 		;
a23(1,26)= a23_b2 		;
a23(1,27)= a23_b3 		;
a23(1,28)= a23_b4 		;

a32_c1 		= a23_c1 		;
a32_c2 		= a23_c2 		;
a32_c3 		= a23_c3 		;
a32_c4 		= a23_c4 		;
a32_c5 		= a23_c5 		;
a32_c6 		= a23_c6 		;
a32_c7 		= a23_c7 		;
a32_c8 		= a23_c8 		;
a32_c9 		= a23_c9 		;
a32_c10		= a23_c10		;
a32_c11		= a23_c11		;
a32_c12		= a23_c12		;
a32_c13		= a23_c13		;
a32_c14		= a23_c14		;
a32_c15		= a23_c15		;
a32_c16		= a23_c16		;
a32_c17		= a23_c17		;
a32_c18		= a23_c18		;
a32_c19		= a23_c19		;
a32_c20		= a23_c20		;
a32_p1 		= a23_p1 		;
a32_p2 		= a23_p2 		;
a32_p3 		= a23_p3 		;
a32_p4 		= a23_p4 		;
a32_b1 		= a23_b1 		;
a32_b2 		= a23_b2 		;
a32_b3 		= a23_b3 		;
a32_b4 		= a23_b4 		;

a32(1,1) = a32_c1 		;
a32(1,2) = a32_c2 		;
a32(1,3) = a32_c3 		;
a32(1,4) = a32_c4 		;
a32(1,5) = a32_c5 		;
a32(1,6) = a32_c6 		;
a32(1,7) = a32_c7 		;
a32(1,8) = a32_c8 		;
a32(1,9) = a32_c9 		;
a32(1,10)= a32_c10		;
a32(1,11)= a32_c11		;
a32(1,12)= a32_c12		;
a32(1,13)= a32_c13		;
a32(1,14)= a32_c14		;
a32(1,15)= a32_c15		;
a32(1,16)= a32_c16		;
a32(1,17)= a32_c17		;
a32(1,18)= a32_c18		;
a32(1,19)= a32_c19		;
a32(1,20)= a32_c20		;
a32(1,21)= a32_p1 		;
a32(1,22)= a32_p2 		;
a32(1,23)= a32_p3 		;
a32(1,24)= a32_p4 		;
a32(1,25)= a32_b1 		;
a32(1,26)= a32_b2 		;
a32(1,27)= a32_b3 		;
a32(1,28)= a32_b4 		;

a33_c1 		= a_c1  * ( cos( gamma_c1  ) )^2	;
a33_c2 		= a_c2  * ( cos( gamma_c2  ) )^2	;
a33_c3 		= a_c3  * ( cos( gamma_c3  ) )^2	;
a33_c4 		= a_c4  * ( cos( gamma_c4  ) )^2	;
a33_c5 		= a_c5  * ( cos( gamma_c5  ) )^2	;
a33_c6 		= a_c6  * ( cos( gamma_c6  ) )^2	;
a33_c7 		= a_c7  * ( cos( gamma_c7  ) )^2	;
a33_c8 		= a_c8  * ( cos( gamma_c8  ) )^2	;
a33_c9 		= a_c9  * ( cos( gamma_c9  ) )^2	;
a33_c10		= a_c10 * ( cos( gamma_c10 ) )^2	;
a33_c11		= a_c11 * ( cos( gamma_c11 ) )^2	;
a33_c12		= a_c12 * ( cos( gamma_c12 ) )^2	;
a33_c13		= a_c13 * ( cos( gamma_c13 ) )^2	;
a33_c14		= a_c14 * ( cos( gamma_c14 ) )^2	;
a33_c15		= a_c15 * ( cos( gamma_c15 ) )^2	;
a33_c16		= a_c16 * ( cos( gamma_c16 ) )^2	;
a33_c17		= a_c17 * ( cos( gamma_c17 ) )^2	;
a33_c18		= a_c18 * ( cos( gamma_c18 ) )^2	;
a33_c19		= a_c19 * ( cos( gamma_c19 ) )^2	;
a33_c20		= a_c20 * ( cos( gamma_c20 ) )^2	;
a33_p1 		= a_p1	* ( cos( gamma_p1  ) )^2	;
a33_p2 		= a_p2	* ( cos( gamma_p2  ) )^2	;
a33_p3 		= a_p3	* ( cos( gamma_p3  ) )^2	;
a33_p4 		= a_p4	* ( cos( gamma_p4  ) )^2	;
a33_b1 		= a_b1	* ( cos( gamma_b1  ) )^2	;
a33_b2 		= a_b2	* ( cos( gamma_b2  ) )^2	;
a33_b3 		= a_b3	* ( cos( gamma_b3  ) )^2	;
a33_b4 		= a_b4	* ( cos( gamma_b4  ) )^2	;
 
a33(1,1) = a33_c1 		;
a33(1,2) = a33_c2 		;
a33(1,3) = a33_c3 		;
a33(1,4) = a33_c4 		;
a33(1,5) = a33_c5 		;
a33(1,6) = a33_c6 		;
a33(1,7) = a33_c7 		;
a33(1,8) = a33_c8 		;
a33(1,9) = a33_c9 		;
a33(1,10)= a33_c10		;
a33(1,11)= a33_c11		;
a33(1,12)= a33_c12		;
a33(1,13)= a33_c13		;
a33(1,14)= a33_c14		;
a33(1,15)= a33_c15		;
a33(1,16)= a33_c16		;
a33(1,17)= a33_c17		;
a33(1,18)= a33_c18		;
a33(1,19)= a33_c19		;
a33(1,20)= a33_c20		;
a33(1,21)= a33_p1 		;
a33(1,22)= a33_p2 		;
a33(1,23)= a33_p3 		;
a33(1,24)= a33_p4 		;
a33(1,25)= a33_b1 		;
a33(1,26)= a33_b2 		;
a33(1,27)= a33_b3 		;
a33(1,28)= a33_b4 		;

a14_c1 		= -a_c1  * (c1h(1,2)  - y0)*sin( gamma_c1  )	* sin( alpha_c1  )	+a_c1  * sin( beta_c1  )	* sin( alpha_c1  )* (c1h(1,3)  - z0)	;
a14_c2 		= -a_c2  * (c2h(1,2)  - y0)*sin( gamma_c2  )	* sin( alpha_c2  )	+a_c2  * sin( beta_c2  )	* sin( alpha_c2  )* (c2h(1,3)  - z0)	;
a14_c3 		= -a_c3  * (c3h(1,2)  - y0)*sin( gamma_c3  )	* sin( alpha_c3  )	+a_c3  * sin( beta_c3  )	* sin( alpha_c3  )* (c3h(1,3)  - z0)	;
a14_c4 		= -a_c4  * (c4h(1,2)  - y0)*sin( gamma_c4  )	* sin( alpha_c4  )	+a_c4  * sin( beta_c4  )	* sin( alpha_c4  )* (c4h(1,3)  - z0)	;
a14_c5 		= -a_c5  * (c5h(1,2)  - y0)*sin( gamma_c5  )	* sin( alpha_c5  )	+a_c5  * sin( beta_c5  )	* sin( alpha_c5  )* (c5h(1,3)  - z0)	;
a14_c6 		= -a_c6  * (c6h(1,2)  - y0)*sin( gamma_c6  )	* sin( alpha_c6  )	+a_c6  * sin( beta_c6  )	* sin( alpha_c6  )* (c6h(1,3)  - z0)	;
a14_c7 		= -a_c7  * (c7h(1,2)  - y0)*sin( gamma_c7  )	* sin( alpha_c7  )	+a_c7  * sin( beta_c7  )	* sin( alpha_c7  )* (c7h(1,3)  - z0)	;
a14_c8 		= -a_c8  * (c8h(1,2)  - y0)*sin( gamma_c8  )	* sin( alpha_c8  )	+a_c8  * sin( beta_c8  )	* sin( alpha_c8  )* (c8h(1,3)  - z0)	;
a14_c9 		= -a_c9  * (c9h(1,2)  - y0)*sin( gamma_c9  )	* sin( alpha_c9  )	+a_c9  * sin( beta_c9  )	* sin( alpha_c9  )* (c9h(1,3)  - z0)	;
a14_c10		= -a_c10 * (c10h(1,2) - y0)*sin( gamma_c10 )	* sin( alpha_c10 )	+a_c10 * sin( beta_c10 )	* sin( alpha_c10 )* (c10h(1,3) - z0)	;
a14_c11		= -a_c11 * (c11h(1,2) - y0)*sin( gamma_c11 )	* sin( alpha_c11 )	+a_c11 * sin( beta_c11 )	* sin( alpha_c11 )* (c11h(1,3) - z0)	;
a14_c12		= -a_c12 * (c12h(1,2) - y0)*sin( gamma_c12 )	* sin( alpha_c12 )	+a_c12 * sin( beta_c12 )	* sin( alpha_c12 )* (c12h(1,3) - z0)	;
a14_c13		= -a_c13 * (c13h(1,2) - y0)*sin( gamma_c13 )	* sin( alpha_c13 )	+a_c13 * sin( beta_c13 )	* sin( alpha_c13 )* (c13h(1,3) - z0)	;
a14_c14		= -a_c14 * (c14h(1,2) - y0)*sin( gamma_c14 )	* sin( alpha_c14 )	+a_c14 * sin( beta_c14 )	* sin( alpha_c14 )* (c14h(1,3) - z0)	;
a14_c15		= -a_c15 * (c15h(1,2) - y0)*sin( gamma_c15 )	* sin( alpha_c15 )	+a_c15 * sin( beta_c15 )	* sin( alpha_c15 )* (c15h(1,3) - z0)	;
a14_c16		= -a_c16 * (c16h(1,2) - y0)*sin( gamma_c16 )	* sin( alpha_c16 )	+a_c16 * sin( beta_c16 )	* sin( alpha_c16 )* (c16h(1,3) - z0)	;
a14_c17		= -a_c17 * (c17h(1,2) - y0)*sin( gamma_c17 )	* sin( alpha_c17 )	+a_c17 * sin( beta_c17 )	* sin( alpha_c17 )* (c17h(1,3) - z0)	;
a14_c18		= -a_c18 * (c18h(1,2) - y0)*sin( gamma_c18 )	* sin( alpha_c18 )	+a_c18 * sin( beta_c18 )	* sin( alpha_c18 )* (c18h(1,3) - z0)	;
a14_c19		= -a_c19 * (c19h(1,2) - y0)*sin( gamma_c19 )	* sin( alpha_c19 )	+a_c19 * sin( beta_c19 )	* sin( alpha_c19 )* (c19h(1,3) - z0)	;
a14_c20		= -a_c20 * (c20h(1,2) - y0)*sin( gamma_c20 )	* sin( alpha_c20 )	+a_c20 * sin( beta_c20 )	* sin( alpha_c20 )* (c20h(1,3) - z0)	;
a14_p1 		= -a_p1	 * (p1a(1,2)  - y0)*sin( gamma_p1  )	* sin( alpha_p1  )	+a_p1  * sin( beta_p1  )	* sin( alpha_p1  )* (p1a(1,3)  - z0) 	;
a14_p2 		= -a_p2	 * (p2a(1,2)  - y0)*sin( gamma_p2  )	* sin( alpha_p2  )	+a_p2  * sin( beta_p2  )	* sin( alpha_p2  )* (p2a(1,3)  - z0) 	;
a14_p3 		= -a_p3	 * (p3a(1,2)  - y0)*sin( gamma_p3  )	* sin( alpha_p3  )	+a_p3  * sin( beta_p3  )	* sin( alpha_p3  )* (p3a(1,3)  - z0) 	;
a14_p4 		= -a_p4	 * (p4a(1,2)  - y0)*sin( gamma_p4  )	* sin( alpha_p4  )	+a_p4  * sin( beta_p4  )	* sin( alpha_p4  )* (p4a(1,3)  - z0) 	;
a14_b1 		= -a_b1	 * (b1s(1,2)  - y0)*sin( gamma_b1  )	* sin( alpha_b1  )	+a_b1  * sin( beta_b1  )	* sin( alpha_b1  )* (b1s(1,3)  - z0) 	;
a14_b2 		= -a_b2	 * (b2s(1,2)  - y0)*sin( gamma_b2  )	* sin( alpha_b2  )	+a_b2  * sin( beta_b2  )	* sin( alpha_b2  )* (b2s(1,3)  - z0) 	;
a14_b3 		= -a_b3	 * (b3s(1,2)  - y0)*sin( gamma_b3  )	* sin( alpha_b3  )	+a_b3  * sin( beta_b3  )	* sin( alpha_b3  )* (b3s(1,3)  - z0) 	;
a14_b4 		= -a_b4	 * (b4s(1,2)  - y0)*sin( gamma_b4  )	* sin( alpha_b4  )	+a_b4  * sin( beta_b4  )	* sin( alpha_b4  )* (b4s(1,3)  - z0) 	;
  
a14(1,1) = a14_c1 		;
a14(1,2) = a14_c2 		;
a14(1,3) = a14_c3 		;
a14(1,4) = a14_c4 		;
a14(1,5) = a14_c5 		;
a14(1,6) = a14_c6 		;
a14(1,7) = a14_c7 		;
a14(1,8) = a14_c8 		;
a14(1,9) = a14_c9 		;
a14(1,10)= a14_c10		;
a14(1,11)= a14_c11		;
a14(1,12)= a14_c12		;
a14(1,13)= a14_c13		;
a14(1,14)= a14_c14		;
a14(1,15)= a14_c15		;
a14(1,16)= a14_c16		;
a14(1,17)= a14_c17		;
a14(1,18)= a14_c18		;
a14(1,19)= a14_c19		;
a14(1,20)= a14_c20		;
a14(1,21)= a14_p1 		;
a14(1,22)= a14_p2 		;
a14(1,23)= a14_p3 		;
a14(1,24)= a14_p4 		;
a14(1,25)= a14_b1 		;
a14(1,26)= a14_b2 		;
a14(1,27)= a14_b3 		;
a14(1,28)= a14_b4 		;

a41_c1 		= a14_c1 		;
a41_c2 		= a14_c2 		;
a41_c3 		= a14_c3 		;
a41_c4 		= a14_c4 		;
a41_c5 		= a14_c5 		;
a41_c6 		= a14_c6 		;
a41_c7 		= a14_c7 		;
a41_c8 		= a14_c8 		;
a41_c9 		= a14_c9 		;
a41_c10		= a14_c10		;
a41_c11		= a14_c11		;
a41_c12		= a14_c12		;
a41_c13		= a14_c13		;
a41_c14		= a14_c14		;
a41_c15		= a14_c15		;
a41_c16		= a14_c16		;
a41_c17		= a14_c17		;
a41_c18		= a14_c18		;
a41_c19		= a14_c19		;
a41_c20		= a14_c20		;
a41_p1 		= a14_p1 		;
a41_p2 		= a14_p2 		;
a41_p3 		= a14_p3 		;
a41_p4 		= a14_p4 		;
a41_b1 		= a14_b1 		;
a41_b2 		= a14_b2 		;
a41_b3 		= a14_b3 		;
a41_b4 		= a14_b4 		;
 
a41(1,1) = a41_c1 		;
a41(1,2) = a41_c2 		;
a41(1,3) = a41_c3 		;
a41(1,4) = a41_c4 		;
a41(1,5) = a41_c5 		;
a41(1,6) = a41_c6 		;
a41(1,7) = a41_c7 		;
a41(1,8) = a41_c8 		;
a41(1,9) = a41_c9 		;
a41(1,10)= a41_c10		;
a41(1,11)= a41_c11		;
a41(1,12)= a41_c12		;
a41(1,13)= a41_c13		;
a41(1,14)= a41_c14		;
a41(1,15)= a41_c15		;
a41(1,16)= a41_c16		;
a41(1,17)= a41_c17		;
a41(1,18)= a41_c18		;
a41(1,19)= a41_c19		;
a41(1,20)= a41_c20		;
a41(1,21)= a41_p1 		;
a41(1,22)= a41_p2 		;
a41(1,23)= a41_p3 		;
a41(1,24)= a41_p4 		;
a41(1,25)= a41_b1 		;
a41(1,26)= a41_b2 		;
a41(1,27)= a41_b3 		;
a41(1,28)= a41_b4 		;

a15_c1 		= s_c1  * sin( gamma_c1  ) + m_c1  	* (c1h(1,1)  - x0)*sin( gamma_c1  )	* sin( alpha_c1  )	+a_c1  * ( cos( alpha_c1  ) )^2* (c1h(1,3)  - z0)	;
a15_c2 		= s_c2  * sin( gamma_c2  ) + m_c2  	* (c2h(1,1)  - x0)*sin( gamma_c2  )	* sin( alpha_c2  )	+a_c2  * ( cos( alpha_c2  ) )^2* (c2h(1,3)  - z0)	;
a15_c3 		= s_c3  * sin( gamma_c3  ) + m_c3  	* (c3h(1,1)  - x0)*sin( gamma_c3  )	* sin( alpha_c3  )	+a_c3  * ( cos( alpha_c3  ) )^2* (c3h(1,3)  - z0)	;
a15_c4 		= s_c4  * sin( gamma_c4  ) + m_c4  	* (c4h(1,1)  - x0)*sin( gamma_c4  )	* sin( alpha_c4  )	+a_c4  * ( cos( alpha_c4  ) )^2* (c4h(1,3)  - z0)	;
a15_c5 		= s_c5  * sin( gamma_c5  ) + m_c5  	* (c5h(1,1)  - x0)*sin( gamma_c5  )	* sin( alpha_c5  )	+a_c5  * ( cos( alpha_c5  ) )^2* (c5h(1,3)  - z0)	;
a15_c6 		= s_c6  * sin( gamma_c6  ) + m_c6  	* (c6h(1,1)  - x0)*sin( gamma_c6  )	* sin( alpha_c6  )	+a_c6  * ( cos( alpha_c6  ) )^2* (c6h(1,3)  - z0)	;
a15_c7 		= s_c7  * sin( gamma_c7  ) + m_c7  	* (c7h(1,1)  - x0)*sin( gamma_c7  )	* sin( alpha_c7  )	+a_c7  * ( cos( alpha_c7  ) )^2* (c7h(1,3)  - z0)	;
a15_c8 		= s_c8  * sin( gamma_c8  ) + m_c8  	* (c8h(1,1)  - x0)*sin( gamma_c8  )	* sin( alpha_c8  )	+a_c8  * ( cos( alpha_c8  ) )^2* (c8h(1,3)  - z0)	;
a15_c9 		= s_c9  * sin( gamma_c9  ) + m_c9  	* (c9h(1,1)  - x0)*sin( gamma_c9  )	* sin( alpha_c9  )	+a_c9  * ( cos( alpha_c9  ) )^2* (c9h(1,3)  - z0)	;
a15_c10		= s_c10 * sin( gamma_c10 ) + m_c10 	* (c10h(1,1) - x0)*sin( gamma_c10 )	* sin( alpha_c10 )	+a_c10 * ( cos( alpha_c10 ) )^2* (c10h(1,3) - z0)	;
a15_c11		= s_c11 * sin( gamma_c11 ) + m_c11 	* (c11h(1,1) - x0)*sin( gamma_c11 )	* sin( alpha_c11 )	+a_c11 * ( cos( alpha_c11 ) )^2* (c11h(1,3) - z0)	;
a15_c12		= s_c12 * sin( gamma_c12 ) + m_c12 	* (c12h(1,1) - x0)*sin( gamma_c12 )	* sin( alpha_c12 )	+a_c12 * ( cos( alpha_c12 ) )^2* (c12h(1,3) - z0)	;
a15_c13		= s_c13 * sin( gamma_c13 ) + m_c13 	* (c13h(1,1) - x0)*sin( gamma_c13 )	* sin( alpha_c13 )	+a_c13 * ( cos( alpha_c13 ) )^2* (c13h(1,3) - z0)	;
a15_c14		= s_c14 * sin( gamma_c14 ) + m_c14 	* (c14h(1,1) - x0)*sin( gamma_c14 )	* sin( alpha_c14 )	+a_c14 * ( cos( alpha_c14 ) )^2* (c14h(1,3) - z0)	;
a15_c15		= s_c15 * sin( gamma_c15 ) + m_c15 	* (c15h(1,1) - x0)*sin( gamma_c15 )	* sin( alpha_c15 )	+a_c15 * ( cos( alpha_c15 ) )^2* (c15h(1,3) - z0)	;
a15_c16		= s_c16 * sin( gamma_c16 ) + m_c16 	* (c16h(1,1) - x0)*sin( gamma_c16 )	* sin( alpha_c16 )	+a_c16 * ( cos( alpha_c16 ) )^2* (c16h(1,3) - z0)	;
a15_c17		= s_c17 * sin( gamma_c17 ) + m_c17 	* (c17h(1,1) - x0)*sin( gamma_c17 )	* sin( alpha_c17 )	+a_c17 * ( cos( alpha_c17 ) )^2* (c17h(1,3) - z0)	;
a15_c18		= s_c18 * sin( gamma_c18 ) + m_c18 	* (c18h(1,1) - x0)*sin( gamma_c18 )	* sin( alpha_c18 )	+a_c18 * ( cos( alpha_c18 ) )^2* (c18h(1,3) - z0)	;
a15_c19		= s_c19 * sin( gamma_c19 ) + m_c19 	* (c19h(1,1) - x0)*sin( gamma_c19 )	* sin( alpha_c19 )	+a_c19 * ( cos( alpha_c19 ) )^2* (c19h(1,3) - z0)	;
a15_c20		= s_c20 * sin( gamma_c20 ) + m_c20 	* (c20h(1,1) - x0)*sin( gamma_c20 )	* sin( alpha_c20 )	+a_c20 * ( cos( alpha_c20 ) )^2* (c20h(1,3) - z0)	;
a15_p1 		= s_p1	* sin( gamma_p1  ) + m_p1	* (p1a(1,1)  - x0)*sin( gamma_p1  )	* sin( alpha_p1  )	+a_p1  * ( cos( alpha_p1  ) )^2* (p1a(1,3)  - z0) 	;
a15_p2 		= s_p2	* sin( gamma_p2  ) + m_p2	* (p2a(1,1)  - x0)*sin( gamma_p2  )	* sin( alpha_p2  )	+a_p2  * ( cos( alpha_p2  ) )^2* (p2a(1,3)  - z0) 	;
a15_p3 		= s_p3	* sin( gamma_p3  ) + m_p3	* (p3a(1,1)  - x0)*sin( gamma_p3  )	* sin( alpha_p3  )	+a_p3  * ( cos( alpha_p3  ) )^2* (p3a(1,3)  - z0) 	;
a15_p4 		= s_p4	* sin( gamma_p4  ) + m_p4	* (p4a(1,1)  - x0)*sin( gamma_p4  )	* sin( alpha_p4  )	+a_p4  * ( cos( alpha_p4  ) )^2* (p4a(1,3)  - z0) 	;
a15_b1 		= s_b1	* sin( gamma_b1  ) + m_b1	* (b1s(1,1)  - x0)*sin( gamma_b1  )	* sin( alpha_b1  )	+a_b1  * ( cos( alpha_b1  ) )^2* (b1s(1,3)  - z0) 	;
a15_b2 		= s_b2	* sin( gamma_b2  ) + m_b2	* (b2s(1,1)  - x0)*sin( gamma_b2  )	* sin( alpha_b2  )	+a_b2  * ( cos( alpha_b2  ) )^2* (b2s(1,3)  - z0) 	;
a15_b3 		= s_b3	* sin( gamma_b3  ) + m_b3	* (b3s(1,1)  - x0)*sin( gamma_b3  )	* sin( alpha_b3  )	+a_b3  * ( cos( alpha_b3  ) )^2* (b3s(1,3)  - z0) 	;
a15_b4 		= s_b4	* sin( gamma_b4  ) + m_b4	* (b4s(1,1)  - x0)*sin( gamma_b4  )	* sin( alpha_b4  )	+a_b4  * ( cos( alpha_b4  ) )^2* (b4s(1,3)  - z0) 	;
  
a15(1,1) = a15_c1 		;
a15(1,2) = a15_c2 		;
a15(1,3) = a15_c3 		;
a15(1,4) = a15_c4 		;
a15(1,5) = a15_c5 		;
a15(1,6) = a15_c6 		;
a15(1,7) = a15_c7 		;
a15(1,8) = a15_c8 		;
a15(1,9) = a15_c9 		;
a15(1,10)= a15_c10		;
a15(1,11)= a15_c11		;
a15(1,12)= a15_c12		;
a15(1,13)= a15_c13		;
a15(1,14)= a15_c14		;
a15(1,15)= a15_c15		;
a15(1,16)= a15_c16		;
a15(1,17)= a15_c17		;
a15(1,18)= a15_c18		;
a15(1,19)= a15_c19		;
a15(1,20)= a15_c20		;
a15(1,21)= a15_p1 		;
a15(1,22)= a15_p2 		;
a15(1,23)= a15_p3 		;
a15(1,24)= a15_p4 		;
a15(1,25)= a15_b1 		;
a15(1,26)= a15_b2 		;
a15(1,27)= a15_b3 		;
a15(1,28)= a15_b4 		;

a51_c1 		= a15_c1 		;
a51_c2 		= a15_c2 		;
a51_c3 		= a15_c3 		;
a51_c4 		= a15_c4 		;
a51_c5 		= a15_c5 		;
a51_c6 		= a15_c6 		;
a51_c7 		= a15_c7 		;
a51_c8 		= a15_c8 		;
a51_c9 		= a15_c9 		;
a51_c10		= a15_c10		;
a51_c11		= a15_c11		;
a51_c12		= a15_c12		;
a51_c13		= a15_c13		;
a51_c14		= a15_c14		;
a51_c15		= a15_c15		;
a51_c16		= a15_c16		;
a51_c17		= a15_c17		;
a51_c18		= a15_c18		;
a51_c19		= a15_c19		;
a51_c20		= a15_c20		;
a51_p1 		= a15_p1 		;
a51_p2 		= a15_p2 		;
a51_p3 		= a15_p3 		;
a51_p4 		= a15_p4 		;
a51_b1 		= a15_b1 		;
a51_b2 		= a15_b2 		;
a51_b3 		= a15_b3 		;
a51_b4 		= a15_b4 		;
 
a51(1,1) = a51_c1 		;
a51(1,2) = a51_c2 		;
a51(1,3) = a51_c3 		;
a51(1,4) = a51_c4 		;
a51(1,5) = a51_c5 		;
a51(1,6) = a51_c6 		;
a51(1,7) = a51_c7 		;
a51(1,8) = a51_c8 		;
a51(1,9) = a51_c9 		;
a51(1,10)= a51_c10		;
a51(1,11)= a51_c11		;
a51(1,12)= a51_c12		;
a51(1,13)= a51_c13		;
a51(1,14)= a51_c14		;
a51(1,15)= a51_c15		;
a51(1,16)= a51_c16		;
a51(1,17)= a51_c17		;
a51(1,18)= a51_c18		;
a51(1,19)= a51_c19		;
a51(1,20)= a51_c20		;
a51(1,21)= a51_p1 		;
a51(1,22)= a51_p2 		;
a51(1,23)= a51_p3 		;
a51(1,24)= a51_p4 		;
a51(1,25)= a51_b1 		;
a51(1,26)= a51_b2 		;
a51(1,27)= a51_b3 		;
a51(1,28)= a51_b4 		;

a16_c1 		= -s_c1  * sin( beta_c1  ) - m_c1  	* (c1h(1,1)  - x0)*sin( beta_c1  )	* sin( alpha_c1  )	-a_c1  * ( cos( alpha_c1  ) )^2* (c1h(1,2)  - y0)	;
a16_c2 		= -s_c2  * sin( beta_c2  ) - m_c2  	* (c2h(1,1)  - x0)*sin( beta_c2  )	* sin( alpha_c2  )	-a_c2  * ( cos( alpha_c2  ) )^2* (c2h(1,2)  - y0)	;
a16_c3 		= -s_c3  * sin( beta_c3  ) - m_c3  	* (c3h(1,1)  - x0)*sin( beta_c3  )	* sin( alpha_c3  )	-a_c3  * ( cos( alpha_c3  ) )^2* (c3h(1,2)  - y0)	;
a16_c4 		= -s_c4  * sin( beta_c4  ) - m_c4  	* (c4h(1,1)  - x0)*sin( beta_c4  )	* sin( alpha_c4  )	-a_c4  * ( cos( alpha_c4  ) )^2* (c4h(1,2)  - y0)	;
a16_c5 		= -s_c5  * sin( beta_c5  ) - m_c5  	* (c5h(1,1)  - x0)*sin( beta_c5  )	* sin( alpha_c5  )	-a_c5  * ( cos( alpha_c5  ) )^2* (c5h(1,2)  - y0)	;
a16_c6 		= -s_c6  * sin( beta_c6  ) - m_c6  	* (c6h(1,1)  - x0)*sin( beta_c6  )	* sin( alpha_c6  )	-a_c6  * ( cos( alpha_c6  ) )^2* (c6h(1,2)  - y0)	;
a16_c7 		= -s_c7  * sin( beta_c7  ) - m_c7  	* (c7h(1,1)  - x0)*sin( beta_c7  )	* sin( alpha_c7  )	-a_c7  * ( cos( alpha_c7  ) )^2* (c7h(1,2)  - y0)	;
a16_c8 		= -s_c8  * sin( beta_c8  ) - m_c8  	* (c8h(1,1)  - x0)*sin( beta_c8  )	* sin( alpha_c8  )	-a_c8  * ( cos( alpha_c8  ) )^2* (c8h(1,2)  - y0)	;
a16_c9 		= -s_c9  * sin( beta_c9  ) - m_c9  	* (c9h(1,1)  - x0)*sin( beta_c9  )	* sin( alpha_c9  )	-a_c9  * ( cos( alpha_c9  ) )^2* (c9h(1,2)  - y0)	;
a16_c10		= -s_c10 * sin( beta_c10 ) - m_c10 	* (c10h(1,1) - x0)*sin( beta_c10 )	* sin( alpha_c10 )	-a_c10 * ( cos( alpha_c10 ) )^2* (c10h(1,2) - y0)	;
a16_c11		= -s_c11 * sin( beta_c11 ) - m_c11 	* (c11h(1,1) - x0)*sin( beta_c11 )	* sin( alpha_c11 )	-a_c11 * ( cos( alpha_c11 ) )^2* (c11h(1,2) - y0)	;
a16_c12		= -s_c12 * sin( beta_c12 ) - m_c12 	* (c12h(1,1) - x0)*sin( beta_c12 )	* sin( alpha_c12 )	-a_c12 * ( cos( alpha_c12 ) )^2* (c12h(1,2) - y0)	;
a16_c13		= -s_c13 * sin( beta_c13 ) - m_c13 	* (c13h(1,1) - x0)*sin( beta_c13 )	* sin( alpha_c13 )	-a_c13 * ( cos( alpha_c13 ) )^2* (c13h(1,2) - y0)	;
a16_c14		= -s_c14 * sin( beta_c14 ) - m_c14 	* (c14h(1,1) - x0)*sin( beta_c14 )	* sin( alpha_c14 )	-a_c14 * ( cos( alpha_c14 ) )^2* (c14h(1,2) - y0)	;
a16_c15		= -s_c15 * sin( beta_c15 ) - m_c15 	* (c15h(1,1) - x0)*sin( beta_c15 )	* sin( alpha_c15 )	-a_c15 * ( cos( alpha_c15 ) )^2* (c15h(1,2) - y0)	;
a16_c16		= -s_c16 * sin( beta_c16 ) - m_c16 	* (c16h(1,1) - x0)*sin( beta_c16 )	* sin( alpha_c16 )	-a_c16 * ( cos( alpha_c16 ) )^2* (c16h(1,2) - y0)	;
a16_c17		= -s_c17 * sin( beta_c17 ) - m_c17 	* (c17h(1,1) - x0)*sin( beta_c17 )	* sin( alpha_c17 )	-a_c17 * ( cos( alpha_c17 ) )^2* (c17h(1,2) - y0)	;
a16_c18		= -s_c18 * sin( beta_c18 ) - m_c18 	* (c18h(1,1) - x0)*sin( beta_c18 )	* sin( alpha_c18 )	-a_c18 * ( cos( alpha_c18 ) )^2* (c18h(1,2) - y0)	;
a16_c19		= -s_c19 * sin( beta_c19 ) - m_c19 	* (c19h(1,1) - x0)*sin( beta_c19 )	* sin( alpha_c19 )	-a_c19 * ( cos( alpha_c19 ) )^2* (c19h(1,2) - y0)	;
a16_c20		= -s_c20 * sin( beta_c20 ) - m_c20 	* (c20h(1,1) - x0)*sin( beta_c20 )	* sin( alpha_c20 )	-a_c20 * ( cos( alpha_c20 ) )^2* (c20h(1,2) - y0)	;
a16_p1 		= -s_p1	*  sin( beta_p1  ) - m_p1	* (p1a(1,1)  - x0)*sin( beta_p1  )	* sin( alpha_p1  )	-a_p1  * ( cos( alpha_p1  ) )^2* (p1a(1,2)  - y0) 	;
a16_p2 		= -s_p2	*  sin( beta_p2  ) - m_p2	* (p2a(1,1)  - x0)*sin( beta_p2  )	* sin( alpha_p2  )	-a_p2  * ( cos( alpha_p2  ) )^2* (p2a(1,2)  - y0) 	;
a16_p3 		= -s_p3	*  sin( beta_p3  ) - m_p3	* (p3a(1,1)  - x0)*sin( beta_p3  )	* sin( alpha_p3  )	-a_p3  * ( cos( alpha_p3  ) )^2* (p3a(1,2)  - y0) 	;
a16_p4 		= -s_p4	*  sin( beta_p4  ) - m_p4	* (p4a(1,1)  - x0)*sin( beta_p4  )	* sin( alpha_p4  )	-a_p4  * ( cos( alpha_p4  ) )^2* (p4a(1,2)  - y0) 	;
a16_b1 		= -s_b1	*  sin( beta_b1  ) - m_b1	* (b1s(1,1)  - x0)*sin( beta_b1  )	* sin( alpha_b1  )	-a_b1  * ( cos( alpha_b1  ) )^2* (b1s(1,2)  - y0) 	;
a16_b2 		= -s_b2	*  sin( beta_b2  ) - m_b2	* (b2s(1,1)  - x0)*sin( beta_b2  )	* sin( alpha_b2  )	-a_b2  * ( cos( alpha_b2  ) )^2* (b2s(1,2)  - y0) 	;
a16_b3 		= -s_b3	*  sin( beta_b3  ) - m_b3	* (b3s(1,1)  - x0)*sin( beta_b3  )	* sin( alpha_b3  )	-a_b3  * ( cos( alpha_b3  ) )^2* (b3s(1,2)  - y0) 	;
a16_b4 		= -s_b4	*  sin( beta_b4  ) - m_b4	* (b4s(1,1)  - x0)*sin( beta_b4  )	* sin( alpha_b4  )	-a_b4  * ( cos( alpha_b4  ) )^2* (b4s(1,2)  - y0) 	;
  
a16(1,1) = a16_c1 		;
a16(1,2) = a16_c2 		;
a16(1,3) = a16_c3 		;
a16(1,4) = a16_c4 		;
a16(1,5) = a16_c5 		;
a16(1,6) = a16_c6 		;
a16(1,7) = a16_c7 		;
a16(1,8) = a16_c8 		;
a16(1,9) = a16_c9 		;
a16(1,10)= a16_c10		;
a16(1,11)= a16_c11		;
a16(1,12)= a16_c12		;
a16(1,13)= a16_c13		;
a16(1,14)= a16_c14		;
a16(1,15)= a16_c15		;
a16(1,16)= a16_c16		;
a16(1,17)= a16_c17		;
a16(1,18)= a16_c18		;
a16(1,19)= a16_c19		;
a16(1,20)= a16_c20		;
a16(1,21)= a16_p1 		;
a16(1,22)= a16_p2 		;
a16(1,23)= a16_p3 		;
a16(1,24)= a16_p4 		;
a16(1,25)= a16_b1 		;
a16(1,26)= a16_b2 		;
a16(1,27)= a16_b3 		;
a16(1,28)= a16_b4 		;

a61_c1 		= a16_c1 		;
a61_c2 		= a16_c2 		;
a61_c3 		= a16_c3 		;
a61_c4 		= a16_c4 		;
a61_c5 		= a16_c5 		;
a61_c6 		= a16_c6 		;
a61_c7 		= a16_c7 		;
a61_c8 		= a16_c8 		;
a61_c9 		= a16_c9 		;
a61_c10		= a16_c10		;
a61_c11		= a16_c11		;
a61_c12		= a16_c12		;
a61_c13		= a16_c13		;
a61_c14		= a16_c14		;
a61_c15		= a16_c15		;
a61_c16		= a16_c16		;
a61_c17		= a16_c17		;
a61_c18		= a16_c18		;
a61_c19		= a16_c19		;
a61_c20		= a16_c20		;
a61_p1 		= a16_p1 		;
a61_p2 		= a16_p2 		;
a61_p3 		= a16_p3 		;
a61_p4 		= a16_p4 		;
a61_b1 		= a16_b1 		;
a61_b2 		= a16_b2 		;
a61_b3 		= a16_b3 		;
a61_b4 		= a16_b4 		;
 
a61(1,1) = a61_c1 		;
a61(1,2) = a61_c2 		;
a61(1,3) = a61_c3 		;
a61(1,4) = a61_c4 		;
a61(1,5) = a61_c5 		;
a61(1,6) = a61_c6 		;
a61(1,7) = a61_c7 		;
a61(1,8) = a61_c8 		;
a61(1,9) = a61_c9 		;
a61(1,10)= a61_c10		;
a61(1,11)= a61_c11		;
a61(1,12)= a61_c12		;
a61(1,13)= a61_c13		;
a61(1,14)= a61_c14		;
a61(1,15)= a61_c15		;
a61(1,16)= a61_c16		;
a61(1,17)= a61_c17		;
a61(1,18)= a61_c18		;
a61(1,19)= a61_c19		;
a61(1,20)= a61_c20		;
a61(1,21)= a61_p1 		;
a61(1,22)= a61_p2 		;
a61(1,23)= a61_p3 		;
a61(1,24)= a61_p4 		;
a61(1,25)= a61_b1 		;
a61(1,26)= a61_b2 		;
a61(1,27)= a61_b3 		;
a61(1,28)= a61_b4 		;

a24_c1 		= -s_c1  * sin( gamma_c1  ) - m_c1  * (c1h(1,3)  - z0)*( cos( beta_c1  ) )^2 -a_c1  * sin( gamma_c1  )	* sin( beta_c1  ) *(c1h(1,2)  - y0)	;
a24_c2 		= -s_c2  * sin( gamma_c2  ) - m_c2  * (c2h(1,3)  - z0)*( cos( beta_c2  ) )^2 -a_c2  * sin( gamma_c2  )	* sin( beta_c2  ) *(c2h(1,2)  - y0)	;
a24_c3 		= -s_c3  * sin( gamma_c3  ) - m_c3  * (c3h(1,3)  - z0)*( cos( beta_c3  ) )^2 -a_c3  * sin( gamma_c3  )	* sin( beta_c3  ) *(c3h(1,2)  - y0)	;
a24_c4 		= -s_c4  * sin( gamma_c4  ) - m_c4  * (c4h(1,3)  - z0)*( cos( beta_c4  ) )^2 -a_c4  * sin( gamma_c4  )	* sin( beta_c4  ) *(c4h(1,2)  - y0)	;
a24_c5 		= -s_c5  * sin( gamma_c5  ) - m_c5  * (c5h(1,3)  - z0)*( cos( beta_c5  ) )^2 -a_c5  * sin( gamma_c5  )	* sin( beta_c5  ) *(c5h(1,2)  - y0)	;
a24_c6 		= -s_c6  * sin( gamma_c6  ) - m_c6  * (c6h(1,3)  - z0)*( cos( beta_c6  ) )^2 -a_c6  * sin( gamma_c6  )	* sin( beta_c6  ) *(c6h(1,2)  - y0)	;
a24_c7 		= -s_c7  * sin( gamma_c7  ) - m_c7  * (c7h(1,3)  - z0)*( cos( beta_c7  ) )^2 -a_c7  * sin( gamma_c7  )	* sin( beta_c7  ) *(c7h(1,2)  - y0)	;
a24_c8 		= -s_c8  * sin( gamma_c8  ) - m_c8  * (c8h(1,3)  - z0)*( cos( beta_c8  ) )^2 -a_c8  * sin( gamma_c8  )	* sin( beta_c8  ) *(c8h(1,2)  - y0)	;
a24_c9 		= -s_c9  * sin( gamma_c9  ) - m_c9  * (c9h(1,3)  - z0)*( cos( beta_c9  ) )^2 -a_c9  * sin( gamma_c9  )	* sin( beta_c9  ) *(c9h(1,2)  - y0)	;
a24_c10		= -s_c10 * sin( gamma_c10 ) - m_c10 * (c10h(1,3) - z0)*( cos( beta_c10 ) )^2 -a_c10 * sin( gamma_c10 )	* sin( beta_c10 ) *(c10h(1,2) - y0)	;
a24_c11		= -s_c11 * sin( gamma_c11 ) - m_c11 * (c11h(1,3) - z0)*( cos( beta_c11 ) )^2 -a_c11 * sin( gamma_c11 )	* sin( beta_c11 ) *(c11h(1,2) - y0)	;
a24_c12		= -s_c12 * sin( gamma_c12 ) - m_c12 * (c12h(1,3) - z0)*( cos( beta_c12 ) )^2 -a_c12 * sin( gamma_c12 )	* sin( beta_c12 ) *(c12h(1,2) - y0)	;
a24_c13		= -s_c13 * sin( gamma_c13 ) - m_c13 * (c13h(1,3) - z0)*( cos( beta_c13 ) )^2 -a_c13 * sin( gamma_c13 )	* sin( beta_c13 ) *(c13h(1,2) - y0)	;
a24_c14		= -s_c14 * sin( gamma_c14 ) - m_c14 * (c14h(1,3) - z0)*( cos( beta_c14 ) )^2 -a_c14 * sin( gamma_c14 )	* sin( beta_c14 ) *(c14h(1,2) - y0)	;
a24_c15		= -s_c15 * sin( gamma_c15 ) - m_c15 * (c15h(1,3) - z0)*( cos( beta_c15 ) )^2 -a_c15 * sin( gamma_c15 )	* sin( beta_c15 ) *(c15h(1,2) - y0)	;
a24_c16		= -s_c16 * sin( gamma_c16 ) - m_c16 * (c16h(1,3) - z0)*( cos( beta_c16 ) )^2 -a_c16 * sin( gamma_c16 )	* sin( beta_c16 ) *(c16h(1,2) - y0)	;
a24_c17		= -s_c17 * sin( gamma_c17 ) - m_c17 * (c17h(1,3) - z0)*( cos( beta_c17 ) )^2 -a_c17 * sin( gamma_c17 )	* sin( beta_c17 ) *(c17h(1,2) - y0)	;
a24_c18		= -s_c18 * sin( gamma_c18 ) - m_c18 * (c18h(1,3) - z0)*( cos( beta_c18 ) )^2 -a_c18 * sin( gamma_c18 )	* sin( beta_c18 ) *(c18h(1,2) - y0)	;
a24_c19		= -s_c19 * sin( gamma_c19 ) - m_c19 * (c19h(1,3) - z0)*( cos( beta_c19 ) )^2 -a_c19 * sin( gamma_c19 )	* sin( beta_c19 ) *(c19h(1,2) - y0)	;
a24_c20		= -s_c20 * sin( gamma_c20 ) - m_c20 * (c20h(1,3) - z0)*( cos( beta_c20 ) )^2 -a_c20 * sin( gamma_c20 )	* sin( beta_c20 ) *(c20h(1,2) - y0)	;
a24_p1 		= -s_p1	*  sin( gamma_p1  ) - m_p1	* (p1a(1,3)  - z0)*( cos( beta_p1  ) )^2 -a_p1  * sin( gamma_p1  )	* sin( beta_p1  ) *(p1a(1,2)  - y0) ;
a24_p2 		= -s_p2	*  sin( gamma_p2  ) - m_p2	* (p2a(1,3)  - z0)*( cos( beta_p2  ) )^2 -a_p2  * sin( gamma_p2  )	* sin( beta_p2  ) *(p2a(1,2)  - y0) ;
a24_p3 		= -s_p3	*  sin( gamma_p3  ) - m_p3	* (p3a(1,3)  - z0)*( cos( beta_p3  ) )^2 -a_p3  * sin( gamma_p3  )	* sin( beta_p3  ) *(p3a(1,2)  - y0) ;
a24_p4 		= -s_p4	*  sin( gamma_p4  ) - m_p4	* (p4a(1,3)  - z0)*( cos( beta_p4  ) )^2 -a_p4  * sin( gamma_p4  )	* sin( beta_p4  ) *(p4a(1,2)  - y0) ;
a24_b1 		= -s_b1	*  sin( gamma_b1  ) - m_b1	* (b1s(1,3)  - z0)*( cos( beta_b1  ) )^2 -a_b1  * sin( gamma_b1  )	* sin( beta_b1  ) *(b1s(1,2)  - y0) ;
a24_b2 		= -s_b2	*  sin( gamma_b2  ) - m_b2	* (b2s(1,3)  - z0)*( cos( beta_b2  ) )^2 -a_b2  * sin( gamma_b2  )	* sin( beta_b2  ) *(b2s(1,2)  - y0) ;
a24_b3 		= -s_b3	*  sin( gamma_b3  ) - m_b3	* (b3s(1,3)  - z0)*( cos( beta_b3  ) )^2 -a_b3  * sin( gamma_b3  )	* sin( beta_b3  ) *(b3s(1,2)  - y0) ;
a24_b4 		= -s_b4	*  sin( gamma_b4  ) - m_b4	* (b4s(1,3)  - z0)*( cos( beta_b4  ) )^2 -a_b4  * sin( gamma_b4  )	* sin( beta_b4  ) *(b4s(1,2)  - y0) ;
  
a24(1,1) = a24_c1 		;
a24(1,2) = a24_c2 		;
a24(1,3) = a24_c3 		;
a24(1,4) = a24_c4 		;
a24(1,5) = a24_c5 		;
a24(1,6) = a24_c6 		;
a24(1,7) = a24_c7 		;
a24(1,8) = a24_c8 		;
a24(1,9) = a24_c9 		;
a24(1,10)= a24_c10		;
a24(1,11)= a24_c11		;
a24(1,12)= a24_c12		;
a24(1,13)= a24_c13		;
a24(1,14)= a24_c14		;
a24(1,15)= a24_c15		;
a24(1,16)= a24_c16		;
a24(1,17)= a24_c17		;
a24(1,18)= a24_c18		;
a24(1,19)= a24_c19		;
a24(1,20)= a24_c20		;
a24(1,21)= a24_p1 		;
a24(1,22)= a24_p2 		;
a24(1,23)= a24_p3 		;
a24(1,24)= a24_p4 		;
a24(1,25)= a24_b1 		;
a24(1,26)= a24_b2 		;
a24(1,27)= a24_b3 		;
a24(1,28)= a24_b4 		;

a42_c1 		= a24_c1 		;
a42_c2 		= a24_c2 		;
a42_c3 		= a24_c3 		;
a42_c4 		= a24_c4 		;
a42_c5 		= a24_c5 		;
a42_c6 		= a24_c6 		;
a42_c7 		= a24_c7 		;
a42_c8 		= a24_c8 		;
a42_c9 		= a24_c9 		;
a42_c10		= a24_c10		;
a42_c11		= a24_c11		;
a42_c12		= a24_c12		;
a42_c13		= a24_c13		;
a42_c14		= a24_c14		;
a42_c15		= a24_c15		;
a42_c16		= a24_c16		;
a42_c17		= a24_c17		;
a42_c18		= a24_c18		;
a42_c19		= a24_c19		;
a42_c20		= a24_c20		;
a42_p1 		= a24_p1 		;
a42_p2 		= a24_p2 		;
a42_p3 		= a24_p3 		;
a42_p4 		= a24_p4 		;
a42_b1 		= a24_b1 		;
a42_b2 		= a24_b2 		;
a42_b3 		= a24_b3 		;
a42_b4 		= a24_b4 		;
 
a42(1,1) = a42_c1 		;
a42(1,2) = a42_c2 		;
a42(1,3) = a42_c3 		;
a42(1,4) = a42_c4 		;
a42(1,5) = a42_c5 		;
a42(1,6) = a42_c6 		;
a42(1,7) = a42_c7 		;
a42(1,8) = a42_c8 		;
a42(1,9) = a42_c9 		;
a42(1,10)= a42_c10		;
a42(1,11)= a42_c11		;
a42(1,12)= a42_c12		;
a42(1,13)= a42_c13		;
a42(1,14)= a42_c14		;
a42(1,15)= a42_c15		;
a42(1,16)= a42_c16		;
a42(1,17)= a42_c17		;
a42(1,18)= a42_c18		;
a42(1,19)= a42_c19		;
a42(1,20)= a42_c20		;
a42(1,21)= a42_p1 		;
a42(1,22)= a42_p2 		;
a42(1,23)= a42_p3 		;
a42(1,24)= a42_p4 		;
a42(1,25)= a42_b1 		;
a42(1,26)= a42_b2 		;
a42(1,27)= a42_b3 		;
a42(1,28)= a42_b4 		;

a25_c1 		= -a_c1  * sin( beta_c1  )	* sin( alpha_c1  )* (c1h(1,3)  - z0) + a_c1  * (c1h(1,1)  - x0) * sin( gamma_c1  )	* sin( beta_c1  )	;
a25_c2 		= -a_c2  * sin( beta_c2  )	* sin( alpha_c2  )* (c2h(1,3)  - z0) + a_c2  * (c2h(1,1)  - x0) * sin( gamma_c2  )	* sin( beta_c2  )	;
a25_c3 		= -a_c3  * sin( beta_c3  )	* sin( alpha_c3  )* (c3h(1,3)  - z0) + a_c3  * (c3h(1,1)  - x0) * sin( gamma_c3  )	* sin( beta_c3  )	;
a25_c4 		= -a_c4  * sin( beta_c4  )	* sin( alpha_c4  )* (c4h(1,3)  - z0) + a_c4  * (c4h(1,1)  - x0) * sin( gamma_c4  )	* sin( beta_c4  )	;
a25_c5 		= -a_c5  * sin( beta_c5  )	* sin( alpha_c5  )* (c5h(1,3)  - z0) + a_c5  * (c5h(1,1)  - x0) * sin( gamma_c5  )	* sin( beta_c5  )	;
a25_c6 		= -a_c6  * sin( beta_c6  )	* sin( alpha_c6  )* (c6h(1,3)  - z0) + a_c6  * (c6h(1,1)  - x0) * sin( gamma_c6  )	* sin( beta_c6  )	;
a25_c7 		= -a_c7  * sin( beta_c7  )	* sin( alpha_c7  )* (c7h(1,3)  - z0) + a_c7  * (c7h(1,1)  - x0) * sin( gamma_c7  )	* sin( beta_c7  )	;
a25_c8 		= -a_c8  * sin( beta_c8  )	* sin( alpha_c8  )* (c8h(1,3)  - z0) + a_c8  * (c8h(1,1)  - x0) * sin( gamma_c8  )	* sin( beta_c8  )	;
a25_c9 		= -a_c9  * sin( beta_c9  )	* sin( alpha_c9  )* (c9h(1,3)  - z0) + a_c9  * (c9h(1,1)  - x0) * sin( gamma_c9  )	* sin( beta_c9  )	;
a25_c10		= -a_c10 * sin( beta_c10 )	* sin( alpha_c10 )* (c10h(1,3) - z0) + a_c10 * (c10h(1,1) - x0) * sin( gamma_c10 )	* sin( beta_c10 )	;
a25_c11		= -a_c11 * sin( beta_c11 )	* sin( alpha_c11 )* (c11h(1,3) - z0) + a_c11 * (c11h(1,1) - x0) * sin( gamma_c11 )	* sin( beta_c11 )	;
a25_c12		= -a_c12 * sin( beta_c12 )	* sin( alpha_c12 )* (c12h(1,3) - z0) + a_c12 * (c12h(1,1) - x0) * sin( gamma_c12 )	* sin( beta_c12 )	;
a25_c13		= -a_c13 * sin( beta_c13 )	* sin( alpha_c13 )* (c13h(1,3) - z0) + a_c13 * (c13h(1,1) - x0) * sin( gamma_c13 )	* sin( beta_c13 )	;
a25_c14		= -a_c14 * sin( beta_c14 )	* sin( alpha_c14 )* (c14h(1,3) - z0) + a_c14 * (c14h(1,1) - x0) * sin( gamma_c14 )	* sin( beta_c14 )	;
a25_c15		= -a_c15 * sin( beta_c15 )	* sin( alpha_c15 )* (c15h(1,3) - z0) + a_c15 * (c15h(1,1) - x0) * sin( gamma_c15 )	* sin( beta_c15 )	;
a25_c16		= -a_c16 * sin( beta_c16 )	* sin( alpha_c16 )* (c16h(1,3) - z0) + a_c16 * (c16h(1,1) - x0) * sin( gamma_c16 )	* sin( beta_c16 )	;
a25_c17		= -a_c17 * sin( beta_c17 )	* sin( alpha_c17 )* (c17h(1,3) - z0) + a_c17 * (c17h(1,1) - x0) * sin( gamma_c17 )	* sin( beta_c17 )	;
a25_c18		= -a_c18 * sin( beta_c18 )	* sin( alpha_c18 )* (c18h(1,3) - z0) + a_c18 * (c18h(1,1) - x0) * sin( gamma_c18 )	* sin( beta_c18 )	;
a25_c19		= -a_c19 * sin( beta_c19 )	* sin( alpha_c19 )* (c19h(1,3) - z0) + a_c19 * (c19h(1,1) - x0) * sin( gamma_c19 )	* sin( beta_c19 )	;
a25_c20		= -a_c20 * sin( beta_c20 )	* sin( alpha_c20 )* (c20h(1,3) - z0) + a_c20 * (c20h(1,1) - x0) * sin( gamma_c20 )	* sin( beta_c20 )	;
a25_p1 		= -a_p1  * sin( beta_p1  )	* sin( alpha_p1  )* (p1a(1,3)  - z0) + a_p1	 * (p1a(1,1)  - x0) * sin( gamma_p1  )	* sin( beta_p1  )	;
a25_p2 		= -a_p2  * sin( beta_p2  )	* sin( alpha_p2  )* (p2a(1,3)  - z0) + a_p2	 * (p2a(1,1)  - x0) * sin( gamma_p2  )	* sin( beta_p2  )	;
a25_p3 		= -a_p3  * sin( beta_p3  )	* sin( alpha_p3  )* (p3a(1,3)  - z0) + a_p3	 * (p3a(1,1)  - x0) * sin( gamma_p3  )	* sin( beta_p3  )	;
a25_p4 		= -a_p4  * sin( beta_p4  )	* sin( alpha_p4  )* (p4a(1,3)  - z0) + a_p4	 * (p4a(1,1)  - x0) * sin( gamma_p4  )	* sin( beta_p4  )	;
a25_b1 		= -a_b1  * sin( beta_b1  )	* sin( alpha_b1  )* (b1s(1,3)  - z0) + a_b1	 * (b1s(1,1)  - x0) * sin( gamma_b1  )	* sin( beta_b1  )	;
a25_b2 		= -a_b2  * sin( beta_b2  )	* sin( alpha_b2  )* (b2s(1,3)  - z0) + a_b2	 * (b2s(1,1)  - x0) * sin( gamma_b2  )	* sin( beta_b2  )	;
a25_b3 		= -a_b3  * sin( beta_b3  )	* sin( alpha_b3  )* (b3s(1,3)  - z0) + a_b3	 * (b3s(1,1)  - x0) * sin( gamma_b3  )	* sin( beta_b3  )	;
a25_b4 		= -a_b4  * sin( beta_b4  )	* sin( alpha_b4  )* (b4s(1,3)  - z0) + a_b4	 * (b4s(1,1)  - x0) * sin( gamma_b4  )	* sin( beta_b4  )	;
  
a25(1,1) = a25_c1 		;
a25(1,2) = a25_c2 		;
a25(1,3) = a25_c3 		;
a25(1,4) = a25_c4 		;
a25(1,5) = a25_c5 		;
a25(1,6) = a25_c6 		;
a25(1,7) = a25_c7 		;
a25(1,8) = a25_c8 		;
a25(1,9) = a25_c9 		;
a25(1,10)= a25_c10		;
a25(1,11)= a25_c11		;
a25(1,12)= a25_c12		;
a25(1,13)= a25_c13		;
a25(1,14)= a25_c14		;
a25(1,15)= a25_c15		;
a25(1,16)= a25_c16		;
a25(1,17)= a25_c17		;
a25(1,18)= a25_c18		;
a25(1,19)= a25_c19		;
a25(1,20)= a25_c20		;
a25(1,21)= a25_p1 		;
a25(1,22)= a25_p2 		;
a25(1,23)= a25_p3 		;
a25(1,24)= a25_p4 		;
a25(1,25)= a25_b1 		;
a25(1,26)= a25_b2 		;
a25(1,27)= a25_b3 		;
a25(1,28)= a25_b4 		;

a52_c1 		= a25_c1 		;
a52_c2 		= a25_c2 		;
a52_c3 		= a25_c3 		;
a52_c4 		= a25_c4 		;
a52_c5 		= a25_c5 		;
a52_c6 		= a25_c6 		;
a52_c7 		= a25_c7 		;
a52_c8 		= a25_c8 		;
a52_c9 		= a25_c9 		;
a52_c10		= a25_c10		;
a52_c11		= a25_c11		;
a52_c12		= a25_c12		;
a52_c13		= a25_c13		;
a52_c14		= a25_c14		;
a52_c15		= a25_c15		;
a52_c16		= a25_c16		;
a52_c17		= a25_c17		;
a52_c18		= a25_c18		;
a52_c19		= a25_c19		;
a52_c20		= a25_c20		;
a52_p1 		= a25_p1 		;
a52_p2 		= a25_p2 		;
a52_p3 		= a25_p3 		;
a52_p4 		= a25_p4 		;
a52_b1 		= a25_b1 		;
a52_b2 		= a25_b2 		;
a52_b3 		= a25_b3 		;
a52_b4 		= a25_b4 		;
 
a52(1,1) = a52_c1 		;
a52(1,2) = a52_c2 		;
a52(1,3) = a52_c3 		;
a52(1,4) = a52_c4 		;
a52(1,5) = a52_c5 		;
a52(1,6) = a52_c6 		;
a52(1,7) = a52_c7 		;
a52(1,8) = a52_c8 		;
a52(1,9) = a52_c9 		;
a52(1,10)= a52_c10		;
a52(1,11)= a52_c11		;
a52(1,12)= a52_c12		;
a52(1,13)= a52_c13		;
a52(1,14)= a52_c14		;
a52(1,15)= a52_c15		;
a52(1,16)= a52_c16		;
a52(1,17)= a52_c17		;
a52(1,18)= a52_c18		;
a52(1,19)= a52_c19		;
a52(1,20)= a52_c20		;
a52(1,21)= a52_p1 		;
a52(1,22)= a52_p2 		;
a52(1,23)= a52_p3 		;
a52(1,24)= a52_p4 		;
a52(1,25)= a52_b1 		;
a52(1,26)= a52_b2 		;
a52(1,27)= a52_b3 		;
a52(1,28)= a52_b4 		;

a26_c1 		= s_c1  * sin( alpha_c1  ) - m_c1  	*(c1h(1,2)  - y0)* sin( alpha_c1  )	* sin( beta_c1  )+a_c1  * ( cos( beta_c1  ) )^2	 *(c1h(1,1)  - x0)	;
a26_c2 		= s_c2  * sin( alpha_c2  ) - m_c2  	*(c2h(1,2)  - y0)* sin( alpha_c2  )	* sin( beta_c2  )+a_c2  * ( cos( beta_c2  ) )^2	 *(c2h(1,1)  - x0)	;
a26_c3 		= s_c3  * sin( alpha_c3  ) - m_c3  	*(c3h(1,2)  - y0)* sin( alpha_c3  )	* sin( beta_c3  )+a_c3  * ( cos( beta_c3  ) )^2	 *(c3h(1,1)  - x0)	;
a26_c4 		= s_c4  * sin( alpha_c4  ) - m_c4  	*(c4h(1,2)  - y0)* sin( alpha_c4  )	* sin( beta_c4  )+a_c4  * ( cos( beta_c4  ) )^2	 *(c4h(1,1)  - x0)	;
a26_c5 		= s_c5  * sin( alpha_c5  ) - m_c5  	*(c5h(1,2)  - y0)* sin( alpha_c5  )	* sin( beta_c5  )+a_c5  * ( cos( beta_c5  ) )^2	 *(c5h(1,1)  - x0)	;
a26_c6 		= s_c6  * sin( alpha_c6  ) - m_c6  	*(c6h(1,2)  - y0)* sin( alpha_c6  )	* sin( beta_c6  )+a_c6  * ( cos( beta_c6  ) )^2	 *(c6h(1,1)  - x0)	;
a26_c7 		= s_c7  * sin( alpha_c7  ) - m_c7  	*(c7h(1,2)  - y0)* sin( alpha_c7  )	* sin( beta_c7  )+a_c7  * ( cos( beta_c7  ) )^2	 *(c7h(1,1)  - x0)	;
a26_c8 		= s_c8  * sin( alpha_c8  ) - m_c8  	*(c8h(1,2)  - y0)* sin( alpha_c8  )	* sin( beta_c8  )+a_c8  * ( cos( beta_c8  ) )^2	 *(c8h(1,1)  - x0)	;
a26_c9 		= s_c9  * sin( alpha_c9  ) - m_c9  	*(c9h(1,2)  - y0)* sin( alpha_c9  )	* sin( beta_c9  )+a_c9  * ( cos( beta_c9  ) )^2	 *(c9h(1,1)  - x0)	;
a26_c10		= s_c10 * sin( alpha_c10 ) - m_c10 	*(c10h(1,2) - y0)* sin( alpha_c10 )	* sin( beta_c10 )+a_c10 * ( cos( beta_c10 ) )^2	 *(c10h(1,1) - x0)	;
a26_c11		= s_c11 * sin( alpha_c11 ) - m_c11 	*(c11h(1,2) - y0)* sin( alpha_c11 )	* sin( beta_c11 )+a_c11 * ( cos( beta_c11 ) )^2	 *(c11h(1,1) - x0)	;
a26_c12		= s_c12 * sin( alpha_c12 ) - m_c12 	*(c12h(1,2) - y0)* sin( alpha_c12 )	* sin( beta_c12 )+a_c12 * ( cos( beta_c12 ) )^2	 *(c12h(1,1) - x0)	;
a26_c13		= s_c13 * sin( alpha_c13 ) - m_c13 	*(c13h(1,2) - y0)* sin( alpha_c13 )	* sin( beta_c13 )+a_c13 * ( cos( beta_c13 ) )^2	 *(c13h(1,1) - x0)	;
a26_c14		= s_c14 * sin( alpha_c14 ) - m_c14 	*(c14h(1,2) - y0)* sin( alpha_c14 )	* sin( beta_c14 )+a_c14 * ( cos( beta_c14 ) )^2	 *(c14h(1,1) - x0)	;
a26_c15		= s_c15 * sin( alpha_c15 ) - m_c15 	*(c15h(1,2) - y0)* sin( alpha_c15 )	* sin( beta_c15 )+a_c15 * ( cos( beta_c15 ) )^2	 *(c15h(1,1) - x0)	;
a26_c16		= s_c16 * sin( alpha_c16 ) - m_c16 	*(c16h(1,2) - y0)* sin( alpha_c16 )	* sin( beta_c16 )+a_c16 * ( cos( beta_c16 ) )^2	 *(c16h(1,1) - x0)	;
a26_c17		= s_c17 * sin( alpha_c17 ) - m_c17 	*(c17h(1,2) - y0)* sin( alpha_c17 )	* sin( beta_c17 )+a_c17 * ( cos( beta_c17 ) )^2	 *(c17h(1,1) - x0)	;
a26_c18		= s_c18 * sin( alpha_c18 ) - m_c18 	*(c18h(1,2) - y0)* sin( alpha_c18 )	* sin( beta_c18 )+a_c18 * ( cos( beta_c18 ) )^2	 *(c18h(1,1) - x0)	;
a26_c19		= s_c19 * sin( alpha_c19 ) - m_c19 	*(c19h(1,2) - y0)* sin( alpha_c19 )	* sin( beta_c19 )+a_c19 * ( cos( beta_c19 ) )^2	 *(c19h(1,1) - x0)	;
a26_c20		= s_c20 * sin( alpha_c20 ) - m_c20 	*(c20h(1,2) - y0)* sin( alpha_c20 )	* sin( beta_c20 )+a_c20 * ( cos( beta_c20 ) )^2	 *(c20h(1,1) - x0)	;
a26_p1 		= s_p1	* sin( alpha_p1  ) - m_p1	*(p1a(1,2)  - y0)* sin( alpha_p1  )	* sin( beta_p1  )+a_p1  * ( cos( beta_p1  ) )^2	 *(p1a(1,1)  - x0) 	;
a26_p2 		= s_p2	* sin( alpha_p2  ) - m_p2	*(p2a(1,2)  - y0)* sin( alpha_p2  )	* sin( beta_p2  )+a_p2  * ( cos( beta_p2  ) )^2	 *(p2a(1,1)  - x0) 	;
a26_p3 		= s_p3	* sin( alpha_p3  ) - m_p3	*(p3a(1,2)  - y0)* sin( alpha_p3  )	* sin( beta_p3  )+a_p3  * ( cos( beta_p3  ) )^2	 *(p3a(1,1)  - x0) 	;
a26_p4 		= s_p4	* sin( alpha_p4  ) - m_p4	*(p4a(1,2)  - y0)* sin( alpha_p4  )	* sin( beta_p4  )+a_p4  * ( cos( beta_p4  ) )^2	 *(p4a(1,1)  - x0) 	;
a26_b1 		= s_b1	* sin( alpha_b1  ) - m_b1	*(b1s(1,2)  - y0)* sin( alpha_b1  )	* sin( beta_b1  )+a_b1  * ( cos( beta_b1  ) )^2	 *(b1s(1,1)  - x0) 	;
a26_b2 		= s_b2	* sin( alpha_b2  ) - m_b2	*(b2s(1,2)  - y0)* sin( alpha_b2  )	* sin( beta_b2  )+a_b2  * ( cos( beta_b2  ) )^2	 *(b2s(1,1)  - x0) 	;
a26_b3 		= s_b3	* sin( alpha_b3  ) - m_b3	*(b3s(1,2)  - y0)* sin( alpha_b3  )	* sin( beta_b3  )+a_b3  * ( cos( beta_b3  ) )^2	 *(b3s(1,1)  - x0) 	;
a26_b4 		= s_b4	* sin( alpha_b4  ) - m_b4	*(b4s(1,2)  - y0)* sin( alpha_b4  )	* sin( beta_b4  )+a_b4  * ( cos( beta_b4  ) )^2	 *(b4s(1,1)  - x0) 	;
  
a26(1,1) = a26_c1 		;
a26(1,2) = a26_c2 		;
a26(1,3) = a26_c3 		;
a26(1,4) = a26_c4 		;
a26(1,5) = a26_c5 		;
a26(1,6) = a26_c6 		;
a26(1,7) = a26_c7 		;
a26(1,8) = a26_c8 		;
a26(1,9) = a26_c9 		;
a26(1,10)= a26_c10		;
a26(1,11)= a26_c11		;
a26(1,12)= a26_c12		;
a26(1,13)= a26_c13		;
a26(1,14)= a26_c14		;
a26(1,15)= a26_c15		;
a26(1,16)= a26_c16		;
a26(1,17)= a26_c17		;
a26(1,18)= a26_c18		;
a26(1,19)= a26_c19		;
a26(1,20)= a26_c20		;
a26(1,21)= a26_p1 		;
a26(1,22)= a26_p2 		;
a26(1,23)= a26_p3 		;
a26(1,24)= a26_p4 		;
a26(1,25)= a26_b1 		;
a26(1,26)= a26_b2 		;
a26(1,27)= a26_b3 		;
a26(1,28)= a26_b4 		;

a62_c1 		= a26_c1 		;
a62_c2 		= a26_c2 		;
a62_c3 		= a26_c3 		;
a62_c4 		= a26_c4 		;
a62_c5 		= a26_c5 		;
a62_c6 		= a26_c6 		;
a62_c7 		= a26_c7 		;
a62_c8 		= a26_c8 		;
a62_c9 		= a26_c9 		;
a62_c10		= a26_c10		;
a62_c11		= a26_c11		;
a62_c12		= a26_c12		;
a62_c13		= a26_c13		;
a62_c14		= a26_c14		;
a62_c15		= a26_c15		;
a62_c16		= a26_c16		;
a62_c17		= a26_c17		;
a62_c18		= a26_c18		;
a62_c19		= a26_c19		;
a62_c20		= a26_c20		;
a62_p1 		= a26_p1 		;
a62_p2 		= a26_p2 		;
a62_p3 		= a26_p3 		;
a62_p4 		= a26_p4 		;
a62_b1 		= a26_b1 		;
a62_b2 		= a26_b2 		;
a62_b3 		= a26_b3 		;
a62_b4 		= a26_b4 		;
 
a62(1,1) = a62_c1 		;
a62(1,2) = a62_c2 		;
a62(1,3) = a62_c3 		;
a62(1,4) = a62_c4 		;
a62(1,5) = a62_c5 		;
a62(1,6) = a62_c6 		;
a62(1,7) = a62_c7 		;
a62(1,8) = a62_c8 		;
a62(1,9) = a62_c9 		;
a62(1,10)= a62_c10		;
a62(1,11)= a62_c11		;
a62(1,12)= a62_c12		;
a62(1,13)= a62_c13		;
a62(1,14)= a62_c14		;
a62(1,15)= a62_c15		;
a62(1,16)= a62_c16		;
a62(1,17)= a62_c17		;
a62(1,18)= a62_c18		;
a62(1,19)= a62_c19		;
a62(1,20)= a62_c20		;
a62(1,21)= a62_p1 		;
a62(1,22)= a62_p2 		;
a62(1,23)= a62_p3 		;
a62(1,24)= a62_p4 		;
a62(1,25)= a62_b1 		;
a62(1,26)= a62_b2 		;
a62(1,27)= a62_b3 		;
a62(1,28)= a62_b4 		;

a34_c1 		= s_c1  * sin( beta_c1  ) - m_c1  	*(c1h(1,3)  - z0)* sin( gamma_c1  )	* sin( beta_c1  )	+a_c1  * ( cos( gamma_c1  ) )^2	 *(c1h(1,2)  - y0)	;
a34_c2 		= s_c2  * sin( beta_c2  ) - m_c2  	*(c2h(1,3)  - z0)* sin( gamma_c2  )	* sin( beta_c2  )	+a_c2  * ( cos( gamma_c2  ) )^2	 *(c2h(1,2)  - y0)	;
a34_c3 		= s_c3  * sin( beta_c3  ) - m_c3  	*(c3h(1,3)  - z0)* sin( gamma_c3  )	* sin( beta_c3  )	+a_c3  * ( cos( gamma_c3  ) )^2	 *(c3h(1,2)  - y0)	;
a34_c4 		= s_c4  * sin( beta_c4  ) - m_c4  	*(c4h(1,3)  - z0)* sin( gamma_c4  )	* sin( beta_c4  )	+a_c4  * ( cos( gamma_c4  ) )^2	 *(c4h(1,2)  - y0)	;
a34_c5 		= s_c5  * sin( beta_c5  ) - m_c5  	*(c5h(1,3)  - z0)* sin( gamma_c5  )	* sin( beta_c5  )	+a_c5  * ( cos( gamma_c5  ) )^2	 *(c5h(1,2)  - y0)	;
a34_c6 		= s_c6  * sin( beta_c6  ) - m_c6  	*(c6h(1,3)  - z0)* sin( gamma_c6  )	* sin( beta_c6  )	+a_c6  * ( cos( gamma_c6  ) )^2	 *(c6h(1,2)  - y0)	;
a34_c7 		= s_c7  * sin( beta_c7  ) - m_c7  	*(c7h(1,3)  - z0)* sin( gamma_c7  )	* sin( beta_c7  )	+a_c7  * ( cos( gamma_c7  ) )^2	 *(c7h(1,2)  - y0)	;
a34_c8 		= s_c8  * sin( beta_c8  ) - m_c8  	*(c8h(1,3)  - z0)* sin( gamma_c8  )	* sin( beta_c8  )	+a_c8  * ( cos( gamma_c8  ) )^2	 *(c8h(1,2)  - y0)	;
a34_c9 		= s_c9  * sin( beta_c9  ) - m_c9  	*(c9h(1,3)  - z0)* sin( gamma_c9  )	* sin( beta_c9  )	+a_c9  * ( cos( gamma_c9  ) )^2	 *(c9h(1,2)  - y0)	;
a34_c10		= s_c10 * sin( beta_c10 ) - m_c10 	*(c10h(1,3) - z0)* sin( gamma_c10 )	* sin( beta_c10 )	+a_c10 * ( cos( gamma_c10 ) )^2	 *(c10h(1,2) - y0)	;
a34_c11		= s_c11 * sin( beta_c11 ) - m_c11 	*(c11h(1,3) - z0)* sin( gamma_c11 )	* sin( beta_c11 )	+a_c11 * ( cos( gamma_c11 ) )^2	 *(c11h(1,2) - y0)	;
a34_c12		= s_c12 * sin( beta_c12 ) - m_c12 	*(c12h(1,3) - z0)* sin( gamma_c12 )	* sin( beta_c12 )	+a_c12 * ( cos( gamma_c12 ) )^2	 *(c12h(1,2) - y0)	;
a34_c13		= s_c13 * sin( beta_c13 ) - m_c13 	*(c13h(1,3) - z0)* sin( gamma_c13 )	* sin( beta_c13 )	+a_c13 * ( cos( gamma_c13 ) )^2	 *(c13h(1,2) - y0)	;
a34_c14		= s_c14 * sin( beta_c14 ) - m_c14 	*(c14h(1,3) - z0)* sin( gamma_c14 )	* sin( beta_c14 )	+a_c14 * ( cos( gamma_c14 ) )^2	 *(c14h(1,2) - y0)	;
a34_c15		= s_c15 * sin( beta_c15 ) - m_c15 	*(c15h(1,3) - z0)* sin( gamma_c15 )	* sin( beta_c15 )	+a_c15 * ( cos( gamma_c15 ) )^2	 *(c15h(1,2) - y0)	;
a34_c16		= s_c16 * sin( beta_c16 ) - m_c16 	*(c16h(1,3) - z0)* sin( gamma_c16 )	* sin( beta_c16 )	+a_c16 * ( cos( gamma_c16 ) )^2	 *(c16h(1,2) - y0)	;
a34_c17		= s_c17 * sin( beta_c17 ) - m_c17 	*(c17h(1,3) - z0)* sin( gamma_c17 )	* sin( beta_c17 )	+a_c17 * ( cos( gamma_c17 ) )^2	 *(c17h(1,2) - y0)	;
a34_c18		= s_c18 * sin( beta_c18 ) - m_c18 	*(c18h(1,3) - z0)* sin( gamma_c18 )	* sin( beta_c18 )	+a_c18 * ( cos( gamma_c18 ) )^2	 *(c18h(1,2) - y0)	;
a34_c19		= s_c19 * sin( beta_c19 ) - m_c19 	*(c19h(1,3) - z0)* sin( gamma_c19 )	* sin( beta_c19 )	+a_c19 * ( cos( gamma_c19 ) )^2	 *(c19h(1,2) - y0)	;
a34_c20		= s_c20 * sin( beta_c20 ) - m_c20 	*(c20h(1,3) - z0)* sin( gamma_c20 )	* sin( beta_c20 )	+a_c20 * ( cos( gamma_c20 ) )^2	 *(c20h(1,2) - y0)	;
a34_p1 		= s_p1	* sin( beta_p1  ) - m_p1	*(p1a(1,3)  - z0)* sin( gamma_p1  )	* sin( beta_p1  )	+a_p1  * ( cos( gamma_p1  ) )^2	 *(p1a(1,2)  - y0) 	;
a34_p2 		= s_p2	* sin( beta_p2  ) - m_p2	*(p2a(1,3)  - z0)* sin( gamma_p2  )	* sin( beta_p2  )	+a_p2  * ( cos( gamma_p2  ) )^2	 *(p2a(1,2)  - y0) 	;
a34_p3 		= s_p3	* sin( beta_p3  ) - m_p3	*(p3a(1,3)  - z0)* sin( gamma_p3  )	* sin( beta_p3  )	+a_p3  * ( cos( gamma_p3  ) )^2	 *(p3a(1,2)  - y0) 	;
a34_p4 		= s_p4	* sin( beta_p4  ) - m_p4	*(p4a(1,3)  - z0)* sin( gamma_p4  )	* sin( beta_p4  )	+a_p4  * ( cos( gamma_p4  ) )^2	 *(p4a(1,2)  - y0) 	;
a34_b1 		= s_b1	* sin( beta_b1  ) - m_b1	*(b1s(1,3)  - z0)* sin( gamma_b1  )	* sin( beta_b1  )	+a_b1  * ( cos( gamma_b1  ) )^2	 *(b1s(1,2)  - y0) 	;
a34_b2 		= s_b2	* sin( beta_b2  ) - m_b2	*(b2s(1,3)  - z0)* sin( gamma_b2  )	* sin( beta_b2  )	+a_b2  * ( cos( gamma_b2  ) )^2	 *(b2s(1,2)  - y0) 	;
a34_b3 		= s_b3	* sin( beta_b3  ) - m_b3	*(b3s(1,3)  - z0)* sin( gamma_b3  )	* sin( beta_b3  )	+a_b3  * ( cos( gamma_b3  ) )^2	 *(b3s(1,2)  - y0) 	;
a34_b4 		= s_b4	* sin( beta_b4  ) - m_b4	*(b4s(1,3)  - z0)* sin( gamma_b4  )	* sin( beta_b4  )	+a_b4  * ( cos( gamma_b4  ) )^2	 *(b4s(1,2)  - y0) 	;
  
a34(1,1) = a34_c1 		;
a34(1,2) = a34_c2 		;
a34(1,3) = a34_c3 		;
a34(1,4) = a34_c4 		;
a34(1,5) = a34_c5 		;
a34(1,6) = a34_c6 		;
a34(1,7) = a34_c7 		;
a34(1,8) = a34_c8 		;
a34(1,9) = a34_c9 		;
a34(1,10)= a34_c10		;
a34(1,11)= a34_c11		;
a34(1,12)= a34_c12		;
a34(1,13)= a34_c13		;
a34(1,14)= a34_c14		;
a34(1,15)= a34_c15		;
a34(1,16)= a34_c16		;
a34(1,17)= a34_c17		;
a34(1,18)= a34_c18		;
a34(1,19)= a34_c19		;
a34(1,20)= a34_c20		;
a34(1,21)= a34_p1 		;
a34(1,22)= a34_p2 		;
a34(1,23)= a34_p3 		;
a34(1,24)= a34_p4 		;
a34(1,25)= a34_b1 		;
a34(1,26)= a34_b2 		;
a34(1,27)= a34_b3 		;
a34(1,28)= a34_b4 		;

a43_c1 		= a34_c1 		;
a43_c2 		= a34_c2 		;
a43_c3 		= a34_c3 		;
a43_c4 		= a34_c4 		;
a43_c5 		= a34_c5 		;
a43_c6 		= a34_c6 		;
a43_c7 		= a34_c7 		;
a43_c8 		= a34_c8 		;
a43_c9 		= a34_c9 		;
a43_c10		= a34_c10		;
a43_c11		= a34_c11		;
a43_c12		= a34_c12		;
a43_c13		= a34_c13		;
a43_c14		= a34_c14		;
a43_c15		= a34_c15		;
a43_c16		= a34_c16		;
a43_c17		= a34_c17		;
a43_c18		= a34_c18		;
a43_c19		= a34_c19		;
a43_c20		= a34_c20		;
a43_p1 		= a34_p1 		;
a43_p2 		= a34_p2 		;
a43_p3 		= a34_p3 		;
a43_p4 		= a34_p4 		;
a43_b1 		= a34_b1 		;
a43_b2 		= a34_b2 		;
a43_b3 		= a34_b3 		;
a43_b4 		= a34_b4 		;
 
a43(1,1) = a43_c1 		;
a43(1,2) = a43_c2 		;
a43(1,3) = a43_c3 		;
a43(1,4) = a43_c4 		;
a43(1,5) = a43_c5 		;
a43(1,6) = a43_c6 		;
a43(1,7) = a43_c7 		;
a43(1,8) = a43_c8 		;
a43(1,9) = a43_c9 		;
a43(1,10)= a43_c10		;
a43(1,11)= a43_c11		;
a43(1,12)= a43_c12		;
a43(1,13)= a43_c13		;
a43(1,14)= a43_c14		;
a43(1,15)= a43_c15		;
a43(1,16)= a43_c16		;
a43(1,17)= a43_c17		;
a43(1,18)= a43_c18		;
a43(1,19)= a43_c19		;
a43(1,20)= a43_c20		;
a43(1,21)= a43_p1 		;
a43(1,22)= a43_p2 		;
a43(1,23)= a43_p3 		;
a43(1,24)= a43_p4 		;
a43(1,25)= a43_b1 		;
a43(1,26)= a43_b2 		;
a43(1,27)= a43_b3 		;
a43(1,28)= a43_b4 		;

a35_c1 		= -s_c1  * sin( alpha_c1  ) - m_c1  *(c1h(1,3)  - z0)* sin( gamma_c1  )	* sin( alpha_c1  )-a_c1  * ( cos( gamma_c1  ) )^2	 *(c1h(1,1)  - x0)	;
a35_c2 		= -s_c2  * sin( alpha_c2  ) - m_c2  *(c2h(1,3)  - z0)* sin( gamma_c2  )	* sin( alpha_c2  )-a_c2  * ( cos( gamma_c2  ) )^2	 *(c2h(1,1)  - x0)	;
a35_c3 		= -s_c3  * sin( alpha_c3  ) - m_c3  *(c3h(1,3)  - z0)* sin( gamma_c3  )	* sin( alpha_c3  )-a_c3  * ( cos( gamma_c3  ) )^2	 *(c3h(1,1)  - x0)	;
a35_c4 		= -s_c4  * sin( alpha_c4  ) - m_c4  *(c4h(1,3)  - z0)* sin( gamma_c4  )	* sin( alpha_c4  )-a_c4  * ( cos( gamma_c4  ) )^2	 *(c4h(1,1)  - x0)	;
a35_c5 		= -s_c5  * sin( alpha_c5  ) - m_c5  *(c5h(1,3)  - z0)* sin( gamma_c5  )	* sin( alpha_c5  )-a_c5  * ( cos( gamma_c5  ) )^2	 *(c5h(1,1)  - x0)	;
a35_c6 		= -s_c6  * sin( alpha_c6  ) - m_c6  *(c6h(1,3)  - z0)* sin( gamma_c6  )	* sin( alpha_c6  )-a_c6  * ( cos( gamma_c6  ) )^2	 *(c6h(1,1)  - x0)	;
a35_c7 		= -s_c7  * sin( alpha_c7  ) - m_c7  *(c7h(1,3)  - z0)* sin( gamma_c7  )	* sin( alpha_c7  )-a_c7  * ( cos( gamma_c7  ) )^2	 *(c7h(1,1)  - x0)	;
a35_c8 		= -s_c8  * sin( alpha_c8  ) - m_c8  *(c8h(1,3)  - z0)* sin( gamma_c8  )	* sin( alpha_c8  )-a_c8  * ( cos( gamma_c8  ) )^2	 *(c8h(1,1)  - x0)	;
a35_c9 		= -s_c9  * sin( alpha_c9  ) - m_c9  *(c9h(1,3)  - z0)* sin( gamma_c9  )	* sin( alpha_c9  )-a_c9  * ( cos( gamma_c9  ) )^2	 *(c9h(1,1)  - x0)	;
a35_c10		= -s_c10 * sin( alpha_c10 ) - m_c10 *(c10h(1,3) - z0)* sin( gamma_c10 )	* sin( alpha_c10 )-a_c10 * ( cos( gamma_c10 ) )^2	 *(c10h(1,1) - x0)	;
a35_c11		= -s_c11 * sin( alpha_c11 ) - m_c11 *(c11h(1,3) - z0)* sin( gamma_c11 )	* sin( alpha_c11 )-a_c11 * ( cos( gamma_c11 ) )^2	 *(c11h(1,1) - x0)	;
a35_c12		= -s_c12 * sin( alpha_c12 ) - m_c12 *(c12h(1,3) - z0)* sin( gamma_c12 )	* sin( alpha_c12 )-a_c12 * ( cos( gamma_c12 ) )^2	 *(c12h(1,1) - x0)	;
a35_c13		= -s_c13 * sin( alpha_c13 ) - m_c13 *(c13h(1,3) - z0)* sin( gamma_c13 )	* sin( alpha_c13 )-a_c13 * ( cos( gamma_c13 ) )^2	 *(c13h(1,1) - x0)	;
a35_c14		= -s_c14 * sin( alpha_c14 ) - m_c14 *(c14h(1,3) - z0)* sin( gamma_c14 )	* sin( alpha_c14 )-a_c14 * ( cos( gamma_c14 ) )^2	 *(c14h(1,1) - x0)	;
a35_c15		= -s_c15 * sin( alpha_c15 ) - m_c15 *(c15h(1,3) - z0)* sin( gamma_c15 )	* sin( alpha_c15 )-a_c15 * ( cos( gamma_c15 ) )^2	 *(c15h(1,1) - x0)	;
a35_c16		= -s_c16 * sin( alpha_c16 ) - m_c16 *(c16h(1,3) - z0)* sin( gamma_c16 )	* sin( alpha_c16 )-a_c16 * ( cos( gamma_c16 ) )^2	 *(c16h(1,1) - x0)	;
a35_c17		= -s_c17 * sin( alpha_c17 ) - m_c17 *(c17h(1,3) - z0)* sin( gamma_c17 )	* sin( alpha_c17 )-a_c17 * ( cos( gamma_c17 ) )^2	 *(c17h(1,1) - x0)	;
a35_c18		= -s_c18 * sin( alpha_c18 ) - m_c18 *(c18h(1,3) - z0)* sin( gamma_c18 )	* sin( alpha_c18 )-a_c18 * ( cos( gamma_c18 ) )^2	 *(c18h(1,1) - x0)	;
a35_c19		= -s_c19 * sin( alpha_c19 ) - m_c19 *(c19h(1,3) - z0)* sin( gamma_c19 )	* sin( alpha_c19 )-a_c19 * ( cos( gamma_c19 ) )^2	 *(c19h(1,1) - x0)	;
a35_c20		= -s_c20 * sin( alpha_c20 ) - m_c20 *(c20h(1,3) - z0)* sin( gamma_c20 )	* sin( alpha_c20 )-a_c20 * ( cos( gamma_c20 ) )^2	 *(c20h(1,1) - x0)	;
a35_p1 		= -s_p1	* sin( alpha_p1  )  - m_p1	*(p1a(1,3)  - z0)* sin( gamma_p1  )	* sin( alpha_p1  )-a_p1  * ( cos( gamma_p1  ) )^2	 *(p1a(1,1)  - x0) 	;
a35_p2 		= -s_p2	* sin( alpha_p2  )  - m_p2	*(p2a(1,3)  - z0)* sin( gamma_p2  )	* sin( alpha_p2  )-a_p2  * ( cos( gamma_p2  ) )^2	 *(p2a(1,1)  - x0) 	;
a35_p3 		= -s_p3	* sin( alpha_p3  )  - m_p3	*(p3a(1,3)  - z0)* sin( gamma_p3  )	* sin( alpha_p3  )-a_p3  * ( cos( gamma_p3  ) )^2	 *(p3a(1,1)  - x0) 	;
a35_p4 		= -s_p4	* sin( alpha_p4  )  - m_p4	*(p4a(1,3)  - z0)* sin( gamma_p4  )	* sin( alpha_p4  )-a_p4  * ( cos( gamma_p4  ) )^2	 *(p4a(1,1)  - x0) 	;
a35_b1 		= -s_b1	* sin( alpha_b1  )  - m_b1	*(b1s(1,3)  - z0)* sin( gamma_b1  )	* sin( alpha_b1  )-a_b1  * ( cos( gamma_b1  ) )^2	 *(b1s(1,1)  - x0) 	;
a35_b2 		= -s_b2	* sin( alpha_b2  )  - m_b2	*(b2s(1,3)  - z0)* sin( gamma_b2  )	* sin( alpha_b2  )-a_b2  * ( cos( gamma_b2  ) )^2	 *(b2s(1,1)  - x0) 	;
a35_b3 		= -s_b3	* sin( alpha_b3  )  - m_b3	*(b3s(1,3)  - z0)* sin( gamma_b3  )	* sin( alpha_b3  )-a_b3  * ( cos( gamma_b3  ) )^2	 *(b3s(1,1)  - x0) 	;
a35_b4 		= -s_b4	* sin( alpha_b4  )  - m_b4	*(b4s(1,3)  - z0)* sin( gamma_b4  )	* sin( alpha_b4  )-a_b4  * ( cos( gamma_b4  ) )^2	 *(b4s(1,1)  - x0) 	;
  
a35(1,1) = a35_c1 		;
a35(1,2) = a35_c2 		;
a35(1,3) = a35_c3 		;
a35(1,4) = a35_c4 		;
a35(1,5) = a35_c5 		;
a35(1,6) = a35_c6 		;
a35(1,7) = a35_c7 		;
a35(1,8) = a35_c8 		;
a35(1,9) = a35_c9 		;
a35(1,10)= a35_c10		;
a35(1,11)= a35_c11		;
a35(1,12)= a35_c12		;
a35(1,13)= a35_c13		;
a35(1,14)= a35_c14		;
a35(1,15)= a35_c15		;
a35(1,16)= a35_c16		;
a35(1,17)= a35_c17		;
a35(1,18)= a35_c18		;
a35(1,19)= a35_c19		;
a35(1,20)= a35_c20		;
a35(1,21)= a35_p1 		;
a35(1,22)= a35_p2 		;
a35(1,23)= a35_p3 		;
a35(1,24)= a35_p4 		;
a35(1,25)= a35_b1 		;
a35(1,26)= a35_b2 		;
a35(1,27)= a35_b3 		;
a35(1,28)= a35_b4 		;

a53_c1 		= a35_c1 		;
a53_c2 		= a35_c2 		;
a53_c3 		= a35_c3 		;
a53_c4 		= a35_c4 		;
a53_c5 		= a35_c5 		;
a53_c6 		= a35_c6 		;
a53_c7 		= a35_c7 		;
a53_c8 		= a35_c8 		;
a53_c9 		= a35_c9 		;
a53_c10		= a35_c10		;
a53_c11		= a35_c11		;
a53_c12		= a35_c12		;
a53_c13		= a35_c13		;
a53_c14		= a35_c14		;
a53_c15		= a35_c15		;
a53_c16		= a35_c16		;
a53_c17		= a35_c17		;
a53_c18		= a35_c18		;
a53_c19		= a35_c19		;
a53_c20		= a35_c20		;
a53_p1 		= a35_p1 		;
a53_p2 		= a35_p2 		;
a53_p3 		= a35_p3 		;
a53_p4 		= a35_p4 		;
a53_b1 		= a35_b1 		;
a53_b2 		= a35_b2 		;
a53_b3 		= a35_b3 		;
a53_b4 		= a35_b4 		;
 
a53(1,1) = a53_c1 		;
a53(1,2) = a53_c2 		;
a53(1,3) = a53_c3 		;
a53(1,4) = a53_c4 		;
a53(1,5) = a53_c5 		;
a53(1,6) = a53_c6 		;
a53(1,7) = a53_c7 		;
a53(1,8) = a53_c8 		;
a53(1,9) = a53_c9 		;
a53(1,10)= a53_c10		;
a53(1,11)= a53_c11		;
a53(1,12)= a53_c12		;
a53(1,13)= a53_c13		;
a53(1,14)= a53_c14		;
a53(1,15)= a53_c15		;
a53(1,16)= a53_c16		;
a53(1,17)= a53_c17		;
a53(1,18)= a53_c18		;
a53(1,19)= a53_c19		;
a53(1,20)= a53_c20		;
a53(1,21)= a53_p1 		;
a53(1,22)= a53_p2 		;
a53(1,23)= a53_p3 		;
a53(1,24)= a53_p4 		;
a53(1,25)= a53_b1 		;
a53(1,26)= a53_b2 		;
a53(1,27)= a53_b3 		;
a53(1,28)= a53_b4 		;

a36_c1 		= -a_c1  * sin( beta_c1  )	* sin( gamma_c1  )* (c1h(1,1)  - x0) + a_c1  * sin( gamma_c1  )	* sin( alpha_c1  ) * (c1h(1,2)  - y0);
a36_c2 		= -a_c2  * sin( beta_c2  )	* sin( gamma_c2  )* (c2h(1,1)  - x0) + a_c2  * sin( gamma_c2  )	* sin( alpha_c2  ) * (c2h(1,2)  - y0);
a36_c3 		= -a_c3  * sin( beta_c3  )	* sin( gamma_c3  )* (c3h(1,1)  - x0) + a_c3  * sin( gamma_c3  )	* sin( alpha_c3  ) * (c3h(1,2)  - y0);
a36_c4 		= -a_c4  * sin( beta_c4  )	* sin( gamma_c4  )* (c4h(1,1)  - x0) + a_c4  * sin( gamma_c4  )	* sin( alpha_c4  ) * (c4h(1,2)  - y0);
a36_c5 		= -a_c5  * sin( beta_c5  )	* sin( gamma_c5  )* (c5h(1,1)  - x0) + a_c5  * sin( gamma_c5  )	* sin( alpha_c5  ) * (c5h(1,2)  - y0);
a36_c6 		= -a_c6  * sin( beta_c6  )	* sin( gamma_c6  )* (c6h(1,1)  - x0) + a_c6  * sin( gamma_c6  )	* sin( alpha_c6  ) * (c6h(1,2)  - y0);
a36_c7 		= -a_c7  * sin( beta_c7  )	* sin( gamma_c7  )* (c7h(1,1)  - x0) + a_c7  * sin( gamma_c7  )	* sin( alpha_c7  ) * (c7h(1,2)  - y0);
a36_c8 		= -a_c8  * sin( beta_c8  )	* sin( gamma_c8  )* (c8h(1,1)  - x0) + a_c8  * sin( gamma_c8  )	* sin( alpha_c8  ) * (c8h(1,2)  - y0);
a36_c9 		= -a_c9  * sin( beta_c9  )	* sin( gamma_c9  )* (c9h(1,1)  - x0) + a_c9  * sin( gamma_c9  )	* sin( alpha_c9  ) * (c9h(1,2)  - y0);
a36_c10		= -a_c10 * sin( beta_c10 )	* sin( gamma_c10 )* (c10h(1,1) - x0) + a_c10 * sin( gamma_c10 )	* sin( alpha_c10 ) * (c10h(1,2) - y0);
a36_c11		= -a_c11 * sin( beta_c11 )	* sin( gamma_c11 )* (c11h(1,1) - x0) + a_c11 * sin( gamma_c11 )	* sin( alpha_c11 ) * (c11h(1,2) - y0);
a36_c12		= -a_c12 * sin( beta_c12 )	* sin( gamma_c12 )* (c12h(1,1) - x0) + a_c12 * sin( gamma_c12 )	* sin( alpha_c12 ) * (c12h(1,2) - y0);
a36_c13		= -a_c13 * sin( beta_c13 )	* sin( gamma_c13 )* (c13h(1,1) - x0) + a_c13 * sin( gamma_c13 )	* sin( alpha_c13 ) * (c13h(1,2) - y0);
a36_c14		= -a_c14 * sin( beta_c14 )	* sin( gamma_c14 )* (c14h(1,1) - x0) + a_c14 * sin( gamma_c14 )	* sin( alpha_c14 ) * (c14h(1,2) - y0);
a36_c15		= -a_c15 * sin( beta_c15 )	* sin( gamma_c15 )* (c15h(1,1) - x0) + a_c15 * sin( gamma_c15 )	* sin( alpha_c15 ) * (c15h(1,2) - y0);
a36_c16		= -a_c16 * sin( beta_c16 )	* sin( gamma_c16 )* (c16h(1,1) - x0) + a_c16 * sin( gamma_c16 )	* sin( alpha_c16 ) * (c16h(1,2) - y0);
a36_c17		= -a_c17 * sin( beta_c17 )	* sin( gamma_c17 )* (c17h(1,1) - x0) + a_c17 * sin( gamma_c17 )	* sin( alpha_c17 ) * (c17h(1,2) - y0);
a36_c18		= -a_c18 * sin( beta_c18 )	* sin( gamma_c18 )* (c18h(1,1) - x0) + a_c18 * sin( gamma_c18 )	* sin( alpha_c18 ) * (c18h(1,2) - y0);
a36_c19		= -a_c19 * sin( beta_c19 )	* sin( gamma_c19 )* (c19h(1,1) - x0) + a_c19 * sin( gamma_c19 )	* sin( alpha_c19 ) * (c19h(1,2) - y0);
a36_c20		= -a_c20 * sin( beta_c20 )	* sin( gamma_c20 )* (c20h(1,1) - x0) + a_c20 * sin( gamma_c20 )	* sin( alpha_c20 ) * (c20h(1,2) - y0);
a36_p1 		= -a_p1  * sin( beta_p1  )	* sin( gamma_p1  )* (p1a(1,1)  - x0) + a_p1	 * sin( gamma_p1  )	* sin( alpha_p1  ) * (p1a(1,2)  - y0);
a36_p2 		= -a_p2  * sin( beta_p2  )	* sin( gamma_p2  )* (p2a(1,1)  - x0) + a_p2	 * sin( gamma_p2  )	* sin( alpha_p2  ) * (p2a(1,2)  - y0);
a36_p3 		= -a_p3  * sin( beta_p3  )	* sin( gamma_p3  )* (p3a(1,1)  - x0) + a_p3	 * sin( gamma_p3  )	* sin( alpha_p3  ) * (p3a(1,2)  - y0);
a36_p4 		= -a_p4  * sin( beta_p4  )	* sin( gamma_p4  )* (p4a(1,1)  - x0) + a_p4	 * sin( gamma_p4  )	* sin( alpha_p4  ) * (p4a(1,2)  - y0);
a36_b1 		= -a_b1  * sin( beta_b1  )	* sin( gamma_b1  )* (b1s(1,1)  - x0) + a_b1	 * sin( gamma_b1  )	* sin( alpha_b1  ) * (b1s(1,2)  - y0);
a36_b2 		= -a_b2  * sin( beta_b2  )	* sin( gamma_b2  )* (b2s(1,1)  - x0) + a_b2	 * sin( gamma_b2  )	* sin( alpha_b2  ) * (b2s(1,2)  - y0);
a36_b3 		= -a_b3  * sin( beta_b3  )	* sin( gamma_b3  )* (b3s(1,1)  - x0) + a_b3	 * sin( gamma_b3  )	* sin( alpha_b3  ) * (b3s(1,2)  - y0);
a36_b4 		= -a_b4  * sin( beta_b4  )	* sin( gamma_b4  )* (b4s(1,1)  - x0) + a_b4	 * sin( gamma_b4  )	* sin( alpha_b4  ) * (b4s(1,2)  - y0);
  
a36(1,1) = a36_c1 		;
a36(1,2) = a36_c2 		;
a36(1,3) = a36_c3 		;
a36(1,4) = a36_c4 		;
a36(1,5) = a36_c5 		;
a36(1,6) = a36_c6 		;
a36(1,7) = a36_c7 		;
a36(1,8) = a36_c8 		;
a36(1,9) = a36_c9 		;
a36(1,10)= a36_c10		;
a36(1,11)= a36_c11		;
a36(1,12)= a36_c12		;
a36(1,13)= a36_c13		;
a36(1,14)= a36_c14		;
a36(1,15)= a36_c15		;
a36(1,16)= a36_c16		;
a36(1,17)= a36_c17		;
a36(1,18)= a36_c18		;
a36(1,19)= a36_c19		;
a36(1,20)= a36_c20		;
a36(1,21)= a36_p1 		;
a36(1,22)= a36_p2 		;
a36(1,23)= a36_p3 		;
a36(1,24)= a36_p4 		;
a36(1,25)= a36_b1 		;
a36(1,26)= a36_b2 		;
a36(1,27)= a36_b3 		;
a36(1,28)= a36_b4 		;

a63_c1 		= a36_c1 		;
a63_c2 		= a36_c2 		;
a63_c3 		= a36_c3 		;
a63_c4 		= a36_c4 		;
a63_c5 		= a36_c5 		;
a63_c6 		= a36_c6 		;
a63_c7 		= a36_c7 		;
a63_c8 		= a36_c8 		;
a63_c9 		= a36_c9 		;
a63_c10		= a36_c10		;
a63_c11		= a36_c11		;
a63_c12		= a36_c12		;
a63_c13		= a36_c13		;
a63_c14		= a36_c14		;
a63_c15		= a36_c15		;
a63_c16		= a36_c16		;
a63_c17		= a36_c17		;
a63_c18		= a36_c18		;
a63_c19		= a36_c19		;
a63_c20		= a36_c20		;
a63_p1 		= a36_p1 		;
a63_p2 		= a36_p2 		;
a63_p3 		= a36_p3 		;
a63_p4 		= a36_p4 		;
a63_b1 		= a36_b1 		;
a63_b2 		= a36_b2 		;
a63_b3 		= a36_b3 		;
a63_b4 		= a36_b4 		;
 
a63(1,1) = a63_c1 		;
a63(1,2) = a63_c2 		;
a63(1,3) = a63_c3 		;
a63(1,4) = a63_c4 		;
a63(1,5) = a63_c5 		;
a63(1,6) = a63_c6 		;
a63(1,7) = a63_c7 		;
a63(1,8) = a63_c8 		;
a63(1,9) = a63_c9 		;
a63(1,10)= a63_c10		;
a63(1,11)= a63_c11		;
a63(1,12)= a63_c12		;
a63(1,13)= a63_c13		;
a63(1,14)= a63_c14		;
a63(1,15)= a63_c15		;
a63(1,16)= a63_c16		;
a63(1,17)= a63_c17		;
a63(1,18)= a63_c18		;
a63(1,19)= a63_c19		;
a63(1,20)= a63_c20		;
a63(1,21)= a63_p1 		;
a63(1,22)= a63_p2 		;
a63(1,23)= a63_p3 		;
a63(1,24)= a63_p4 		;
a63(1,25)= a63_b1 		;
a63(1,26)= a63_b2 		;
a63(1,27)= a63_b3 		;
a63(1,28)= a63_b4 		;

a44_c1 		= I_c1  *( ( sin( beta_c1  ) )^2 + ( sin( gamma_c1  ) )^2) +2*s_c1  * ( sin( beta_c1  )*(c1h(1,2)  - y0)+(c1h(1,3)  - z0)* sin( gamma_c1  ) ) + a_c1  * ( ((c1h(1,2)  - y0))^2 * ( cos( gamma_c1  ) )^2+ ((c1h(1,3)  - z0))^2 * ( cos( beta_c1  ) )^2+ 2*(c1h(1,2)  - y0)*(c1h(1,3)  - z0)* sin( beta_c1  )	* sin( gamma_c1  )  )  ;
a44_c2 		= I_c2  *( ( sin( beta_c2  ) )^2 + ( sin( gamma_c2  ) )^2) +2*s_c2  * ( sin( beta_c2  )*(c2h(1,2)  - y0)+(c2h(1,3)  - z0)* sin( gamma_c2  ) ) + a_c2  * ( ((c2h(1,2)  - y0))^2 * ( cos( gamma_c2  ) )^2+ ((c2h(1,3)  - z0))^2 * ( cos( beta_c2  ) )^2+ 2*(c2h(1,2)  - y0)*(c2h(1,3)  - z0)* sin( beta_c2  )	* sin( gamma_c2  )  )  ;
a44_c3 		= I_c3  *( ( sin( beta_c3  ) )^2 + ( sin( gamma_c3  ) )^2) +2*s_c3  * ( sin( beta_c3  )*(c3h(1,2)  - y0)+(c3h(1,3)  - z0)* sin( gamma_c3  ) ) + a_c3  * ( ((c3h(1,2)  - y0))^2 * ( cos( gamma_c3  ) )^2+ ((c3h(1,3)  - z0))^2 * ( cos( beta_c3  ) )^2+ 2*(c3h(1,2)  - y0)*(c3h(1,3)  - z0)* sin( beta_c3  )	* sin( gamma_c3  )  )  ;
a44_c4 		= I_c4  *( ( sin( beta_c4  ) )^2 + ( sin( gamma_c4  ) )^2) +2*s_c4  * ( sin( beta_c4  )*(c4h(1,2)  - y0)+(c4h(1,3)  - z0)* sin( gamma_c4  ) ) + a_c4  * ( ((c4h(1,2)  - y0))^2 * ( cos( gamma_c4  ) )^2+ ((c4h(1,3)  - z0))^2 * ( cos( beta_c4  ) )^2+ 2*(c4h(1,2)  - y0)*(c4h(1,3)  - z0)* sin( beta_c4  )	* sin( gamma_c4  )  )  ;
a44_c5 		= I_c5  *( ( sin( beta_c5  ) )^2 + ( sin( gamma_c5  ) )^2) +2*s_c5  * ( sin( beta_c5  )*(c5h(1,2)  - y0)+(c5h(1,3)  - z0)* sin( gamma_c5  ) ) + a_c5  * ( ((c5h(1,2)  - y0))^2 * ( cos( gamma_c5  ) )^2+ ((c5h(1,3)  - z0))^2 * ( cos( beta_c5  ) )^2+ 2*(c5h(1,2)  - y0)*(c5h(1,3)  - z0)* sin( beta_c5  )	* sin( gamma_c5  )  )  ;
a44_c6 		= I_c6  *( ( sin( beta_c6  ) )^2 + ( sin( gamma_c6  ) )^2) +2*s_c6  * ( sin( beta_c6  )*(c6h(1,2)  - y0)+(c6h(1,3)  - z0)* sin( gamma_c6  ) ) + a_c6  * ( ((c6h(1,2)  - y0))^2 * ( cos( gamma_c6  ) )^2+ ((c6h(1,3)  - z0))^2 * ( cos( beta_c6  ) )^2+ 2*(c6h(1,2)  - y0)*(c6h(1,3)  - z0)* sin( beta_c6  )	* sin( gamma_c6  )  )  ;
a44_c7 		= I_c7  *( ( sin( beta_c7  ) )^2 + ( sin( gamma_c7  ) )^2) +2*s_c7  * ( sin( beta_c7  )*(c7h(1,2)  - y0)+(c7h(1,3)  - z0)* sin( gamma_c7  ) ) + a_c7  * ( ((c7h(1,2)  - y0))^2 * ( cos( gamma_c7  ) )^2+ ((c7h(1,3)  - z0))^2 * ( cos( beta_c7  ) )^2+ 2*(c7h(1,2)  - y0)*(c7h(1,3)  - z0)* sin( beta_c7  )	* sin( gamma_c7  )  )  ;
a44_c8 		= I_c8  *( ( sin( beta_c8  ) )^2 + ( sin( gamma_c8  ) )^2) +2*s_c8  * ( sin( beta_c8  )*(c8h(1,2)  - y0)+(c8h(1,3)  - z0)* sin( gamma_c8  ) ) + a_c8  * ( ((c8h(1,2)  - y0))^2 * ( cos( gamma_c8  ) )^2+ ((c8h(1,3)  - z0))^2 * ( cos( beta_c8  ) )^2+ 2*(c8h(1,2)  - y0)*(c8h(1,3)  - z0)* sin( beta_c8  )	* sin( gamma_c8  )  )  ;
a44_c9 		= I_c9  *( ( sin( beta_c9  ) )^2 + ( sin( gamma_c9  ) )^2) +2*s_c9  * ( sin( beta_c9  )*(c9h(1,2)  - y0)+(c9h(1,3)  - z0)* sin( gamma_c9  ) ) + a_c9  * ( ((c9h(1,2)  - y0))^2 * ( cos( gamma_c9  ) )^2+ ((c9h(1,3)  - z0))^2 * ( cos( beta_c9  ) )^2+ 2*(c9h(1,2)  - y0)*(c9h(1,3)  - z0)* sin( beta_c9  )	* sin( gamma_c9  )  )  ;
a44_c10		= I_c10 *( ( sin( beta_c10 ) )^2 + ( sin( gamma_c10 ) )^2) +2*s_c10 * ( sin( beta_c10 )*(c10h(1,2) - y0)+(c10h(1,3) - z0)* sin( gamma_c10 ) ) + a_c10 * ( ((c10h(1,2) - y0))^2 * ( cos( gamma_c10 ) )^2+ ((c10h(1,3) - z0))^2 * ( cos( beta_c10 ) )^2+ 2*(c10h(1,2) - y0)*(c10h(1,3) - z0)* sin( beta_c10 )	* sin( gamma_c10 )  )  ;
a44_c11		= I_c11 *( ( sin( beta_c11 ) )^2 + ( sin( gamma_c11 ) )^2) +2*s_c11 * ( sin( beta_c11 )*(c11h(1,2) - y0)+(c11h(1,3) - z0)* sin( gamma_c11 ) ) + a_c11 * ( ((c11h(1,2) - y0))^2 * ( cos( gamma_c11 ) )^2+ ((c11h(1,3) - z0))^2 * ( cos( beta_c11 ) )^2+ 2*(c11h(1,2) - y0)*(c11h(1,3) - z0)* sin( beta_c11 )	* sin( gamma_c11 )  )  ;
a44_c12		= I_c12 *( ( sin( beta_c12 ) )^2 + ( sin( gamma_c12 ) )^2) +2*s_c12 * ( sin( beta_c12 )*(c12h(1,2) - y0)+(c12h(1,3) - z0)* sin( gamma_c12 ) ) + a_c12 * ( ((c12h(1,2) - y0))^2 * ( cos( gamma_c12 ) )^2+ ((c12h(1,3) - z0))^2 * ( cos( beta_c12 ) )^2+ 2*(c12h(1,2) - y0)*(c12h(1,3) - z0)* sin( beta_c12 )	* sin( gamma_c12 )  )  ;
a44_c13		= I_c13 *( ( sin( beta_c13 ) )^2 + ( sin( gamma_c13 ) )^2) +2*s_c13 * ( sin( beta_c13 )*(c13h(1,2) - y0)+(c13h(1,3) - z0)* sin( gamma_c13 ) ) + a_c13 * ( ((c13h(1,2) - y0))^2 * ( cos( gamma_c13 ) )^2+ ((c13h(1,3) - z0))^2 * ( cos( beta_c13 ) )^2+ 2*(c13h(1,2) - y0)*(c13h(1,3) - z0)* sin( beta_c13 )	* sin( gamma_c13 )  )  ;
a44_c14		= I_c14 *( ( sin( beta_c14 ) )^2 + ( sin( gamma_c14 ) )^2) +2*s_c14 * ( sin( beta_c14 )*(c14h(1,2) - y0)+(c14h(1,3) - z0)* sin( gamma_c14 ) ) + a_c14 * ( ((c14h(1,2) - y0))^2 * ( cos( gamma_c14 ) )^2+ ((c14h(1,3) - z0))^2 * ( cos( beta_c14 ) )^2+ 2*(c14h(1,2) - y0)*(c14h(1,3) - z0)* sin( beta_c14 )	* sin( gamma_c14 )  )  ;
a44_c15		= I_c15 *( ( sin( beta_c15 ) )^2 + ( sin( gamma_c15 ) )^2) +2*s_c15 * ( sin( beta_c15 )*(c15h(1,2) - y0)+(c15h(1,3) - z0)* sin( gamma_c15 ) ) + a_c15 * ( ((c15h(1,2) - y0))^2 * ( cos( gamma_c15 ) )^2+ ((c15h(1,3) - z0))^2 * ( cos( beta_c15 ) )^2+ 2*(c15h(1,2) - y0)*(c15h(1,3) - z0)* sin( beta_c15 )	* sin( gamma_c15 )  )  ;
a44_c16		= I_c16 *( ( sin( beta_c16 ) )^2 + ( sin( gamma_c16 ) )^2) +2*s_c16 * ( sin( beta_c16 )*(c16h(1,2) - y0)+(c16h(1,3) - z0)* sin( gamma_c16 ) ) + a_c16 * ( ((c16h(1,2) - y0))^2 * ( cos( gamma_c16 ) )^2+ ((c16h(1,3) - z0))^2 * ( cos( beta_c16 ) )^2+ 2*(c16h(1,2) - y0)*(c16h(1,3) - z0)* sin( beta_c16 )	* sin( gamma_c16 )  )  ;
a44_c17		= I_c17 *( ( sin( beta_c17 ) )^2 + ( sin( gamma_c17 ) )^2) +2*s_c17 * ( sin( beta_c17 )*(c17h(1,2) - y0)+(c17h(1,3) - z0)* sin( gamma_c17 ) ) + a_c17 * ( ((c17h(1,2) - y0))^2 * ( cos( gamma_c17 ) )^2+ ((c17h(1,3) - z0))^2 * ( cos( beta_c17 ) )^2+ 2*(c17h(1,2) - y0)*(c17h(1,3) - z0)* sin( beta_c17 )	* sin( gamma_c17 )  )  ;
a44_c18		= I_c18 *( ( sin( beta_c18 ) )^2 + ( sin( gamma_c18 ) )^2) +2*s_c18 * ( sin( beta_c18 )*(c18h(1,2) - y0)+(c18h(1,3) - z0)* sin( gamma_c18 ) ) + a_c18 * ( ((c18h(1,2) - y0))^2 * ( cos( gamma_c18 ) )^2+ ((c18h(1,3) - z0))^2 * ( cos( beta_c18 ) )^2+ 2*(c18h(1,2) - y0)*(c18h(1,3) - z0)* sin( beta_c18 )	* sin( gamma_c18 )  )  ;
a44_c19		= I_c19 *( ( sin( beta_c19 ) )^2 + ( sin( gamma_c19 ) )^2) +2*s_c19 * ( sin( beta_c19 )*(c19h(1,2) - y0)+(c19h(1,3) - z0)* sin( gamma_c19 ) ) + a_c19 * ( ((c19h(1,2) - y0))^2 * ( cos( gamma_c19 ) )^2+ ((c19h(1,3) - z0))^2 * ( cos( beta_c19 ) )^2+ 2*(c19h(1,2) - y0)*(c19h(1,3) - z0)* sin( beta_c19 )	* sin( gamma_c19 )  )  ;
a44_c20		= I_c20 *( ( sin( beta_c20 ) )^2 + ( sin( gamma_c20 ) )^2) +2*s_c20 * ( sin( beta_c20 )*(c20h(1,2) - y0)+(c20h(1,3) - z0)* sin( gamma_c20 ) ) + a_c20 * ( ((c20h(1,2) - y0))^2 * ( cos( gamma_c20 ) )^2+ ((c20h(1,3) - z0))^2 * ( cos( beta_c20 ) )^2+ 2*(c20h(1,2) - y0)*(c20h(1,3) - z0)* sin( beta_c20 )	* sin( gamma_c20 )  )  ;
a44_p1 		= I_p1	*( ( sin( beta_p1  ) )^2 + ( sin( gamma_p1  ) )^2) +2*s_p1  * ( sin( beta_p1  )*(p1a(1,2)  - y0)+(p1a(1,3)  - z0)* sin( gamma_p1  ) ) + a_p1  * ( ((p1a(1,2)  - y0))^2 * ( cos( gamma_p1  ) )^2+ ((p1a(1,3)  - z0))^2 * ( cos( beta_p1  ) )^2+ 2*(p1a(1,2)  - y0)*(p1a(1,3)  - z0)* sin( beta_p1  )	* sin( gamma_p1  )  )  ;
a44_p2 		= I_p2	*( ( sin( beta_p2  ) )^2 + ( sin( gamma_p2  ) )^2) +2*s_p2  * ( sin( beta_p2  )*(p2a(1,2)  - y0)+(p2a(1,3)  - z0)* sin( gamma_p2  ) ) + a_p2  * ( ((p2a(1,2)  - y0))^2 * ( cos( gamma_p2  ) )^2+ ((p2a(1,3)  - z0))^2 * ( cos( beta_p2  ) )^2+ 2*(p2a(1,2)  - y0)*(p2a(1,3)  - z0)* sin( beta_p2  )	* sin( gamma_p2  )  )  ;
a44_p3 		= I_p3	*( ( sin( beta_p3  ) )^2 + ( sin( gamma_p3  ) )^2) +2*s_p3  * ( sin( beta_p3  )*(p3a(1,2)  - y0)+(p3a(1,3)  - z0)* sin( gamma_p3  ) ) + a_p3  * ( ((p3a(1,2)  - y0))^2 * ( cos( gamma_p3  ) )^2+ ((p3a(1,3)  - z0))^2 * ( cos( beta_p3  ) )^2+ 2*(p3a(1,2)  - y0)*(p3a(1,3)  - z0)* sin( beta_p3  )	* sin( gamma_p3  )  )  ;
a44_p4 		= I_p4	*( ( sin( beta_p4  ) )^2 + ( sin( gamma_p4  ) )^2) +2*s_p4  * ( sin( beta_p4  )*(p4a(1,2)  - y0)+(p4a(1,3)  - z0)* sin( gamma_p4  ) ) + a_p4  * ( ((p4a(1,2)  - y0))^2 * ( cos( gamma_p4  ) )^2+ ((p4a(1,3)  - z0))^2 * ( cos( beta_p4  ) )^2+ 2*(p4a(1,2)  - y0)*(p4a(1,3)  - z0)* sin( beta_p4  )	* sin( gamma_p4  )  )  ;
a44_b1 		= I_b1	*( ( sin( beta_b1  ) )^2 + ( sin( gamma_b1  ) )^2) +2*s_b1  * ( sin( beta_b1  )*(b1s(1,2)  - y0)+(b1s(1,3)  - z0)* sin( gamma_b1  ) ) + a_b1  * ( ((b1s(1,2)  - y0))^2 * ( cos( gamma_b1  ) )^2+ ((b1s(1,3)  - z0))^2 * ( cos( beta_b1  ) )^2+ 2*(b1s(1,2)  - y0)*(b1s(1,3)  - z0)* sin( beta_b1  )	* sin( gamma_b1  )  )  ;
a44_b2 		= I_b2	*( ( sin( beta_b2  ) )^2 + ( sin( gamma_b2  ) )^2) +2*s_b2  * ( sin( beta_b2  )*(b2s(1,2)  - y0)+(b2s(1,3)  - z0)* sin( gamma_b2  ) ) + a_b2  * ( ((b2s(1,2)  - y0))^2 * ( cos( gamma_b2  ) )^2+ ((b2s(1,3)  - z0))^2 * ( cos( beta_b2  ) )^2+ 2*(b2s(1,2)  - y0)*(b2s(1,3)  - z0)* sin( beta_b2  )	* sin( gamma_b2  )  )  ;
a44_b3 		= I_b3	*( ( sin( beta_b3  ) )^2 + ( sin( gamma_b3  ) )^2) +2*s_b3  * ( sin( beta_b3  )*(b3s(1,2)  - y0)+(b3s(1,3)  - z0)* sin( gamma_b3  ) ) + a_b3  * ( ((b3s(1,2)  - y0))^2 * ( cos( gamma_b3  ) )^2+ ((b3s(1,3)  - z0))^2 * ( cos( beta_b3  ) )^2+ 2*(b3s(1,2)  - y0)*(b3s(1,3)  - z0)* sin( beta_b3  )	* sin( gamma_b3  )  )  ;
a44_b4 		= I_b4	*( ( sin( beta_b4  ) )^2 + ( sin( gamma_b4  ) )^2) +2*s_b4  * ( sin( beta_b4  )*(b4s(1,2)  - y0)+(b4s(1,3)  - z0)* sin( gamma_b4  ) ) + a_b4  * ( ((b4s(1,2)  - y0))^2 * ( cos( gamma_b4  ) )^2+ ((b4s(1,3)  - z0))^2 * ( cos( beta_b4  ) )^2+ 2*(b4s(1,2)  - y0)*(b4s(1,3)  - z0)* sin( beta_b4  )	* sin( gamma_b4  )  )  ;
  
a44(1,1) = a44_c1 		;
a44(1,2) = a44_c2 		;
a44(1,3) = a44_c3 		;
a44(1,4) = a44_c4 		;
a44(1,5) = a44_c5 		;
a44(1,6) = a44_c6 		;
a44(1,7) = a44_c7 		;
a44(1,8) = a44_c8 		;
a44(1,9) = a44_c9 		;
a44(1,10)= a44_c10		;
a44(1,11)= a44_c11		;
a44(1,12)= a44_c12		;
a44(1,13)= a44_c13		;
a44(1,14)= a44_c14		;
a44(1,15)= a44_c15		;
a44(1,16)= a44_c16		;
a44(1,17)= a44_c17		;
a44(1,18)= a44_c18		;
a44(1,19)= a44_c19		;
a44(1,20)= a44_c20		;
a44(1,21)= a44_p1 		;
a44(1,22)= a44_p2 		;
a44(1,23)= a44_p3 		;
a44(1,24)= a44_p4 		;
a44(1,25)= a44_b1 		;
a44(1,26)= a44_b2 		;
a44(1,27)= a44_b3 		;
a44(1,28)= a44_b4 		;

a55_c1 		= I_c1  *( ( sin( alpha_c1  ) )^2 + ( sin( gamma_c1  ) )^2) +2*s_c1  * ( sin( alpha_c1  )*(c1h(1,1)  - x0)+(c1h(1,3)  - z0)* sin( gamma_c1  ) ) + a_c1  * ( ((c1h(1,1)  - x0))^2 * ( cos( gamma_c1  ) )^2+ ((c1h(1,3)  - z0))^2 * ( cos( alpha_c1  ) )^2+ 2*(c1h(1,1)  - x0)*(c1h(1,3)  - z0)* sin( alpha_c1  )	* sin( gamma_c1  )  )  ;
a55_c2 		= I_c2  *( ( sin( alpha_c2  ) )^2 + ( sin( gamma_c2  ) )^2) +2*s_c2  * ( sin( alpha_c2  )*(c2h(1,1)  - x0)+(c2h(1,3)  - z0)* sin( gamma_c2  ) ) + a_c2  * ( ((c2h(1,1)  - x0))^2 * ( cos( gamma_c2  ) )^2+ ((c2h(1,3)  - z0))^2 * ( cos( alpha_c2  ) )^2+ 2*(c2h(1,1)  - x0)*(c2h(1,3)  - z0)* sin( alpha_c2  )	* sin( gamma_c2  )  )  ;
a55_c3 		= I_c3  *( ( sin( alpha_c3  ) )^2 + ( sin( gamma_c3  ) )^2) +2*s_c3  * ( sin( alpha_c3  )*(c3h(1,1)  - x0)+(c3h(1,3)  - z0)* sin( gamma_c3  ) ) + a_c3  * ( ((c3h(1,1)  - x0))^2 * ( cos( gamma_c3  ) )^2+ ((c3h(1,3)  - z0))^2 * ( cos( alpha_c3  ) )^2+ 2*(c3h(1,1)  - x0)*(c3h(1,3)  - z0)* sin( alpha_c3  )	* sin( gamma_c3  )  )  ;
a55_c4 		= I_c4  *( ( sin( alpha_c4  ) )^2 + ( sin( gamma_c4  ) )^2) +2*s_c4  * ( sin( alpha_c4  )*(c4h(1,1)  - x0)+(c4h(1,3)  - z0)* sin( gamma_c4  ) ) + a_c4  * ( ((c4h(1,1)  - x0))^2 * ( cos( gamma_c4  ) )^2+ ((c4h(1,3)  - z0))^2 * ( cos( alpha_c4  ) )^2+ 2*(c4h(1,1)  - x0)*(c4h(1,3)  - z0)* sin( alpha_c4  )	* sin( gamma_c4  )  )  ;
a55_c5 		= I_c5  *( ( sin( alpha_c5  ) )^2 + ( sin( gamma_c5  ) )^2) +2*s_c5  * ( sin( alpha_c5  )*(c5h(1,1)  - x0)+(c5h(1,3)  - z0)* sin( gamma_c5  ) ) + a_c5  * ( ((c5h(1,1)  - x0))^2 * ( cos( gamma_c5  ) )^2+ ((c5h(1,3)  - z0))^2 * ( cos( alpha_c5  ) )^2+ 2*(c5h(1,1)  - x0)*(c5h(1,3)  - z0)* sin( alpha_c5  )	* sin( gamma_c5  )  )  ;
a55_c6 		= I_c6  *( ( sin( alpha_c6  ) )^2 + ( sin( gamma_c6  ) )^2) +2*s_c6  * ( sin( alpha_c6  )*(c6h(1,1)  - x0)+(c6h(1,3)  - z0)* sin( gamma_c6  ) ) + a_c6  * ( ((c6h(1,1)  - x0))^2 * ( cos( gamma_c6  ) )^2+ ((c6h(1,3)  - z0))^2 * ( cos( alpha_c6  ) )^2+ 2*(c6h(1,1)  - x0)*(c6h(1,3)  - z0)* sin( alpha_c6  )	* sin( gamma_c6  )  )  ;
a55_c7 		= I_c7  *( ( sin( alpha_c7  ) )^2 + ( sin( gamma_c7  ) )^2) +2*s_c7  * ( sin( alpha_c7  )*(c7h(1,1)  - x0)+(c7h(1,3)  - z0)* sin( gamma_c7  ) ) + a_c7  * ( ((c7h(1,1)  - x0))^2 * ( cos( gamma_c7  ) )^2+ ((c7h(1,3)  - z0))^2 * ( cos( alpha_c7  ) )^2+ 2*(c7h(1,1)  - x0)*(c7h(1,3)  - z0)* sin( alpha_c7  )	* sin( gamma_c7  )  )  ;
a55_c8 		= I_c8  *( ( sin( alpha_c8  ) )^2 + ( sin( gamma_c8  ) )^2) +2*s_c8  * ( sin( alpha_c8  )*(c8h(1,1)  - x0)+(c8h(1,3)  - z0)* sin( gamma_c8  ) ) + a_c8  * ( ((c8h(1,1)  - x0))^2 * ( cos( gamma_c8  ) )^2+ ((c8h(1,3)  - z0))^2 * ( cos( alpha_c8  ) )^2+ 2*(c8h(1,1)  - x0)*(c8h(1,3)  - z0)* sin( alpha_c8  )	* sin( gamma_c8  )  )  ;
a55_c9 		= I_c9  *( ( sin( alpha_c9  ) )^2 + ( sin( gamma_c9  ) )^2) +2*s_c9  * ( sin( alpha_c9  )*(c9h(1,1)  - x0)+(c9h(1,3)  - z0)* sin( gamma_c9  ) ) + a_c9  * ( ((c9h(1,1)  - x0))^2 * ( cos( gamma_c9  ) )^2+ ((c9h(1,3)  - z0))^2 * ( cos( alpha_c9  ) )^2+ 2*(c9h(1,1)  - x0)*(c9h(1,3)  - z0)* sin( alpha_c9  )	* sin( gamma_c9  )  )  ;
a55_c10		= I_c10 *( ( sin( alpha_c10 ) )^2 + ( sin( gamma_c10 ) )^2) +2*s_c10 * ( sin( alpha_c10 )*(c10h(1,1) - x0)+(c10h(1,3) - z0)* sin( gamma_c10 ) ) + a_c10 * ( ((c10h(1,1) - x0))^2 * ( cos( gamma_c10 ) )^2+ ((c10h(1,3) - z0))^2 * ( cos( alpha_c10 ) )^2+ 2*(c10h(1,1) - x0)*(c10h(1,3) - z0)* sin( alpha_c10 )	* sin( gamma_c10 )  )  ;
a55_c11		= I_c11 *( ( sin( alpha_c11 ) )^2 + ( sin( gamma_c11 ) )^2) +2*s_c11 * ( sin( alpha_c11 )*(c11h(1,1) - x0)+(c11h(1,3) - z0)* sin( gamma_c11 ) ) + a_c11 * ( ((c11h(1,1) - x0))^2 * ( cos( gamma_c11 ) )^2+ ((c11h(1,3) - z0))^2 * ( cos( alpha_c11 ) )^2+ 2*(c11h(1,1) - x0)*(c11h(1,3) - z0)* sin( alpha_c11 )	* sin( gamma_c11 )  )  ;
a55_c12		= I_c12 *( ( sin( alpha_c12 ) )^2 + ( sin( gamma_c12 ) )^2) +2*s_c12 * ( sin( alpha_c12 )*(c12h(1,1) - x0)+(c12h(1,3) - z0)* sin( gamma_c12 ) ) + a_c12 * ( ((c12h(1,1) - x0))^2 * ( cos( gamma_c12 ) )^2+ ((c12h(1,3) - z0))^2 * ( cos( alpha_c12 ) )^2+ 2*(c12h(1,1) - x0)*(c12h(1,3) - z0)* sin( alpha_c12 )	* sin( gamma_c12 )  )  ;
a55_c13		= I_c13 *( ( sin( alpha_c13 ) )^2 + ( sin( gamma_c13 ) )^2) +2*s_c13 * ( sin( alpha_c13 )*(c13h(1,1) - x0)+(c13h(1,3) - z0)* sin( gamma_c13 ) ) + a_c13 * ( ((c13h(1,1) - x0))^2 * ( cos( gamma_c13 ) )^2+ ((c13h(1,3) - z0))^2 * ( cos( alpha_c13 ) )^2+ 2*(c13h(1,1) - x0)*(c13h(1,3) - z0)* sin( alpha_c13 )	* sin( gamma_c13 )  )  ;
a55_c14		= I_c14 *( ( sin( alpha_c14 ) )^2 + ( sin( gamma_c14 ) )^2) +2*s_c14 * ( sin( alpha_c14 )*(c14h(1,1) - x0)+(c14h(1,3) - z0)* sin( gamma_c14 ) ) + a_c14 * ( ((c14h(1,1) - x0))^2 * ( cos( gamma_c14 ) )^2+ ((c14h(1,3) - z0))^2 * ( cos( alpha_c14 ) )^2+ 2*(c14h(1,1) - x0)*(c14h(1,3) - z0)* sin( alpha_c14 )	* sin( gamma_c14 )  )  ;
a55_c15		= I_c15 *( ( sin( alpha_c15 ) )^2 + ( sin( gamma_c15 ) )^2) +2*s_c15 * ( sin( alpha_c15 )*(c15h(1,1) - x0)+(c15h(1,3) - z0)* sin( gamma_c15 ) ) + a_c15 * ( ((c15h(1,1) - x0))^2 * ( cos( gamma_c15 ) )^2+ ((c15h(1,3) - z0))^2 * ( cos( alpha_c15 ) )^2+ 2*(c15h(1,1) - x0)*(c15h(1,3) - z0)* sin( alpha_c15 )	* sin( gamma_c15 )  )  ;
a55_c16		= I_c16 *( ( sin( alpha_c16 ) )^2 + ( sin( gamma_c16 ) )^2) +2*s_c16 * ( sin( alpha_c16 )*(c16h(1,1) - x0)+(c16h(1,3) - z0)* sin( gamma_c16 ) ) + a_c16 * ( ((c16h(1,1) - x0))^2 * ( cos( gamma_c16 ) )^2+ ((c16h(1,3) - z0))^2 * ( cos( alpha_c16 ) )^2+ 2*(c16h(1,1) - x0)*(c16h(1,3) - z0)* sin( alpha_c16 )	* sin( gamma_c16 )  )  ;
a55_c17		= I_c17 *( ( sin( alpha_c17 ) )^2 + ( sin( gamma_c17 ) )^2) +2*s_c17 * ( sin( alpha_c17 )*(c17h(1,1) - x0)+(c17h(1,3) - z0)* sin( gamma_c17 ) ) + a_c17 * ( ((c17h(1,1) - x0))^2 * ( cos( gamma_c17 ) )^2+ ((c17h(1,3) - z0))^2 * ( cos( alpha_c17 ) )^2+ 2*(c17h(1,1) - x0)*(c17h(1,3) - z0)* sin( alpha_c17 )	* sin( gamma_c17 )  )  ;
a55_c18		= I_c18 *( ( sin( alpha_c18 ) )^2 + ( sin( gamma_c18 ) )^2) +2*s_c18 * ( sin( alpha_c18 )*(c18h(1,1) - x0)+(c18h(1,3) - z0)* sin( gamma_c18 ) ) + a_c18 * ( ((c18h(1,1) - x0))^2 * ( cos( gamma_c18 ) )^2+ ((c18h(1,3) - z0))^2 * ( cos( alpha_c18 ) )^2+ 2*(c18h(1,1) - x0)*(c18h(1,3) - z0)* sin( alpha_c18 )	* sin( gamma_c18 )  )  ;
a55_c19		= I_c19 *( ( sin( alpha_c19 ) )^2 + ( sin( gamma_c19 ) )^2) +2*s_c19 * ( sin( alpha_c19 )*(c19h(1,1) - x0)+(c19h(1,3) - z0)* sin( gamma_c19 ) ) + a_c19 * ( ((c19h(1,1) - x0))^2 * ( cos( gamma_c19 ) )^2+ ((c19h(1,3) - z0))^2 * ( cos( alpha_c19 ) )^2+ 2*(c19h(1,1) - x0)*(c19h(1,3) - z0)* sin( alpha_c19 )	* sin( gamma_c19 )  )  ;
a55_c20		= I_c20 *( ( sin( alpha_c20 ) )^2 + ( sin( gamma_c20 ) )^2) +2*s_c20 * ( sin( alpha_c20 )*(c20h(1,1) - x0)+(c20h(1,3) - z0)* sin( gamma_c20 ) ) + a_c20 * ( ((c20h(1,1) - x0))^2 * ( cos( gamma_c20 ) )^2+ ((c20h(1,3) - z0))^2 * ( cos( alpha_c20 ) )^2+ 2*(c20h(1,1) - x0)*(c20h(1,3) - z0)* sin( alpha_c20 )	* sin( gamma_c20 )  )  ;
a55_p1 		= I_p1	*( ( sin( alpha_p1  ) )^2 + ( sin( gamma_p1  ) )^2) +2*s_p1  * ( sin( alpha_p1  )*(p1a(1,1)  - x0)+(p1a(1,3)  - z0)* sin( gamma_p1  ) ) + a_p1  * ( ((p1a(1,1)  - x0))^2 * ( cos( gamma_p1  ) )^2+ ((p1a(1,3)  - z0))^2 * ( cos( alpha_p1  ) )^2+ 2*(p1a(1,1)  - x0)*(p1a(1,3)  - z0)* sin( alpha_p1  )	* sin( gamma_p1  )  )  ;
a55_p2 		= I_p2	*( ( sin( alpha_p2  ) )^2 + ( sin( gamma_p2  ) )^2) +2*s_p2  * ( sin( alpha_p2  )*(p2a(1,1)  - x0)+(p2a(1,3)  - z0)* sin( gamma_p2  ) ) + a_p2  * ( ((p2a(1,1)  - x0))^2 * ( cos( gamma_p2  ) )^2+ ((p2a(1,3)  - z0))^2 * ( cos( alpha_p2  ) )^2+ 2*(p2a(1,1)  - x0)*(p2a(1,3)  - z0)* sin( alpha_p2  )	* sin( gamma_p2  )  )  ;
a55_p3 		= I_p3	*( ( sin( alpha_p3  ) )^2 + ( sin( gamma_p3  ) )^2) +2*s_p3  * ( sin( alpha_p3  )*(p3a(1,1)  - x0)+(p3a(1,3)  - z0)* sin( gamma_p3  ) ) + a_p3  * ( ((p3a(1,1)  - x0))^2 * ( cos( gamma_p3  ) )^2+ ((p3a(1,3)  - z0))^2 * ( cos( alpha_p3  ) )^2+ 2*(p3a(1,1)  - x0)*(p3a(1,3)  - z0)* sin( alpha_p3  )	* sin( gamma_p3  )  )  ;
a55_p4 		= I_p4	*( ( sin( alpha_p4  ) )^2 + ( sin( gamma_p4  ) )^2) +2*s_p4  * ( sin( alpha_p4  )*(p4a(1,1)  - x0)+(p4a(1,3)  - z0)* sin( gamma_p4  ) ) + a_p4  * ( ((p4a(1,1)  - x0))^2 * ( cos( gamma_p4  ) )^2+ ((p4a(1,3)  - z0))^2 * ( cos( alpha_p4  ) )^2+ 2*(p4a(1,1)  - x0)*(p4a(1,3)  - z0)* sin( alpha_p4  )	* sin( gamma_p4  )  )  ;
a55_b1 		= I_b1	*( ( sin( alpha_b1  ) )^2 + ( sin( gamma_b1  ) )^2) +2*s_b1  * ( sin( alpha_b1  )*(b1s(1,1)  - x0)+(b1s(1,3)  - z0)* sin( gamma_b1  ) ) + a_b1  * ( ((b1s(1,1)  - x0))^2 * ( cos( gamma_b1  ) )^2+ ((b1s(1,3)  - z0))^2 * ( cos( alpha_b1  ) )^2+ 2*(b1s(1,1)  - x0)*(b1s(1,3)  - z0)* sin( alpha_b1  )	* sin( gamma_b1  )  )  ;
a55_b2 		= I_b2	*( ( sin( alpha_b2  ) )^2 + ( sin( gamma_b2  ) )^2) +2*s_b2  * ( sin( alpha_b2  )*(b2s(1,1)  - x0)+(b2s(1,3)  - z0)* sin( gamma_b2  ) ) + a_b2  * ( ((b2s(1,1)  - x0))^2 * ( cos( gamma_b2  ) )^2+ ((b2s(1,3)  - z0))^2 * ( cos( alpha_b2  ) )^2+ 2*(b2s(1,1)  - x0)*(b2s(1,3)  - z0)* sin( alpha_b2  )	* sin( gamma_b2  )  )  ;
a55_b3 		= I_b3	*( ( sin( alpha_b3  ) )^2 + ( sin( gamma_b3  ) )^2) +2*s_b3  * ( sin( alpha_b3  )*(b3s(1,1)  - x0)+(b3s(1,3)  - z0)* sin( gamma_b3  ) ) + a_b3  * ( ((b3s(1,1)  - x0))^2 * ( cos( gamma_b3  ) )^2+ ((b3s(1,3)  - z0))^2 * ( cos( alpha_b3  ) )^2+ 2*(b3s(1,1)  - x0)*(b3s(1,3)  - z0)* sin( alpha_b3  )	* sin( gamma_b3  )  )  ;
a55_b4 		= I_b4	*( ( sin( alpha_b4  ) )^2 + ( sin( gamma_b4  ) )^2) +2*s_b4  * ( sin( alpha_b4  )*(b4s(1,1)  - x0)+(b4s(1,3)  - z0)* sin( gamma_b4  ) ) + a_b4  * ( ((b4s(1,1)  - x0))^2 * ( cos( gamma_b4  ) )^2+ ((b4s(1,3)  - z0))^2 * ( cos( alpha_b4  ) )^2+ 2*(b4s(1,1)  - x0)*(b4s(1,3)  - z0)* sin( alpha_b4  )	* sin( gamma_b4  )  )  ;
  
a55(1,1) = a55_c1 		;
a55(1,2) = a55_c2 		;
a55(1,3) = a55_c3 		;
a55(1,4) = a55_c4 		;
a55(1,5) = a55_c5 		;
a55(1,6) = a55_c6 		;
a55(1,7) = a55_c7 		;
a55(1,8) = a55_c8 		;
a55(1,9) = a55_c9 		;
a55(1,10)= a55_c10		;
a55(1,11)= a55_c11		;
a55(1,12)= a55_c12		;
a55(1,13)= a55_c13		;
a55(1,14)= a55_c14		;
a55(1,15)= a55_c15		;
a55(1,16)= a55_c16		;
a55(1,17)= a55_c17		;
a55(1,18)= a55_c18		;
a55(1,19)= a55_c19		;
a55(1,20)= a55_c20		;
a55(1,21)= a55_p1 		;
a55(1,22)= a55_p2 		;
a55(1,23)= a55_p3 		;
a55(1,24)= a55_p4 		;
a55(1,25)= a55_b1 		;
a55(1,26)= a55_b2 		;
a55(1,27)= a55_b3 		;
a55(1,28)= a55_b4 		;

a66_c1 		= I_c1  *( ( sin( alpha_c1  ) )^2 + ( sin( beta_c1  ) )^2) +2*s_c1  * ( sin( alpha_c1  )*(c1h(1,1)  - x0)+(c1h(1,2)  - y0)* sin( beta_c1  ) ) + a_c1  * ( ((c1h(1,1)  - x0))^2 * ( cos( beta_c1  ) )^2+ ((c1h(1,2)  - y0))^2 * ( cos( alpha_c1  ) )^2+ 2*(c1h(1,1)  - x0)*(c1h(1,2)  - y0)* sin( alpha_c1  )	* sin( beta_c1  )  )  ;
a66_c2 		= I_c2  *( ( sin( alpha_c2  ) )^2 + ( sin( beta_c2  ) )^2) +2*s_c2  * ( sin( alpha_c2  )*(c2h(1,1)  - x0)+(c2h(1,2)  - y0)* sin( beta_c2  ) ) + a_c2  * ( ((c2h(1,1)  - x0))^2 * ( cos( beta_c2  ) )^2+ ((c2h(1,2)  - y0))^2 * ( cos( alpha_c2  ) )^2+ 2*(c2h(1,1)  - x0)*(c2h(1,2)  - y0)* sin( alpha_c2  )	* sin( beta_c2  )  )  ;
a66_c3 		= I_c3  *( ( sin( alpha_c3  ) )^2 + ( sin( beta_c3  ) )^2) +2*s_c3  * ( sin( alpha_c3  )*(c3h(1,1)  - x0)+(c3h(1,2)  - y0)* sin( beta_c3  ) ) + a_c3  * ( ((c3h(1,1)  - x0))^2 * ( cos( beta_c3  ) )^2+ ((c3h(1,2)  - y0))^2 * ( cos( alpha_c3  ) )^2+ 2*(c3h(1,1)  - x0)*(c3h(1,2)  - y0)* sin( alpha_c3  )	* sin( beta_c3  )  )  ;
a66_c4 		= I_c4  *( ( sin( alpha_c4  ) )^2 + ( sin( beta_c4  ) )^2) +2*s_c4  * ( sin( alpha_c4  )*(c4h(1,1)  - x0)+(c4h(1,2)  - y0)* sin( beta_c4  ) ) + a_c4  * ( ((c4h(1,1)  - x0))^2 * ( cos( beta_c4  ) )^2+ ((c4h(1,2)  - y0))^2 * ( cos( alpha_c4  ) )^2+ 2*(c4h(1,1)  - x0)*(c4h(1,2)  - y0)* sin( alpha_c4  )	* sin( beta_c4  )  )  ;
a66_c5 		= I_c5  *( ( sin( alpha_c5  ) )^2 + ( sin( beta_c5  ) )^2) +2*s_c5  * ( sin( alpha_c5  )*(c5h(1,1)  - x0)+(c5h(1,2)  - y0)* sin( beta_c5  ) ) + a_c5  * ( ((c5h(1,1)  - x0))^2 * ( cos( beta_c5  ) )^2+ ((c5h(1,2)  - y0))^2 * ( cos( alpha_c5  ) )^2+ 2*(c5h(1,1)  - x0)*(c5h(1,2)  - y0)* sin( alpha_c5  )	* sin( beta_c5  )  )  ;
a66_c6 		= I_c6  *( ( sin( alpha_c6  ) )^2 + ( sin( beta_c6  ) )^2) +2*s_c6  * ( sin( alpha_c6  )*(c6h(1,1)  - x0)+(c6h(1,2)  - y0)* sin( beta_c6  ) ) + a_c6  * ( ((c6h(1,1)  - x0))^2 * ( cos( beta_c6  ) )^2+ ((c6h(1,2)  - y0))^2 * ( cos( alpha_c6  ) )^2+ 2*(c6h(1,1)  - x0)*(c6h(1,2)  - y0)* sin( alpha_c6  )	* sin( beta_c6  )  )  ;
a66_c7 		= I_c7  *( ( sin( alpha_c7  ) )^2 + ( sin( beta_c7  ) )^2) +2*s_c7  * ( sin( alpha_c7  )*(c7h(1,1)  - x0)+(c7h(1,2)  - y0)* sin( beta_c7  ) ) + a_c7  * ( ((c7h(1,1)  - x0))^2 * ( cos( beta_c7  ) )^2+ ((c7h(1,2)  - y0))^2 * ( cos( alpha_c7  ) )^2+ 2*(c7h(1,1)  - x0)*(c7h(1,2)  - y0)* sin( alpha_c7  )	* sin( beta_c7  )  )  ;
a66_c8 		= I_c8  *( ( sin( alpha_c8  ) )^2 + ( sin( beta_c8  ) )^2) +2*s_c8  * ( sin( alpha_c8  )*(c8h(1,1)  - x0)+(c8h(1,2)  - y0)* sin( beta_c8  ) ) + a_c8  * ( ((c8h(1,1)  - x0))^2 * ( cos( beta_c8  ) )^2+ ((c8h(1,2)  - y0))^2 * ( cos( alpha_c8  ) )^2+ 2*(c8h(1,1)  - x0)*(c8h(1,2)  - y0)* sin( alpha_c8  )	* sin( beta_c8  )  )  ;
a66_c9 		= I_c9  *( ( sin( alpha_c9  ) )^2 + ( sin( beta_c9  ) )^2) +2*s_c9  * ( sin( alpha_c9  )*(c9h(1,1)  - x0)+(c9h(1,2)  - y0)* sin( beta_c9  ) ) + a_c9  * ( ((c9h(1,1)  - x0))^2 * ( cos( beta_c9  ) )^2+ ((c9h(1,2)  - y0))^2 * ( cos( alpha_c9  ) )^2+ 2*(c9h(1,1)  - x0)*(c9h(1,2)  - y0)* sin( alpha_c9  )	* sin( beta_c9  )  )  ;
a66_c10		= I_c10 *( ( sin( alpha_c10 ) )^2 + ( sin( beta_c10 ) )^2) +2*s_c10 * ( sin( alpha_c10 )*(c10h(1,1) - x0)+(c10h(1,2) - y0)* sin( beta_c10 ) ) + a_c10 * ( ((c10h(1,1) - x0))^2 * ( cos( beta_c10 ) )^2+ ((c10h(1,2) - y0))^2 * ( cos( alpha_c10 ) )^2+ 2*(c10h(1,1) - x0)*(c10h(1,2) - y0)* sin( alpha_c10 )	* sin( beta_c10 )  )  ;
a66_c11		= I_c11 *( ( sin( alpha_c11 ) )^2 + ( sin( beta_c11 ) )^2) +2*s_c11 * ( sin( alpha_c11 )*(c11h(1,1) - x0)+(c11h(1,2) - y0)* sin( beta_c11 ) ) + a_c11 * ( ((c11h(1,1) - x0))^2 * ( cos( beta_c11 ) )^2+ ((c11h(1,2) - y0))^2 * ( cos( alpha_c11 ) )^2+ 2*(c11h(1,1) - x0)*(c11h(1,2) - y0)* sin( alpha_c11 )	* sin( beta_c11 )  )  ;
a66_c12		= I_c12 *( ( sin( alpha_c12 ) )^2 + ( sin( beta_c12 ) )^2) +2*s_c12 * ( sin( alpha_c12 )*(c12h(1,1) - x0)+(c12h(1,2) - y0)* sin( beta_c12 ) ) + a_c12 * ( ((c12h(1,1) - x0))^2 * ( cos( beta_c12 ) )^2+ ((c12h(1,2) - y0))^2 * ( cos( alpha_c12 ) )^2+ 2*(c12h(1,1) - x0)*(c12h(1,2) - y0)* sin( alpha_c12 )	* sin( beta_c12 )  )  ;
a66_c13		= I_c13 *( ( sin( alpha_c13 ) )^2 + ( sin( beta_c13 ) )^2) +2*s_c13 * ( sin( alpha_c13 )*(c13h(1,1) - x0)+(c13h(1,2) - y0)* sin( beta_c13 ) ) + a_c13 * ( ((c13h(1,1) - x0))^2 * ( cos( beta_c13 ) )^2+ ((c13h(1,2) - y0))^2 * ( cos( alpha_c13 ) )^2+ 2*(c13h(1,1) - x0)*(c13h(1,2) - y0)* sin( alpha_c13 )	* sin( beta_c13 )  )  ;
a66_c14		= I_c14 *( ( sin( alpha_c14 ) )^2 + ( sin( beta_c14 ) )^2) +2*s_c14 * ( sin( alpha_c14 )*(c14h(1,1) - x0)+(c14h(1,2) - y0)* sin( beta_c14 ) ) + a_c14 * ( ((c14h(1,1) - x0))^2 * ( cos( beta_c14 ) )^2+ ((c14h(1,2) - y0))^2 * ( cos( alpha_c14 ) )^2+ 2*(c14h(1,1) - x0)*(c14h(1,2) - y0)* sin( alpha_c14 )	* sin( beta_c14 )  )  ;
a66_c15		= I_c15 *( ( sin( alpha_c15 ) )^2 + ( sin( beta_c15 ) )^2) +2*s_c15 * ( sin( alpha_c15 )*(c15h(1,1) - x0)+(c15h(1,2) - y0)* sin( beta_c15 ) ) + a_c15 * ( ((c15h(1,1) - x0))^2 * ( cos( beta_c15 ) )^2+ ((c15h(1,2) - y0))^2 * ( cos( alpha_c15 ) )^2+ 2*(c15h(1,1) - x0)*(c15h(1,2) - y0)* sin( alpha_c15 )	* sin( beta_c15 )  )  ;
a66_c16		= I_c16 *( ( sin( alpha_c16 ) )^2 + ( sin( beta_c16 ) )^2) +2*s_c16 * ( sin( alpha_c16 )*(c16h(1,1) - x0)+(c16h(1,2) - y0)* sin( beta_c16 ) ) + a_c16 * ( ((c16h(1,1) - x0))^2 * ( cos( beta_c16 ) )^2+ ((c16h(1,2) - y0))^2 * ( cos( alpha_c16 ) )^2+ 2*(c16h(1,1) - x0)*(c16h(1,2) - y0)* sin( alpha_c16 )	* sin( beta_c16 )  )  ;
a66_c17		= I_c17 *( ( sin( alpha_c17 ) )^2 + ( sin( beta_c17 ) )^2) +2*s_c17 * ( sin( alpha_c17 )*(c17h(1,1) - x0)+(c17h(1,2) - y0)* sin( beta_c17 ) ) + a_c17 * ( ((c17h(1,1) - x0))^2 * ( cos( beta_c17 ) )^2+ ((c17h(1,2) - y0))^2 * ( cos( alpha_c17 ) )^2+ 2*(c17h(1,1) - x0)*(c17h(1,2) - y0)* sin( alpha_c17 )	* sin( beta_c17 )  )  ;
a66_c18		= I_c18 *( ( sin( alpha_c18 ) )^2 + ( sin( beta_c18 ) )^2) +2*s_c18 * ( sin( alpha_c18 )*(c18h(1,1) - x0)+(c18h(1,2) - y0)* sin( beta_c18 ) ) + a_c18 * ( ((c18h(1,1) - x0))^2 * ( cos( beta_c18 ) )^2+ ((c18h(1,2) - y0))^2 * ( cos( alpha_c18 ) )^2+ 2*(c18h(1,1) - x0)*(c18h(1,2) - y0)* sin( alpha_c18 )	* sin( beta_c18 )  )  ;
a66_c19		= I_c19 *( ( sin( alpha_c19 ) )^2 + ( sin( beta_c19 ) )^2) +2*s_c19 * ( sin( alpha_c19 )*(c19h(1,1) - x0)+(c19h(1,2) - y0)* sin( beta_c19 ) ) + a_c19 * ( ((c19h(1,1) - x0))^2 * ( cos( beta_c19 ) )^2+ ((c19h(1,2) - y0))^2 * ( cos( alpha_c19 ) )^2+ 2*(c19h(1,1) - x0)*(c19h(1,2) - y0)* sin( alpha_c19 )	* sin( beta_c19 )  )  ;
a66_c20		= I_c20 *( ( sin( alpha_c20 ) )^2 + ( sin( beta_c20 ) )^2) +2*s_c20 * ( sin( alpha_c20 )*(c20h(1,1) - x0)+(c20h(1,2) - y0)* sin( beta_c20 ) ) + a_c20 * ( ((c20h(1,1) - x0))^2 * ( cos( beta_c20 ) )^2+ ((c20h(1,2) - y0))^2 * ( cos( alpha_c20 ) )^2+ 2*(c20h(1,1) - x0)*(c20h(1,2) - y0)* sin( alpha_c20 )	* sin( beta_c20 )  )  ;
a66_p1 		= I_p1	*( ( sin( alpha_p1  ) )^2 + ( sin( beta_p1  ) )^2) +2*s_p1  * ( sin( alpha_p1  )*(p1a(1,1)  - x0)+(p1a(1,2)  - y0)* sin( beta_p1  ) ) + a_p1  * ( ((p1a(1,1)  - x0))^2 * ( cos( beta_p1  ) )^2+ ((p1a(1,2)  - y0))^2 * ( cos( alpha_p1  ) )^2+ 2*(p1a(1,1)  - x0)*(p1a(1,2)  - y0)* sin( alpha_p1  )	* sin( beta_p1  )  )  ;
a66_p2 		= I_p2	*( ( sin( alpha_p2  ) )^2 + ( sin( beta_p2  ) )^2) +2*s_p2  * ( sin( alpha_p2  )*(p2a(1,1)  - x0)+(p2a(1,2)  - y0)* sin( beta_p2  ) ) + a_p2  * ( ((p2a(1,1)  - x0))^2 * ( cos( beta_p2  ) )^2+ ((p2a(1,2)  - y0))^2 * ( cos( alpha_p2  ) )^2+ 2*(p2a(1,1)  - x0)*(p2a(1,2)  - y0)* sin( alpha_p2  )	* sin( beta_p2  )  )  ;
a66_p3 		= I_p3	*( ( sin( alpha_p3  ) )^2 + ( sin( beta_p3  ) )^2) +2*s_p3  * ( sin( alpha_p3  )*(p3a(1,1)  - x0)+(p3a(1,2)  - y0)* sin( beta_p3  ) ) + a_p3  * ( ((p3a(1,1)  - x0))^2 * ( cos( beta_p3  ) )^2+ ((p3a(1,2)  - y0))^2 * ( cos( alpha_p3  ) )^2+ 2*(p3a(1,1)  - x0)*(p3a(1,2)  - y0)* sin( alpha_p3  )	* sin( beta_p3  )  )  ;
a66_p4 		= I_p4	*( ( sin( alpha_p4  ) )^2 + ( sin( beta_p4  ) )^2) +2*s_p4  * ( sin( alpha_p4  )*(p4a(1,1)  - x0)+(p4a(1,2)  - y0)* sin( beta_p4  ) ) + a_p4  * ( ((p4a(1,1)  - x0))^2 * ( cos( beta_p4  ) )^2+ ((p4a(1,2)  - y0))^2 * ( cos( alpha_p4  ) )^2+ 2*(p4a(1,1)  - x0)*(p4a(1,2)  - y0)* sin( alpha_p4  )	* sin( beta_p4  )  )  ;
a66_b1 		= I_b1	*( ( sin( alpha_b1  ) )^2 + ( sin( beta_b1  ) )^2) +2*s_b1  * ( sin( alpha_b1  )*(b1s(1,1)  - x0)+(b1s(1,2)  - y0)* sin( beta_b1  ) ) + a_b1  * ( ((b1s(1,1)  - x0))^2 * ( cos( beta_b1  ) )^2+ ((b1s(1,2)  - y0))^2 * ( cos( alpha_b1  ) )^2+ 2*(b1s(1,1)  - x0)*(b1s(1,2)  - y0)* sin( alpha_b1  )	* sin( beta_b1  )  )  ;
a66_b2 		= I_b2	*( ( sin( alpha_b2  ) )^2 + ( sin( beta_b2  ) )^2) +2*s_b2  * ( sin( alpha_b2  )*(b2s(1,1)  - x0)+(b2s(1,2)  - y0)* sin( beta_b2  ) ) + a_b2  * ( ((b2s(1,1)  - x0))^2 * ( cos( beta_b2  ) )^2+ ((b2s(1,2)  - y0))^2 * ( cos( alpha_b2  ) )^2+ 2*(b2s(1,1)  - x0)*(b2s(1,2)  - y0)* sin( alpha_b2  )	* sin( beta_b2  )  )  ;
a66_b3 		= I_b3	*( ( sin( alpha_b3  ) )^2 + ( sin( beta_b3  ) )^2) +2*s_b3  * ( sin( alpha_b3  )*(b3s(1,1)  - x0)+(b3s(1,2)  - y0)* sin( beta_b3  ) ) + a_b3  * ( ((b3s(1,1)  - x0))^2 * ( cos( beta_b3  ) )^2+ ((b3s(1,2)  - y0))^2 * ( cos( alpha_b3  ) )^2+ 2*(b3s(1,1)  - x0)*(b3s(1,2)  - y0)* sin( alpha_b3  )	* sin( beta_b3  )  )  ;
a66_b4 		= I_b4	*( ( sin( alpha_b4  ) )^2 + ( sin( beta_b4  ) )^2) +2*s_b4  * ( sin( alpha_b4  )*(b4s(1,1)  - x0)+(b4s(1,2)  - y0)* sin( beta_b4  ) ) + a_b4  * ( ((b4s(1,1)  - x0))^2 * ( cos( beta_b4  ) )^2+ ((b4s(1,2)  - y0))^2 * ( cos( alpha_b4  ) )^2+ 2*(b4s(1,1)  - x0)*(b4s(1,2)  - y0)* sin( alpha_b4  )	* sin( beta_b4  )  )  ;
  
a66(1,1) = a66_c1 		;
a66(1,2) = a66_c2 		;
a66(1,3) = a66_c3 		;
a66(1,4) = a66_c4 		;
a66(1,5) = a66_c5 		;
a66(1,6) = a66_c6 		;
a66(1,7) = a66_c7 		;
a66(1,8) = a66_c8 		;
a66(1,9) = a66_c9 		;
a66(1,10)= a66_c10		;
a66(1,11)= a66_c11		;
a66(1,12)= a66_c12		;
a66(1,13)= a66_c13		;
a66(1,14)= a66_c14		;
a66(1,15)= a66_c15		;
a66(1,16)= a66_c16		;
a66(1,17)= a66_c17		;
a66(1,18)= a66_c18		;
a66(1,19)= a66_c19		;
a66(1,20)= a66_c20		;
a66(1,21)= a66_p1 		;
a66(1,22)= a66_p2 		;
a66(1,23)= a66_p3 		;
a66(1,24)= a66_p4 		;
a66(1,25)= a66_b1 		;
a66(1,26)= a66_b2 		;
a66(1,27)= a66_b3 		;
a66(1,28)= a66_b4 		;

a54_c1 		= -I_c1  * sin( alpha_c1  )	* sin( beta_c1  )-s_c1  * ( sin( alpha_c1  )*(c1h(1,2)  - y0)+(c1h(1,1)  - x0)* sin( beta_c1  ) ) - a_c1  * ( (c1h(1,3)  - z0)* sin( alpha_c1  )-(c1h(1,1)  - x0)*sin( gamma_c1  ) ) * ( (c1h(1,2)  - y0)* sin( gamma_c1  )	-(c1h(1,3)  - z0)* sin( beta_c1  ) ) -a_c1  * (c1h(1,1)  - x0)*(c1h(1,2)  - y0) ;
a54_c2 		= -I_c2  * sin( alpha_c2  )	* sin( beta_c2  )-s_c2  * ( sin( alpha_c2  )*(c2h(1,2)  - y0)+(c2h(1,1)  - x0)* sin( beta_c2  ) ) - a_c2  * ( (c2h(1,3)  - z0)* sin( alpha_c2  )-(c2h(1,1)  - x0)*sin( gamma_c2  ) ) * ( (c2h(1,2)  - y0)* sin( gamma_c2  )	-(c2h(1,3)  - z0)* sin( beta_c2  ) ) -a_c2  * (c2h(1,1)  - x0)*(c2h(1,2)  - y0) ;
a54_c3 		= -I_c3  * sin( alpha_c3  )	* sin( beta_c3  )-s_c3  * ( sin( alpha_c3  )*(c3h(1,2)  - y0)+(c3h(1,1)  - x0)* sin( beta_c3  ) ) - a_c3  * ( (c3h(1,3)  - z0)* sin( alpha_c3  )-(c3h(1,1)  - x0)*sin( gamma_c3  ) ) * ( (c3h(1,2)  - y0)* sin( gamma_c3  )	-(c3h(1,3)  - z0)* sin( beta_c3  ) ) -a_c3  * (c3h(1,1)  - x0)*(c3h(1,2)  - y0) ;
a54_c4 		= -I_c4  * sin( alpha_c4  )	* sin( beta_c4  )-s_c4  * ( sin( alpha_c4  )*(c4h(1,2)  - y0)+(c4h(1,1)  - x0)* sin( beta_c4  ) ) - a_c4  * ( (c4h(1,3)  - z0)* sin( alpha_c4  )-(c4h(1,1)  - x0)*sin( gamma_c4  ) ) * ( (c4h(1,2)  - y0)* sin( gamma_c4  )	-(c4h(1,3)  - z0)* sin( beta_c4  ) ) -a_c4  * (c4h(1,1)  - x0)*(c4h(1,2)  - y0) ;
a54_c5 		= -I_c5  * sin( alpha_c5  )	* sin( beta_c5  )-s_c5  * ( sin( alpha_c5  )*(c5h(1,2)  - y0)+(c5h(1,1)  - x0)* sin( beta_c5  ) ) - a_c5  * ( (c5h(1,3)  - z0)* sin( alpha_c5  )-(c5h(1,1)  - x0)*sin( gamma_c5  ) ) * ( (c5h(1,2)  - y0)* sin( gamma_c5  )	-(c5h(1,3)  - z0)* sin( beta_c5  ) ) -a_c5  * (c5h(1,1)  - x0)*(c5h(1,2)  - y0) ;
a54_c6 		= -I_c6  * sin( alpha_c6  )	* sin( beta_c6  )-s_c6  * ( sin( alpha_c6  )*(c6h(1,2)  - y0)+(c6h(1,1)  - x0)* sin( beta_c6  ) ) - a_c6  * ( (c6h(1,3)  - z0)* sin( alpha_c6  )-(c6h(1,1)  - x0)*sin( gamma_c6  ) ) * ( (c6h(1,2)  - y0)* sin( gamma_c6  )	-(c6h(1,3)  - z0)* sin( beta_c6  ) ) -a_c6  * (c6h(1,1)  - x0)*(c6h(1,2)  - y0) ;
a54_c7 		= -I_c7  * sin( alpha_c7  )	* sin( beta_c7  )-s_c7  * ( sin( alpha_c7  )*(c7h(1,2)  - y0)+(c7h(1,1)  - x0)* sin( beta_c7  ) ) - a_c7  * ( (c7h(1,3)  - z0)* sin( alpha_c7  )-(c7h(1,1)  - x0)*sin( gamma_c7  ) ) * ( (c7h(1,2)  - y0)* sin( gamma_c7  )	-(c7h(1,3)  - z0)* sin( beta_c7  ) ) -a_c7  * (c7h(1,1)  - x0)*(c7h(1,2)  - y0) ;
a54_c8 		= -I_c8  * sin( alpha_c8  )	* sin( beta_c8  )-s_c8  * ( sin( alpha_c8  )*(c8h(1,2)  - y0)+(c8h(1,1)  - x0)* sin( beta_c8  ) ) - a_c8  * ( (c8h(1,3)  - z0)* sin( alpha_c8  )-(c8h(1,1)  - x0)*sin( gamma_c8  ) ) * ( (c8h(1,2)  - y0)* sin( gamma_c8  )	-(c8h(1,3)  - z0)* sin( beta_c8  ) ) -a_c8  * (c8h(1,1)  - x0)*(c8h(1,2)  - y0) ;
a54_c9 		= -I_c9  * sin( alpha_c9  )	* sin( beta_c9  )-s_c9  * ( sin( alpha_c9  )*(c9h(1,2)  - y0)+(c9h(1,1)  - x0)* sin( beta_c9  ) ) - a_c9  * ( (c9h(1,3)  - z0)* sin( alpha_c9  )-(c9h(1,1)  - x0)*sin( gamma_c9  ) ) * ( (c9h(1,2)  - y0)* sin( gamma_c9  )	-(c9h(1,3)  - z0)* sin( beta_c9  ) ) -a_c9  * (c9h(1,1)  - x0)*(c9h(1,2)  - y0) ;
a54_c10		= -I_c10 * sin( alpha_c10 )	* sin( beta_c10 )-s_c10 * ( sin( alpha_c10 )*(c10h(1,2) - y0)+(c10h(1,1) - x0)* sin( beta_c10 ) ) - a_c10 * ( (c10h(1,3) - z0)* sin( alpha_c10 )-(c10h(1,1) - x0)*sin( gamma_c10 ) ) * ( (c10h(1,2) - y0)* sin( gamma_c10 )	-(c10h(1,3) - z0)* sin( beta_c10 ) ) -a_c10 * (c10h(1,1) - x0)*(c10h(1,2) - y0) ;
a54_c11		= -I_c11 * sin( alpha_c11 )	* sin( beta_c11 )-s_c11 * ( sin( alpha_c11 )*(c11h(1,2) - y0)+(c11h(1,1) - x0)* sin( beta_c11 ) ) - a_c11 * ( (c11h(1,3) - z0)* sin( alpha_c11 )-(c11h(1,1) - x0)*sin( gamma_c11 ) ) * ( (c11h(1,2) - y0)* sin( gamma_c11 )	-(c11h(1,3) - z0)* sin( beta_c11 ) ) -a_c11 * (c11h(1,1) - x0)*(c11h(1,2) - y0) ;
a54_c12		= -I_c12 * sin( alpha_c12 )	* sin( beta_c12 )-s_c12 * ( sin( alpha_c12 )*(c12h(1,2) - y0)+(c12h(1,1) - x0)* sin( beta_c12 ) ) - a_c12 * ( (c12h(1,3) - z0)* sin( alpha_c12 )-(c12h(1,1) - x0)*sin( gamma_c12 ) ) * ( (c12h(1,2) - y0)* sin( gamma_c12 )	-(c12h(1,3) - z0)* sin( beta_c12 ) ) -a_c12 * (c12h(1,1) - x0)*(c12h(1,2) - y0) ;
a54_c13		= -I_c13 * sin( alpha_c13 )	* sin( beta_c13 )-s_c13 * ( sin( alpha_c13 )*(c13h(1,2) - y0)+(c13h(1,1) - x0)* sin( beta_c13 ) ) - a_c13 * ( (c13h(1,3) - z0)* sin( alpha_c13 )-(c13h(1,1) - x0)*sin( gamma_c13 ) ) * ( (c13h(1,2) - y0)* sin( gamma_c13 )	-(c13h(1,3) - z0)* sin( beta_c13 ) ) -a_c13 * (c13h(1,1) - x0)*(c13h(1,2) - y0) ;
a54_c14		= -I_c14 * sin( alpha_c14 )	* sin( beta_c14 )-s_c14 * ( sin( alpha_c14 )*(c14h(1,2) - y0)+(c14h(1,1) - x0)* sin( beta_c14 ) ) - a_c14 * ( (c14h(1,3) - z0)* sin( alpha_c14 )-(c14h(1,1) - x0)*sin( gamma_c14 ) ) * ( (c14h(1,2) - y0)* sin( gamma_c14 )	-(c14h(1,3) - z0)* sin( beta_c14 ) ) -a_c14 * (c14h(1,1) - x0)*(c14h(1,2) - y0) ;
a54_c15		= -I_c15 * sin( alpha_c15 )	* sin( beta_c15 )-s_c15 * ( sin( alpha_c15 )*(c15h(1,2) - y0)+(c15h(1,1) - x0)* sin( beta_c15 ) ) - a_c15 * ( (c15h(1,3) - z0)* sin( alpha_c15 )-(c15h(1,1) - x0)*sin( gamma_c15 ) ) * ( (c15h(1,2) - y0)* sin( gamma_c15 )	-(c15h(1,3) - z0)* sin( beta_c15 ) ) -a_c15 * (c15h(1,1) - x0)*(c15h(1,2) - y0) ;
a54_c16		= -I_c16 * sin( alpha_c16 )	* sin( beta_c16 )-s_c16 * ( sin( alpha_c16 )*(c16h(1,2) - y0)+(c16h(1,1) - x0)* sin( beta_c16 ) ) - a_c16 * ( (c16h(1,3) - z0)* sin( alpha_c16 )-(c16h(1,1) - x0)*sin( gamma_c16 ) ) * ( (c16h(1,2) - y0)* sin( gamma_c16 )	-(c16h(1,3) - z0)* sin( beta_c16 ) ) -a_c16 * (c16h(1,1) - x0)*(c16h(1,2) - y0) ;
a54_c17		= -I_c17 * sin( alpha_c17 )	* sin( beta_c17 )-s_c17 * ( sin( alpha_c17 )*(c17h(1,2) - y0)+(c17h(1,1) - x0)* sin( beta_c17 ) ) - a_c17 * ( (c17h(1,3) - z0)* sin( alpha_c17 )-(c17h(1,1) - x0)*sin( gamma_c17 ) ) * ( (c17h(1,2) - y0)* sin( gamma_c17 )	-(c17h(1,3) - z0)* sin( beta_c17 ) ) -a_c17 * (c17h(1,1) - x0)*(c17h(1,2) - y0) ;
a54_c18		= -I_c18 * sin( alpha_c18 )	* sin( beta_c18 )-s_c18 * ( sin( alpha_c18 )*(c18h(1,2) - y0)+(c18h(1,1) - x0)* sin( beta_c18 ) ) - a_c18 * ( (c18h(1,3) - z0)* sin( alpha_c18 )-(c18h(1,1) - x0)*sin( gamma_c18 ) ) * ( (c18h(1,2) - y0)* sin( gamma_c18 )	-(c18h(1,3) - z0)* sin( beta_c18 ) ) -a_c18 * (c18h(1,1) - x0)*(c18h(1,2) - y0) ;
a54_c19		= -I_c19 * sin( alpha_c19 )	* sin( beta_c19 )-s_c19 * ( sin( alpha_c19 )*(c19h(1,2) - y0)+(c19h(1,1) - x0)* sin( beta_c19 ) ) - a_c19 * ( (c19h(1,3) - z0)* sin( alpha_c19 )-(c19h(1,1) - x0)*sin( gamma_c19 ) ) * ( (c19h(1,2) - y0)* sin( gamma_c19 )	-(c19h(1,3) - z0)* sin( beta_c19 ) ) -a_c19 * (c19h(1,1) - x0)*(c19h(1,2) - y0) ;
a54_c20		= -I_c20 * sin( alpha_c20 )	* sin( beta_c20 )-s_c20 * ( sin( alpha_c20 )*(c20h(1,2) - y0)+(c20h(1,1) - x0)* sin( beta_c20 ) ) - a_c20 * ( (c20h(1,3) - z0)* sin( alpha_c20 )-(c20h(1,1) - x0)*sin( gamma_c20 ) ) * ( (c20h(1,2) - y0)* sin( gamma_c20 )	-(c20h(1,3) - z0)* sin( beta_c20 ) ) -a_c20 * (c20h(1,1) - x0)*(c20h(1,2) - y0) ;
a54_p1 		= -I_p1	 * sin( alpha_p1  )	* sin( beta_p1  )-s_p1  * ( sin( alpha_p1  )*(p1a(1,2)  - y0)+(p1a(1,1)  - x0)* sin( beta_p1  ) ) - a_p1  * ( (p1a(1,3)  - z0)* sin( alpha_p1  )-(p1a(1,1)  - x0)*sin( gamma_p1  ) ) * ( (p1a(1,2)  - y0)* sin( gamma_p1  )	-(p1a(1,3)  - z0)* sin( beta_p1  ) ) -a_p1	* (p1a(1,1)  - x0)*(p1a(1,2)  - y0) ;
a54_p2 		= -I_p2	 * sin( alpha_p2  )	* sin( beta_p2  )-s_p2  * ( sin( alpha_p2  )*(p2a(1,2)  - y0)+(p2a(1,1)  - x0)* sin( beta_p2  ) ) - a_p2  * ( (p2a(1,3)  - z0)* sin( alpha_p2  )-(p2a(1,1)  - x0)*sin( gamma_p2  ) ) * ( (p2a(1,2)  - y0)* sin( gamma_p2  )	-(p2a(1,3)  - z0)* sin( beta_p2  ) ) -a_p2	* (p2a(1,1)  - x0)*(p2a(1,2)  - y0) ;
a54_p3 		= -I_p3	 * sin( alpha_p3  )	* sin( beta_p3  )-s_p3  * ( sin( alpha_p3  )*(p3a(1,2)  - y0)+(p3a(1,1)  - x0)* sin( beta_p3  ) ) - a_p3  * ( (p3a(1,3)  - z0)* sin( alpha_p3  )-(p3a(1,1)  - x0)*sin( gamma_p3  ) ) * ( (p3a(1,2)  - y0)* sin( gamma_p3  )	-(p3a(1,3)  - z0)* sin( beta_p3  ) ) -a_p3	* (p3a(1,1)  - x0)*(p3a(1,2)  - y0) ;
a54_p4 		= -I_p4	 * sin( alpha_p4  )	* sin( beta_p4  )-s_p4  * ( sin( alpha_p4  )*(p4a(1,2)  - y0)+(p4a(1,1)  - x0)* sin( beta_p4  ) ) - a_p4  * ( (p4a(1,3)  - z0)* sin( alpha_p4  )-(p4a(1,1)  - x0)*sin( gamma_p4  ) ) * ( (p4a(1,2)  - y0)* sin( gamma_p4  )	-(p4a(1,3)  - z0)* sin( beta_p4  ) ) -a_p4	* (p4a(1,1)  - x0)*(p4a(1,2)  - y0) ;
a54_b1 		= -I_b1	 * sin( alpha_b1  )	* sin( beta_b1  )-s_b1  * ( sin( alpha_b1  )*(b1s(1,2)  - y0)+(b1s(1,1)  - x0)* sin( beta_b1  ) ) - a_b1  * ( (b1s(1,3)  - z0)* sin( alpha_b1  )-(b1s(1,1)  - x0)*sin( gamma_b1  ) ) * ( (b1s(1,2)  - y0)* sin( gamma_b1  )	-(b1s(1,3)  - z0)* sin( beta_b1  ) ) -a_b1	* (b1s(1,1)  - x0)*(b1s(1,2)  - y0) ;
a54_b2 		= -I_b2	 * sin( alpha_b2  )	* sin( beta_b2  )-s_b2  * ( sin( alpha_b2  )*(b2s(1,2)  - y0)+(b2s(1,1)  - x0)* sin( beta_b2  ) ) - a_b2  * ( (b2s(1,3)  - z0)* sin( alpha_b2  )-(b2s(1,1)  - x0)*sin( gamma_b2  ) ) * ( (b2s(1,2)  - y0)* sin( gamma_b2  )	-(b2s(1,3)  - z0)* sin( beta_b2  ) ) -a_b2	* (b2s(1,1)  - x0)*(b2s(1,2)  - y0) ;
a54_b3 		= -I_b3	 * sin( alpha_b3  )	* sin( beta_b3  )-s_b3  * ( sin( alpha_b3  )*(b3s(1,2)  - y0)+(b3s(1,1)  - x0)* sin( beta_b3  ) ) - a_b3  * ( (b3s(1,3)  - z0)* sin( alpha_b3  )-(b3s(1,1)  - x0)*sin( gamma_b3  ) ) * ( (b3s(1,2)  - y0)* sin( gamma_b3  )	-(b3s(1,3)  - z0)* sin( beta_b3  ) ) -a_b3	* (b3s(1,1)  - x0)*(b3s(1,2)  - y0) ;
a54_b4 		= -I_b4	 * sin( alpha_b4  )	* sin( beta_b4  )-s_b4  * ( sin( alpha_b4  )*(b4s(1,2)  - y0)+(b4s(1,1)  - x0)* sin( beta_b4  ) ) - a_b4  * ( (b4s(1,3)  - z0)* sin( alpha_b4  )-(b4s(1,1)  - x0)*sin( gamma_b4  ) ) * ( (b4s(1,2)  - y0)* sin( gamma_b4  )	-(b4s(1,3)  - z0)* sin( beta_b4  ) ) -a_b4	* (b4s(1,1)  - x0)*(b4s(1,2)  - y0) ;
  
a54(1,1) = a54_c1 		;
a54(1,2) = a54_c2 		;
a54(1,3) = a54_c3 		;
a54(1,4) = a54_c4 		;
a54(1,5) = a54_c5 		;
a54(1,6) = a54_c6 		;
a54(1,7) = a54_c7 		;
a54(1,8) = a54_c8 		;
a54(1,9) = a54_c9 		;
a54(1,10)= a54_c10		;
a54(1,11)= a54_c11		;
a54(1,12)= a54_c12		;
a54(1,13)= a54_c13		;
a54(1,14)= a54_c14		;
a54(1,15)= a54_c15		;
a54(1,16)= a54_c16		;
a54(1,17)= a54_c17		;
a54(1,18)= a54_c18		;
a54(1,19)= a54_c19		;
a54(1,20)= a54_c20		;
a54(1,21)= a54_p1 		;
a54(1,22)= a54_p2 		;
a54(1,23)= a54_p3 		;
a54(1,24)= a54_p4 		;
a54(1,25)= a54_b1 		;
a54(1,26)= a54_b2 		;
a54(1,27)= a54_b3 		;
a54(1,28)= a54_b4 		;

a45_c1 		= a54_c1 		;
a45_c2 		= a54_c2 		;
a45_c3 		= a54_c3 		;
a45_c4 		= a54_c4 		;
a45_c5 		= a54_c5 		;
a45_c6 		= a54_c6 		;
a45_c7 		= a54_c7 		;
a45_c8 		= a54_c8 		;
a45_c9 		= a54_c9 		;
a45_c10		= a54_c10		;
a45_c11		= a54_c11		;
a45_c12		= a54_c12		;
a45_c13		= a54_c13		;
a45_c14		= a54_c14		;
a45_c15		= a54_c15		;
a45_c16		= a54_c16		;
a45_c17		= a54_c17		;
a45_c18		= a54_c18		;
a45_c19		= a54_c19		;
a45_c20		= a54_c20		;
a45_p1 		= a54_p1 		;
a45_p2 		= a54_p2 		;
a45_p3 		= a54_p3 		;
a45_p4 		= a54_p4 		;
a45_b1 		= a54_b1 		;
a45_b2 		= a54_b2 		;
a45_b3 		= a54_b3 		;
a45_b4 		= a54_b4 		;
 
a45(1,1) = a45_c1 		;
a45(1,2) = a45_c2 		;
a45(1,3) = a45_c3 		;
a45(1,4) = a45_c4 		;
a45(1,5) = a45_c5 		;
a45(1,6) = a45_c6 		;
a45(1,7) = a45_c7 		;
a45(1,8) = a45_c8 		;
a45(1,9) = a45_c9 		;
a45(1,10)= a45_c10		;
a45(1,11)= a45_c11		;
a45(1,12)= a45_c12		;
a45(1,13)= a45_c13		;
a45(1,14)= a45_c14		;
a45(1,15)= a45_c15		;
a45(1,16)= a45_c16		;
a45(1,17)= a45_c17		;
a45(1,18)= a45_c18		;
a45(1,19)= a45_c19		;
a45(1,20)= a45_c20		;
a45(1,21)= a45_p1 		;
a45(1,22)= a45_p2 		;
a45(1,23)= a45_p3 		;
a45(1,24)= a45_p4 		;
a45(1,25)= a45_b1 		;
a45(1,26)= a45_b2 		;
a45(1,27)= a45_b3 		;
a45(1,28)= a45_b4 		;

a56_c1 		= -I_c1  * sin( gamma_c1  )	* sin( beta_c1  )-s_c1  * ( sin( gamma_c1  )*(c1h(1,2)  - y0)+(c1h(1,3)  - z0)* sin( beta_c1  ) ) - a_c1  * ( (c1h(1,1)  - x0)*sin( beta_c1  )-(c1h(1,2)  - y0)* sin( alpha_c1  ) ) * ( (c1h(1,3)  - z0)* sin( alpha_c1  )-(c1h(1,1)  - x0)* sin( gamma_c1  ) ) -a_c1  * (c1h(1,3)  - z0)*(c1h(1,2)  - y0) ;
a56_c2 		= -I_c2  * sin( gamma_c2  )	* sin( beta_c2  )-s_c2  * ( sin( gamma_c2  )*(c2h(1,2)  - y0)+(c2h(1,3)  - z0)* sin( beta_c2  ) ) - a_c2  * ( (c2h(1,1)  - x0)*sin( beta_c2  )-(c2h(1,2)  - y0)* sin( alpha_c2  ) ) * ( (c2h(1,3)  - z0)* sin( alpha_c2  )-(c2h(1,1)  - x0)* sin( gamma_c2  ) ) -a_c2  * (c2h(1,3)  - z0)*(c2h(1,2)  - y0) ;
a56_c3 		= -I_c3  * sin( gamma_c3  )	* sin( beta_c3  )-s_c3  * ( sin( gamma_c3  )*(c3h(1,2)  - y0)+(c3h(1,3)  - z0)* sin( beta_c3  ) ) - a_c3  * ( (c3h(1,1)  - x0)*sin( beta_c3  )-(c3h(1,2)  - y0)* sin( alpha_c3  ) ) * ( (c3h(1,3)  - z0)* sin( alpha_c3  )-(c3h(1,1)  - x0)* sin( gamma_c3  ) ) -a_c3  * (c3h(1,3)  - z0)*(c3h(1,2)  - y0) ;
a56_c4 		= -I_c4  * sin( gamma_c4  )	* sin( beta_c4  )-s_c4  * ( sin( gamma_c4  )*(c4h(1,2)  - y0)+(c4h(1,3)  - z0)* sin( beta_c4  ) ) - a_c4  * ( (c4h(1,1)  - x0)*sin( beta_c4  )-(c4h(1,2)  - y0)* sin( alpha_c4  ) ) * ( (c4h(1,3)  - z0)* sin( alpha_c4  )-(c4h(1,1)  - x0)* sin( gamma_c4  ) ) -a_c4  * (c4h(1,3)  - z0)*(c4h(1,2)  - y0) ;
a56_c5 		= -I_c5  * sin( gamma_c5  )	* sin( beta_c5  )-s_c5  * ( sin( gamma_c5  )*(c5h(1,2)  - y0)+(c5h(1,3)  - z0)* sin( beta_c5  ) ) - a_c5  * ( (c5h(1,1)  - x0)*sin( beta_c5  )-(c5h(1,2)  - y0)* sin( alpha_c5  ) ) * ( (c5h(1,3)  - z0)* sin( alpha_c5  )-(c5h(1,1)  - x0)* sin( gamma_c5  ) ) -a_c5  * (c5h(1,3)  - z0)*(c5h(1,2)  - y0) ;
a56_c6 		= -I_c6  * sin( gamma_c6  )	* sin( beta_c6  )-s_c6  * ( sin( gamma_c6  )*(c6h(1,2)  - y0)+(c6h(1,3)  - z0)* sin( beta_c6  ) ) - a_c6  * ( (c6h(1,1)  - x0)*sin( beta_c6  )-(c6h(1,2)  - y0)* sin( alpha_c6  ) ) * ( (c6h(1,3)  - z0)* sin( alpha_c6  )-(c6h(1,1)  - x0)* sin( gamma_c6  ) ) -a_c6  * (c6h(1,3)  - z0)*(c6h(1,2)  - y0) ;
a56_c7 		= -I_c7  * sin( gamma_c7  )	* sin( beta_c7  )-s_c7  * ( sin( gamma_c7  )*(c7h(1,2)  - y0)+(c7h(1,3)  - z0)* sin( beta_c7  ) ) - a_c7  * ( (c7h(1,1)  - x0)*sin( beta_c7  )-(c7h(1,2)  - y0)* sin( alpha_c7  ) ) * ( (c7h(1,3)  - z0)* sin( alpha_c7  )-(c7h(1,1)  - x0)* sin( gamma_c7  ) ) -a_c7  * (c7h(1,3)  - z0)*(c7h(1,2)  - y0) ;
a56_c8 		= -I_c8  * sin( gamma_c8  )	* sin( beta_c8  )-s_c8  * ( sin( gamma_c8  )*(c8h(1,2)  - y0)+(c8h(1,3)  - z0)* sin( beta_c8  ) ) - a_c8  * ( (c8h(1,1)  - x0)*sin( beta_c8  )-(c8h(1,2)  - y0)* sin( alpha_c8  ) ) * ( (c8h(1,3)  - z0)* sin( alpha_c8  )-(c8h(1,1)  - x0)* sin( gamma_c8  ) ) -a_c8  * (c8h(1,3)  - z0)*(c8h(1,2)  - y0) ;
a56_c9 		= -I_c9  * sin( gamma_c9  )	* sin( beta_c9  )-s_c9  * ( sin( gamma_c9  )*(c9h(1,2)  - y0)+(c9h(1,3)  - z0)* sin( beta_c9  ) ) - a_c9  * ( (c9h(1,1)  - x0)*sin( beta_c9  )-(c9h(1,2)  - y0)* sin( alpha_c9  ) ) * ( (c9h(1,3)  - z0)* sin( alpha_c9  )-(c9h(1,1)  - x0)* sin( gamma_c9  ) ) -a_c9  * (c9h(1,3)  - z0)*(c9h(1,2)  - y0) ;
a56_c10		= -I_c10 * sin( gamma_c10 )	* sin( beta_c10 )-s_c10 * ( sin( gamma_c10 )*(c10h(1,2) - y0)+(c10h(1,3) - z0)* sin( beta_c10 ) ) - a_c10 * ( (c10h(1,1) - x0)*sin( beta_c10 )-(c10h(1,2) - y0)* sin( alpha_c10 ) ) * ( (c10h(1,3) - z0)* sin( alpha_c10 )-(c10h(1,1) - x0)* sin( gamma_c10 ) ) -a_c10 * (c10h(1,3) - z0)*(c10h(1,2) - y0) ;
a56_c11		= -I_c11 * sin( gamma_c11 )	* sin( beta_c11 )-s_c11 * ( sin( gamma_c11 )*(c11h(1,2) - y0)+(c11h(1,3) - z0)* sin( beta_c11 ) ) - a_c11 * ( (c11h(1,1) - x0)*sin( beta_c11 )-(c11h(1,2) - y0)* sin( alpha_c11 ) ) * ( (c11h(1,3) - z0)* sin( alpha_c11 )-(c11h(1,1) - x0)* sin( gamma_c11 ) ) -a_c11 * (c11h(1,3) - z0)*(c11h(1,2) - y0) ;
a56_c12		= -I_c12 * sin( gamma_c12 )	* sin( beta_c12 )-s_c12 * ( sin( gamma_c12 )*(c12h(1,2) - y0)+(c12h(1,3) - z0)* sin( beta_c12 ) ) - a_c12 * ( (c12h(1,1) - x0)*sin( beta_c12 )-(c12h(1,2) - y0)* sin( alpha_c12 ) ) * ( (c12h(1,3) - z0)* sin( alpha_c12 )-(c12h(1,1) - x0)* sin( gamma_c12 ) ) -a_c12 * (c12h(1,3) - z0)*(c12h(1,2) - y0) ;
a56_c13		= -I_c13 * sin( gamma_c13 )	* sin( beta_c13 )-s_c13 * ( sin( gamma_c13 )*(c13h(1,2) - y0)+(c13h(1,3) - z0)* sin( beta_c13 ) ) - a_c13 * ( (c13h(1,1) - x0)*sin( beta_c13 )-(c13h(1,2) - y0)* sin( alpha_c13 ) ) * ( (c13h(1,3) - z0)* sin( alpha_c13 )-(c13h(1,1) - x0)* sin( gamma_c13 ) ) -a_c13 * (c13h(1,3) - z0)*(c13h(1,2) - y0) ;
a56_c14		= -I_c14 * sin( gamma_c14 )	* sin( beta_c14 )-s_c14 * ( sin( gamma_c14 )*(c14h(1,2) - y0)+(c14h(1,3) - z0)* sin( beta_c14 ) ) - a_c14 * ( (c14h(1,1) - x0)*sin( beta_c14 )-(c14h(1,2) - y0)* sin( alpha_c14 ) ) * ( (c14h(1,3) - z0)* sin( alpha_c14 )-(c14h(1,1) - x0)* sin( gamma_c14 ) ) -a_c14 * (c14h(1,3) - z0)*(c14h(1,2) - y0) ;
a56_c15		= -I_c15 * sin( gamma_c15 )	* sin( beta_c15 )-s_c15 * ( sin( gamma_c15 )*(c15h(1,2) - y0)+(c15h(1,3) - z0)* sin( beta_c15 ) ) - a_c15 * ( (c15h(1,1) - x0)*sin( beta_c15 )-(c15h(1,2) - y0)* sin( alpha_c15 ) ) * ( (c15h(1,3) - z0)* sin( alpha_c15 )-(c15h(1,1) - x0)* sin( gamma_c15 ) ) -a_c15 * (c15h(1,3) - z0)*(c15h(1,2) - y0) ;
a56_c16		= -I_c16 * sin( gamma_c16 )	* sin( beta_c16 )-s_c16 * ( sin( gamma_c16 )*(c16h(1,2) - y0)+(c16h(1,3) - z0)* sin( beta_c16 ) ) - a_c16 * ( (c16h(1,1) - x0)*sin( beta_c16 )-(c16h(1,2) - y0)* sin( alpha_c16 ) ) * ( (c16h(1,3) - z0)* sin( alpha_c16 )-(c16h(1,1) - x0)* sin( gamma_c16 ) ) -a_c16 * (c16h(1,3) - z0)*(c16h(1,2) - y0) ;
a56_c17		= -I_c17 * sin( gamma_c17 )	* sin( beta_c17 )-s_c17 * ( sin( gamma_c17 )*(c17h(1,2) - y0)+(c17h(1,3) - z0)* sin( beta_c17 ) ) - a_c17 * ( (c17h(1,1) - x0)*sin( beta_c17 )-(c17h(1,2) - y0)* sin( alpha_c17 ) ) * ( (c17h(1,3) - z0)* sin( alpha_c17 )-(c17h(1,1) - x0)* sin( gamma_c17 ) ) -a_c17 * (c17h(1,3) - z0)*(c17h(1,2) - y0) ;
a56_c18		= -I_c18 * sin( gamma_c18 )	* sin( beta_c18 )-s_c18 * ( sin( gamma_c18 )*(c18h(1,2) - y0)+(c18h(1,3) - z0)* sin( beta_c18 ) ) - a_c18 * ( (c18h(1,1) - x0)*sin( beta_c18 )-(c18h(1,2) - y0)* sin( alpha_c18 ) ) * ( (c18h(1,3) - z0)* sin( alpha_c18 )-(c18h(1,1) - x0)* sin( gamma_c18 ) ) -a_c18 * (c18h(1,3) - z0)*(c18h(1,2) - y0) ;
a56_c19		= -I_c19 * sin( gamma_c19 )	* sin( beta_c19 )-s_c19 * ( sin( gamma_c19 )*(c19h(1,2) - y0)+(c19h(1,3) - z0)* sin( beta_c19 ) ) - a_c19 * ( (c19h(1,1) - x0)*sin( beta_c19 )-(c19h(1,2) - y0)* sin( alpha_c19 ) ) * ( (c19h(1,3) - z0)* sin( alpha_c19 )-(c19h(1,1) - x0)* sin( gamma_c19 ) ) -a_c19 * (c19h(1,3) - z0)*(c19h(1,2) - y0) ;
a56_c20		= -I_c20 * sin( gamma_c20 )	* sin( beta_c20 )-s_c20 * ( sin( gamma_c20 )*(c20h(1,2) - y0)+(c20h(1,3) - z0)* sin( beta_c20 ) ) - a_c20 * ( (c20h(1,1) - x0)*sin( beta_c20 )-(c20h(1,2) - y0)* sin( alpha_c20 ) ) * ( (c20h(1,3) - z0)* sin( alpha_c20 )-(c20h(1,1) - x0)* sin( gamma_c20 ) ) -a_c20 * (c20h(1,3) - z0)*(c20h(1,2) - y0) ;
a56_p1 		= -I_p1	 * sin( gamma_p1  )	* sin( beta_p1  )-s_p1  * ( sin( gamma_p1  )*(p1a(1,2)  - y0)+(p1a(1,3)  - z0)* sin( beta_p1  ) ) - a_p1  * ( (p1a(1,1)  - x0)*sin( beta_p1  )-(p1a(1,2)  - y0)* sin( alpha_p1  ) ) * ( (p1a(1,3)  - z0)* sin( alpha_p1  )-(p1a(1,1)  - x0)* sin( gamma_p1  ) ) -a_p1  * (p1a(1,3)  - z0)*(p1a(1,2)  - y0) ;
a56_p2 		= -I_p2	 * sin( gamma_p2  )	* sin( beta_p2  )-s_p2  * ( sin( gamma_p2  )*(p2a(1,2)  - y0)+(p2a(1,3)  - z0)* sin( beta_p2  ) ) - a_p2  * ( (p2a(1,1)  - x0)*sin( beta_p2  )-(p2a(1,2)  - y0)* sin( alpha_p2  ) ) * ( (p2a(1,3)  - z0)* sin( alpha_p2  )-(p2a(1,1)  - x0)* sin( gamma_p2  ) ) -a_p2  * (p2a(1,3)  - z0)*(p2a(1,2)  - y0) ;
a56_p3 		= -I_p3	 * sin( gamma_p3  )	* sin( beta_p3  )-s_p3  * ( sin( gamma_p3  )*(p3a(1,2)  - y0)+(p3a(1,3)  - z0)* sin( beta_p3  ) ) - a_p3  * ( (p3a(1,1)  - x0)*sin( beta_p3  )-(p3a(1,2)  - y0)* sin( alpha_p3  ) ) * ( (p3a(1,3)  - z0)* sin( alpha_p3  )-(p3a(1,1)  - x0)* sin( gamma_p3  ) ) -a_p3  * (p3a(1,3)  - z0)*(p3a(1,2)  - y0) ;
a56_p4 		= -I_p4	 * sin( gamma_p4  )	* sin( beta_p4  )-s_p4  * ( sin( gamma_p4  )*(p4a(1,2)  - y0)+(p4a(1,3)  - z0)* sin( beta_p4  ) ) - a_p4  * ( (p4a(1,1)  - x0)*sin( beta_p4  )-(p4a(1,2)  - y0)* sin( alpha_p4  ) ) * ( (p4a(1,3)  - z0)* sin( alpha_p4  )-(p4a(1,1)  - x0)* sin( gamma_p4  ) ) -a_p4  * (p4a(1,3)  - z0)*(p4a(1,2)  - y0) ;
a56_b1 		= -I_b1	 * sin( gamma_b1  )	* sin( beta_b1  )-s_b1  * ( sin( gamma_b1  )*(b1s(1,2)  - y0)+(b1s(1,3)  - z0)* sin( beta_b1  ) ) - a_b1  * ( (b1s(1,1)  - x0)*sin( beta_b1  )-(b1s(1,2)  - y0)* sin( alpha_b1  ) ) * ( (b1s(1,3)  - z0)* sin( alpha_b1  )-(b1s(1,1)  - x0)* sin( gamma_b1  ) ) -a_b1  * (b1s(1,3)  - z0)*(b1s(1,2)  - y0) ;
a56_b2 		= -I_b2	 * sin( gamma_b2  )	* sin( beta_b2  )-s_b2  * ( sin( gamma_b2  )*(b2s(1,2)  - y0)+(b2s(1,3)  - z0)* sin( beta_b2  ) ) - a_b2  * ( (b2s(1,1)  - x0)*sin( beta_b2  )-(b2s(1,2)  - y0)* sin( alpha_b2  ) ) * ( (b2s(1,3)  - z0)* sin( alpha_b2  )-(b2s(1,1)  - x0)* sin( gamma_b2  ) ) -a_b2  * (b2s(1,3)  - z0)*(b2s(1,2)  - y0) ;
a56_b3 		= -I_b3	 * sin( gamma_b3  )	* sin( beta_b3  )-s_b3  * ( sin( gamma_b3  )*(b3s(1,2)  - y0)+(b3s(1,3)  - z0)* sin( beta_b3  ) ) - a_b3  * ( (b3s(1,1)  - x0)*sin( beta_b3  )-(b3s(1,2)  - y0)* sin( alpha_b3  ) ) * ( (b3s(1,3)  - z0)* sin( alpha_b3  )-(b3s(1,1)  - x0)* sin( gamma_b3  ) ) -a_b3  * (b3s(1,3)  - z0)*(b3s(1,2)  - y0) ;
a56_b4 		= -I_b4	 * sin( gamma_b4  )	* sin( beta_b4  )-s_b4  * ( sin( gamma_b4  )*(b4s(1,2)  - y0)+(b4s(1,3)  - z0)* sin( beta_b4  ) ) - a_b4  * ( (b4s(1,1)  - x0)*sin( beta_b4  )-(b4s(1,2)  - y0)* sin( alpha_b4  ) ) * ( (b4s(1,3)  - z0)* sin( alpha_b4  )-(b4s(1,1)  - x0)* sin( gamma_b4  ) ) -a_b4  * (b4s(1,3)  - z0)*(b4s(1,2)  - y0) ;
  
a56(1,1) = a56_c1 		;
a56(1,2) = a56_c2 		;
a56(1,3) = a56_c3 		;
a56(1,4) = a56_c4 		;
a56(1,5) = a56_c5 		;
a56(1,6) = a56_c6 		;
a56(1,7) = a56_c7 		;
a56(1,8) = a56_c8 		;
a56(1,9) = a56_c9 		;
a56(1,10)= a56_c10		;
a56(1,11)= a56_c11		;
a56(1,12)= a56_c12		;
a56(1,13)= a56_c13		;
a56(1,14)= a56_c14		;
a56(1,15)= a56_c15		;
a56(1,16)= a56_c16		;
a56(1,17)= a56_c17		;
a56(1,18)= a56_c18		;
a56(1,19)= a56_c19		;
a56(1,20)= a56_c20		;
a56(1,21)= a56_p1 		;
a56(1,22)= a56_p2 		;
a56(1,23)= a56_p3 		;
a56(1,24)= a56_p4 		;
a56(1,25)= a56_b1 		;
a56(1,26)= a56_b2 		;
a56(1,27)= a56_b3 		;
a56(1,28)= a56_b4 		;

a65_c1 		= a56_c1 		;
a65_c2 		= a56_c2 		;
a65_c3 		= a56_c3 		;
a65_c4 		= a56_c4 		;
a65_c5 		= a56_c5 		;
a65_c6 		= a56_c6 		;
a65_c7 		= a56_c7 		;
a65_c8 		= a56_c8 		;
a65_c9 		= a56_c9 		;
a65_c10		= a56_c10		;
a65_c11		= a56_c11		;
a65_c12		= a56_c12		;
a65_c13		= a56_c13		;
a65_c14		= a56_c14		;
a65_c15		= a56_c15		;
a65_c16		= a56_c16		;
a65_c17		= a56_c17		;
a65_c18		= a56_c18		;
a65_c19		= a56_c19		;
a65_c20		= a56_c20		;
a65_p1 		= a56_p1 		;
a65_p2 		= a56_p2 		;
a65_p3 		= a56_p3 		;
a65_p4 		= a56_p4 		;
a65_b1 		= a56_b1 		;
a65_b2 		= a56_b2 		;
a65_b3 		= a56_b3 		;
a65_b4 		= a56_b4 		;
 
a65(1,1) = a65_c1 		;
a65(1,2) = a65_c2 		;
a65(1,3) = a65_c3 		;
a65(1,4) = a65_c4 		;
a65(1,5) = a65_c5 		;
a65(1,6) = a65_c6 		;
a65(1,7) = a65_c7 		;
a65(1,8) = a65_c8 		;
a65(1,9) = a65_c9 		;
a65(1,10)= a65_c10		;
a65(1,11)= a65_c11		;
a65(1,12)= a65_c12		;
a65(1,13)= a65_c13		;
a65(1,14)= a65_c14		;
a65(1,15)= a65_c15		;
a65(1,16)= a65_c16		;
a65(1,17)= a65_c17		;
a65(1,18)= a65_c18		;
a65(1,19)= a65_c19		;
a65(1,20)= a65_c20		;
a65(1,21)= a65_p1 		;
a65(1,22)= a65_p2 		;
a65(1,23)= a65_p3 		;
a65(1,24)= a65_p4 		;
a65(1,25)= a65_b1 		;
a65(1,26)= a65_b2 		;
a65(1,27)= a65_b3 		;
a65(1,28)= a65_b4 		;

a46_c1 		= -I_c1  * sin( gamma_c1  )	* sin( alpha_c1  )-s_c1  * ( (c1h(1,3)  - z0)* sin( alpha_c1  )+(c1h(1,1)  - x0) ) - a_c1  * ( (c1h(1,2)  - y0)* sin( gamma_c1  )-(c1h(1,3)  - z0)*sin( beta_c1  ) ) * ( (c1h(1,1)  - x0)*sin( beta_c1  )-(c1h(1,2)  - y0)* sin( alpha_c1  )) -a_c1  * (c1h(1,3)  - z0)*(c1h(1,1)  - x0) ;
a46_c2 		= -I_c2  * sin( gamma_c2  )	* sin( alpha_c2  )-s_c2  * ( (c2h(1,3)  - z0)* sin( alpha_c2  )+(c2h(1,1)  - x0) ) - a_c2  * ( (c2h(1,2)  - y0)* sin( gamma_c2  )-(c2h(1,3)  - z0)*sin( beta_c2  ) ) * ( (c2h(1,1)  - x0)*sin( beta_c2  )-(c2h(1,2)  - y0)* sin( alpha_c2  )) -a_c2  * (c2h(1,3)  - z0)*(c2h(1,1)  - x0) ;
a46_c3 		= -I_c3  * sin( gamma_c3  )	* sin( alpha_c3  )-s_c3  * ( (c3h(1,3)  - z0)* sin( alpha_c3  )+(c3h(1,1)  - x0) ) - a_c3  * ( (c3h(1,2)  - y0)* sin( gamma_c3  )-(c3h(1,3)  - z0)*sin( beta_c3  ) ) * ( (c3h(1,1)  - x0)*sin( beta_c3  )-(c3h(1,2)  - y0)* sin( alpha_c3  )) -a_c3  * (c3h(1,3)  - z0)*(c3h(1,1)  - x0) ;
a46_c4 		= -I_c4  * sin( gamma_c4  )	* sin( alpha_c4  )-s_c4  * ( (c4h(1,3)  - z0)* sin( alpha_c4  )+(c4h(1,1)  - x0) ) - a_c4  * ( (c4h(1,2)  - y0)* sin( gamma_c4  )-(c4h(1,3)  - z0)*sin( beta_c4  ) ) * ( (c4h(1,1)  - x0)*sin( beta_c4  )-(c4h(1,2)  - y0)* sin( alpha_c4  )) -a_c4  * (c4h(1,3)  - z0)*(c4h(1,1)  - x0) ;
a46_c5 		= -I_c5  * sin( gamma_c5  )	* sin( alpha_c5  )-s_c5  * ( (c5h(1,3)  - z0)* sin( alpha_c5  )+(c5h(1,1)  - x0) ) - a_c5  * ( (c5h(1,2)  - y0)* sin( gamma_c5  )-(c5h(1,3)  - z0)*sin( beta_c5  ) ) * ( (c5h(1,1)  - x0)*sin( beta_c5  )-(c5h(1,2)  - y0)* sin( alpha_c5  )) -a_c5  * (c5h(1,3)  - z0)*(c5h(1,1)  - x0) ;
a46_c6 		= -I_c6  * sin( gamma_c6  )	* sin( alpha_c6  )-s_c6  * ( (c6h(1,3)  - z0)* sin( alpha_c6  )+(c6h(1,1)  - x0) ) - a_c6  * ( (c6h(1,2)  - y0)* sin( gamma_c6  )-(c6h(1,3)  - z0)*sin( beta_c6  ) ) * ( (c6h(1,1)  - x0)*sin( beta_c6  )-(c6h(1,2)  - y0)* sin( alpha_c6  )) -a_c6  * (c6h(1,3)  - z0)*(c6h(1,1)  - x0) ;
a46_c7 		= -I_c7  * sin( gamma_c7  )	* sin( alpha_c7  )-s_c7  * ( (c7h(1,3)  - z0)* sin( alpha_c7  )+(c7h(1,1)  - x0) ) - a_c7  * ( (c7h(1,2)  - y0)* sin( gamma_c7  )-(c7h(1,3)  - z0)*sin( beta_c7  ) ) * ( (c7h(1,1)  - x0)*sin( beta_c7  )-(c7h(1,2)  - y0)* sin( alpha_c7  )) -a_c7  * (c7h(1,3)  - z0)*(c7h(1,1)  - x0) ;
a46_c8 		= -I_c8  * sin( gamma_c8  )	* sin( alpha_c8  )-s_c8  * ( (c8h(1,3)  - z0)* sin( alpha_c8  )+(c8h(1,1)  - x0) ) - a_c8  * ( (c8h(1,2)  - y0)* sin( gamma_c8  )-(c8h(1,3)  - z0)*sin( beta_c8  ) ) * ( (c8h(1,1)  - x0)*sin( beta_c8  )-(c8h(1,2)  - y0)* sin( alpha_c8  )) -a_c8  * (c8h(1,3)  - z0)*(c8h(1,1)  - x0) ;
a46_c9 		= -I_c9  * sin( gamma_c9  )	* sin( alpha_c9  )-s_c9  * ( (c9h(1,3)  - z0)* sin( alpha_c9  )+(c9h(1,1)  - x0) ) - a_c9  * ( (c9h(1,2)  - y0)* sin( gamma_c9  )-(c9h(1,3)  - z0)*sin( beta_c9  ) ) * ( (c9h(1,1)  - x0)*sin( beta_c9  )-(c9h(1,2)  - y0)* sin( alpha_c9  )) -a_c9  * (c9h(1,3)  - z0)*(c9h(1,1)  - x0) ;
a46_c10		= -I_c10 * sin( gamma_c10 )	* sin( alpha_c10 )-s_c10 * ( (c10h(1,3) - z0)* sin( alpha_c10 )+(c10h(1,1) - x0) ) - a_c10 * ( (c10h(1,2) - y0)* sin( gamma_c10 )-(c10h(1,3) - z0)*sin( beta_c10 ) ) * ( (c10h(1,1) - x0)*sin( beta_c10 )-(c10h(1,2) - y0)* sin( alpha_c10 )) -a_c10 * (c10h(1,3) - z0)*(c10h(1,1) - x0) ;
a46_c11		= -I_c11 * sin( gamma_c11 )	* sin( alpha_c11 )-s_c11 * ( (c11h(1,3) - z0)* sin( alpha_c11 )+(c11h(1,1) - x0) ) - a_c11 * ( (c11h(1,2) - y0)* sin( gamma_c11 )-(c11h(1,3) - z0)*sin( beta_c11 ) ) * ( (c11h(1,1) - x0)*sin( beta_c11 )-(c11h(1,2) - y0)* sin( alpha_c11 )) -a_c11 * (c11h(1,3) - z0)*(c11h(1,1) - x0) ;
a46_c12		= -I_c12 * sin( gamma_c12 )	* sin( alpha_c12 )-s_c12 * ( (c12h(1,3) - z0)* sin( alpha_c12 )+(c12h(1,1) - x0) ) - a_c12 * ( (c12h(1,2) - y0)* sin( gamma_c12 )-(c12h(1,3) - z0)*sin( beta_c12 ) ) * ( (c12h(1,1) - x0)*sin( beta_c12 )-(c12h(1,2) - y0)* sin( alpha_c12 )) -a_c12 * (c12h(1,3) - z0)*(c12h(1,1) - x0) ;
a46_c13		= -I_c13 * sin( gamma_c13 )	* sin( alpha_c13 )-s_c13 * ( (c13h(1,3) - z0)* sin( alpha_c13 )+(c13h(1,1) - x0) ) - a_c13 * ( (c13h(1,2) - y0)* sin( gamma_c13 )-(c13h(1,3) - z0)*sin( beta_c13 ) ) * ( (c13h(1,1) - x0)*sin( beta_c13 )-(c13h(1,2) - y0)* sin( alpha_c13 )) -a_c13 * (c13h(1,3) - z0)*(c13h(1,1) - x0) ;
a46_c14		= -I_c14 * sin( gamma_c14 )	* sin( alpha_c14 )-s_c14 * ( (c14h(1,3) - z0)* sin( alpha_c14 )+(c14h(1,1) - x0) ) - a_c14 * ( (c14h(1,2) - y0)* sin( gamma_c14 )-(c14h(1,3) - z0)*sin( beta_c14 ) ) * ( (c14h(1,1) - x0)*sin( beta_c14 )-(c14h(1,2) - y0)* sin( alpha_c14 )) -a_c14 * (c14h(1,3) - z0)*(c14h(1,1) - x0) ;
a46_c15		= -I_c15 * sin( gamma_c15 )	* sin( alpha_c15 )-s_c15 * ( (c15h(1,3) - z0)* sin( alpha_c15 )+(c15h(1,1) - x0) ) - a_c15 * ( (c15h(1,2) - y0)* sin( gamma_c15 )-(c15h(1,3) - z0)*sin( beta_c15 ) ) * ( (c15h(1,1) - x0)*sin( beta_c15 )-(c15h(1,2) - y0)* sin( alpha_c15 )) -a_c15 * (c15h(1,3) - z0)*(c15h(1,1) - x0) ;
a46_c16		= -I_c16 * sin( gamma_c16 )	* sin( alpha_c16 )-s_c16 * ( (c16h(1,3) - z0)* sin( alpha_c16 )+(c16h(1,1) - x0) ) - a_c16 * ( (c16h(1,2) - y0)* sin( gamma_c16 )-(c16h(1,3) - z0)*sin( beta_c16 ) ) * ( (c16h(1,1) - x0)*sin( beta_c16 )-(c16h(1,2) - y0)* sin( alpha_c16 )) -a_c16 * (c16h(1,3) - z0)*(c16h(1,1) - x0) ;
a46_c17		= -I_c17 * sin( gamma_c17 )	* sin( alpha_c17 )-s_c17 * ( (c17h(1,3) - z0)* sin( alpha_c17 )+(c17h(1,1) - x0) ) - a_c17 * ( (c17h(1,2) - y0)* sin( gamma_c17 )-(c17h(1,3) - z0)*sin( beta_c17 ) ) * ( (c17h(1,1) - x0)*sin( beta_c17 )-(c17h(1,2) - y0)* sin( alpha_c17 )) -a_c17 * (c17h(1,3) - z0)*(c17h(1,1) - x0) ;
a46_c18		= -I_c18 * sin( gamma_c18 )	* sin( alpha_c18 )-s_c18 * ( (c18h(1,3) - z0)* sin( alpha_c18 )+(c18h(1,1) - x0) ) - a_c18 * ( (c18h(1,2) - y0)* sin( gamma_c18 )-(c18h(1,3) - z0)*sin( beta_c18 ) ) * ( (c18h(1,1) - x0)*sin( beta_c18 )-(c18h(1,2) - y0)* sin( alpha_c18 )) -a_c18 * (c18h(1,3) - z0)*(c18h(1,1) - x0) ;
a46_c19		= -I_c19 * sin( gamma_c19 )	* sin( alpha_c19 )-s_c19 * ( (c19h(1,3) - z0)* sin( alpha_c19 )+(c19h(1,1) - x0) ) - a_c19 * ( (c19h(1,2) - y0)* sin( gamma_c19 )-(c19h(1,3) - z0)*sin( beta_c19 ) ) * ( (c19h(1,1) - x0)*sin( beta_c19 )-(c19h(1,2) - y0)* sin( alpha_c19 )) -a_c19 * (c19h(1,3) - z0)*(c19h(1,1) - x0) ;
a46_c20		= -I_c20 * sin( gamma_c20 )	* sin( alpha_c20 )-s_c20 * ( (c20h(1,3) - z0)* sin( alpha_c20 )+(c20h(1,1) - x0) ) - a_c20 * ( (c20h(1,2) - y0)* sin( gamma_c20 )-(c20h(1,3) - z0)*sin( beta_c20 ) ) * ( (c20h(1,1) - x0)*sin( beta_c20 )-(c20h(1,2) - y0)* sin( alpha_c20 )) -a_c20 * (c20h(1,3) - z0)*(c20h(1,1) - x0) ;
a46_p1 		= -I_p1	 * sin( gamma_p1  )	* sin( alpha_p1  )-s_p1  * ( (p1a(1,3)  - z0)* sin( alpha_p1  )+(p1a(1,1)  - x0) ) - a_p1  * ( (p1a(1,2)  - y0)* sin( gamma_p1  )-(p1a(1,3)  - z0)*sin( beta_p1  ) ) * ( (p1a(1,1)  - x0)*sin( beta_p1  )-(p1a(1,2)  - y0)* sin( alpha_p1  )) -a_p1  * (p1a(1,3)  - z0)*(p1a(1,1)  - x0) ;
a46_p2 		= -I_p2	 * sin( gamma_p2  )	* sin( alpha_p2  )-s_p2  * ( (p2a(1,3)  - z0)* sin( alpha_p2  )+(p2a(1,1)  - x0) ) - a_p2  * ( (p2a(1,2)  - y0)* sin( gamma_p2  )-(p2a(1,3)  - z0)*sin( beta_p2  ) ) * ( (p2a(1,1)  - x0)*sin( beta_p2  )-(p2a(1,2)  - y0)* sin( alpha_p2  )) -a_p2  * (p2a(1,3)  - z0)*(p2a(1,1)  - x0) ;
a46_p3 		= -I_p3	 * sin( gamma_p3  )	* sin( alpha_p3  )-s_p3  * ( (p3a(1,3)  - z0)* sin( alpha_p3  )+(p3a(1,1)  - x0) ) - a_p3  * ( (p3a(1,2)  - y0)* sin( gamma_p3  )-(p3a(1,3)  - z0)*sin( beta_p3  ) ) * ( (p3a(1,1)  - x0)*sin( beta_p3  )-(p3a(1,2)  - y0)* sin( alpha_p3  )) -a_p3  * (p3a(1,3)  - z0)*(p3a(1,1)  - x0) ;
a46_p4 		= -I_p4	 * sin( gamma_p4  )	* sin( alpha_p4  )-s_p4  * ( (p4a(1,3)  - z0)* sin( alpha_p4  )+(p4a(1,1)  - x0) ) - a_p4  * ( (p4a(1,2)  - y0)* sin( gamma_p4  )-(p4a(1,3)  - z0)*sin( beta_p4  ) ) * ( (p4a(1,1)  - x0)*sin( beta_p4  )-(p4a(1,2)  - y0)* sin( alpha_p4  )) -a_p4  * (p4a(1,3)  - z0)*(p4a(1,1)  - x0) ;
a46_b1 		= -I_b1	 * sin( gamma_b1  )	* sin( alpha_b1  )-s_b1  * ( (b1s(1,3)  - z0)* sin( alpha_b1  )+(b1s(1,1)  - x0) ) - a_b1  * ( (b1s(1,2)  - y0)* sin( gamma_b1  )-(b1s(1,3)  - z0)*sin( beta_b1  ) ) * ( (b1s(1,1)  - x0)*sin( beta_b1  )-(b1s(1,2)  - y0)* sin( alpha_b1  )) -a_b1  * (b1s(1,3)  - z0)*(b1s(1,1)  - x0) ;
a46_b2 		= -I_b2	 * sin( gamma_b2  )	* sin( alpha_b2  )-s_b2  * ( (b2s(1,3)  - z0)* sin( alpha_b2  )+(b2s(1,1)  - x0) ) - a_b2  * ( (b2s(1,2)  - y0)* sin( gamma_b2  )-(b2s(1,3)  - z0)*sin( beta_b2  ) ) * ( (b2s(1,1)  - x0)*sin( beta_b2  )-(b2s(1,2)  - y0)* sin( alpha_b2  )) -a_b2  * (b2s(1,3)  - z0)*(b2s(1,1)  - x0) ;
a46_b3 		= -I_b3	 * sin( gamma_b3  )	* sin( alpha_b3  )-s_b3  * ( (b3s(1,3)  - z0)* sin( alpha_b3  )+(b3s(1,1)  - x0) ) - a_b3  * ( (b3s(1,2)  - y0)* sin( gamma_b3  )-(b3s(1,3)  - z0)*sin( beta_b3  ) ) * ( (b3s(1,1)  - x0)*sin( beta_b3  )-(b3s(1,2)  - y0)* sin( alpha_b3  )) -a_b3  * (b3s(1,3)  - z0)*(b3s(1,1)  - x0) ;
a46_b4 		= -I_b4	 * sin( gamma_b4  )	* sin( alpha_b4  )-s_b4  * ( (b4s(1,3)  - z0)* sin( alpha_b4  )+(b4s(1,1)  - x0) ) - a_b4  * ( (b4s(1,2)  - y0)* sin( gamma_b4  )-(b4s(1,3)  - z0)*sin( beta_b4  ) ) * ( (b4s(1,1)  - x0)*sin( beta_b4  )-(b4s(1,2)  - y0)* sin( alpha_b4  )) -a_b4  * (b4s(1,3)  - z0)*(b4s(1,1)  - x0) ;
  
a46(1,1) = a46_c1 		;
a46(1,2) = a46_c2 		;
a46(1,3) = a46_c3 		;
a46(1,4) = a46_c4 		;
a46(1,5) = a46_c5 		;
a46(1,6) = a46_c6 		;
a46(1,7) = a46_c7 		;
a46(1,8) = a46_c8 		;
a46(1,9) = a46_c9 		;
a46(1,10)= a46_c10		;
a46(1,11)= a46_c11		;
a46(1,12)= a46_c12		;
a46(1,13)= a46_c13		;
a46(1,14)= a46_c14		;
a46(1,15)= a46_c15		;
a46(1,16)= a46_c16		;
a46(1,17)= a46_c17		;
a46(1,18)= a46_c18		;
a46(1,19)= a46_c19		;
a46(1,20)= a46_c20		;
a46(1,21)= a46_p1 		;
a46(1,22)= a46_p2 		;
a46(1,23)= a46_p3 		;
a46(1,24)= a46_p4 		;
a46(1,25)= a46_b1 		;
a46(1,26)= a46_b2 		;
a46(1,27)= a46_b3 		;
a46(1,28)= a46_b4 		;

a64_c1 		= a46_c1 		;
a64_c2 		= a46_c2 		;
a64_c3 		= a46_c3 		;
a64_c4 		= a46_c4 		;
a64_c5 		= a46_c5 		;
a64_c6 		= a46_c6 		;
a64_c7 		= a46_c7 		;
a64_c8 		= a46_c8 		;
a64_c9 		= a46_c9 		;
a64_c10		= a46_c10		;
a64_c11		= a46_c11		;
a64_c12		= a46_c12		;
a64_c13		= a46_c13		;
a64_c14		= a46_c14		;
a64_c15		= a46_c15		;
a64_c16		= a46_c16		;
a64_c17		= a46_c17		;
a64_c18		= a46_c18		;
a64_c19		= a46_c19		;
a64_c20		= a46_c20		;
a64_p1 		= a46_p1 		;
a64_p2 		= a46_p2 		;
a64_p3 		= a46_p3 		;
a64_p4 		= a46_p4 		;
a64_b1 		= a46_b1 		;
a64_b2 		= a46_b2 		;
a64_b3 		= a46_b3 		;
a64_b4 		= a46_b4 		;
 
a64(1,1) = a64_c1 		;
a64(1,2) = a64_c2 		;
a64(1,3) = a64_c3 		;
a64(1,4) = a64_c4 		;
a64(1,5) = a64_c5 		;
a64(1,6) = a64_c6 		;
a64(1,7) = a64_c7 		;
a64(1,8) = a64_c8 		;
a64(1,9) = a64_c9 		;
a64(1,10)= a64_c10		;
a64(1,11)= a64_c11		;
a64(1,12)= a64_c12		;
a64(1,13)= a64_c13		;
a64(1,14)= a64_c14		;
a64(1,15)= a64_c15		;
a64(1,16)= a64_c16		;
a64(1,17)= a64_c17		;
a64(1,18)= a64_c18		;
a64(1,19)= a64_c19		;
a64(1,20)= a64_c20		;
a64(1,21)= a64_p1 		;
a64(1,22)= a64_p2 		;
a64(1,23)= a64_p3 		;
a64(1,24)= a64_p4 		;
a64(1,25)= a64_b1 		;
a64(1,26)= a64_b2 		;
a64(1,27)= a64_b3 		;
a64(1,28)= a64_b4 		;

% A nested for loop is written below to fill elements of added mass coefficients to the added mass matrix. Does not work
%	for i = 1:1:6
%      for j = 1:1:6
%           a(i,j) = sum(a);
%      end
%   end

a(1,1)=	sum	(a11)	;
a(1,2)=	sum	(a12)	;
a(1,3)=	sum	(a13)	;
a(1,4)=	sum	(a14)	;
a(1,5)=	sum	(a15)	;
a(1,6)=	sum	(a16)	;
a(2,1)=	sum	(a21)	;
a(2,2)=	sum	(a22)	;
a(2,3)=	sum	(a23)	;
a(2,4)=	sum	(a24)	;
a(2,5)=	sum	(a25)	;
a(2,6)=	sum	(a26)	;
a(3,1)=	sum	(a31)	;
a(3,2)=	sum	(a32)	;
a(3,3)=	sum	(a33)	;
a(3,4)=	sum	(a34)	;
a(3,5)=	sum	(a35)	;
a(3,6)=	sum	(a36)	;
a(4,1)=	sum	(a41)	;
a(4,2)=	sum	(a42)	;
a(4,3)=	sum	(a43)	;
a(4,4)=	sum	(a44)	;
a(4,5)=	sum	(a45)	;
a(4,6)=	sum	(a46)	;
a(5,1)=	sum	(a51)	;
a(5,2)=	sum	(a52)	;
a(5,3)=	sum	(a53)	;
a(5,4)=	sum	(a54)	;
a(5,5)=	sum	(a55)	;
a(5,6)=	sum	(a56)	;
a(6,1)=	sum	(a61)	;
a(6,2)=	sum	(a62)	;
a(6,3)=	sum	(a63)	;
a(6,4)=	sum	(a64)	;
a(6,5)=	sum	(a65)	;
a(6,6)=	sum	(a66)	;

compare_a(1,1)=	sum	(a11) / g;
compare_a(1,2)=	sum	(a12) / g;
compare_a(1,3)=	sum	(a13) / g;
compare_a(1,4)=	sum	(a14) / g;
compare_a(1,5)=	sum	(a15) / g;
compare_a(1,6)=	sum	(a16) / g;
compare_a(2,1)=	sum	(a21) / g;
compare_a(2,2)=	sum	(a22) / g;
compare_a(2,3)=	sum	(a23) / g;
compare_a(2,4)=	sum	(a24) / g;
compare_a(2,5)=	sum	(a25) / g;
compare_a(2,6)=	sum	(a26) / g;
compare_a(3,1)=	sum	(a31) / g;
compare_a(3,2)=	sum	(a32) / g;
compare_a(3,3)=	sum	(a33) / g;
compare_a(3,4)=	sum	(a34) / g;
compare_a(3,5)=	sum	(a35) / g;
compare_a(3,6)=	sum	(a36) / g;
compare_a(4,1)=	sum	(a41) / g;
compare_a(4,2)=	sum	(a42) / g;
compare_a(4,3)=	sum	(a43) / g;
compare_a(4,4)=	sum	(a44) / g;
compare_a(4,5)=	sum	(a45) / g;
compare_a(4,6)=	sum	(a46) / g;
compare_a(5,1)=	sum	(a51) / g;
compare_a(5,2)=	sum	(a52) / g;
compare_a(5,3)=	sum	(a53) / g;
compare_a(5,4)=	sum	(a54) / g;
compare_a(5,5)=	sum	(a55) / g;
compare_a(5,6)=	sum	(a56) / g;
compare_a(6,1)=	sum	(a61) / g;
compare_a(6,2)=	sum	(a62) / g;
compare_a(6,3)=	sum	(a63) / g;
compare_a(6,4)=	sum	(a64) / g;
compare_a(6,5)=	sum	(a65) / g;
compare_a(6,6)=	sum	(a66) / g;

for i = 1:1:36
    temp1 = abs( compare_a(i) );
	    if  temp1 < 1
	    compare_a(i) = 0;
    end
end


% Below the [restoring coefficients matrix] will be defined.
c = zeros(6,6)	;
c(3,3) = rho*g*pi*(Rc1^2+Rc2^2+Rc3^2+Rc4^2+Rc5^2+Rc6^2+Rc7^2+Rc8^2+...
Rc9^2+Rc10^2+Rc11^2+Rc12^2+Rc13^2+Rc14^2+Rc15^2+Rc16^2+Rc17^2+Rc18^2+Rc19^2+Rc20^2);
%	For calculation of c(3,3), the formula requires input of the 'endpoint' on the water surface.
%	Only the terms piercing the water surface (columns) is included.
%	And in my case, the point on the water surface is denoted like 
c(3,4) = rho*g*pi*( Rc1^2*(c1h(1,2)-y0)+Rc2^2*(c2h(1,2)-y0)+Rc3^2*(c3h(1,2)-y0)+Rc4^2*(c4h(1,2)-y0)+Rc5^2*(c5h(1,2)-y0)+...
					Rc6^2*(c6h(1,2)-y0)+Rc7^2*(c7h(1,2)-y0)+Rc8^2*(c8h(1,2)-y0)+Rc9^2*(c9h(1,2)-y0)+Rc10^2*(c10h(1,2)-y0)+...
					Rc11^2*(c11h(1,2)-y0)+Rc12^2*(c12h(1,2)-y0)+Rc13^2*(c13h(1,2)-y0)+Rc14^2*(c14h(1,2)-y0)+Rc15^2*(c15h(1,2)-y0)+...
					Rc16^2*(c16h(1,2)-y0)+Rc17^2*(c17h(1,2)-y0)+Rc18^2*(c18h(1,2)-y0)+Rc19^2*(c19h(1,2)-y0)+Rc20^2*(c20h(1,2)-y0) );
c(3,5) = -rho*g*pi*( Rc1^2*(c1h(1,1)-x0)+Rc2^2*(c2h(1,1)-x0)+Rc3^2*(c3h(1,1)-x0)+Rc4^2*(c4h(1,1)-x0)+Rc5^2*(c5h(1,1)-x0)+...
					 Rc6^2*(c6h(1,1)-x0)+Rc7^2*(c7h(1,1)-x0)+Rc8^2*(c8h(1,1)-x0)+Rc9^2*(c9h(1,1)-x0)+Rc10^2*(c10h(1,1)-x0)+...
					 Rc11^2*(c11h(1,1)-x0)+Rc12^2*(c12h(1,1)-x0)+Rc13^2*(c13h(1,1)-x0)+Rc14^2*(c14h(1,1)-x0)+Rc15^2*(c15h(1,1)-x0)+...
					 Rc16^2*(c16h(1,1)-x0)+Rc17^2*(c17h(1,1)-x0)+Rc18^2*(c18h(1,1)-x0)+Rc19^2*(c19h(1,1)-x0)+Rc20^2*(c20h(1,1)-x0) )	;
%	The calculation of c(4,4) is partly finished when I calculated gyration of 
%	rotation (located at line 239, or use search function for the keyword 'I4xx').
c(4,4) = -( z0 - cob(1,3) )*dis*rho*g + rho*g*I4xx;
c(5,5) = -( z0 - cob(1,3) )*dis*rho*g + rho*g*I4yy;
c(4,5) = -rho*g*( Rc1^2*(c1h(1,2)-y0)*(c1h(1,1)-x0)+Rc2^2*(c2h(1,2)-y0)*(c2h(1,1)-x0)+Rc3^2*(c3h(1,2)-y0)*(c3h(1,1)-x0)+...
				  Rc4^2*(c4h(1,2)-y0)*(c4h(1,1)-x0)+Rc5^2*(c5h(1,2)-y0)*(c5h(1,1)-x0)+Rc6^2*(c6h(1,2)-y0)*(c6h(1,1)-x0)+...
				  Rc7^2*(c7h(1,2)-y0)*(c7h(1,1)-x0)+Rc8^2*(c8h(1,2)-y0)*(c8h(1,1)-x0)+Rc9^2*(c9h(1,2)-y0)*(c9h(1,1)-x0)+...
				  Rc10^2*(c10h(1,2)-y0)*(c10h(1,1)-x0)+Rc11^2*(c11h(1,2)-y0)*(c11h(1,1)-x0)+Rc12^2*(c12h(1,2)-y0)*(c12h(1,1)-x0)+...
				  Rc13^2*(c13h(1,2)-y0)*(c13h(1,1)-x0)+Rc14^2*(c14h(1,2)-y0)*(c14h(1,1)-x0)+Rc15^2*(c15h(1,2)-y0)*(c15h(1,1)-x0)+...
				  Rc16^2*(c16h(1,2)-y0)*(c16h(1,1)-x0)+Rc17^2*(c17h(1,2)-y0)*(c17h(1,1)-x0)+Rc18^2*(c18h(1,2)-y0)*(c18h(1,1)-x0)+...
				  Rc19^2*(c19h(1,2)-y0)*(c19h(1,1)-x0)+Rc20^2*(c20h(1,2)-y0)*(c20h(1,1)-x0) );
c(5,4) = c(4,5);

%for j = 1:1:36
%    temp2 = abs( compare_c(j) );
%	    if  temp2 < 1
%	    compare_c(i) = 0;
%    end
%end



% [Wave Excited Forces]
% Note: waves are long-crested, regualr, have the shape of a cos function.
% Note: System of coordiantes of water particle motions is defined with (p,q,r); waves could attack the hull at an angle of miu
miu=pi/6;
% miu = 0:pi/20:pi; 	miu could also be definedd with an array, for 
% simplicity, the wave angle is only one number, randomed at 30 degrees


% This follwoing document/part is defined to solve transendental functions.
% In particular to solve the wave length given frequency of waves using the
% dispersion relationship.

omega = 0.1:0.03:1.3;
i2 = size(omega);
h = 2000;	%water depth
lamda = zeros(1,i2(1,2));
for i = 1:1:i2(1,2)
% 	syms wavelength
	eqn =tanh(2*pi*h/wavelength) * (2*pi*g/wavelength) == (omega(1,i))^2;
	solx = solve(eqn,wavelength) 
eqn =@(wavelength)tanh(2*pi*h/wavelength) * (2*pi*g/wavelength) - (omega(1,i))^2;
 	solx = fsolve(eqn,2);
	
	if solx < 0
		solx = -solx;
	end
	lamda(1,i) = solx;
	
end
 
% lamda = int16( lamda );
% Thought process as follows, 1-wave frequencies are defined with
% omega; 2-i2 is a temporary variable designed to hold the size information
% of omega; 3-wavelength is defined as a variable;4-if statement is defined
% to make sure all wavelength values are positive; 5-lamda is the target
% matrix to hold all wavelength values, should be the same size as the
% omega matrix; 6-solve the transendental function, ignore the warnings  


% Wave force calculations are further divided into three parts. 
% 	1-parts that are small;
% 	2-parts that are big relative to wave length;
% 	3-surfaces

% [page 112, part that are small] are not present in this model
% type two and three are present, so the following part of code will calculate forces on cylinders without the ends
% and the ends, which are surfaces.

% [Wave Force Part One]-forces on cylinders without the ends (surfaces):
% For cylinders(columns, braces and pontoons), their length are generally big compared with wave length. 
% Below, the endpoints will be defined in the wave reference coordinate system (p,q,r), in the previous section where
% endpoint positions are defined, a typical point is called c14l(1,1), now a new set of points location will be defined
% and this time, the c will be made upper case for ease of differentiation.

% 	End points for [pontoons] [4 pontoons in total]
%	The points denoted 'a' takes the place of 'd1' for (xyz) coordinate system
%	The points denoted 'b' takes the place of 'd2' for (xyz) coordinate system
%	While the points for (pqr) coordinate system will continue using the a&b notation.
P1a(1,1) = p1a(1,1)*cos(miu) + p1a(1,2)*sin(miu)	;	P1a(1,2) = -p1a(1,1)*sin(miu) + p1a(1,2)*cos(miu);	P1a(1,3) = p1a(1,3)		;
P2a(1,1) = p2a(1,1)*cos(miu) + p2a(1,2)*sin(miu)	;	P2a(1,2) = -p2a(1,1)*sin(miu) + p2a(1,2)*cos(miu);	P2a(1,3) = p2a(1,3)		;
P3a(1,1) = p3a(1,1)*cos(miu) + p3a(1,2)*sin(miu)	;	P3a(1,2) = -p3a(1,1)*sin(miu) + p3a(1,2)*cos(miu);	P3a(1,3) = p3a(1,3)		;
P4a(1,1) = p4a(1,1)*cos(miu) + p4a(1,2)*sin(miu)	;	P4a(1,2) = -p4a(1,1)*sin(miu) + p4a(1,2)*cos(miu);	P4a(1,3) = p4a(1,3)		;
P1b(1,1) = p1b(1,1)*cos(miu) + p1b(1,2)*sin(miu)	;	P1b(1,2) = -p1b(1,1)*sin(miu) + p1b(1,2)*cos(miu);	P1b(1,3) = p1b(1,3)		;
P2b(1,1) = p2b(1,1)*cos(miu) + p2b(1,2)*sin(miu)	;	P2b(1,2) = -p2b(1,1)*sin(miu) + p2b(1,2)*cos(miu);	P2b(1,3) = p2b(1,3)		;
P3b(1,1) = p3b(1,1)*cos(miu) + p3b(1,2)*sin(miu)	;	P3b(1,2) = -p3b(1,1)*sin(miu) + p3b(1,2)*cos(miu);	P3b(1,3) = p3b(1,3)		;
P4b(1,1) = p4b(1,1)*cos(miu) + p4b(1,2)*sin(miu)	;	P4b(1,2) = -p4b(1,1)*sin(miu) + p4b(1,2)*cos(miu);	P4b(1,3) = p4b(1,3)		;
% 	End points for [columns] [20 columns in total]
%	The 'h' in the denotation means high, closer to the deck. 
C1h(1,1)	= c1h(1,1)	*cos(miu) + c1h(1,2)*sin(miu)	;	 C1h(1,2)	= -c1h(1,1)	*sin(miu) + c1h(1,2)*cos(miu)	;   C1h(1,3)	= c1h(1,3);
C2h(1,1)	= c2h(1,1)	*cos(miu) + c2h(1,2)*sin(miu)	;	 C2h(1,2)	= -c2h(1,1)	*sin(miu) + c2h(1,2)*cos(miu)	;   C2h(1,3)	= c2h(1,3);
C3h(1,1)	= c3h(1,1)	*cos(miu) + c3h(1,2)*sin(miu)	;	 C3h(1,2)	= -c3h(1,1)	*sin(miu) + c3h(1,2)*cos(miu)	;   C3h(1,3)	= c3h(1,3);
C4h(1,1)	= c4h(1,1)	*cos(miu) + c4h(1,2)*sin(miu)   ;	 C4h(1,2)	= -c4h(1,1)	*sin(miu) + c4h(1,2)*cos(miu)   ;   C4h(1,3)	= c4h(1,3);
C5h(1,1)	= c5h(1,1)	*cos(miu) + c5h(1,2)*sin(miu)   ;	 C5h(1,2)	= -c5h(1,1)	*sin(miu) + c5h(1,2)*cos(miu)   ;   C5h(1,3)	= c5h(1,3);
C6h(1,1)	= c6h(1,1)	*cos(miu) + c6h(1,2)*sin(miu)	;	 C6h(1,2)	= -c6h(1,1)	*sin(miu) + c6h(1,2)*cos(miu)	;   C6h(1,3)	= c6h(1,3);
C7h(1,1)	= c7h(1,1)	*cos(miu) + c7h(1,2)*sin(miu)	;	 C7h(1,2)	= -c7h(1,1)	*sin(miu) + c7h(1,2)*cos(miu)	;   C7h(1,3)	= c7h(1,3);
C8h(1,1)	= c8h(1,1)	*cos(miu) + c8h(1,2)*sin(miu)	;	 C8h(1,2)	= -c8h(1,1)	*sin(miu) + c8h(1,2)*cos(miu)	;   C8h(1,3)	= c8h(1,3);
C9h(1,1)	= c9h(1,1)	*cos(miu) + c9h(1,2)*sin(miu)	;	 C9h(1,2)	= -c9h(1,1)	*sin(miu) + c9h(1,2)*cos(miu)	;   C9h(1,3)	= c9h(1,3);
C10h(1,1)	= c10h(1,1)	*cos(miu) + c10h(1,2)*sin(miu)	;	 C10h(1,2)	= -c10h(1,1)*sin(miu) + c10h(1,2)*cos(miu)	;   C10h(1,3)	= c10h(1,3);
C11h(1,1)	= c11h(1,1)	*cos(miu) + c11h(1,2)*sin(miu)	;	 C11h(1,2)	= -c11h(1,1)*sin(miu) + c11h(1,2)*cos(miu)	;   C11h(1,3)	= c11h(1,3);
C12h(1,1)	= c12h(1,1)	*cos(miu) + c12h(1,2)*sin(miu)	;	 C12h(1,2)	= -c12h(1,1)*sin(miu) + c12h(1,2)*cos(miu)	;   C12h(1,3)	= c12h(1,3);
C13h(1,1)	= c13h(1,1)	*cos(miu) + c13h(1,2)*sin(miu)	;	 C13h(1,2)	= -c13h(1,1)*sin(miu) + c13h(1,2)*cos(miu)	;   C13h(1,3)	= c13h(1,3);
C14h(1,1)	= c14h(1,1)	*cos(miu) + c14h(1,2)*sin(miu) 	;	 C14h(1,2)	= -c14h(1,1)*sin(miu) + c14h(1,2)*cos(miu) 	;   C14h(1,3)	= c14h(1,3);
C15h(1,1)	= c15h(1,1)	*cos(miu) + c15h(1,2)*sin(miu) 	;	 C15h(1,2)	= -c15h(1,1)*sin(miu) + c15h(1,2)*cos(miu) 	;   C15h(1,3)	= c15h(1,3);
C16h(1,1)	= c16h(1,1)	*cos(miu) + c16h(1,2)*sin(miu)	;	 C16h(1,2)	= -c16h(1,1)*sin(miu) + c16h(1,2)*cos(miu)	;   C16h(1,3)	= c16h(1,3);
C17h(1,1)	= c17h(1,1)	*cos(miu) + c17h(1,2)*sin(miu)	;	 C17h(1,2)	= -c17h(1,1)*sin(miu) + c17h(1,2)*cos(miu)	;   C17h(1,3)	= c17h(1,3);
C18h(1,1)	= c18h(1,1)	*cos(miu) + c18h(1,2)*sin(miu)	;	 C18h(1,2)	= -c18h(1,1)*sin(miu) + c18h(1,2)*cos(miu)	;   C18h(1,3)	= c18h(1,3);
C19h(1,1)	= c19h(1,1)	*cos(miu) + c19h(1,2)*sin(miu)	;	 C19h(1,2)	= -c19h(1,1)*sin(miu) + c19h(1,2)*cos(miu)	;   C19h(1,3)	= c19h(1,3);
C20h(1,1)	= c20h(1,1)	*cos(miu) + c20h(1,2)*sin(miu)	;	 C20h(1,2)	= -c20h(1,1)*sin(miu) + c20h(1,2)*cos(miu)	;   C20h(1,3)	= c20h(1,3);                                                   
C1l(1,1)	= c1l(1,1)	*cos(miu) + c1l(1,2)*sin(miu)	;	 C1l(1,2)	= -c1l(1,1)	*sin(miu) + c1l(1,2)*cos(miu)	;   C1l(1,3)	= c1l(1,3);
C2l(1,1)	= c2l(1,1)	*cos(miu) + c2l(1,2)*sin(miu)	;	 C2l(1,2)	= -c2l(1,1)	*sin(miu) + c2l(1,2)*cos(miu)	;   C2l(1,3)	= c2l(1,3);
C3l(1,1)	= c3l(1,1)	*cos(miu) + c3l(1,2)*sin(miu)	;	 C3l(1,2)	= -c3l(1,1)	*sin(miu) + c3l(1,2)*cos(miu)	;   C3l(1,3)	= c3l(1,3);
C4l(1,1)	= c4l(1,1)	*cos(miu) + c4l(1,2)*sin(miu)   ;	 C4l(1,2)	= -c4l(1,1)	*sin(miu) + c4l(1,2)*cos(miu)   ;   C4l(1,3)	= c4l(1,3);
C5l(1,1)	= c5l(1,1)	*cos(miu) + c5l(1,2)*sin(miu)   ;	 C5l(1,2)	= -c5l(1,1)	*sin(miu) + c5l(1,2)*cos(miu)   ;   C5l(1,3)	= c5l(1,3);
C6l(1,1)	= c6l(1,1)	*cos(miu) + c6l(1,2)*sin(miu)	;	 C6l(1,2)	= -c6l(1,1)	*sin(miu) + c6l(1,2)*cos(miu)	;  	C6l(1,3)	= c6l(1,3);
C7l(1,1)	= c7l(1,1)	*cos(miu) + c7l(1,2)*sin(miu)	;	 C7l(1,2)	= -c7l(1,1)	*sin(miu) + c7l(1,2)*cos(miu)	;  	C7l(1,3)	= c7l(1,3);
C8l(1,1)	= c8l(1,1)	*cos(miu) + c8l(1,2)*sin(miu)	;	 C8l(1,2)	= -c8l(1,1)	*sin(miu) + c8l(1,2)*cos(miu)	;  	C8l(1,3)	= c8l(1,3);
C9l(1,1)	= c9l(1,1)	*cos(miu) + c9l(1,2)*sin(miu)	;	 C9l(1,2)	= -c9l(1,1)	*sin(miu) + c9l(1,2)*cos(miu)	;  	C9l(1,3)	= c9l(1,3);
C10l(1,1)	= c10l(1,1)	*cos(miu) + c10l(1,2)*sin(miu)	;	 C10l(1,2)	= -c10l(1,1)*sin(miu) + c10l(1,2)*cos(miu)	;  	C10l(1,3)	= c10l(1,3);
C11l(1,1)	= c11l(1,1)	*cos(miu) + c11l(1,2)*sin(miu)	;	 C11l(1,2)	= -c11l(1,1)*sin(miu) + c11l(1,2)*cos(miu)	;   C11l(1,3)	= c11l(1,3);
C12l(1,1)	= c12l(1,1)	*cos(miu) + c12l(1,2)*sin(miu)	;	 C12l(1,2)	= -c12l(1,1)*sin(miu) + c12l(1,2)*cos(miu)	;   C12l(1,3)	= c12l(1,3);
C13l(1,1)	= c13l(1,1)	*cos(miu) + c13l(1,2)*sin(miu)	;	 C13l(1,2)	= -c13l(1,1)*sin(miu) + c13l(1,2)*cos(miu)	;   C13l(1,3)	= c13l(1,3);
C14l(1,1)	= c14l(1,1)	*cos(miu) + c14l(1,2)*sin(miu) 	;	 C14l(1,2)	= -c14l(1,1)*sin(miu) + c14l(1,2)*cos(miu) 	;   C14l(1,3)	= c14l(1,3);
C15l(1,1)	= c15l(1,1)	*cos(miu) + c15l(1,2)*sin(miu) 	;	 C15l(1,2)	= -c15l(1,1)*sin(miu) + c15l(1,2)*cos(miu) 	;   C15l(1,3)	= c15l(1,3);
C16l(1,1)	= c16l(1,1)	*cos(miu) + c16l(1,2)*sin(miu)	;	 C16l(1,2)	= -c16l(1,1)*sin(miu) + c16l(1,2)*cos(miu)	;   C16l(1,3)	= c16l(1,3);
C17l(1,1)	= c17l(1,1)	*cos(miu) + c17l(1,2)*sin(miu)	;	 C17l(1,2)	= -c17l(1,1)*sin(miu) + c17l(1,2)*cos(miu)	;   C17l(1,3)	= c17l(1,3);
C18l(1,1)	= c18l(1,1)	*cos(miu) + c18l(1,2)*sin(miu)	;	 C18l(1,2)	= -c18l(1,1)*sin(miu) + c18l(1,2)*cos(miu)	;   C18l(1,3)	= c18l(1,3);
C19l(1,1)	= c19l(1,1)	*cos(miu) + c19l(1,2)*sin(miu)	;	 C19l(1,2)	= -c19l(1,1)*sin(miu) + c19l(1,2)*cos(miu)	;   C19l(1,3)	= c19l(1,3);
C20l(1,1)	= c20l(1,1)	*cos(miu) + c20l(1,2)*sin(miu)	;	 C20l(1,2)	= -c20l(1,1)*sin(miu) + c20l(1,2)*cos(miu)	;   C20l(1,3)	= c20l(1,3);

% 	End points for [Braces] [4 braces in total]
B1p(1,1) = b1p(1,1)*cos(miu) + b1p(1,2)*sin(miu)	;	B1p(1,2) = -b1p(1,1)*sin(miu) + b1p(1,2)*cos(miu)	;	B1p(1,3) = b1p(1,3)		;
B2p(1,1) = b2p(1,1)*cos(miu) + b2p(1,2)*sin(miu)	;	B2p(1,2) = -b2p(1,1)*sin(miu) + b2p(1,2)*cos(miu)	;	B2p(1,3) = b2p(1,3)		;
B3p(1,1) = b3p(1,1)*cos(miu) + b3p(1,2)*sin(miu)	;	B3p(1,2) = -b3p(1,1)*sin(miu) + b3p(1,2)*cos(miu)	;	B3p(1,3) = b3p(1,3)		;
B4p(1,1) = b4p(1,1)*cos(miu) + b4p(1,2)*sin(miu)	;	B4p(1,2) = -b4p(1,1)*sin(miu) + b4p(1,2)*cos(miu)	;	B4p(1,3) = b4p(1,3)		;                                                         
B1s(1,1) = b1s(1,1)*cos(miu) + b1s(1,2)*sin(miu)	;	B1s(1,2) = -b1s(1,1)*sin(miu) + b1s(1,2)*cos(miu)	;	B1s(1,3) = b1s(1,3)		;
B2s(1,1) = b2s(1,1)*cos(miu) + b2s(1,2)*sin(miu)	;	B2s(1,2) = -b2s(1,1)*sin(miu) + b2s(1,2)*cos(miu)	;	B2s(1,3) = b2s(1,3)		;
B3s(1,1) = b3s(1,1)*cos(miu) + b3s(1,2)*sin(miu)	;	B3s(1,2) = -b3s(1,1)*sin(miu) + b3s(1,2)*cos(miu)	;	B3s(1,3) = b3s(1,3)		;
B4s(1,1) = b4s(1,1)*cos(miu) + b4s(1,2)*sin(miu)	;	B4s(1,2) = -b4s(1,1)*sin(miu) + b4s(1,2)*cos(miu)	;	B4s(1,3) = b4s(1,3)		;

% After defining all points it is useful to put them into one matrix for ease of use later. 
% Note that this matrix holds the location information for end points in the (pqr) coordinate system.
% [Matrix name] [P,Q,R]
% Note that P(1,1)  till P(1,20) holds column high points; P(1,21) till P(1,24) holds pontoon aft  points; P(1,25) till P(1,28) holds brace port      points;  
%			P(1,29) till P(1,48) holds column low  points; P(1,49) till P(1,52) holds pontoon fore points; P(1,53) till P(1,56) holds brace starboard points;
P(1,1)		=	C1h(1,1)		;  Q(1,1)		=	C1h(1,2)  		;  R(1,1)		=	C1h(1,3) 		;
P(1,2)		=	C2h(1,1)		;  Q(1,2)		=	C2h(1,2)  		;  R(1,2)		=	C2h(1,3) 		;
P(1,3)		=	C3h(1,1)		;  Q(1,3)		=	C3h(1,2)  		;  R(1,3)		=	C3h(1,3) 		;
P(1,4)		=	C4h(1,1)		;  Q(1,4)		=	C4h(1,2)  		;  R(1,4)		=	C4h(1,3) 		;
P(1,5)		=	C5h(1,1)		;  Q(1,5)		=	C5h(1,2)  		;  R(1,5)		=	C5h(1,3) 		;
P(1,6)		=	C6h(1,1)		;  Q(1,6)		=	C6h(1,2)  		;  R(1,6)		=	C6h(1,3) 		;
P(1,7)		=	C7h(1,1)		;  Q(1,7)		=	C7h(1,2)  		;  R(1,7)		=	C7h(1,3) 		;
P(1,8)		=	C8h(1,1)		;  Q(1,8)		=	C8h(1,2)  		;  R(1,8)		=	C8h(1,3) 		;
P(1,9)		=	C9h(1,1)		;  Q(1,9)		=	C9h(1,2)  		;  R(1,9)		=	C9h(1,3) 		;
P(1,10)		=	C10h(1,1)		;  Q(1,10)		=	C10h(1,2)  		;  R(1,10)		=	C10h(1,3) 			;
P(1,11)		=	C11h(1,1)		;  Q(1,11)		=	C11h(1,2)  		;  R(1,11)		=	C11h(1,3) 			;
P(1,12)		=	C12h(1,1)		;  Q(1,12)		=	C12h(1,2)  		;  R(1,12)		=	C12h(1,3) 			;
P(1,13)		=	C13h(1,1)		;  Q(1,13)		=	C13h(1,2)  		;  R(1,13)		=	C13h(1,3) 			;
P(1,14)		=	C14h(1,1)		;  Q(1,14)		=	C14h(1,2)  		;  R(1,14)		=	C14h(1,3) 			;
P(1,15)		=	C15h(1,1)		;  Q(1,15)		=	C15h(1,2)  		;  R(1,15)		=	C15h(1,3) 			;
P(1,16)		=	C16h(1,1)		;  Q(1,16)		=	C16h(1,2)  		;  R(1,16)		=	C16h(1,3) 			;
P(1,17)		=	C17h(1,1)		;  Q(1,17)		=	C17h(1,2)  		;  R(1,17)		=	C17h(1,3) 			;
P(1,18)		=	C18h(1,1)		;  Q(1,18)		=	C18h(1,2)  		;  R(1,18)		=	C18h(1,3) 			;
P(1,19)		=	C19h(1,1)		;  Q(1,19)		=	C19h(1,2)  		;  R(1,19)		=	C19h(1,3) 			;
P(1,20)		=	C20h(1,1)		;  Q(1,20)		=	C20h(1,2)  		;  R(1,20)		=	C20h(1,3) 			;
P(1,21)		=	P1a(1,1)  		;  Q(1,21)		=	P1a(1,2)   		;  R(1,21)		=	P1a(1,3)  			;
P(1,22)		=	P2a(1,1)  		;  Q(1,22)		=	P2a(1,2)   		;  R(1,22)		=	P2a(1,3)  			;
P(1,23)		=	P3a(1,1)  		;  Q(1,23)		=	P3a(1,2)   		;  R(1,23)		=	P3a(1,3)  			;
P(1,24)		=	P4a(1,1)  		;  Q(1,24)		=	P4a(1,2)   		;  R(1,24)		=	P4a(1,3)  			;
P(1,25)		=	B1p(1,1)  		;  Q(1,25)		=	B1p(1,2)  		;  R(1,25)		=	B1p(1,3)  			;
P(1,26)		=	B2p(1,1)  		;  Q(1,26)		=	B2p(1,2)  		;  R(1,26)		=	B2p(1,3)  			;
P(1,27)		=	B3p(1,1)  		;  Q(1,27)		=	B3p(1,2)  		;  R(1,27)		=	B3p(1,3)  			;
P(1,28)		=	B4p(1,1)  		;  Q(1,28)		=	B4p(1,2)  		;  R(1,28)		=	B4p(1,3)  			;
P(1,29)		=	C1l(1,1) 		;  Q(1,29)		=	C1l(1,2)  		;  R(1,29)		=	C1l(1,3) 		;
P(1,30)		=	C2l(1,1) 		;  Q(1,30)		=	C2l(1,2)  		;  R(1,30)		=	C2l(1,3) 		;
P(1,31)		=	C3l(1,1) 		;  Q(1,31)		=	C3l(1,2)  		;  R(1,31)		=	C3l(1,3) 		;
P(1,32)		=	C4l(1,1) 		;  Q(1,32)		=	C4l(1,2)  		;  R(1,32)		=	C4l(1,3) 		;
P(1,33)		=	C5l(1,1) 		;  Q(1,33)		=	C5l(1,2)  		;  R(1,33)		=	C5l(1,3) 		;
P(1,34)		=	C6l(1,1) 		;  Q(1,34)		=	C6l(1,2)  		;  R(1,34)		=	C6l(1,3) 		;
P(1,35)		=	C7l(1,1) 		;  Q(1,35)		=	C7l(1,2)  		;  R(1,35)		=	C7l(1,3) 		;
P(1,36)		=	C8l(1,1) 		;  Q(1,36)		=	C8l(1,2)  		;  R(1,36)		=	C8l(1,3) 		;
P(1,37)		=	C9l(1,1) 		;  Q(1,37)		=	C9l(1,2)  		;  R(1,37)		=	C9l(1,3) 		;
P(1,38)		=	C10l(1,1) 		;  Q(1,38)		=	C10l(1,2)  		;  R(1,38)		=	C10l(1,3) 			;
P(1,39)		=	C11l(1,1) 		;  Q(1,39)		=	C11l(1,2)  		;  R(1,39)		=	C11l(1,3) 			;
P(1,40)		=	C12l(1,1) 		;  Q(1,40)		=	C12l(1,2)  		;  R(1,40)		=	C12l(1,3) 			;
P(1,41)		=	C13l(1,1) 		;  Q(1,41)		=	C13l(1,2)  		;  R(1,41)		=	C13l(1,3) 			;
P(1,42)		=	C14l(1,1) 		;  Q(1,42)		=	C14l(1,2)  		;  R(1,42)		=	C14l(1,3) 			;
P(1,43)		=	C15l(1,1) 		;  Q(1,43)		=	C15l(1,2)  		;  R(1,43)		=	C15l(1,3) 			;
P(1,44)		=	C16l(1,1) 		;  Q(1,44)		=	C16l(1,2)  		;  R(1,44)		=	C16l(1,3) 			;
P(1,45)		=	C17l(1,1) 		;  Q(1,45)		=	C17l(1,2)  		;  R(1,45)		=	C17l(1,3) 			;
P(1,46)		=	C18l(1,1) 		;  Q(1,46)		=	C18l(1,2)  		;  R(1,46)		=	C18l(1,3) 			;
P(1,47)		=	C19l(1,1) 		;  Q(1,47)		=	C19l(1,2)  		;  R(1,47)		=	C19l(1,3) 			;
P(1,48)		=	C20l(1,1) 		;  Q(1,48)		=	C20l(1,2)  		;  R(1,48)		=	C20l(1,3) 			;
P(1,49)		=	P1b(1,1)  		;  Q(1,49)		=	P1b(1,2)   		;  R(1,49)		=	P1b(1,3)  			;
P(1,50)		=	P2b(1,1)  		;  Q(1,50)		=	P2b(1,2)   		;  R(1,50)		=	P2b(1,3)  			;
P(1,51)		=	P3b(1,1)  		;  Q(1,51)		=	P3b(1,2)   		;  R(1,51)		=	P3b(1,3)  			;
P(1,52)		=	P4b(1,1)  		;  Q(1,52)		=	P4b(1,2)   		;  R(1,52)		=	P4b(1,3)  			;
P(1,53)		=	B1s(1,1)  		;  Q(1,53)		=	B1s(1,2)   		;  R(1,53)		=	B1s(1,3)  			;
P(1,54)		=	B2s(1,1)  		;  Q(1,54)		=	B2s(1,2)   		;  R(1,54)		=	B2s(1,3)  			;
P(1,55)		=	B3s(1,1)  		;  Q(1,55)		=	B3s(1,2)   		;  R(1,55)		=	B3s(1,3)  			;
P(1,56)		=	B4s(1,1)  		;  Q(1,56)		=	B4s(1,2)   		;  R(1,56)		=	B4s(1,3)  			;
 
%	All Angles for (xyz) coordinate system are Listed Below.
%	all angles are renamed with a 'w' before the old angle names (which means wave).
walpha_c1  = asin( (C1l(1,1) -  C1h(1,1) ) / Lc1 )	;
walpha_c2  = asin( (C2l(1,1) -  C2h(1,1) ) / Lc2 )	;
walpha_c3  = asin( (C3l(1,1) -  C3h(1,1) ) / Lc3 )	;
walpha_c4  = asin( (C4l(1,1) -  C4h(1,1) ) / Lc4 )	;
walpha_c5  = asin( (C5l(1,1) -  C5h(1,1) ) / Lc5 )	;
walpha_c6  = asin( (C6l(1,1) -  C6h(1,1) ) / Lc6 )	;
walpha_c7  = asin( (C7l(1,1) -  C7h(1,1) ) / Lc7 )	;
walpha_c8  = asin( (C8l(1,1) -  C8h(1,1) ) / Lc8 )	;
walpha_c9  = asin( (C9l(1,1) -  C9h(1,1) ) / Lc9 )	;
walpha_c10 = asin( (C10l(1,1) - C10h(1,1) ) / Lc10)	;
walpha_c11 = asin( (C11l(1,1) - C11h(1,1) ) / Lc11)	;
walpha_c12 = asin( (C12l(1,1) - C12h(1,1) ) / Lc12)	;
walpha_c13 = asin( (C13l(1,1) - C13h(1,1) ) / Lc13)	;
walpha_c14 = asin( (C14l(1,1) - C14h(1,1) ) / Lc14)	;
walpha_c15 = asin( (C15l(1,1) - C15h(1,1) ) / Lc15)	;
walpha_c16 = asin( (C16l(1,1) - C16h(1,1) ) / Lc16)	;
walpha_c17 = asin( (C17l(1,1) - C17h(1,1) ) / Lc17)	;
walpha_c18 = asin( (C18l(1,1) - C18h(1,1) ) / Lc18)	;
walpha_c19 = asin( (C19l(1,1) - C19h(1,1) ) / Lc19)	;
walpha_c20 = asin( (C20l(1,1) - C20h(1,1) ) / Lc20)	;
walpha_p1  = asin( (P1b(1,1) -  P1a(1,1) ) / Lp1 )	;
walpha_p2  = asin( (P2b(1,1) -  P2a(1,1) ) / Lp2 )	;
walpha_p3  = asin( (P3b(1,1) -  P3a(1,1) ) / Lp3 )	;
walpha_p4  = asin( (P4b(1,1) -  P4a(1,1) ) / Lp4 )	;
walpha_b1  = asin( (B1p(1,1) -  B1s(1,1) ) / Lb1 )	;
walpha_b2  = asin( (B2p(1,1) -  B2s(1,1) ) / Lb2 )	;
walpha_b3  = asin( (B3p(1,1) -  B3s(1,1) ) / Lb3 )	;
walpha_b4  = asin( (B4p(1,1) -  B4s(1,1) ) / Lb4 )	;
wbeta_c1   = asin( (C1l(1,2) -  C1h(1,2) ) / Lc1 )	;
wbeta_c2   = asin( (C2l(1,2) -  C2h(1,2) ) / Lc2 )	;
wbeta_c3   = asin( (C3l(1,2) -  C3h(1,2) ) / Lc3 )	;
wbeta_c4   = asin( (C4l(1,2) -  C4h(1,2) ) / Lc4 )	;
wbeta_c5   = asin( (C5l(1,2) -  C5h(1,2) ) / Lc5 )	;
wbeta_c6   = asin( (C6l(1,2) -  C6h(1,2) ) / Lc6 )	;
wbeta_c7   = asin( (C7l(1,2) -  C7h(1,2) ) / Lc7 )	;
wbeta_c8   = asin( (C8l(1,2) -  C8h(1,2) ) / Lc8 )	;
wbeta_c9   = asin( (C9l(1,2) -  C9h(1,2) ) / Lc9 )	;
wbeta_c10  = asin( (C10l(1,2) - C10h(1,2) ) / Lc10)	;
wbeta_c11  = asin( (C11l(1,2) - C11h(1,2) ) / Lc11)	;
wbeta_c12  = asin( (C12l(1,2) - C12h(1,2) ) / Lc12)	;
wbeta_c13  = asin( (C13l(1,2) - C13h(1,2) ) / Lc13)	;
wbeta_c14  = asin( (C14l(1,2) - C14h(1,2) ) / Lc14)	;
wbeta_c15  = asin( (C15l(1,2) - C15h(1,2) ) / Lc15)	;
wbeta_c16  = asin( (C16l(1,2) - C16h(1,2) ) / Lc16)	;
wbeta_c17  = asin( (C17l(1,2) - C17h(1,2) ) / Lc17)	;
wbeta_c18  = asin( (C18l(1,2) - C18h(1,2) ) / Lc18)	;
wbeta_c19  = asin( (C19l(1,2) - C19h(1,2) ) / Lc19)	;
wbeta_c20  = asin( (C20l(1,2) - C20h(1,2) ) / Lc20)	;
wbeta_p1   = asin( (P1b(1,2) -  P1a(1,2) ) / Lp1 )	;
wbeta_p2   = asin( (P2b(1,2) -  P2a(1,2) ) / Lp2 )	;
wbeta_p3   = asin( (P3b(1,2) -  P3a(1,2) ) / Lp3 )	;
wbeta_p4   = asin( (P4b(1,2) -  P4a(1,2) ) / Lp4 )	;
wbeta_b1   = asin( (B1p(1,2) -  B1s(1,2) ) / Lb1 )	;
wbeta_b2   = asin( (B2p(1,2) -  B2s(1,2) ) / Lb2 )	;
wbeta_b3   = asin( (B3p(1,2) -  B3s(1,2) ) / Lb3 )	;
wbeta_b4   = asin( (B4p(1,2) -  B4s(1,2) ) / Lb4 )	;      
wgamma_c1  = asin( (C1l(1,3) -  C1h(1,3) ) / Lc1 )	;
wgamma_c2  = asin( (C2l(1,3) -  C2h(1,3) ) / Lc2 )	;
wgamma_c3  = asin( (C3l(1,3) -  C3h(1,3) ) / Lc3 )	;
wgamma_c4  = asin( (C4l(1,3) -  C4h(1,3) ) / Lc4 )	;
wgamma_c5  = asin( (C5l(1,3) -  C5h(1,3) ) / Lc5 )	;
wgamma_c6  = asin( (C6l(1,3) -  C6h(1,3) ) / Lc6 )	;
wgamma_c7  = asin( (C7l(1,3) -  C7h(1,3) ) / Lc7 )	;
wgamma_c8  = asin( (C8l(1,3) -  C8h(1,3) ) / Lc8 )	;
wgamma_c9  = asin( (C9l(1,3) -  C9h(1,3) ) / Lc9 )	;
wgamma_c10 = asin( (C10l(1,3) - C10h(1,3) ) / Lc10)	;
wgamma_c11 = asin( (C11l(1,3) - C11h(1,3) ) / Lc11)	;
wgamma_c12 = asin( (C12l(1,3) - C12h(1,3) ) / Lc12)	;
wgamma_c13 = asin( (C13l(1,3) - C13h(1,3) ) / Lc13)	;
wgamma_c14 = asin( (C14l(1,3) - C14h(1,3) ) / Lc14)	;
wgamma_c15 = asin( (C15l(1,3) - C15h(1,3) ) / Lc15)	;
wgamma_c16 = asin( (C16l(1,3) - C16h(1,3) ) / Lc16)	;
wgamma_c17 = asin( (C17l(1,3) - C17h(1,3) ) / Lc17)	;
wgamma_c18 = asin( (C18l(1,3) - C18h(1,3) ) / Lc18)	;
wgamma_c19 = asin( (C19l(1,3) - C19h(1,3) ) / Lc19)	;
wgamma_c20 = asin( (C20l(1,3) - C20h(1,3) ) / Lc20)	;
wgamma_p1  = asin( (P1b(1,3) -  P1a(1,3) ) / Lp1 )	;
wgamma_p2  = asin( (P2b(1,3) -  P2a(1,3) ) / Lp2 )	;
wgamma_p3  = asin( (P3b(1,3) -  P3a(1,3) ) / Lp3 )	;
wgamma_p4  = asin( (P4b(1,3) -  P4a(1,3) ) / Lp4 )	;
wgamma_b1  = asin( (B1p(1,3) -  B1s(1,3) ) / Lb1 )	;
wgamma_b2  = asin( (B2p(1,3) -  B2s(1,3) ) / Lb2 )	;
wgamma_b3  = asin( (B3p(1,3) -  B3s(1,3) ) / Lb3 )	;
wgamma_b4  = asin( (B4p(1,3) -  B4s(1,3) ) / Lb4 )	;

% After defining all angles, they are put together into one matrix named wgamma, walpha and wbeta. The size of these matrices are
% 1 * 28.
walpha(1,1)		=	walpha_c1  		;  wbeta(1,1)		=	wbeta_c1   		;  wgamma(1,1)		=	wgamma_c1  		;
walpha(1,2)		=	walpha_c2  		;  wbeta(1,2)		=	wbeta_c2   		;  wgamma(1,2)		=	wgamma_c2  		;
walpha(1,3)		=	walpha_c3  		;  wbeta(1,3)		=	wbeta_c3   		;  wgamma(1,3)		=	wgamma_c3  		;
walpha(1,4)		=	walpha_c4  		;  wbeta(1,4)		=	wbeta_c4   		;  wgamma(1,4)		=	wgamma_c4  		;
walpha(1,5)		=	walpha_c5  		;  wbeta(1,5)		=	wbeta_c5   		;  wgamma(1,5)		=	wgamma_c5  		;
walpha(1,6)		=	walpha_c6  		;  wbeta(1,6)		=	wbeta_c6   		;  wgamma(1,6)		=	wgamma_c6  		;
walpha(1,7)		=	walpha_c7  		;  wbeta(1,7)		=	wbeta_c7   		;  wgamma(1,7)		=	wgamma_c7  		;
walpha(1,8)		=	walpha_c8  		;  wbeta(1,8)		=	wbeta_c8   		;  wgamma(1,8)		=	wgamma_c8  		;
walpha(1,9)		=	walpha_c9  		;  wbeta(1,9)		=	wbeta_c9   		;  wgamma(1,9)		=	wgamma_c9  		;
walpha(1,10)	=	walpha_c10 		;  wbeta(1,10)		=	wbeta_c10  		;  wgamma(1,10)		=	wgamma_c10 			;
walpha(1,11)	=	walpha_c11 		;  wbeta(1,11)		=	wbeta_c11  		;  wgamma(1,11)		=	wgamma_c11 			;
walpha(1,12)	=	walpha_c12 		;  wbeta(1,12)		=	wbeta_c12  		;  wgamma(1,12)		=	wgamma_c12 			;
walpha(1,13)	=	walpha_c13 		;  wbeta(1,13)		=	wbeta_c13  		;  wgamma(1,13)		=	wgamma_c13 			;
walpha(1,14)	=	walpha_c14 		;  wbeta(1,14)		=	wbeta_c14  		;  wgamma(1,14)		=	wgamma_c14 			;
walpha(1,15)	=	walpha_c15 		;  wbeta(1,15)		=	wbeta_c15  		;  wgamma(1,15)		=	wgamma_c15 			;
walpha(1,16)	=	walpha_c16 		;  wbeta(1,16)		=	wbeta_c16  		;  wgamma(1,16)		=	wgamma_c16 			;
walpha(1,17)	=	walpha_c17 		;  wbeta(1,17)		=	wbeta_c17  		;  wgamma(1,17)		=	wgamma_c17 			;
walpha(1,18)	=	walpha_c18 		;  wbeta(1,18)		=	wbeta_c18  		;  wgamma(1,18)		=	wgamma_c18 			;
walpha(1,19)	=	walpha_c19 		;  wbeta(1,19)		=	wbeta_c19  		;  wgamma(1,19)		=	wgamma_c19 			;
walpha(1,20)	=	walpha_c20 		;  wbeta(1,20)		=	wbeta_c20  		;  wgamma(1,20)		=	wgamma_c20 			;
walpha(1,21)	=	walpha_p1  		;  wbeta(1,21)		=	wbeta_p1   		;  wgamma(1,21)		=	wgamma_p1  			;
walpha(1,22)	=	walpha_p2  		;  wbeta(1,22)		=	wbeta_p2   		;  wgamma(1,22)		=	wgamma_p2  			;
walpha(1,23)	=	walpha_p3  		;  wbeta(1,23)		=	wbeta_p3   		;  wgamma(1,23)		=	wgamma_p3  			;
walpha(1,24)	=	walpha_p4  		;  wbeta(1,24)		=	wbeta_p4   		;  wgamma(1,24)		=	wgamma_p4  			;
walpha(1,25)	=	walpha_b1  		;  wbeta(1,25)		=	wbeta_b1   		;  wgamma(1,25)		=	wgamma_b1  			;
walpha(1,26)	=	walpha_b2  		;  wbeta(1,26)		=	wbeta_b2   		;  wgamma(1,26)		=	wgamma_b2  			;
walpha(1,27)	=	walpha_b3  		;  wbeta(1,27)		=	wbeta_b3   		;  wgamma(1,27)		=	wgamma_b3  			;
walpha(1,28)	=	walpha_b4  		;  wbeta(1,28)		=	wbeta_b4   		;  wgamma(1,28)		=	wgamma_b4  			;


% Below, the Id numbers are defined.
% Id numbers are obtained with Matlab integration function - int. 
% Sth to say about these integratoins: 1-variable of integration is r; 2-integration upper and lower limits are both real numbers;
% 3-all other variables are already determined previously;4-different elements may have a series of different Id values depending on 
% the location of the cylinder as well as their dimensions. Another thing to consider is the kappa number, a.k.a. the wave number,
% wave number are pre-defined with an array. And in theory, each wave number will produce one Id number. 
% 
% The name 'I_c1' was already used during definition of added mass coefficients, so, the variable names of this part will be
% lower case i_c1, i_p1, etc; and i_c1 will be an array the size of wave number (wave frequency) array.
% I also need to calculate all Id1 through Id8 for the element d, so, the naming part will be changed to i1_e01 to i8_e01.
% The 'e' in the denotation means 'element number 1', in the light of this idea, there will be i1_e01 to i1_e28(which in total 
% of 28*8*number of wave frequencies of variables will be calculated).
% 
% [Elements?] [which are which]The first 20 elements are the columns, the 21st to 24th elements are for pontoon number 1 to 4, and the last four are for 
% barces 1 to 4, the same as the angle calculations. 

%[i1 to i8 arraies] are defined below. To think in a multi-dimensional way, I used another method to do this part of the calculation.
%For details, see the 3d-matrix file in the void folder.
%3d-matrix
%	I need to create a whole slew of numbers, and a three dimensional matrix is the perfecr choice for this. 
%	To visiualize this matrix, I can use the metaphor of a book, on page one, is a normal two dimensional 
%	matrix and on page two is another 2d matrix. More pages could be added later. 
%	
%	Creating a 3d matrix could sovle the problem of painstakingly defining all of those matrix names. A normal
%	loop is good at controlling the elements of a matrix instead of creating variable names; once the matrix 
%	is constructed, it would be easier to access those variables as well compared with a series of 2d matrices. 
%	
%	Below is the test I did trying to get a hold of 3d matrix
%	
%	test(:,:,1) =									%first page defined
%	
%	     1     2
%	     3     4
%	
%	
%	test(:,:,2) =									%second page defined
%	
%	     5     6
%	     7     8
%	
%	>> tttt=[11 1;3 8]
%	
%	tttt =											%another matrix defined
%	
%	    11     1
%	     3     8
%	
%	>> tttt.*test(:,:,2)							%here we see it is possible to operate on one page
%	
%	ans =
%	
%	    55     6
%	    21    64

% The name of this matrix is 'II'
% Each page represent one wave number/frequency (these two are linked via the dispersion relationship), 
% in another word, on each page, all i numbers for all 28 elements will be find; and on another page, another 
% wave frequency will be used. 

% On page one, matrix II features the first column of the first row of [matrix lamda] wave number
% II(2,4,1) - element number 4, calculated according to Id4.




syms r
tic
for page = 1:1:i2(1,2)	%i2 is the size of [omag e] array
	kappa(1,page) = 2*pi/lamda(1,page);
	for element_number = 1:1:28
		II(1,element_number,page) = 1/sinh(kappa(1,page)*h)*int( cosh(kappa(1,page)*(h - R(1,element_number) - r * sin(wgamma(1,element_number)))) * cosh(kappa(1,page)*(P(1,element_number) + r * sin(walpha(1,element_number)))), 0, L(1,element_number) );
% 		II(1,element_number,page) = 1/sinh(kappa(1,page)*h)*quadgk(@(r)int_fcn(r,kappa(1,page),h,R(1,element_number),wgamma(1,element_number),walpha(1,element_number),P(1,element_number)),0,L(1,element_number));
% 	    II(2,element_number,page) = 1/sinh(kappa(1,page)*h)*int( cosh(kappa(1,page)*(h - R(1,element_number) - r * sin(wgamma(1,element_number)))) * sinh(kappa(1,page)*(P(1,element_number) + r * sin(walpha(1,element_number)))), 0, L(1,element_number) );
% 	    II(3,element_number,page) = 1/sinh(kappa(1,page)*h)*int( sinh(kappa(1,page)*(h - R(1,element_number) - r * sin(wgamma(1,element_number)))) * cosh(kappa(1,page)*(P(1,element_number) + r * sin(walpha(1,element_number)))), 0, L(1,element_number) );
% 	    II(4,element_number,page) = 1/sinh(kappa(1,page)*h)*int( sinh(kappa(1,page)*(h - R(1,element_number) - r * sin(wgamma(1,element_number)))) * sinh(kappa(1,page)*(P(1,element_number) + r * sin(walpha(1,element_number)))), 0, L(1,element_number) );
% 	    II(5,element_number,page) = 1/sinh(kappa(1,page)*h)*int( cosh(kappa(1,page)*(h - R(1,element_number) - r * sin(wgamma(1,element_number)))) * r * cosh(kappa(1,page)*(P(1,element_number) + r * sin(walpha(1,element_number)))), 0, L(1,element_number) );
% 	    II(6,element_number,page) = 1/sinh(kappa(1,page)*h)*int( cosh(kappa(1,page)*(h - R(1,element_number) - r * sin(wgamma(1,element_number)))) * r * sinh(kappa(1,page)*(P(1,element_number) + r * sin(walpha(1,element_number)))), 0, L(1,element_number) );
% 	    II(7,element_number,page) = 1/sinh(kappa(1,page)*h)*int( sinh(kappa(1,page)*(h - R(1,element_number) - r * sin(wgamma(1,element_number)))) * r * cosh(kappa(1,page)*(P(1,element_number) + r * sin(walpha(1,element_number)))), 0, L(1,element_number) );
% 	    II(8,element_number,page) = 1/sinh(kappa(1,page)*h)*int( sinh(kappa(1,page)*(h - R(1,element_number) - r * sin(wgamma(1,element_number)))) * r * sinh(kappa(1,page)*(P(1,element_number) + r * sin(walpha(1,element_number)))), 0, L(1,element_number) );
	end
end	
toc
% According to page 115, a few variables yet need to be defined.
% On top of page 115, there's definition of Ad, which is defined below.
%	note that to distinguish from variables like a_b1, here they are named with aa_e01 to aa_e28.
%	The "e" here means element. wave height is defined as 1.
%	Added mass for the whole cylinder is the same as the dispalced volume. So, aa_ecc could be calculated as 
%	2 * omega^2 * wave height * rho * (pi*R^2) (among which, O_ecc and omega are stored in previous arraies)
%	This array will be named as "unit_a"
%wave_height = 2;
%for page = 1:1:i2(1,2)	%i2 is the size of [omega] array
%	for element_number = 1:1:28
%		unit_a(element_number,page) = 2*omega(1,page).^2*wave_height*rho*(pi*R(1,element_number).^2)	;
%	end
%end























	
	
	













	
	
	
	
	
	









% [Wave Forces on Cylinder End Planes]