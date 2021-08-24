function [Sxx,Sxz,Szz]=EdgeStress(m,x,z,nu,mu)
%EdgeStress     [Sxx,Sxz,Szz]=EdgeStress(m,x,z,nu,{mu})
%
%Computes stresses at (x,z) due to an edge dislocation specified in m.   
%All length inputs should be in same units; otherwise stresses will be off   
%by a constant.  Nu is Poisson's ratio and mu is shear modulus.  If mu is   
%not specified, it defaults to one.   
% 
%The form of m:   
%   
%     m(1) = Depth of updip end    
%     m(2) = Horizontal position of updip end    
%     m(3) = Length downdip   
%     m(4) = Dip (degrees)   
%     m(5) = Slip   
%   
%The vertical axis is reckoned NEGATIVE down.


if nargin <5	mu=1;end   
%Make sure vertical coordinates are negative	
	z=-abs(z);	Z=-abs(m(1));
   
%Define Constants	
	dip=m(4)*pi/180;	
	sd=sin(dip);	
	cd=cos(dip);	
	bx=m(5)*sd;	
	bz=m(5)*cd;	
	a1=[bx -bx]*mu/(2*pi*(1-nu));	
	a2=[bz -bz]*mu/(2*pi*(1-nu));	
	X=[m(2) m(2)+m(3)*cd];	
	Z=[Z Z-m(3)*sd];	
	Sxx=zeros(size(x));	
	Sxz=zeros(size(x));	
	Szz=zeros(size(x));   

%Compute stresses	
for i=1:2		
		zmZ=z-Z(i);
		zmZ2=zmZ.*zmZ;
		zpZ=z+Z(i);
		zpZ2=zpZ.*zpZ;
		xmX=x-X(i);
		xmX2=xmX.*xmX;		
		r12=zmZ2+xmX2;		
		r22=zpZ2+xmX2;		
		r122=r12.*r12;		
		r222=r22.*r22;		
		r223=r222.*r22;		
		b1=xmX.*(3*zmZ2+xmX2)./(r122);
		b2=xmX.*(3*zpZ2+xmX2)./(r222);	
		b3=4.0*Z(i)*z.*xmX.*(3.0*zpZ2-xmX2)./(r223);
		b4=zmZ.*(zmZ2-xmX2)./(r122);
		b5=zpZ.*(zpZ2-xmX2)./(r222);
		a3=(3.0*z+Z(i)).*zpZ2.*zpZ-6.0*z.*zpZ.*xmX2-xmX2.*xmX2;
		b6=2.0*Z(i).*a3./(r223);

		c1=xmX.*(zmZ2-xmX2)./(r122);
		c2=xmX.*(zpZ2-xmX2)./(r222);
		a4=(2.0*Z(i)-z).*zpZ2+(3.0*z+2.0*Z(i)).*xmX2;
		c3=4.0*Z(i).*xmX.*a4./(r223);
		c4=zmZ.*(zmZ2+3.0*xmX2)./(r122);
		c5=zpZ.*(zpZ2+3.0*xmX2)./(r222);
		a5=zmZ.*zpZ2.*zpZ-6.0*z.*zpZ.*xmX2+xmX2.*xmX2;
		c6=2.0*Z(i)*a5./(r223);	

		Sxx=Sxx+a1(i)*(c1-c2+c3)+a2(i)*(-c4+c5+c6);
		Sxz=Sxz+a1(i)*(b4-b5-c6)+a2(i)*(-c1+c2-b3);
		Szz=Szz+a1(i)*(-b1+b2+b3)+a2(i)*(-b4+b5-b6);	
end
