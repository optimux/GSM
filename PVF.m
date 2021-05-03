%Particle volume fraction in 2D in the region x<a, y<b (NOTE: a^2+b^2<1.0), from Eq.(24-27), Suckale et al., 2012
%2020-10-16

a=0.5;
b=0.0;

%For a*b>0.0
if(a*b>=0.0)
phiabp=2.0*a*b*sqrt(1.0-a*a-b*b)/(4.0*pi)+...
    b*(3.0-b*b)*(asin(a/sqrt(1.0-b*b))+pi/2.0)/(4.0*pi)+...
    a*(3.0-a*a)*(asin(b/sqrt(1.0-a*a))+pi/2.0)/(4.0*pi)+...
    (asin((1.0-a*a-b*b*(1.0+a*a))/((1.0-b*b)*(1.0-a*a)))+pi/2.0)/(4.0*pi)

else
%For a*b<0.0
phiabn=2.0*a*b*sqrt(1.0-a*a-b*b)/(4.0*pi)+...
    b*(3.0-b*b)*(asin(a/sqrt(1.0-b*b))+pi/2.0)/(4.0*pi)+...
    a*(3.0-a*a)*(asin(b/sqrt(1.0-a*a))+pi/2.0)/(4.0*pi)+...
    (-asin((1.0-a*a-b*b*(1.0+a*a))/((1.0-b*b)*(1.0-a*a)))+1.5*pi)/(4.0*pi)
end
