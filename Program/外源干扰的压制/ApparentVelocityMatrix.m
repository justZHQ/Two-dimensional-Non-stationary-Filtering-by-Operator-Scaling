function [ApparentVelocityMatrix]=ApparentVelocityMatrix(t,x,tnot,xnot,v,L,ApparentVelocityMatrix)
v=v/2;
[nsamp,nc]=size(ApparentVelocityMatrix);
dt=t(2)-t(1);
tmin=t(1);
for k=1:nc
	xoff=x(k)-xnot;
		tk = sqrt(tnot^2+(xoff/v)^2);
		ik=(tk-tmin)/dt+1;
		if( between(1,nsamp,ik) )
			ik1=floor(ik);
			ik2=ceil(ik);
            if xoff==0
            ApparentVelocityMatrix(ik1-L:ik2+L,k)=10000;
            else
            ApparentVelocityMatrix(ik1-L:ik2+L,k)=abs(v^2*tk/(xoff));
            end
		end
	end
end
