function ApparentVelocityMatrix=ApparentVelocityMatrix(ApparentVelocityMatrix,t,x,tlims,xlims,L)
nc= between(xlims(1),xlims(2),x,2);
if(nc~=0)
	tmin=t(1);
	dt=t(2)-t(1);
	for k=nc
		tk = interp1(xlims,tlims,x(k));
        tt=floor((tk-tmin)/dt);
        ApparentVelocityMatrix(tt-L:tt+L,k)=abs((xlims(2)-xlims(1))/(tlims(2)-tlims(1)));
	end
end


