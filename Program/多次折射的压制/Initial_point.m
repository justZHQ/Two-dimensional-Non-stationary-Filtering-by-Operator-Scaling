function t_start=Initial_point(t,x,x_end,x_point)
t_start=[];
nc= between(0,x_end,x,2);
if(nc~=0)
    tlims=[t(1),t(end)];
    xlims=[0,x_end];
	for k=x_point
		tk = interp1(xlims,tlims,k);
        t_start=[t_start;tk];
	end
end
