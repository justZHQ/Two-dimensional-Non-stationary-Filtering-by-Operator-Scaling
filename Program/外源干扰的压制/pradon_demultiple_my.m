function [prim,m,tau,q] = pradon_demultiple_my(d,dt,h,qmin,qmax,nq,flow,fhigh,mu,q_cutstart,q_cutend);
%PRADON_DEMULTIPLE: Multiple removal using parabolic Radon Transform.
%
%  [prim,m,tau,q] = pradon_demultiple(d,dt,h,qmin,qmax,nq,flow,fhigh,mu,q_cut);
%
%  IN   d:     data   (d(nt,nh))
%       dt:    sampling interval in secs
%       h:     offset in meters  (h(nh))
%       qmin:  min residual moveout at far offset (secs)
%       qmax:  max residual moveout at far offset (secs)
%       nq:    number of samples of the residual moveout axis
%       flow:  min freq to process (Hz)
%       flow:  max freq to process (Hz)
%       mu:    trade-off parameter for LS solution used
%              to retrieve the Radon panel
%       q_cut: keep contributions from q_cut to qmax,
%              this is to estimate the multiples that 
%              are removed from the primaries
%
%  OUT  prim: primaries obtained by removing from the data the
%             multiples modelled  with the Radon transform
%       m:    panel with the Radon transform  (m(nt,nq))
%       tau:  vertical axis for the Radon panel (tau(nt))
%       q:    horizontal axis for the Radon panel (q(nq))
%
%  Reference: Hampson, D., 1986, Inverse velocity stacking for multiple elimination,
%             Journal of the CSEG, vol 22, no 1., 44-55.
%
%  Example: see radon_demo.m 
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://www-geo.phys.ualberta.ca/saig/SeismicLab
%  Author: M.D.Sacchi
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%


  [nt,nx] = size(d);
  dq = (qmax-qmin)/(nq-1);
  q = qmin+dq*[0:1:nq-1];

  N = 2;                % parabolic trasnsform

% Transform from t-x to tau-q

  [m] = inverse_radon_freq(d,dt,h,q,N,flow,fhigh,mu,'ls');
  nq = length(q);
  
  tt=1;
  n_end=length(q_cutend(:,1)) ;
  vector_end=zeros(nt,1);
  for i=1:n_end-1
      k=(q_cutend(i+1,1)-q_cutend(i,1))/(q_cutend(i+1,2)-q_cutend(i,2));
      for j=0:dt:(q_cutend(i+1,2)-q_cutend(i,2)-dt)
          vector_end(tt)=floor((q_cutend(i,1)-qmin+k*j)/dq);
          tt=tt+1;
      end
  end
  tt=1;
  n_start=length(q_cutstart(:,1)) ;
  vector_start=zeros(nt,1);
  for i=1:n_start-1
      k=(q_cutstart(i+1,1)-q_cutstart(i,1))/(q_cutstart(i+1,2)-q_cutstart(i,2));
      for j=0:dt:(q_cutstart(i+1,2)-q_cutstart(i,2)-dt)
          vector_start(tt)=floor((q_cutstart(i,1)-qmin+k*j)/dq);
          tt=tt+1;
      end
  end
  
  mc = m;   
%   iq_cut_start=floor((q_cutstart-qmin)/dq)+1;
  for i=1:nt
      mc(i, vector_start(i):vector_end(i)) = 0;
  end
    
%   iq_cut_end = floor((q_cutend-qmin)/dq)+1;
%   iq_cut_start=floor((q_cutstart-qmin)/dq)+1;
%   mc = m;
%   mc(:, iq_cut_start:iq_cut_end) = 0;   % Keep multiples in the Radon panel

% Transform from tau-q to t-x

  [dm] = forward_radon_freq(mc,dt,h,q,N,flow,fhigh);

  prim = d-dm;

  tau = (0:1:nt-1)*dt;

  return;
