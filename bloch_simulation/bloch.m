%
%	Bloch simulator.
%
%	[mx,my,mz] = bloch(b1,gr,tp,t1,t2,df,dp,mode,mx,my,mz)
%
%	Bloch simulation of rotations due to B1, gradient and
%	off-resonance, including relaxation effects.  At each time
%	point, the rotation matrix and decay matrix are calculated.
%	Simulation can simulate the steady-state if the sequence
%	is applied repeatedly, or the magnetization starting at m0.
%
%	INPUT:
%		b1 = (Mx1) RF pulse in G.  Can be complex.
%		gr = (Mx1,2,or 3) 1,2 or 3-dimensional gradient in G/cm.
%		tp = (Mx1) time duration of each b1 and gr point, in seconds,
%				or 1x1 time step if constant for all points
%				or monotonically INCREASING endtime of each
%				interval..
%		t1 = T1 relaxation time in seconds.
%		t2 = T2 relaxation time in seconds.
%		df = (Nx1) Array of off-resonance frequencies (Hz)
%		dp = (Px1,2,or 3) Array of spatial positions (cm).  
%			Width should match width of gr.
%		mode= Bitmask mode:
%			Bit 0:  0-Simulate from start or M0, 1-Steady State
%			Bit 1:  1-Record m at time points.  0-just end time.
%
%	(optional)
%		mx,my,mz (PxN) arrays of starting magnetization, where N
%			is the number of frequencies and P is the number
%			of spatial positions.
%
%	OUTPUT:
%		mx,my,mz = PxN arrays of the resulting magnetization
%				components at each position and frequency.
%
%	B. Hargreaves.	Nov 2003.
%






