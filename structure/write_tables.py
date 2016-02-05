#!/usr/bin/env python

import sys,os,math

#rdeb = 1.0	# nm
#ionic_strength=0.4 # M
#rrep = 0.5 # nm
#ratt = 0.7 # nm
rrep = float(sys.argv[1])
ratt = float(sys.argv[2])
#ionic_strength=float(sys.argv[1])
#ionic_strength_number = ionic_strength * 6.02e23/(1.e-3)
eps0 = 8.854e-12 
eps_r = 80.
kB = 8.31/6.02e23
T=300.

#rdeb = (eps_r*eps0*kB*T/(ionic_strength_number*(1.602e-19)**2))**0.5 * 1.e9

rlo = 0.	# nm
rhi = 11.0	# nm
eps = 80.
dr = 0.002	
#dr = 0.0001
nr = int(math.floor((rhi-rlo)/dr))+1

#for i in range(nr):
#	r = rlo+float(i)*dr
#	rx = r
#	if r == 0:
#		r = dr
#	expr = math.exp(-r/rdeb)
#	f = expr/(eps*r)
#	df = f * ( 1./r + 1./rdeb )
#	sys.stdout.write("%15.10e   %15.10e %15.10e   %15.10e %15.10e   %15.10e %15.10e\n" \
#			% (rx,f,df,f,df,f,df))

# write generic 6-9 (test ...)... seems to work
#for i in range(nr):
#	r = rlo+float(i)*dr
#	rx = r
#        if r < 4.e-2:
#            # avoid region where things get out of control
#            f = 0; df = 0; 
#            g = 0; dg = 0; 
#            h = 0; dh = 0; 
#        else:
#	    f = 1./r
#	    df = 1./r**2 # -f'
#            g= -1./r**6
#            dg = -6./r**7
#            h= 1./r**9
#            dh = 9./r**10
#	sys.stdout.write("%15.10e   %15.10e %15.10e   %15.10e %15.10e   %15.10e %15.10e\n" \
#			% (rx,f,df,g,dg,h,dh))

# write Weeks-Chandler-Anderson type repulsive and attractive parts...
for i in range(nr):
	r = rlo+float(i)*dr
	rx = r
        if r < 4.e-2:
            # avoid region where things get out of control
            f = 0; df = 0; 
            g = 0; dg = 0; 
            h = 0; dh = 0; 
        else:
            # elect fwiw
	    f = 1./r
	    df = 1./r**2 # -f''
            # attractive arm: 
            if r < ratt:
                g = -1.
                dg = 0.
            else:
                rattr = ratt/r
                rattr6 = rattr**6
                g = rattr6**2 - 2.*rattr6
                dg = 12.*rattr6**2/r - 12.*rattr6/r # = -g'
            if r > rrep:
                h = 0.
                dh = 0.
            else:
                rrepr = rrep/r
                rrepr6 = rrepr**6
                h = rrepr6**2 - 2.*rrepr6 + 1.
                dh = 12.*rrepr6**2/r - 12.*rrepr6/r # = -g'
	sys.stdout.write("%15.10e   %15.10e %15.10e   %15.10e %15.10e   %15.10e %15.10e\n" \
			% (rx,f,df,g,dg,h,dh))
