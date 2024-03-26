import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.constants import g


primitive_variables = ('Depth', 'Velocity')


def primitive_to_conservative(h, u):
    hu = h*u
    return h, hu


def pospart(x):
    return np.maximum(1.e-15,x)


def conservative_to_primitive(h, hu):
    # Check that h>=0 where it is not np.nan:
    assert np.nanmin(h)>=0
    # We should instead check that hu is zero everywhere that h is
    u = hu/pospart(h)
    return h, u


def exact_riemann_solution(q_l, q_r, force_waves=None, include_contact=False):
    """Return the exact solution to the Riemann problem with initial states q_l, q_r.
       The solution is given in terms of a list of states, a list of speeds (each of which
       may be a pair in case of a rarefaction fan), and a function reval(xi) that gives the
       solution at a point xi=x/t.

       The input and output vectors are the conserved quantities.
    """
    h_l, u_l = q_l
    h_r, u_r = q_r
    hu_l = h_l*u_l
    hu_r = h_r*u_r
    # if primitive_inputs:
    #     h_l, u_l = q_l
    #     h_r, u_r = q_r
    #     hu_l = h_l*u_l
    #     hu_r = h_r*u_r
    # else:
    #     h_l, u_l = conservative_to_primitive(*q_l)
    #     h_r, u_r = conservative_to_primitive(*q_r)
    #     hu_l = q_l[1]
    #     hu_r = q_r[1]

    # Compute left and right state sound speeds
    c_l = np.sqrt(g*h_l)
    c_r = np.sqrt(g*h_r)

    # Define the integral curves and hugoniot loci
    # Avoid warnings due to negative depths in fsolve calls
    integral_curve_1   = lambda h: u_l + 2*(np.sqrt(g*h_l) -
                                            np.sqrt(g*np.maximum(h,0)))
    integral_curve_2   = lambda h: u_r - 2*(np.sqrt(g*h_r) -
                                            np.sqrt(g*np.maximum(h,0)))
    hugoniot_locus_1 = lambda h: (h_l*u_l + (h-h_l)*(u_l -
                                  np.sqrt(g*h_l*(1 + (h-h_l)/h_l) * (1 + (h-h_l)/(2*h_l)))))/h
    hugoniot_locus_2 = lambda h: (h_r*u_r + (h-h_r)*(u_r +
                                  np.sqrt(g*h_r*(1 + (h-h_r)/h_r) * (1 + (h-h_r)/(2*h_r)))))/h

    # Check whether the 1-wave is a shock or rarefaction
    def phi_l(h):
        if (h>=h_l and force_waves!='raref') or force_waves=='shock':
            return hugoniot_locus_1(h)
        else:
            return integral_curve_1(h)

    # Check whether the 2-wave is a shock or rarefaction
    def phi_r(h):
        if (h>=h_r and force_waves!='raref') or force_waves=='shock':
            return hugoniot_locus_2(h)
        else:
            return integral_curve_2(h)

    ws = np.zeros(4)
    wave_types = ['', '']

    dry_velocity_l = u_l + 2*np.sqrt(g*h_l)
    dry_velocity_r = u_r - 2*np.sqrt(g*h_r)
    if dry_velocity_l < dry_velocity_r:
        # Dry middle state
        h_m = 0
        # This is a bit arbitrary:
        u_m = 0.5*(dry_velocity_l + dry_velocity_r)

        hu_m = u_m * h_m
        ws[0] = u_l - c_l
        ws[1] = dry_velocity_l
        ws[2] = dry_velocity_r
        ws[3] = u_r + c_r

    elif h_l == 0:
        # Dry left state; 2-rarefaction only
        h_m = 0
        u_m = dry_velocity_r
        hu_m = u_m * h_m
        ws[0] = 0
        ws[1] = 0
        ws[2] = dry_velocity_r
        ws[3] = u_r + c_r

    elif h_r == 0:
        # Dry right state; 1-rarefaction only
        h_m = 0
        u_m = dry_velocity_l
        hu_m = u_m * h_m
        ws[0] = u_l - c_l
        ws[1] = dry_velocity_l
        ws[2] = 0
        ws[3] = 0

    else:
        phi = lambda h: phi_l(h)-phi_r(h)

        # Compute middle state h, hu by finding curve intersection
        guess = (u_l-u_r+2.*np.sqrt(g)*(np.sqrt(h_l)+np.sqrt(h_r)))**2./16./g
        h_m, _, ier, msg = fsolve(phi, guess, full_output=True, xtol=1.e-14)
        # For strong rarefactions, sometimes fsolve needs help
        if ier!=1:
            h_m, _, ier, msg = fsolve(phi, guess,full_output=True,factor=0.1,xtol=1.e-10)
            # This should not happen:
            if ier!=1:
                print('Warning: fsolve did not converge.')
                print(msg)

        u_m = phi_l(h_m)
        hu_m = u_m * h_m

        # Find shock and rarefaction speeds
        if (h_m>h_l and force_waves!='raref') or force_waves=='shock':
            wave_types[0] = 'shock'
            ws[0] = (hu_l - hu_m) / (h_l - h_m)
            ws[1] = ws[0]
        else:
            wave_types[0] = 'raref'
            c_m = np.sqrt(g * h_m)
            ws[0] = u_l - c_l
            ws[1] = u_m - c_m

        if (h_m>h_r and force_waves!='raref') or force_waves=='shock':
            wave_types[1] = 'shock'
            ws[2] = (hu_r - hu_m) / (h_r - h_m)
            ws[3] = ws[2]
        else:
            wave_types[1] = 'raref'
            c_m = np.sqrt(g * h_m)
            ws[2] = u_m + c_m
            ws[3] = u_r + c_r

    # Find solution inside rarefaction fans (in primitive variables)
    def raref1(xi):
        RiemannInvariant = u_l + 2*np.sqrt(g*h_l)
        h = ((RiemannInvariant - xi)**2 / (9*g))
        u = (xi + np.sqrt(g*h))
        hu = h*u
        return h, hu

    def raref2(xi):
        RiemannInvariant = u_r - 2*np.sqrt(g*h_r)
        h = ((RiemannInvariant - xi)**2 / (9*g))
        u = (xi - np.sqrt(g*h))
        hu = h*u
        return h, hu

    q_m = np.squeeze(np.array((h_m, hu_m)))

    states = np.column_stack([q_l,q_m,q_r])
    speeds = [[], []]
    if wave_types[0] == 'shock':
        speeds[0] = ws[0]
    else:
        speeds[0] = (ws[0],ws[1])
    if wave_types[1] == 'shock':
        speeds[1] = ws[2]
    else:
        speeds[1] = (ws[2],ws[3])

    if include_contact:
        states = np.column_stack([q_l,q_m,q_m,q_r])
        um = q_m[1]/q_m[0]
        speeds = [speeds[0],um,speeds[1]]
        wave_types = [wave_types[0],'contact',wave_types[1]]
        
    def reval(xi):
        """
        Function that evaluates the Riemann solution for arbitrary xi = x/t.
        Sets the solution to nan in an over-turning rarefaction wave
        for illustration purposes of this non-physical solution.
        """
        rar1 = raref1(xi)
        rar2 = raref2(xi)
        # h_out = (xi<=ws[0])*h_l + \
        #     (xi>ws[0])*(xi<=ws[1])*rar1[0] + \
        #     (xi>ws[1])*(xi<=ws[0])*1e9 +  \
        #     (xi>ws[1])*(xi<=ws[2])*h_m +  \
        #     (xi>ws[2])*(xi<=ws[3])*rar2[0] +  \
        #     (xi>ws[3])*(xi<=ws[2])*1e9 +  \
        #     (xi>ws[3])*h_r
        # h_out[h_out>1e8] = np.nan
        # hu_out = (xi<=ws[0])*hu_l + \
        #     (xi>ws[0])*(xi<=ws[1])*rar1[1] + \
        #     (xi>ws[1])*(xi<=ws[0])*1e9 +  \
        #     (xi>ws[1])*(xi<=ws[2])*hu_m +  \
        #     (xi>ws[2])*(xi<=ws[3])*rar2[1] +  \
        #     (xi>ws[3])*(xi<=ws[2])*1e9 +  \
        #     (xi>ws[3])*hu_r
        # hu_out[hu_out>1e8] = np.nan
        h_out = (xi<=ws[0])*h_l + \
            (xi>ws[0])*(xi<=ws[1])*rar1[0] + \
            (xi>ws[1])*(xi<=ws[2])*h_m +  \
            (xi>ws[2])*(xi<=ws[3])*rar2[0] +  \
            (xi>ws[3])*h_r
        hu_out = (xi<=ws[0])*hu_l + \
            (xi>ws[0])*(xi<=ws[1])*rar1[1] + \
            (xi>ws[1])*(xi<=ws[2])*hu_m +  \
            (xi>ws[2])*(xi<=ws[3])*rar2[1] +  \
            (xi>ws[3])*hu_r
        return h_out, hu_out

    return states, speeds, reval, wave_types


def get_solution_func(hl=10.0, hr=1.0, ul=0.0, ur=0.0):
    states, speeds, reval, wave_types = exact_riemann_solution([hl, ul], [hr, ur])

    def get_exact_solution(xs, t):
        h, u = None, None
        if t == 0:
            h = hl*(xs<=0) + hr*(xs>0)
            u = ul*(xs<=0) + ur*(xs>0)
        else:
            q = np.array(reval(xs/t))
            h, u = conservative_to_primitive(q[0],q[1])
        return h, u
    
    return get_exact_solution


if __name__ == '__main__':
    pass
