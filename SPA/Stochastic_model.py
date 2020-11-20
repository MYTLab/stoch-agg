# Author: Ace Shen, 2019-9-18
"""This source file contains all the ODEs (different closing orders)
    needed for generating numerical solutions used in the paper.
    All functions are constructed for the use of the scipy.integrate.odeint module:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html#scipy-integrate-odeint
    
    model(z, t, para):
        _ = parameters(*para)  # the name here has to be consistent with that assigned in the Mathematica module
        p1m0, ... = z  # list of variables to be solved

        # higher moments in terms of lower moments
        ... ...
        
        # set of ODEs
        ... ...
    
    The variables in z are expected values of p and m, < p^a * m^b >, denoted as pamb.
    e.g, p1m0 == <p>, p2m1 == <p**2 * m>
    
    If the needed stochastic model is not in this script,
    please update it from the result of the Mathematica module "moment_closure_generate_pythonscript.nb",
    following the same format as described above.
    """
from collections import namedtuple


parameters = namedtuple('Parameters', ['kn', 'kp', 'km', 'kf', 'mt', 'k2', 'nc', 'n2'])

# =============================
# ====== Deterministic ========
# =============================
# Note: need to specify nc & n2 in para


def deterministic(z, t, para):
    _ = parameters(*para)
    p1m0, p0m1 = z
    # ODEs
    dp1m0dt = _.kn * pow(p0m1, _.nc) + _.kf * (_.mt - p0m1 - (2 * -_.nc - 1) * p1m0) +\
              _.k2 * pow(p0m1, _.n2) * (_.mt - p0m1)
    dp0m1dt = -_.nc * _.kn * pow(p0m1, -_.nc) - 2 * _.kp * p0m1 * p1m0 + 2 * _.km * p1m0 +\
              _.kf * -_.nc * (-_.nc - 1) * p1m0 - _.n2 * _.k2 * pow(p0m1, _.n2) * (_.mt - p0m1)
    dzdt = [dp1m0dt, dp0m1dt]
    return dzdt

# =============================
# ====== nc = 2 ===============
# =============================


def stochastic_close2_nc2_n20(z, t, para):
    _ = parameters(*para)
    p1m0, p0m1, p2m0, p1m1, p0m2, = z
    # higher moments in terms of lower moments
    # p3m0 = pow(p1m0, -3) * pow(p0m1, 0) * pow(p2m0, 3) * pow(p1m1, 0) * pow(p0m2, 0)
    p2m1 = pow(p1m0, -2) * pow(p0m1, -1) * pow(p2m0, 1) * pow(p1m1, 2) * pow(p0m2, 0)
    p1m2 = pow(p1m0, -1) * pow(p0m1, -2) * pow(p2m0, 0) * pow(p1m1, 2) * pow(p0m2, 1)
    p0m3 = pow(p1m0, 0) * pow(p0m1, -3) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, 3)
    # p4m0 = pow(p1m0, -8) * pow(p0m1, 0) * pow(p2m0, 6) * pow(p1m1, 0) * pow(p0m2, 0)
    # p3m1 = pow(p1m0, -6) * pow(p0m1, -2) * pow(p2m0, 3) * pow(p1m1, 3) * pow(p0m2, 0)
    # p2m2 = pow(p1m0, -4) * pow(p0m1, -4) * pow(p2m0, 1) * pow(p1m1, 4) * pow(p0m2, 1)
    # p1m3 = pow(p1m0, -2) * pow(p0m1, -6) * pow(p2m0, 0) * pow(p1m1, 3) * pow(p0m2, 3)
    # p0m4 = pow(p1m0, 0) * pow(p0m1, -8) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, 6)
    # ODEs
    dp1m0dt = _.kn*p0m2 - _.kf*p0m1 - (2*2-1)*_.kf*p1m0 + _.kf*_.mt + _.k2*_.mt*1 - _.k2*p0m1
    dp0m1dt = (-2)*_.kn*p0m2 - 2*_.kp*p1m1 + 2*_.km*p1m0 + 2*(2-1)*_.kf*p1m0 - 0*_.k2*_.mt*1 + 0*_.k2*p0m1
    dp2m0dt = 2*(_.kn*p1m2 - _.kf*p1m1 - (2*2-1)*_.kf*p2m0 + _.mt*_.kf*p1m0 + _.k2*_.mt*p1m0 - _.k2*p1m1) + 1*(_.kn*p0m2 - _.kf*p0m1 - (2*2-1)*_.kf*p1m0 + _.mt*_.kf*1 + _.k2*_.mt*1 - _.k2*p0m1)
    dp1m1dt = ((-2)*_.kn*p1m2 - 2*_.kp*p2m1 + 2*_.km*p2m0 + _.kf*2*(2-1)*p2m0 - 0*_.k2*_.mt*p1m0 + 0*_.k2*p1m1) + (_.kn*p0m3 - _.kf*p0m2 - (2*2-1)*_.kf*p1m1 + _.mt*_.kf*p0m1 + _.k2*_.mt*p0m1 - _.k2*p0m2) + ((-2)*_.kn * p0m2 - 0*_.k2*_.mt*1 + 0*_.k2*p0m1)
    dp0m2dt = 2*((-2)*_.kn*p0m3 - 2*_.kp*p1m2 + 2*_.km*p1m1 + _.kf*2*(2-1)*p1m1 - 0*_.k2*_.mt*p0m1 + 0*_.k2*p0m2) + 1*(pow( 2, 2) * _.kn*p0m2 + 2*_.kp*p1m1 + 2*_.km*p1m0 + 2*(2-1)*(2*2-1)*_.kf*p1m0/3. + pow(0, 2) *_.k2*_.mt*1 - pow(0, 2)*_.k2*p0m1)
    dzdt = [dp1m0dt, dp0m1dt, dp2m0dt, dp1m1dt, dp0m2dt]
    return dzdt


def stochastic_close2_nc2_n22(z, t, para):
    _ = parameters(*para)
    p1m0, p0m1, p2m0, p1m1, p0m2, = z
    # higher moments in terms of lower moments
    # p3m0 = pow(p1m0, -3) * pow(p0m1, 0) * pow(p2m0, 3) * pow(p1m1, 0) * pow(p0m2, 0)
    p2m1 = pow(p1m0, -2) * pow(p0m1, -1) * pow(p2m0, 1) * pow(p1m1, 2) * pow(p0m2, 0)
    p1m2 = pow(p1m0, -1) * pow(p0m1, -2) * pow(p2m0, 0) * pow(p1m1, 2) * pow(p0m2, 1)
    p0m3 = pow(p1m0, 0) * pow(p0m1, -3) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, 3)
    # p4m0 = pow(p1m0, -8) * pow(p0m1, 0) * pow(p2m0, 6) * pow(p1m1, 0) * pow(p0m2, 0)
    # p3m1 = pow(p1m0, -6) * pow(p0m1, -2) * pow(p2m0, 3) * pow(p1m1, 3) * pow(p0m2, 0)
    # p2m2 = pow(p1m0, -4) * pow(p0m1, -4) * pow(p2m0, 1) * pow(p1m1, 4) * pow(p0m2, 1)
    p1m3 = pow(p1m0, -2) * pow(p0m1, -6) * pow(p2m0, 0) * pow(p1m1, 3) * pow(p0m2, 3)
    p0m4 = pow(p1m0, 0) * pow(p0m1, -8) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, 6)
    # ODEs
    dp1m0dt = _.kn*p0m2 - _.kf*p0m1 - (2*2-1)*_.kf*p1m0 + _.kf*_.mt + _.k2*_.mt*p0m2 - _.k2*p0m3
    dp0m1dt = (-2)*_.kn*p0m2 - 2*_.kp*p1m1 + 2*_.km*p1m0 + 2*(2-1)*_.kf*p1m0 - 2*_.k2*_.mt*p0m2 + 2*_.k2*p0m3
    dp2m0dt = 2*(_.kn*p1m2 - _.kf*p1m1 - (2*2-1)*_.kf*p2m0 + _.mt*_.kf*p1m0 + _.k2*_.mt*p1m2 - _.k2*p1m3) +\
              1*(_.kn*p0m2 - _.kf*p0m1 - (2*2-1)*_.kf*p1m0 + _.mt*_.kf*1 + _.k2*_.mt*p0m2 - _.k2*p0m3)
    dp1m1dt = ((-2)*_.kn*p1m2 - 2*_.kp*p2m1 + 2*_.km*p2m0 + _.kf*2*(2-1)*p2m0 - 2*_.k2*_.mt*p1m2 + 2*_.k2*p1m3) +\
              (_.kn*p0m3 - _.kf*p0m2 - (2*2-1)*_.kf*p1m1 + _.mt*_.kf*p0m1 + _.k2*_.mt*p0m3 - _.k2*p0m4) + ((-2)*_.kn * p0m2 - 2*_.k2*_.mt*p0m2 + 2*_.k2*p0m3)
    dp0m2dt = 2*((-2)*_.kn*p0m3 - 2*_.kp*p1m2 + 2*_.km*p1m1 + _.kf*2*(2-1)*p1m1 - 2*_.k2*_.mt*p0m3 + 2*_.k2*p0m4) +\
              1*(pow( 2, 2) * _.kn*p0m2 + 2*_.kp*p1m1 + 2*_.km*p1m0 + 2*(2-1)*(2*2-1)*_.kf*p1m0/3. + pow(2, 2) *_.k2*_.mt*p0m2 - pow(2, 2)*_.k2*p0m3)
    dzdt = [dp1m0dt, dp0m1dt, dp2m0dt, dp1m1dt, dp0m2dt]
    return dzdt

# Close at 3rd order. Use z0_close3 as initial condition


def stochastic_close3_nc2_n20(z, t, para):
    _ = parameters(*para)
    p1m0, p0m1, p2m0, p1m1, p0m2, p3m0, p2m1, p1m2, p0m3, = z
    # higher moments in terms of lower moments
    # p4m0 = pow(p1m0, 4) * pow(p0m1, 0) * pow(p2m0, -6) * pow(p1m1, 0) * pow(p0m2, 0) * pow(p3m0, 4) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 0)
    p3m1 = pow(p1m0, 3) * pow(p0m1, 1) * pow(p2m0, -3) * pow(p1m1, -3) * pow(p0m2, 0) * pow(p3m0, 1) * pow(p2m1, 3) * pow(p1m2, 0) * pow(p0m3, 0)
    p2m2 = pow(p1m0, 2) * pow(p0m1, 2) * pow(p2m0, -1) * pow(p1m1, -4) * pow(p0m2, -1) * pow(p3m0, 0) * pow(p2m1, 2) * pow(p1m2, 2) * pow(p0m3, 0)
    p1m3 = pow(p1m0, 1) * pow(p0m1, 3) * pow(p2m0, 0) * pow(p1m1, -3) * pow(p0m2, -3) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 3) * pow(p0m3, 1)
    p0m4 = pow(p1m0, 0) * pow(p0m1, 4) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, -6) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 4)
    # p5m0 = pow(p1m0, 15) * pow(p0m1, 0) * pow(p2m0, -20) * pow(p1m1, 0) * pow(p0m2, 0) * pow(p3m0, 10) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 0)
    # p4m1 = pow(p1m0, 12) * pow(p0m1, 3) * pow(p2m0, -12) * pow(p1m1, -8) * pow(p0m2, 0) * pow(p3m0, 4) * pow(p2m1, 6) * pow(p1m2, 0) * pow(p0m3, 0)
    # p3m2 = pow(p1m0, 9) * pow(p0m1, 6) * pow(p2m0, -6) * pow(p1m1, -12) * pow(p0m2, -2) * pow(p3m0, 1) * pow(p2m1, 6) * pow(p1m2, 3) * pow(p0m3, 0)
    # p2m3 = pow(p1m0, 6) * pow(p0m1, 9) * pow(p2m0, -2) * pow(p1m1, -12) * pow(p0m2, -6) * pow(p3m0, 0) * pow(p2m1, 3) * pow(p1m2, 6) * pow(p0m3, 1)
    # p1m4 = pow(p1m0, 3) * pow(p0m1, 12) * pow(p2m0, 0) * pow(p1m1, -8) * pow(p0m2, -12) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 6) * pow(p0m3, 4)
    # p0m5 = pow(p1m0, 0) * pow(p0m1, 15) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, -20) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 10)
    # ODEs
    dp1m0dt = _.kn*p0m2 - _.kf*p0m1 - (2*2-1)*_.kf*p1m0 + _.kf*_.mt + _.k2*_.mt*1 - _.k2*p0m1
    dp0m1dt = (-2)*_.kn*p0m2 - 2*_.kp*p1m1 + 2*_.km*p1m0 + 2*(2-1)*_.kf*p1m0 - 0*_.k2*_.mt*1 + 0*_.k2*p0m1
    dp2m0dt = 2*(_.kn*p1m2 - _.kf*p1m1 - (2*2-1)*_.kf*p2m0 + _.mt*_.kf*p1m0 + _.k2*_.mt*p1m0 - _.k2*p1m1) +\
              1*(_.kn*p0m2 - _.kf*p0m1 - (2*2-1)*_.kf*p1m0 + _.mt*_.kf*1 + _.k2*_.mt*1 - _.k2*p0m1)
    dp1m1dt = ((-2)*_.kn*p1m2 - 2*_.kp*p2m1 + 2*_.km*p2m0 + _.kf*2*(2-1)*p2m0 - 0*_.k2*_.mt*p1m0 + 0*_.k2*p1m1) +\
              (_.kn*p0m3 - _.kf*p0m2 - (2*2-1)*_.kf*p1m1 + _.mt*_.kf*p0m1 + _.k2*_.mt*p0m1 - _.k2*p0m2) + ((-2)*_.kn * p0m2 - 0*_.k2*_.mt*1 + 0*_.k2*p0m1)
    dp0m2dt = 2*((-2)*_.kn*p0m3 - 2*_.kp*p1m2 + 2*_.km*p1m1 + _.kf*2*(2-1)*p1m1 - 0*_.k2*_.mt*p0m1 + 0*_.k2*p0m2) +\
              1*(pow( 2, 2) * _.kn*p0m2 + 2*_.kp*p1m1 + 2*_.km*p1m0 + 2*(2-1)*(2*2-1)*_.kf*p1m0/3. + pow(0, 2) *_.k2*_.mt*1 - pow(0, 2)*_.k2*p0m1)
    dp3m0dt = 3*(_.kn*p2m2 - _.kf*p2m1 - (2*2-1)*_.kf*p3m0 + _.mt*_.kf*p2m0 + _.k2*_.mt*p2m0 - _.k2*p2m1) +\
              3*(_.kn*p1m2 - _.kf*p1m1 - (2*2-1)*_.kf*p2m0 + _.mt*_.kf*p1m0 + _.k2*_.mt*p1m0 - _.k2*p1m1)
    dp2m1dt = 1*((-2)*_.kn*p2m2 - 2*_.kp*p3m1 + 2*_.km*p3m0 + _.kf*2*(2-1)*p3m0 - 0*_.k2*_.mt*p2m0 + 0*_.k2*p2m1) +\
              2*(_.kn*p1m3 - _.kf*p1m2 - (2*2-1)*_.kf*p2m1 + _.mt*_.kf*p1m1 + _.k2*_.mt*p1m1 - _.k2*p1m2) +\
              1*(_.kn*p0m3 - _.kf*p0m2 - (2*2-1)*_.kf*p1m1 + _.mt*_.kf*p0m1 + _.k2*_.mt*p0m1 - _.k2*p0m2) + 2*1*((-2)*_.kn * p1m2 - 0*_.k2*_.mt*p1m0 + 0*_.k2*p1m1)
    dp1m2dt = 2*((-2)*_.kn*p1m3 - 2*_.kp*p2m2 + 2*_.km*p2m1 + _.kf*2*(2-1)*p2m1 - 0*_.k2*_.mt*p1m1 + 0*_.k2*p1m2) +\
              1*(_.kn*p0m4 - _.kf*p0m3 - (2*2-1)*_.kf*p1m2 + _.mt*_.kf*p0m2 + _.k2*_.mt*p0m2 - _.k2*p0m3) +\
              1*(pow( 2, 2) * _.kn*p1m2 + 2*_.kp*p2m1 + 2*_.km*p2m0 + 2*(2-1)*(2*2-1)*_.kf*p2m0/3. +\
                 pow(0, 2) *_.k2*_.mt*p1m0 - pow(0, 2)*_.k2*p1m1) + 1*2*((-2)*_.kn * p0m3 - 0*_.k2*_.mt*p0m1 + 0*_.k2*p0m2)
    dp0m3dt = 3*((-2)*_.kn*p0m4 - 2*_.kp*p1m3 + 2*_.km*p1m2 + _.kf*2*(2-1)*p1m2 - 0*_.k2*_.mt*p0m2 + 0*_.k2*p0m3) +\
              3*(pow( 2, 2) * _.kn*p0m3 + 2*_.kp*p1m2 + 2*_.km*p1m1 + 2*(2-1)*(2*2-1)*_.kf*p1m1/3. + pow(0, 2) *_.k2*_.mt*p0m1 - pow(0, 2)*_.k2*p0m2)
    dzdt = [dp1m0dt, dp0m1dt, dp2m0dt, dp1m1dt, dp0m2dt, dp3m0dt, dp2m1dt, dp1m2dt, dp0m3dt]
    return dzdt


def stochastic_close3_nc2_n22(z, t, para):
    _ = parameters(*para)
    p1m0, p0m1, p2m0, p1m1, p0m2, p3m0, p2m1, p1m2, p0m3, = z
    # higher moments in terms of lower moments
    # p4m0 = pow(p1m0, 4) * pow(p0m1, 0) * pow(p2m0, -6) * pow(p1m1, 0) * pow(p0m2, 0) * pow(p3m0, 4) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 0)
    p3m1 = pow(p1m0, 3) * pow(p0m1, 1) * pow(p2m0, -3) * pow(p1m1, -3) * pow(p0m2, 0) * pow(p3m0, 1) * pow(p2m1, 3) * pow(p1m2, 0) * pow(p0m3, 0)
    p2m2 = pow(p1m0, 2) * pow(p0m1, 2) * pow(p2m0, -1) * pow(p1m1, -4) * pow(p0m2, -1) * pow(p3m0, 0) * pow(p2m1, 2) * pow(p1m2, 2) * pow(p0m3, 0)
    p1m3 = pow(p1m0, 1) * pow(p0m1, 3) * pow(p2m0, 0) * pow(p1m1, -3) * pow(p0m2, -3) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 3) * pow(p0m3, 1)
    p0m4 = pow(p1m0, 0) * pow(p0m1, 4) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, -6) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 4)
    # p5m0 = pow(p1m0, 15) * pow(p0m1, 0) * pow(p2m0, -20) * pow(p1m1, 0) * pow(p0m2, 0) * pow(p3m0, 10) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 0)
    # p4m1 = pow(p1m0, 12) * pow(p0m1, 3) * pow(p2m0, -12) * pow(p1m1, -8) * pow(p0m2, 0) * pow(p3m0, 4) * pow(p2m1, 6) * pow(p1m2, 0) * pow(p0m3, 0)
    # p3m2 = pow(p1m0, 9) * pow(p0m1, 6) * pow(p2m0, -6) * pow(p1m1, -12) * pow(p0m2, -2) * pow(p3m0, 1) * pow(p2m1, 6) * pow(p1m2, 3) * pow(p0m3, 0)
    p2m3 = pow(p1m0, 6) * pow(p0m1, 9) * pow(p2m0, -2) * pow(p1m1, -12) * pow(p0m2, -6) * pow(p3m0, 0) * pow(p2m1, 3) * pow(p1m2, 6) * pow(p0m3, 1)
    p1m4 = pow(p1m0, 3) * pow(p0m1, 12) * pow(p2m0, 0) * pow(p1m1, -8) * pow(p0m2, -12) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 6) * pow(p0m3, 4)
    p0m5 = pow(p1m0, 0) * pow(p0m1, 15) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, -20) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 10)
    # ODEs
    dp1m0dt = _.kn*p0m2 - _.kf*p0m1 - (2*2-1)*_.kf*p1m0 + _.kf*_.mt + _.k2*_.mt*p0m2 - _.k2*p0m3
    dp0m1dt = (-2)*_.kn*p0m2 - 2*_.kp*p1m1 + 2*_.km*p1m0 + 2*(2-1)*_.kf*p1m0 - 2*_.k2*_.mt*p0m2 + 2*_.k2*p0m3
    dp2m0dt = 2*(_.kn*p1m2 - _.kf*p1m1 - (2*2-1)*_.kf*p2m0 + _.mt*_.kf*p1m0 + _.k2*_.mt*p1m2 - _.k2*p1m3) +\
              1*(_.kn*p0m2 - _.kf*p0m1 - (2*2-1)*_.kf*p1m0 + _.mt*_.kf*1 + _.k2*_.mt*p0m2 - _.k2*p0m3)
    dp1m1dt = ((-2)*_.kn*p1m2 - 2*_.kp*p2m1 + 2*_.km*p2m0 + _.kf*2*(2-1)*p2m0 - 2*_.k2*_.mt*p1m2 + 2*_.k2*p1m3) +\
              (_.kn*p0m3 - _.kf*p0m2 - (2*2-1)*_.kf*p1m1 + _.mt*_.kf*p0m1 + _.k2*_.mt*p0m3 - _.k2*p0m4) + ((-2)*_.kn * p0m2 - 2*_.k2*_.mt*p0m2 + 2*_.k2*p0m3)
    dp0m2dt = 2*((-2)*_.kn*p0m3 - 2*_.kp*p1m2 + 2*_.km*p1m1 + _.kf*2*(2-1)*p1m1 - 2*_.k2*_.mt*p0m3 + 2*_.k2*p0m4) +\
              1*(pow( 2, 2) * _.kn*p0m2 + 2*_.kp*p1m1 + 2*_.km*p1m0 + 2*(2-1)*(2*2-1)*_.kf*p1m0/3. + pow(2, 2) *_.k2*_.mt*p0m2 - pow(2, 2)*_.k2*p0m3)
    dp3m0dt = 3*(_.kn*p2m2 - _.kf*p2m1 - (2*2-1)*_.kf*p3m0 + _.mt*_.kf*p2m0 + _.k2*_.mt*p2m2 - _.k2*p2m3) +\
              3*(_.kn*p1m2 - _.kf*p1m1 - (2*2-1)*_.kf*p2m0 + _.mt*_.kf*p1m0 + _.k2*_.mt*p1m2 - _.k2*p1m3)
    dp2m1dt = 1*((-2)*_.kn*p2m2 - 2*_.kp*p3m1 + 2*_.km*p3m0 + _.kf*2*(2-1)*p3m0 - 2*_.k2*_.mt*p2m2 + 2*_.k2*p2m3) +\
              2*(_.kn*p1m3 - _.kf*p1m2 - (2*2-1)*_.kf*p2m1 + _.mt*_.kf*p1m1 + _.k2*_.mt*p1m3 - _.k2*p1m4) +\
              1*(_.kn*p0m3 - _.kf*p0m2 - (2*2-1)*_.kf*p1m1 + _.mt*_.kf*p0m1 + _.k2*_.mt*p0m3 - _.k2*p0m4) + 2*1*((-2)*_.kn * p1m2 - 2*_.k2*_.mt*p1m2 + 2*_.k2*p1m3)
    dp1m2dt = 2*((-2)*_.kn*p1m3 - 2*_.kp*p2m2 + 2*_.km*p2m1 + _.kf*2*(2-1)*p2m1 - 2*_.k2*_.mt*p1m3 + 2*_.k2*p1m4) +\
              1*(_.kn*p0m4 - _.kf*p0m3 - (2*2-1)*_.kf*p1m2 + _.mt*_.kf*p0m2 + _.k2*_.mt*p0m4 - _.k2*p0m5) +\
              1*(pow( 2, 2) * _.kn*p1m2 + 2*_.kp*p2m1 + 2*_.km*p2m0 + 2*(2-1)*(2*2-1)*_.kf*p2m0/3. +\
                 pow(2, 2) *_.k2*_.mt*p1m2 - pow(2, 2)*_.k2*p1m3) + 1*2*((-2)*_.kn * p0m3 - 2*_.k2*_.mt*p0m3 + 2*_.k2*p0m4)
    dp0m3dt = 3*((-2)*_.kn*p0m4 - 2*_.kp*p1m3 + 2*_.km*p1m2 + _.kf*2*(2-1)*p1m2 - 2*_.k2*_.mt*p0m4 + 2*_.k2*p0m5) +\
              3*(pow( 2, 2) * _.kn*p0m3 + 2*_.kp*p1m2 + 2*_.km*p1m1 + 2*(2-1)*(2*2-1)*_.kf*p1m1/3. + pow(2, 2) *_.k2*_.mt*p0m3 - pow(2, 2)*_.k2*p0m4)
    dzdt = [dp1m0dt, dp0m1dt, dp2m0dt, dp1m1dt, dp0m2dt, dp3m0dt, dp2m1dt, dp1m2dt, dp0m3dt]
    return dzdt

# =============================
# ====== nc = 3 ===============
# =============================
# Close at 2nd order. Use z0_close2 as initial condition:


def stochastic_close2_nc3_n20(z, t, para):
    _ = parameters(*para)
    p1m0, p0m1, p2m0, p1m1, p0m2, = z
    # higher moments in terms of lower moments
    # p3m0 = pow(p1m0, -3) * pow(p0m1, 0) * pow(p2m0, 3) * pow(p1m1, 0) * pow(p0m2, 0)
    p2m1 = pow(p1m0, -2) * pow(p0m1, -1) * pow(p2m0, 1) * pow(p1m1, 2) * pow(p0m2, 0)
    p1m2 = pow(p1m0, -1) * pow(p0m1, -2) * pow(p2m0, 0) * pow(p1m1, 2) * pow(p0m2, 1)
    p0m3 = pow(p1m0, 0) * pow(p0m1, -3) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, 3)
    # p4m0 = pow(p1m0, -8) * pow(p0m1, 0) * pow(p2m0, 6) * pow(p1m1, 0) * pow(p0m2, 0)
    # p3m1 = pow(p1m0, -6) * pow(p0m1, -2) * pow(p2m0, 3) * pow(p1m1, 3) * pow(p0m2, 0)
    # p2m2 = pow(p1m0, -4) * pow(p0m1, -4) * pow(p2m0, 1) * pow(p1m1, 4) * pow(p0m2, 1)
    p1m3 = pow(p1m0, -2) * pow(p0m1, -6) * pow(p2m0, 0) * pow(p1m1, 3) * pow(p0m2, 3)
    p0m4 = pow(p1m0, 0) * pow(p0m1, -8) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, 6)
    # p5m0 = pow(p1m0, -15) * pow(p0m1, 0) * pow(p2m0, 10) * pow(p1m1, 0) * pow(p0m2, 0)
    # p4m1 = pow(p1m0, -12) * pow(p0m1, -3) * pow(p2m0, 6) * pow(p1m1, 4) * pow(p0m2, 0)
    # p3m2 = pow(p1m0, -9) * pow(p0m1, -6) * pow(p2m0, 3) * pow(p1m1, 6) * pow(p0m2, 1)
    # p2m3 = pow(p1m0, -6) * pow(p0m1, -9) * pow(p2m0, 1) * pow(p1m1, 6) * pow(p0m2, 3)
    # p1m4 = pow(p1m0, -3) * pow(p0m1, -12) * pow(p2m0, 0) * pow(p1m1, 4) * pow(p0m2, 6)
    # p0m5 = pow(p1m0, 0) * pow(p0m1, -15) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, 10)
    # ODEs
    dp1m0dt = _.kn*p0m3 - _.kf*p0m1 - (2*3-1)*_.kf*p1m0 + _.kf*_.mt + _.k2*_.mt*1 - _.k2*p0m1
    dp0m1dt = (-3)*_.kn*p0m3 - 2*_.kp*p1m1 + 2*_.km*p1m0 + 3*(3-1)*_.kf*p1m0 - 0*_.k2*_.mt*1 + 0*_.k2*p0m1
    dp2m0dt = 2*(_.kn*p1m3 - _.kf*p1m1 - (2*3-1)*_.kf*p2m0 + _.mt*_.kf*p1m0 + _.k2*_.mt*p1m0 - _.k2*p1m1) +\
              1*(_.kn*p0m3 - _.kf*p0m1 - (2*3-1)*_.kf*p1m0 + _.mt*_.kf*1 + _.k2*_.mt*1 - _.k2*p0m1)
    dp1m1dt = ((-3)*_.kn*p1m3 - 2*_.kp*p2m1 + 2*_.km*p2m0 + _.kf*3*(3-1)*p2m0 - 0*_.k2*_.mt*p1m0 + 0*_.k2*p1m1) +\
              (_.kn*p0m4 - _.kf*p0m2 - (2*3-1)*_.kf*p1m1 + _.mt*_.kf*p0m1 + _.k2*_.mt*p0m1 - _.k2*p0m2) + ((-3)*_.kn * p0m3 - 0*_.k2*_.mt*1 + 0*_.k2*p0m1)
    dp0m2dt = 2*((-3)*_.kn*p0m4 - 2*_.kp*p1m2 + 2*_.km*p1m1 + _.kf*3*(3-1)*p1m1 - 0*_.k2*_.mt*p0m1 + 0*_.k2*p0m2) +\
              1*(pow( 3, 2) * _.kn*p0m3 + 2*_.kp*p1m1 + 2*_.km*p1m0 + 3*(3-1)*(2*3-1)*_.kf*p1m0/3. + pow(0, 2) *_.k2*_.mt*1 - pow(0, 2)*_.k2*p0m1)
    dzdt = [dp1m0dt, dp0m1dt, dp2m0dt, dp1m1dt, dp0m2dt]
    return dzdt


def stochastic_close2_nc3_n22(z, t, para):
    _ = parameters(*para)
    p1m0, p0m1, p2m0, p1m1, p0m2, = z
    # higher moments in terms of lower moments
    # p3m0 = pow(p1m0, -3) * pow(p0m1, 0) * pow(p2m0, 3) * pow(p1m1, 0) * pow(p0m2, 0)
    p2m1 = pow(p1m0, -2) * pow(p0m1, -1) * pow(p2m0, 1) * pow(p1m1, 2) * pow(p0m2, 0)
    p1m2 = pow(p1m0, -1) * pow(p0m1, -2) * pow(p2m0, 0) * pow(p1m1, 2) * pow(p0m2, 1)
    p0m3 = pow(p1m0, 0) * pow(p0m1, -3) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, 3)
    # p4m0 = pow(p1m0, -8) * pow(p0m1, 0) * pow(p2m0, 6) * pow(p1m1, 0) * pow(p0m2, 0)
    # p3m1 = pow(p1m0, -6) * pow(p0m1, -2) * pow(p2m0, 3) * pow(p1m1, 3) * pow(p0m2, 0)
    # p2m2 = pow(p1m0, -4) * pow(p0m1, -4) * pow(p2m0, 1) * pow(p1m1, 4) * pow(p0m2, 1)
    p1m3 = pow(p1m0, -2) * pow(p0m1, -6) * pow(p2m0, 0) * pow(p1m1, 3) * pow(p0m2, 3)
    p0m4 = pow(p1m0, 0) * pow(p0m1, -8) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, 6)
    # p5m0 = pow(p1m0, -15) * pow(p0m1, 0) * pow(p2m0, 10) * pow(p1m1, 0) * pow(p0m2, 0)
    # p4m1 = pow(p1m0, -12) * pow(p0m1, -3) * pow(p2m0, 6) * pow(p1m1, 4) * pow(p0m2, 0)
    # p3m2 = pow(p1m0, -9) * pow(p0m1, -6) * pow(p2m0, 3) * pow(p1m1, 6) * pow(p0m2, 1)
    # p2m3 = pow(p1m0, -6) * pow(p0m1, -9) * pow(p2m0, 1) * pow(p1m1, 6) * pow(p0m2, 3)
    # p1m4 = pow(p1m0, -3) * pow(p0m1, -12) * pow(p2m0, 0) * pow(p1m1, 4) * pow(p0m2, 6)
    # p0m5 = pow(p1m0, 0) * pow(p0m1, -15) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, 10)
    # ODEs
    dp1m0dt = _.kn*p0m3 - _.kf*p0m1 - (2*3-1)*_.kf*p1m0 + _.kf*_.mt + _.k2*_.mt*p0m2 - _.k2*p0m3
    dp0m1dt = (-3)*_.kn*p0m3 - 2*_.kp*p1m1 + 2*_.km*p1m0 + 3*(3-1)*_.kf*p1m0 - 2*_.k2*_.mt*p0m2 + 2*_.k2*p0m3
    dp2m0dt = 2*(_.kn*p1m3 - _.kf*p1m1 - (2*3-1)*_.kf*p2m0 + _.mt*_.kf*p1m0 + _.k2*_.mt*p1m2 - _.k2*p1m3) +\
              1*(_.kn*p0m3 - _.kf*p0m1 - (2*3-1)*_.kf*p1m0 + _.mt*_.kf*1 + _.k2*_.mt*p0m2 - _.k2*p0m3)
    dp1m1dt = ((-3)*_.kn*p1m3 - 2*_.kp*p2m1 + 2*_.km*p2m0 + _.kf*3*(3-1)*p2m0 - 2*_.k2*_.mt*p1m2 + 2*_.k2*p1m3) +\
              (_.kn*p0m4 - _.kf*p0m2 - (2*3-1)*_.kf*p1m1 + _.mt*_.kf*p0m1 + _.k2*_.mt*p0m3 - _.k2*p0m4) + ((-3)*_.kn * p0m3 - 2*_.k2*_.mt*p0m2 + 2*_.k2*p0m3)
    dp0m2dt = 2*((-3)*_.kn*p0m4 - 2*_.kp*p1m2 + 2*_.km*p1m1 + _.kf*3*(3-1)*p1m1 - 2*_.k2*_.mt*p0m3 + 2*_.k2*p0m4) +\
              1*(pow( 3, 2) * _.kn*p0m3 + 2*_.kp*p1m1 + 2*_.km*p1m0 + 3*(3-1)*(2*3-1)*_.kf*p1m0/3. + pow(2, 2) *_.k2*_.mt*p0m2 - pow(2, 2)*_.k2*p0m3)
    dzdt = [dp1m0dt, dp0m1dt, dp2m0dt, dp1m1dt, dp0m2dt]
    return dzdt


# Close at 3rd order. Use z0_close3 as initial condition:


def stochastic_close3_nc3_n20(z, t, para):
    _ = parameters(*para)
    p1m0, p0m1, p2m0, p1m1, p0m2, p3m0, p2m1, p1m2, p0m3, = z
    # higher moments in terms of lower moments
    # p4m0 = pow(p1m0, 4) * pow(p0m1, 0) * pow(p2m0, -6) * pow(p1m1, 0) * pow(p0m2, 0) * pow(p3m0, 4) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 0)
    p3m1 = pow(p1m0, 3) * pow(p0m1, 1) * pow(p2m0, -3) * pow(p1m1, -3) * pow(p0m2, 0) * pow(p3m0, 1) * pow(p2m1, 3) * pow(p1m2, 0) * pow(p0m3, 0)
    p2m2 = pow(p1m0, 2) * pow(p0m1, 2) * pow(p2m0, -1) * pow(p1m1, -4) * pow(p0m2, -1) * pow(p3m0, 0) * pow(p2m1, 2) * pow(p1m2, 2) * pow(p0m3, 0)
    p1m3 = pow(p1m0, 1) * pow(p0m1, 3) * pow(p2m0, 0) * pow(p1m1, -3) * pow(p0m2, -3) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 3) * pow(p0m3, 1)
    p0m4 = pow(p1m0, 0) * pow(p0m1, 4) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, -6) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 4)
    # p5m0 = pow(p1m0, 15) * pow(p0m1, 0) * pow(p2m0, -20) * pow(p1m1, 0) * pow(p0m2, 0) * pow(p3m0, 10) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 0)
    # p4m1 = pow(p1m0, 12) * pow(p0m1, 3) * pow(p2m0, -12) * pow(p1m1, -8) * pow(p0m2, 0) * pow(p3m0, 4) * pow(p2m1, 6) * pow(p1m2, 0) * pow(p0m3, 0)
    # p3m2 = pow(p1m0, 9) * pow(p0m1, 6) * pow(p2m0, -6) * pow(p1m1, -12) * pow(p0m2, -2) * pow(p3m0, 1) * pow(p2m1, 6) * pow(p1m2, 3) * pow(p0m3, 0)
    p2m3 = pow(p1m0, 6) * pow(p0m1, 9) * pow(p2m0, -2) * pow(p1m1, -12) * pow(p0m2, -6) * pow(p3m0, 0) * pow(p2m1, 3) * pow(p1m2, 6) * pow(p0m3, 1)
    p1m4 = pow(p1m0, 3) * pow(p0m1, 12) * pow(p2m0, 0) * pow(p1m1, -8) * pow(p0m2, -12) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 6) * pow(p0m3, 4)
    p0m5 = pow(p1m0, 0) * pow(p0m1, 15) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, -20) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 10)
    # p6m0 = pow(p1m0, 36) * pow(p0m1, 0) * pow(p2m0, -45) * pow(p1m1, 0) * pow(p0m2, 0) * pow(p3m0, 20) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 0)
    # p5m1 = pow(p1m0, 30) * pow(p0m1, 6) * pow(p2m0, -30) * pow(p1m1, -15) * pow(p0m2, 0) * pow(p3m0, 10) * pow(p2m1, 10) * pow(p1m2, 0) * pow(p0m3, 0)
    # p4m2 = pow(p1m0, 24) * pow(p0m1, 12) * pow(p2m0, -18) * pow(p1m1, -24) * pow(p0m2, -3) * pow(p3m0, 4) * pow(p2m1, 12) * pow(p1m2, 4) * pow(p0m3, 0)
    # p3m3 = pow(p1m0, 18) * pow(p0m1, 18) * pow(p2m0, -9) * pow(p1m1, -27) * pow(p0m2, -9) * pow(p3m0, 1) * pow(p2m1, 9) * pow(p1m2, 9) * pow(p0m3, 1)
    # p2m4 = pow(p1m0, 12) * pow(p0m1, 24) * pow(p2m0, -3) * pow(p1m1, -24) * pow(p0m2, -18) * pow(p3m0, 0) * pow(p2m1, 4) * pow(p1m2, 12) * pow(p0m3, 4)
    # p1m5 = pow(p1m0, 6) * pow(p0m1, 30) * pow(p2m0, 0) * pow(p1m1, -15) * pow(p0m2, -30) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 10) * pow(p0m3, 10)
    # p0m6 = pow(p1m0, 0) * pow(p0m1, 36) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, -45) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 20)
    # ODEs
    dp1m0dt = _.kn*p0m3 - _.kf*p0m1 - (2*3-1)*_.kf*p1m0 + _.kf*_.mt + _.k2*_.mt*1 - _.k2*p0m1
    dp0m1dt = (-3)*_.kn*p0m3 - 2*_.kp*p1m1 + 2*_.km*p1m0 + 3*(3-1)*_.kf*p1m0 - 0*_.k2*_.mt*1 + 0*_.k2*p0m1
    dp2m0dt = 2*(_.kn*p1m3 - _.kf*p1m1 - (2*3-1)*_.kf*p2m0 + _.mt*_.kf*p1m0 + _.k2*_.mt*p1m0 - _.k2*p1m1) +\
              1*(_.kn*p0m3 - _.kf*p0m1 - (2*3-1)*_.kf*p1m0 + _.mt*_.kf*1 + _.k2*_.mt*1 - _.k2*p0m1)
    dp1m1dt = ((-3)*_.kn*p1m3 - 2*_.kp*p2m1 + 2*_.km*p2m0 + _.kf*3*(3-1)*p2m0 - 0*_.k2*_.mt*p1m0 + 0*_.k2*p1m1) +\
              (_.kn*p0m4 - _.kf*p0m2 - (2*3-1)*_.kf*p1m1 + _.mt*_.kf*p0m1 + _.k2*_.mt*p0m1 - _.k2*p0m2) + ((-3)*_.kn * p0m3 - 0*_.k2*_.mt*1 + 0*_.k2*p0m1)
    dp0m2dt = 2*((-3)*_.kn*p0m4 - 2*_.kp*p1m2 + 2*_.km*p1m1 + _.kf*3*(3-1)*p1m1 - 0*_.k2*_.mt*p0m1 + 0*_.k2*p0m2) +\
              1*(pow( 3, 2) * _.kn*p0m3 + 2*_.kp*p1m1 + 2*_.km*p1m0 + 3*(3-1)*(2*3-1)*_.kf*p1m0/3. + pow(0, 2) *_.k2*_.mt*1 - pow(0, 2)*_.k2*p0m1)
    dp3m0dt = 3*(_.kn*p2m3 - _.kf*p2m1 - (2*3-1)*_.kf*p3m0 + _.mt*_.kf*p2m0 + _.k2*_.mt*p2m0 - _.k2*p2m1) +\
              3*(_.kn*p1m3 - _.kf*p1m1 - (2*3-1)*_.kf*p2m0 + _.mt*_.kf*p1m0 + _.k2*_.mt*p1m0 - _.k2*p1m1)
    dp2m1dt = 1*((-3)*_.kn*p2m3 - 2*_.kp*p3m1 + 2*_.km*p3m0 + _.kf*3*(3-1)*p3m0 - 0*_.k2*_.mt*p2m0 + 0*_.k2*p2m1) +\
              2*(_.kn*p1m4 - _.kf*p1m2 - (2*3-1)*_.kf*p2m1 + _.mt*_.kf*p1m1 + _.k2*_.mt*p1m1 - _.k2*p1m2) +\
              1*(_.kn*p0m4 - _.kf*p0m2 - (2*3-1)*_.kf*p1m1 + _.mt*_.kf*p0m1 + _.k2*_.mt*p0m1 - _.k2*p0m2) + 2*1*((-3)*_.kn * p1m3 - 0*_.k2*_.mt*p1m0 + 0*_.k2*p1m1)
    dp1m2dt = 2*((-3)*_.kn*p1m4 - 2*_.kp*p2m2 + 2*_.km*p2m1 + _.kf*3*(3-1)*p2m1 - 0*_.k2*_.mt*p1m1 + 0*_.k2*p1m2) +\
              1*(_.kn*p0m5 - _.kf*p0m3 - (2*3-1)*_.kf*p1m2 + _.mt*_.kf*p0m2 + _.k2*_.mt*p0m2 - _.k2*p0m3) +\
              1*(pow( 3, 2) * _.kn*p1m3 + 2*_.kp*p2m1 + 2*_.km*p2m0 + 3*(3-1)*(2*3-1)*_.kf*p2m0/3. +\
                 pow(0, 2) *_.k2*_.mt*p1m0 - pow(0, 2)*_.k2*p1m1) + 1*2*((-3)*_.kn * p0m4 - 0*_.k2*_.mt*p0m1 + 0*_.k2*p0m2)
    dp0m3dt = 3*((-3)*_.kn*p0m5 - 2*_.kp*p1m3 + 2*_.km*p1m2 + _.kf*3*(3-1)*p1m2 - 0*_.k2*_.mt*p0m2 + 0*_.k2*p0m3) +\
              3*(pow( 3, 2) * _.kn*p0m4 + 2*_.kp*p1m2 + 2*_.km*p1m1 + 3*(3-1)*(2*3-1)*_.kf*p1m1/3. + pow(0, 2) *_.k2*_.mt*p0m1 - pow(0, 2)*_.k2*p0m2)
    dzdt = [dp1m0dt, dp0m1dt, dp2m0dt, dp1m1dt, dp0m2dt, dp3m0dt, dp2m1dt, dp1m2dt, dp0m3dt]
    return dzdt


def stochastic_close3_nc3_n22(z, t, para):
    _ = parameters(*para)
    p1m0, p0m1, p2m0, p1m1, p0m2, p3m0, p2m1, p1m2, p0m3, = z
    # higher moments in terms of lower moments
    # p4m0 = pow(p1m0, 4) * pow(p0m1, 0) * pow(p2m0, -6) * pow(p1m1, 0) * pow(p0m2, 0) * pow(p3m0, 4) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 0)
    p3m1 = pow(p1m0, 3) * pow(p0m1, 1) * pow(p2m0, -3) * pow(p1m1, -3) * pow(p0m2, 0) * pow(p3m0, 1) * pow(p2m1, 3) * pow(p1m2, 0) * pow(p0m3, 0)
    p2m2 = pow(p1m0, 2) * pow(p0m1, 2) * pow(p2m0, -1) * pow(p1m1, -4) * pow(p0m2, -1) * pow(p3m0, 0) * pow(p2m1, 2) * pow(p1m2, 2) * pow(p0m3, 0)
    p1m3 = pow(p1m0, 1) * pow(p0m1, 3) * pow(p2m0, 0) * pow(p1m1, -3) * pow(p0m2, -3) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 3) * pow(p0m3, 1)
    p0m4 = pow(p1m0, 0) * pow(p0m1, 4) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, -6) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 4)
    # p5m0 = pow(p1m0, 15) * pow(p0m1, 0) * pow(p2m0, -20) * pow(p1m1, 0) * pow(p0m2, 0) * pow(p3m0, 10) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 0)
    # p4m1 = pow(p1m0, 12) * pow(p0m1, 3) * pow(p2m0, -12) * pow(p1m1, -8) * pow(p0m2, 0) * pow(p3m0, 4) * pow(p2m1, 6) * pow(p1m2, 0) * pow(p0m3, 0)
    # p3m2 = pow(p1m0, 9) * pow(p0m1, 6) * pow(p2m0, -6) * pow(p1m1, -12) * pow(p0m2, -2) * pow(p3m0, 1) * pow(p2m1, 6) * pow(p1m2, 3) * pow(p0m3, 0)
    p2m3 = pow(p1m0, 6) * pow(p0m1, 9) * pow(p2m0, -2) * pow(p1m1, -12) * pow(p0m2, -6) * pow(p3m0, 0) * pow(p2m1, 3) * pow(p1m2, 6) * pow(p0m3, 1)
    p1m4 = pow(p1m0, 3) * pow(p0m1, 12) * pow(p2m0, 0) * pow(p1m1, -8) * pow(p0m2, -12) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 6) * pow(p0m3, 4)
    p0m5 = pow(p1m0, 0) * pow(p0m1, 15) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, -20) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 10)
    # p6m0 = pow(p1m0, 36) * pow(p0m1, 0) * pow(p2m0, -45) * pow(p1m1, 0) * pow(p0m2, 0) * pow(p3m0, 20) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 0)
    # p5m1 = pow(p1m0, 30) * pow(p0m1, 6) * pow(p2m0, -30) * pow(p1m1, -15) * pow(p0m2, 0) * pow(p3m0, 10) * pow(p2m1, 10) * pow(p1m2, 0) * pow(p0m3, 0)
    # p4m2 = pow(p1m0, 24) * pow(p0m1, 12) * pow(p2m0, -18) * pow(p1m1, -24) * pow(p0m2, -3) * pow(p3m0, 4) * pow(p2m1, 12) * pow(p1m2, 4) * pow(p0m3, 0)
    # p3m3 = pow(p1m0, 18) * pow(p0m1, 18) * pow(p2m0, -9) * pow(p1m1, -27) * pow(p0m2, -9) * pow(p3m0, 1) * pow(p2m1, 9) * pow(p1m2, 9) * pow(p0m3, 1)
    # p2m4 = pow(p1m0, 12) * pow(p0m1, 24) * pow(p2m0, -3) * pow(p1m1, -24) * pow(p0m2, -18) * pow(p3m0, 0) * pow(p2m1, 4) * pow(p1m2, 12) * pow(p0m3, 4)
    # p1m5 = pow(p1m0, 6) * pow(p0m1, 30) * pow(p2m0, 0) * pow(p1m1, -15) * pow(p0m2, -30) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 10) * pow(p0m3, 10)
    # p0m6 = pow(p1m0, 0) * pow(p0m1, 36) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, -45) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 20)
    # ODEs
    dp1m0dt = _.kn*p0m3 - _.kf*p0m1 - (2*3-1)*_.kf*p1m0 + _.kf*_.mt + _.k2*_.mt*p0m2 - _.k2*p0m3
    dp0m1dt = (-3)*_.kn*p0m3 - 2*_.kp*p1m1 + 2*_.km*p1m0 + 3*(3-1)*_.kf*p1m0 - 2*_.k2*_.mt*p0m2 + 2*_.k2*p0m3
    dp2m0dt = 2*(_.kn*p1m3 - _.kf*p1m1 - (2*3-1)*_.kf*p2m0 + _.mt*_.kf*p1m0 + _.k2*_.mt*p1m2 - _.k2*p1m3) +\
              1*(_.kn*p0m3 - _.kf*p0m1 - (2*3-1)*_.kf*p1m0 + _.mt*_.kf*1 + _.k2*_.mt*p0m2 - _.k2*p0m3)
    dp1m1dt = ((-3)*_.kn*p1m3 - 2*_.kp*p2m1 + 2*_.km*p2m0 + _.kf*3*(3-1)*p2m0 - 2*_.k2*_.mt*p1m2 + 2*_.k2*p1m3) +\
              (_.kn*p0m4 - _.kf*p0m2 - (2*3-1)*_.kf*p1m1 + _.mt*_.kf*p0m1 + _.k2*_.mt*p0m3 - _.k2*p0m4) + ((-3)*_.kn * p0m3 - 2*_.k2*_.mt*p0m2 + 2*_.k2*p0m3)
    dp0m2dt = 2*((-3)*_.kn*p0m4 - 2*_.kp*p1m2 + 2*_.km*p1m1 + _.kf*3*(3-1)*p1m1 - 2*_.k2*_.mt*p0m3 + 2*_.k2*p0m4) +\
              1*(pow( 3, 2) * _.kn*p0m3 + 2*_.kp*p1m1 + 2*_.km*p1m0 + 3*(3-1)*(2*3-1)*_.kf*p1m0/3. + pow(2, 2) *_.k2*_.mt*p0m2 - pow(2, 2)*_.k2*p0m3)
    dp3m0dt = 3*(_.kn*p2m3 - _.kf*p2m1 - (2*3-1)*_.kf*p3m0 + _.mt*_.kf*p2m0 + _.k2*_.mt*p2m2 - _.k2*p2m3) +\
              3*(_.kn*p1m3 - _.kf*p1m1 - (2*3-1)*_.kf*p2m0 + _.mt*_.kf*p1m0 + _.k2*_.mt*p1m2 - _.k2*p1m3)
    dp2m1dt = 1*((-3)*_.kn*p2m3 - 2*_.kp*p3m1 + 2*_.km*p3m0 + _.kf*3*(3-1)*p3m0 - 2*_.k2*_.mt*p2m2 + 2*_.k2*p2m3) +\
              2*(_.kn*p1m4 - _.kf*p1m2 - (2*3-1)*_.kf*p2m1 + _.mt*_.kf*p1m1 + _.k2*_.mt*p1m3 - _.k2*p1m4) +\
              1*(_.kn*p0m4 - _.kf*p0m2 - (2*3-1)*_.kf*p1m1 + _.mt*_.kf*p0m1 + _.k2*_.mt*p0m3 - _.k2*p0m4) + 2*1*((-3)*_.kn * p1m3 - 2*_.k2*_.mt*p1m2 + 2*_.k2*p1m3)
    dp1m2dt = 2*((-3)*_.kn*p1m4 - 2*_.kp*p2m2 + 2*_.km*p2m1 + _.kf*3*(3-1)*p2m1 - 2*_.k2*_.mt*p1m3 + 2*_.k2*p1m4) +\
              1*(_.kn*p0m5 - _.kf*p0m3 - (2*3-1)*_.kf*p1m2 + _.mt*_.kf*p0m2 + _.k2*_.mt*p0m4 - _.k2*p0m5) +\
              1*(pow( 3, 2) * _.kn*p1m3 + 2*_.kp*p2m1 + 2*_.km*p2m0 + 3*(3-1)*(2*3-1)*_.kf*p2m0/3. +\
                 pow(2, 2) *_.k2*_.mt*p1m2 - pow(2, 2)*_.k2*p1m3) + 1*2*((-3)*_.kn * p0m4 - 2*_.k2*_.mt*p0m3 + 2*_.k2*p0m4)
    dp0m3dt = 3*((-3)*_.kn*p0m5 - 2*_.kp*p1m3 + 2*_.km*p1m2 + _.kf*3*(3-1)*p1m2 - 2*_.k2*_.mt*p0m4 + 2*_.k2*p0m5) +\
              3*(pow( 3, 2) * _.kn*p0m4 + 2*_.kp*p1m2 + 2*_.km*p1m1 + 3*(3-1)*(2*3-1)*_.kf*p1m1/3. + pow(2, 2) *_.k2*_.mt*p0m3 - pow(2, 2)*_.k2*p0m4)
    dzdt = [dp1m0dt, dp0m1dt, dp2m0dt, dp1m1dt, dp0m2dt, dp3m0dt, dp2m1dt, dp1m2dt, dp0m3dt]
    return dzdt


# =============================
# ====== nc = 4 ===============
# =============================
# Close at 2nd order. Use z0_close2 as initial condition:


def stochastic_close2_nc4_n20(z, t, para):
    _ = parameters(*para)
    p1m0, p0m1, p2m0, p1m1, p0m2, = z
    # higher moments in terms of lower moments
    # p3m0 = pow(p1m0, -3) * pow(p0m1, 0) * pow(p2m0, 3) * pow(p1m1, 0) * pow(p0m2, 0)
    p2m1 = pow(p1m0, -2) * pow(p0m1, -1) * pow(p2m0, 1) * pow(p1m1, 2) * pow(p0m2, 0)
    p1m2 = pow(p1m0, -1) * pow(p0m1, -2) * pow(p2m0, 0) * pow(p1m1, 2) * pow(p0m2, 1)
    # p0m3 = pow(p1m0, 0) * pow(p0m1, -3) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, 3)
    # p4m0 = pow(p1m0, -8) * pow(p0m1, 0) * pow(p2m0, 6) * pow(p1m1, 0) * pow(p0m2, 0)
    # p3m1 = pow(p1m0, -6) * pow(p0m1, -2) * pow(p2m0, 3) * pow(p1m1, 3) * pow(p0m2, 0)
    # p2m2 = pow(p1m0, -4) * pow(p0m1, -4) * pow(p2m0, 1) * pow(p1m1, 4) * pow(p0m2, 1)
    # p1m3 = pow(p1m0, -2) * pow(p0m1, -6) * pow(p2m0, 0) * pow(p1m1, 3) * pow(p0m2, 3)
    p0m4 = pow(p1m0, 0) * pow(p0m1, -8) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, 6)
    # p5m0 = pow(p1m0, -15) * pow(p0m1, 0) * pow(p2m0, 10) * pow(p1m1, 0) * pow(p0m2, 0)
    # p4m1 = pow(p1m0, -12) * pow(p0m1, -3) * pow(p2m0, 6) * pow(p1m1, 4) * pow(p0m2, 0)
    # p3m2 = pow(p1m0, -9) * pow(p0m1, -6) * pow(p2m0, 3) * pow(p1m1, 6) * pow(p0m2, 1)
    # p2m3 = pow(p1m0, -6) * pow(p0m1, -9) * pow(p2m0, 1) * pow(p1m1, 6) * pow(p0m2, 3)
    p1m4 = pow(p1m0, -3) * pow(p0m1, -12) * pow(p2m0, 0) * pow(p1m1, 4) * pow(p0m2, 6)
    p0m5 = pow(p1m0, 0) * pow(p0m1, -15) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, 10)
    # p6m0 = pow(p1m0, -24) * pow(p0m1, 0) * pow(p2m0, 15) * pow(p1m1, 0) * pow(p0m2, 0)
    # p5m1 = pow(p1m0, -20) * pow(p0m1, -4) * pow(p2m0, 10) * pow(p1m1, 5) * pow(p0m2, 0)
    # p4m2 = pow(p1m0, -16) * pow(p0m1, -8) * pow(p2m0, 6) * pow(p1m1, 8) * pow(p0m2, 1)
    # p3m3 = pow(p1m0, -12) * pow(p0m1, -12) * pow(p2m0, 3) * pow(p1m1, 9) * pow(p0m2, 3)
    # p2m4 = pow(p1m0, -8) * pow(p0m1, -16) * pow(p2m0, 1) * pow(p1m1, 8) * pow(p0m2, 6)
    # p1m5 = pow(p1m0, -4) * pow(p0m1, -20) * pow(p2m0, 0) * pow(p1m1, 5) * pow(p0m2, 10)
    # p0m6 = pow(p1m0, 0) * pow(p0m1, -24) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, 15)
    # ODEs
    dp1m0dt = _.kn*p0m4 - _.kf*p0m1 - (2*4-1)*_.kf*p1m0 + _.kf*_.mt + _.k2*_.mt*1 - _.k2*p0m1
    dp0m1dt = (-4)*_.kn*p0m4 - 2*_.kp*p1m1 + 2*_.km*p1m0 + 4*(4-1)*_.kf*p1m0 - 0*_.k2*_.mt*1 + 0*_.k2*p0m1
    dp2m0dt = 2*(_.kn*p1m4 - _.kf*p1m1 - (2*4-1)*_.kf*p2m0 + _.mt*_.kf*p1m0 + _.k2*_.mt*p1m0 - _.k2*p1m1) +\
              1*(_.kn*p0m4 - _.kf*p0m1 - (2*4-1)*_.kf*p1m0 + _.mt*_.kf*1 + _.k2*_.mt*1 - _.k2*p0m1)
    dp1m1dt = ((-4)*_.kn*p1m4 - 2*_.kp*p2m1 + 2*_.km*p2m0 + _.kf*4*(4-1)*p2m0 - 0*_.k2*_.mt*p1m0 + 0*_.k2*p1m1) +\
              (_.kn*p0m5 - _.kf*p0m2 - (2*4-1)*_.kf*p1m1 + _.mt*_.kf*p0m1 + _.k2*_.mt*p0m1 - _.k2*p0m2) + ((-4)*_.kn * p0m4 - 0*_.k2*_.mt*1 + 0*_.k2*p0m1)
    dp0m2dt = 2*((-4)*_.kn*p0m5 - 2*_.kp*p1m2 + 2*_.km*p1m1 + _.kf*4*(4-1)*p1m1 - 0*_.k2*_.mt*p0m1 + 0*_.k2*p0m2) +\
              1*(pow( 4, 2) * _.kn*p0m4 + 2*_.kp*p1m1 + 2*_.km*p1m0 + 4*(4-1)*(2*4-1)*_.kf*p1m0/3. + pow(0, 2) *_.k2*_.mt*1 - pow(0, 2)*_.k2*p0m1)
    dzdt = [dp1m0dt, dp0m1dt, dp2m0dt, dp1m1dt, dp0m2dt]
    return dzdt


# Close at 3rd order. Use z0_close3 as initial condition:


def stochastic_close3_nc4_n20(z, t, para):
    _ = parameters(*para)
    p1m0, p0m1, p2m0, p1m1, p0m2, p3m0, p2m1, p1m2, p0m3, = z
    # higher moments in terms of lower moments
    # p4m0 = pow(p1m0, 4) * pow(p0m1, 0) * pow(p2m0, -6) * pow(p1m1, 0) * pow(p0m2, 0) * pow(p3m0, 4) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 0)
    p3m1 = pow(p1m0, 3) * pow(p0m1, 1) * pow(p2m0, -3) * pow(p1m1, -3) * pow(p0m2, 0) * pow(p3m0, 1) * pow(p2m1, 3) * pow(p1m2, 0) * pow(p0m3, 0)
    p2m2 = pow(p1m0, 2) * pow(p0m1, 2) * pow(p2m0, -1) * pow(p1m1, -4) * pow(p0m2, -1) * pow(p3m0, 0) * pow(p2m1, 2) * pow(p1m2, 2) * pow(p0m3, 0)
    p1m3 = pow(p1m0, 1) * pow(p0m1, 3) * pow(p2m0, 0) * pow(p1m1, -3) * pow(p0m2, -3) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 3) * pow(p0m3, 1)
    p0m4 = pow(p1m0, 0) * pow(p0m1, 4) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, -6) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 4)
    # p5m0 = pow(p1m0, 15) * pow(p0m1, 0) * pow(p2m0, -20) * pow(p1m1, 0) * pow(p0m2, 0) * pow(p3m0, 10) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 0)
    # p4m1 = pow(p1m0, 12) * pow(p0m1, 3) * pow(p2m0, -12) * pow(p1m1, -8) * pow(p0m2, 0) * pow(p3m0, 4) * pow(p2m1, 6) * pow(p1m2, 0) * pow(p0m3, 0)
    # p3m2 = pow(p1m0, 9) * pow(p0m1, 6) * pow(p2m0, -6) * pow(p1m1, -12) * pow(p0m2, -2) * pow(p3m0, 1) * pow(p2m1, 6) * pow(p1m2, 3) * pow(p0m3, 0)
    # p2m3 = pow(p1m0, 6) * pow(p0m1, 9) * pow(p2m0, -2) * pow(p1m1, -12) * pow(p0m2, -6) * pow(p3m0, 0) * pow(p2m1, 3) * pow(p1m2, 6) * pow(p0m3, 1)
    p1m4 = pow(p1m0, 3) * pow(p0m1, 12) * pow(p2m0, 0) * pow(p1m1, -8) * pow(p0m2, -12) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 6) * pow(p0m3, 4)
    p0m5 = pow(p1m0, 0) * pow(p0m1, 15) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, -20) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 10)
    # p6m0 = pow(p1m0, 36) * pow(p0m1, 0) * pow(p2m0, -45) * pow(p1m1, 0) * pow(p0m2, 0) * pow(p3m0, 20) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 0)
    # p5m1 = pow(p1m0, 30) * pow(p0m1, 6) * pow(p2m0, -30) * pow(p1m1, -15) * pow(p0m2, 0) * pow(p3m0, 10) * pow(p2m1, 10) * pow(p1m2, 0) * pow(p0m3, 0)
    # p4m2 = pow(p1m0, 24) * pow(p0m1, 12) * pow(p2m0, -18) * pow(p1m1, -24) * pow(p0m2, -3) * pow(p3m0, 4) * pow(p2m1, 12) * pow(p1m2, 4) * pow(p0m3, 0)
    # p3m3 = pow(p1m0, 18) * pow(p0m1, 18) * pow(p2m0, -9) * pow(p1m1, -27) * pow(p0m2, -9) * pow(p3m0, 1) * pow(p2m1, 9) * pow(p1m2, 9) * pow(p0m3, 1)
    p2m4 = pow(p1m0, 12) * pow(p0m1, 24) * pow(p2m0, -3) * pow(p1m1, -24) * pow(p0m2, -18) * pow(p3m0, 0) * pow(p2m1, 4) * pow(p1m2, 12) * pow(p0m3, 4)
    p1m5 = pow(p1m0, 6) * pow(p0m1, 30) * pow(p2m0, 0) * pow(p1m1, -15) * pow(p0m2, -30) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 10) * pow(p0m3, 10)
    p0m6 = pow(p1m0, 0) * pow(p0m1, 36) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, -45) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 20)
    # p7m0 = pow(p1m0, 70) * pow(p0m1, 0) * pow(p2m0, -84) * pow(p1m1, 0) * pow(p0m2, 0) * pow(p3m0, 35) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 0)
    # p6m1 = pow(p1m0, 60) * pow(p0m1, 10) * pow(p2m0, -60) * pow(p1m1, -24) * pow(p0m2, 0) * pow(p3m0, 20) * pow(p2m1, 15) * pow(p1m2, 0) * pow(p0m3, 0)
    # p5m2 = pow(p1m0, 50) * pow(p0m1, 20) * pow(p2m0, -40) * pow(p1m1, -40) * pow(p0m2, -4) * pow(p3m0, 10) * pow(p2m1, 20) * pow(p1m2, 5) * pow(p0m3, 0)
    # p4m3 = pow(p1m0, 40) * pow(p0m1, 30) * pow(p2m0, -24) * pow(p1m1, -48) * pow(p0m2, -12) * pow(p3m0, 4) * pow(p2m1, 18) * pow(p1m2, 12) * pow(p0m3, 1)
    # p3m4 = pow(p1m0, 30) * pow(p0m1, 40) * pow(p2m0, -12) * pow(p1m1, -48) * pow(p0m2, -24) * pow(p3m0, 1) * pow(p2m1, 12) * pow(p1m2, 18) * pow(p0m3, 4)
    # p2m5 = pow(p1m0, 20) * pow(p0m1, 50) * pow(p2m0, -4) * pow(p1m1, -40) * pow(p0m2, -40) * pow(p3m0, 0) * pow(p2m1, 5) * pow(p1m2, 20) * pow(p0m3, 10)
    # p1m6 = pow(p1m0, 10) * pow(p0m1, 60) * pow(p2m0, 0) * pow(p1m1, -24) * pow(p0m2, -60) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 15) * pow(p0m3, 20)
    # p0m7 = pow(p1m0, 0) * pow(p0m1, 70) * pow(p2m0, 0) * pow(p1m1, 0) * pow(p0m2, -84) * pow(p3m0, 0) * pow(p2m1, 0) * pow(p1m2, 0) * pow(p0m3, 35)
    # ODEs
    dp1m0dt = _.kn*p0m4 - _.kf*p0m1 - (2*4-1)*_.kf*p1m0 + _.kf*_.mt + _.k2*_.mt*1 - _.k2*p0m1
    dp0m1dt = (-4)*_.kn*p0m4 - 2*_.kp*p1m1 + 2*_.km*p1m0 + 4*(4-1)*_.kf*p1m0 - 0*_.k2*_.mt*1 + 0*_.k2*p0m1
    dp2m0dt = 2*(_.kn*p1m4 - _.kf*p1m1 - (2*4-1)*_.kf*p2m0 + _.mt*_.kf*p1m0 + _.k2*_.mt*p1m0 - _.k2*p1m1) +\
              1*(_.kn*p0m4 - _.kf*p0m1 - (2*4-1)*_.kf*p1m0 + _.mt*_.kf*1 + _.k2*_.mt*1 - _.k2*p0m1)
    dp1m1dt = ((-4)*_.kn*p1m4 - 2*_.kp*p2m1 + 2*_.km*p2m0 + _.kf*4*(4-1)*p2m0 - 0*_.k2*_.mt*p1m0 + 0*_.k2*p1m1) +\
              (_.kn*p0m5 - _.kf*p0m2 - (2*4-1)*_.kf*p1m1 + _.mt*_.kf*p0m1 + _.k2*_.mt*p0m1 - _.k2*p0m2) + ((-4)*_.kn * p0m4 - 0*_.k2*_.mt*1 + 0*_.k2*p0m1)
    dp0m2dt = 2*((-4)*_.kn*p0m5 - 2*_.kp*p1m2 + 2*_.km*p1m1 + _.kf*4*(4-1)*p1m1 - 0*_.k2*_.mt*p0m1 +\
                 0*_.k2*p0m2) + 1*(pow( 4, 2) * _.kn*p0m4 + 2*_.kp*p1m1 + 2*_.km*p1m0 + 4*(4-1)*(2*4-1)*_.kf*p1m0/3. + pow(0, 2) *_.k2*_.mt*1 - pow(0, 2)*_.k2*p0m1)
    dp3m0dt = 3*(_.kn*p2m4 - _.kf*p2m1 - (2*4-1)*_.kf*p3m0 + _.mt*_.kf*p2m0 + _.k2*_.mt*p2m0 - _.k2*p2m1) +\
              3*(_.kn*p1m4 - _.kf*p1m1 - (2*4-1)*_.kf*p2m0 + _.mt*_.kf*p1m0 + _.k2*_.mt*p1m0 - _.k2*p1m1)
    dp2m1dt = 1*((-4)*_.kn*p2m4 - 2*_.kp*p3m1 + 2*_.km*p3m0 + _.kf*4*(4-1)*p3m0 - 0*_.k2*_.mt*p2m0 + 0*_.k2*p2m1) +\
              2*(_.kn*p1m5 - _.kf*p1m2 - (2*4-1)*_.kf*p2m1 + _.mt*_.kf*p1m1 + _.k2*_.mt*p1m1 - _.k2*p1m2) +\
              1*(_.kn*p0m5 - _.kf*p0m2 - (2*4-1)*_.kf*p1m1 + _.mt*_.kf*p0m1 + _.k2*_.mt*p0m1 - _.k2*p0m2) + 2*1*((-4)*_.kn * p1m4 - 0*_.k2*_.mt*p1m0 + 0*_.k2*p1m1)
    dp1m2dt = 2*((-4)*_.kn*p1m5 - 2*_.kp*p2m2 + 2*_.km*p2m1 + _.kf*4*(4-1)*p2m1 - 0*_.k2*_.mt*p1m1 + 0*_.k2*p1m2) +\
              1*(_.kn*p0m6 - _.kf*p0m3 - (2*4-1)*_.kf*p1m2 + _.mt*_.kf*p0m2 + _.k2*_.mt*p0m2 - _.k2*p0m3) +\
              1*(pow( 4, 2) * _.kn*p1m4 + 2*_.kp*p2m1 + 2*_.km*p2m0 + 4*(4-1)*(2*4-1)*_.kf*p2m0/3. +\
                 pow(0, 2) *_.k2*_.mt*p1m0 - pow(0, 2)*_.k2*p1m1) + 1*2*((-4)*_.kn * p0m5 - 0*_.k2*_.mt*p0m1 + 0*_.k2*p0m2)
    dp0m3dt = 3*((-4)*_.kn*p0m6 - 2*_.kp*p1m3 + 2*_.km*p1m2 + _.kf*4*(4-1)*p1m2 - 0*_.k2*_.mt*p0m2 + 0*_.k2*p0m3) +\
              3*(pow( 4, 2) * _.kn*p0m5 + 2*_.kp*p1m2 + 2*_.km*p1m1 + 4*(4-1)*(2*4-1)*_.kf*p1m1/3. + pow(0, 2) *_.k2*_.mt*p0m1 - pow(0, 2)*_.k2*p0m2)
    dzdt = [dp1m0dt, dp0m1dt, dp2m0dt, dp1m1dt, dp0m2dt, dp3m0dt, dp2m1dt, dp1m2dt, dp0m3dt]
    return dzdt


