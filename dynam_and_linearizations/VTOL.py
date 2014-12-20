from nlsymb import np, sym, matmult, tensor
import nlsymb.tensor as tn

from sympy import Symbol as S, sin, cos

from Cdynamics import WriteDynamics, WriteLinearizations

if __name__ == "__main__":
    # x, u
    x = np.array([S('x1'), S('x2'), S('x3'), S('x4'), S('x5'), S('x6')]) # x,y,theta, x',y',theta'
    u = np.array([S('u1'),S('u2')])                            # x''

    # constants
    g = S('g_')
    eps = S('eps_')

    # dynamics
    f = np.array([x[3], x[4], x[5], -u[0]*sin(x[2])+eps*u[1]*cos(x[2]), u[0]*cos(x[2])+eps*u[1]*sin(x[2])-g, u[1]])

    # linearizations
    A = tn.diff(f, x)
    B = tn.diff(f, u)

    # write to file
    WriteDynamics( f, x, u )
    WriteLinearizations( A, B, x, u )
