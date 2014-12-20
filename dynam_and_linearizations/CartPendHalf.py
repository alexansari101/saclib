from nlsymb import np, sym, matmult, tensor
import nlsymb.tensor as tn

from sympy import Symbol as S, sin, cos

from Cdynamics import WriteDynamics, WriteLinearizations

if __name__ == "__main__":
    # x, u
    x = np.array([S('x1'), S('x2')])            # theta, theta'
    u = np.array([S('u1')])                     # x''

    # constants
    g = S('m_g')
    h = S('h')

    # dynamics
    f = np.array([x[1], (g*sin(x[0])+cos(x[0])*u[0])/h])

    # linearizations
    A = tn.diff(f, x)
    B = tn.diff(f, u)

    # write to file
    WriteDynamics( f, x, u )
    WriteLinearizations( A, B, x, u )
