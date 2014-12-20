from nlsymb import np, sym, matmult, tensor
import nlsymb.tensor as tn

from sympy import Symbol as S, sin, cos

from Cdynamics import WriteDynamics, WriteLinearizations

if __name__ == "__main__":
    x = np.array([S('x1'), S('x2'), S('x3')]) # x, y, theta
    u = np.array([S('u1'), S('u2')])          # v, omega
    f = np.array([cos(x[2])*u[0], sin(x[2])*u[0], u[1]])

    A = tn.diff(f, x)
    B = tn.diff(f, u)

    WriteDynamics( f, x, u )
    WriteLinearizations( A, B, x, u )
