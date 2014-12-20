from nlsymb import np, sym, matmult, tensor
import nlsymb.tensor as tn

from sympy import Symbol as S, sin, cos, ccode
import re


def ReplaceXU( s, x, u ):
    '''Takes symbolic vectors x and u and replaces each of these elements
    in string s with suitable C/C++ implementable versions.
    e.g: x1 -> x[0], u1 -> u[0]'''
    for i in range(0,len(x)):
        xistr = str(x[i])[0:1] + '[' + str(i) + ']'
        s = re.sub(str(x[i]), xistr, s)
    for i in range(0,len(u)):
        uistr = str(u[i])[0:1] + '[' + str(i) + ']'
        s = re.sub(str(u[i]), uistr, s)
    return s


def WriteDynamics( dyn, x, u ):
    xlen = dyn.size
    f = open('./output/sys_dynam.hpp', 'w')
    f.write('#ifndef SYS_DYNAM_HPP\r\n')
    f.write('#define SYS_DYNAM_HPP\r\n\r\n')
    f.write('//[ The rhs of xdot = f(x) defined as a class\r\n')
    f.write('// USER SPECIFIED:\r\n')
    f.write('class sys_dynam {\r\n')
    f.write('  b_control & m_u;\r\n')
    f.write('  state_type u;\r\n\r\n')
    f.write('public:\r\n')
    f.write('  sys_dynam( b_control & uu ) : m_u(uu) , u(ulen) {  }\r\n\r\n')
    f.write('  void operator() (const state_type &x, state_type &dxdt,')
    f.write(' const double t)\r\n')
    f.write('  {\r\n')
    f.write('    m_u(t, u);\r\n')
    f.write('    //\r\n')
    for i in range (0,xlen):
        s = ccode(dyn[i])
        s = ReplaceXU( s, x, u )
        f.write('    dxdt[' + str(i) + '] = ' + s + ';\r\n')
    f.write('  }\r\n')
    f.write('};\r\n')
    f.write('//]\r\n\r\n')
    f.write('#endif  // SYS_DYNAM_HPP')
    f.close()


def WriteLinearizations( A, B, x, u ):
    xlen = (A.shape)[0]
    ulen = (B.shape)[1]
    f = open('./output/sys_lin.hpp', 'w')
    f.write('#ifndef SYS_LIN_HPP\r\n')
    f.write('#define SYS_LIN_HPP\r\n\r\n')
    f.write('//[ Linearizations of the system defined as a class\r\n')
    f.write('// USER SPECIFIED:\r\n')
    f.write('class sys_lin {\r\n\r\n')
    f.write('public:\r\n')
    f.write('  sys_lin( ) {  }\r\n\r\n')
    f.write('  void A( const state_type & x, const state_type & u,\r\n')
    f.write('          Eigen::Matrix< double, xlen, xlen > & Amat ) {\r\n')
    for i in range (0,xlen):
        for j in range (0,xlen):
            s = ccode(A[i,j])
            s = ReplaceXU( s, x, u )
            f.write('    Amat(' + str(i) + ',' + str(j) + ') = '
                    + s + ';\r\n')
    f.write('  }\r\n\r\n')
    f.write('  void B( const state_type & x, const state_type & u,\r\n')
    f.write('          Eigen::Matrix< double, xlen, ulen > & Bmat ) {\r\n')
    for i in range (0,xlen):
        for j in range (0,ulen):
            s = ccode(B[i,j])
            s = ReplaceXU( s, x, u )
            f.write('    Bmat(' + str(i) + ',' + str(j) + ') = '
                    + s + ';\r\n')
    f.write('  }\r\n')
    f.write('};\r\n')
    f.write('//]\r\n\r\n')
    f.write('#endif  // SYS_LIN_HPP')
    f.close()
