import numpy as np
import math

np.set_printoptions(precision=10, suppress=True)


def Euler2A(fi, teta, psi):
    '''
            Forumule preuzete sa prezentacije,
            strana 6.
    '''
    Rz = [[np.cos(psi), -np.sin(psi), 0],
          [np.sin(psi), np.cos(psi), 0],
          [0, 0, 1]]

    Ry = [[np.cos(teta), 0, np.sin(teta)],
          [0, 1, 0],
          [-np.sin(teta), 0, np.cos(teta)]]

    Rx = [[1, 0, 0],
          [0, np.cos(fi), -np.sin(fi)],
          [0, np.sin(fi), np.cos(fi)]]

    '''
            Vraca matricu
            A = Rz(psi)Ry(teta)Rx(fi)
    '''
    return np.matmul(np.matmul(Rz, Ry), Rx)


from numpy import linalg as la
def AxisAngle(A):

    '''
            Odredjujemo sopstvene vrednosti i vektore.
    '''
    w, v = la.eig(np.array(A))
    '''
            Odrediti jedinicni sopstveni vektor p
            za lambda = 1.
    '''
    '''
            Nalazimo lambda = 1.
    '''
    i_lambda = np.where(w.real.round() == 1)[0][0]
    p = v[:, i_lambda].real

    '''
            Odrediti proizvoljan jedinicni vektor u _|_ p.
            
            Posto je matrica ulaz matrica A koja je simetricna,
            Sopstveni vektori su normalni jedni na druge, 
            pa cemo uzeti proizvoljan sopstveni vektor koji nismo iskoristili.
    '''
    u = v[:, len(p) - i_lambda -1].real

    '''
            u' = Au,
            gde je u' jedinicni.
    '''
    u_p = np.dot(A, u)
    u_p = u_p / la.norm(u_p, 2)

    '''
            fi = arccos(u, u')
    '''
    fi = math.degrees(math.acos(np.dot(u, u_p)))

    '''
            Vraca jedinicni vektor p 
            i ugao fi 
            tako da A = Rp(fi).
    '''
    return (p, fi)

def main():
    A = Euler2A(-np.arctan(1/4), -np.arcsin(8/9), np.arctan(4))
    print(AxisAngle(A))

if __name__ == '__main__':
	main()


