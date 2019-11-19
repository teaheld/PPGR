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
    fi = math.acos(np.dot(u, u_p))

    '''
            Da bi rotacija bila u pozitivnom smeru,
            ako je [u, u', p] < 0
            p = -p. 
    '''
    if (np.cross(u, u_p)).dot(p):
        p = -p

    '''
            Vraca jedinicni vektor p 
            i ugao fi 
            tako da A = Rp(fi).
    '''
    return (p, fi)

def Rodriqez(p, fi):
    '''
            Formule preuzete sa prezentacije, strana 7.
    '''
    ppt = np.matmul(p[np.newaxis].transpose(), p[np.newaxis])
    px = np.matrix([[0, -p[2], p[1]],
          [p[2], 0, -p[0]],
          [-p[1], p[0], 0]])

    return ppt + math.cos(fi) * (np.eye(3) - ppt) + math.sin(fi) * px

def A2Euler(A):
    '''
            Slajdovi, strana 22.
    '''
    fi = 0
    teta = math.pi / 2
    if A[2][0] < 1:
        if A[2][0] > -1:
            '''
                    Jedinstveno resenje.
            '''
            psi = math.atan2(A[1][0], A[0][0])
            teta = math.asin(-A[2][0])
            fi = math.atan2(A[2][1], A[2][2])
        else:
            '''
                    Nije jedinstveno. 
            '''
            psi = math.atan2(-A[0][1], A[1][1])
    else:
        psi = math.atan2(-A[0][1], A[1][1])

    return (fi, teta, psi)

def AxisAngle2Q(p, fi):
    '''
            Slajdovi, strana 34.
    '''
    w = math.cos(fi / 2)
    p = p / la.norm(p)
    x, y, z = math.sin(fi / 2) * p

    return (x, y, z, w)

def main():
    A = Euler2A(-np.arctan(1/4), -np.arcsin(8/9), np.arctan(4))
    p, fi = AxisAngle(A)
    Rp = Rodriqez(p, fi)
    E_fi, E_teta, E_psi = A2Euler(A)
    quaternions = AxisAngle2Q(p, fi)

if __name__ == '__main__':
	main()


