import numpy as np
import math

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

def main():
    print(Euler2A(-np.arctan(1/4), -np.arcsin(8/9), np.arctan(4)))

if __name__ == '__main__':
	main()


