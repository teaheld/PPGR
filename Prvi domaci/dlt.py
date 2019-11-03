'''
        DLT algoritam uradjen za dat test primer
        sa 6 tacaka.
'''

'''
        Ukljucivanje potrebnih biblioteka
'''
import numpy as np

def matrica_2_9(X, Xp):
    '''
            Ulaz: tacka i njena slika

            Izlaz: 2 x 9 matrica data u lemi (prezentacija 1, str 31)
    '''
    return [[0, 0, 0, -Xp[2] * X[0], -Xp[2] * X[1], -Xp[2] * X[2], Xp[1] * X[0], Xp[1] * X[1], Xp[1] * X[2]],
            [Xp[2] * X[0], Xp[2] * X[1], Xp[2] * X[2], 0, 0, 0, -Xp[0] * X[0], -Xp[0] * X[1], -Xp[0] * X[2]]]

def dlt(tacke):
    '''
            I korak:
            Za svaku korespodenciju tacaka odredjujemo 2 x 9 matricu
    '''
    M1 = matrica_2_9(tacke['A'], tacke['Ap'])
    M2 = matrica_2_9(tacke['B'], tacke['Bp'])
    M3 = matrica_2_9(tacke['C'], tacke['Cp'])
    M4 = matrica_2_9(tacke['D'], tacke['Dp'])
    M5 = matrica_2_9(tacke['E'], tacke['Ep'])
    M6 = matrica_2_9(tacke['F'], tacke['Fp'])

    '''
            II korak:
            Spajamo ove matrice u jednu matricu A formata 2n x 9    
    '''
    A = np.concatenate((M1, M2, M3, M4, M5, M6), axis=0)
    print("A = \n", A.round(6))

    '''
            III korak:
            Odredjujemo SVD dekompoziciju matrice A = UDV^T
    '''
    U, D, V_T = np.linalg.svd(A)
    print("U = \n", U)
    print("D = \n", D)
    print("V_T = \n", V_T)

    '''
            Izlaz:
            Resenje je poslednja kolona(vrsta) matrice V^T 
            Zaokruzena na 6 decimala
    '''
    P = np.reshape(V_T[-1,:], (3,3)).round(6)
    return P

def main():
    tacke = {}

    tacke['A'] = (-3, -1, 1)
    tacke['B'] = (3, -1, 1)
    tacke['C'] = (1, 1, 1)
    tacke['D'] = (-1, 1, 1)
    tacke['E'] = (1, 2, 3)
    tacke['F'] = (-8, -2, 1)

    tacke['Ap'] = (-2, -1, 1)
    tacke['Bp'] = (2, -1, 1)
    tacke['Cp'] = (2, 1, 1)
    tacke['Dp'] = (-2, 1, 1)
    tacke['Ep'] = (2, 1, 4)
    tacke['Fp'] = (-16, -5, 4)

    print("P = \n", dlt(tacke))

if __name__ == '__main__':
	main()