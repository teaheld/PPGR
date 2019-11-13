'''
        Normalizovani DLT algoritam uradjen za dat test primer
        sa 6 tacaka.
'''

'''
        Inicijalizacija potrebnih biblioteka
'''
import numpy as np

def normalizacija_tacaka(tacke):
    '''
            I korak:
            Izracunati teziste sistema tacaka (afino)
    '''

    '''
            Delimo sa 3.koordinatom, posto se nalazimo u afinim koordinata
    '''
    Tx = sum([t[0] / t[2] for t in tacke]) / len(tacke)
    Ty = sum([t[1] / t[2] for t in tacke]) / len(tacke)

    '''
            II korak:
            Translirati teziste u koordinatni pocetak
    '''

    '''
            Matrica translacije G
    '''
    G = np.array([[1, 0, -Tx],
                  [0, 1, -Ty],
                  [0, 0,  1]])

    '''
            III korak:
            Skalirati tacke 
            tako da je prosecna udaljenost tacke od koordinatnog pocetka bude sqrt(2)
    '''
    '''
            Prvo trazimo prosecnu udaljenost od koordinatnog pocetka l
            Opet smo delili sa 3.koordinatom, jer se nalazimo u afinim koordinatama
    '''
    l = sum(np.sqrt((t[0] / t[2] - Tx) ** 2 + (t[1] / t[2] - Ty) ** 2) for t in tacke) / len(tacke)

    '''
            Odredjujemo koeficijent homotetije k
    '''
    k = np.sqrt(2) / l

    '''
            Na kraju je matrica homotetije S
    '''
    S = np.array([[k, 0, 0],
                  [0, k, 0],
                  [0, 0, 1]])

    '''
            IV korak:
            Trazena matrica normalizacije je T = S*G
    '''
    T = np.matmul(S, G)
    return T

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

    '''
            III korak:
            Odredjujemo SVD dekompoziciju matrice A = UDV^T
    '''
    U, D, V_T = np.linalg.svd(A)

    '''
            Izlaz:
            Resenje je poslednja kolona(vrsta) matrice V^T
            Zaokruzena na 6 decimala
    '''
    P = np.reshape(V_T[-1,:], (3,3))
    return P

def normalizovani_DLT(tacke):
    '''
            Ulaz:
            Homogene koordinate
            4 ili vise orignialnih tacaka
            i njihovih slika
    '''

    normalizovane = {}
    '''
            I korak:
            Normalizacija originalnih tacaka transformacijom T
    '''
    T = normalizacija_tacaka([tacke['A'], tacke['B'], tacke['C'], tacke['D'], tacke['E'], tacke['F']])
    print("T = \n", T.round(6))

    normalizovane['A'] = np.matmul(T, tacke['A'])
    print("\n Normalizovane ['A'] = ", normalizovane['A'])
    normalizovane['B'] = np.matmul(T, tacke['B'])
    normalizovane['C'] = np.matmul(T, tacke['C'])
    normalizovane['D'] = np.matmul(T, tacke['D'])
    normalizovane['E'] = np.matmul(T, tacke['E'])
    normalizovane['F'] = np.matmul(T, tacke['F'])

    '''
            II korak:
            Normalizacija slika transformacijom Tp
    '''
    Tp = normalizacija_tacaka( [tacke['Ap'], tacke['Bp'], tacke['Cp'], tacke['Dp'], tacke['Ep'], tacke['Fp']])
    print("Tp = \n", Tp.round(6))

    normalizovane['Ap'] = np.matmul(Tp, tacke['Ap'])
    normalizovane['Bp'] = np.matmul(Tp, tacke['Bp'])
    normalizovane['Cp'] = np.matmul(Tp, tacke['Cp'])
    normalizovane['Dp'] = np.matmul(Tp, tacke['Dp'])
    normalizovane['Ep'] = np.matmul(Tp, tacke['Ep'])
    normalizovane['Fp'] = np.matmul(Tp, tacke['Fp'])

    '''
            III korak:
            DLT algoritmom odrediti matricu transformacije Pn
            koja slika normalizovane tacke
            u normalizovane slike
    '''
    Pn = dlt(normalizovane)
    print("Pn = \n", Pn.round(6))

    '''
            IV korak:
            Trazena matrica transformacije je P = Tp^-1 * Pn * T
    '''
    Tp_inv = np.linalg.inv(Tp)
    P = Tp_inv.dot(Pn).dot(T)

    '''
            Izlaz:
            3 x 3 matrica P projektivnog preslikavanja
    '''
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

    print("P = \n", normalizovani_DLT(tacke).round(6))

if __name__ == '__main__':
    main()