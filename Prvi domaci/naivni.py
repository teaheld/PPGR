'''
        Naivni algoritam uradjen za dati test primer
        za preslikavanje trapeza u pravougaonik.
'''

'''
        Ucitavanje potrebnih biblioteka
'''
import numpy as np

'''
        Naivni algoritam
'''
def naivni(tacke):
    '''
            Ulaz:
                    homogene koordinate
                    4 originalne tacke A,  B,  C,  D
                    i 4 njihove slike  Ap, Bp, Cp, Dp
    '''

    '''
            I korak:
            Odredjujemo alfa, beta i gama =/= 0
            D = alfa * A + beta * B + gama * C
            P1 je matrica sa kolonama alfa * A, beta * B, gama * C
    '''

    '''
            Matrica sa kolonama originalnih tacaka
            potrebna i za resavanje sistema 
            i kasnije za pravljenje P1
    '''
    or_kol = np.transpose([tacke['A'], tacke['B'], tacke['C']])

    ''' 
            D = alfa * A + beta * B + gama * C
                        <=>
            -1 = alfa * -3 + beta *  3 + gama * 1
             1 = alfa * -1 + beta * -1 + gama * 1 
             1 = alfa *  1 + beta *  1 + gama * 1
    '''
    alfa_beta_gama = np.round(np.linalg.solve(or_kol, np.transpose(tacke['D'])), 6)
    print("alfa = {0}, beta = {1}, gama = {2}".format(alfa_beta_gama[0], alfa_beta_gama[1], alfa_beta_gama[2]))

    '''
            P1 je matrica sa kolonama alfa * A, beta * B, gama * C
    '''
    P1 = alfa_beta_gama * or_kol
    print("P1 = \n", P1)

    '''
            II korak:
            Isto radimo i za Ap, Bp, Cp i Dp
    '''
    '''
            Matrica sa kolonama slika
    '''
    sl_kol = np.transpose([tacke['Ap'], tacke['Bp'], tacke['Cp']])

    '''
            Dp = alfap * Ap + betap * Bp + gamap * Cp
    '''
    alfap_betap_gamap = np.round(np.linalg.solve(sl_kol, np.transpose(tacke['Dp'])), 6)
    print("alfap = {0}, betap = {1}, gamap = {2}".format(alfap_betap_gamap[0], alfap_betap_gamap[1], alfap_betap_gamap[2]))

    '''
            P2 je matrica sa kolonama alfap * Ap, betap * Bp, gamap * Cp
    '''
    P2 = alfap_betap_gamap * sl_kol
    print("P2 = \n", P2.round(6))

    '''
            III korak:
            P = P2 * P1^-1 je trazena matrica preslikavanja
    '''

    ''' 
            Racunamo P1^-1
    '''
    P1_inv = np.linalg.inv(P1)
    print("P1^-1 = \n", P1_inv.round(6))

    ''' 
            P = P2 * P1^-1
    '''
    P = np.matmul(P2, P1_inv)

    '''      Izlaz:
                    3 x 3 matrica P projektivnog preslikavanja ravni
                    koje slika
                    A,  B,  C,  D  u
                    Ap, Bp, Cp, Dp
    '''
    return P


def main():
    tacke = {}

    tacke['A'] = (-3, -1, 1)
    tacke['B'] = (3, -1, 1)
    tacke['C'] = (1, 1, 1)
    tacke['D'] = (-1, 1, 1)

    tacke['Ap'] = (-2, -1, 1)
    tacke['Bp'] = (2, -1, 1)
    tacke['Cp'] = (2, 1, 1)
    tacke['Dp'] = (-2, 1, 1)

    print("P = \n", naivni(tacke).round(6))

if __name__ == '__main__':
	main()