'''
        Ucitavamo potrebna zaglavlja.
'''
import numpy as np

def ucitaj_tacke(data):
    '''
            Fja prima putanju do fajla.
    '''

    '''
            Otvaramo fajl u kome se nalaze tacke.
    '''
    f = open(data, "r")

    '''
            Cuvamo tacke kao listu.
    '''
    t = []
    for line in f:
        if line.strip():
            '''
                    U svakoj liniji koja nije prazna se nalazi po jedna tacka.
                    u obliku:
                    -3, -1, 1
                    Cuvamo ih kao trojku.
            '''
            x1, x2, x3 = (float(x) for x in line.split(', '))
            t.append((x1, x2, x3))

    f.close()
    '''
            Fja vraca listu tacaka iz prosledjenog fajla.
    '''
    return t

def naivni(tacke):
    '''
            Ulaz:
                    homogene koordinate
                    4 originalne tacke A,  B,  C,  D
                    i 4 njihove slike  Ap, Bp, Cp, Dp.
            Mi prosledjujemo sve tacke, ali cemo manipulisati samo sa 4.
    '''

    '''
            I korak:
                    Odredjujemo alfa, beta i gama =/= 0
                    D = alfa * A + beta * B + gama * C.
                    P1 je matrica sa kolonama alfa * A, beta * B, gama * C.
    '''

    '''
            Matrica sa kolonama originalnih tacaka
            potrebna i za resavanje sistema 
            i kasnije za pravljenje P1.
    '''
    or_kol = np.transpose([tacke['A'], tacke['B'], tacke['C']])

    ''' 
            D = alfa * A + beta * B + gama * C
                        <=>
            -1 = alfa * -3 + beta *  3 + gama * 1
             1 = alfa * -1 + beta * -1 + gama * 1 
             1 = alfa *  1 + beta *  1 + gama * 1
    '''
    alfa_beta_gama = np.linalg.solve(or_kol, np.transpose(tacke['D']))

    '''
            P1 je matrica sa kolonama alfa * A, beta * B, gama * C.
    '''
    P1 = alfa_beta_gama * or_kol

    '''
            II korak:
                    Isto radimo i za Ap, Bp, Cp i Dp.
    '''
    '''
            Matrica sa kolonama slika.
    '''
    sl_kol = np.transpose([tacke['Ap'], tacke['Bp'], tacke['Cp']])

    '''
            Dp = alfap * Ap + betap * Bp + gamap * Cp
    '''
    alfap_betap_gamap = np.linalg.solve(sl_kol, np.transpose(tacke['Dp']))

    '''
            P2 je matrica sa kolonama alfap * Ap, betap * Bp, gamap * Cp
    '''
    P2 = alfap_betap_gamap * sl_kol

    '''
            III korak:
                    P = P2 * P1^-1 je trazena matrica preslikavanja
    '''

    ''' 
            Racunamo P1^-1
    '''
    P1_inv = np.linalg.inv(P1)

    ''' 
            P = P2 * P1^-1
    '''
    P = np.matmul(P2, P1_inv)

    '''     
            Izlaz:
                    3 x 3 matrica P projektivnog preslikavanja ravni
                    koje slika
                    A,  B,  C,  D  u
                    Ap, Bp, Cp, Dp.
    '''
    return P

def matrica_2_9(T, Tp):
    '''
            Ulaz: 
                    tacka i njena slika

            Izlaz: 
                    2 x 9 matrica data u lemi (prezentacija 1, str 31)
    '''
    return [[0, 0, 0, -Tp[2] * T[0], -Tp[2] * T[1], -Tp[2] * T[2], Tp[1] * T[0], Tp[1] * T[1], Tp[1] * T[2]],
            [Tp[2] * T[0], Tp[2] * T[1], Tp[2] * T[2], 0, 0, 0, -Tp[0] * T[0], -Tp[0] * T[1], -Tp[0] * T[2]]]

def DLT(tacke, br_tacaka):
    '''
            Ulaz:
                    homogene koordinate
                    n >= 4 originalnih tacaka A,  B,  C,  D, ..
                    i n njihovih slika  Ap, Bp, Cp, Dp, ..
            Mi cemo proslediti i broj tacaka.
    '''

    '''
            I korak:
                    Za svaku korespodenciju tacaka odredjujemo 2 x 9 matricu.
                    
            II korak:
                    Spajamo ove matrice u jednu matricu A formata 2n x 9.
    '''
    A = np.reshape([matrica_2_9(tacke[chr(ord('A') + x)], tacke[(chr(ord('A') + x) + "p")]) for x in range(0, br_tacaka)],
                   (2 * br_tacaka, 9))

    '''
            III korak:
                    Odredjujemo SVD dekompoziciju matrice A = UDV^T.
    '''
    U, D, V_T = np.linalg.svd(A)

    '''
                Izlaz:
                Resenje je poslednja kolona(vrsta) matrice V^T.
    '''
    P = np.reshape(V_T[-1, :], (3, 3))
    return P

def homogene_u_afine_koordinate(tacke):
    return [(x[0] / x[2], x[1] / x[2]) for x in tacke]

def normalizacija_tacaka(tacke):
    '''
            I korak:
                    Izracunati teziste sistema tacaka (afino).
    '''
    '''
            Sva izracunavanja se vrse sa afinim koordinatama,
            zato tacke iz homogenih 
            prebacujemo u afine koordinate.
    '''
    tacke_a = homogene_u_afine_koordinate(tacke)
    Tx = sum([t[0] for t in tacke_a]) / len(tacke_a)
    Ty = sum([t[1] for t in tacke_a]) / len(tacke_a)

    '''
            II korak:
                    Translirati teziste u koordinatni pocetak.
    '''
    '''
            Matrica translacije G:
    '''
    G = np.array([[1, 0, -Tx],
                  [0, 1, -Ty],
                  [0, 0, 1]])

    '''
            III korak:
                    Skalirati tacke 
                    tako da je prosecna udaljenost tacke od koordinatnog pocetka bude sqrt(2).
    '''
    '''
            Prvo trazimo prosecnu udaljenost od koordinatnog pocetka l.
    '''
    l = sum(np.sqrt((t[0] - Tx) ** 2 + (t[1] - Ty) ** 2) for t in tacke_a) / len(tacke_a)

    '''
            Odredjujemo koeficijent homotetije k:
    '''
    k = np.sqrt(2) / l

    '''
            Na kraju je matrica homotetije S:
    '''
    S = np.array([[k, 0, 0],
                  [0, k, 0],
                  [0, 0, 1]])

    '''
            IV korak:
                    Trazena matrica normalizacije je T = S*G.
    '''
    T = np.matmul(S, G)
    return T

def modifikovani_DLT(tacke, br_tacaka):
    '''
            Ulaz:
                    homogene koordinate
                    n >= 4 originalnih tacaka A,  B,  C,  D, ..
                    i n njihovih slika  Ap, Bp, Cp, Dp, ..
            Mi cemo proslediti i broj tacaka.
    '''

    '''
            I korak:
                    Normalizacija originalnih tacaka transformacijom T.
    '''
    '''
            II korak:
                    Normalizacija slika transformacijom Tp.
    '''
    normalizovane = {}
    T = normalizacija_tacaka([tacke[chr(ord('A') + x)] for x in range(0, br_tacaka)])
    Tp = normalizacija_tacaka([tacke[(chr(ord('A') + x) + "p")] for x in range(0, br_tacaka)])

    normalizovane = dict(zip([chr(ord('A') + x) for x in range(0, br_tacaka)] + [(chr(ord('A') + x) + "p") for x in range(0, br_tacaka)], \
                             [tuple(np.matmul(T, tacke[chr(ord('A') + x)]).flatten()) for x in range(0, br_tacaka)] +
                             [tuple(np.matmul(Tp, tacke[(chr(ord('A') + x) + "p")]).flatten()) for x in range(0, br_tacaka)]))

    '''
           III korak:
                   DLT algoritmom odrediti matricu transformacije Pn
                   koja slika normalizovane tacke
                   u normalizovane slike.
   '''
    Pn = DLT(normalizovane, br_tacaka)

    '''
            IV korak:
                    Trazena matrica transformacije je P = Tp^-1 * Pn * T.
    '''
    Tp_inv = np.linalg.inv(Tp)
    P = Tp_inv.dot(Pn).dot(T)

    '''
            Izlaz:
                    3 x 3 matrica P projektivnog preslikavanja.
    '''
    return P


def main():
    '''
            Prvo uzimamo listu tacaka iz zeljenog primera.
    '''
    t = ucitaj_tacke("data/primer1.txt")
    br_tacaka = int(len(t) / 2)

    '''
            Pravimo mapu oblika 
            (tacka, koordinate) 
            odnosno (A, (-3, -1, 1))
            i (slika, koordinate)
            odnosno (Ap, (-2, -1, 1)).
    '''
    tacke = dict(zip([chr(ord('A') + x) for x in range(0, br_tacaka)] + [(chr(ord('A') + x) + "p") for x in range(0, br_tacaka)], t))
    print("\nTacke i njihove koordinate: ")
    for tacka, koordinate in tacke.items():
        print(tacka, ": ", koordinate)

    '''
            Za drugi primer, za poredjenje sa matricom sa pdf-a
            moze se izvrsiti sledeci kod:
                print(np.reshape(np.array(naivni(tacke).flatten()) * 0.142857238, (3, 3)))
    '''
    print("\nNaivni: \n", naivni(tacke).round(7))

    print("\nDLT: \n", DLT(tacke, br_tacaka).round(5))

    print("\nModifikovani DLT:\n", modifikovani_DLT(tacke, br_tacaka).round(5))

if __name__ == '__main__':
	main()