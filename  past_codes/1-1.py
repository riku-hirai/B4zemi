
C = 0.5
nmax = 3

def init(nmax):
    for n in range (1, nmax+1):
        x = 0.1*(n-1.0)
        A = 1 + 2.2 * (x - 1.5) ** 2
        rho = 1 - 0.3146 * x
        T = 1 - 0.2314 * x
        V = (0.1 + 1.09 * x) * T ** 0.5
        print(A,V)

init(nmax)