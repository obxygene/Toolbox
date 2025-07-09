
import kwant

def square_Hamiltonian(a = 1, t = 1, t2 = 0.5, X = 5, Y = 5):
    syst = kwant.Builder()
    lat = kwant.lattice.honeycomb(a, norbs=1, name=['a', 'b'])
    def square(pos):
        x, y = pos
        return (abs(x) < X) & (abs(y) < Y)
    syst[lat.shape(square,(0,0))] = 0.
    syst[lat.neighbors()] = t
    # add second neighbours hoppings
    syst[lat.a.neighbors()] = 1j * t2
    syst[lat.b.neighbors()] = -1j * t2
    syst.eradicate_dangling()
    kwant.plot(syst)
    fsyst = syst.finalized()
    hamiltonian_Matrix = fsyst.hamiltonian_submatrix(sparse=True)

    return hamiltonian_Matrix.toarray()