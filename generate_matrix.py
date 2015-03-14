from ngsolve import *
pde = comp.PDE("high_conductivity.pde")
mesh = pde.Mesh()

v = pde.spaces["v"]
v.Update (heapsize=1000000)

a = pde.bilinearforms["a"]
a.Assemble()

archive = ngstd.Archive("matrix.out",True, True)
mat = a.mat
la.DoArchive(archive,mat)

