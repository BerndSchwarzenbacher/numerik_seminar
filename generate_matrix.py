from ngsolve import *
pde = comp.PDE("pde_file_eintragen") 
mesh = pde.Mesh()

v = pde.spaces["fespace_eintragen"]
v.Update (heapsize=1000000)

a = pde.bilinearforms["bilinearform eintragen"]
a.Assemble()

archive = ngstd.Archive("matrix.out",True, True)
mat = a.mat
la.DoArchive(archive,mat)


