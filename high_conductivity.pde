geometry = square_conductivity.in2d
mesh = square_conductivity_fine.vol

constant numthreads = 1

define coefficient lam
1, 1000

define coefficient penalty
1e5, 0, 0

define coefficient coef_source
1,

define fespace v -type=h1ho -order=4 -dirichlet=[1]

define gridfunction u -fespace=v -nested

define bilinearform a -fespace=v -symmetric -printelmat -fespace=v
laplace lam
# robin penalty

define linearform f -fespace=v
source coef_source

#define preconditioner c -type=myamg -bilinearform=a
#define preconditioner c -type=amg -bilinearform=a -test
#define preconditioner c -type=local -bilinearform=a -test
define preconditioner c -type=bddc -bilinearform=a

numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=u  -preconditioner=c -maxsteps=1000

numproc visualization npvis -scalarfunction=u -subdivision=1 -nolineartexture

