// /micNfs/netgen/bin/ngscxx test_matrix_mul.cpp -lngstd -lsolve -lngcomp -lngfem -linterface -limf -lsvml -lirng -lintlc -lmkl_rt -L/micNfs/intel/mkl/lib/mic/
#include <solve.hpp>
#include "mysparsematrix.hpp"

using namespace ngsolve;

int main( int argc, char *argv[] ) {
    if(argc<2) {
        cout << "usage: " << argv[0] << " matfile" << endl;
        return -1;
    }
    char *matfile = argv[1];
    cout << "Load Matrix from file " << matfile << endl;
    SparseMatrix<double> tmat(0,0);
    LoadMatrix(matfile, tmat);
    MySparseMatrix mymat(tmat);

    AutoVector vecx = mymat.CreateVector();
    AutoVector vecy = mymat.CreateVector();

    mymat.MultAdd(1.0, vecx, vecy);
}
