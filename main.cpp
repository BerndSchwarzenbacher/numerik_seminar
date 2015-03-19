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
    vecx = 2.0;
    AutoVector vecy = mymat.CreateVector();
    vecy = 3.0;
    AutoVector vecy2 = mymat.CreateVector();
    vecy2 = 3.0;
    AutoVector vecy3 = mymat.CreateVector();
    vecy3 = 3.0;
    AutoVector vecy4 = mymat.CreateVector();
    vecy4 = 3.0;
    AutoVector vecy5 = mymat.CreateVector();
    vecy5 = 3.0;
    AutoVector vecy6 = mymat.CreateVector();
    vecy5 = 3.0;

    mymat.MultAdd1(1.0, vecx, vecy);
    mymat.MultAdd2(1.0, vecx, vecy2);
    mymat.MultAdd3(1.0, vecx, vecy3);
    mymat.TranMultAdd1(1.0, vecx, vecy4);
    mymat.TranMultAdd2(1.0, vecx, vecy5);
    mymat.TranMultAdd3(1.0, vecx, vecy6);

    diff1 = vecy - vecy2;
    cout << "LNORM2(vecy - vecy2)" << L2Norm(diff1) << endl;
}

