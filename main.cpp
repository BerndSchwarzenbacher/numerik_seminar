// /micNfs/netgen/bin/ngscxx test_matrix_mul.cpp -lngstd -lsolve -lngcomp -lngfem -linterface -limf -lsvml -lirng -lintlc -lmkl_rt -L/micNfs/intel/mkl/lib/mic/
#include <solve.hpp>
#include "mysparsematrix.hpp"

using namespace ngsolve;

void LoadMatrix(string, SparseMatrix<double> &);

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
    AutoVector vecy1 = mymat.CreateVector();
    vecy1 = 3.0;
    AutoVector vecy2 = mymat.CreateVector();
    vecy2 = 3.0;
    AutoVector vecy3 = mymat.CreateVector();
    vecy3 = 3.0;
    AutoVector vecy4 = mymat.CreateVector();
    vecy4 = 3.0;
    AutoVector vecy5 = mymat.CreateVector();
    vecy5 = 3.0;
    AutoVector vecy6 = mymat.CreateVector();
    vecy6 = 3.0;
    AutoVector vecy7 = mymat.CreateVector();
    vecy7 = 3.0;
    AutoVector vecy8 = mymat.CreateVector();
    vecy8 = 3.0;

    mymat.MultAdd1(2.0, vecx, vecy);
    mymat.TranMultAdd1(2.0, vecx, vecy1);

    mymat.MultAdd2(2.0, vecx, vecy2);
    mymat.MultAdd3(2.0, vecx, vecy3);
    mymat.TranMultAdd2(2.0, vecx, vecy5);
    mymat.TranMultAdd3(2.0, vecx, vecy6);
    mymat.TranMultAdd4(2.0, vecx, vecy7);
    mymat.TranMultAdd5(2.0, vecx, vecy8);

    AutoVector diff2= mymat.CreateVector();
    AutoVector diff3= mymat.CreateVector();
    AutoVector diff5= mymat.CreateVector();
    AutoVector diff6= mymat.CreateVector();
    AutoVector diff7= mymat.CreateVector();
    AutoVector diff8= mymat.CreateVector();
    diff2 = vecy - vecy2;
    diff3 = vecy - vecy3;
    diff5 = vecy1 - vecy5;
    diff6 = vecy1 - vecy6;
    diff7 = vecy1 - vecy7;
    diff8 = vecy1 - vecy8;

    cout << "LNORM2(vecy - vecy2) -sequential  " << L2Norm(diff2) << endl;
    cout << "LNORM2(vecy - vecy3) -parallel  " << L2Norm(diff3) << endl;
    cout << "LNORM2(vecy1 - vecy5) -Trans parallel for static  " << L2Norm(diff5) << endl;
    cout << "LNORM2(vecy1 - vecy6) -Trans balancing  " << L2Norm(diff6) << endl;
    cout << "LNORM2(vecy1 - vecy7) -Trans coloring  " << L2Norm(diff7) << endl;
    cout << "LNORM2(vecy1 - vecy8) -Trans parallel for dynamic 100  " << L2Norm(diff8) << endl;
}

void LoadMatrix(string filename, SparseMatrix<double> & mat) {
    BinaryInArchive ar(filename.c_str());
    Timer bl("binary - load");
    bl.Start();
    try {
        mat.DoArchive(ar);
    }
    catch (std::bad_alloc& ba) {
        std::cerr << "bad_alloc caught: " << ba.what() << '\n';
    }
    bl.Stop();
    cout << mat.Size() << endl;
    cout << mat.Width() << endl;
    cout << mat.NZE() << endl;
}


