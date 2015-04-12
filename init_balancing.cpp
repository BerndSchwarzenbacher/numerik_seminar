#include <solve.hpp>
#include "mysparsematrix_init_balancing.hpp"

using namespace ngsolve;

void LoadMatrix(string, MySparseMatrix &);

int main( int argc, char *argv[] ) {
    if(argc<2) {
        cout << "usage: " << argv[0] << " matfile" << endl;
        return -1;
    }
    char *matfile = argv[1];
    cout << "Load Matrix from file " << matfile << endl;
    SparseMatrix<double> tmat(0,0);
    MySparseMatrix mymat(tmat);
    LoadMatrix(matfile, mymat);



    AutoVector vecx = mymat.CreateVector();
    vecx = 2.0;
    AutoVector vecy1 = mymat.CreateVector();
    vecy1 = 3.0;
    AutoVector vecy5 = mymat.CreateVector();
    vecy5 = 3.0;

    mymat.TranMultAdd1(2.0, vecx, vecy1);
    mymat.TranMultAdd2(2.0, vecx, vecy5);

    AutoVector diff5= mymat.CreateVector();
    diff5 = vecy1 - vecy5;

    cout << "LNORM2(vecy1 - vecy5) - Trans balancing "
         << L2Norm(diff5) << endl;
}

void LoadMatrix(string filename, MySparseMatrix & mat) {
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

