#ifndef MYSPARSEMATRIX_HPP___
#define MYSPARSEMATRIX_HPP___
#include <solve.hpp>
#include <map>


using namespace ngsolve;

class MySparseMatrix : public SparseMatrix<double,double,double> {
public:
    MySparseMatrix(const SparseMatrixTM<double> & mat) : SparseMatrix<double,double,double>(mat), SparseMatrixTM<double>(mat) {
//         your own constructor...
    }

        void MultAdd1 (double s, const BaseVector & x, BaseVector & y) const {
            static Timer timer("SparseMatrix::MultAdd1");
            RegionTimer reg (timer);
            FlatVector<double> fx = x.FV<double> ();
            FlatVector<double> fy = y.FV<double> ();
            for (int i = 0; i < this->Height(); ++i)
              fy(i) += s * RowTimesVector (i, fx);
        }

        ///////////////////////////////////////////////////////////////////////
        void MultAdd2 (double s, const BaseVector & x, BaseVector & y) const {
            static Timer timer("SparseMatrix::MultAdd2");
            RegionTimer reg (timer);
            FlatVector<double> fx = x.FV<double> ();
            FlatVector<double> fy = y.FV<double> ();

//#pragma omp parallel for
            for (int i = 0; i < this->Height(); ++i)
            {
              int first = firsti [i];
              int last  = firsti [i+1];

              for (int j = first; j < last; ++j)
              {
                fy(i) += data[j] * fx(colnr[i]);
              }

              fy(i) *= s;
            }

            // your special implementation of MultAdd...
        }
};

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

#endif // MYSPARSEMATRIX_HPP___
