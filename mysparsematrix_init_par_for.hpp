#ifndef MYSPARSEMATRIX_HPP___
#define MYSPARSEMATRIX_HPP___
#include <solve.hpp>
#include <map>


using namespace ngsolve;

class MySparseMatrix : public SparseMatrix<double,double,double>
{
public:
  MySparseMatrix(const SparseMatrixTM<double> & mat)
    : SparseMatrix<double,double,double>(mat), SparseMatrixTM<double>(mat)
  { }

	////////////////////////////////////////////////////////////////////
	void TranMultAdd1 (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer timer("SparseMatrix::TranMultAdd-Sequential");
    RegionTimer reg (timer);
    timer.AddFlops(this->nze);
    FlatVector<double> fx = x.FV<double> ();
    FlatVector<double> fy = y.FV<double> ();

    for (int i = 0; i < this->Height(); ++i)
    {
      int first = firsti [i];
      int last  = firsti [i+1];

      for (int j = first; j < last; ++j)
      {
        fy(colnr[j]) += s * data[j] * fx(i);
      }

    }
  }

	////////////////////////////////////////////////////////////////////
	void TranMultAdd2 (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer timer("SparseMatrix::TranMultAdd-Parallel-For");
    RegionTimer reg (timer);
    timer.AddFlops(this->nze);
    FlatVector<double> fx = x.FV<double> ();
    FlatVector<double> fy = y.FV<double> ();

#pragma omp parallel
    {
#pragma omp for
      for (int i = 0; i < this->Height(); ++i) {
        int first = firsti [i];
        int last  = firsti [i+1];

        for (int j = first; j < last; ++j)
        {
#pragma omp atomic
          fy(colnr[j]) += s * data[j] * fx(i);
        }

      }
    }
  }

  void DoArchive (Archive & ar)
  {
    ar & this->size;
    ar & this->width;
    ar & this->nze;

    int * colnr = new int[nze];
    double * data = new double[nze];

#pragma omp parallel for
    for (int i = 0; i < nze; ++i)
    {
      colnr[i] = 0;
      data[i] = 0;
    }

    this->colnr = Array<int>(nze, colnr);
    this->data = Array<double>(nze, data);

    ar & this->firsti;
    ar & this->colnr;
    ar & this->data;
    cout << "sparsemat, doarch, sizeof (firstint) = " << firsti.Size() << endl;
  }

};

#endif // MYSPARSEMATRIX_HPP___
