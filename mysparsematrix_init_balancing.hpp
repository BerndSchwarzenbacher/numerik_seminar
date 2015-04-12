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
    static Timer timer("SparseMatrix::TranMultAdd-Balancing");
    RegionTimer reg (timer);
    timer.AddFlops(this->nze);
    FlatVector<double> fx = x.FV<double> ();
    FlatVector<double> fy = y.FV<double> ();

    int height = this->Height();

    Array<int> thread_seperation;
#pragma omp parallel
    {
      int num_threads = omp_get_num_threads();

#pragma omp single
      {
        thread_seperation = Array<int>(num_threads+1);
        int seperation_step = ceil(this->nze / num_threads);

        int thread_i = 1;
        for (int row = 0; row < height; ++row)
        {
          if (firsti[row] >= thread_i * seperation_step)
          {
            thread_seperation[thread_i] = row;
            ++thread_i;
          }
        }
        thread_seperation[0] = 0;
        thread_seperation[num_threads] = height;
      }

#pragma omp for
      for (int thread_i = 1; thread_i <= num_threads; ++thread_i)
      {
        for (int i = thread_seperation[thread_i-1];
            i < thread_seperation[thread_i]; ++i)
        {
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
  }

  void DoArchive (Archive & ar)
  {
    ar & this->size;
    ar & this->width;
    ar & this->nze;

    int * colnr = new int[nze];
    double * data = new double[nze];

    int seperation_step = 0;
#pragma omp parallel
    {
      int num_threads = omp_get_num_threads();

#pragma omp single
      {
        seperation_step = ceil(this->nze / num_threads);
      }

#pragma omp for
      for (int thread_i = 0; thread_i < num_threads; ++thread_i)
      {
        for (int i = thread_i * seperation_step;
             i < (thread_i + 1) * seperation_step && i < nze;
             ++i)
        {
          colnr[i] = 0;
          data[i] = 0;
        }
      }
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
