#ifndef MYSPARSEMATRIX_HPP___
#define MYSPARSEMATRIX_HPP___
#include <solve.hpp>
#include <map>


using namespace ngsolve;

class MySparseMatrix : public SparseMatrix<double,double,double>
{
private:
  Table<int> coloring_;

public:
  MySparseMatrix(const SparseMatrixTM<double> & mat)
    : SparseMatrix<double,double,double>(mat), SparseMatrixTM<double>(mat)
  { }

  void Coloring ()
  {
    static Timer timer("SparseMatrix::Coloring");
    RegionTimer reg (timer);
    timer.AddFlops(this->nze);

    int height = this->Height();
    int width = this->Width();

    Array<int> row_color(height);
    row_color = -1;

    int maxcolor = 0;
    int basecol = 0;

    Array<unsigned int> mask(width);

    int found = 0;

    do
    {
      mask = 0;

      for (int row = 0; row < height; ++row)
      {
        if (row_color[row] >= 0) continue;

        int first = firsti [row];
        int last  = firsti [row+1];

        unsigned check = 0;
        for (int i = first; i < last; ++i)
          check |= mask[colnr[i]];

        if (check != UINT_MAX) // 0xFFFFFFFF)
        {
          found++;
          unsigned checkbit = 1;
          int color = basecol;
          while (check & checkbit)
          {
            color++;
            checkbit *= 2;
          }

          row_color[row] = color;
          if (color > maxcolor) maxcolor = color;

          for (int i = first; i < last; ++i)
            mask[colnr[i]] |= checkbit;
        }
      }

      basecol += 8*sizeof(unsigned int); // 32;
    }
    while (found < height);

    Array<int> cntcol(maxcolor+1);
    cntcol = 0;
    for (int row = 0; row < height; ++row)
      ++cntcol[row_color[row]];

    coloring_ = Table<int>(cntcol);

    cntcol = 0;
    for (int row = 0; row < height; ++row)
      coloring_[row_color[row]][cntcol[row_color[row]]++] = row;

    std::cout << "needed " << maxcolor+1 << " colors" << std::endl;
  }

  void MultAdd1 (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer timer("SparseMatrix::MultAdd-Original");
    RegionTimer reg (timer);
    timer.AddFlops(this->nze);
    FlatVector<double> fx = x.FV<double> ();
    FlatVector<double> fy = y.FV<double> ();
    for (int i = 0; i < this->Height(); ++i)
      fy(i) += s * RowTimesVector (i, fx);
  }

  ///////////////////////////////////////////////////////////////////////
  void MultAdd2 (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer timer("SparseMatrix::MultAdd-Sequential");
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
        fy(i) += s * data[j] * fx(colnr[j]);
      }
    }
  }

  /////////////////////////////////////////////////////////////////////
  void MultAdd3 (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer timer("SparseMatrix::MultAdd-Parallel");
    RegionTimer reg (timer);
    timer.AddFlops(this->nze);
    FlatVector<double> fx = x.FV<double> ();
    FlatVector<double> fy = y.FV<double> ();

#pragma omp parallel for
    for (int i = 0; i < this->Height(); ++i)
    {
      int first = firsti [i];
      int last  = firsti [i+1];

      for (int j = first; j < last; ++j)
      {
        fy(i) += s * data[j] * fx(colnr[j]);
      }
    }
  }

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

  ////////////////////////////////////////////////////////////////////
  void TranMultAdd3 (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer timer("SparseMatrix::TranMultAdd-Balancing");
    RegionTimer reg (timer);
    timer.AddFlops(this->nze);
    FlatVector<double> fx = x.FV<double> ();
    FlatVector<double> fy = y.FV<double> ();

    int height = this->Height();

    Array<int> thread_seperation;
#pragma omp parallel //shared(thread_seperation)
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

  void TranMultAdd4 (double s, const BaseVector & x, BaseVector & y)
  {
    static Timer timer("SparseMatrix::TranMultAdd-Coloring");
    RegionTimer reg (timer);
    timer.AddFlops(this->nze);
    FlatVector<double> fx = x.FV<double> ();
    FlatVector<double> fy = y.FV<double> ();

    Coloring();

#pragma omp parallel
    {
      for (auto color : coloring_)
      {
#pragma omp for
        for (int i = 0; i < color.Size(); ++i)
        {
          int first = firsti [color[i]];
          int last  = firsti [color[i]+1];

          for (int j = first; j < last; ++j)
          {
            fy(colnr[j]) += s * data[j] * fx(color[i]);
          }
        }
      }

    }
  }

	////////////////////////////////////////////////////////////////////
	void TranMultAdd5 (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer timer("SparseMatrix::TranMultAdd-Parallel-For-Dyn");
    RegionTimer reg (timer);
    timer.AddFlops(this->nze);
    FlatVector<double> fx = x.FV<double> ();
    FlatVector<double> fy = y.FV<double> ();

#pragma omp parallel
    {
#pragma omp for schedule (dynamic, 100)
      for (int i = 0; i < this->Height(); ++i)
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
