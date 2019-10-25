#ifndef matrix_h
#define matrix_h

#include <stdexcept>
#include <vector>
#include <function>

// undefine to disable range checking
#define RANGE_CHECK

class overdetermined : public std::domain_error
{
   public:
     overdetermined() 
        : std::domain_error("solution is over-determined")
     {}
};

class underdetermined : public std::domain_error
{
   public:
      underdetermined() 
         : std::domain_error("solution is under-determined")
      {}
};

template<class T>
class kahn_sum {     // implements Kahn Summation method
   public:
      kahn_sum() : sum(0.0), cor(0.0) {}
      kahn_sum<T>& operator+=( const T& val ) {
         T old_sum = sum;
         T next = val-cor;
         cor = ( (sum += next) - old_sum ) - next;
         return *this;
      }
      kahn_sum<T>& operator-=( const T& val ) {
         T old_sum = sum;
         T next = val+cor;
         cor = ( (sum -= val) - old_sum ) + next;
         return *this;
      }
      operator T&() { return sum; }
   private:
      T sum;  // running sum
      T cor;  // correction term
};

template<class T>
class matrix {
   private:
      std::vector<T> elements; // array of elements
   public:
      const unsigned rows;  // number of rows
      const unsigned cols;  // number of columns
      
   protected:
      // range check function for matrix access
      void range_check( unsigned i, unsigned j ) const;
   public:
      T& operator()( unsigned i, unsigned j ) {
         #ifdef RANGE_CHECK
         range_check(i,j);
         #endif
         return elements[i*cols+j];
      }
      const T& operator()( unsigned i, unsigned j ) const {
         #ifdef RANGE_CHECK
         range_check(i,j);
         #endif
         return elements[i*cols+j];
      }
      const T& element(unsigned i, unsigned j) const {
         #ifdef RANGE_CHECK
         range_check(i,j);
         #endif
         return elements[i*cols+j];
      }
      T& element(unsigned i, unsigned j) {
         #ifdef RANGE_CHECK
         range_check(i,j);
         #endif
         return elements[i*cols+j];
      }
   public:
      // constructors
      matrix( unsigned rows, unsigned cols, const T* elements = 0 );
      matrix( const matrix<T>& );
      // destructor
      ~matrix();

      // assignment
      matrix<T>& operator=( const matrix<T>& );

      // comparison
      bool operator==( const matrix<T>& ) const;
      bool iszero() const;
      bool operator!() const {
        return iszero();
      }

      // scalar multiplication/division
      matrix<T>& operator*=( const T& a );
      matrix<T> operator*( const T& a ) const {
         return matrix<T>(*this).operator*=(a);
      }
      matrix<T>& operator/=( const T& a );
      matrix<T> operator/( const T& a ) {
        return matrix<T>(*this).operator/=(a);
      }
      matrix<T> operator-() const;
      matrix<T> operator+() const;

      // addition/subtraction
      matrix<T>& operator+=( const matrix<T>& );
      matrix<T>& operator-=( const matrix<T>& );
      matrix<T> operator+( const matrix<T>& M ) const {
         return matrix<T>(*this).operator+=(M);
      }
      matrix<T> operator-( const matrix<T>& M ) const {
         return matrix<T>(*this).operator-=(M);
      }

      // matrix multiplication
      matrix<T> operator*( const matrix<T>& ) const;
      matrix<T>& operator*=( const matrix<T>& M ) {
         return *this = *this * M;
      }
      // matrix division
      matrix<T> leftdiv( const matrix<T>& ) const;
      matrix<T> rightdiv( const matrix<T>& D ) const {
         return transpose().leftdiv(D.transpose()).transpose();
      }
      matrix<T> operator/( const matrix<T>& D ) const {
         return rightdiv(D);
      }
      matrix<T>& operator/=( const matrix<T>& M ) {
         return *this = *this/M;
      }

      // determinants
      matrix<T> minor( unsigned i, unsigned j ) const;
      T det() const;
      T minor_det( unsigned i, unsigned j ) const;

      // these member functions are only valid for squares
      matrix<T> inverse() const;
      matrix<T> pow( int exp ) const;
      matrix<T> identity() const;
      bool isidentity() const;

      // vector operations
      matrix<T> getrow( unsigned j ) const;
      matrix<T> getcol( unsigned i ) const;
      matrix<T>& setcol( unsigned j, const matrix<T>& C );
      matrix<T>& setrow( unsigned i, const matrix<T>& R );
      matrix<T> delrow( unsigned i ) const;
      matrix<T> delcol( unsigned j ) const;

      matrix<T> transpose() const;
      matrix<T> operator~() const {
         return transpose();
      }
};

template<class T>
matrix<T>::matrix( unsigned rows, unsigned cols, const T* elements = 0 )
  : rows(rows), cols(cols), elements(rows*cols,T(0.0))
{
   if( rows == 0 || cols == 0 )
      throw std::range_error("attempt to create a degenerate matrix");
   // initialze from array
   if(elements)
      for(unsigned i=0;i<rows*cols;i++)
         this->elements[i] = elements[i];
};

template<class T>
matrix<T>::matrix( const matrix<T>& cp )
  : rows(cp.rows), cols(cp.cols), elements(cp.elements)
{
}

template<class T>
matrix<T>::~matrix()
{
}

template<class T>
matrix<T>& matrix<T>::operator=( const matrix<T>& cp )
{
   if(cp.rows != rows && cp.cols != cols )
      throw std::domain_error("matrix op= not of same order");
   for(unsigned i=0;i<rows*cols;i++)
      elements[i] = cp.elements[i];
   return *this;
}
template<class T>
void matrix<T>::range_check( unsigned i, unsigned j ) const
{
   if( rows <= i )
      throw std::range_error("matrix access row out of range");
   if( cols <= j )
      throw std::range_error("matrix access col out of range");
}

//not shown: operator==(const matrix<T>& A) const, iszero() const,
//operator*=(const T& a), operator/=(const T& a), operator-() const,
//operator+() const, operator*(const T& a, const matrix<T>& M),
//operator+=(const matrix<T>& M), operator-=(const matrix<T>& M),
//minor(unsigned i, unsigned j) const,
//minor_det(unsigned i, unsigned j) const, det() const -- mb
//operator*( const matrix<T>& B) const ... mb

template<class T>
matrix<T> matrix<T>::leftdiv( const matrix<T>& D ) const
{
   const matrix<T>& N = *this;

   if( N.rows != D.rows )
      throw std::domain_error("matrix divide: incompatible orders");

   matrix<T> Q(D.cols,N.cols);  // quotient matrix

   if( N.cols > 1 ) {
       // break this down for each column of the numerator
       for(unsigned j=0;j<Q.cols;j++)
          Q.setcol( j, N.getcol(j).leftdiv(D) );  // note: recursive
       return Q;
   }

   // from here on, N.col == 1

   if( D.rows < D.cols )
      throw underdetermined();

   if( D.rows > D.cols ) {
      bool solution = false;
      for(unsigned i=0;i<D.rows;i++) {
         matrix<T> D2 = D.delrow(i);  // delete a row from the matrix
         matrix<T> N2 = N.delrow(i);
         matrix<T> Q2(Q);  
         try {
            Q2 = N2.leftdiv(D2);
         } 
         catch( underdetermined x ) {
            continue;  // try again with next row
         }
         if( !solution ) {
            // this is our possible solution
            solution = true;
            Q = Q2;
         } else {
            // do the solutions agree?
            if( Q != Q2 )
               throw overdetermined();
         }
      }
      if( !solution )
         throw underdetermined();
      return Q;
   }

   // D.rows == D.cols && N.cols == 1
   // use Kramer's Rule
   //
   // D is a square matrix of order N x N
   // N is a matrix of order N x 1

   const T T0(0.0); // additive identity

   if( D.cols<=3 ) {
      T ddet = D.det();
      if( ddet == T0 )
         throw underdetermined();
      for(unsigned j=0;j<D.cols;j++) {
         matrix<T> A(D); // make a copy of the D matrix
         // replace column with numerator vector
         A.setcol(j,N);
         Q(j,0) = A.det()/ddet;
      }
   } else {
      // this method optimizes the determinant calculations
      // by saving a cofactors used in calculating the
      // denominator determinant.
      kahn_sum<T> sum;
      vector<T> cofactors(D.cols);  // save cofactors
      for(unsigned j=0;j<D.cols;j++) {
         T c = D.minor_det(0,j);
         cofactors[j] = c;
         T a = D(0,j);
         if( a != T0 ) {
            a *= c;
            if(j%2)
               sum -= a;
            else
               sum += a;
         }
      }
      T ddet = sum;
      if( ddet == T0 )
         throw underdetermined();
      for(unsigned j = 0;j<D.cols;j++) {
         matrix<T> A(D);
         A.setcol(j,N);
         kahn_sum<T> ndet;
         for(unsigned k=0;k<D.cols;k++) {
            T a = A(0,k);
            if( a != T0 ) {
               if(k==j)
                  a *= cofactors[k];  // use previously calculated
                                      // cofactor
               else
                  a *= A.minor_det(0,k); // calculate minor's determinant
               if(k%2)
                  ndet -= a;
               else
                  ndet += a;
            }
         }
         Q(j,0) = T(ndet)/ddet;
      }
   }
   return Q;
}

//not shown: inverse() const, getrow(unsigned i) const,
//getcol(unsigned j) const, delcol(unsigned j) const,
//delrow(unsigned i) const, setcol(unsigned j, const matrix<T>& C)
//setrow(unsigned i, const matrix<T>& R), identity() const,
//isidentity() const, pow(int exp) const,
//pow(const matrix<T>& M, int exp) --mb
// ...
#endif
