#ifndef MVN_HPP__
#define MVN_HPP__

#include<vector>
#include<string>

#include<gsl/gsl_matrix.h>
#include<cmath>

#include<cassert>

namespace mvn {

	struct VecCol { // A column vector
		private:
		gsl_vector * m_a;
		public:
		explicit VecCol(const size_t D);
		inline size_t size() const {return m_a->size;}
		inline int    len () const {return m_a->size;}
		inline gsl_vector * get() { return m_a; }
		inline const gsl_vector * get() const { return m_a; }
		       bool               isnull() const { return gsl_vector_isnull(get()); }
		void set (size_t c, double d) {
			return gsl_vector_set(m_a, c, d);
		}

		~VecCol() { // Rule of Five
			if(m_a)
				gsl_vector_free(m_a);
			// otherwise, it was moved from
		}
		VecCol(const VecCol &in); // Rule Of Five
		VecCol& operator=(const VecCol &in); // Rule Of Five
		VecCol(VecCol &&in); // Rule Of Five
		VecCol& operator=(VecCol &&in); // Rule Of Five

		VecCol& operator+=(const VecCol &in);
		VecCol& operator-=(const VecCol &in);
		VecCol& operator*=(const long double scale);
		VecCol& operator/=(const long double scale);

		double operator() (size_t r) const {
			return gsl_vector_get(m_a, r);
		}
	};
	static_assert(sizeof(VecCol) == sizeof( gsl_vector *), "extra elements in VecCol?");
	struct SquareMatrix { // A square matrix
		gsl_matrix * a;
		struct IdentityMatrixRequested {};
		explicit SquareMatrix(const size_t D);
		explicit SquareMatrix(const size_t D, IdentityMatrixRequested);
		~SquareMatrix(); // Rule Of Five
		void set_identity();
		size_t size() const { assert(a); return a->size1;}
		int    len () const { assert(a); return a->size1;}
		gsl_matrix * get() {
			assert(a); return a;
		}
		const gsl_matrix * get() const {  assert(a); return a; }
		       bool               isnull() const { return gsl_matrix_isnull(get()); }

		SquareMatrix(const SquareMatrix &in); // Rule Of Five
		SquareMatrix& operator=(const SquareMatrix &in); // Rule Of Five
		SquareMatrix(SquareMatrix &&in); // C++11 rule of Five
		SquareMatrix& operator=(SquareMatrix &&in); // C++11 rule of Five
		SquareMatrix& operator+=(const SquareMatrix &in);
		SquareMatrix& operator-=(const SquareMatrix &in);
		SquareMatrix& operator*=(long double);
		SquareMatrix& operator/=(long double);
		double operator() (size_t r, size_t c) const {
			return gsl_matrix_get(a, r, c);
		}
		void set (size_t r, size_t c, double d) {
			return gsl_matrix_set(a, r, c, d);
		}
	};
	static_assert(sizeof(SquareMatrix) == sizeof(SquareMatrix :: a), "extra elements in SquareMatrix?");
	struct Matrix { // A matrix, need not be square
		gsl_matrix * a;
		explicit Matrix(const size_t rows, const size_t columns);
		~Matrix(); // Rule Of Five
		size_t size1() const { assert(a); return a->size1;}
		size_t size2() const { assert(a); return a->size2;}
		gsl_matrix * get()             { assert(a); return a; }
		const gsl_matrix * get() const {  assert(a); return a; }

		Matrix(const Matrix &in); // Rule Of Three/Five
		Matrix& operator=(const Matrix &in); // Rule Of Three/Five
		Matrix(Matrix &&in); // C++11 rule of Five
		Matrix& operator=(Matrix &&in); // C++11 rule of Five
		//Matrix(Matrix &&in); // C++11 rule of Five
		//Matrix& operator=(Matrix &&in); // C++11 rule of Five
		//Matrix& operator+=(const Matrix &in);
		//Matrix& operator-=(const Matrix &in);
		//Matrix& operator*=(long double);
		//Matrix& operator/=(long double);
		double operator() (size_t r, size_t c) const {
			return gsl_matrix_get(a, r, c);
		}
		void set (size_t r, size_t c, double d) {
			return gsl_matrix_set(a, r, c, d);
		}
	};
	typedef std :: vector< VecCol* > mvn_data_t;
	void load_mvn_data (mvn_data_t &data, const std :: string & dataFileName);

	VecCol operator-(const VecCol &lhs, const VecCol &rhs);
	VecCol operator+(const VecCol &lhs,       VecCol  rhs);
	VecCol operator*(long double lhs, VecCol rhs);
	VecCol operator/(VecCol rhs, long double lhs);
	VecCol multiply_rowvec_by_matrix_giving_rowvec(const VecCol &rowvec, const SquareMatrix &matrix);
	VecCol multiply_matrix_by_colvec_giving_colvec(const SquareMatrix &matrix, const VecCol &colvec);
	VecCol multiply_matrix_by_colvec_giving_colvec(const       Matrix &matrix, const VecCol &colvec);
	VecCol multiply_rowvec_by_colvec_giving_scalar(const VecCol &lhs, const VecCol &rhs);
	SquareMatrix multiply_colvec_by_rowvec_giving_matrix(const VecCol &lhs, const VecCol &rhs);
	SquareMatrix operator*(const SquareMatrix &lhs, const SquareMatrix &rhs);
	SquareMatrix operator+(const SquareMatrix &lhs, const SquareMatrix &rhs);
	SquareMatrix operator-(const SquareMatrix &lhs, const SquareMatrix &rhs);
	SquareMatrix operator*(long double, SquareMatrix rhs);
	inline
	SquareMatrix operator*(SquareMatrix rhs, long double x) {
		return x * std::move(rhs);
	}
	SquareMatrix operator/(SquareMatrix lhs, long double);
	SquareMatrix		invert_a_matrix_impl(SquareMatrix copy_with_which_to_do_LU, const char * file, const size_t line) ;
#define invert_a_matrix(copy) invert_a_matrix_impl(copy, __FILE__, __LINE__)
	VecCol      		solve_a_matrix(SquareMatrix copy_with_which_to_do_LU, VecCol const &b) ;
	SquareMatrix		cholesky_upper(SquareMatrix copy_with_which_to_do_LL) ;

	Matrix operator*(const Matrix &lhs, const SquareMatrix &rhs);


bool veccol_is_close_to_zero(const VecCol &x);
bool matrix_is_close_to_zero(const SquareMatrix &x);
double norm_2(const SquareMatrix &x);
double norm_2(const VecCol       &x);

inline
bool    operator==(SquareMatrix const& lhs, SquareMatrix const &rhs) {
    assert(lhs.size() == rhs.size());
    return gsl_matrix_equal(lhs.get(), rhs.get());
}
inline
bool    operator==(VecCol const& lhs, VecCol const &rhs) {
    assert(lhs.size() == rhs.size());
    return gsl_vector_equal(lhs.get(), rhs.get());
}

long double
calculate_lndeterminant(mvn:: SquareMatrix copy_with_which_to_do_LU);

mvn:: VecCol make_VecCol(std::vector<double> const & in);

#ifdef USE_GPERF_CPUPROFILE
	void dump_final_destruct_Square_counters();
#endif
} // namespace mvn

namespace std {
	std :: ostream& operator << (std :: ostream & o, const mvn :: VecCol &);
	std :: ostream& operator << (std :: ostream & o, const mvn :: SquareMatrix &);
	std :: ostream& operator << (std :: ostream & o, const mvn :: Matrix &);
}

#endif
