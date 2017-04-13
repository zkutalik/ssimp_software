#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <limits>
#include <iostream>

#include "mvn.hh"

#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>

#include "PP.hh"

using namespace std;
namespace mvn {
	void load_mvn_data (mvn_data_t &out_data, const std :: string & dataFileName) {
		// First load the data into a vector of vectors, then convert it to a vector< VecCol* >
		vector< vector<double> > data;
		assert(data.empty());
		ifstream datastream(dataFileName);
		string line;
		while(getline(datastream, line)) {
			data.push_back( vector<double>() ); // start an empty vector to store this data point
			istringstream fieldstream(line);
			double field;
			while(fieldstream >> field) {
				data.back().push_back(field);
			}
			assert(fieldstream.eof());
			assert(!fieldstream.bad());
		}
		assert(datastream.eof());
		assert(!datastream.bad());

		assert(!data.empty());
		const size_t D = data.front().size();

		PP(data.size(), D);
		for( auto& point : data ) {
			assert(D == point.size());
		}

		assert(out_data.empty());
		for( auto& point : data ) {
			VecCol *new_point = new VecCol(D);
			for (size_t d = 0; d<D; ++d) {
				gsl_vector_set(new_point->get(),d, point.at(d));
			}
			out_data.push_back(new_point);
		}
		assert(out_data.size() == data.size());
	}
	VecCol :: VecCol(const size_t D) {
		m_a = gsl_vector_alloc(D);
		gsl_vector_set_zero(m_a);
	}
				VecCol :: VecCol(const VecCol &in) {
					m_a = gsl_vector_alloc(in.size());
					gsl_vector_memcpy(m_a, in.m_a);
	}
	VecCol&			VecCol :: operator=(const VecCol &in) {
				if(this->size() == in.size()) {
					if(this != &in) {
						gsl_vector_memcpy(m_a, in.get());
					}
				} else {
					// size mismatch. Must be funky with swap
					VecCol in_copy(in);
					swap(this->m_a, in_copy.m_a);
				}
				return *this;
	}
				VecCol :: VecCol(VecCol &&in) {
					m_a = in.m_a;
					in.m_a = 0;
	}
	VecCol&			VecCol :: operator=(VecCol &&in) {
					swap(m_a, in.m_a);
					return *this;
	}
	VecCol&			VecCol :: operator+=(const VecCol &in) {
					assert(this->size() == in.size());
					assert(this != &in);
					gsl_vector_add(m_a, in.get());
					return *this;
	}
	VecCol&			VecCol :: operator-=(const VecCol &in) {
					assert(this->size() == in.size());
					assert(this != &in);
					gsl_vector_sub(m_a, in.get());
					return *this;
	}
	VecCol&			VecCol :: operator*=(const long double scale) {
					gsl_vector_scale(m_a, (double)scale);
					return *this;
	}
	VecCol&			VecCol :: operator/=(const long double scale) {
					gsl_vector_scale(m_a, 1.0/double(scale));
					return *this;
	}
	void SquareMatrix :: set_identity() {
		gsl_matrix_set_identity(a);
	}
	SquareMatrix :: SquareMatrix(const size_t D) {
		a = gsl_matrix_alloc(D, D);
		gsl_matrix_set_zero(a);
	}
	SquareMatrix :: SquareMatrix(const size_t D, IdentityMatrixRequested) {
		a = gsl_matrix_alloc(D, D);
		gsl_matrix_set_zero(a);
		this->set_identity();
	}

	Matrix :: Matrix(const size_t rows, const size_t columns) {
		a = gsl_matrix_alloc(rows, columns);
		gsl_matrix_set_zero(a);
	}
	Matrix :: Matrix(const Matrix &in) {
					a = gsl_matrix_alloc(in.size1(), in.size2());
					gsl_matrix_memcpy(a, in.get());
	}
	Matrix&			Matrix :: operator=(const Matrix &in) {
					assert(this->size1() == in.size1());
					assert(this->size2() == in.size2());
					if(this != &in) {
						gsl_matrix_memcpy(a, in.get());
					}
					return *this;
	}
	Matrix :: Matrix(Matrix &&in) {
		a = in.a;
		in.a = 0;
	}
	Matrix&			Matrix :: operator=(Matrix &&in) {
		swap(a,in.a);
		return *this;
	}
	Matrix :: ~Matrix() {
		if(a) {
			gsl_matrix_free(a);
		} // otherwise, this has been moved-from
	}

#ifdef USE_GPERF_CPUPROFILE
	static size_t destruct_Square_counter = 0;
	static size_t destruct_Square_counter_free = 0;
	void dump_final_destruct_Square_counters() {
		PP(destruct_Square_counter, destruct_Square_counter_free);
	}
#endif
	/* Destructors and constructors */
	SquareMatrix :: ~SquareMatrix() {
		#ifdef USE_GPERF_CPUPROFILE
		++destruct_Square_counter;
		#endif
		if(a) {
			#ifdef USE_GPERF_CPUPROFILE
			++destruct_Square_counter_free;
			#endif
			gsl_matrix_free(a);
		} // otherwise, this has been moved-from
	}
	SquareMatrix :: SquareMatrix(const SquareMatrix &in) {
					a = gsl_matrix_alloc(in.size(), in.size());
					gsl_matrix_memcpy(a, in.get());
	}
	SquareMatrix&			SquareMatrix :: operator=(const SquareMatrix &in) {
				if(this->size() == in.size()) {
					if(this != &in) {
						gsl_matrix_memcpy(a, in.get());
					}
				} else {
					SquareMatrix copy(in);
					swap(this->a,copy.a);
				}
				return *this;
	}
	SquareMatrix :: SquareMatrix(SquareMatrix &&in) {
		a = in.a;
		in.a = 0;
	}
	SquareMatrix&			SquareMatrix :: operator=(SquareMatrix &&in) {
		swap(a,in.a);
		return *this;
	}




	SquareMatrix&		SquareMatrix :: operator+=(const SquareMatrix &in) {
					assert(this->size() == in.size());
					assert(this != &in);
					gsl_matrix_add(a, in.get());
					return *this;
	}
	SquareMatrix&		SquareMatrix :: operator-=(const SquareMatrix &in) {
					assert(this->size() == in.size());
					assert(this != &in);
					gsl_matrix_sub(a, in.get());
					return *this;
	}
	SquareMatrix&		SquareMatrix :: operator*=(const long double in) {
					gsl_matrix_scale(a, double(in));
					return *this;
	}
	SquareMatrix&		SquareMatrix :: operator/=(const long double in) {
					gsl_matrix_scale(a, 1.0 / double(in));
					return *this;
	}
	VecCol operator-(const VecCol &lhs, const VecCol &rhs) {
		assert(lhs.size() == rhs.size());
		VecCol diff(lhs.size());
		diff = lhs; // copy the data in
		gsl_vector_sub(diff.get(), rhs.get());
		return diff;
	}
	VecCol operator+(const VecCol &lhs, VecCol rhs) {
		// these two "if"s aren't beautiful. Just a hack for when
		// trying to add to a  single-element zero-vector
		if(rhs.len() == 1 && rhs(0) == 0) {
			return lhs;
		}
		if(lhs.len() == 1 && lhs(0) == 0) {
			return rhs;
		}
		assert(lhs.size() == rhs.size());
		gsl_vector_add(rhs.get(), lhs.get());
		return rhs;
	}
	VecCol operator*(const long double lhs, VecCol rhs) {
		gsl_vector_scale(rhs.get(), (double)lhs);
		return rhs;
	}
	VecCol operator/(VecCol rhs, const long double lhs) {
		gsl_vector_scale(rhs.get(), 1.0/(double)lhs);
		return rhs;
	}
VecCol			multiply_rowvec_by_matrix_giving_rowvec(const VecCol &rowvec, const SquareMatrix &matrix) {

	const size_t D = rowvec.size();
	assert(D == matrix.size());

	// somewhere to store the answer
	VecCol answer(D);

	// We need a matrix view of the lhs, and of the answer, because GSL refuses to do vector-matrix
	gsl_matrix_const_view mv_lhs    = gsl_matrix_const_view_vector(rowvec.get(), 1, D); // make it look like a row-vector
	gsl_matrix_view mv_answer = gsl_matrix_view_vector(answer.get(), 1, D); // make it look like a row-vector

	const int res_0 = gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, &(mv_lhs.matrix), matrix.get(), 0, &(mv_answer.matrix));
	assert(res_0 == 0);

	return answer;
}
VecCol			multiply_matrix_by_colvec_giving_colvec(const SquareMatrix &matrix, const VecCol &colvec) {
	const size_t D = matrix.size();
	assert(D == colvec.size());

	VecCol answer(D);
	const int res_0 = gsl_blas_dgemv (CblasNoTrans, 1.0, matrix.get(), colvec.get(), 0, answer.get());
	assert(res_0 == 0);
	return answer;
}
VecCol			multiply_matrix_by_colvec_giving_colvec(const       Matrix &matrix, const VecCol &colvec) {
	const size_t rows = matrix.get()->size1;
	const size_t cols = matrix.get()->size2;
	assert(cols == colvec.size());

	VecCol answer(rows);
	const int res_0 = gsl_blas_dgemv (CblasNoTrans, 1.0, matrix.get(), colvec.get(), 0, answer.get());
	assert(res_0 == 0);
	return answer;
}
VecCol			multiply_rowvec_by_colvec_giving_scalar(const VecCol &lhs, const VecCol &rhs) {
	const size_t D = lhs.size();
	assert(D == rhs.size());

	gsl_matrix_const_view mv_lhs    = gsl_matrix_const_view_vector(lhs.get(), 1, D); // make it look like a row-vector

	VecCol answer(1);
	const int res_0 = gsl_blas_dgemv (CblasNoTrans, 1.0, &(mv_lhs.matrix), rhs.get(), 0, answer.get());
	assert(res_0 == 0);
	return answer;
}
SquareMatrix		multiply_colvec_by_rowvec_giving_matrix(const VecCol &lhs, const VecCol &rhs) {
	const size_t D = lhs.size();
	assert(D == rhs.size());

	gsl_matrix_const_view mv_lhs    = gsl_matrix_const_view_vector(lhs.get(), D, 1); // make it look like a col-vector
	gsl_matrix_const_view mv_rhs    = gsl_matrix_const_view_vector(rhs.get(), 1, D); // make it look like a row-vector

	SquareMatrix answer(D);
	const int res_0 = gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, &(mv_lhs.matrix), &(mv_rhs.matrix), 0.0, answer.get());
	assert(res_0 == 0);
	return answer;
}
SquareMatrix		operator* (const SquareMatrix &lhs, const SquareMatrix &rhs) {
	const size_t D = lhs.size();
	assert(D == rhs.size());

	SquareMatrix answer(D);
	const int res_0 = gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, lhs.get(), rhs.get(), 0, answer.get());
	assert(res_0 == 0);
	return answer;
}
SquareMatrix		operator+ (const SquareMatrix &lhs, const SquareMatrix &rhs) {
	if(lhs.len() == 1 && lhs(0,0) == 0)
		return rhs;
	if(rhs.len() == 1 && rhs(0,0) == 0)
		return lhs;

	assert(lhs.size() == rhs.size());

	SquareMatrix answer(lhs);
	gsl_matrix_add(answer.get(), rhs.get());
	return answer;
}
SquareMatrix		operator- (const SquareMatrix &lhs, const SquareMatrix &rhs) {
	assert(lhs.size() == rhs.size());

	SquareMatrix answer(lhs);
	gsl_matrix_sub(answer.get(), rhs.get());
	return answer;
}
SquareMatrix		operator* (const long double lhs, SquareMatrix rhs) {
	/* Note, this works as desired. It copies the input,
	 * (taking it by value), then multiplies in place
	 * If you doubt me, see:  http://stackoverflow.com/a/7592741/146041
	 * :-)
	 */
	gsl_matrix_scale(rhs.get(), (double)lhs);
	return rhs;
}
SquareMatrix		operator/ (SquareMatrix lhs, long double rhs) {
	gsl_matrix_scale(lhs.get(), 1.0 / (double)rhs);
	return lhs;
}
static inline bool is_symmetric(const SquareMatrix &m) {
	const int D = (int)m.size();
	for(int r=0;r<D;++r) {
		for(int c=r+1;c<D;++c) {
			if ( gsl_matrix_get(m.get(), r,c)
			      !=gsl_matrix_get(m.get(), c,r)
			      )
				return false;
		}
	}
	return true;
}
#if 0
int64_t l2_determinant_of_PosDefMatrix_CountOfErrors = 0;
static long double			l2_determinant_of_PosDefMatrix(SquareMatrix LU) {
	//assert(is_symmetric(LU));
	const size_t D = LU.size();

	long double l2_deter_cholesky=0.0L;
	{
		gsl_error_handler_t * previous_handler = gsl_set_error_handler_off ();
		const int gsl_err = gsl_linalg_cholesky_decomp(LU.get());
		unless(gsl_err == 0) {
			// Not a positive definite matrix. Just set the determinant to as close to zero as possibble
			assert(gsl_err == GSL_EDOM);
			++l2_determinant_of_PosDefMatrix_CountOfErrors;
			return -1e+2000L;
		}
		assert(gsl_err == 0);
		gsl_set_error_handler(previous_handler);
	}

	for(size_t d = 0; d<D; ++d) {
		const double diagonal_entry = gsl_matrix_get(LU.a, d, d);
		l2_deter_cholesky += log2l(diagonal_entry);
	}
	return 2.0L*l2_deter_cholesky;
}
#endif
SquareMatrix			invert_a_matrix_impl(SquareMatrix another_copy_for_cholesky, const char * file, const size_t line) {
	assert(is_symmetric(another_copy_for_cholesky));
	(void)file; (void)line;
#if 0
	{
		SquareMatrix copy_with_which_to_do_LU(another_copy_for_cholesky);
		const size_t D = copy_with_which_to_do_LU.size();
		gsl_permutation * p = gsl_permutation_calloc(D);
		int signum;
		const int res_0 = gsl_linalg_LU_decomp(copy_with_which_to_do_LU.get(), p, &signum);
		assert(res_0 == 0);
		const double determinant = gsl_linalg_LU_det(copy_with_which_to_do_LU.get(), signum);
		if(determinant<=0) {
			PP(determinant, file, line);
			cerr << "Error: The input data is singular. Please change it. Exiting." << endl;
			exit(1);
		}
		assert(determinant > 0);
		SquareMatrix inverse(D);
		gsl_linalg_LU_invert(copy_with_which_to_do_LU.get(), p, inverse.get());
		gsl_permutation_free(p);
	}
#endif

	{ // calculate via Cholesky
		gsl_linalg_cholesky_decomp(another_copy_for_cholesky.get());
		gsl_linalg_cholesky_invert(another_copy_for_cholesky.get());
	}

	// Finally, return one of them!
	//inverse -= another_copy_for_cholesky;
	//assert(matrix_is_close_to_zero(inverse));
	return another_copy_for_cholesky;
}
SquareMatrix		cholesky_upper(SquareMatrix thematrix) {
		gsl_linalg_cholesky_decomp(thematrix.get());
		const size_t D = thematrix.size();
		for(size_t r = 0; r<D; ++r) {
		for(size_t d = 0; d<r; ++d) {
			gsl_matrix_set(thematrix.get(),r,d,0); // upper
			//gsl_matrix_set(thematrix.get(),d,r,0); // lower
		}
		}
		return thematrix;
}

#if 0
long double		l2_P_under_Normal		(const VecCol &value, const VecCol &mean, SquareMatrix precision) {
	const size_t D = value.size();
	assert(D == mean.size());
	assert(D == precision.size());

	VecCol offset = value - mean;
	//cout << "offset = (" << *value << "-" << *mean << ") = " << offset << endl;

	VecCol left_times_middle = multiply_rowvec_by_matrix_giving_rowvec(offset, precision);
	//cout << "left_time_middle = " << left_times_middle << endl;

	VecCol product_of_three = multiply_rowvec_by_colvec_giving_scalar(left_times_middle, offset);
	//cout << "product_of_three = " << product_of_three << endl;
	assert(product_of_three.size()==1);

	const long double l2_deter_of_COVARIANCE = - l2_determinant_of_PosDefMatrix(move(precision));
	assert(precision.get()==0);

	const long double l2_density =
		-(D/2.0) * log2(2*M_PI)
		- 0.5 * l2_deter_of_COVARIANCE
		+ (-0.5) * gsl_vector_get(product_of_three.get(),0) * M_LOG2E
		;
	//PP(l2_density);
	assert(isfinite(l2_density));
	return l2_density;
}
#endif


bool veccol_is_close_to_zero(const VecCol &x) {
	const size_t D = x.size();
	for(size_t d=0; d<D; ++d) {
		//PP(d, gsl_vector_get(x.a, d));
		const double x_d = gsl_vector_get(x.get(), d);
		assert(isfinite(x_d));
		if(abs(x_d) > 1.0) {
			return false;
		}
	}
	return true;
}
bool matrix_is_close_to_zero(const SquareMatrix &x) {
	const size_t D = x.size();
	for(size_t d=0; d<D; ++d) {
	for(size_t d2=0; d2<D; ++d2) {
		//PP(d, gsl_vector_get(x.a, d));
		const double x_d = gsl_matrix_get(x.a, d, d2);
		assert(isfinite(x_d));
		if(abs(x_d) > 1.0) {
			return false;
		}
	}
	}
	return true;
}
double norm_2(const SquareMatrix &x) {
	const size_t D = x.size();
	double total_x2 = 0.0;
	for(size_t d=0; d<D; ++d) {
	for(size_t d2=0; d2<D; ++d2) {
		const double x_d = gsl_matrix_get(x.a, d, d2);
		assert(isfinite(x_d));
		total_x2 += x_d*x_d;
	}
	}
	return sqrt(total_x2) / double(D*D);
}
double norm_2(const VecCol &x) {
	const size_t D = x.size();
	double total_x2 = 0.0;
	for(size_t d=0; d<D; ++d) {
		const double x_d = gsl_vector_get(x.get(), d);
		assert(isfinite(x_d));
		total_x2 += x_d*x_d;
	}
	return sqrt(total_x2) / double(D);
}

} // namespace mvn

	std :: ostream& std :: operator << (std :: ostream & o, const mvn :: VecCol &v) {
		std :: cout << "[";
		for(size_t d=0; d < v.size(); ++d) {
			if(d!=0)
				std :: cout << ',';
			std :: cout << gsl_vector_get(v.get(), d);
		}
		std :: cout << "]";
		return o;
	}
	std :: ostream& std :: operator << (std :: ostream & o, const mvn :: SquareMatrix &v) {
		for(size_t r=0; r < v.size(); ++r) {
			if(r==0)
				std :: cout << "[[";
			else
				std :: cout << " [";
			for(size_t c=0; c < v.size(); ++c) {
				std :: cout << '\t' << gsl_matrix_get(v.a, r, c);
			}
			if(r+1 == v.size())
				cout << "]]";
			else
				cout << "]" << endl;
		}
		return o;
	}
long double
mvn:: calculate_lndeterminant(mvn:: SquareMatrix copy_with_which_to_do_LU) {
    const size_t D = copy_with_which_to_do_LU.size();
    gsl_permutation * p = gsl_permutation_calloc(D);
    int signum;
    const int res_0 = gsl_linalg_LU_decomp(copy_with_which_to_do_LU.get(), p, &signum);
    assert(res_0 == 0);

    const double determinant = gsl_linalg_LU_lndet(copy_with_which_to_do_LU.get());
    gsl_permutation_free(p);
    assert(std::isfinite(determinant));
    return determinant;
}
