/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _petsc_wrapper_H
#define _petsc_wrapper_H

#include <inttypes.h>
#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"

#include "file_utils.h"
#include "general_utils.h"

class cOwnership;//forward declaration only
class cPetscDistVector;//forward declaration only
class cPetscDistMatrix;//forward declaration only
class cPetscDistShellMatrix;//forward declaration only

#if PETSC_VERSION_GE(3,5,0)
#define __SDIR__ "src\\"
#endif

#define CHKERR(ierr) cPetscObject::chkerrabort(ierr, __LINE__, __FUNCTION__,__FILE__,__SDIR__)
#define dbgprint cPetscObject::debugprint(__LINE__, __FUNCTION__,__FILE__,__SDIR__)

class cSparseMatrixEntries{

public:
	std::vector<PetscInt> rowind;
	std::vector<PetscInt> colind;
	std::vector<double> value;

	cSparseMatrixEntries(){ resize(0); }
	cSparseMatrixEntries(const PetscInt nnz){ resize(nnz); }

	void resize(const PetscInt nnz){
		rowind.resize(nnz);
		colind.resize(nnz);
		value.resize(nnz);
	};	
};

class cOwnership{
	
public:

	PetscInt start;//index of first 
	PetscInt end;//index of last + 1

	cOwnership(){
		start = -1;
		end   = -1;
	};

	cOwnership(const PetscInt _start, const PetscInt _end){
		start = _start;
		end   = _end;
	}

	cOwnership(const int size, const int rank, const PetscInt nglobal){
		set_petsc_default(size, rank, nglobal);
	}

	void set_petsc_default(const int size, const int rank, const PetscInt nglobal)
	{		
		end = 0;
		for (PetscInt i = 0; i <= rank; i++){
			PetscInt n = nglobal / (PetscInt)size + ((nglobal % (PetscInt)size) > i);
			end += n;
			start = end - n;
		}
	}

	PetscInt nlocal() const
	{
		return end - start;
	};  //number of local items
	
	PetscInt globalind(const PetscInt& ilocal) const 
	{
		return start + ilocal; 
	}

	PetscInt localind(const PetscInt& iglobal) const
	{
		return iglobal - start; 
	}

	bool owns(const PetscInt& iglobal) const
	{
		if (iglobal < start || iglobal >= end) return false;
		return true;
	}

};

class cPetscObject{
	
private:
	
public:	
		
	cPetscObject(){ };

	virtual const PetscObject* pobj() const = 0;

	void setname(const std::string& inname) {
		if (pobj()){
			PetscErrorCode ierr = PetscObjectSetName(*pobj(), inname.c_str()); CHKERR(ierr);
		}
	}

	const char* cname() const
	{
		if (pobj()){
			const char* name;
			PetscErrorCode ierr = PetscObjectGetName(*pobj(), &name);
			CHKERR(ierr);
			return name;
		}
		else return (char*)NULL;
	}

	const std::string name() const
	{
		if (cname()){
			return std::string(cname());
		}
		return std::string("unnamed");
	}

	MPI_Comm mpicomm() const
	{
		MPI_Comm comm;
		PetscErrorCode ierr = PetscObjectGetComm(*pobj(), &comm);
		CHKERR(ierr);
		return comm;
	}

	int mpisize() const 
	{
		int size;
		PetscErrorCode ierr = MPI_Comm_size(mpicomm(), &size);
		CHKERR(ierr);
		return size;
	}

	int mpirank() const 
	{
		int rank;
		PetscErrorCode ierr = MPI_Comm_rank(mpicomm(), &rank); CHKERR(ierr);
		return rank;
	}

	PetscErrorCode mpibarrier() const 
	{
		PetscErrorCode ierr = MPI_Barrier(mpicomm()); CHKERR(ierr);
		return ierr;
	}

	std::string mpiprocname() const 
	{
		char name[MPI_MAX_PROCESSOR_NAME + 1];
		int len;
		PetscErrorCode ierr = MPI_Get_processor_name(name, &len); CHKERR(ierr);
		return std::string(name);
	}

	static void printsizeofs()
	{
		printf("sizeof(MPI_Int) = %zu\n", sizeof(MPI_INT));
		printf("sizeof(bool) = %zu\n", sizeof(bool));
		printf("sizeof(int) = %zu\n", sizeof(int));
		printf("sizeof(int32_t) = %zu\n", sizeof(int32_t));
		printf("sizeof(int64_t) = %zu\n", sizeof(int64_t));
		printf("sizeof(size_t) = %zu\n", sizeof(size_t));
		printf("sizeof(PetscInt) = %zu\n", sizeof(PetscInt));
		printf("sizeof(PetscScalar) = %zu\n", sizeof(PetscScalar));
#if PETSC_VERSION_LT(3,8,0)
		printf("sizeof(Petsc64bitInt) = %zu\n", sizeof(Petsc64bitInt));
#else
		printf("sizeof(PetscInt64) = %zu\n", sizeof(PetscInt64));
#endif
		printf("sizeof(PetscBLASInt) = %zu\n", sizeof(PetscBLASInt));
		printf("sizeof(PetscMPIInt) = %zu\n", sizeof(PetscMPIInt));
	};

	static void chkerrabort(PetscErrorCode ierr, const int linenumber, const char* functionname, const char* srcfile, const char* srcdirectory){
		do {
			if (PetscUnlikely(ierr)){
#if PETSC_VERSION_LT(3,5,0)
				PetscError(PETSC_COMM_SELF,linenumber,functionname,srcfile,srcdirectory,ierr,PETSC_ERROR_REPEAT," ");
#else
				PetscError(PETSC_COMM_SELF,linenumber,functionname,srcfile,ierr,PETSC_ERROR_REPEAT," ");
#endif
				MPI_Abort(PETSC_COMM_SELF, ierr);
			}
		} while (0);		
	}

	static void debugprint(const int linenumber, const char* functionname, const char* srcfile, const char* srcdirectory){
		PetscPrintf(PETSC_COMM_WORLD, "**Debug: at line %d : function %s : file %s : directory %s\n", linenumber, functionname, srcfile, srcdirectory);
	}

};



class cPetscDistVector : public cPetscObject {

	friend class cPetscDistMatrix;
	friend class cPetscDistShellMatrix;

private:
	bool isowner = true;
	Vec* pVec = PETSC_NULL;	

	const PetscObject* pobj() const{
		if (pVec) return (PetscObject*)pVec;
		else return (PetscObject*)PETSC_NULL;
	}

	void deallocate()
	{
		if (allocated()){
			if (isowner){
				//printf("Destroying Vector %s\n", name().c_str());
				PetscErrorCode ierr = VecDestroy(pVec);	CHKERR(ierr);
				pVec = PETSC_NULL;				
			}
		}
	}

	void resize(const MPI_Comm comm, const PetscInt _localsize, const PetscInt _globalsize)
	{
		std::string myname = name();
		if (allocated()){
			if (localsize() == _localsize && globalsize() == _globalsize){
				return;
			}
			else deallocate();			
		}
		create(myname, comm, _localsize, _globalsize);
	}

public:
	
	cPetscDistVector() {};

	cPetscDistVector(const std::string inname, const MPI_Comm comm, const PetscInt localsize, const PetscInt globalsize, const double value = 0.0)
	{				
		create(inname, comm, localsize, globalsize, value);
	}
	
	cPetscDistVector(const MPI_Comm comm, const PetscInt localsize, const PetscInt globalsize)
	{				
		create(comm, localsize, globalsize);		
	}

	cPetscDistVector(const cPetscDistVector& rhs)
	{				
		PetscErrorCode ierr;		
		pVec = new Vec;
		ierr = VecDuplicate(rhs.vec(), ncvecptr()); CHKERR(ierr);
		ierr = VecCopy(rhs.vec(), vec()); CHKERR(ierr);
		//printf("Copied %s to %s\n", rhs.cname(), cname());
	}

	cPetscDistVector(Vec& vec)
	{		
		isowner = false;
		pVec = &vec;
	}

	cPetscDistVector(cPetscDistVector&& rhs)
	{				
		//std::swap(isowner, rhs.isowner);
		std::swap(pVec, rhs.pVec);
	}

	~cPetscDistVector()
	{
		//printf("Destructor for Vector %s\n", name().c_str());
		deallocate();
	}

	const Vec& vec() const { return *pVec; }
	Vec& ncvec() const { return *pVec; }
	const Vec* vecptr() const { return pVec; }
	Vec* ncvecptr() const { return pVec; }

	PetscInt gi(const PetscInt& localind) const
	{
		return ownership().globalind(localind);
	};

	PetscInt li(const PetscInt& globalind) const
	{
		return ownership().localind(globalind);
	};

	bool allocated() const {
		if (pVec && *pVec)return true;
		else return false;
	}

	void create(const MPI_Comm comm, const PetscInt _localsize, const PetscInt _globalsize)
	{				
		pVec = new Vec;
		PetscErrorCode ierr = VecCreateMPI(comm, _localsize, _globalsize, ncvecptr()); CHKERR(ierr);
	}

	void create(const std::string& inname, const MPI_Comm comm, const PetscInt _localsize, const PetscInt _globalsize)
	{
		create(comm, _localsize, _globalsize);
		setname(inname);		
	}

	void create(const std::string& inname, const MPI_Comm comm, const PetscInt _localsize, const PetscInt _globalsize, const double value)
	{
		create(comm, _localsize, _globalsize);				
		setname(inname);		
		set(value);		
	}

	cOwnership ownership() const
	{
		PetscInt start;
		PetscInt end;
		PetscErrorCode ierr = VecGetOwnershipRange(vec(), &start, &end); CHKERR(ierr);
		return cOwnership(start, end);
	}

	PetscInt globalsize() const
	{
		if (vec()){
			PetscInt N;
			PetscErrorCode ierr = VecGetSize(vec(), &N); CHKERR(ierr);
			return N;
		}
		else return 0;
	}

	PetscInt localsize() const
	{
		if (pVec){
			PetscInt n;
			PetscErrorCode ierr = VecGetLocalSize(vec(), &n); CHKERR(ierr);
			return n;
		}
		else return 0;
	}

	const double* getlocalreadonlyarray() const
	{
		const double* localarray;
		if (allocated()==false){
			PetscPrintf(mpicomm(), "Error cPetscDistVector::getlocalarray() (%s) attempting to access NULL vec() pointer\n", name().c_str());
			throw(strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__));
		}
		PetscErrorCode ierr = VecGetArrayRead(vec(), &localarray); CHKERR(ierr);
		return localarray;
	}

	void restorelocalreadonlyarray(const double* localarray) const
	{
		PetscErrorCode ierr = VecRestoreArrayRead(vec(), &localarray); CHKERR(ierr);
	}

	double* getlocalarray()
	{
		double* localarray;
		if (vec() == 0){
			PetscPrintf(mpicomm(), "Error cPetscDistVector::getlocalarray() (%s) attempting to access NULL vec() pointer\n", cname());
			throw(strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__));
		}
		PetscErrorCode ierr = VecGetArray(vec(), &localarray); CHKERR(ierr);
		return localarray;
	}

	void restorelocalarray(double* localarray)
	{
		PetscErrorCode ierr = VecRestoreArray(vec(), &localarray); CHKERR(ierr);
	}

	void set(const double& value)
	{
		PetscErrorCode ierr = VecSet(vec(), value);
		CHKERR(ierr);
		return;
	}

	void set_local(const std::vector<double>& v)
	{
		double* a = getlocalarray();
		for (PetscInt i = 0; i<localsize(); ++i){
			a[i] = v[i];
		}
		restorelocalarray(a);
	}

	std::vector<double> getglobal() const
	{
		std::vector<double> v(globalsize());
		double* a;
		VecScatter ctx;
		Vec V_SEQ;
		VecScatterCreateToAll(vec(), &ctx, &V_SEQ);
		VecScatterBegin(ctx, vec(), V_SEQ, INSERT_VALUES, SCATTER_FORWARD);
		VecScatterEnd(ctx, vec(), V_SEQ, INSERT_VALUES, SCATTER_FORWARD);
		VecGetArray(V_SEQ, &a);
		for (size_t i = 0; i<(size_t)globalsize(); i++) v[i] = a[i];
		VecRestoreArray(V_SEQ, &a);
		VecScatterDestroy(&ctx);
		VecDestroy(&V_SEQ);
		return v;
	}

	double dot(const cPetscDistVector& b) const
	{
		double x;
		PetscErrorCode ierr = VecDot(vec(), b.vec(), &x); CHKERR(ierr);
		return x;
	}

	double norm2() const
	{
		double x;
		PetscErrorCode ierr = VecNorm(vec(), NORM_2, &x); CHKERR(ierr);
		return x;
	}

	double l2_norm() const
	{
		double x;
		PetscErrorCode ierr = VecDot(vec(), vec(), &x); CHKERR(ierr);
		return x;
	}

	friend double dot(const cPetscDistVector& v1, const cPetscDistVector& v2);
	friend cPetscDistVector normalized_residual(const cPetscDistVector& v0, const cPetscDistVector& v, const cPetscDistVector& ev);
	friend double data_misfit(const cPetscDistVector& v0, cPetscDistVector& v, const cPetscDistVector& ev);
	friend double data_misfit_normalized(const cPetscDistVector& v0, const cPetscDistVector& v, const cPetscDistVector& ev);


	double localmemory() const
	{
		return (double)localsize()*sizeof(PetscReal);	
	}

	double globalmemory() const
	{
		return (double)globalsize()*sizeof(PetscReal);
	}

	void printsize() const
	{
		if (!vec()){
			PetscPrintf(mpicomm(), "Vector %s (not valid)\n", cname());
			return;
		}
		PetscPrintf(mpicomm(), "Vector (%s) size = %ld\n", cname(), (int64_t)globalsize());
	}
	void printinfo()const
	{
		VecType type;
		PetscErrorCode ierr = VecGetType(vec(), &type); CHKERR(ierr);
		PetscPrintf(mpicomm(), "\n");
		PetscPrintf(mpicomm(), "==== Vector %s (%s) ================\n", cname(), type);
		PetscPrintf(mpicomm(), "Global rows          \t%7ld\n", (int64_t)globalsize());
		PetscPrintf(mpicomm(), "=======================================\n");
		PetscPrintf(mpicomm(), "\n");

	}
	void printdistributions() const
	{
		for (int i = 0; i < mpisize(); i++){
			mpibarrier();
			if (i == mpirank()){
				printf("%s rank %d owns %d to %d  - total of %d rows\n",
					mpiprocname().c_str(),
					mpirank(),
					ownership().start,
					ownership().end - 1,
					localsize());
				fflush(stdout);
			}
			mpibarrier();
		}
	}
	void printvalues(const std::string& format)
	{
		double* a = getlocalarray();
		for (int p = 0; p < mpisize(); p++){
			mpibarrier();
			if (p == mpirank()){
				cOwnership r = ownership();
				for (PetscInt i = 0; i<r.nlocal(); i++){
					printf(format.c_str(), a[i]);
					fflush(stdout);
				}
				if (p == (mpisize() - 1)){
					printf("=============\n");
					fflush(stdout);
				}
			}
			mpibarrier();
		}
		restorelocalarray(a);
		mpibarrier();
	}
	void writetextfile(const std::string& filename)
	{
		double* a = getlocalarray();
		for (int p = 0; p < mpisize(); p++){
			mpibarrier();
			if (p == mpirank()){
				FILE* fp;
				if (p == 0){
					fp = fileopen(filename, "w");
				}
				else{
					fp = fileopen(filename, "a");
				}

				cOwnership r = ownership();
				for (PetscInt i = 0; i<r.nlocal(); i++){
					fprintf(fp, "%20.16le\n", a[i]);
				}
				fclose(fp);
			}
			mpibarrier();
		}
		restorelocalarray(a);
		mpibarrier();
	}

	cPetscDistVector& operator=(const cPetscDistVector& rhs)
	{				
		PetscErrorCode ierr;				
		if (localsize() != rhs.localsize() || globalsize() != rhs.globalsize() ){			
			deallocate();
			pVec = new Vec;
			ierr = VecDuplicate(rhs.vec(), ncvecptr()); CHKERR(ierr);
		}
		ierr = VecCopy(rhs.vec(), ncvec()); CHKERR(ierr);					
		//printf("Assigned %s = %s\n", cname(), rhs.cname());
		return *this;
	}

	cPetscDistVector& operator<<(const cPetscDistVector& rhs)
	{
		PetscErrorCode ierr;		
		ierr = VecCopy(rhs.vec(), ncvec()); CHKERR(ierr);		
		return *this;
	}

	cPetscDistVector& operator+=(const cPetscDistVector& rhs)
	{
		PetscErrorCode ierr = VecAXPY(vec(), 1.0, rhs.vec()); CHKERR(ierr);
		return *this;
	}
	
	cPetscDistVector& operator-=(const cPetscDistVector& rhs)
	{
		PetscErrorCode ierr = VecAXPY(vec(), -1.0, rhs.vec()); CHKERR(ierr);
		return *this;
	}

	cPetscDistVector& operator*=(const cPetscDistVector& rhs)
	{
		PetscErrorCode ierr = VecPointwiseMult(vec(), vec(), rhs.vec()); CHKERR(ierr);
		return *this;
	}

	cPetscDistVector& operator*=(const double& s)
	{
		PetscErrorCode ierr = VecScale(vec(), s); CHKERR(ierr);
		return *this;
	}

	cPetscDistVector& operator/=(const cPetscDistVector& rhs)
	{
		PetscErrorCode ierr = VecPointwiseDivide(vec(), vec(), rhs.vec()); CHKERR(ierr);
		return *this;
	}

	cPetscDistVector& operator=(const double& value)
	{
		set(value);		
		return *this;
	}

	cPetscDistVector& operator+=(const double& value)
	{
		double* a = getlocalarray();
		cOwnership r = ownership();
		for (PetscInt i = 0; i<r.nlocal(); ++i){
			a[i] += value;
		}
		restorelocalarray(a);
		return *this;
	}

	cPetscDistVector operator+(const cPetscDistVector& b) const
	{
		cPetscDistVector c(*this);
		c += b;
		return c;
	}

	cPetscDistVector operator-(const cPetscDistVector& b) const
	{
		cPetscDistVector c(*this);
		c -= b;
		return c;
	}

	cPetscDistVector operator*(const cPetscDistVector& b) const
	{
		cPetscDistVector c(*this);
		c *= b;
		return c;
	}

	cPetscDistVector operator/(const cPetscDistVector& b) const
	{
		cPetscDistVector c(*this);
		c /= b;
		return c;
	}

	cPetscDistVector operator*(const double& s) const
	{
		cPetscDistVector c(*this);
		PetscErrorCode ierr = VecScale(c.vec(), s);	CHKERR(ierr);
		return c;
	}

	//Friends	
	friend double vtAv(const cPetscDistMatrix& A, const cPetscDistVector& v);
	friend cPetscDistVector operator*(const double& s, const cPetscDistVector& v)
	{
		cPetscDistVector b(v);
		PetscErrorCode ierr = VecScale(b.vec(), s); CHKERR(ierr);
		return b;
	}
	
	void pow(const double& p)
	{		
		double* v = getlocalarray();		
		for (size_t i = 0; i < (size_t)localsize(); i++){
			v[i] = std::pow(v[i],p);
		}
		restorelocalarray(v);
		return;
	}
		
};

class cPetscDistVectorLocalView{

private:
	cPetscDistVector* pv;
	double* data;

public:

	cPetscDistVectorLocalView(cPetscDistVector& v){
		pv = &v;
		data = pv->getlocalarray();
	}

	~cPetscDistVectorLocalView(){
		pv->restorelocalarray(data);
	}

	double& operator[](const size_t& i){
		return data[i];
	}

};

class cPetscDistMatrix : public cPetscObject {

	friend class cPetscDistShellMatrix;

private:
		
	Mat* pMat = PETSC_NULL;
	bool isowner = true;
	std::vector<PetscInt> mRowInd;
	std::vector<PetscInt> mColInd;

	const PetscObject* pobj() const{
		if (pMat) return (PetscObject*)pMat;
		else return (PetscObject*)PETSC_NULL;
	}

	bool allocated() const {
		if (pMat && *pMat)return true;
		else return false;
	}

	void deallocate()
	{
		if (allocated()){
			if (isowner){
				//printf("Destroying Matrix %s\n", name().c_str());
				PetscErrorCode ierr = MatDestroy(pMat); CHKERR(ierr);
				pMat = PETSC_NULL;
			}
		}
	}

	bool init_sparse(const MPI_Comm comm, const PetscInt& _nlocalrows, const PetscInt& _nlocalcols, const PetscInt& _nglobalrows, const PetscInt& _nglobalcols){
		deallocate();
		pMat = new Mat;
		PetscErrorCode ierr;
		ierr = MatCreate(comm, ncmatptr()); CHKERR(ierr);
		ierr = MatSetSizes(mat(), _nlocalrows, _nlocalcols, _nglobalrows, _nglobalcols); CHKERR(ierr);
		ierr = MatSetType(mat(), MATMPIAIJ); CHKERR(ierr);
		ierr = MatSetUp(mat()); CHKERR(ierr);
		return true;
	}

public:
	
	//Constructors & Destructors
	cPetscDistMatrix() {};
	
	cPetscDistMatrix(const cPetscDistMatrix& rhs)
	{					
		PetscErrorCode ierr;
		pMat = new Mat;
		ierr = MatDuplicate(rhs.mat(), MAT_COPY_VALUES, pMat); CHKERR(ierr);
		printf("Copied Mat %s to Mat %s\n", rhs.cname(), cname());
	}

	cPetscDistMatrix(cPetscDistMatrix&& rhs)
	{
		std::swap(pMat, rhs.pMat);
	}

	~cPetscDistMatrix()
	{
		deallocate();
	}

	//Member Functions				
	const Mat& mat() const { return *pMat; }
	Mat& ncmat() const { return *pMat; }
	const Mat* matptr() const { return pMat; }
	Mat* ncmatptr() const { return pMat; }
	
	PetscInt nglobalrows() const
	{
		PetscErrorCode ierr;
		PetscInt M, N;
		ierr = MatGetSize(mat(), &M, &N); CHKERR(ierr);
		return M;
	};

	PetscInt nglobalcols() const
	{
		PetscErrorCode ierr;
		PetscInt M, N;
ierr = MatGetSize(mat(), &M, &N); CHKERR(ierr);
		return N;
	};

	PetscInt nlocalrows() const
	{
		PetscErrorCode ierr;
		PetscInt m, n;
ierr = MatGetLocalSize(mat(), &m, &n); CHKERR(ierr);
		return m;
	};

	PetscInt nlocalcols() const
	{
		//All cols are actually local but this is about the 
		//local cols of a matrix vector multiply
		//see http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatMPIAIJSetPreallocation.html		
		PetscErrorCode ierr;
		PetscInt m, n;
ierr = MatGetLocalSize(mat(), &m, &n); CHKERR(ierr);
		return n;
	};

	const cOwnership rowownership() const
	{
		PetscInt start;
		PetscInt end;
		PetscErrorCode ierr = MatGetOwnershipRange(mat(), &start, &end); CHKERR(ierr);
		return cOwnership(start, end);
	}

	const cOwnership columnownership() const
	{
		PetscInt start;
		PetscInt end;
		PetscErrorCode ierr = MatGetOwnershipRangeColumn(mat(), &start, &end); CHKERR(ierr);
		return cOwnership(start, end);
	}

	PetscInt gri(const PetscInt& localrow) const
	{
		return rowownership().globalind(localrow);
	};

	PetscInt lri(const PetscInt& globalrow) const
	{
		return rowownership().localind(globalrow);
	};

	void create_sparse(const std::string& name, const MPI_Comm comm, const PetscInt& _nlocalrows, const PetscInt& _nlocalcols, const PetscInt& _nglobalrows, const PetscInt& _nglobalcols){
		init_sparse(comm, _nlocalrows, _nlocalcols, _nglobalrows, _nglobalcols);
		setname(name);		
	}

	bool create_diagonal_structure(const std::string& name, const MPI_Comm comm, const PetscInt& localsize, const PetscInt& globalsize)
	{
		create_sparse(name, comm, localsize, localsize, globalsize, globalsize);
		PetscErrorCode ierr;
		ierr = MatMPIAIJSetPreallocation(mat(), 1, PETSC_NULL, 0, PETSC_NULL); CHKERR(ierr);		
		const cOwnership r = rowownership();		
		for (PetscInt gi = r.start; gi<r.end; gi++){
			ierr = MatSetValue(mat(), gi, gi, 0.0, INSERT_VALUES); CHKERR(ierr);
		}
		assemble();		
		return true;
	}

	bool create_identity(const std::string& name, const MPI_Comm comm, const PetscInt localsize, const PetscInt globalsize)
	{
		create_diagonal_structure(name, comm, localsize, globalsize);
		std::vector<double> localdiagonal(localsize, 1.0);
		set_diagonal_local(localdiagonal);
		return true;
	}

	bool create_diagonal(const std::string& name, const cPetscDistVector& diagonal)
	{
		//W_ii = d_i
		PetscErrorCode ierr;
		create_diagonal_structure(name, diagonal.mpicomm(), diagonal.localsize(), diagonal.globalsize());
		ierr = MatDiagonalSet(mat(), diagonal.vec(), INSERT_VALUES); CHKERR(ierr);
		return true;
	}

	bool create_diagonal_to_power(const std::string& name, const cPetscDistVector& v, const double& power)
	{
		//W_ii = (1/d_i)^2		
		cPetscDistVector d = v;
		d.pow(power);
		create_diagonal(name, d);		
		return true;
	}

	bool create_diagonal_local(const std::string& name, const MPI_Comm comm, const std::vector<double>& localdiagonal)
	{				
		create_diagonal_structure(name,comm,(PetscInt)localdiagonal.size(),PETSC_DETERMINE);
		set_diagonal_local(localdiagonal);				
		return true;
	}

	bool create_dense(const std::string& name, const MPI_Comm comm, const PetscInt numrows, const PetscInt numcols, const std::string mattype)
	{
		PetscErrorCode ierr;

		//Create and set sizes
		ierr = MatCreate(comm, ncmatptr()); CHKERR(ierr);
		ierr = MatSetSizes(mat(), PETSC_DECIDE, PETSC_DECIDE, numrows, numcols); CHKERR(ierr);

		//Set the matrix type from the petsc options
		//PetscBool issetflg;
		//const MatType defmtype = MATMPIDENSE;
		//char mtype[256];
		//ierr = PetscOptionsList("-mat_type", "Matrix type", "MatSetType", MatList, defmtype, mtype, 256, &issetflg); CHKERR(ierr);
		//if (issetflg){
		//	ierr = MatSetType(mat(),mtype); CHKERR(ierr);
		//}
		//else{
		//	ierr = MatSetType(mat(),defmtype); CHKERR(ierr);
		//}

		//Set the matrix type
		if (mattype != "mpidense" && mattype != "elemental"){
			PetscPrintf(comm, "Error: cPetscDistMatrix::create_dense() the matrix type should be either be 'mpidense' or 'elemental' \n");
			throw(strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__));
		}
		ierr = MatSetType(mat(), mattype.c_str()); CHKERR(ierr);

		//Setup and preallocation
		ierr = MatSetUp(mat()); CHKERR(ierr);

		//Create arrays that hold the local processe's owned indices
		PetscPrintf(PETSC_COMM_WORLD, "Entering setindexset\n");
		setindexset();
		PetscPrintf(PETSC_COMM_WORLD, "Finished setindexset\n");
		return true;
	}

	bool create_fromglobalindices(const std::string& name, const MPI_Comm comm, const std::vector<PetscInt>& rowind, const std::vector<PetscInt>& colind, const std::vector<double>& values, const PetscInt _nlocalrows, const PetscInt _nlocalcols, PetscInt _nglobalrows, PetscInt _nglobalcols)
	{
		PetscErrorCode ierr;
		auto prmm = std::minmax_element(rowind.begin(), rowind.end());
		auto pcmm = std::minmax_element(colind.begin(), colind.end());
		PetscInt rmin = *prmm.first;
		PetscInt rmax = *prmm.second;
		PetscInt cmin = *pcmm.first;
		PetscInt cmax = *pcmm.second;

		if (_nglobalrows <= 0) _nglobalrows = rmax + 1;
		if (_nglobalcols <= 0) _nglobalcols = cmax + 1;

		PetscInt nnz = (PetscInt)rowind.size();
		if ((int)colind.size() != nnz || (int)values.size() != nnz){
			printf("cPetscDistMatrix::create_globalindices(...) rowind, colind, and vals must be the same size\n");
			throw(strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__));
		}
		if (rmin<0 || rmax >= _nglobalrows){
			printf("cPetscDistMatrix::create_globalindices(...) row index out of range\n");
			throw(strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__));
		}
		if (cmin<0 || cmax >= _nglobalcols){
			printf("cPetscDistMatrix::create_globalindices(...) column index out of range\n");
			throw(strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__));
		}

		create_sparse(name, comm, _nlocalrows, _nlocalcols, _nglobalrows, _nglobalcols);

		PetscInt nlr = nlocalrows();
		std::vector<PetscInt> d_nnz(nlr, 0);
		std::vector<PetscInt> o_nnz(nlr, 0);

		cOwnership rown = rowownership();
		cOwnership cown = columnownership();

		//Preallocate
		for (PetscInt i = 0; i<nnz; ++i){
			if (rown.owns(rowind[i])){
				PetscInt lr = lri(rowind[i]);
				if (cown.owns(colind[i])){
					d_nnz[lr]++;
				}
				else{
					o_nnz[lr]++;
				}
			}
		}

		//Cater for duplicates in one location
		PetscInt nlc = nlocalcols();
		for (PetscInt k = 0; k < nlr; k++){
			if (d_nnz[k] > nlc)d_nnz[k] = nlc;
			if (o_nnz[k] > _nglobalcols - nlc)o_nnz[k] = _nglobalcols - nlc;
		}
		ierr = MatMPIAIJSetPreallocation(mat(), 0, d_nnz.data(), 0, o_nnz.data()); CHKERR(ierr);

		//Now set the values
		for (PetscInt i = 0; i < nnz; ++i){
			if (rown.owns(rowind[i])){
				ierr = MatSetValue(mat(), rowind[i], colind[i], values[i], ADD_VALUES); CHKERR(ierr);
			}
		}
		ierr = MatAssemblyBegin(mat(), MAT_FINAL_ASSEMBLY);	CHKERR(ierr);
		ierr = MatAssemblyEnd(mat(), MAT_FINAL_ASSEMBLY); CHKERR(ierr);
		return true;
	}

	bool preallocate(const PetscInt d_nz, const PetscInt o_nz) const
	{
		PetscErrorCode ierr = MatMPIAIJSetPreallocation(mat(), d_nz, PETSC_NULL, o_nz, PETSC_NULL); CHKERR(ierr);
		return true;
	}

	bool preallocate(std::vector<PetscInt>& d_nnz, std::vector<PetscInt>& o_nnz) const
	{
		PetscErrorCode ierr = MatMPIAIJSetPreallocation(mat(), (PetscInt)PETSC_NULL, d_nnz.data(), (PetscInt)PETSC_NULL, o_nnz.data());
		CHKERR(ierr);
		return true;
	}

	bool inownerdiagonalblock(const PetscInt& ri, const PetscInt& ci) const
	{		
		if (rowownership().owns(ri) && columnownership().owns(ci)) return true;
		else return false;
	}

	PetscInt ownsrow(const PetscInt& globalrow) const
	{
		return rowownership().owns(globalrow);
	}

	bool isassembled() const {
		if (mat() == PETSC_NULL)return false;
		PetscBool status;
		PetscErrorCode ierr = MatAssembled(mat(), &status); CHKERR(ierr);
		if (status == PETSC_TRUE)return true;
		return false;
	}
	
	void set(const PetscInt globalrow, const PetscInt globalcol, double value)
	{
		PetscErrorCode  ierr;
		ierr = MatSetValue(mat(), globalrow, globalcol, value, INSERT_VALUES);
		CHKERR(ierr);
	};

	bool set_diagonal_local(const std::vector<double>& localdiagonal){
		
		if (nglobalrows() != nglobalcols()){
			printf("cPetscDistMatrix::set_diagonal_local() not a square matrix\n");
			throw(strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__));
		}
		else if((int)localdiagonal.size() != nlocalrows()){
			printf("cPetscDistMatrix::set_diagonal_local() incompatible size\n");
			throw(strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__));
		}
		
		PetscErrorCode ierr;
		const cOwnership r = rowownership();
		PetscInt k = 0;
		for (PetscInt gi = r.start; gi < r.end; gi++){
			ierr = MatSetValue(mat(), gi, gi, localdiagonal[k], INSERT_VALUES); CHKERR(ierr);
			k++;
		}
		assemble();		
		return true;
	}

	PetscErrorCode setrow(const PetscInt globalrow,  const std::vector<PetscInt>& globalcols, const std::vector<double>& values)
	{				
		PetscErrorCode  ierr;
		ierr = MatSetValues(mat(),1,&globalrow, (PetscInt)globalcols.size(),globalcols.data(),values.data(),INSERT_VALUES); CHKERR(ierr);
		return ierr;
	}

	PetscErrorCode setrowdense(const PetscInt globalrow,  const std::vector<PetscInt>& globalcols, const std::vector<double>& rowvalues)
	{		
		PetscErrorCode  ierr;		
		PetscInt nlr = nlocalrows();
		PetscInt nc  = nglobalcols();
		double* a=NULL;
		Mat lMat;		
		PetscInt localrow = lri(globalrow);
		ierr = MatDenseGetLocalMatrix(mat(), &lMat); CHKERR(ierr);
		ierr = MatDenseGetArray(lMat, &a); CHKERR(ierr);
		for(PetscInt j=0; j<nc; ++j){ 			
			a[localrow + j*nlr] = rowvalues[j];
		}
		ierr = MatDenseRestoreArray(lMat, &a); CHKERR(ierr);
		return ierr;
	};

	double get(const PetscInt globalrow, const PetscInt globalcol) const
	{
		double value;		
		PetscErrorCode ierr = MatGetValues(mat(), 1, &globalrow, 1, &globalcol, &value); CHKERR(ierr);
		return value;
	};

	PetscErrorCode assemble(){
		PetscErrorCode ierr;
		ierr = MatAssemblyBegin(mat(), MAT_FINAL_ASSEMBLY); CHKERR(ierr);
		ierr = MatAssemblyEnd(mat(), MAT_FINAL_ASSEMBLY); CHKERR(ierr);
		return ierr;
	};

	void add(const cPetscDistMatrix& A, const double s=1.0)
	{					
		PetscErrorCode ierr = MatAXPY(mat(), s, A.mat(), DIFFERENT_NONZERO_PATTERN); CHKERR(ierr);
	}

	void setindexset(){
		IS isrows, iscols;
		PetscInt lnrow, lncol;
		const PetscInt *rind, *cind;
		MatGetOwnershipIS(mat(), &isrows, &iscols);
		ISGetLocalSize(isrows, &lnrow);
		ISGetLocalSize(iscols, &lncol);
		ISGetIndices(isrows, &rind);
		ISGetIndices(iscols, &cind);
		mRowInd = std::vector<PetscInt>(rind, rind + lnrow);
		mColInd = std::vector<PetscInt>(cind, cind + lncol);
		ISRestoreIndices(isrows, &rind);
		ISRestoreIndices(iscols, &cind);
		ISDestroy(&isrows);
		ISDestroy(&iscols);
	};

	std::string type() const
	{
		MatType type;
		PetscErrorCode ierr = MatGetType(mat(), &type); CHKERR(ierr);
		return std::string(type);
	}

	double localmemory() const
	{
		MatInfo info;
		PetscErrorCode ierr = MatGetInfo(mat(), MAT_LOCAL, &info);
		CHKERR(ierr);
		return info.memory;
	}

	double globalmemory() const
	{
		MatInfo info;
		PetscErrorCode ierr = MatGetInfo(mat(), MAT_GLOBAL_SUM, &info);
		CHKERR(ierr);
		return info.memory;
	}

	std::string infostring() const
	{
		PetscErrorCode ierr;
		MatInfo info;
		MatType type;
		PetscInt rows_global, cols_global;
		PetscInt rows_local, cols_local;
		PetscInt rfrom, rto;

		ierr = MatGetType(mat(), &type); CHKERR(ierr);
		ierr = MatGetOwnershipRange(mat(), &rfrom, &rto); CHKERR(ierr);
		ierr = MatGetInfo(mat(), MAT_LOCAL, &info); CHKERR(ierr);
		ierr = MatGetSize(mat(), &rows_global, &cols_global); CHKERR(ierr);
		ierr = MatGetLocalSize(mat(), &rows_local, &cols_local); CHKERR(ierr);

		std::string s = "";
		s += strprint("\n");
		s += strprint("==== Matrix %s (%s) ================\n", cname(), type);
		s += strprint("Global rows x columns\t%7ld  x  %7ld\n", (int64_t)rows_global, (int64_t)cols_global);
		s += strprint("Local  rows x columns\t%7ld  x  %7ld\n", (int64_t)rows_local, (int64_t)cols_local);
		s += strprint("Local row range      \t%7ld  to %7ld\n", (int64_t)rfrom, (int64_t)rto - 1);
		s += strprint("block_size           \t%7ld\n", (int64_t)info.block_size);
		s += strprint("nz_allocated         \t%7ld\n", (int64_t)info.nz_allocated);
		s += strprint("nz_used              \t%7ld\n", (int64_t)info.nz_used);
		s += strprint("nz_unneeded          \t%7ld\n", (int64_t)info.nz_unneeded);
		s += strprint("memory               \t%7lf\n", info.memory);
		s += strprint("assemblies           \t%7ld\n", (int64_t)info.assemblies);
		s += strprint("mallocs              \t%7ld\n", (int64_t)info.mallocs);
		s += strprint("fill_ratio_given     \t%7lf\n", info.fill_ratio_given);
		s += strprint("fill_ratio_needed    \t%7lf\n", info.fill_ratio_needed);
		s += strprint("factor_mallocs       \t%7ld\n", (int64_t)info.factor_mallocs);
		s += strprint("=======================================\n");

		ierr = MatGetInfo(mat(), MAT_GLOBAL_SUM, &info); CHKERR(ierr);

		s += strprint("global nz_allocated  \t%7ld\n", (int64_t)info.nz_allocated);
		s += strprint("global nz_used       \t%7ld\n", (int64_t)info.nz_used);
		s += strprint("global nz_unneeded   \t%7ld\n", (int64_t)info.nz_unneeded);
		s += strprint("global memory        \t%7lf\n", info.memory);
		s += strprint("=======================================\n");
		s += strprint("\n");
		return s;
	}

	void printinfostring(int rank=0) const
	{
		std::string s = infostring();
		if (mpirank() == rank){
			printf(s.c_str());
			fflush(stdout);
		}
	}

	void printsize() const
	{
		printf("Matrix (%s) size = %d x %d\n", cname(), nglobalrows(), nglobalcols());
	}
	
	void printdistributions() const
	{
		for (int i = 0; i < mpisize(); i++){
			mpibarrier();
			if (i == mpirank()){
				printf("%s rank %d owns rows %d to %d  - total of %d rows\n", mpiprocname().c_str(), mpirank(), rowownership().start, rowownership().end - 1, nlocalrows());
				fflush(stdout);
			}			
		}
	}
	
	void printvalues(const std::string& format = "%7.3lf\t") const
	{
		mpibarrier();
		for (int p = 0; p<mpisize(); p++){
			mpibarrier();
			if (p == mpirank()){
				for (PetscInt i = rowownership().start; i<rowownership().end; i++){
					for (PetscInt j = 0; j<nglobalcols(); j++){
						double val;
						MatGetValue(mat(), i, j, &val);
						printf(format.c_str(), val);
					}
					printf("\n"); fflush(stdout);
				}
			}
		}
	}

	PetscErrorCode writetextfile(const std::string& filename)const
	{
		std::string mtype = type();
		if (mtype == "mpiaij" || mtype == "seqaij"){
			writetextfilesparse(filename);
		}
		else if (mtype == "mpidense" || mtype == "seqdense" || mtype == "elemental"){
			PetscViewer viewer;
			PetscViewerASCIIOpen(mpicomm(), filename.c_str(), &viewer);
			MatView(mat(), viewer);
			PetscViewerDestroy(&viewer);
		}
		return 0;
	}

	PetscErrorCode writetextfilesparse(const std::string& filename) const
	{
		for (int p = 0; p<mpisize(); p++){
			mpibarrier();
			if (p == mpirank()){
				FILE* fp;
				if (p == 0){
					fp = fileopen(filename, "w");
				}
				else{
					fp = fileopen(filename, "a");
				}

				bool writezeroatlowerright = true;
				cOwnership r = rowownership();
				for (PetscInt gi = r.start; gi<r.end; gi++){
					PetscInt nnz;
					const PetscInt* colind;
					const double* val;

					PetscErrorCode ierr = MatGetRow(mat(), gi, &nnz, &colind, &val); CHKERR(ierr);
					for (PetscInt j = 0; j<nnz; j++){
						fprintf(fp, "%d\t%d\t%20.16le\n", gi, colind[j], val[j]);
					}

					if (gi == nglobalrows() - 1){
						if (colind[nnz - 1] == nglobalcols() - 1){
							writezeroatlowerright = false;
						}
					}
					ierr = MatRestoreRow(mat(), gi, &nnz, &colind, &val); CHKERR(ierr);
				}
				//Always write a zero into last row/col if it is empty
				if (p == mpisize() - 1 && writezeroatlowerright){
					fprintf(fp, "%d\t%d\t%20.16le\n", nglobalrows() - 1, nglobalcols() - 1, 0.0);
				}

				fclose(fp);

			}
		}
		mpibarrier();
		return 0;
	}

	void vecmult(const cPetscDistVector& x, cPetscDistVector& b) const
	{
		if(nglobalcols() != x.globalsize()){
			printf("cPetscDistMatrix::vecmult() matrix %s and vector %s dimensions are not compatible\n", cname(), x.cname());
			throw(strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__));
		}

		if (nglobalrows() != b.globalsize()){
			printf("cPetscDistMatrix::vecmult() matrix %s and vector %s dimensions are not compatible\n", cname(), b.cname());
			throw(strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__));
		}
		
		PetscErrorCode ierr = MatMult(mat(), x.vec(), b.vec()); CHKERR(ierr);
	}

	cPetscDistVector vecmult(const cPetscDistVector& x) const
	{
	cPetscDistVector b(mpicomm(), nlocalrows(), nglobalrows());
	vecmult(x, b);
	return b;
	}

	void matmult(const cPetscDistMatrix& B, cPetscDistMatrix& C) const
	{
		PetscErrorCode ierr;
		if (nglobalcols() != B.nglobalrows()){
			printf("cPetscDistMatrix::matmult() matrix %s and matrix %s dimensions are not compatible\n", cname(), B.cname());
			throw(strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__));
		}
		ierr = MatMatMult(mat(), B.mat(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, C.ncmatptr()); CHKERR(ierr);
	}

	cPetscDistMatrix matmult(const cPetscDistMatrix& B) const
	{				
		cPetscDistMatrix C;
		matmult(B, C);
		return C;
	}
	
	void transvecmult(const cPetscDistVector& x, cPetscDistVector& b) const
	{		
	PetscErrorCode ierr = MatMultTranspose(mat(), x.vec(), b.ncvec()); CHKERR(ierr);
	}

	cPetscDistVector transvecmult(const cPetscDistVector& x) const
	{		
cPetscDistVector b(mpicomm(), nlocalcols(), nglobalcols());
transvecmult(x, b);
return b;
	}

	void transmatmult(const cPetscDistMatrix& B, cPetscDistMatrix& C) const
	{
		PetscErrorCode ierr;
		if (nglobalrows() != B.nglobalrows()){
			printf("cPetscDistMatrix::transmatmult() matrix %s and matrix %s dimensions are not compatible\n", cname(), B.cname());
			throw(strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__));
		}
		ierr = MatTransposeMatMult(mat(), B.mat(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, C.ncmatptr()); CHKERR(ierr);
	}
	
	cPetscDistMatrix transmatmult(const cPetscDistMatrix& B) const
	{		
		cPetscDistMatrix C;
		transmatmult(B,C);
		return C;		
	}

	void leftdiagonalscale(const cPetscDistVector& v) const
	{
		PetscErrorCode ierr = MatDiagonalScale(mat(), v.vec(),PETSC_NULL);	CHKERR(ierr);
		return;
	}

	double vtAv(const cPetscDistVector& v)
	{
		return v.dot(vecmult(v));
	}

	cPetscDistVector solve_CG(const cPetscDistVector& b, const cPetscDistVector& initialguess = cPetscDistVector())
	{
		//A x = b
		cPetscDistVector x;//solution vector that will be returned
		if (initialguess.globalsize() == 0){
			x.create("x",mpicomm(), nlocalcols(), nglobalcols(), 0.0);
		}
		else{
			x = initialguess;
		}

		cPetscDistVector ir = vecmult(x) - b;
		double inorm = ir.norm2();

		KSP ksp;
		PC  pc;
		PetscErrorCode ierr;
		ierr = KSPCreate(mpicomm(), &ksp); CHKERR(ierr);
#if PETSC_VERSION_LT(3,5,0)
		ierr = KSPSetOperators(ksp, mat(), mat(), SAME_NONZERO_PATTERN); CHKERR(ierr);		
#else
		ierr = KSPSetOperators(ksp, mat(), mat()); CHKERR(ierr);
#endif
		ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERR(ierr);
		ierr = KSPGetPC(ksp, &pc); CHKERR(ierr);
		ierr = KSPSetType(ksp, KSPCG); CHKERR(ierr);
		ierr = PCSetType(pc, PCJACOBI); CHKERR(ierr);

		//rtol    - the relative convergence tolerance (relative decrease in the residual norm)  
		//abstol  - the absolute convergence tolerance (absolute size of the residual norm)  
		//dtol    - the divergence tolerance (amount residual can increase before KSPDefaultConverged() concludes that the method is diverging)  
		//maxits  - maximum number of iterations to use  

		double reltol = 1e-6;
		double abstol = 1e-12;
		double divtol = 1e2;
		PetscInt maxits = nglobalcols();
		//PetscInt maxits = 400;


		ierr = KSPSetTolerances(ksp, reltol, abstol, divtol, maxits); CHKERR(ierr);

		//Point to function that monitors/outputs convergence progress at each CG iteration
		//ierr = KSPMonitorSet(ksp, kspmonitor, this, PETSC_NULL); CHKERR(ierr);
		//See also PETSC functions KSPMonitorDefault, KSPMonitorTrueResidualNorm, KSPMonitorTrueResidualMaxNorm
		//ierr = KSPMonitorSet(ksp, KSPMonitorDefault, PETSC_NULL, PETSC_NULL); CHKERR(ierr);

		//Setup and run the Conjugate Gradients algorithm
		ierr = KSPSetUp(ksp); CHKERR(ierr);
		ierr = KSPSolve(ksp, b.vec(), x.vec()); CHKERR(ierr);

		cPetscDistVector r = vecmult(x) - b;
		double anorm = r.norm2();

		//Output reason for convergence/divergence
		PetscInt niterations;
		double residnorm;
		ierr = KSPGetIterationNumber(ksp, &niterations); CHKERR(ierr);
		ierr = KSPGetResidualNorm(ksp, &residnorm); CHKERR(ierr);
		KSPConvergedReason reason;
		ierr = KSPGetConvergedReason(ksp, &reason); CHKERR(ierr);
		std::string reasonstring = std::string(KSPConvergedReasons[reason]);
		ierr = KSPDestroy(&ksp);

		if (mpirank() == 0){
			std::string s;
			s += strprint("CG %s ", reasonstring.c_str());
			s += strprint("%d Its ", niterations);
			s += strprint("rnorm = %3.1le ", residnorm);
			s += strprint("inorm = %3.1le ", inorm);
			s += strprint("anorm = %3.1le ", anorm);
			s += strprint("reduc = %3.1le ", anorm / inorm);
			s += strprint("\n");
			glog.logmsg(s.c_str());
		}

		return x;
	}

	//Operators
	cPetscDistMatrix& operator=(const cPetscDistMatrix& rhs)
	{
		printf("Entering cPetscDistMatrix& operator=(const cPetscDistMatrix& rhs)\n");
		PetscErrorCode ierr;
		deallocate();		
		pMat = new Mat;		
		ierr = MatDuplicate(rhs.mat(), MAT_COPY_VALUES, ncmatptr()); CHKERR(ierr);
		printf("Leaving cPetscDistMatrix& operator=(const cPetscDistMatrix& rhs)\n");
		return *this;
	}

	cPetscDistMatrix& operator+=(const cPetscDistMatrix& rhs)
	{
		PetscErrorCode ierr = MatAXPY(mat(), 1.0, rhs.mat(), DIFFERENT_NONZERO_PATTERN); CHKERR(ierr);
		return *this;
	}

	cPetscDistMatrix& operator*=(const double s)
	{
		PetscErrorCode ierr = MatScale(mat(), s); CHKERR(ierr);
		return *this;
	}

	cPetscDistVector operator*(const cPetscDistVector& x) const 
	{
	cPetscDistVector b = vecmult(x);
	return b;
	}
	
	cPetscDistVector operator^(const cPetscDistVector& x) const
	{
		//b = A'b;
	cPetscDistVector b = transvecmult(x);
	return b;
	}

	cPetscDistMatrix operator+(const cPetscDistMatrix& B) const 
	{		
		cPetscDistMatrix C(*this);
		C += B;
		return C;
	}

	cPetscDistMatrix operator*(const double& s) const 
	{
		cPetscDistMatrix B(*this);
		PetscErrorCode ierr = MatScale(B.mat(), s); CHKERR(ierr);
		return B;		
	}
	
	cPetscDistMatrix operator*(const cPetscDistMatrix& B) const
	{
		cPetscDistMatrix C;
		matmult(B, C);
		return C;
	}

	//Friends	
	friend cPetscDistMatrix operator*(const double& s, const cPetscDistMatrix& A)
	{
		cPetscDistMatrix B(A);
		PetscErrorCode ierr = MatScale(B.mat(), s); CHKERR(ierr);
		return B;
	}

};

class cPetscDistShellMatrix : public cPetscObject {

private:

	//Private data members		
	Mat* pMat = PETSC_NULL;

	const PetscObject* pobj() const{
		if (pMat) return (PetscObject*) pMat;
		else return (PetscObject*) PETSC_NULL;
	}

public:		
	bool _converged;
	KSPConvergedReason _conv_reason_code;
	std::string _conv_reason_str;
	PetscInt _niterations;	
	double _rnorm_start;
	double _rnorm_end;	
	double _solve_time;

	//Constructors & Destructors		
	cPetscDistShellMatrix(const std::string& name, const MPI_Comm& comm, const PetscInt& nlocalrows, const PetscInt& nlocalcols, const PetscInt& nglobalrows, const PetscInt& nglobalcols, void* context){
		create(name, comm, nlocalrows, nlocalcols, nglobalrows, nglobalcols, context);
	}

	~cPetscDistShellMatrix()
	{
		PetscErrorCode ierr = MatDestroy(ncmatptr()); CHKERR(ierr);
	}

	const Mat& mat() const { return *pMat; }
	Mat& ncmat() const { return *pMat; }
	const Mat* matptr() const { return pMat; }
	Mat* ncmatptr() const { return pMat; }

	//Member Functions			
	void create(const std::string& name, const MPI_Comm& comm, const PetscInt& _nlocalrows, const PetscInt& _nlocalcols, const PetscInt& _nglobalrows, const PetscInt& _nglobalcols, void* context)
	{		
		PetscErrorCode ierr;		
		pMat = new Mat;
		ierr = MatCreateShell(comm,
			_nlocalrows, _nlocalcols,
			_nglobalrows, _nglobalcols,
			context,ncmatptr()); CHKERR(ierr);		
		setname(name);
	}

	bool set_multiply_function_vec(void* func){
		PetscErrorCode ierr;
		ierr = MatShellSetOperation(mat(), MATOP_MULT,(void(*)(void))func); CHKERR(ierr);		
		return true;
	}

	
	PetscInt nglobalrows() const
	{
		PetscErrorCode ierr;
		PetscInt M, N;
		ierr = MatGetSize(mat(), &M, &N); CHKERR(ierr);
		return M;
	};

	PetscInt nglobalcols() const
	{
		PetscErrorCode ierr;
		PetscInt M, N;
		ierr = MatGetSize(mat(), &M, &N); CHKERR(ierr);
		return N;
	};

	PetscInt nlocalrows() const
	{
		PetscErrorCode ierr;
		PetscInt m, n;
		ierr = MatGetLocalSize(mat(), &m, &n); CHKERR(ierr);
		return m;
	};

	PetscInt nlocalcols() const
	{		
		PetscErrorCode ierr;
		PetscInt m, n;
		ierr = MatGetLocalSize(mat(), &m, &n); CHKERR(ierr);
		return n;
	};

	cPetscDistVector vecmult(const cPetscDistVector& x) const
	{		
		cPetscDistVector b(mpicomm(), nlocalcols(), nglobalcols());
		vecmult(x, b);
		return b;
	}

	void vecmult(const cPetscDistVector& x, cPetscDistVector& b) const
	{
	PetscErrorCode ierr = MatMult(mat(), x.vec(), b.vec()); CHKERR(ierr);
	}

	static PetscErrorCode kspmonitor(KSP ksp, PetscInt n, double rnorm, void *context)
	{				
		int k;
		if (n>1000)k = 1000;
		else if (n>100)k = 100;
		else if (n>10)k = 10;
		else k = 1;

		if (n%k == 0){
			cPetscDistShellMatrix* pA = (cPetscDistShellMatrix*)context;
			if (pA->mpirank() == 0){
				//printf("CG Iteration %6d residual norm = %8.6le\n", n, rnorm);
			}
			if (n == 0) pA->_rnorm_start = rnorm;
		}
		return 0;
	}

	cPetscDistVector solve_CG(const cPetscDistMatrix& P, const cPetscDistVector& b, const cPetscDistVector& initialguess = cPetscDistVector())
	{		
		//A x = b
		cPetscDistVector x;//solution vector that will be returned
		if (initialguess.globalsize() == 0){			
			x.create("x",mpicomm(), nlocalcols(), nglobalcols(), 0.0);
		}
		else{	
			x = initialguess;
		}

		KSP ksp;
		PC  pc;
		PetscErrorCode ierr;
		ierr = KSPCreate(mpicomm(), &ksp); CHKERR(ierr);
#if PETSC_VERSION_LT(3,5,0)
		ierr = KSPSetOperators(ksp, mat(), P.mat(), SAME_NONZERO_PATTERN); CHKERR(ierr);
#else
		ierr = KSPSetOperators(ksp, mat(), P.mat()); CHKERR(ierr);
#endif

		ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERR(ierr);
		ierr = KSPGetPC(ksp, &pc); CHKERR(ierr);
		ierr = KSPSetType(ksp, KSPCG); CHKERR(ierr);
		ierr = PCSetType(pc, PCJACOBI); CHKERR(ierr);

		//rtol    - the relative convergence tolerance (relative decrease in the residual norm)  
		//abstol  - the absolute convergence tolerance (absolute size of the residual norm)  
		//dtol    - the divergence tolerance (amount residual can increase before KSPDefaultConverged() concludes that the method is diverging)  
		//maxits  - maximum number of iterations to use  

		double reltol, abstol, divtol;
		PetscInt maxits;
		ierr = KSPGetTolerances(ksp, &reltol, &abstol, &divtol, &maxits); CHKERR(ierr);
		reltol = 1e-10;
		abstol = 1e-50;
		divtol = 1e4;
		maxits = nglobalcols();
		ierr = KSPSetTolerances(ksp, reltol, abstol, divtol, maxits); CHKERR(ierr);
		//ierr = KSPSetTolerances(ksp, PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT, maxits); CHKERR(ierr);

		//Point to function that monitors/outputs convergence progress at each CG iteration
		ierr = KSPMonitorSet(ksp, kspmonitor, this, PETSC_NULL); CHKERR(ierr);
		//See also PETSC functions KSPMonitorDefault, KSPMonitorTrueResidualNorm, KSPMonitorTrueResidualMaxNorm
		//ierr = KSPMonitorSet(ksp, KSPMonitorDefault, PETSC_NULL, PETSC_NULL); CHKERR(ierr);
				
		//Setup and run the Conjugate Gradients algorithm
		ierr = KSPSetUp(ksp); CHKERR(ierr);

		double t1 = gettime();
		ierr = KSPSolve(ksp, b.vec(), x.vec()); CHKERR(ierr);
		double t2 = gettime();
				
		ierr = KSPGetIterationNumber(ksp, &_niterations); CHKERR(ierr);
		ierr = KSPGetResidualNorm(ksp, &_rnorm_end); CHKERR(ierr);
		ierr = KSPGetConvergedReason(ksp, &_conv_reason_code); CHKERR(ierr);
		
		_solve_time = t2 - t1;
		_conv_reason_str = std::string(KSPConvergedReasons[_conv_reason_code]);
		if (_conv_reason_code <= 0){
			_converged = false;
		}
		else{
			_converged = true;
		}

		ierr = KSPDestroy(&ksp);		
		return x;
	}
	
	std::string convergence_summary(){		
		std::string s;
		if (_converged) s += strprint("CG Converged\n");
		else  s += strprint("CG Diverged\n");		
		s += strprint("\tCG Iterations=%d\n", _niterations);
		s += strprint("\tReason = %d (%s)\n", _conv_reason_code, _conv_reason_str.c_str());
		s += strprint("\tResidual norm at start = %8.6le\n", _rnorm_start);
		s += strprint("\tResidual norm at end = %8.6le\n", _rnorm_end);
		s += strprint("\tResidual norm reduction = %8.6le\n", _rnorm_end / _rnorm_start);
		s += strprint("\tSolve time = %8.5lf\n", _solve_time);
		s += strprint("\tTime per CG iteration = %8.5lf\n", _solve_time / (double)_niterations);
		return s;
	}

	//Operators
	cPetscDistVector operator*(const cPetscDistVector& x) const 
	{				
	cPetscDistVector b(mpicomm(), nlocalcols(), nglobalcols());
	vecmult(x, b);
	return b;
	}

};

#endif

