/*  fwdmodel_cest_devel.h - Development resting state ASL model

    Michael Chappell, IBME PUMMA & FMRIB Image Analysis Group

    Copyright (C) 2010 University of Oxford  */

/*  CCOPYRIGHT */

#include "fabber_core/fwdmodel.h"

#include <string>
using namespace std;

class CESTFwdModel : public FwdModel
{
public:
    static FwdModel *NewInstance();

    // Virtual function overrides
    virtual void Initialize(ArgsType &args);
    virtual void Evaluate(const NEWMAT::ColumnVector &params,
        NEWMAT::ColumnVector &result) const;
    virtual vector<string> GetUsage() const;
    virtual string ModelVersion() const;
    virtual void GetOptions(std::vector<OptionSpec> &opts) const;
    virtual std::string GetDescription() const;

    virtual void DumpParameters(const NEWMAT::ColumnVector &vec,
        const string &indents = "") const;

    virtual void NameParams(vector<string> &names) const;
    virtual int NumParams() const
    {
        return (3 * npool - 1) + 1 + (inferdrift ? 1 : 0) + (t12soft ? (2 * npool) : 0) + (3 * nexpool);
    }

    virtual ~CESTFwdModel() { return; }
    virtual void HardcodedInitialDists(MVNDist &prior, MVNDist &posterior) const;
    virtual void InitParams(MVNDist &posterior) const;

protected:
    //specific functions
    void Mz_spectrum(NEWMAT::ColumnVector &Mz, const NEWMAT::ColumnVector &wvec, const NEWMAT::ColumnVector &w1, const NEWMAT::ColumnVector &t, const NEWMAT::ColumnVector &M0, const NEWMAT::Matrix &wi, const NEWMAT::Matrix &kij, const NEWMAT::Matrix &T12) const;
    NEWMAT::ReturnMatrix Mz_spectrum_lorentz(const NEWMAT::ColumnVector &wvec, const NEWMAT::ColumnVector &w1, const NEWMAT::ColumnVector &t, const NEWMAT::ColumnVector &M0, const NEWMAT::Matrix &wi, const NEWMAT::Matrix &kij, const NEWMAT::Matrix &T12) const;

    void Mz_contribution_lorentz_simple(NEWMAT::ColumnVector &Mzc, const NEWMAT::ColumnVector &wvec, const double &I, const NEWMAT::ColumnVector &wi, const double &R) const;

    //maths functions
    void Ainverse(const NEWMAT::Matrix A, NEWMAT::RowVector &Ai) const;
    NEWMAT::ReturnMatrix expm(NEWMAT::Matrix inmatrix) const;
    NEWMAT::ReturnMatrix expm_eig(NEWMAT::Matrix inmatrix) const;
    NEWMAT::ReturnMatrix expm_pade(NEWMAT::Matrix inmatrix) const;
    NEWMAT::ReturnMatrix PadeApproximant(NEWMAT::Matrix inmatrix, int m) const;
    NEWMAT::ReturnMatrix PadeCoeffs(int m) const;

    // Constants

    // Lookup the starting indices of the parameters

    // vector indices for the parameters to expereicne ARD
    vector<int> ard_index;

    //flags
    bool t12soft;
    bool inferdrift;
    //bool pvcorr;
    bool lorentz;
    bool steadystate;
    bool setconcprior;

    // scan parameters
    vector<float> t;

    //model parameters
    int npool;
    float wlam;
    float B1set;
    NEWMAT::ColumnVector poolppm;
    NEWMAT::ColumnVector poolk;
    NEWMAT::Matrix T12master;
    NEWMAT::ColumnVector poolcon;
    NEWMAT::ColumnVector poolconprec;
    // NEWMAT::Matrix T12WMmaster;
    // NEWMAT::Matrix T12CSFmaster;

    // MT pool
    //bool mtpool;

    // extra pools
    int nexpool;
    NEWMAT::ColumnVector expoolppm;
    NEWMAT::ColumnVector expoolR;

    // Data specification
    NEWMAT::ColumnVector wvec;
    NEWMAT::ColumnVector w1vec;
    NEWMAT::ColumnVector tsatvec;

    //pulse specificaiton
    NEWMAT::ColumnVector pmagvec;
    NEWMAT::ColumnVector ptvec;
    int nseg;

    // processing flags
    mutable bool fastgrad; //use a fast approximation to the expm because we are caculating the gradient

    // ard flags
    bool doard;

private:
    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, CESTFwdModel> registration;
};
