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
    void GetOutputs(std::vector<std::string> &outputs) const;
    void EvaluateModel(
        const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result, const std::string &key = "") const;

    virtual string ModelVersion() const;
    virtual void GetOptions(std::vector<OptionSpec> &opts) const;
    virtual std::string GetDescription() const;

    virtual void DumpParameters(const NEWMAT::ColumnVector &vec, const string &indents = "") const;

    virtual void NameParams(vector<string> &names) const;
    virtual int NumParams() const
    {
        return (3 * npool - 1) + 1 + (inferdrift ? 1 : 0) + (t12soft ? (2 * npool) : 0) + (3 * nexpool)
            + (m_pvcorr ? 2 : 0);
    }

    virtual ~CESTFwdModel()
    {
        return;
    }
    virtual void HardcodedInitialDists(MVNDist &prior, MVNDist &posterior) const;
    virtual void InitParams(MVNDist &posterior) const;

protected:
    // specific functions

    void EvaluateCestRstar(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result, int pool_num) const;
    void RestrictPools(
        NEWMAT::ColumnVector &M0, NEWMAT::Matrix &wimat, NEWMAT::Matrix &kij, NEWMAT::Matrix &T12, int pool) const;
    void Evaluate(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result, int restrict_pool = -1) const;
    void Mz_spectrum(NEWMAT::ColumnVector &Mz, const NEWMAT::ColumnVector &wvec, const NEWMAT::ColumnVector &w1,
        const NEWMAT::ColumnVector &t, const NEWMAT::ColumnVector &M0, const NEWMAT::Matrix &wi,
        const NEWMAT::Matrix &kij, const NEWMAT::Matrix &T12) const;
    NEWMAT::ReturnMatrix Mz_spectrum_lorentz(const NEWMAT::ColumnVector &wvec, const NEWMAT::ColumnVector &w1,
        const NEWMAT::ColumnVector &t, const NEWMAT::ColumnVector &M0, const NEWMAT::Matrix &wi,
        const NEWMAT::Matrix &kij, const NEWMAT::Matrix &T12) const;

    void Mz_contribution_lorentz_simple(NEWMAT::ColumnVector &Mzc, const NEWMAT::ColumnVector &wvec, const double &I,
        const NEWMAT::ColumnVector &wi, const double &R) const;

    // Steady State Mz Spectrum
    void Mz_spectrum_SS(NEWMAT::ColumnVector &Mz, const NEWMAT::ColumnVector &wvec, const NEWMAT::ColumnVector &w1,
        const NEWMAT::ColumnVector &t, const NEWMAT::ColumnVector &M0, const NEWMAT::Matrix &wi,
        const NEWMAT::Matrix &kij, const NEWMAT::Matrix &T12, double w1EX, int pool_num) const;

    // Function to raise a matrix to a power
    inline NEWMAT::ReturnMatrix mpower(const NEWMAT::Matrix &Mat_Base, int Power) const;
    template <typename T> vector<T> spower(const vector<T> &Mat_Base, int Power) const;
    NEWMAT::ReturnMatrix spower_Mat(const NEWMAT::Matrix &, int) const;

    // Function to create a lineshape if an MT pool is present
    NEWMAT::ReturnMatrix absLineShape(const NEWMAT::ColumnVector &wvec, double T2) const;
    vector<double> SuperLorentzianGenerator(vector<double> &deltac, double T2) const;

    // maths functions
    void Ainverse(const NEWMAT::Matrix A, NEWMAT::RowVector &Ai) const;
    NEWMAT::ReturnMatrix expm(NEWMAT::Matrix inmatrix) const;
    NEWMAT::ReturnMatrix expm_eig(NEWMAT::Matrix inmatrix) const;
    NEWMAT::ReturnMatrix expm_pade(NEWMAT::Matrix inmatrix) const;
    NEWMAT::ReturnMatrix PadeApproximant(NEWMAT::Matrix inmatrix, int m, int &s) const;
    NEWMAT::ReturnMatrix PadeCoeffs(int m) const;

    // flags
    bool t12soft;
    bool inferdrift;
    bool use_b1off;

    // bool pvcorr;
    bool lorentz;
    bool steadystate;
    bool setconcprior;

    // scan parameters
    vector<float> t;

    // model parameters
    int npool;
    float wlam;
    float B1set;
    NEWMAT::ColumnVector poolppm;
    NEWMAT::ColumnVector poolk;
    NEWMAT::Matrix T12master;
    NEWMAT::ColumnVector poolcon;
    NEWMAT::ColumnVector poolconprec;

    // extra pools
    int nexpool;
    NEWMAT::ColumnVector expoolppm;
    NEWMAT::ColumnVector expoolR;

    // Data specification
    NEWMAT::ColumnVector wvec;
    NEWMAT::ColumnVector w1vec;
    NEWMAT::ColumnVector tsatvec;

    // pulse specificaiton
    NEWMAT::ColumnVector pmagvec;
    NEWMAT::ColumnVector ptvec;
    int nseg;

    // Readout specifications
    double m_tr;
    double m_fa;

    // Flags for New Steady State Sequence
    bool m_new_ss;           // Flag to use new steady state sequence
    bool m_inter_spoil;      // Flag to use interpulse saturation spoiling
    std::string m_lineshape; // String that describes which lineshape to use

    // Variables used to increase processing speed of SuperLorentzian Lineshape
    mutable double m_t2m;
    mutable std::vector<double> m_gc;

    // Partial volume correction
    bool m_pvcorr;
    NEWMAT::ColumnVector m_pv_img;
    double m_pv_threshold;
    double m_pv_csf_tiss_m0ratio;

private:
    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, CESTFwdModel> registration;
};
