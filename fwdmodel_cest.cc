/*  fwdmodel_cest_devel.cc - Developement CEST APT model

    Michael Chappell, IBME PUMMA & FMRIB Image Analysis Group

    Copyright (C) 2010 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_cest.h"

#include "miscmaths/miscprob.h"
#include "newimage/newimageall.h"
#include <iostream>
#include <newmatio.h>
#include <stdexcept>
using namespace NEWIMAGE;
#include "fabber_core/easylog.h"

FactoryRegistration<FwdModelFactory, CESTFwdModel> CESTFwdModel::registration("cest");

static OptionSpec OPTIONS[] = {
    { "spec", OPT_MATRIX, "Data specification", OPT_REQ, "" },
    { "pools", OPT_MATRIX, "Pool specification", OPT_REQ, "" },
    { "expools", OPT_MATRIX, "Extra pool specification", OPT_NONREQ, "" },
    { "ptrain", OPT_MATRIX, "Pulse specification", OPT_REQ, "" },
    { "t12prior", OPT_BOOL, "Include uncertainty in T1 and T2 values", OPT_NONREQ, "" },
    { "inferdrift", OPT_BOOL, "Infer drift parameter", OPT_NONREQ, "" },
    { "lorentz", OPT_BOOL, "Alternative solution to Bloch equations", OPT_NONREQ, "" },
    { "steadystate", OPT_BOOL, "Steady state", OPT_NONREQ, "" },
    { "pvimg", OPT_IMAGE,
        "Tissue partial volume image. Should be 3D image containing tissue partial volumes, i.e. sum of GM and WM "
        "partial volumes",
        OPT_NONREQ, "" },
    { "pv-threshold", OPT_FLOAT, "Partial volume threshold for including tissue contribution", OPT_NONREQ, "0.5" },
    { "csf-tiss-m0ratio", OPT_FLOAT, "Used for fixing CSF M0", OPT_NONREQ, "0.5269" },
    { "" },
};

void CESTFwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

std::string CESTFwdModel::GetDescription() const
{
    return "Model for Chemical Exchange Saturation transfer, with correction for partial volume effects using a tissue "
           "PV map";
}

string CESTFwdModel::ModelVersion() const
{
    string version = "fwdmodel_cest.cc";
#ifdef GIT_SHA1
    version += string(" Revision ") + GIT_SHA1;
#endif
#ifdef GIT_DATE
    version += string(" Last commit ") + GIT_DATE;
#endif
    return version;
}

void CESTFwdModel::HardcodedInitialDists(MVNDist &prior, MVNDist &posterior) const
{
    assert(prior.means.Nrows() == NumParams());

    SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;

    // Set priors

    int place = 1;
    // M0
    prior.means.Rows(1, npool) = 0.0;
    precisions(place, place) = 1e-12;
    place++;

    if (npool > 1)
    {
        for (int i = 2; i <= npool; i++)
        {
            if (setconcprior)
            {
                // priors have been specified via the poolmat
                prior.means(place)
                    = poolcon(place - 1); // NB poolcon vec doesn't have water in so entry 1 is pool 2 etc
                precisions(place, place) = poolconprec(place - 1);
            }
            else
            {
                // hardcoded default of zero mean and (relatively) uniformative precision
                // prior mean of zero set above
                precisions(place, place) = 1e2;
                place++;
            }
        }
    }

    // exchnage consts (these are log_e)
    if (npool > 1)
    {
        for (int i = 2; i <= npool; i++)
        {
            prior.means(place) = poolk(i - 1); // NOTE these shoudl already be log in poolk
            precisions(place, place) = 1;
            place++;
        }
    }

    // frequency offsets (ppm)
    prior.means(place) = 0; // water centre offset
    precisions(place, place) = 10;
    place++;

    if (npool > 1)
    {
        for (int i = 2; i <= npool; i++)
        {
            prior.means(place) = poolppm(i - 1);
            precisions(place, place) = 1e12; // 1e6;
            place++;
        }
    }

    // B1 offset (fractional)
    prior.means(place) = 0;
    precisions(place, place) = 1e99; // 1e20; //1e6;
    place++;

    if (inferdrift)
    {
        // Drift (ppm/sample)
        prior.means(place) = 0;
        precisions(place, place) = 1e17;
        place++;
    }

    if (t12soft)
    {
        // T1 values
        for (int i = 1; i <= npool; i++)
        {
            prior.means(place) = T12master(1, i);
            precisions(place, place) = 44.4; // all T1s have same prior uncertainty
            place++;
        }

        // T12 values
        for (int i = 1; i <= npool; i++)
        {
            float T12 = T12master(2, i);
            prior.means(place) = T12;
            precisions(place, place) = 1 / std::pow(T12 / 5, 2); // prior has std dev of 1/5 of the
                                                                 // value, to try and get the
                                                                 // scaling about right
            place++;
        }
    }

    if (use_pvcorr)
    {
        // CSF pool T1 and T2 priors
        prior.means(place) = 1.9;
        precisions(place, place) = 44.4;
        place++;
        prior.means(place) = 0.25;
        precisions(place, place) = 1 / std::pow(prior.means(place) / 5, 2);
        place++;
    }

    // Extra ('indepdnent') pools
    if (nexpool > 0)
    {
        for (int i = 1; i <= nexpool; i++)
        {
            prior.means(place) = 1;
            precisions(place, place) = 1;
            place++;
            prior.means(place) = expoolppm(i);
            precisions(place, place) = 1e12;
            place++;
            prior.means(place) = expoolR(i);
            precisions(place, place) = 1;
        }
    }

    // Set precsions on priors
    prior.SetPrecisions(precisions);

    // Set initial posterior
    posterior = prior;

    // For parameters with uniformative prior chosoe more sensible inital posterior
    // now done in initialise!
    //  posterior.means(1) = 1000;
    //  precisions(1,1) = 10;
    // posterior.SetPrecisions(precisions);
}

void CESTFwdModel::InitParams(MVNDist &posterior) const
{
    // Check out dataspec is the right size
    int nt = data.Nrows();
    if (wvec.Nrows() != nt)
    {
        throw InvalidOptionValue(
            "dataspec", "", "Incorrect number of frequencies - should match number of data volumes");
    }

    // load the existing precisions as the basis for any update
    SymmetricMatrix precisions;
    precisions = posterior.GetPrecisions();

    // init the M0a value  - to max value in the z-spectrum
    posterior.means(1) = data.Maximum();
    precisions(1, 1) = 10;

    // init the pool concentraitons
    if (npool > 1)
    {
        for (int i = 2; i <= npool; i++)
        {
            // NB poolcon vec doesn't have water in so entry 1 is pool 2 etc
            posterior.means(i) = poolcon(i - 1);
            precisions(i, i) = 1e6;
        }
    }

    // init the ppmoff value - by finding the freq where the min of z-spectrum is
    int ind;
    float val;
    val = data.Minimum1(ind);     // find the minimum in the z-spectrum
    val = wvec(ind) * 1e6 / wlam; // frequency of the minimum in ppm
    if (val > 0.5)
        val = 0.5; // put a limit on the value
    if (val < -0.5)
        val = -0.5;
    int ppmind = 2 * npool;
    posterior.means(ppmind) = val;

    posterior.SetPrecisions(precisions);
}

void CESTFwdModel::GetOutputs(std::vector<std::string> &outputs) const
{
    for (int p = 1; p < npool; p++)
    {
        char pool_char = char(int('a') + p);
        outputs.push_back(string("cest_rstar_") + pool_char);
    }
}

double lin_interp(const ColumnVector &x, const ColumnVector &y, double pos)
{
    // Quick-and-dumb linear interpolation. Assume x, pos > 0
    // and function starts at 0, 0
    double prev_val = 0;
    double prev_x = 0;
    for (int i = 1; i <= x.Nrows(); i++)
    {
        if (pos < x(i))
        {
            double frac = (pos - prev_x) / (x(i) - prev_x);
            return prev_val + frac * (y(i) - prev_val);
        }
        prev_val = y(i);
        prev_x = x(i);
    }

    // pos beyond last x value - just return last y value
    return y(y.Nrows());
}

void CESTFwdModel::EvaluateModel(const ColumnVector &params, ColumnVector &result, const std::string &key) const
{
    if (key == "")
    {
        Evaluate(params, result);
    }
    else
    {
        // Outputting  CEST R* - need to know pool number which is encoded by the letter
        // at the end of the key (a=water, b=pool 2, etc). This will always be > 1
        int pool_num = 1 + int(key[key.length() - 1]) - int('a');
        assert(pool_num > 1);
        EvaluateCestRstar(params, result, pool_num);
    }
}

void CESTFwdModel::EvaluateCestRstar(const ColumnVector &params, ColumnVector &result, int pool_num) const
{
    assert(pool_num > 1);
    // LOG << "Outputting CESTRstar for pool " << pool_num << endl;
    // For this calculation, we use default T1/T2, NOT any inferred values or
    // image priors
    ColumnVector mod_params = params;
    if (t12soft)
    {
        int t12_idx = (npool - 1) * 3 + 3 + (inferdrift ? 1 : 0);
        for (int i = 1; i <= npool; i++)
        {
            mod_params(t12_idx + i) = T12master(1, i);
            mod_params(t12_idx + npool + i) = T12master(2, i);
        }
    }
    // We also do not use a water ppm offset
    int ppm_off_idx = (npool - 1) * 2 + 2;
    mod_params(ppm_off_idx) = 0;

    ColumnVector water_only, with_pool;
    Evaluate(mod_params, water_only, 1);
    Evaluate(mod_params, with_pool, pool_num);
    // LOG << "freq: " << wvec.t();
    // LOG << "water only: " << water_only.t();
    // LOG << "With pool: " << with_pool.t();
    // We evaluate the spectrum at a fixed PPM from the poolmat file, unless
    // the value is 0 in which case we use 50ppm (works for semisolid pool)
    double ppm_eval = 50;
    if (poolppm(pool_num - 1) != 0)
    {
        ppm_eval = poolppm(pool_num - 1);
    }
    // LOG << "Evaluating at " << ppm_eval << ", " << (ppm_eval* wlam / 1e6) << endl;
    // Evaluate at fixed PPM by linear interpolation. Note freq transformation
    // same as transformation applied to wvec
    double water = lin_interp(wvec, water_only, ppm_eval * wlam / 1e6);
    double pool = lin_interp(wvec, with_pool, ppm_eval * wlam / 1e6);
    // LOG << "water " << water << endl;
    // LOG << "pool " << pool << endl;
    // LOG << "frac " << pool/water << endl;
    result.ReSize(1);
    result(1) = 100 * (water - pool) / params(1);
    // LOG << "res " << result.t() << endl;
}

/**
 * To evaluate only with single pool or 2-pool, need to restrict matrices to the relevant pools
 *
 * This is used in order to calculate CEST R* in the final output stage
 *
 * wvec - no change, this is vector of size nsamp
 * w1   - no change, this is vector of size nsamp
 * tsatvec - no change, this is vector of size nsamp
 * M0 - need to restrict to rows 1 and poolnum
 * wimat - need to restric to rows 1 and poolnum
 * kij - need to restict to rows/columns 1 and poolnum
 * T12 - need to restrict to rows 1 and poolnum
 */
void CESTFwdModel::RestrictPools(ColumnVector &M0, Matrix &wimat, Matrix &kij, Matrix &T12, int pool) const
{
    if (pool == 1)
    {
        // Restrict solution to water pool only
        ColumnVector M0_res(1);
        M0_res(1) = M0(1);
        M0 = M0_res;

        Matrix wimat_res(1, wimat.Ncols());
        wimat_res.Row(1) = wimat.Row(1);
        wimat = wimat_res;

        Matrix kij_res(1, 1);
        kij_res(1, 1) = kij(1, 1);
        kij = kij_res;

        Matrix T12_res(2, 1);
        T12_res.Column(1) = T12.Column(1);
        T12 = T12_res;
    }
    else if (pool > 1)
    {
        // Restrict solution to water pool + other specified pool
        ColumnVector M0_res(2);
        M0_res(1) = M0(1);
        M0_res(2) = M0(pool);
        M0 = M0_res;

        Matrix wimat_res(2, wimat.Ncols());
        wimat_res.Row(1) = wimat.Row(1);
        wimat_res.Row(2) = wimat.Row(pool);
        wimat = wimat_res;

        Matrix kij_res(2, 2);
        kij_res(1, 1) = kij(1, 1);
        kij_res(1, 2) = kij(1, pool);
        kij_res(2, 1) = kij(pool, 1);
        kij_res(2, 2) = kij(pool, pool);
        kij = kij_res;

        Matrix T12_res(2, 2);
        T12_res.Column(1) = T12.Column(1);
        T12_res.Column(2) = T12.Column(pool);
        T12 = T12_res;
    }
}

void CESTFwdModel::Evaluate(const ColumnVector &params, ColumnVector &result, int restrict_pool) const
{
    // ensure that values are reasonable
    // negative check
    ColumnVector paramcpy = params;
    for (int i = 1; i <= NumParams(); i++)
    {
        if (params(i) < 0)
        {
            paramcpy(i) = 0;
        }
    }

    // model matrices
    ColumnVector M0(npool);
    Matrix kij(npool, npool);
    Matrix T12(2, npool);
    ColumnVector ppmvec(npool);
    // MT pool
    ColumnVector mtparams(5);
    // Extra pools
    ColumnVector exI(nexpool);
    ColumnVector exppmvec(nexpool);
    ColumnVector exR(nexpool);

    // extract values from params
    // M0 comes first
    int place = 1;
    M0(1) = paramcpy(place); // this is the 'master' M0 value of water
    if (M0(1) < 1e-4)
        M0(1) = 1e-4; // M0 of water cannot disapear all together
    place++;

    // values in the parameters are ratios of M0_water
    float M0ratio;
    for (int j = 2; j <= npool; j++)
    {
        M0ratio = paramcpy(place);
        if (M0ratio > 1.0)
        { // dont expect large ratios
            M0ratio = 1.0;
        }
        M0(j) = M0ratio * M0(1);

        place++;
    }

    // now exchange - we assume that only significant exchnage occurs with water
    kij = 0.0; // float ktemp;
    for (int j = 2; j <= npool; j++)
    {
        kij(j, 1) = exp(params(place)); // non-linear transformation
        if (kij(j, 1) > 1e6)
            kij(j, 1) = 1e6; // exclude really extreme values
        kij(1, j) = kij(j, 1) * M0(j) / M0(1);
        place++;
    }

    // frequency offset next
    ppmvec = params.Rows(place, place + npool);
    place += npool;

    // now B1 offset
    // scale the B1_off parameter to achieve proper updating of this parameter
    float B1off = params(place) * 1e6;
    // (otherwise it is sufficently small that numerical diff is ineffective)
    place++;

    // Drift
    float drift = 0;
    if (inferdrift)
    {
        drift = params(place) * 1e6; // scale this parameters like B1 offset
        place++;
    }

    // T12 parameter values
    if (t12soft)
    {
        // T12 values
        for (int i = 1; i <= npool; i++)
        {
            T12(1, i) = paramcpy(place);
            if (T12(1, i) < 1e-12)
                T12(1, i) = 1e-12; // 0 is no good for a T1 value
            if (T12(1, i) > 10)
                T12(1, i) = 10; // Prevent convergence issues causing T1 to blow up
            place++;
        }
        for (int i = 1; i <= npool; i++)
        {
            T12(2, i) = paramcpy(place);
            if (T12(2, i) < 1e-12)
                T12(2, i) = 1e-12; // 0 is no good for a T2 value

            if (T12(2, i) > 1)
                T12(2, i) = 1; // Prevent convergence issues causing T2 to blow up
            place++;
        }
    }
    else
    {
        T12 = T12master;
    }

    // PVC parameters
    double t1csf, t2csf;
    if (use_pvcorr)
    {
        t1csf = paramcpy(place++);
        t2csf = paramcpy(place++);
    }

    // Extra pools
    if (nexpool > 0)
    {
        for (int i = 1; i <= nexpool; i++)
        {
            exI(i) = params(place);
            place++;
            exppmvec(i) = params(place);
            place++;
            exR(i) = params(place);
            place++;
        }
    }

    // MODEL - CEST

    // deal with frequencies
    int nsamp = wvec.Nrows();
    Matrix wimat(npool, nsamp);

    for (int j = 1; j <= nsamp; j++)
    {
        wimat(1, j) = wlam * (ppmvec(1) + (j - 1) * drift) / 1e6;
    }
    if (npool > 1)
    {
        for (int i = 2; i <= npool; i++)
        {
            for (int j = 1; j <= nsamp; j++)
            {
                wimat(i, j)
                    = wlam * ppmvec(i) / 1e6 + wlam * (ppmvec(1) + (j - 1) * drift) / 1e6 * (1 + ppmvec(i) / 1e6);
            }
        }
    }
    // species b is at ppm*wlam, but also include offset of main field

    // deal with B1
    // B1 cannot go too small
    if (B1off < -0.5)
        B1off = -0.5;
    // unlikely to get this big (hardlimit in case of convergence problems)
    if (B1off > 10)
        B1off = 10;
    ColumnVector w1 = w1vec * (1 + B1off); // w1 in radians!

    // frequencies for the extra pools
    Matrix exwimat(nexpool, nsamp);
    if (nexpool > 0)
    {
        for (int i = 1; i <= nexpool; i++)
        {
            for (int j = 1; j <= nsamp; j++)
            {
                exwimat(i, j)
                    = wlam * ppmvec(i) / 1e6 + wlam * (ppmvec(1) + (j - 1) * drift) / 1e6 * (1 + exppmvec(i) / 1e6);
            }
        }
    }

    if (restrict_pool > 0)
    {
        // We are restriction the solution to the water pool + possibly one other pool
        RestrictPools(M0, wimat, kij, T12, restrict_pool);
    }
    else
    {
        float pv_val = 1.0;
        if (use_pvcorr)
        {
            pv_val = tissue_pv(voxel);
        }

        // Generate tissue z spectrum if pv > threshold (otherwise we're doing csf-only fit)
        ColumnVector Mztissue(wvec);

        if (pv_val >= pv_threshold)
        {
            // PV correction supports lorentz option for 'tissue' spectrum
            if (lorentz)
            {
                Mztissue = Mz_spectrum_lorentz(wvec, w1, tsatvec, M0, wimat, kij, T12);
            }
            else
            {
                // The complete tissue spectrum
                Mz_spectrum(Mztissue, wvec, w1, tsatvec, M0, wimat, kij, T12);
            }
        }

        if (use_pvcorr)
        {
            // Now repeat but for one pool to get CSF Lorentzian
            ColumnVector M0csf(1);
            M0csf(1) = M0(1) * csf_tiss_m0ratio;

            if (M0csf(1) < 1e-4)
            {
                M0csf(1) = 1e-4;
            }

            Matrix T12csf(2, 1);
            T12csf(1, 1) = t1csf;
            T12csf(2, 1) = t2csf;
            Matrix kij_csf(1, 1);
            kij_csf(1, 1) = 0.0;

            Matrix wimat_csf(1, nsamp);
            for (int j = 1; j <= nsamp; j++)
            {
                wimat_csf(1, j) = wlam * (ppmvec(1) + (j - 1) * drift) / 1e6;
            }

            // a csf spectrum is always generated
            ColumnVector Mzcsf(wvec);
            Mz_spectrum(Mzcsf, wvec, w1, tsatvec, M0csf, wimat_csf, kij_csf, T12csf);

            if (pv_val >= pv_threshold)
            {
                result = Mztissue * pv_val + Mzcsf * (1.0 - pv_val);
            }
            else
            {
                result = Mzcsf;
            }
        }
        else
        {
            result = Mztissue;
        }
    }

    // extra pools
    ColumnVector Mz_extrapools(nsamp);
    Mz_extrapools = 0.0;
    ColumnVector Mzc(nsamp);
    Mzc = 0.0;
    if (nexpool > 0)
    {
        for (int i = 1; i <= nexpool; i++)
        {
            Mz_contribution_lorentz_simple(
                Mzc, wvec, exI(i), (exwimat.SubMatrix(i, i, 1, nsamp)).t(), exR(i));
            Mz_extrapools += Mzc;
        }
    }

    result = result - M0(1) * Mz_extrapools;
}

FwdModel *CESTFwdModel::NewInstance()
{
    return new CESTFwdModel();
}

void CESTFwdModel::Initialize(ArgsType &args)
{
    // Matrix containing spec for each datapoint
    // 3 columns: Freq (ppm), B1 (T), tsat (s)
    // Nrows = number data points
    Matrix dataspec;

    // Matrix containing setup information for the pools
    // 4 columns: res. freq (rel. to water) (ppm), rate (species-> water), T1, T2.
    // 1st row is water, col 1 interpreted as the actaul centre freq of water if > 0, col 2 ignored.
    // Further rows for pools to be modelled.
    Matrix poolmat;

    // Matrix for the 'extra' pools
    Matrix expoolmat;

    string expoolmatfile;
    string pulsematfile;

    t12soft = false;

    // read data specification from file
    dataspec = read_ascii_matrix(args.Read("spec"));

    // read pool specification from file
    poolmat = read_ascii_matrix(args.Read("pools"));

    // read extra pool specification from file
    expoolmatfile = args.ReadWithDefault("expools", "none");

    // read pulsed saturation specification
    pulsematfile = args.ReadWithDefault("ptrain", "none");

    t12soft = args.ReadBool("t12prior");
    inferdrift = args.ReadBool("inferdrift");

    // alternatives to Matrix exponential solution to Bloch equations
    lorentz = args.ReadBool("lorentz"); // NB only compatible with single pool
    steadystate = args.ReadBool("steadystate");

    // Deal with the specification of the pools
    npool = poolmat.Nrows();
    if (poolmat.Ncols() < 4 || poolmat.Ncols() > 6)
    {
        throw invalid_argument("Incorrect number of columns in pool spefication file");
    }
    // water centre
    float wdefault = 42.58e6 * 3 * 2 * M_PI; // the default centre freq. (3T)
    if (poolmat(1, 1) > 0)
        wlam = poolmat(1, 1) * 2 * M_PI;
    else
        wlam = wdefault;
    // ppm ofsets
    poolppm = poolmat.SubMatrix(2, npool, 1, 1);
    // exchange rate
    poolk = log(poolmat.SubMatrix(2, npool, 2, 2)); // NOTE the log_e transformation
    // T1 and T2 values
    T12master = (poolmat.SubMatrix(1, npool, 3, 4)).t();

    if (poolmat.Ncols() > 4)
    {
        // pool ratios (concetrations) - these are used for initialisaiton, but may also be used as
        // prior means if precisions are provided
        poolcon = poolmat.SubMatrix(2, npool, 5, 5); // ignore water pool
    }
    else
    {
        poolcon.ReSize(npool);
        poolcon = 0;
    }

    setconcprior = false;
    if (poolmat.Ncols() > 5)
    {
        // pool ratio precisions - this overrides the inbuilt priors for the concentrations using
        // these precision and the means in poolcon
        poolconprec = poolmat.SubMatrix(2, npool, 6, 6); // ingore water pool
        setconcprior = true;
    }

    // check that the method chosen is possible
    if ((npool > 1) & lorentz)
        throw invalid_argument("Lorentzian (analytic) solution only compatible with single pool");

    if (expoolmatfile != "none")
    {
        // process the Extra pool specification matrix
        expoolmat = read_ascii_matrix(expoolmatfile);
        nexpool = expoolmat.Nrows();
        // cout << "nexpool: " << nexpool << endl;
        if (expoolmat.Ncols() != 2)
            throw invalid_argument("Incorrect number of columns in extra pool spefication file");
        // ppm ofsets
        expoolppm = expoolmat.SubMatrix(1, nexpool, 1, 1);
        // (prior) R values
        expoolR = expoolmat.SubMatrix(1, nexpool, 2, 2);
    }
    else
    {
        nexpool = 0;
    }

    // setup vectors that specify deatils of each data point
    // sampling frequency
    wvec = dataspec.Column(1) * wlam / 1e6;
    // B1 value, convert to radians equivalent
    w1vec = dataspec.Column(2) * 42.58e6 * 2 * M_PI;
    // Saturation time
    tsatvec = dataspec.Column(3);

    LOG << " Model parameters: " << endl;
    LOG << " Water - freq. (MHz) = " << wlam / 2 / M_PI << endl;
    LOG << "         T1    (s)   = " << T12master(1, 1) << endl;
    LOG << "         T2    (s)   = " << T12master(2, 1) << endl;
    ;
    for (int i = 2; i <= npool; i++)
    {
        LOG << " Pool " << i << " - freq. (ppm)  = " << poolppm(i - 1) << endl;
        LOG << "       - kiw   (s^-1) = " << exp(poolk(i - 1)) << endl;
        LOG << "       - T1    (s)    = " << T12master(1, i) << endl;
        LOG << "       - T2    (s)    = " << T12master(2, i) << endl;
    }

    LOG << "Sampling frequencies (ppm):" << endl
        << (dataspec.Column(1)).t() << endl
        << "                     (rad/s):" << endl
        << wvec.t() << endl;
    LOG << "B1 values (uT):" << endl
        << (dataspec.Column(2)).t() * 1e6 << endl
        << "          (rad/s):" << endl
        << w1vec.t() << endl;

    // Pulsed saturation modelling
    // For pulsed saturation the B1 values are taken as peak values
    if (pulsematfile == "none")
    {
        cout << "WARNGING! - you should supply a pulsemat file, in future the calling of this "
                "model without one (to get continuous saturation) will be depreceated"
             << endl;
        nseg = 1;
        pmagvec.ReSize(1);
        ptvec.ReSize(1);
        pmagvec = 1.0;
        ptvec = 1e12; // a very long value for continuous saturation

        LOG << "Saturation times (s):" << endl << tsatvec.t() << endl;
    }
    else
    {
        Matrix pulsemat;
        pulsemat = read_ascii_matrix(pulsematfile);
        // vector of (relative) magnitude values for each segment
        pmagvec = pulsemat.Column(1);
        // vector of time durations for each segment
        // note that we are loding in a vector of times for the end of each segment
        ColumnVector pttemp = pulsemat.Column(2);
        ColumnVector nought(1);
        nought = 0.0;
        ptvec = pttemp - (nought & pttemp.Rows(1, pttemp.Nrows() - 1));
        nseg = pulsemat.Nrows();

        LOG << "Pulse repeats:" << endl << tsatvec.t() << endl;

        LOG << "Pulse shape:" << endl << "Number of segments: " << nseg << endl;
        LOG << " Magnitudes (relative): " << pmagvec.t() << endl;
        LOG << " Durations (s): " << ptvec.t() << endl;
    }

    // Partial volume correction
    try
    {
        // Partial volume image should be 3D map
        Matrix pvimg = args.GetVoxelData("pvimg");
        LOG << "Tissue partial volume image found - using partial volume correction" << endl;
        use_pvcorr = true;
        if (pvimg.Nrows() > 1)
        {
            throw InvalidOptionValue("pvimg", "4D data set", "3D data set");
        }
        tissue_pv = pvimg.Row(1).t();

        pv_threshold = args.GetDoubleDefault("pv-threshold", 0.5, 0, 1);
        csf_tiss_m0ratio = args.GetDoubleDefault("csf-tiss-m0ratio", 0.5269);
    }
    catch (DataNotFound &e)
    {
        use_pvcorr = false;
        LOG << "Tissue partial volume image not found - will not use partial volume correction" << endl;
    }
}

void CESTFwdModel::DumpParameters(const ColumnVector &vec, const string &indent) const
{
    // cout << vec.t() << endl;
}

void CESTFwdModel::NameParams(vector<string> &names) const
{
    names.clear();

    // name the parameters for the pools using letters
    string lettervec[] = { "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r",
        "s", "t", "u", "v", "w", "x", "y", "z" };

    names.push_back("M0a");
    if (npool > 1)
    {
        for (int i = 2; i <= npool; i++)
        {
            names.push_back("M0" + lettervec[i - 1] + "_r");
        }
    }
    if (npool > 1)
    {
        for (int i = 2; i <= npool; i++)
        {
            names.push_back("k" + lettervec[i - 1] + "a");
        }
    }
    names.push_back("ppm_off");
    if (npool > 1)
    {
        for (int i = 2; i <= npool; i++)
        {
            names.push_back("ppm_" + lettervec[i - 1]);
        }
    }
    names.push_back("B1_off");
    if (inferdrift)
    {
        names.push_back("drift");
    }

    if (t12soft)
    {
        for (int i = 1; i <= npool; i++)
        {
            names.push_back("T1" + lettervec[i - 1]);
        }
        for (int i = 1; i <= npool; i++)
        {
            names.push_back("T2" + lettervec[i - 1]);
        }
    }

    if (use_pvcorr)
    {
        // Additional PV paramters
        names.push_back("T1csf");
        names.push_back("T2csf");
    }

    if (nexpool > 0)
    {
        for (int i = 1; i <= nexpool; i++)
        {
            names.push_back("I_" + lettervec[npool - 1 + i]);
            names.push_back("ppm_" + lettervec[npool - 1 + i]);
            names.push_back("R_" + lettervec[npool - 1 + i]);
        }
    }
}

ReturnMatrix CESTFwdModel::expm(Matrix inmatrix) const
{
    // to set the routine we use to do expm
    return expm_pade(inmatrix);
}

/**
 * Do matrix exponential using eigen decomposition of the matrix
 */
ReturnMatrix CESTFwdModel::expm_eig(Matrix inmatrix) const
{
    // a bit poor - the matrix coming in should be symmetric, but I haven't
    // implemented this elsewhere in the code (yet!)
    SymmetricMatrix A;
    A << inmatrix;

    // eigen decomposition
    DiagonalMatrix D;
    Matrix V;
    Jacobi(A, D, V);

    // exponent of D
    for (int i = 1; i <= D.Nrows(); i++)
    {
        // if (D(i)>100) D(i)=100;
        D(i) = exp(D(i));
    }

    // now matrix exponential
    Matrix EXP(A);
    EXP = V * D * V.t();
    return EXP;
}

/**
 * Do matrix exponential
 *
 * Algorithm from Higham, SIAM J. Matrix Analysis App. 24(4) 2005, 1179-1193
 */
ReturnMatrix CESTFwdModel::expm_pade(Matrix inmatrix) const
{
    Matrix A = inmatrix;
    Matrix X(A.Nrows(), A.Ncols());

    int metflag = 0;
    int mvals[5] = { 3, 5, 7, 9, 13 };
    float thetam[5] = { 1.495585217958292e-002, 2.539398330063230e-001, 9.504178996162932e-001, 2.097847961257068e+000,
        5.371920351148152e+000 };
    int i = 0;

    while (metflag < 1)
    {
        if (i < 4)
        {
            if (A.Norm1() <= thetam[i])
            {
                X = PadeApproximant(A, mvals[i]);
                metflag = 1;
            }
        }
        else
        {
            // Using Pade Approximant degree 13 and scaling and squaring
            int s = ceil(log2(A.Norm1() / thetam[4]));
            float half = 0.5;
            A *= MISCMATHS::pow(half, s);
            X = PadeApproximant(A, 13);
            for (int i = 1; i <= s; i++)
            {
                X *= X;
            }
            metflag = 1;
        }
        i++;
    }

    return X;
}

ReturnMatrix CESTFwdModel::PadeApproximant(Matrix inmatrix, int m) const
{
    assert(inmatrix.Nrows() == inmatrix.Ncols());
    int n = inmatrix.Nrows();
    IdentityMatrix I(n);
    ColumnVector coeff = PadeCoeffs(m);
    Matrix X(n, n);
    Matrix U(n, n);
    Matrix V(n, n);
    vector<Matrix> Apowers;

    switch (m)
    {
    case 3:
    case 5:
    case 7:
    case 9:
        Apowers.push_back(I);
        Apowers.push_back(inmatrix * inmatrix);
        for (int j = 3; j <= ceil((m + 1 / 2)); j++)
        {
            Apowers.push_back(Apowers[j - 2] * Apowers[1]);
        }
        U = 0.0;
        V = 0.0;
        for (int j = m + 1; j >= 2; j -= 2)
        {
            U += coeff(j) * Apowers[j / 2 - 1];
        }
        U = inmatrix * U;
        for (int j = m; j >= 1; j -= 2)
        {
            V += coeff(j) * Apowers[(j + 1) / 2 - 1];
        }
        X = (U + V) * (-U + V).i();
        break;
    case 13:
        Matrix A2(n, n);
        Matrix A4(n, n);
        Matrix A6(n, n);
        A2 = inmatrix * inmatrix;
        A4 = A2 * A2;
        A6 = A2 * A4;
        U = inmatrix
            * (A6 * (coeff(14) * A6 + coeff(12) * A4 + coeff(10) * A2) + coeff(8) * A6 + coeff(6) * A4 + coeff(4) * A2
                  + coeff(2) * I);
        V = A6 * (coeff(13) * A6 + coeff(11) * A4 + coeff(9) * A2) + coeff(7) * A6 + coeff(5) * A4 + coeff(3) * A2
            + coeff(1) * I;
        X = (U + V) * (-U + V).i();
        break;
    }

    return X;
}

ReturnMatrix CESTFwdModel::PadeCoeffs(int m) const
{
    ColumnVector C;
    C.ReSize(m + 1);

    switch (m)
    {
    case 3:
        C << 120 << 60 << 12 << 1;
        break;
    case 5:
        C << 30240 << 15120 << 3360 << 420 << 30 << 1;
        break;
    case 7:
        C << 17297280 << 8648640 << 1995840 << 277200 << 25200 << 1512 << 56 << 1;
        break;
    case 9:
        C << 1.7643225600e10 << 8.821612800e9 << 2.075673600e9 << 3.02702400e8 << 30270240 << 2162160 << 110880 << 3960
          << 90 << 1;
        break;
    case 13:
        C << 6.4764752532480000e16 << 3.2382376266240000e16 << 7.771770303897600e15 << 1.187353796428800e15
          << 1.29060195264000e14 << 1.0559470521600e13 << 6.70442572800e11 << 3.3522128640e10 << 1323241920 << 40840800
          << 960960 << 16380 << 182 << 1;
        break;
    }
    return C;
}

void CESTFwdModel::Mz_spectrum(ColumnVector &Mz, const ColumnVector &wvec, const ColumnVector &w1,
    const ColumnVector &t, const ColumnVector &M0, const Matrix &wi, const Matrix &kij, const Matrix &T12) const
{
    int nfreq = wvec.Nrows();
    int mpool = M0.Nrows();

    Mz.ReSize(nfreq);
    Mz = 0.0;

    // assmeble model matrices
    ColumnVector k1i(mpool);
    ColumnVector k2i(mpool);
    for (int i = 1; i <= mpool; i++)
    {
        k1i(i) = 1 / T12(1, i) + (kij.Row(i)).Sum();
        k2i(i) = 1 / T12(2, i) + (kij.Row(i)).Sum();
    }
    
    // first population of A (with w=0)
    Matrix A(mpool * 3, mpool * 3);
    A = 0.0;
    int st = 0;
    for (int i = 1; i <= mpool; i++)
    {
        Matrix D(3, 3);
        D = 0.0;
        D(1, 1) = -k2i(i);
        D(2, 2) = -k2i(i);
        D(1, 2) = -(wi(i, 1));
        D(2, 1) = wi(i, 1);
        // D(2,3) = -w1; D(3,2) = w1;
        D(3, 3) = -k1i(i);
        st = (i - 1) * 3;
        A.SubMatrix(st + 1, st + 3, st + 1, st + 3) = D;
    }

    int st2 = 0;
    IdentityMatrix I(3);
    for (int i = 1; i <= mpool; i++)
    {
        for (int j = 1; j <= mpool; j++)
        {
            if (i != j)
            {
                st = (i - 1) * 3;
                st2 = (j - 1) * 3;
                // NB 'reversal' of indices is correct here
                A.SubMatrix(st + 1, st + 3, st2 + 1, st2 + 3) = I * kij(j, i);
            }
        }
    }

    ColumnVector M0i(mpool * 3);
    ColumnVector B(mpool * 3);
    B = 0.0;
    M0i = 0.0;
    for (int i = 1; i <= mpool; i++)
    {
        M0i(i * 3) = M0(i);
        B(i * 3) = M0(i) / T12(1, i);
    }

    Matrix M(mpool * 3, nfreq);
    M = 0.0;

    Matrix AinvB(mpool, mpool);
    AinvB = 0.0;
    Matrix ExpmAt;

    for (int k = 1; k <= nfreq; k++)
    {
        if (w1(k) == 0.0)
        {
            // no saturation image - the z water magnetization is just M0
            // (save doing the expm calculation here)
            M(3, k) = M0(1);
        }
        else
        {
            // Calculate new A matrix for this sample
            for (int i = 1; i <= mpool; i++)
            {
                st = (i - 1) * 3;
                // terms that involve the saturation frequency (wvec)
                A(st + 1, st + 2) = -(wi(i, k) - wvec(k));
                A(st + 2, st + 1) = wi(i, k) - wvec(k);

                if (steadystate)
                {
                    // terms that involve the B1 (w1) value - now done below for the pulsed case
                    A(st + 2, st + 3) = -w1(k);
                    A(st + 3, st + 2) = w1(k);
                }
            }

            if (abs(A.Determinant()) < 1e-12)
            {
                cout << A << endl;
            }

            M.Column(k) = M0i; // set result as the intial conditions

            if (steadystate)
            {
                Matrix Ai;
                Ai = A.i();
                M.Column(k) = -Ai * B;
            }
            else
            {
                // work through the segments of the saturation
                int npulse;

                // first: go once through the full pulse and calcualte all the required matrices
                vector<Matrix> AiBseg;
                vector<Matrix> expmAseg;
                for (int s = 1; s <= nseg; s++)
                {
                    // assemble the appropriate A matrix
                    Matrix Atemp(A);
                    for (int i = 1; i <= mpool; i++)
                    {
                        st = (i - 1) * 3;
                        Atemp(st + 2, st + 3) = -w1(k) * pmagvec(s);
                        Atemp(st + 3, st + 2) = w1(k) * pmagvec(s);
                    }

                    // Make the AinvB term
                    Matrix AiBtemp;
                    AiBtemp = Atemp.i() * B;
                    AiBseg.push_back(AiBtemp);

                    // make matrix exponential term
                    Matrix expmAtemp;
                    float tseg;
                    // sort out the duration of the pulse
                    if (ptvec(s) > 1e6)
                    {
                        // a 'continuous' pulse - find duration from dataspec
                        tseg = t(k);
                        npulse = 1;
                    }
                    else
                    {
                        // short duration pulse segment - duration from the pulse specification
                        tseg = ptvec(s);

                        // t now contains the number of times this pulse is repeated
                        npulse = t(k);
                    }

                    expmAtemp = expm(Atemp * tseg);
                    expmAseg.push_back(expmAtemp);
                }

                // now we step through all the pulses
                ColumnVector Mtemp;
                Mtemp = M.Column(k);
                for (int p = 1; p <= npulse; p++)
                {
                    for (int s = 1; s <= nseg; s++)
                    {
                        // Crusher gradients - force transverse magentizations to zero at the very
                        // start of a new pulse
                        if (s == 1)
                        {
                            for (int i = 1; i <= mpool; i++)
                            {
                                Mtemp((i - 1) * 3 + 1) = 0.0; // x
                                Mtemp((i - 1) * 3 + 2) = 0.0; // y
                            }
                        }
                        Mtemp = expmAseg[s - 1] * (Mtemp + AiBseg[s - 1]) - AiBseg[s - 1];
                    }
                }
                M.Column(k) = Mtemp;
            }
        }
    }

    Mz = (M.Row(3)).AsColumn();
}

/**
 * Analytic *steady state* solution to the *one pool* Bloch equations
 * NB t is ignored becasue it is steady state
 */
ReturnMatrix CESTFwdModel::Mz_spectrum_lorentz(const ColumnVector &wvec, const ColumnVector &w1, const ColumnVector &t,
    const ColumnVector &M0, const Matrix &wi, const Matrix &kij, const Matrix &T12) const
{
    int nfreq = wvec.Nrows();

    double R1 = 1. / T12(1, 1);
    double R2 = 1. / T12(2, 1);

    ColumnVector result(wvec);
    double delw;
    for (int k = 1; k <= nfreq; k++)
    {
        delw = wi(1, k) - wvec(k);
        result(k) = (M0(1) * R1 * (R2 * R2 + delw * delw))
            / (R1 * (R2 * R2 + delw * delw) + w1(k) * w1(k) * R2);
    }

    return result;
}

/**
 * More efficicent matrix inversion using the block structure of the problem
 * Implicitly assumes no exchange between pools (aside from water)
 */
void CESTFwdModel::Ainverse(const Matrix A, RowVector &Ai) const
{
    int npool = A.Nrows() / 3;
    int subsz = (npool - 1) * 3;

    Ai.ReSize(npool * 3);

    Matrix DDi(subsz, subsz);
    DDi = 0.0;
    // DDi needs to be the inverse of the lower right blocks refering to the CEST pools
    int st;
    int en;
    for (int i = 1; i < npool; i++)
    {
        st = (i - 1) * 3 + 1;
        en = st + 2;
        DDi.SubMatrix(st, en, st, en) = (A.SubMatrix(st + 3, en + 3, st + 3, en + 3)).i();
    }

    Matrix ABDCi(3, 3);
    ABDCi = (A.SubMatrix(1, 3, 1, 3)
                - A.SubMatrix(1, 3, 4, npool * 3) * DDi * A.SubMatrix(4, npool * 3, 1, 3))
                .i();

    Ai.Columns(1, 3) = ABDCi.Row(3);
    Matrix UR;
    UR = -ABDCi * A.SubMatrix(1, 3, 4, npool * 3) * DDi;
    Ai.Columns(4, npool * 3) = UR.Row(3);
}

/**
 * Contirbution to the spectrum from a  single ('independent') pool
 * A simply parameterised Lorentzian (represeting a single pool solution to the Bloch equation)
 * Following Jones MRM 2012
 */
void CESTFwdModel::Mz_contribution_lorentz_simple(
    ColumnVector &Mzc, const ColumnVector &wvec, const double &I, const ColumnVector &wi, const double &R) const
{
    int nfreq = wvec.Nrows();

    double delw;
    for (int k = 1; k <= nfreq; k++)
    {
        delw = wi(k) - wvec(k);
        Mzc(k) = I * ((R * R) / (4 * delw * delw + R * R));
    }
}