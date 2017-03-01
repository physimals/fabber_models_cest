/*  fwdmodel_cest_devel.cc - Developement CEST APT model

 Michael Chappell, IBME PUMMA & FMRIB Image Analysis Group

 Copyright (C) 2010 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_cest.h"

#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include "newimage/newimageall.h"
#include "miscmaths/miscprob.h"
using namespace NEWIMAGE;
#include "fabber_core/easylog.h"

FactoryRegistration<FwdModelFactory, CESTFwdModel> CESTFwdModel::registration(
		"cest");

static OptionSpec OPTIONS[] =
{
		{ "spec", OPT_MATRIX, "ASCII matrix containing data specification",
				OPT_REQ, "" },
		{ "pools", OPT_MATRIX, "ASCII matrix containing pool specification",
				OPT_REQ, "" },
		{ "expools", OPT_MATRIX,
				"ASCII matrix containing extra pool specification", OPT_NONREQ,
				"" },
		{ "ptrain", OPT_MATRIX,
				"ASCII matrix containing pulsed saturation specification",
				OPT_NONREQ, "" },
		{ "t12prior", OPT_BOOL, "Include uncertainty in T1 and T2 values",
				OPT_NONREQ, "" },
		{ "inferdrift", OPT_BOOL, "", OPT_NONREQ, "" },
		{ "lorentz", OPT_BOOL,
				"Alternative to Matrix exponential solution to Bloch equations",
				OPT_NONREQ, "" },
		{ "steadystate", OPT_BOOL,
				"Alternative to Matrix exponential solution to Bloch equations",
				OPT_NONREQ, "" },
		{ "" }, };

void CESTFwdModel::GetOptions(vector<OptionSpec> &opts) const
{
	for (int i = 0; OPTIONS[i].name != ""; i++)
	{
		opts.push_back(OPTIONS[i]);
	}
}

std::string CESTFwdModel::GetDescription() const
{
	return "Model for Chemical Exchange Saturation transfer";
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

void CESTFwdModel::HardcodedInitialDists(MVNDist& prior,
		MVNDist& posterior) const
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
				prior.means(place) = poolcon(place - 1); //NB poolcon vec doesn't have water in so entry 1 is pool 2 etc
				precisions(place, place) = poolconprec(place - 1);
			}
			else
			{
				//hardcoded default of zero mean and (relatively) uniformative precision
				//prior mean of zero set above
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
			prior.means(place) = poolk(i - 1); //NOTE these shoudl already be log in poolk
			precisions(place, place) = 1;
			place++;
		}
	}

	// frequency offsets (ppm)
	prior.means(place) = 0; //water centre offset
	precisions(place, place) = 10;
	place++;

	if (npool > 1)
	{
		for (int i = 2; i <= npool; i++)
		{
			prior.means(place) = poolppm(i - 1);
			precisions(place, place) = 1e12; //1e6;
			place++;
		}
	}

	// B1 offset (fractional)
	prior.means(place) = 0;
	precisions(place, place) = 1e99; //1e20; //1e6;
	place++;

	if (inferdrift)
	{
		// Drift (ppm/sample)
		prior.means(place) = 0;
		precisions(place, place) = 1e17;
		place++;
	}

	// trough thresh
	//prior.means(place) = 0;
	//precisions(place,place) = 1e-12;
	//place++;

	/*
	 int idx = 10; //should be next entry in priors
	 if (pvcorr) {
	 // M0
	 prior.means(idx) = 0.0;
	 prior.means(idx+1) = 0.0;
	 prior.means(idx+2) = 0.15;

	 precisions(idx,idx) = 1e-12;
	 precisions(idx+1,idx+1) = 100;
	 precisions(idx+2,idx+2) = 2500;

	 // exchnage consts (these are log_e)
	 prior.means(idx+3) = 3;
	 precisions(idx+3,idx+3) = 1;
	 prior.means(idx+4) = 3.7;
	 precisions(idx+4,idx+4) = 1;

	 prior.means(idx+5) = 0.0;
	 precisions(idx+5,idx+5) = 1e-12;

	 prior.means(idx+6) = 1.0;
	 precisions(idx+6,idx+6) = 1e12;
	 prior.means(idx+7) = 0.0;
	 precisions(idx+7,idx+7) = 1e12;
	 prior.means(idx+8) = 0.0;
	 precisions(idx+8,idx+8) = 1e12;
	 idx += 9;
	 }
	 */

	if (t12soft)
	{
		//T1 values
		for (int i = 1; i <= npool; i++)
		{
			prior.means(place) = T12master(1, i);
			precisions(place, place) = 44.4; //all T1s have same prior uncertainty
			place++;
		}

		//T12 values
		for (int i = 1; i <= npool; i++)
		{
			float T12 = T12master(2, i);
			prior.means(place) = T12;
			precisions(place, place) = 1 / std::pow(T12 / 5, 2); //prior has std dev of 1/5 of the value, to try and get the scaling about right
			place++;
		}
	}

	//else { //no T12 inference
	//  for (int i=0; i<npool; i++) {
	//	 precisions(idx+i,idx+i) = 1e99;
	//  }
	//}

	/*
	 // MT pool
	 if (mtpool) {
	 prior.means(place) = 0;
	 precisions(place,place) = 1e12;
	 place++;
	 prior.means(place) = 0;
	 precisions(place,place) = 1e12;
	 place++;
	 prior.means(place) = 0;
	 precisions(place,place) = 1e12;
	 place++;
	 prior.means(place) = 0;
	 precisions(place,place) = 1e12;
	 place++;
	 prior.means(place) = 0;
	 precisions(place,place) = 1e12;
	 place++;
	 }
	 */

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
	//posterior.SetPrecisions(precisions);

	//cout << prior.means << endl;
	//cout << posterior.means << endl;

}

void CESTFwdModel::InitParams(MVNDist& posterior) const
{
	//load the existing precisions as the basis for any update
	SymmetricMatrix precisions;
	precisions = posterior.GetPrecisions();

	//init the M0a value  - to max value in the z-spectrum
	posterior.means(1) = data.Maximum();
	precisions(1, 1) = 10;

	// init the pool concentraitons
	if (npool > 1)
	{
		for (int i = 2; i <= npool; i++)
		{
			posterior.means(i) = poolcon(i - 1); //NB poolcon vec doesn't have water in so entry 1 is pool 2 etc
			precisions(i, i) = 1e6;
		}
	}

	//init the ppmoff value - by finding the freq where the min of z-spectrum is
	int ind;
	float val;
	val = data.Minimum1(ind); // find the minimum in the z-spectrum
	val = wvec(ind) * 1e6 / wlam; //frequency of the minimum in ppm
	if (val > 0.5)
		val = 0.5; //put a limit on the value
	if (val < -0.5)
		val = -0.5;
	int ppmind = 2 * npool;
	posterior.means(ppmind) = val;

	posterior.SetPrecisions(precisions);

}

void CESTFwdModel::Evaluate(const ColumnVector& params,
		ColumnVector& result) const
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

	//model matrices
	ColumnVector M0(npool);
	Matrix kij(npool, npool);
	Matrix T12(2, npool);
	ColumnVector ppmvec(npool);
	//MT pool
	ColumnVector mtparams(5);
	//Extra pools
	ColumnVector exI(nexpool);
	ColumnVector exppmvec(nexpool);
	ColumnVector exR(nexpool);

	// extract values from params
	// M0 comes first
	int place = 1;
	M0(1) = paramcpy(place); // this is the 'master' M0 value of water
	if (M0(1) < 1e-4)
		M0(1) = 1e-4; //M0 of water cannot disapear all together
	place++;

	// values in the parameters are ratios of M0_water
	float M0ratio;
	for (int j = 2; j <= npool; j++)
	{
		M0ratio = paramcpy(place);
		if (M0ratio > 1.0)
		{ //dont expect large ratios
			M0ratio = 1.0;
		}
		M0(j) = M0ratio * M0(1);

		place++;
	}

	// now exchange - we assume that only significant exchnage occurs with water
	kij = 0.0; //float ktemp;
	for (int j = 2; j <= npool; j++)
	{
		kij(j, 1) = exp(params(place));  //non-linear transformation
		if (kij(j, 1) > 1e6)
			kij(j, 1) = 1e6; //exclude really extreme values
		kij(1, j) = kij(j, 1) * M0(j) / M0(1);

		place++;

	}

	// frequency offset next
	ppmvec = params.Rows(place, place + npool);
	place += npool;
	// Frist entry is offset due to field, rest are res freq. of the pools rel. to water

	/* OLD
	 float ppm_off = params(place); // frequncy offset due to field
	 place++;

	 float ppm_apt = params(place);
	 place++;
	 float ppm_mt = params(place);
	 place++;
	 END OLD */

	// now B1 offset
	float B1off = params(place) * 1e6; //scale the B1_off parameter to achieve proper updating of this parameter
									   // (otherwise it is sufficently small that numerical diff is ineffective)
	place++;

	// Drift
	float drift = 0;
	if (inferdrift)
	{
		drift = params(place) * 1e6; //scale this parameters like B1 offset
		place++;
	}

	//float floorval = params(place);

	//place++;

	/*
	 // PV correction section (untested)
	 ColumnVector M0_WM(npool);
	 Matrix kij_WM(npool,npool);
	 Matrix T12_WM(2,npool);
	 Matrix T12_CSF(2,1);
	 Matrix M0_CSF(1,1);
	 float PV_GM;
	 float PV_WM;
	 float PV_CSF;
	 if (pvcorr)
	 {

	 // now parameters for WM
	 M0_WM(1) = paramcpy(place); // this is the M0 value of water
	 if (M0_WM(1)<1e-4) M0_WM(1)=1e-4; //M0 of water cannot disapear all together
	 place++;
	 // this code here if parameters contains actual M0 values
	 //M0.Rows(2,npool) = paramcpy.Rows(place,place+npool-1-1);
	 // place +=npool-1;

	 // values in the parameters are ratios of M0_water
	 float M0WMratio;
	 for (int j=2; j<=npool; j++) {
	 M0WMratio = paramcpy(place);
	 if (M0WMratio > 0.1) { //dont expect large ratios
	 M0WMratio = 0.1;
	 }
	 M0_WM(j) = M0WMratio*M0_WM(1);

	 place++;
	 }

	 // now exchange - we assume that only significant exchnage occurs with water
	 kij_WM=0.0; //float ktemp;
	 for (int j=2; j<=npool; j++) {
	 //ktemp = exp(params(place)); //this is the 'fundamental rate const - non-linear transformation
	 //kij(j,1) = ktemp*M0(1);
	 //kij(1,j) = ktemp*M0(j);
	 kij_WM(j,1) = exp(params(place));  //non-linear transformation
	 kij_WM(1,j) =kij_WM(j,1)*M0_WM(j)/M0_WM(1);

	 place++;
	 }

	 // CSF
	 M0_CSF = paramcpy(place);
	 place++;

	 //partial volume estimates
	 PV_GM = paramcpy(place);
	 place++;
	 PV_WM = paramcpy(place);
	 place++;
	 PV_CSF = paramcpy(place);
	 place++;

	 //T12 fixed
	 T12_WM = T12WMmaster;
	 T12_CSF = T12CSFmaster;
	 }
	 */

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
				T12(2, i) = 1e-12; //0 is no good for a T2 value

			if (T12(2, i) > 1)
				T12(2, i) = 1; // Prevent convergence issues causing T2 to blow up
			place++;
		}
	}
	else
	{
		T12 = T12master;
	}

	// MT pool
	/*
	 if (mtpool) {
	 mtparams = params.SubMatrix(place,place+5,1,1);
	 place += 5;
	 }
	 */

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

	//cout << "Parameters set up" << endl;
	//cout << "M0: " << M0.t() << endl << "wi: " << wi.t() << endl << "kij: " << kij << endl;

	//MODEL - CEST

	//no saturation image first
	//ColumnVector nosat(1);
	/*if (pvcorr) {
	 nosat = PV_GM*M0.Row(1) + PV_WM*M0_WM.Row(1) + PV_CSF*M0_CSF.Row(1);}
	 */
	//nosat = M0.Row(1);
	// deal with frequencies
	//ColumnVector wi(npool);
	int nsamp = wvec.Nrows();
	Matrix wimat(npool, nsamp);

	//float wlocal = wlam*ppmvec(1)/1e6; //local water frequency
	//cout << wlocal << endl;
	//wi(1) = wlocal;
	for (int j = 1; j <= nsamp; j++)
	{
		wimat(1, j) = wlam * (ppmvec(1) + (j - 1) * drift) / 1e6;
	}
	if (npool > 1)
	{
		for (int i = 2; i <= npool; i++)
		{
			//wi(i) = wlam*ppmvec(i)/1e6 + wlocal*(1+ppmvec(i)/1e6);
			for (int j = 1; j <= nsamp; j++)
			{
				wimat(i, j) = wlam * ppmvec(i) / 1e6
						+ wlam * (ppmvec(1) + (j - 1) * drift) / 1e6
								* (1 + ppmvec(i) / 1e6);
			}
		}
	}
	// species b is at ppm*wlam, but also include offset of main field
	//cout << wi << endl;

	// frequencies for the extra pools
	Matrix exwimat(nexpool, nsamp);
	if (nexpool > 0)
	{
		for (int i = 1; i <= nexpool; i++)
		{
			for (int j = 1; j <= nsamp; j++)
			{
				exwimat(i, j) = wlam * ppmvec(i) / 1e6
						+ wlam * (ppmvec(1) + (j - 1) * drift) / 1e6
								* (1 + exppmvec(i) / 1e6);
			}
		}
	}

	//deal with B1
	if (B1off < -0.5)
		B1off = -0.5; // B1 cannot go too small
	if (B1off > 10)
		B1off = 10; //unlikely to get this big (hardlimit in case of convergence problems)
	ColumnVector w1 = w1vec * (1 + B1off); // w1 in radians!

	/*
	 if (pvcorr) {
	 ColumnVector GM_result;
	 GM_result = Mz_spectrum(wvec[n],w1,t[n],M0,wi,kij,T12);
	 ColumnVector WM_result;
	 WM_result = Mz_spectrum(wvec[n],w1,t[n],M0_WM,wi,kij_WM,T12_WM);
	 ColumnVector CSF_result;
	 CSF_result = Mz_spectrum(wvec[n],w1,t[n],M0_CSF,wi.Row(1),kij.SubMatrix(1,1,1,1),T12_CSF);

	 thisresult = PV_GM*GM_result + PV_WM*WM_result + PV_CSF*CSF_result;
	 }
	 */

	if (lorentz)
	{
		result = Mz_spectrum_lorentz(wvec, w1, tsatvec, M0, wimat, kij, T12);
	}
	else
	{
		// Call to Steady State CEST Model
		if (m_SS)
		{
			// Only need to fix w1EX if using SS CEST
			ColumnVector w1EX = m_EXmagMax * (1+B1off);
			Mz_spectrum_SS(result, wvec, w1, tsatvec, M0, wimat, kij, T12, w1EX);
		}
		else
		{
			Mz_spectrum(result, wvec, w1, tsatvec, M0, wimat, kij, T12);
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
			Mz_contribution_lorentz_simple(Mzc, wvec, exI(i),
					(exwimat.SubMatrix(i, i, 1, nsamp)).t(), exR(i));
			Mz_extrapools += Mzc;
		}
	}

	//cout << Mz_extrapools << endl;
	result = result - M0(1) * Mz_extrapools;

	//for (int i=1; i<=result.Nrows(); i++) {
	//  if(result(i) < floorval) result(i) = floorval;
	//}

	//result &= nosat;

	//cout << cest_result.t() << endl;

	return;
}

FwdModel* CESTFwdModel::NewInstance()
{
	return new CESTFwdModel();
}

void CESTFwdModel::Initialize(ArgsType& args)
{
	string scanParams = args.ReadWithDefault("scan-params", "cmdline");

	Matrix dataspec; //matrix containing spec for each datapoint
	// 3 Columns: Freq (ppm), B1 (T), tsat (s)
	// Nrows = number data points

	Matrix poolmat; // matrix containing setup information for the pools
	//  columns: res. freq (rel. to water) (ppm), rate (species-> water), T1, T2.
	// 1st row is water, col 1 will be interpreted as the actaul centre freq of water if value is >0, col 2 will be ignored.
	// Further rows for pools to be modelled.

	Matrix expoolmat; //Matrix for the 'extra' pools

	string expoolmatfile;
	string pulsematfile;
	string Readmatfile;
	string EXmatfile;

	t12soft = false; //pvcorr=false;

	if (scanParams == "cmdline")
	{
		// read data specification from file
		dataspec = read_ascii_matrix(args.Read("spec"));

		//read pool specification from file
		poolmat = read_ascii_matrix(args.Read("pools"));

		//read extra pool specification from file
		expoolmatfile = args.ReadWithDefault("expools", "none");

		//read pulsed saturation specification
		pulsematfile = args.ReadWithDefault("ptrain", "none");

		//read Excitation Pulse specification
		EXmatfile = args.GetStringDefault("EXtrain", "none");

		//read readout specification
		Readmatfile = args.ReadWithDefault("readspec", "none");

		//basic  = args.ReadBool("basic");
		//      pvcorr = args.ReadBool("pvcorr");
		t12soft = args.ReadBool("t12prior");
		inferdrift = args.ReadBool("inferdrift");

		// alternatives to Matrix exponential solution to Bloch equations
		lorentz = args.ReadBool("lorentz"); //NB only compatible with single pool
		steadystate = args.ReadBool("steadystate");

	}

	else
		throw invalid_argument(
				"Only --scan-params=cmdline is accepted at the moment");

	//Deal with the specification of the pools
	npool = poolmat.Nrows();
	if (poolmat.Ncols() < 4 || poolmat.Ncols() > 6)
	{
		throw invalid_argument(
				"Incorrect number of columns in pool specification file");
	}
	// water centre
	float wdefault = 42.58e6 * 3 * 2 * M_PI; //the default centre freq. (3T)
	if (poolmat(1, 1) > 0)
		wlam = -poolmat(1, 1) * 2 * M_PI; // Should be negative to reflect conversion between ppm and Hz
	else
		wlam = -wdefault;
	// ppm offsets
	poolppm = poolmat.SubMatrix(2, npool, 1, 1);
	// exchange rate
	poolk = log(poolmat.SubMatrix(2, npool, 2, 2)); //NOTE the log_e transformation
	// T1 and T2 values
	T12master = (poolmat.SubMatrix(1, npool, 3, 4)).t();

	if (poolmat.Ncols() > 4)
	{
		// pool ratios (concentrations) - these are used for initialisaiton, but may also be used as prior means if precisions are provided
		poolcon = poolmat.SubMatrix(2, npool, 5, 5); //ignore water pool
	}
	else
	{
		poolcon.ReSize(npool);
		poolcon = 0;
	}

	setconcprior = false;
	if (poolmat.Ncols() > 5)
	{
		// pool ratio precisions - this overrides the inbuilt priors for the concentrations using these precision and the means in poolcon
		poolconprec = poolmat.SubMatrix(2, npool, 6, 6); //ingore water pool
		setconcprior = true;
	}

	// check that the method chosen is possible
	if ((npool > 1) & lorentz)
		throw invalid_argument(
				"Lorentzian (analytic) solution only compatible with single pool");

	if (expoolmatfile != "none")
	{
		//process the Extra pool specification matrix
		expoolmat = read_ascii_matrix(expoolmatfile);
		nexpool = expoolmat.Nrows();
		//cout << "nexpool: " << nexpool << endl;
		if (expoolmat.Ncols() != 2)
			throw invalid_argument(
					"Incorrect number of columns in extra pool spefication file");
		// ppm ofsets
		expoolppm = expoolmat.SubMatrix(1, nexpool, 1, 1);
		// (prior) R values
		expoolR = expoolmat.SubMatrix(1, nexpool, 2, 2);
	}

	/* OLD
	 //initialization
	 npool = 3;

	 // some fixed things
	 T12master.ReSize(2,npool);
	 T12master << 1.3 << 0.77  << 1
	 << 0.05 << 0.01 << 1e-5;  //T2
	 END OLD */

	/*T12WMmaster.ReSize(2,npool);
	 T12WMmaster << 1 << 0.77  << 1
	 << 0.02 << 0.01 << 1e-5;  //T2
	 T12CSFmaster.ReSize(2,1);
	 T12CSFmaster << 3.5 << .75;
	 */

	//ppm_apt = -3.5;
	//ppm_mt = -2.41;
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

	LOG << "Sampling frequencies (ppm):" << endl << (dataspec.Column(1)).t()
			<< endl << "                     (rad/s):" << endl << wvec.t()
			<< endl;
	LOG << "B1 values (uT):" << endl << (dataspec.Column(2)).t() * 1e6 << endl
			<< "          (rad/s):" << endl << w1vec.t() << endl;

	//Pulsed saturation modelling
	// For pulsed saturation the B1 values are taken as peak values
	if (pulsematfile == "none")
	{
		cout
				<< "WARNGING! - you should supply a pulsemat file, in future the calling of this model without one (to get continuous saturation) will be depreceated"
				<< endl;
		nseg = 1;
		pmagvec.ReSize(1);
		ptvec.ReSize(1);
		pmagvec = 1.0;
		ptvec = 1e12; //a very long value for continuous saturation

		LOG << "Saturation times (s):" << endl << tsatvec.t() << endl;
	}
	else
	{
		Matrix pulsemat;
		pulsemat = read_ascii_matrix(pulsematfile);
		pmagvec = pulsemat.Column(1); //vector of (relative) magnitude values for each segment
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

	//Steady State Modeling
	// For Modeling the steady state Signal according to Listerud, Magn Reson Med 1997; 37: 693–705.
	if (Readmatfile == "none")
	{
		m_SS = false;
		LOG << "Running Original Fabber CEST Method" << endl;
	}
	else
	{
		// Load Readout Parameters
		Matrix Readmat;
		Readmat = read_ascii_matrix(Readmatfile);
		m_EXmagMax = Readmat.Column(1) * 42.58e6 * 2 * M_PI;
		m_TR = Readmat.Column(2);
		m_SS = true;

		LOG << "\nRunning Steady State CEST Model Based upon Listerud, MRM: 37 (1997)\n" << endl;

		Matrix EXmat;
		EXmat = read_ascii_matrix(EXmatfile);
		m_EXmagvec = EXmat.Column(1); //vector of normalized (to 1) magnitude values for each segment
		// vector of time durations for each segment
		// We are loading a vector of times for the end of each segment + zero
		ColumnVector ettemp = EXmat.Column(2);
		m_EXtvec = ettemp.Rows(2,ettemp.Nrows()) - ettemp.Rows(1,ettemp.Nrows()-1);

		m_nEXseg = EXmat.Nrows();

		LOG << "Excitation Pulse length (s): \n" << endl << ettemp.Row(m_nEXseg) << endl;
		LOG << "TR (s): \n" << m_TR.t() << endl;

		LOG << "Pulse shape:" << endl << "Number of segments: " << m_nEXseg << endl;
		LOG << "Excitation B1 values (uT):" << endl << Readmat.Column(1).t()*1e6 << endl
				<< "          (rad/s):" << endl << m_EXmagMax.t() << endl;
		LOG << " Magnitudes (Normalized to 1): \n" << m_EXmagvec.t() << endl;
		LOG << " Durations (s): \n" << m_EXtvec.t() << endl;

	}


	/* OLD
	 //deal with water centre frequency *in radians!*
	 float wdefault = 42.58e6*3*2*M_PI; // this is the default value (3T)

	 wlam.resize(t.size());
	 for (unsigned int n=0; n<t.size(); n++) {
	 //iterate over all the datasets

	 //water centre frequency
	 if (wcentre[n]>0)     wlam[n] = wcentre[n]*1e6*2*M_PI; //NB input is in MHz
	 else                  wlam[n] = wdefault;

	 // CEST frequency vector
	 //int nfreq=32;
	 ColumnVector ppmvec(nfreq[n]);
	 float ppminc;
	 if (wrange[n] > 0.0)
	 { //deafult case is that wrange is specified
	 ppmvec(1) = -wrange[n];
	 ppminc = (2*wrange[n])/(nfreq[n]-1);
	 }
	 else
	 { //if wrange is zero then we use wstart and wstep parameters
	 ppmvec(1) = wstart[n];
	 ppminc = wstep[n];
	 }
	 //cout << ppminc << endl;
	 for (int i=2; i<=nfreq[n]; i++) {
	 ppmvec(i) = ppmvec(i-1) + ppminc;
	 }

	 //wvec.ReSize(nfreq);
	 ColumnVector freqvec(nfreq[n]);
	 freqvec = ppmvec*wlam[n]/1e6;

	 if (revfrq) freqvec = -freqvec;
	 //cout << freqvec.t() << endl;
	 wvec.push_back( freqvec );
	 //cout << wvec[n] << endl;
	 //cout << wlam[n] << endl;
	 }

	 */

}

vector<string> CESTFwdModel::GetUsage() const
{
	vector<string> usage;
	usage.push_back("\nUsage info for --model=cest:\n");
	usage.push_back("Undefined\n");

	return usage;
}

void CESTFwdModel::DumpParameters(const ColumnVector& vec,
		const string& indent) const
{
	//cout << vec.t() << endl;
}

void CESTFwdModel::NameParams(vector<string>& names) const
{
	names.clear();

	// name the parameters for the pools using letters
	string lettervec[] =
	{ "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o",
			"p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z" };

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
	/*  if (pvcorr) {
	 names.push_back("WM_M0a");
	 names.push_back("WM_M0b_r");
	 names.push_back("WM_M0c_r");
	 names.push_back("WM_kba");
	 names.push_back("WM_kca");
	 names.push_back("CSF_M0a");
	 }*/
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

ReturnMatrix CESTFwdModel::expm_eig(Matrix inmatrix) const
{
	// Do matrix exponential using eigen decomposition of the matrix
	SymmetricMatrix A;
	A << inmatrix; // a bit poor - the matrix coming in should be symmetric, but I haven't implemented this elsewhere in the code (yet!)

	//eigen decomposition
	DiagonalMatrix D;
	Matrix V;
	Jacobi(A, D, V);

	//exponent of D
	for (int i = 1; i <= D.Nrows(); i++)
	{
		//if (D(i)>100) D(i)=100;
		D(i) = exp(D(i));
	}

	// now matrix exponential
	Matrix EXP(A);
	EXP = V * D * V.t();
	return EXP;

}

ReturnMatrix CESTFwdModel::expm_pade(Matrix inmatrix) const
{
	// Do matrix exponential
	// Algorithm from Higham, SIAM J. Matrix Analysis App. 24(4) 2005, 1179-1193

	Matrix A = inmatrix;
	Matrix X(A.Nrows(), A.Ncols());

	//Coeff of degree 13 Pade approximant
	//ColumnVector b(13);
	//  b= PadeCoeffs(13);

	int metflag = 0;
	int mvals[5] =
	{ 3, 5, 7, 9, 13 };
	//ColumnVector thetam(5);
	float thetam[5] =
	{ 1.495585217958292e-002, 2.539398330063230e-001, 9.504178996162932e-001,
			2.097847961257068e+000, 5.371920351148152e+000 };
	int i = 0;

	while (metflag < 1)
	{
		if (i < 4)
		{
			//cout << "Try m = " << mvals[i] << endl;
			if (A.Norm1() <= thetam[i])
			{
				X = PadeApproximant(A, mvals[i]);
				metflag = 1;
			}
		}
		else
		{
			// Using Pade Approximant degree 13 and scaling and squaring
			//cout << "Doing m = 13" << endl;
			int s = ceil(log2(A.Norm1() / thetam[4]));
			//cout << s << endl;
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
	//cout << "PadeApproximant" << endl;
	//cout << inmatrix << endl;
	assert(inmatrix.Nrows() == inmatrix.Ncols());
	int n = inmatrix.Nrows();
	IdentityMatrix I(n);
	ColumnVector coeff = PadeCoeffs(m);
	Matrix X(n, n);
	Matrix U(n, n);
	Matrix V(n, n);
	vector<Matrix> Apowers;

	//cout << coeff << endl;

	switch (m)
	{
	case 3:
	case 5:
	case 7:
	case 9:
		//cout << "case " << m << endl;
		Apowers.push_back(I);
		Apowers.push_back(inmatrix * inmatrix);
		for (int j = 3; j <= ceil((m + 1 / 2)); j++)
		{
			Apowers.push_back(Apowers[j - 2] * Apowers[1]);
		}
		U = 0.0;
		V = 0.0;
		//cout << "Apowers" <<endl;
		for (int j = m + 1; j >= 2; j -= 2)
		{
			U += coeff(j) * Apowers[j / 2 - 1];
		}
		U = inmatrix * U;
		for (int j = m; j >= 1; j -= 2)
		{
			V += coeff(j) * Apowers[(j + 1) / 2 - 1];
		}
		//cout << U << endl << V << endl;
		X = (U + V) * (-U + V).i();
		//cout << X<< endl;
		break;
	case 13:
		//cout << "case 13" << endl;
		Matrix A2(n, n);
		Matrix A4(n, n);
		Matrix A6(n, n);
		A2 = inmatrix * inmatrix;
		A4 = A2 * A2;
		A6 = A2 * A4;
		U = inmatrix
				* (A6 * (coeff(14) * A6 + coeff(12) * A4 + coeff(10) * A2)
						+ coeff(8) * A6 + coeff(6) * A4 + coeff(4) * A2
						+ coeff(2) * I);
		//cout << U;
		V = A6 * (coeff(13) * A6 + coeff(11) * A4 + coeff(9) * A2)
				+ coeff(7) * A6 + coeff(5) * A4 + coeff(3) * A2 + coeff(1) * I;
		//cout << V;
		X = (U + V) * (-U + V).i();
		//cout << X <<endl;
		break;
	}

	return X;
}

ReturnMatrix CESTFwdModel::PadeCoeffs(int m) const
{

	ColumnVector C;
	C.ReSize(m + 1);

	//cout << "PadeCoeffs" << endl;

	switch (m)
	{
	case 3:
		C << 120 << 60 << 12 << 1;
		break;
	case 5:
		C << 30240 << 15120 << 3360 << 420 << 30 << 1;
		break;
	case 7:
		C << 17297280 << 8648640 << 1995840 << 277200 << 25200 << 1512 << 56
				<< 1;
		break;
	case 9:
		C << 1.7643225600e10 << 8.821612800e9 << 2.075673600e9 << 3.02702400e8
				<< 30270240 << 2162160 << 110880 << 3960 << 90 << 1;
		break;
	case 13:
		C << 6.4764752532480000e16 << 3.2382376266240000e16
				<< 7.771770303897600e15 << 1.187353796428800e15
				<< 1.29060195264000e14 << 1.0559470521600e13 << 6.70442572800e11
				<< 3.3522128640e10 << 1323241920 << 40840800 << 960960 << 16380
				<< 182 << 1;
		break;
	}
	return C;
}

void CESTFwdModel::Mz_spectrum(ColumnVector& Mz, const ColumnVector& wvec,
		const ColumnVector& w1, const ColumnVector& t, const ColumnVector& M0,
		const Matrix& wi, const Matrix& kij, const Matrix& T12) const
{

	int nfreq = wvec.Nrows(); // total number of samples collected
	int mpool = M0.Nrows();

	Mz.ReSize(nfreq);
	Mz = 0.0;

	//assmeble model matrices
	ColumnVector k1i(mpool);
	ColumnVector k2i(mpool);
	for (int i = 1; i <= mpool; i++)
	{
		k1i(i) = 1 / T12(1, i) + (kij.Row(i)).Sum();
		k2i(i) = 1 / T12(2, i) + (kij.Row(i)).Sum();
	}

	//cout << "k matrices generated" << endl;
	//cout << k1i << endl << k2i << endl;

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
		//D(2,3) = -w1; D(3,2) = w1;
		D(3, 3) = -k1i(i);
		//cout << D << endl;
		st = (i - 1) * 3;
		A.SubMatrix(st + 1, st + 3, st + 1, st + 3) = D;
	}
	//cout << A << endl;

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
				A.SubMatrix(st + 1, st + 3, st2 + 1, st2 + 3) = I * kij(j, i); //NB 'reversal' of indices is correct here
			}
		}
	}

	//cout << "Inital A matrix populated" << endl;
	//cout << A << endl;

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
			//Mz(k) = M0(1);
		}
		else
		{
			//Calculate new A matrix for this sample
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

			M.Column(k) = M0i; //set result as the intial conditions

			if (steadystate)
			{
				Matrix Ai;
				Ai = A.i();
				//Mz(k) = -(Ai.Row(3)*B).AsScalar(); //Only want the z-component of the water pool
				M.Column(k) = -Ai * B;
			}
			else
			{
				//work through the segments of the saturation
				int npulse;

				// first: go once through the full pulse and calcualte all the required matrices
				vector<Matrix> AiBseg;
				vector<Matrix> expmAseg;
				for (int s = 1; s <= nseg; s++)
				{
					//assemble the appropriate A matrix
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
						//a 'continuous' pulse - find duration from dataspec
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

						//Crusher gradients - force transverse magentizations to zero at the very start of a new pulse
						if (s == 1)
						{
							for (int i = 1; i <= mpool; i++)
							{
								Mtemp((i - 1) * 3 + 1) = 0.0; //x
								Mtemp((i - 1) * 3 + 2) = 0.0; //y
							}
						}
						Mtemp = expmAseg[s - 1] * (Mtemp + AiBseg[s - 1])
								- AiBseg[s - 1];
					}
				}
				M.Column(k) = Mtemp;

				/*
				 //segment 1 - from the initial conditions
				 int s=1;
				 float tcum = 0;

				 for (int i=1; i<=mpool; i++) {
				 st = (i-1)*3;
				 A(st+2,st+3) = -w1(k)*pmagvec(s);
				 A(st+3,st+2) = w1(k)*pmagvec(s);
				 }
				 AinvB = A.i()*B;
				 //ExpmAt = expm(A*t(k,1));
				 //Mz(k) = (ExpmAt.Row(3) * (M0i + AinvB) - AinvB).AsScalar();
				 M.Column(k) = expm(A*ptvec(s)) * (M0i + AinvB) - AinvB;

				 tcum += ptvec(s);
				 cout << "Finished segment 1, saturation time completed is:" <<tcum << endl;
				 s++;

				 while (tcum<t(k)) {
				 // terms in A the involve the B1 (w1) value
				 for (int i=1; i<=mpool; i++) {
				 st = (i-1)*3;
				 A(st+2,st+3) = -w1(k)*pmagvec(s);
				 A(st+3,st+2) = w1(k)*pmagvec(s);
				 }

				 /
				 if (pmagvec(s)<1e-12) { //crush the transverse components during the zero saturation phases
				 for (int i=1; i<=mpool; i++) {
				 st = (i-1)*3;
				 M(i+1,k)=0.0;
				 M(i+2,k)=0.0;
				 }
				 }
				 /

				 AinvB = A.i()*B;
				 //ExpmAt = expm(A*t(k,s));
				 //Mz(k) = (ExpmAt.Row(3) * (M0i + AinvB) - AinvB).AsScalar();
				 M.Column(k) = expm(A*ptvec(s)) * (M.Column(k) + AinvB) - AinvB;

				 tcum += ptvec(s);
				 //cout << "Finished segment " << s << ", saturation time completed is:" <<tcum << endl;
				 s++;
				 if (s>nseg) s=1; //reset back tot he first segment of the pulsetrain
				 }
				 */
			}
		}
	}

	Mz = (M.Row(3)).AsColumn();

	//ColumnVector result;
	//cout << M.Row(3);
	//result = (M.Row(3)).AsColumn();
	//result = Mz;
	//return result;
}

ReturnMatrix CESTFwdModel::Mz_spectrum_lorentz(const ColumnVector& wvec,
		const ColumnVector& w1, const ColumnVector& t, const ColumnVector& M0,
		const Matrix& wi, const Matrix& kij, const Matrix& T12) const
{

	//Analytic *steady state* solution to the *one pool* Bloch equations
	// NB t is ignored becasue it is ss

	int nfreq = wvec.Nrows(); // total number of samples collected

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

void CESTFwdModel::Ainverse(const Matrix A, RowVector& Ai) const
{
	// More efficicent matrix inversion using the block structure of the problem
	// Implicitly assumes no exchange between pools (aside from water)

	int npool = A.Nrows() / 3;
	int subsz = (npool - 1) * 3;

	//Ai.ReSize(npool*3,npool*3);;
	Ai.ReSize(npool * 3);

	// Matrix AA(3,3);
//   AA = A.SubMatrix(1,3,1,3);
//   Matrix BB(3,subsz);
//   BB = A.SubMatrix(1,3,4,npool*3);
//   Matrix CC(subsz,3);
//   CC = A.SubMatrix(4,npool*3,1,3);

	Matrix DDi(subsz, subsz);
	DDi = 0.0;
	//DDi needs to be the inverse of the lower right blocks refering to the CEST pools
	int st;
	int en;
	for (int i = 1; i < npool; i++)
	{
		st = (i - 1) * 3 + 1;
		en = st + 2;
		DDi.SubMatrix(st, en, st, en) = (A.SubMatrix(st + 3, en + 3, st + 3,
				en + 3)).i();
	}

	Matrix ABDCi(3, 3);
	ABDCi = (A.SubMatrix(1, 3, 1, 3)
			- A.SubMatrix(1, 3, 4, npool * 3) * DDi
					* A.SubMatrix(4, npool * 3, 1, 3)).i();

	Ai.Columns(1, 3) = ABDCi.Row(3);
	Matrix UR;
	UR = -ABDCi * A.SubMatrix(1, 3, 4, npool * 3) * DDi;
	Ai.Columns(4, npool * 3) = UR.Row(3);

	//Ai.SubMatrix(1,3,1,3) = ABDCi;
	//Ai.SubMatrix(1,3,4,npool*3) = -ABDCi*A.SubMatrix(1,3,4,npool*3)*DDi;
	//Ai.SubMatrix(4,npool*3,1,3) = -DDi*A.SubMatrix(4,npool*3,1,3)*ABDCi;
	//Ai.SubMatrix(4,npool*3,4,npool*3) = DDi + DDi*A.SubMatrix(4,npool*3,1,3)*ABDCi*A.SubMatrix(1,3,4,npool*3)*DDi;

	//return Ai;

}

void CESTFwdModel::Mz_contribution_lorentz_simple(ColumnVector& Mzc,
		const ColumnVector& wvec, const double& I, const ColumnVector& wi,
		const double& R) const
{

	// Contirbution to the spectrum from a  single ('independent') pool
	// A simply parameterised Lorentzian (represeting a single pool solution to the Bloch equation)
	// Following Jones MRM 2012

	int nfreq = wvec.Nrows(); // total number of samples collected

	double delw;
	for (int k = 1; k <= nfreq; k++)
	{
		delw = wi(k) - wvec(k);
		Mzc(k) = I * ((R * R) / (4 * delw * delw + R * R));
	}

}


// Models the Bloch-McConnell equations based on Listerud, Magn Reson Med 1997; 37: 693–705.
void CESTFwdModel::Mz_spectrum_SS(
		ColumnVector& Mz,			// Vector: Magnetization
		const ColumnVector& wvec, 	// Vector: Saturation Pulse Offset (radians/s = ppm * 42.58*B0*2*pi)
		const ColumnVector& w1,   	// Vector: B1-corrected Saturation Pulse (radians = uT*42.58*2*pi)
		const ColumnVector& t,    	// Vector: Number Pulses
		const ColumnVector& M0,		// Vector: Pool Sizes
		const Matrix& wi,			// Matrix: Pool offsets (radians/s = ppm * 42.58*B0*2*pi)
		const Matrix& kij,			// Matrix: exchange rates for each pool (see below for in depth description)
		const Matrix& T12,			// Matrix: T1's (Row 1) and T2's (Row 2)
		const ColumnVector& w1EX     // Vector: B1-corrected Excitation Flip Angle
) const
{

	/*********************************************************************
	 * Description of Exchange Rate Matrix
	 *********************************************************************
	 3 Pool for example (Water - f, MT - m, Solute - s)
	 	  [  0		kfm		kfs ]
	 kij =[  kmf     0		0	]
	  	  [	 ksf	 0		0	]


	 *********************************************************************/

	// total number of samples collected
	int nfreq = wvec.Nrows();

	// Number of Pools to be solved
	int mpool = M0.Nrows();


	Mz.ReSize(nfreq);
	Mz = 0.0;

	/**********************************************************************
	 *					Assemble model matrices
	 **********************************************************************/

	// Find Diagonals of Relaxation Matrix
	ColumnVector k1i(mpool);
	ColumnVector k2i(mpool);
	for (int i = 1; i <= mpool; i++)
	{
		k1i(i) = 1 / T12(1, i) + (kij.Row(i)).Sum();
		k2i(i) = 1 / T12(2, i) + (kij.Row(i)).Sum();
	}


	// Populate Diagonals of Relaxation Matrix
	Matrix A(mpool * 3, mpool * 3);
	A = 0.0;
	int st = 0;
	for (int i = 1; i <= mpool; i++)
	{
		Matrix D(3, 3);
		D = 0.0;
		D(1, 1) = -k2i(i);
		D(2, 2) = -k2i(i);
		D(3, 3) = -k1i(i);
		st = (i - 1) * 3;
		A.SubMatrix(st + 1, st + 3, st + 1, st + 3) = D;
	}

	// Populate Exchange Parameters of Relaxation Matrix
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
				A.SubMatrix(st + 1, st + 3, st2 + 1, st2 + 3) = I * kij(j, i); //NB 'reversal' of indices is correct here
			}
		}
	}

	// Find Readout Matrix

	double Tr = m_TR(1);
	Tr-= m_EXtvec.Sum();
	Tr -= ptvec.Rows(1,ptvec.Nrows()-1).Sum()*t(1);
	Tr -= ptvec(ptvec.Nrows())*(t(1)-1);

	Matrix Er = expm(A*Tr);

	// Find CEST Saturation Train Pause (for Duty Cycle < 1.0) Matrix
	float Tdc = ptvec(ptvec.Nrows());
	Matrix Edc = expm(A*Tdc);

	// Create Spoiling Matrix Using Square matrix with Transverse elements = 0

	Matrix Spoil (mpool * 3, mpool * 3);
	Spoil = 0.0;
	for (int i = 1; i<= mpool; ++i)
	{
		st = i*3;
		Spoil(st,st) = 1.0;
	}

	// Create End of Saturation Pulse Train Spoiling Delay Matrix with arbitrary 3 ms time
	/* 	NB: This is also the point in time where the magnetization vector, Mz, is solved for.
		Right after this rest, the magnetization is excited into the XY plane.
	*/

	Matrix Es (mpool*3,mpool*3);
	Es = expm(A*3e-3);

	// Create Excitation Matrix

	double ds0 = m_EXtvec.Sum()/m_EXtvec.Nrows();
	Matrix Eps (mpool*3, 1);
	Matrix Ept (mpool * 3, mpool * 3);
	IdentityMatrix Eye(mpool*3);

	for (int ii (1); ii <= m_EXmagvec.Nrows(); ++ii)
	{
		Matrix W0 (mpool*3, mpool *3);
		W0 = 0.0;
		for (int kk (1); kk <= mpool; ++kk)
		{
			Matrix W0s (3, 3);
			W0s = 0.0;
			W0s(1,2) = -wi(kk,1);
			W0s(2,1) = wi(kk,1);
			W0s(2,3) = -m_EXmagvec(ii)*w1EX(1);
			W0s(3,2) = m_EXmagvec(ii)*w1EX(1);
			st = (kk - 1) * 3;
			W0.SubMatrix(st + 1, st + 3, st + 1, st + 3) = W0s;
		}

		Matrix Ep (mpool * 3, mpool * 3);

		Ep = expm((A+W0)*ds0);

		if (ii == 1)
		{
			Ept = Ep;
			Matrix RlW0 (mpool*3,mpool*3);
			RlW0 = A+W0;
			Eps = (Eye-Ep)*RlW0.i();
		}
		else
		{
			Ept = Ep*Ept;
			Matrix RlW0 (mpool*3,mpool*3);
			RlW0 = A+W0;
			Eps = Ep*Eps+(Eye-Ep)*RlW0.i();
		}
	}

	// Create M0i Vector

	ColumnVector M0i(mpool * 3);
	M0i = 0.0;

	for (int i = 1; i <= mpool; i++)
	{
		M0i(i * 3) = M0(i);
	}

	Matrix M(mpool * 3, nfreq);
	M = 0.0;


	/**********************************************************************
	 *					Solve for Mz
	 **********************************************************************/

	for (int k = 1; k <= nfreq; k++)
	{
		if (w1(k) == 0.0)
		{
			// no saturation image - the z water magnetization is just a normal FLASH sequence
			Matrix Er0 (mpool*3, mpool *3);
		    Er0 = expm(A*(m_TR(1)-m_EXtvec.Sum()));
		    ColumnVector Mz0 (mpool * 3);
		    Matrix Er0Temp (mpool*3, mpool*3);
		    Er0Temp = (Eye-Spoil*Er0*Ept);
		    Mz0 = Er0Temp.i()*((Spoil*Er0*Eps*A+Eye-Er0)*M0i);
		    M(3, k) = Mz0(3);
		}
		else
		{
			Matrix Ems (mpool * 3, mpool * 3);
			Matrix Emt (mpool * 3, mpool * 3);

			for (int jj = 1; jj <= nseg-1; ++jj)
			{
				Matrix W (mpool * 3, mpool * 3);
				W = 0.0;
				for (int nn = 1; nn <= mpool; ++nn)
				{
					Matrix Ws (3, 3);
					Ws = 0.0;
					Ws(1,2) = -(wi(nn,k) - wvec(k));
					Ws(2,1) = -Ws(1,2);
					Ws(2,3) = -w1(k)*pmagvec(jj);
					Ws(3,2) = -Ws(2,3);
					st = (nn - 1) * 3;
					W.SubMatrix(st + 1, st + 3, st + 1, st + 3) = Ws;
				}

				Matrix Em (mpool * 3, mpool * 3);

				Em = expm((A+W)*ptvec(jj));

				if (jj == 1)
				{
					Emt = Em;
					Matrix RlW0 (mpool*3,mpool*3);
					RlW0 = A+W;
					Ems = (Eye-Em)*RlW0.i();
				}
				else
				{
					Emt = Em*Emt;
					Matrix RlW0 (mpool*3,mpool*3);
					RlW0 = A+W;
					Ems = Em*Ems+(Eye-Em)*RlW0.i();
				}

			}

			// Build the Pulse Train
			Matrix Emdc (Emt);
			Emdc = mpower(Emt*Edc,t(k)-1);

			Matrix Emm (Emt);
			Emm = Eye;
			Matrix Emb (Emt);
			Emb = Eye;

			for (int jj=1; jj < t(k); ++jj)
			{
				if (jj == t(k)-1)
					Emb = Emm;
				Emm += mpower(Emt*Edc,jj);
			}


			Matrix Mztemp (Emt);
			Mztemp = Eye - Es*Emdc*Emt*Spoil*Er*Ept*Spoil;

			M.Column(k) = Mztemp.i() * (Es*Emdc*Emt*Spoil*Er*Eps*A+Es*Emdc*Emt*Spoil*(Eye-Er)+Es*Emb*Emt*(Eye-Edc)
					+Eye-Es+Es*Emm*Ems*A)*M0i;

		}
	}

	ColumnVector Mtemp = (M.Row(3)).AsColumn();
	Mz = (Mtemp/Mtemp.Maximum())*M0(1);

}

// Function that will raise a matrix to a power Power
inline ReturnMatrix CESTFwdModel::mpower(const Matrix& Mat_Base, int Power) const
{
	Matrix MExp;
	MExp = Mat_Base;
	for (int ii{1}; ii< Power; ++ii)
		MExp *= Mat_Base;

	return MExp;
}
