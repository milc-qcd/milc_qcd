// rationals.h
// rational function approximations for RHMC algorithm
//
// approximation form:  a0 + a1/(x+b1) + ... aN/(x+bN)
//
// fermion action is phi^dagger ( (M^daggerM)^(-NF/4) ) phi
//
// we need x^(-NF/4) for molecular dynamics integration
// we need x^(+NF/8) for heat bath computation of pseudofermions
// we need x^(-NF/8) for fermion action computations
// all of these for both quark masses
//
// The numbers of terms we need and all the coefficients will depend
// on quark masses and weakly on beta as u_0 changes.

#define NFLAVORS1 2
#define NFLAVORS2 1
// These had better match the parameters in the input file!!

#define CG1_ERRORSTRETCH 1.0	// factor to increase CG residual for factor 1 in molecular dynamics
#define MOLDYN_ORDER_1	4	// order of approximation for updating quark type 1 (light, usually)
#define MOLDYN_ORDER_2	6	// order of approximation for updating quark type 2 (strange, usually)
#define GRSOURCE_ORDER_1 6	// order for heat bath pseudofermion type 1
#define GRSOURCE_ORDER_2 8	// type 2
#define ACTION_ORDER_1	6	// order for evaluating fermion action type 1
#define ACTION_ORDER_2	8	// type 2
#define MAX_RAT_ORDER 8	// largest of above orders

#ifdef CONTROL
// Values for m=0.01/.05, nf=2/1  emin=1e=-15 emax=90 
// (x + 4*m1^2)^(-1/2) (x + 4*m2^2 )^(1/2) where x = M^dagger M
// (x + 4*m1^2)^(1/4) (x + 4*m2^2 )^(-1/4) where x = M^dagger M
// (x + 4*m1^2)^(-1/4) (x + 4*m2^2 )^(1/4) where x = M^dagger M
 Real A_MD_1[MOLDYN_ORDER_1+1] =	{ 
1.0000013909229337e+00,
8.5971364812362813e-04,
1.1750830149597894e-03,
1.6111370468384256e-03,
1.1534635785172337e-03
 };
 Real B_MD_1[MOLDYN_ORDER_1+1] =	{
99.9, //DUMMY
4.4462635312990379e-04,
9.1988952691570737e-04,
2.6122364595470749e-03,
6.6966664426365927e-03
 };
 Real A_GR_1[GRSOURCE_ORDER_1+1] =	{
9.9999999870775846e-01,
-1.0465703802991746e-05,
-4.7088487089484707e-05,
-1.5044863737285111e-04,
-4.0578764266913043e-04,
-8.3702518676026158e-04,
-9.4918316418484588e-04
 };
 Real B_GR_1[GRSOURCE_ORDER_1+1] =	{
99.9,  //DUMMY
4.6808443256449259e-04,
7.6466958480072352e-04,
1.5165615643089338e-03,
3.1609298697618791e-03,
6.0591633587739170e-03,
9.2041295998850949e-03
 };
 Real A_FA_1[ACTION_ORDER_1+1] =	{
1.0000000012922416e+00,
1.0020633119693040e-04,
2.0390410350389638e-04,
3.6323500936637476e-04,
5.8505946932939520e-04,
7.2031751800945190e-04,
4.2727639667627308e-04
 };
 Real B_FA_1[ACTION_ORDER_1+1] =	{
99.9, //DUMMY
4.3458385994507538e-04,
6.6013015248349359e-04,
1.2653664261906257e-03,
2.6373623879382994e-03,
5.2307615518001850e-03,
8.5453352067564187e-03
 };

// 3 flavors at the strange mass
// Values for m=0.01/.05, nf=2/1  emin=1e=-15 emax=90 
// (x + 4*m2^2)^(-3/4) where x = M^dagger M
// (x + 4*m2^2)^(3/8) where x = M^dagger M
// (x + 4*m2^2)^(-3/8) where x = M^dagger M
 Real A_MD_2[MOLDYN_ORDER_2+1] =	{
3.8483298363361355e-03,
2.6558360892988070e-01,
1.9569771859133908e-01,
2.7882543496824874e-01,
4.3687850968576963e-01,
7.2108707456346488e-01,
1.5251671480903228e+00
 };
 Real B_MD_2[MOLDYN_ORDER_2+1] =	{
99.9, //DUMMY
1.0978140977078202e-02,
3.5452669733142685e-02,
1.9549894369987475e-01,
1.2064028470198347e+00,
7.8144931675995775e+00,
6.2571222838716828e+01
 };
 Real A_GR_2[GRSOURCE_ORDER_2+1] =	{
1.6016869138318473e+01,
-4.8613389840934673e-04,
-4.2661397572737320e-03,
-3.0415805410802844e-02,
-2.1066081719549276e-01,
-1.4754461005108406e+00,
-1.1122407266520270e+01,
-1.1505975008324116e+02,
-6.0511740559838308e+03
 };
 Real B_GR_2[GRSOURCE_ORDER_2+1] =	{
99.9,  //DUMMY
1.5007952641723991e-02,
4.2787202777922015e-02,
1.5701948223357723e-01,
6.2040079974468609e-01,
2.5165803503208233e+00,
1.0588849834697589e+01,
5.1215275203934681e+01,
5.3595695231128923e+02
 };
 Real A_FA_2[ACTION_ORDER_2+1] =	{
6.2434174329839394e-02,
1.9725488167977902e-02,
4.1089016220267456e-02,
9.3057978626662180e-02,
2.1988308085060027e-01,
5.2940166483909257e-01,
1.3175929708128822e+00,
3.7158491203968143e+00,
1.8149606153354981e+01
 };
 Real B_FA_2[ACTION_ORDER_2+1] =	{
99.9, //DUMMY
1.1679457259936602e-02,
2.7578267012827722e-02,
9.5084864050896223e-02,
3.6909481213514034e-01,
1.4846048831793428e+00,
6.1323178474404925e+00,
2.7462784127293482e+01,
1.7974412777526393e+02
};
 

#else //not control
extern Real A_MD_1[MOLDYN_ORDER_1+1] ;
extern Real B_MD_1[MOLDYN_ORDER_1+1] ;
extern Real A_GR_1[GRSOURCE_ORDER_1+1] ;
extern Real B_GR_1[GRSOURCE_ORDER_1+1] ;
extern Real A_FA_1[ACTION_ORDER_1+1] ;
extern Real B_FA_1[ACTION_ORDER_1+1] ;
extern Real A_MD_2[MOLDYN_ORDER_2+1] ;
extern Real B_MD_2[MOLDYN_ORDER_2+1] ;
extern Real A_GR_2[GRSOURCE_ORDER_2+1] ;
extern Real B_GR_2[GRSOURCE_ORDER_2+1] ;
extern Real A_FA_2[ACTION_ORDER_2+1] ;
extern Real B_FA_2[ACTION_ORDER_2+1] ;
#endif
