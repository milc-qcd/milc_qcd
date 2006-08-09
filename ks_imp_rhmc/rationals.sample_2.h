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
#define MOLDYN_ORDER_1	8	// order of approximation for updating quark type 1 (light, usually)
#define MOLDYN_ORDER_2	0	// order of approximation for updating quark type 2 (strange, usually)
#define GRSOURCE_ORDER_1 12	// order for heat bath pseudofermion type 1
#define GRSOURCE_ORDER_2 0	// type 2
#define ACTION_ORDER_1	12	// order for evaluating fermion action type 1
#define ACTION_ORDER_2	0	// type 2
#define MAX_RAT_ORDER 12	// largest of above orders

#ifdef CONTROL
// Values for m=0.01/.05, nf=2/1  emin=1e=-15 emax=90 
// (x + 4*m1^2)^(-1/2) (x + 4*m2^2 )^(-1/4) where x = M^dagger M
// (x + 4*m1^2)^(1/4) (x + 4*m2^2 )^(1/8) where x = M^dagger M
// (x + 4*m1^2)^(-1/4) (x + 4*m2^2 )^(-1/8) where x = M^dagger M
 Real A_MD_1[MOLDYN_ORDER_1+1] =	{ 
3.7231533757634953e-03,
3.8823067852292226e-02,
7.1660822950868341e-02,
1.8180723733444168e-01,
1.9414456350624171e-01,
2.8213819645237437e-01,
4.4081487715740719e-01,
7.2355902302894326e-01,
1.5320334553134476e+00
 };
 Real B_MD_1[MOLDYN_ORDER_1+1] =	{
99.9, //DUMMY
4.8044115433228059e-04,
1.5592286846126098e-03,
7.4423198488717890e-03,
3.7182394100419500e-02,
2.2485318197011245e-01,
1.3627884462163506e+00,
8.5321493951455203e+00,
6.6308280952380500e+01
 };
 Real A_GR_1[GRSOURCE_ORDER_1+1] =	{
1.7935943463534663e+01,
-4.2346294815369983e-06,
-2.5270331053834703e-05,
-1.1836306543414480e-04,
-6.8529899477899137e-04,
-4.1969376437118880e-03,
-2.2459034790239547e-02,
-1.1833950698934990e-01,
-6.2664182568496885e-01,
-3.4220851811052952e+00,
-2.0823173512912046e+01,
-1.8590501922084979e+02,
-9.2169567420014664e+03
 };
 Real B_GR_1[GRSOURCE_ORDER_1+1] =	{
99.9,  //DUMMY
5.2899758990177350e-04,
1.2063852203952261e-03,
3.5912783720297137e-03,
1.2717949710773490e-02,
4.1183573731478476e-02,
1.3603285915322907e-01,
4.5195832007026948e-01,
1.5096906250425794e+00,
5.1139469876177301e+00,
1.8148993493158216e+01,
7.6326835865934015e+01,
7.3473857224592007e+02
 };
 Real A_FA_1[ACTION_ORDER_1+1] =	{
5.5753967001127495e-02,
9.3171482992475310e-04,
2.3277982531619045e-03,
5.8015756074727594e-03,
1.8082916607904138e-02,
3.8188779635102962e-02,
7.7636887303168853e-02,
1.6302085246807335e-01,
3.4632477082441759e-01,
7.4767386052197859e-01,
1.6924256672549463e+00,
4.4934129170867925e+00,
2.1720781252125374e+01
 };
 Real B_FA_1[ACTION_ORDER_1+1] =	{
99.9, //DUMMY
4.6411483597786962e-04,
9.4805329824817987e-04,
2.6971999057830622e-03,
8.7656333533561344e-03,
2.6398346658857257e-02,
8.6822436946333278e-02,
2.8796316220913870e-01,
9.5956003714193683e-01,
3.2271930117183918e+00,
1.1177329734840960e+01,
4.2910898957414993e+01,
2.5215453589547161e+02
 };

// dummy values - no action for second phi
 Real A_MD_2[MOLDYN_ORDER_2+1] =	{
1.0
 };
 Real B_MD_2[MOLDYN_ORDER_2+1] =	{
99.9, //DUMMY
 };
 Real A_GR_2[GRSOURCE_ORDER_2+1] =	{
1.0
 };
 Real B_GR_2[GRSOURCE_ORDER_2+1] =	{
99.9,  //DUMMY
 };
 Real A_FA_2[ACTION_ORDER_2+1] =	{
1.0
 };
 Real B_FA_2[ACTION_ORDER_2+1] =	{
99.9, //DUMMY
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
