# 1 "hybrid_loop1.c"
 
 
 














 

























 


 

 



typedef struct
{
  char name[12];                   
  int length;                      
  int disp[4];                     
  int dir[4 ];        
} link_path;

# 1 "hvy_qpot_includes.h" 1
 
 



 
# 1 "../include/config.h" 1


 

 


 
 
 

 
 
 



 
 
 

 
 
 
 

 
 
 

 
 


 
 


 


 
   
 



# 7 "hvy_qpot_includes.h" 2

# 1 "/usr/include/stdio.h" 1 3
 
 

 
 
 

 




 






#pragma ident	"@(#)stdio.h	1.69	98/07/13 SMI"	

# 1 "/usr/include/sys/feature_tests.h" 1 3
 




 
 
 




#pragma ident	"@(#)feature_tests.h	1.17	97/12/04 SMI"

# 1 "/usr/include/sys/isa_defs.h" 1 3
 







#pragma ident	"@(#)isa_defs.h	1.16	99/05/25 SMI"

 












































































































































 






# 218 "/usr/include/sys/isa_defs.h" 3


 






 






 









 
















 






 





 







 








 



# 316 "/usr/include/sys/isa_defs.h" 3


 















# 15 "/usr/include/sys/feature_tests.h" 2 3






 












 








































 



























 



















 






















 
































# 22 "/usr/include/stdio.h" 2 3

# 1 "/usr/include/sys/va_list.h" 1 3
 







#pragma ident	"@(#)va_list.h	1.11	97/11/22 SMI"

 














typedef void *__va_list;









# 23 "/usr/include/stdio.h" 2 3

# 1 "/usr/include/stdio_tag.h" 1 3
 







#pragma ident	"@(#)stdio_tag.h	1.3	98/04/20 SMI"











typedef struct __FILE  __FILE;







# 24 "/usr/include/stdio.h" 2 3

# 1 "/usr/include/stdio_impl.h" 1 3
 







#pragma ident	"@(#)stdio_impl.h	1.8	99/06/30 SMI"












typedef int	ssize_t;		 



# 36 "/usr/include/stdio_impl.h" 3


struct __FILE 	 
{




	ssize_t		_cnt;	 
	unsigned char	*_ptr;	 

	unsigned char	*_base;	 
	unsigned char	_flag;	 
	unsigned char	_file;	 
	unsigned	__orientation:2;  
	unsigned	__ionolock:1;	 
	unsigned	__filler:5;
};








# 25 "/usr/include/stdio.h" 2 3


 













typedef	__FILE FILE;







typedef unsigned int	size_t;		 




 





typedef	long long	__longlong_t;













typedef long		off_t;







typedef __longlong_t	off64_t;






typedef long		fpos_t;







typedef __longlong_t	fpos64_t;













 










 


































































extern FILE	__iob[20 ];











extern unsigned char	 _sibuf[], _sobuf[];


 
# 227 "/usr/include/stdio.h" 3



extern unsigned char	*_bufendtab[];
extern FILE		*_lastbuf;


 
# 257 "/usr/include/stdio.h" 3




extern int	remove(const char *);
extern int	rename(const char *, const char *);
extern FILE	*tmpfile(void);
extern char	*tmpnam(char *);



extern int	fclose(FILE *);
extern int	fflush(FILE *);
extern FILE	*fopen(const char *, const char *);
extern FILE	*freopen(const char *, const char *, FILE *);
extern void	setbuf(FILE *, char *);


extern void setbuffer(FILE *, char *, size_t);
extern int setlinebuf(FILE *);

extern int	setvbuf(FILE *, char *, int, size_t);
 
extern int	fprintf(FILE *, const char *, ...);
 
extern int	fscanf(FILE *, const char *, ...);
 
extern int	printf(const char *, ...);
 
extern int	scanf(const char *, ...);



 
extern int	snprintf(char *, size_t, const char *, ...);

 
extern int	sprintf(char *, const char *, ...);
 
extern int	sscanf(const char *, const char *, ...);
extern int	vfprintf(FILE *, const char *, __va_list);
extern int	vprintf(const char *, __va_list);



extern int	vsnprintf(char *, size_t, const char *, __va_list);

extern int	vsprintf(char *, const char *, __va_list);
extern int	fgetc(FILE *);
extern char	*fgets(char *, int, FILE *);
extern int	fputc(int, FILE *);
extern int	fputs(const char *, FILE *);
extern int	getc(FILE *);
extern int	getchar(void);
extern char	*gets(char *);
extern int	putc(int, FILE *);
extern int	putchar(int);
extern int	puts(const char *);
extern int	ungetc(int, FILE *);
extern size_t	fread(void *, size_t, size_t, FILE *);
extern size_t	fwrite(const void *, size_t, size_t, FILE *);
extern int	fgetpos(FILE *, fpos_t *);
extern int	fseek(FILE *, long, int);
extern int	fsetpos(FILE *, const fpos_t *);
extern long	ftell(FILE *);
extern void	rewind(FILE *);
extern void	clearerr(FILE *);
extern int	feof(FILE *);
extern int	ferror(FILE *);
extern void	perror(const char *);


extern int	__filbuf(FILE *);
extern int	__flsbuf(int, FILE *);


 





extern FILE	*fdopen(int, const char *);
extern char	*ctermid(char *);
extern int	fileno(FILE *);



 


# 358 "/usr/include/stdio.h" 3


 




extern FILE	*popen(const char *, const char *);
extern char	*cuserid(char *);
extern char	*tempnam(const char *, const char *);
extern int	getopt(int, char *const *, const char *);

extern int	getsubopt(char **, char *const *, char **);

extern char	*optarg;
extern int	optind, opterr, optopt;
extern int	getw(FILE *);
extern int	putw(int, FILE *);
extern int	pclose(FILE *);



 



extern int	fseeko(FILE *, off_t, int);
extern off_t	ftello(FILE *);


 





extern FILE	*fopen64(const char *, const char *);
extern FILE	*freopen64(const char *, const char *, FILE *);
extern FILE	*tmpfile64(void);
extern int	fgetpos64(FILE *, fpos64_t *);
extern int	fsetpos64(FILE *, const fpos64_t *);
extern int	fseeko64(FILE *, off64_t, int);
extern off64_t	ftello64(FILE *);


# 515 "/usr/include/stdio.h" 3





























# 567 "/usr/include/stdio.h" 3









# 8 "hvy_qpot_includes.h" 2

# 1 "/usr/include/stdlib.h" 1 3
 




 
 

 
 
 




#pragma ident	"@(#)stdlib.h	1.44	98/01/22 SMI"	












typedef	struct {
	int	quot;
	int	rem;
} div_t;

typedef struct {
	long	quot;
	long	rem;
} ldiv_t;


typedef struct {
	long long	quot;
	long long	rem;
} lldiv_t;
















typedef long	uid_t;			 




















typedef long	wchar_t;





 
# 101 "/usr/include/stdlib.h" 3


 
# 117 "/usr/include/stdlib.h" 3


extern unsigned char	__ctype[];



extern double atof(const char *);
extern int atoi(const char *);
extern long int atol(const char *);
extern double strtod(const char *, char **);
extern long int strtol(const char *, char **, int);
extern unsigned long int strtoul(const char *, char **, int);

extern int rand(void);
extern void srand(unsigned int);





extern void *calloc(size_t, size_t);
extern void free(void *);
extern void *malloc(size_t);
extern void *realloc(void *, size_t);

extern void abort(void);
extern int atexit(void (*)(void));
extern void exit(int);
extern void _exithandle(void);
extern char *getenv(const char *);
extern int system(const char *);

extern void *bsearch(const void *, const void *, size_t, size_t,
	int (*)(const void *, const void *));
extern void qsort(void *, size_t, size_t,
	int (*)(const void *, const void *));

extern int abs(int);
extern div_t div(int, int);
extern long int labs(long);
extern ldiv_t ldiv(long, long);

extern int mbtowc(wchar_t *, const char *, size_t);
extern int mblen(const char *, size_t);
extern int wctomb(char *, wchar_t);

extern size_t mbstowcs(wchar_t *, const char *, size_t);
extern size_t wcstombs(char *, const wchar_t *, size_t);




extern double drand48(void);
extern double erand48(unsigned short *);
extern long jrand48(unsigned short *);
extern void lcong48(unsigned short *);
extern long lrand48(void);
extern long mrand48(void);
extern long nrand48(unsigned short *);
extern unsigned short *seed48(unsigned short *);
extern void srand48(long);
extern int putenv(char *);
extern void setkey(const char *);





extern void swab(const char *, char *, int);




extern int	mkstemp(char *);


extern int	mkstemp64(char *);






extern long a64l(const char *);
extern char *ecvt(double, int, int *, int *);
extern char *fcvt(double, int, int *, int *);
extern char *gcvt(double, int, char *);
extern int getsubopt(char **, char *const *, char **);
extern int  grantpt(int);
extern char *initstate(unsigned, char *, size_t);
extern char *l64a(long);
extern char *mktemp(char *);
extern char *ptsname(int);
extern long random(void);
extern char *realpath(const char *, char *);
extern char *setstate(const char *);
extern void srandom(unsigned);
extern int ttyslot(void);
extern int  unlockpt(int);
extern void *valloc(size_t);




extern int dup2(int, int);
extern char *qecvt(long double, int, int *, int *);
extern char *qfcvt(long double, int, int *, int *);
extern char *qgcvt(long double, int, char *);
extern char *getcwd(char *, size_t);
extern const char *getexecname(void);
extern char *getlogin(void);
extern int getopt(int, char *const *, const char *);
extern char *optarg;
extern int optind, opterr, optopt;
extern char *getpass(const char *);
extern char *getpassphrase(const char *);
extern int getpw(uid_t, char *);
extern int isatty(int);
extern void *memalign(size_t, size_t);
extern char *ttyname(int);


extern long long atoll(const char *);
extern long long llabs(long long);
extern lldiv_t lldiv(long long, long long);
extern char *lltostr(long long, char *);
extern long long strtoll(const char *, char **, int);
extern unsigned long long strtoull(const char *, char **, int);
extern char *ulltostr(unsigned long long, char *);




# 380 "/usr/include/stdlib.h" 3







# 9 "hvy_qpot_includes.h" 2

# 1 "/usr/local/lib/gcc-lib/sparc-sun-solaris2.7/2.95.3/include/math.h" 1 3




# 1 "/usr/include/math.h" 1 3
 










#pragma ident	"@(#)math.h	2.7	98/01/27"











 


typedef union _h_val {
  	unsigned long _i[sizeof(double) / sizeof(unsigned long)];
	double _d;
} _h_val;


extern const _h_val __huge_val;








 
















extern int signgam;




 


enum version {libm_ieee = -1, c_issue_4, ansi_1, strict_ansi};


extern const enum version _lib_version;




struct exception {
	int type;
	char *name;
	double arg1;
	double arg2;
	double retval;
};




























 


extern double acos  (double)  ;
extern double asin  (double)  ;
extern double atan  (double)  ;
extern double atan2  (double, double)  ;
extern double cos  (double)  ;
extern double sin  (double)  ;
extern double tan  (double)  ;

extern double cosh  (double)  ;
extern double sinh  (double)  ;
extern double tanh  (double)  ;

extern double exp  (double)  ;
extern double frexp  (double, int *)  ;
extern double ldexp  (double, int)  ;
extern double log  (double)  ;
extern double log10  (double)  ;
extern double modf  (double, double *)  ;

extern double pow  (double, double)  ;
extern double sqrt  (double)  ;

extern double ceil  (double)  ;
extern double fabs  (double)  ;
extern double floor  (double)  ;
extern double fmod  (double, double)  ;



 


extern double erf  (double)  ;
extern double erfc  (double)  ;
extern double gamma  (double)  ;
extern double hypot  (double, double)  ;
extern int isnan  (double)  ;
extern double j0  (double)  ;
extern double j1  (double)  ;
extern double jn  (int, double)  ;
extern double lgamma  (double)  ;
extern double y0  (double)  ;
extern double y1  (double)  ;
extern double yn  (int, double)  ;




 


extern double acosh  (double)  ;
extern double asinh  (double)  ;
extern double atanh  (double)  ;
extern double cbrt  (double)  ;
extern double logb  (double)  ;
extern double nextafter  (double, double)  ;
extern double remainder  (double, double)  ;
extern double scalb  (double, double)  ;

 


extern double expm1  (double)  ;
extern int ilogb  (double)  ;
extern double log1p  (double)  ;
extern double rint  (double)  ;



 


extern int matherr  (struct exception *)  ;

 


extern double significand  (double)  ;

 


extern double copysign  (double, double)  ;
extern double scalbn  (double, int)  ;

 








 


extern float modff  (float, float *)  ;

# 1 "/usr/include/floatingpoint.h" 1 3
 
 

 
 
 








#pragma ident	"@(#)floatingpoint.h	2.4 94/06/09"

 




 












# 1 "/usr/include/sys/ieeefp.h" 1 3
 






#pragma ident	"@(#)ieeefp.h	2.7 94/11/09"





 



enum fp_direction_type {	 
	fp_nearest	= 0,
	fp_tozero	= 1,
	fp_positive	= 2,
	fp_negative	= 3
};

enum fp_precision_type {	 
	fp_extended	= 0,
	fp_single	= 1,
	fp_double	= 2,
	fp_precision_3	= 3
};

enum fp_exception_type {	 
	fp_inexact	= 0,
	fp_division	= 1,
	fp_underflow	= 2,
	fp_overflow	= 3,
	fp_invalid	= 4
};

enum fp_trap_enable_type {	 
	fp_trap_inexact	= 0,
	fp_trap_division	= 1,
	fp_trap_underflow	= 2,
	fp_trap_overflow	= 3,
	fp_trap_invalid	= 4
};


# 81 "/usr/include/sys/ieeefp.h" 3


# 122 "/usr/include/sys/ieeefp.h" 3


enum fp_class_type {		 
	fp_zero		= 0,
	fp_subnormal	= 1,
	fp_normal	= 2,
	fp_infinity   	= 3,
	fp_quiet	= 4,
	fp_signaling	= 5
};






# 35 "/usr/include/floatingpoint.h" 2 3










typedef int sigfpe_code_type;	 

typedef void (*sigfpe_handler_type)();	 





extern sigfpe_handler_type sigfpe  (sigfpe_code_type, sigfpe_handler_type)  ;

 


typedef float single;			



typedef unsigned extended[3];


typedef long double quadruple;	 

typedef unsigned fp_exception_field_type;
				 



 




typedef char decimal_string[512 ];	
				 

typedef struct {
	enum fp_class_type fpclass;
	int	sign;
	int	exponent;
	decimal_string ds;	 


	int	more;		 


	int 	ndigits;	 


} decimal_record;

enum decimal_form {
	fixed_form,		 


	floating_form		 

};

typedef struct {
	enum fp_direction_type rd;	
				 
	enum decimal_form df;	 

	int ndigits;		 
} decimal_mode;

enum decimal_string_form {	 
	invalid_form,		 
	whitespace_form,	 
	fixed_int_form,		 
	fixed_intdot_form,	 
	fixed_dotfrac_form,	 
	fixed_intdotfrac_form,	 
	floating_int_form,	 	
	floating_intdot_form,	 
	floating_dotfrac_form,	 
	floating_intdotfrac_form, 
	inf_form,		 
	infinity_form,		 
	nan_form,		 
	nanstring_form		 
};

extern void single_to_decimal  (single *, decimal_mode *, decimal_record *,
				fp_exception_field_type *)  ;
extern void double_to_decimal  (double *, decimal_mode *, decimal_record *,
				fp_exception_field_type *)  ;
extern void extended_to_decimal  (extended *, decimal_mode *,
				decimal_record *, fp_exception_field_type *)  ;
extern void quadruple_to_decimal  (quadruple *, decimal_mode *,
				decimal_record *, fp_exception_field_type *)  ;

extern void decimal_to_single  (single *, decimal_mode *, decimal_record *,
				fp_exception_field_type *)  ;
extern void decimal_to_double  (double *, decimal_mode *, decimal_record *,
				fp_exception_field_type *)  ;
extern void decimal_to_extended  (extended *, decimal_mode *,
				decimal_record *, fp_exception_field_type *)  ;
extern void decimal_to_quadruple  (quadruple *, decimal_mode *,
				decimal_record *, fp_exception_field_type *)  ;

extern void string_to_decimal  (char **, int, int, decimal_record *,
				enum decimal_string_form *, char **)  ;
extern void func_to_decimal  (char **, int, int, decimal_record *,
				enum decimal_string_form *, char **,
				int (*)(void), int *, int (*)(int))  ;
extern void file_to_decimal  (char **, int, int, decimal_record *,
				enum decimal_string_form *, char **,
				FILE *, int *)  ;

extern char *seconvert  (single *, int, int *, int *, char *)  ;
extern char *sfconvert  (single *, int, int *, int *, char *)  ;
extern char *sgconvert  (single *, int, int, char *)  ;
extern char *econvert  (double, int, int *, int *, char *)  ;
extern char *fconvert  (double, int, int *, int *, char *)  ;
extern char *gconvert  (double, int, int, char *)  ;
extern char *qeconvert  (quadruple *, int, int *, int *, char *)  ;
extern char *qfconvert  (quadruple *, int, int *, int *, char *)  ;
extern char *qgconvert  (quadruple *, int, int, char *)  ;

extern char *ecvt  (double, int, int *, int *)  ;
extern char *fcvt  (double, int, int *, int *)  ;
extern char *gcvt  (double, int, char *)  ;

 



extern double atof  (const char *)  ;
extern double strtod  (const char *, char **)  ;






# 213 "/usr/include/math.h" 2 3









# 5 "/usr/local/lib/gcc-lib/sparc-sun-solaris2.7/2.95.3/include/math.h" 2 3






# 10 "hvy_qpot_includes.h" 2

# 1 "../include/complex.h" 1



 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 

 






 







typedef struct {	    
   float real;		    
   float imag;
} complex;
typedef struct {            
   double real;		    
   double imag;
} double_complex;

 
 






 
 


 
complex cmplx(  float x, float y );
complex cadd( complex *a, complex *b );
complex cmul( complex *a, complex *b );
complex csub( complex *a, complex *b );
complex cdiv( complex *a, complex *b );
complex conjg( complex *a );
complex cexp( complex *a );       
complex clog( complex *a );        
complex csqrt( complex *z );        
complex ce_itheta( float theta );    

double_complex dcmplx( double x, double y );
double_complex dcadd( double_complex *a, double_complex *b );
double_complex dcmul( double_complex *a, double_complex *b );
double_complex dcsub( double_complex *a, double_complex *b );
double_complex dcdiv( double_complex *a, double_complex *b );
double_complex dconjg(  double_complex *a );
double_complex dcexp(  double_complex *a ); 
double_complex dclog(  double_complex *a );  
double_complex dcsqrt( double_complex *z );  
double_complex dce_itheta( double theta );   

 
								 

								 

								 


								 

								 


								 


								 

								 


								 


								 



								 


								 


								 


								 

								 

								 

								 

								 


                                                                


                                                                 




# 11 "hvy_qpot_includes.h" 2

# 1 "../include/su3.h" 1




# 1 "../include/../include/random.h" 1



 

typedef struct {
   
  unsigned long r0,r1,r2,r3,r4,r5,r6;
  unsigned long multiplier,addend,ic_state;
  float scale;
} double_prn;

 

float myrand(double_prn *prn_pt);


# 5 "../include/su3.h" 2


 





 
typedef struct { complex e[3][3]; } su3_matrix;
typedef struct { complex c[3]; } su3_vector;
typedef struct
  { complex m01,m02,m12; float m00im,m11im,m22im; float space; } anti_hermitmat;

 
typedef struct { complex e[2][2]; } su2_matrix;

 

 
 
 
 
 
 
 
 
 
 
 
 
 
 

 






typedef struct { su3_vector d[4]; } wilson_vector;
typedef struct { su3_vector h[2]; } half_wilson_vector;
typedef struct { wilson_vector c[3]; } color_wilson_vector;
typedef struct { wilson_vector d[4]; } spin_wilson_vector;
typedef struct { color_wilson_vector d[4]; } wilson_matrix;
typedef struct { spin_wilson_vector c[3]; } wilson_propagator;




 





 







































































































































































































































































































void mult_su3_nn ( su3_matrix *a, su3_matrix *b, su3_matrix *c );
void mult_su3_na ( su3_matrix *a, su3_matrix *b, su3_matrix *c );
void mult_su3_an ( su3_matrix *a, su3_matrix *b, su3_matrix *c );
float realtrace_su3(  su3_matrix *a, su3_matrix *b );
complex trace_su3(  su3_matrix *a );
complex complextrace_su3( su3_matrix *a, su3_matrix *b );
complex det_su3( su3_matrix *a );
void add_su3_matrix( su3_matrix *a, su3_matrix *b, su3_matrix *c );
void sub_su3_matrix( su3_matrix *a, su3_matrix *b, su3_matrix *c );
void scalar_mult_su3_matrix( su3_matrix *src, float scalar, su3_matrix *dest);
void scalar_mult_add_su3_matrix( su3_matrix *src1, su3_matrix *src2,
	float scalar, su3_matrix *dest);
void scalar_mult_sub_su3_matrix( su3_matrix *src1, su3_matrix *src2,
	float scalar, su3_matrix *dest);
void c_scalar_mult_su3mat( su3_matrix *src, complex *scalar,
	su3_matrix *dest);
void c_scalar_mult_add_su3mat( su3_matrix *src1, su3_matrix *src2,
	complex *scalar, su3_matrix *dest);
void c_scalar_mult_sub_su3mat( su3_matrix *src1, su3_matrix *src2,
	complex *scalar, su3_matrix *dest);
void su3_adjoint( su3_matrix *a, su3_matrix *b );
void make_anti_hermitian( su3_matrix *m3, anti_hermitmat *ah3 );
void random_anti_hermitian( anti_hermitmat *mat_antihermit, double_prn *prn_pt );
void uncompress_anti_hermitian( anti_hermitmat *mat_anti, su3_matrix *mat );
void compress_anti_hermitian( su3_matrix *mat, anti_hermitmat *mat_anti);
void clear_su3mat( su3_matrix *dest );
void su3mat_copy( su3_matrix *a, su3_matrix *b );
void dumpmat( su3_matrix *m );

void su3_projector( su3_vector *a, su3_vector *b, su3_matrix *c );
complex su3_dot( su3_vector *a, su3_vector *b );
float su3_rdot( su3_vector *a, su3_vector *b );
float magsq_su3vec( su3_vector *a );
void su3vec_copy( su3_vector *a, su3_vector *b );
void dumpvec( su3_vector *v );
void clearvec( su3_vector *v );

void mult_su3_mat_vec( su3_matrix *a, su3_vector *b, su3_vector *c );
void mult_su3_mat_vec_sum(  su3_matrix *a, su3_vector *b, su3_vector *c );
void mult_su3_mat_vec_sum_4dir( su3_matrix *a, su3_vector *b0,
	su3_vector *b1, su3_vector *b2, su3_vector *b3, su3_vector *c );
void mult_su3_mat_vec_nsum( su3_matrix *a, su3_vector *b, su3_vector *c );
void mult_adj_su3_mat_vec( su3_matrix *a, su3_vector *b, su3_vector *c );
void mult_adj_su3_mat_vec_4dir( su3_matrix *a, su3_vector *b, su3_vector *c );
void mult_adj_su3_mat_4vec( su3_matrix *mat, su3_vector *src,
			    su3_vector *dest0, su3_vector *dest1, 
			    su3_vector *dest2, su3_vector *dest3  ) ;
void mult_adj_su3_mat_vec_sum( su3_matrix *a, su3_vector *b, su3_vector *c );
void mult_adj_su3_mat_vec_nsum( su3_matrix *a, su3_vector *b, su3_vector *c );

void add_su3_vector( su3_vector *a, su3_vector *b, su3_vector *c );
void sub_su3_vector( su3_vector *a, su3_vector *b, su3_vector *c );
void sub_four_su3_vecs( su3_vector *a, su3_vector *b1, su3_vector *b2,
	su3_vector *b3, su3_vector *b4 );

void scalar_mult_su3_vector(  su3_vector *src, float scalar, 
	su3_vector *dest);
void scalar_mult_add_su3_vector( su3_vector *src1, su3_vector *src2,
	float scalar, su3_vector *dest);
void scalar_mult_sum_su3_vector( su3_vector *src1, su3_vector *src2,
	float scalar);
void scalar_mult_sub_su3_vector( su3_vector *src1, su3_vector *src2,
	float scalar, su3_vector *dest);
void scalar_mult_wvec( wilson_vector *src, float s, wilson_vector *dest );
void scalar_mult_hwvec( half_wilson_vector *src, float s, 
    half_wilson_vector *dest );
void scalar_mult_add_wvec( wilson_vector *src1, wilson_vector *src2,
	float scalar, wilson_vector *dest );
void scalar_mult_addtm_wvec( wilson_vector *src1, wilson_vector *src2,
	float scalar, wilson_vector *dest );
void c_scalar_mult_wvec(wilson_vector *src1, complex *phase,
	wilson_vector *dest );
void c_scalar_mult_add_wvec(wilson_vector *src1, wilson_vector *src2,
	complex *phase, wilson_vector *dest );
void c_scalar_mult_add_wvec2(wilson_vector *src1, wilson_vector *src2,
	complex s, wilson_vector *dest );
void c_scalar_mult_su3vec( su3_vector *src, complex *phase, su3_vector *dest );
void c_scalar_mult_add_su3vec(su3_vector *v1, complex *phase, su3_vector *v2);
void c_scalar_mult_sub_su3vec(su3_vector *v1, complex *phase, su3_vector *v2);

void left_su2_hit_n(su2_matrix *u,int p,int q,su3_matrix *link);
void right_su2_hit_a(su2_matrix *u,int p,int q,su3_matrix *link);
void dumpsu2(su2_matrix *u);
void mult_su2_mat_vec_elem_n(su2_matrix *u,complex *x0,complex *x1);
void mult_su2_mat_vec_elem_a(su2_matrix *u,complex *x0,complex *x1);

void mult_mat_wilson_vec( su3_matrix *mat, wilson_vector *src,
	wilson_vector *dest );
void mult_su3_mat_hwvec( su3_matrix *mat, half_wilson_vector *src,
	half_wilson_vector *dest );
void mult_adj_mat_wilson_vec( su3_matrix *mat, wilson_vector *src,
	wilson_vector *dest);
void mult_adj_su3_mat_hwvec( su3_matrix *mat, half_wilson_vector *src,
	half_wilson_vector *dest );

void add_wilson_vector( wilson_vector *src1, wilson_vector *src2,
	wilson_vector *dest );
void sub_wilson_vector( wilson_vector *src1, wilson_vector *src2,
       wilson_vector *dest );
float magsq_wvec( wilson_vector *src );
complex wvec_dot( wilson_vector *src1, wilson_vector *src2 );
complex wvec2_dot( wilson_vector *src1, wilson_vector *src2 );
float wvec_rdot( wilson_vector *a, wilson_vector *b );

void wp_shrink( wilson_vector *src, half_wilson_vector *dest,
	int dir, int sign );
void wp_shrink_4dir( wilson_vector *a,  half_wilson_vector *b1,
	half_wilson_vector *b2, half_wilson_vector *b3,
	half_wilson_vector *b4,	int sign );
void wp_grow(  half_wilson_vector *src, wilson_vector *dest,
	int dir, int sign );
void wp_grow_add( half_wilson_vector *src, wilson_vector *dest,
	int dir, int sign );
void grow_add_four_wvecs( wilson_vector *a, half_wilson_vector *b1,
	half_wilson_vector *b2, half_wilson_vector *b3,
	half_wilson_vector *b4, int sign, int sum );
void mult_by_gamma( wilson_vector *src, wilson_vector *dest, int dir );
void mult_by_gamma_left( wilson_matrix *src,  wilson_matrix *dest, int dir );
void mult_by_gamma_right( wilson_matrix *src,  wilson_matrix *dest, int dir );
void mult_swv_by_gamma_l(spin_wilson_vector *src, spin_wilson_vector *dest, int dir);
void mult_swv_by_gamma_r(spin_wilson_vector *src, spin_wilson_vector *dest, int dir);
void su3_projector_w( wilson_vector *a, wilson_vector *b, su3_matrix *c );
void clear_wvec( wilson_vector *dest );
void copy_wvec( wilson_vector *src, wilson_vector *dest );
void dump_wilson_vec( wilson_vector *src );

float gaussian_rand_no( double_prn *prn_pt );
# 1 "../include/../include/int32type.h" 1
 



 










 














typedef int int32type;
typedef unsigned int u_int32type;




# 485 "../include/su3.h" 2

void byterevn(int32type w[], int n);

# 12 "hvy_qpot_includes.h" 2

# 1 "../include/macros.h" 1



# 1 "defines.h" 1



 

 	 

# 4 "../include/macros.h" 2


 

 
 



 
 





 
 



 
 





 
 

 
 





typedef int field_offset;




 
 








 
 


 



 
 

 
 

 


 

 

 





















 
# 123 "../include/macros.h"



# 135 "../include/macros.h"








# 13 "hvy_qpot_includes.h" 2

# 1 "lattice.h" 1


 

 





# 1 "../include/io_lat.h" 1


 

 


 
 



 















 







# 1 "/usr/include/unistd.h" 1 3
 
 

 
 
 

 







#pragma ident	"@(#)unistd.h	1.55	98/04/14 SMI"	



# 1 "/usr/include/sys/types.h" 1 3
 
 

 
 
 

 







#pragma ident	"@(#)types.h	1.63	98/04/28 SMI"




 


# 1 "/usr/include/sys/machtypes.h" 1 3
 
 

 
 
 

 







#pragma ident	"@(#)machtypes.h	1.12	96/08/01 SMI"







 








typedef struct  _physadr_t { int r[1]; } *physadr_t;

typedef	struct	_label_t { long val[2]; } label_t;



typedef	unsigned char	lock_t;		 






# 24 "/usr/include/sys/types.h" 2 3


 









# 1 "/usr/include/sys/int_types.h" 1 3
 







#pragma ident	"@(#)int_types.h	1.6	97/08/20 SMI"

 




































 













typedef char			int8_t;





typedef short			int16_t;
typedef int			int32_t;




typedef	long long		int64_t;



typedef unsigned char		uint8_t;
typedef unsigned short		uint16_t;
typedef unsigned int		uint32_t;




typedef unsigned long long	uint64_t;



 




typedef int64_t			intmax_t;
typedef uint64_t		uintmax_t;





 








typedef	int			intptr_t;
typedef	unsigned int		uintptr_t;


 




typedef char			int_least8_t;





typedef short			int_least16_t;
typedef int			int_least32_t;




typedef long long		int_least64_t;



typedef unsigned char		uint_least8_t;
typedef unsigned short		uint_least16_t;
typedef unsigned int		uint_least32_t;




typedef unsigned long long	uint_least64_t;








# 36 "/usr/include/sys/types.h" 2 3











 





typedef	long long		longlong_t;
typedef	unsigned long long	u_longlong_t;
# 65 "/usr/include/sys/types.h" 3


 








typedef long		t_scalar_t;	 
typedef unsigned long	t_uscalar_t;


 


typedef	unsigned char	uchar_t;
typedef	unsigned short	ushort_t;
typedef	unsigned int	uint_t;
typedef	unsigned long	ulong_t;

typedef	char		*caddr_t;	 
typedef	long		daddr_t;	 
typedef	short		cnt_t;		 


typedef	ulong_t		paddr_t;	 







typedef	int	ptrdiff_t;		 



 


typedef	ulong_t		pfn_t;		 
typedef	ulong_t		pgcnt_t;	 
typedef	long		spgcnt_t;	 

typedef	uchar_t		use_t;		 
typedef	short		sysid_t;
typedef	short		index_t;
typedef void		*timeout_id_t;	 
typedef void		*bufcall_id_t;	 

 







# 143 "/usr/include/sys/types.h" 3



typedef ulong_t		ino_t;		 
typedef long		blkcnt_t;	 
typedef ulong_t		fsblkcnt_t;	 
typedef ulong_t		fsfilcnt_t;	 














typedef u_longlong_t	ino64_t;	 
typedef longlong_t	blkcnt64_t;	 
typedef u_longlong_t	fsblkcnt64_t;	 
typedef u_longlong_t	fsfilcnt64_t;	 






typedef	long		blksize_t;	 





typedef enum { B_FALSE, B_TRUE } boolean_t;


 







typedef int64_t		pad64_t;
typedef	uint64_t	upad64_t;
# 202 "/usr/include/sys/types.h" 3


typedef	longlong_t	offset_t;
typedef	u_longlong_t	u_offset_t;
typedef u_longlong_t	len_t;
typedef	longlong_t	diskaddr_t;

 




# 222 "/usr/include/sys/types.h" 3



typedef union {
	offset_t	_f;	 
	struct {
		int32_t	_u;	 
		int32_t	_l;	 
	} _p;
} lloff_t;


# 242 "/usr/include/sys/types.h" 3



typedef union {
	diskaddr_t	_f;	 
	struct {
		int32_t	_u;	 
		int32_t	_l;	 
	} _p;
} lldaddr_t;


typedef uint_t k_fltset_t;	 

 













typedef	long		id_t;		 


 



typedef uint_t		useconds_t;	 



typedef long	suseconds_t;	 


 






typedef ulong_t	major_t;	 
typedef ulong_t	minor_t;	 


 


typedef short	pri_t;

 










typedef	ushort_t o_mode_t;		 
typedef short	o_dev_t;		 
typedef	ushort_t o_uid_t;		 
typedef	o_uid_t	o_gid_t;		 
typedef	short	o_nlink_t;		 
typedef short	o_pid_t;		 
typedef ushort_t o_ino_t;		 


 


typedef	int	key_t;			 



typedef	ulong_t	mode_t;			 











typedef	uid_t	gid_t;			 

 




typedef	unsigned int	pthread_t;	 
typedef	unsigned int	pthread_key_t;	 

typedef	struct	_pthread_mutex {		 
	struct {
		uint8_t		__pthread_mutex_flag[4];
		uint32_t 	__pthread_mutex_type;
	} __pthread_mutex_flags;
	union {
		struct {
			uint8_t	__pthread_mutex_pad[8];
		} __pthread_mutex_lock64;
		upad64_t __pthread_mutex_owner64;
	} __pthread_mutex_lock;
	upad64_t __pthread_mutex_data;
} pthread_mutex_t;

typedef	struct	_pthread_cond {		 
	struct {
		uint8_t		__pthread_cond_flag[4];
		uint32_t 	__pthread_cond_type;
	} __pthread_cond_flags;
	upad64_t __pthread_cond_data;
} pthread_cond_t;

 


typedef	struct _pthread_rwlock {	 
	int32_t		__pthread_rwlock_readers;
	uint16_t	__pthread_rwlock_type;
	uint16_t	__pthread_rwlock_magic;
	upad64_t	__pthread_rwlock_pad1[3];
	upad64_t	__pthread_rwlock_pad2[2];
	upad64_t	__pthread_rwlock_pad3[2];
} pthread_rwlock_t;

 


typedef struct _pthread_attr {
	void	*__pthread_attrp;
} pthread_attr_t;


 


typedef struct _pthread_mutexattr {
	void	*__pthread_mutexattrp;
} pthread_mutexattr_t;


 


typedef struct _pthread_condattr {
	void	*__pthread_condattrp;
} pthread_condattr_t;

 


typedef	struct	_once {
	upad64_t	__pthread_once_pad[4];
} pthread_once_t;

 



typedef struct _pthread_rwlockattr {
	void	*__pthread_rwlockattrp;
} pthread_rwlockattr_t;

typedef ulong_t	dev_t;			 





typedef	ulong_t	nlink_t;		 
typedef	long	pid_t;			 






















typedef	long		time_t;	 




typedef	long		clock_t;  




typedef	int	clockid_t;	 




typedef	int	timer_t;	 





 
typedef	unsigned char	unchar;
typedef	unsigned short	ushort;
typedef	unsigned int	uint;
typedef	unsigned long	ulong;
 

# 501 "/usr/include/sys/types.h" 3




 














 





 
typedef unsigned char	u_char;
typedef unsigned short	u_short;
typedef unsigned int	u_int;
typedef unsigned long	u_long;
typedef struct _quad { int val[2]; } quad_t;	 
typedef quad_t		quad;			 
 

 



# 1 "/usr/include/sys/select.h" 1 3
 
 
 

 
 
 




#pragma ident	"@(#)select.h	1.16	98/04/27 SMI"	




# 1 "/usr/include/sys/time.h" 1 3
 
 

 
 
 

 





 







#pragma ident	"@(#)time.h	2.62	98/02/04 SMI"	



 






















struct timeval {
	time_t		tv_sec;		 
	suseconds_t	tv_usec;	 
};

# 74 "/usr/include/sys/time.h" 3








struct timezone {
	int	tz_minuteswest;	 
	int	tz_dsttime;	 
};








 





# 1 "/usr/include/sys/types.h" 1 3
 
 

 
 
 

 




# 559 "/usr/include/sys/types.h" 3

# 100 "/usr/include/sys/time.h" 2 3






















 


















 






				 
				 

				 
				 


struct	itimerval {
	struct	timeval it_interval;	 
	struct	timeval it_value;	 
};

# 181 "/usr/include/sys/time.h" 3







 











 


typedef	longlong_t	hrtime_t;

# 310 "/usr/include/sys/time.h" 3





int adjtime(struct timeval *, struct timeval *);










int getitimer(int, struct itimerval *);
int utimes(const char *, const struct timeval *);



int setitimer(int, struct itimerval *, struct itimerval *);











 




















int settimeofday(struct timeval *, void *);

hrtime_t	gethrtime(void);
hrtime_t	gethrvtime(void);
















int gettimeofday(struct timeval *, void *);







 












# 1 "/usr/include/time.h" 1 3
 
 

 
 
 

 







#pragma ident	"@(#)time.h	1.36	98/03/27 SMI"	


# 1 "/usr/include/sys/types.h" 1 3
 
 

 
 
 

 




# 559 "/usr/include/sys/types.h" 3

# 19 "/usr/include/time.h" 2 3













































struct	tm {	 
	int	tm_sec;
	int	tm_min;
	int	tm_hour;
	int	tm_mday;
	int	tm_mon;
	int	tm_year;
	int	tm_wday;
	int	tm_yday;
	int	tm_isdst;
};



extern clock_t clock(void);
extern double difftime(time_t, time_t);
extern time_t mktime(struct tm *);
extern time_t time(time_t *);
extern char *asctime(const struct tm *);
extern char *ctime(const time_t *);
extern struct tm *gmtime(const time_t *);
extern struct tm *localtime(const time_t *);
extern size_t strftime(char *, size_t, const char *, const struct tm *);










extern char *strptime(const char *, const char *, struct tm *);





# 1 "/usr/include/sys/time_impl.h" 1 3
 




 








#pragma ident	"@(#)time_impl.h	1.3	98/02/13 SMI"














 





typedef struct  timespec {		 
	time_t		tv_sec;		 
	long		tv_nsec;	 
} timespec_t;

# 58 "/usr/include/sys/time_impl.h" 3


typedef struct timespec timestruc_t;	 



 









 




typedef struct itimerspec {		 
	struct timespec	it_interval;	 
	struct timespec	it_value;	 
} itimerspec_t;

# 100 "/usr/include/sys/time_impl.h" 3













 




















# 103 "/usr/include/time.h" 2 3


 












union sigval {
	int	sival_int;	 
	void	*sival_ptr;	 
};




struct sigevent {
	int		sigev_notify;	 
	int		sigev_signo;	 
	union sigval	sigev_value;	 
	void		(*sigev_notify_function)(union sigval);
	pthread_attr_t	*sigev_notify_attributes;
	int		__sigev_pad2;
};


extern int clock_getres(clockid_t, struct timespec *);
extern int clock_gettime(clockid_t, struct timespec *);
extern int clock_settime(clockid_t, const struct timespec *);
extern int timer_create(clockid_t, struct sigevent *, timer_t *);
extern int timer_delete(timer_t);
extern int timer_getoverrun(timer_t);
extern int timer_gettime(timer_t, struct itimerspec *);
extern int timer_settime(timer_t, int, const struct itimerspec *,
		struct itimerspec *);
extern int nanosleep(const struct timespec *, struct timespec *);




extern void tzset(void);

extern char *tzname[2];


extern long _sysconf(int);	 

				 




extern long timezone;
extern int daylight;





extern int cftime(char *, char *, const time_t *);
extern int ascftime(char *, const char *, const struct tm *);
extern long altzone;




extern struct tm *getdate(const char *);






extern int getdate_err;



# 222 "/usr/include/time.h" 3


 



 





























# 331 "/usr/include/time.h" 3







# 405 "/usr/include/sys/time.h" 2 3



 









# 1 "/usr/include/sys/select.h" 1 3
 
 
 

 
 
 

# 107 "/usr/include/sys/select.h" 3

# 418 "/usr/include/sys/time.h" 2 3










# 17 "/usr/include/sys/select.h" 2 3







 




















typedef	long	fd_mask;

typedef	long	fds_mask;

 























typedef	struct fd_set {



	long	fds_bits[((( 1024  )+((  (sizeof (fds_mask) * 8 )  )-1))/(  (sizeof (fds_mask) * 8 )  )) ];
} fd_set;


















extern int select(int, fd_set *, fd_set *, fd_set *, struct timeval *);










# 539 "/usr/include/sys/types.h" 2 3




 








 







# 20 "/usr/include/unistd.h" 2 3

# 1 "/usr/include/sys/unistd.h" 1 3
 









 




 








#pragma ident	"@(#)unistd.h	1.36	98/07/16 SMI"	







 



 


 




 





 


















 

 








 




 




 



























 





















 



 





 







 


















 




















 








 

 









 




 










































# 21 "/usr/include/unistd.h" 2 3






 













 















 




 





















 

























 









 
# 157 "/usr/include/unistd.h" 3


 
# 196 "/usr/include/unistd.h" 3




extern int access(const char *, int);


extern int acct(const char *);

extern unsigned alarm(unsigned);


extern int brk(void *);

extern int chdir(const char *);
extern int chown(const char *, uid_t, gid_t);


extern int chroot(const char *);

extern int close(int);







extern char *ctermid(char *);




extern char *cuserid(char *);
extern int dup(int);
extern int dup2(int, int);




extern int execl(const char *, const char *, ...);
extern int execle(const char *, const char *, ...);
extern int execlp(const char *, const char *, ...);
extern int execv(const char *, char *const *);
extern int execve(const char *, char *const *, char *const *);
extern int execvp(const char *, char *const *);
extern void _exit(int);
 






extern int fattach(int, const char *);



extern int fchdir(int);
extern int fchown(int, uid_t, gid_t);



extern int fchroot(int);



extern int fdatasync(int);

 






extern int fdetach(const char *);

extern pid_t fork(void);


extern pid_t fork1(void);

extern long fpathconf(int, int);


extern int fsync(int);



extern int ftruncate(int, off_t);

extern char *getcwd(char *, size_t);


extern int getdtablesize(void);

extern gid_t getegid(void);
extern uid_t geteuid(void);
extern gid_t getgid(void);
extern int getgroups(int, gid_t *);


extern long gethostid(void);





extern int gethostname(char *, int);

extern char *getlogin(void);









extern int getpagesize(void);
extern pid_t getpgid(pid_t);

extern pid_t getpid(void);
extern pid_t getppid(void);
extern pid_t getpgrp(void);


char *gettxt(const char *, const char *);



extern pid_t getsid(pid_t);

extern uid_t getuid(void);


extern char *getwd(char *);

 






extern int ioctl(int, int, ...);



extern int isaexec(const char *, char *const *, char *const *);

extern int isatty(int);
extern int link(const char *, const char *);


extern int lchown(const char *, uid_t, gid_t);
extern int lockf(int, int, off_t);
extern int readlink(const char *, char *, size_t);

extern off_t lseek(int, off_t, int);


extern int nice(int);

extern long pathconf(const char *, int);
extern int pause(void);
extern int pipe(int *);


extern offset_t llseek(int, offset_t, int);
extern off_t tell(int);
extern int mincore(caddr_t, size_t, char *);



extern ssize_t pread(int, void *, size_t, off_t);



extern void profil(unsigned short *, size_t, unsigned long, unsigned int);



extern int pthread_atfork(void (*) (void), void (*) (void), void (*) (void));



extern long ptrace(int, pid_t, long, long);



extern ssize_t pwrite(int, const void *, size_t, off_t);

extern ssize_t read(int, void *, size_t);


extern int rename(const char *, const char *);



extern int resolvepath(const char *, char *, size_t);

extern int rmdir(const char *);


extern void *sbrk(intptr_t);

extern int setgid(gid_t);


extern int setegid(gid_t);



extern int setgroups(int, const gid_t *);

extern int setpgid(pid_t, pid_t);


extern pid_t setpgrp(void);
extern int setregid(gid_t, gid_t);
extern int setreuid(uid_t, uid_t);

extern pid_t setsid(void);
extern int setuid(uid_t);


extern int seteuid(uid_t);

extern unsigned sleep(unsigned);


extern int stime(const time_t *);







extern int symlink(const char *, const char *);
extern void sync(void);

extern long sysconf(int);








extern pid_t tcgetpgrp(int);
extern int tcsetpgrp(int, pid_t);


extern off_t tell(int);



extern int truncate(const char *, off_t);

extern char *ttyname(int);


extern useconds_t ualarm(useconds_t, useconds_t);

extern int unlink(const char *);


extern int usleep(useconds_t);



extern pid_t vfork(void);



extern void vhangup(void);

extern ssize_t write(int, const void *, size_t);


extern void yield(void);


 




extern int ftruncate64(int, off64_t);

extern off64_t lseek64(int, off64_t, int);


extern ssize_t	pread64(int, void *, size_t, off64_t);
extern ssize_t	pwrite64(int, const void *, size_t, off64_t);
extern off64_t	tell64(int);
extern int	truncate64(const char *, off64_t);
extern int	lockf64(int, int, off64_t);



# 789 "/usr/include/unistd.h" 3


 







#pragma unknown_control_flow(vfork)



 



 





























# 911 "/usr/include/unistd.h" 3







# 37 "../include/io_lat.h" 2







 
 
 
 









 
 




 

typedef struct {
  int32type magic_number;                
  char   time_stamp[64 ];  



  int32type dims[4];                     
  int32type header_bytes;                

  int32type order;                       





} gauge_header;


 








 

 

typedef struct {
  u_int32type sum31;
  u_int32type sum29;
} gauge_check;

 
 

 












# 158 "../include/io_lat.h"

extern char *gauge_info_keyword[];


 



 
 










 


 

 






 

typedef struct {
  int32type magic_number;           
  int32type dims[4];                
  int32type header_bytes;           


  int32type order;                  



  struct {                       
    int32type n_descript;           
    char   descript[200 ];   
    int32type n_param;              
    float  param[2 ];         
  } gauge_field;
} gauge_header_1996  ;


 

typedef struct {   
  int ntoken;
  char **token;
  char **value;
} QCDheader ;


 














 

 


 
 



 

 

 

 

 
typedef struct {
  FILE *         fp;             
  gauge_header*  header;         
  char *         filename;        
  int            byterevflag;    
  int32type *       rank2rcv;       
 
  int            parallel;       

  gauge_check    check;          
} gauge_file;

 
 
 

extern   char ensemble_id[256 ];
extern   int sequence_number;

 
 

void complete_U(double *u);
int big_endian();

gauge_file *restore_ascii(char *filename);
gauge_file *save_ascii(char *filename);
gauge_file *restore_serial(char *filename);
gauge_file *save_serial(char *filename);
gauge_file *restore_parallel(char *filename);
gauge_file *save_parallel(char *filename);
gauge_file *save_checkpoint(char *filename);
gauge_file *save_serial_archive(char *filename);
gauge_file *save_parallel_archive(char *filename);
int write_gauge_info_item( FILE *fpout,  
		       char *keyword,    
		       char *fmt,        

		       char *src,        
		       int count,        
		       int stride);      


 
gauge_file *save_old_binary(char *filename, float c1, float c2);

 
 
void write_appl_gauge_info(FILE *fp);

 
 
gauge_file *save_lattice( int flag, char *filename );
gauge_file *reload_lattice( int flag, char *filename);
int ask_starting_lattice( int prompt, int *flag, char *filename );
int ask_ending_lattice( int prompt, int *flag, char *filename );
void coldlat();
void funnylat();
int get_f( int prompt, char *variable_name_string, float *value );
int get_i( int prompt, char *variable_name_string, int *value );
int get_s( int prompt, char *variable_name_string, char *value );
int get_prompt( int *value );


 
 


FILE *g_open(const char *filename, const char *mode);
int g_seek(FILE *stream, off_t offset, int whence);
size_t g_write(const void *ptr, size_t size, size_t nmemb,FILE *stream);
size_t g_read(void *ptr, size_t size, size_t nmemb, FILE *stream);
int g_close(FILE *stream);

 
 

void byterevn(int32type w[], int n);
void swrite_data(FILE* fp, void *src, size_t size, 
		 char *myname, char *descrip);
void pwrite_data(FILE* fp, void *src, size_t size, 
		 char *myname, char *descrip);
void pswrite_data(int parallel, FILE* fp, void *src, size_t size, 
		  char *myname, char *descrip);
int sread_data(FILE* fp, void *src, size_t size, 
	       char *myname, char *descrip);
int pread_data(FILE* fp, void *src, size_t size, 
	       char *myname, char *descrip);
int pread_byteorder(int byterevflag, FILE* fp, void *src, size_t size, 
		    char *myname, char *descrip);
int sread_byteorder(int byterevflag, FILE* fp, void *src, size_t size, 
		    char *myname, char *descrip);
int psread_data(int parallel, FILE* fp, void *src, size_t size, 
		char *myname, char *descrip);
int psread_byteorder(int byterevflag, int parallel, FILE* fp, 
		      void *src, size_t size, char *myname, char *descrip);


# 11 "lattice.h" 2


 



 
struct site {
     
	 
	short x,y,z,t;
	 
	char parity;
	 
	int index;

     
	 
	su3_matrix link[4];

   
        su3_matrix diag,staple,tempmat1;
};
typedef struct site site;


 

 







 
extern 	int nx,ny,nz,nt;	 
extern   int volume;			 
extern 	int no_smear_level,smear_num[5],off_axis_flag;
extern   int tot_smear;   
extern 	float smear_fac;
extern 	char startfile[256 ],savefile[256 ];
extern 	int startflag;	 
extern 	int saveflag;	 
extern 	int total_iters;

 
 
extern 	int sites_on_node;		 
extern 	int even_sites_on_node;	 
extern 	int odd_sites_on_node;	 
extern 	int number_of_nodes;	 
extern   int this_node;		 

extern  gauge_file *startlat_p;

 

extern  double_prn node_prn ;


 

extern  site *lattice;

 
 

extern  char ** gen_pt[8 ];


# 14 "hvy_qpot_includes.h" 2

# 1 "../include/comdefs.h" 1


 






 

 












# 1 "../include/../include/dirs.h" 1



 
 

 















# 25 "../include/comdefs.h" 2















 











 
typedef struct {
	 
    int msg_node;
	 
    int msg_size;
	 
    char *msg_buf;
	 



    int msg_id;

} msg_tag;

 
 

void start_handlers();
void initialize_machine(int argc, char **argv);
void make_nn_gathers();
void sort_eight_neighborlists(int index);
void sort_eight_special(void **pt);

void neighbor_coords_special(
 int x,int y,int z,int t,  
 int *dirpt,               
 int fb,                   
 int *x2p,int *y2p,int *z2p,int *t2p);
                           
int make_gather(
 void (*func)(int, int, int, int, int *, int, int *, int *, int *, int *),
         		 
 int *args,		 
 int inverse,		 
 int want_even_odd,	 
 int parity_conserve);	 

void neighbor_coords(
 int x, int y, int z, int t,   
 int dir,	               
 int *x2p, int *y2p, int *z2p, int *t2p);
                              
msg_tag * start_gather(
 
 field_offset field,	 
 int size,		 
 int index,		 

 int parity,		 

 char ** dest);		 

void restart_gather(
 
 field_offset field,	 
 int size,		 
 int index,		 

 int parity,		 

 char ** dest,		 
 msg_tag *mbuf);         

msg_tag * start_gather_from_temp(
 
 void * field,		 
 int size,		 
 int index,		 

 int parity,		 

 char ** dest);		 

void restart_gather_from_temp(
 
 void * field,		 
 int size,		 
 int index,		 

 int parity,		 

 char ** dest,		 
 msg_tag *mbuf);         

void wait_gather(msg_tag *mbuf);
void cleanup_gather(msg_tag *mbuf);

msg_tag * start_general_gather(
 
 field_offset field,	 
 int size,		 
 int *displacement,	 
 int parity,		 

 char ** dest);		 

msg_tag * start_general_gather_from_temp(
 
 void * field,	         
 int size,		 
 int *displacement,	 
 int parity,		 

 char ** dest);		 


void wait_general_gather(msg_tag *mbuf);
void cleanup_general_gather(msg_tag *mbuf);

char * field_pointer_at_coordinates(
 
 int field,	 
 int size,	 
 int x,int y,int z,int t);	 

char * field_pointer_at_direction(
 
 field_offset field,	 
 int size,	 
 site *s,	 
 int direction);	 


void cleanup_field_pointer(char * buf);
void send_field(char *buf, int size, int tonode);
void get_field(char *buf, int size, int fromnode);
char * machine_type();
int mynode();
int numnodes();
void g_sync();
void g_floatsum( float *fpt);
void g_vecfloatsum( float *fpt, int nfloats);
void g_doublesum( double *dpt);
void g_vecdoublesum( double *dpt, int ndoubles);
void g_complexsum( complex *cpt);
void g_dcomplexsum( double_complex *cpt);
void g_veccomplexsum( complex *cpt, int ncomplex);
void g_vecdcomplexsum( double_complex *cpt, int ncomplex);
void g_wvectorsum( wilson_vector *wvpt);
void g_xor32( u_int32type *pt );
void g_floatmax( float *fpt);
void g_doublemax( double *dpt);
void broadcast_float(float *fpt);
void broadcast_double(double *dpt);
void broadcast_complex(complex *cpt);
void broadcast_dcomplex(double_complex *cpt);
void broadcast_bytes(char *buf,int size);
void send_integer(int tonode, int *address);
void receive_integer(int fromnode, int *address);

 

 
double dclock();
void time_stamp(char *msg);

void terminate(int status);

void normal_exit(int status);


# 15 "hvy_qpot_includes.h" 2


# 1 "../include/generic.h" 1


 







 











 
void ape_smear_dir(
  field_offset src,        

  int dir1,                
  field_offset dest,       


  float staple_weight,     
  float link_u0,           

  int space_only,          



  int nhits,               

  float tol                


 
  );

void ape_smear(
  field_offset src,        

  field_offset dest,       

  float staple_weight,     
  float link_u0,           

  int space_only,          



  int nhits,               

  float tol                


 
  );

 
void ax_gauge();

 
int32type bsd_sum (char *data,int32type total_bytes);

 
float check_unitarity( void );

 
void d_plaquette(double *ss_plaq,double *st_plaq);

 
void make_field_strength(
  field_offset link_src,        

  field_offset field_dest       

  );

 
void gaugefix(int gauge_dir,float relax_boost,int max_gauge_iter,
	      float gauge_fix_tol, field_offset diffmat, field_offset sumvec,
	      int nvector, field_offset vector_offset[], int vector_parity[],
	      int nantiherm, field_offset antiherm_offset[], 
	      int antiherm_parity[] );

 
double imp_gauge_action();
void imp_gauge_force( float eps, field_offset mom_off );
void make_loop_table();
void dsdu_qhb_subl(int dir, int subl);

 
void make_glueball_ops();
void measure_glueball_ops();

 
void hvy_pot( field_offset links );

 
void setup_layout( void );
int node_number(int x,int y,int z,int t);
int node_index(int x,int y,int z,int t);
int num_sites(int node);

 
void make_lattice();

 
void make_global_fields();

 
void path_product( const int *dir, const int length);
void path_prod_subl(const int *dir, const int length, const int subl);

 
void plaquette(float *ss_plaq,float *st_plaq);

 
complex ploop( void );

 
complex ploop_staple(float alpha_fuzz);

 
void project_su3(
   su3_matrix *w,          

   su3_matrix *q,          
   int Nhit,               
   float tol               


 
   );

 
void rand_gauge(field_offset G);

 
void ranmom( void );

 
void initialize_prn(double_prn *prn_pt, int seed, int index);
float myrand(double_prn *prn_pt);

 
void setup_restrict_fourier( int *key, int *slice);
void restrict_fourier( 
     field_offset src,	  
     field_offset space,  
     field_offset space2, 
                          
     int size,		  



     int isign);	  

 
void reunitarize( void );
int reunit_su3(su3_matrix *c);


# 17 "hvy_qpot_includes.h" 2


void ax_gauge();
int setup();
int readin(int prompt);
void gball_simp(int tot_smear);
void smearing(void);
void w_loop1(int tot_smear);
void w_loop2(int tot_smear);
void hybrid_loop1(int tot_smear);

# 61 "hybrid_loop1.c" 2







 





void hybrid_loop1(int tot_smear) {

  char myname[] = "hybrid_loop1";
  register int i,j,dir,r,t;
  int dir1,dir2,trans_path1,trans_path2;
  int disp[4];
  int nth,nxh;
  register site *s;
  su3_matrix tmat1,tmat2;
  su3_matrix *tmatp;
  float *wils_loop1,ftmp;
  su3_matrix *s_link, *flux_links_f;
  su3_matrix *trans_links, *trans_links_f; 
  su3_matrix *shape_links;
  
   
   
  enum{ M_S_LINK, M_F_LINKS_F, M_T_LINKS_F, M_STAP_POS1, 
	  M_STAP_NEG1, M_STAP_POS2, M_STAP_NEG2, NMSGS };
  
  msg_tag *mtag[NMSGS];

   

  
  enum{ STAP_POS1, STAP_NEG1, STAP_POS2, STAP_NEG2, NSHAPE };

   

  enum{ XX, YY, ZZ, NTRANS };

   
   
  enum{ TRANS_PATH1_F, TRANS_PATH2_F, T_LINK_F, NTRANS_F };
  
   
  
  const link_path trans_path[NTRANS] =
    {
      { "XX",   2, {2,0,0,0}, {0 , 0 , -1 , -1 } },
      { "YY",   2, {0,2,0,0}, {1 , 1 , -1 , -1 } },
      { "ZZ",   2, {0,0,2,0}, {2 , 2 , -1 , -1 } }
    };

   
   
  enum{ S_LINK_F, STAP_POS1_F, NFLUX_F };

   
   
  enum{ W_LOOP1, STAP_SIG_GP1, STAP_PI_U1, STAP_DELTA_G1, NWLOOP1 };
  
   
  if(NMSGS > 8 ){
    if(this_node == 0)fprintf((&__iob[2]) ,"%s: Aborted. gen_pt array too small.",myname);
    return;
  }

  if( nx != ny || nx != ny){
    if(this_node == 0)fprintf((&__iob[2]) ,"%s: Aborted. assumes nx=ny=nz",myname);
    return;
  }
  
   
  nth = nt/2;  nxh = nx/2;
  wils_loop1 = (float *)malloc(nth*nxh*sizeof(float)*NWLOOP1);
  if(wils_loop1 == 0 ){
    fprintf((&__iob[2]) ,"%s: CAN'T MALLOC wils_loop1\n",myname);
    fflush((&__iob[2]) );
    terminate(1);
  }
  
  for(i=0;i<NWLOOP1;i++) for(t=0;t<nth;t++) for(r=0;r<nxh;r++)
    wils_loop1[ ( i )*nxh*nth + ( t )*nxh +  r  ]  = 0.0;

  
   
  s_link = 
    (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
  if(s_link == 0 ){
    fprintf((&__iob[2]) ,"%s: CAN'T MALLOC s_link\n",myname);
    fflush((&__iob[2]) );
    terminate(1);
  }
  
   
  flux_links_f = 
    (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix)*NFLUX_F);
  if(flux_links_f == 0 ){
    fprintf((&__iob[2]) ,"%s: CAN'T MALLOC flux_links_f\n",myname);
    fflush((&__iob[2]) );
    terminate(1);
  }
  
   
  trans_links = 
    (su3_matrix *)malloc(NTRANS*sites_on_node*sizeof(su3_matrix));
  if(trans_links == 0 ){
    fprintf((&__iob[2]) ,"%s: CAN'T MALLOC trans_links\n",myname);
    fflush((&__iob[2]) );
    terminate(1);
  }
  
   
  trans_links_f = 
    (su3_matrix *)malloc(NTRANS_F*sites_on_node*sizeof(su3_matrix));
  if(trans_links_f == 0 ){
    fprintf((&__iob[2]) ,"%s: CAN'T MALLOC trans_links_f\n",myname);
    fflush((&__iob[2]) );
    terminate(1);
  }
     
   
   

  for(j = 0; j < NTRANS; j++){
    path_product(trans_path[j].dir, trans_path[j].length);
    for( i =0, s =lattice; i <sites_on_node; i ++, s ++) {
      su3mat_copy(&(s->tempmat1), trans_links + ( i )*NTRANS + ( j ) );
    }
  }
  
   
  shape_links = 
    (su3_matrix *)malloc(NSHAPE*sites_on_node*sizeof(su3_matrix));
  if(shape_links == 0 ){
    fprintf((&__iob[2]) ,"%s: CAN'T MALLOC shape_links\n",myname);
    fflush((&__iob[2]) );
    terminate(1);
  }
  
  
   
  for(dir= 0 ;dir<= 2 ;dir++){
    
       
      switch(dir){
      case 0 :
	dir1 = 1 ; dir2 = 2 ; 
	trans_path1 = YY; trans_path2 = ZZ;
	break;
      case 1 :
	dir1 = 2 ; dir2 = 0 ; 
	trans_path1 = ZZ; trans_path2 = XX;
	break;
      case 2 :
	dir1 = 0 ; dir2 = 1 ; 
	trans_path1 = XX; trans_path2 = YY;
	break;
      default:
	if(this_node == 0)fprintf((&__iob[2]) ,"%s unknown direction %d\n", myname, dir);
	break;
      }
      
       
      
       
       
      for( i =0, s =lattice; i <sites_on_node; i ++, s ++) {
	su3mat_copy( &(s->link[dir]), s_link +  i   );
	su3mat_copy(trans_links + ( i )*NTRANS + ( trans_path1 ) , 
		    trans_links + ( i )*NTRANS_F + ( TRANS_PATH1_F ) );
	su3mat_copy(trans_links + ( i )*NTRANS + ( trans_path2 ) , 
		    trans_links + ( i )*NTRANS_F + ( TRANS_PATH2_F ) );
	su3mat_copy( &(s->link[3 ]), trans_links + ( i )*NTRANS_F + ( T_LINK_F ) );
      }
      
       
       

      mtag[M_T_LINKS_F] = start_gather_from_temp( 
		       (void *)(trans_links_f), 
		       NTRANS_F*sizeof(su3_matrix),
		       dir, 0x03 , gen_pt[M_T_LINKS_F] );
      
       

      
      for(r=0;r<nxh;r++){
	
	if( r>0 ){
	   

	   

	  wait_gather( mtag[M_S_LINK]);
	  for( i =0, s =lattice; i <sites_on_node; i ++, s ++) {
	    su3mat_copy( (su3_matrix *)(gen_pt[M_S_LINK][i]), &(s->staple));
	  }
	  for( i =0, s =lattice; i <sites_on_node; i ++, s ++) {
	    mult_su3_nn( &(s->link[dir]), &(s->staple), 
			 s_link +  i   );
	  }
	}
	
	 
	 
	if( r==0 ){
	  mtag[M_S_LINK] = start_gather_from_temp( 
			   (void *)s_link, sizeof(su3_matrix),
			   dir, 0x03 , gen_pt[M_S_LINK] );
	}
	else if( r<(nxh-1) ){
	  restart_gather_from_temp( 
	   (void *)s_link, sizeof(su3_matrix),
	   dir, 0x03 , gen_pt[M_S_LINK], mtag[M_S_LINK] );
	}
	else{
	  cleanup_gather( mtag[M_S_LINK]);
	}
	
	 
	 
	wait_gather( mtag[M_T_LINKS_F]);
	for(j = 0; j < NTRANS_F; j++){ 
	  for( i =0, s =lattice; i <sites_on_node; i ++, s ++) {
	    su3mat_copy( (su3_matrix *)(gen_pt[M_T_LINKS_F][i]) + j, 
			 &(s->staple));
	  }
	  for( i =0, s =lattice; i <sites_on_node; i ++, s ++) {
	    su3mat_copy( &(s->staple), trans_links + ( i )*NTRANS_F + ( j )  );
	  }
	}
	
	 

	for( j = 0 ;  j <= 3 ;  j ++) {disp[j] = -trans_path[dir1].disp[j];}
	mtag[M_STAP_NEG1] = start_general_gather_from_temp(
			 s_link, sizeof(su3_matrix),
			 disp, 0x03 , gen_pt[M_STAP_NEG1] );

	 
	for( i =0, s =lattice; i <sites_on_node; i ++, s ++) {
	  mult_su3_nn( trans_links + ( i )*NTRANS + ( trans_path1 ) ,
		       s_link +  i  , &tmat1 );
	  mult_su3_na( &tmat1, 
		       trans_links + ( i )*NTRANS_F + ( TRANS_PATH1_F ) , 
		       shape_links + ( i ) + ( STAP_POS1 )*sites_on_node  );
	}
	
	 
	wait_general_gather(mtag[M_STAP_NEG1]);
	
	
	 

	for( j = 0 ;  j <= 3 ;  j ++) {disp[j] = -trans_path[dir2].disp[j];}
	mtag[M_STAP_NEG2] = start_general_gather_from_temp(
			 s_link, sizeof(su3_matrix),
			 disp, 0x03 , gen_pt[M_STAP_NEG2] );
	
	
	 
	for( i =0, s =lattice; i <sites_on_node; i ++, s ++) {
	  mult_su3_nn( trans_links + ( i )*NTRANS + ( trans_path2 ) , 
		       s_link +  i  , &tmat1 );
	  mult_su3_na( &tmat1, 
		       trans_links + ( i )*NTRANS_F + ( TRANS_PATH2_F ) , 
		       shape_links + ( i ) + ( STAP_POS2 )*sites_on_node  );
	}
	
	wait_general_gather(mtag[M_STAP_NEG2]);

	 
	for( j = 0 ;  j <= 3 ;  j ++) {disp[j] = trans_path[dir1].disp[j];}
	mtag[M_STAP_POS1] = start_general_gather_from_temp(
			 shape_links + ( 0 ) + ( STAP_POS1 )*sites_on_node , sizeof(su3_matrix),
			 disp, 0x03 , gen_pt[M_STAP_POS1] );
	
	 
	for( i =0, s =lattice; i <sites_on_node; i ++, s ++) {
	  mult_su3_an( trans_links + ( i )*NTRANS + ( trans_path1 ) , 
		       (su3_matrix *)gen_pt[M_STAP_NEG1][i], 
		       &tmat1 );
	  mult_su3_nn( &tmat1, 
		       trans_links + ( i )*NTRANS_F + ( TRANS_PATH1_F ) , 
		       shape_links + ( i ) + ( STAP_NEG1 )*sites_on_node  );
	}
	
	 
	for( i =0, s =lattice; i <sites_on_node; i ++, s ++) {
	  mult_su3_an( trans_links + ( i )*NTRANS + ( trans_path2 ) , 
		       (su3_matrix *)gen_pt[M_STAP_NEG2][i], 
		       &tmat1 );
	  mult_su3_nn( &tmat1, 
		       trans_links + ( i )*NTRANS_F + ( TRANS_PATH2_F ) , 
		       shape_links + ( i ) + ( STAP_NEG2 )*sites_on_node  );
	}
	
	wait_general_gather(mtag[M_STAP_POS1]);

	 
	for( j = 0 ;  j <= 3 ;  j ++) {disp[j] = trans_path[dir2].disp[j];}
	mtag[M_STAP_POS2] = start_general_gather_from_temp(
			 shape_links + ( 0 ) + ( STAP_POS2 )*sites_on_node , sizeof(su3_matrix),
			 disp, 0x03 , gen_pt[M_STAP_POS2] );
	
	wait_general_gather(mtag[M_STAP_POS2]);

	 

	 
	 
	 

	for( i =0, s =lattice; i <sites_on_node; i ++, s ++) {
	  su3mat_copy( s_link +  i  , flux_links_f + ( i )*NFLUX_F + ( S_LINK_F )  );
	  su3mat_copy( shape_links + ( i ) + ( STAP_POS1 )*sites_on_node , 
		       flux_links_f + ( i )*NFLUX_F + ( STAP_POS1_F )  );
	}
	
	 

	mtag[M_F_LINKS_F] = start_gather_from_temp( 
			   (void *)flux_links_f, 
			   sizeof(su3_matrix)*NFLUX_F,
			   3 , 0x03 , gen_pt[M_F_LINKS_F] );
	
	
	 

	for(t=0;t<nth;t++){
	  
	   
	  wait_gather( mtag[M_F_LINKS_F]);
	  for( i =0, s =lattice; i <sites_on_node; i ++, s ++) {
	    su3mat_copy( (su3_matrix *)(gen_pt[M_F_LINKS_F][i])+S_LINK_F, 
			 &(s->staple));
	    su3mat_copy( (su3_matrix *)(gen_pt[M_F_LINKS_F][i])+STAP_POS1_F, 
			 &(s->diag));
	  }
	  for( i =0, s =lattice; i <sites_on_node; i ++, s ++) {
	    su3mat_copy( &(s->staple), flux_links_f + ( i )*NFLUX_F + ( S_LINK_F )  );
	    su3mat_copy( &(s->diag), flux_links_f + ( i )*NFLUX_F + ( STAP_POS1_F )  );
	  }
	  
	   
	  if( t<(nth-1) ){
	    restart_gather_from_temp( 
		     (void *)flux_links_f, 
		     sizeof(su3_matrix)*NFLUX_F,
		     3 , 0x03 , gen_pt[M_F_LINKS_F], 
		     mtag[M_F_LINKS_F] );
	  }
	  else{
	    cleanup_gather( mtag[M_F_LINKS_F]);
	  }
	  
	   
	  
	   
	  for( i =0, s =lattice; i <sites_on_node; i ++, s ++) {
	     

	    if( ((s->t)+t+1)>=nt ){
	      mult_su3_nn( &(s->link[3 ]), 
			   flux_links_f + ( i )*NFLUX_F + ( S_LINK_F ) , &tmat1);
	      mult_su3_na( &tmat1, 
			   trans_links + ( i )*NTRANS_F + ( T_LINK_F ) , &tmat2);
	      tmatp = &tmat2;
	    }
	    else
	      tmatp = flux_links_f + ( i )*NFLUX_F + ( S_LINK_F ) ;

	    ftmp = realtrace_su3( tmatp, s_link +  i   );
	    wils_loop1[ ( W_LOOP1 )*nxh*nth + ( t )*nxh +  r  ]  += ftmp;

	    if(s->x ==0 && s->y ==0 && s->z ==0 && s->t ==0){
	      printf("WL1 %d %d %d %d  %d %d  %f %f\n",s->x,s->y,s->z,s->t,t,r,
		     wils_loop1[ ( W_LOOP1 )*nxh*nth + ( t )*nxh +  r  ] , ftmp);
	      dumpmat(flux_links_f + ( i )*NFLUX_F + ( S_LINK_F ) );
	      dumpmat(s_link +  i  );
	    }
	  }
	  
	   
	  
	   
	  
	  for( i =0, s =lattice; i <sites_on_node; i ++, s ++) {
	     

	    if( ((s->t)+t+1)>=nt ){
	      mult_su3_nn( &(s->link[3 ]), 
			   flux_links_f + ( i )*NFLUX_F + ( STAP_POS1_F ) , &tmat1);
	      mult_su3_na( &tmat1, 
			   trans_links + ( i )*NTRANS_F + ( T_LINK_F ) , &tmat2);
	      tmatp = &tmat2;   
	    }
	    else
	       
	      tmatp = flux_links_f + ( i )*NFLUX_F + ( STAP_POS1_F ) ;
	    
	     
	    
	    wils_loop1[ ( STAP_SIG_GP1 )*nxh*nth + ( t )*nxh +  r  ]  +=
	      realtrace_su3( tmatp, shape_links + ( i ) + ( STAP_POS1 )*sites_on_node  ) +
	      realtrace_su3( tmatp, shape_links + ( i ) + ( STAP_NEG1 )*sites_on_node  ) +
	      realtrace_su3( tmatp, shape_links + ( i ) + ( STAP_POS2 )*sites_on_node  ) +
	      realtrace_su3( tmatp, shape_links + ( i ) + ( STAP_NEG2 )*sites_on_node  );
	    
	    wils_loop1[ ( STAP_PI_U1 )*nxh*nth + ( t )*nxh +  r  ]  +=
	      realtrace_su3( tmatp, shape_links + ( i ) + ( STAP_POS1 )*sites_on_node  ) -
	      realtrace_su3( tmatp, shape_links + ( i ) + ( STAP_NEG1 )*sites_on_node  );
	    
	    wils_loop1[ ( STAP_DELTA_G1 )*nxh*nth + ( t )*nxh +  r  ]  +=
	      realtrace_su3( tmatp, shape_links + ( i ) + ( STAP_POS1 )*sites_on_node  ) -
	      realtrace_su3( tmatp, shape_links + ( i ) + ( STAP_NEG1 )*sites_on_node  ) +
	      realtrace_su3( tmatp, shape_links + ( i ) + ( STAP_POS2 )*sites_on_node  ) -
	      realtrace_su3( tmatp, shape_links + ( i ) + ( STAP_NEG2 )*sites_on_node  );

	    if(s->x ==0 && s->y ==0 && s->z ==0 && s->t ==0){
	      printf("WL1p %d %d %d %d  %d %d  %f %f\n",s->x,s->y,s->z,s->t,t,r,
		     wils_loop1[ ( W_LOOP1 )*nxh*nth + ( t )*nxh +  r  ] , ftmp);
	    }
	  }
	  
	  
	  
	}  
	
	 
	if( r<(nxh-1) ){
	  restart_gather_from_temp( 
		   (void *)trans_links_f, 
		   NTRANS_F*sizeof(su3_matrix),
		   dir, 0x03 , gen_pt[M_T_LINKS_F], 
		   mtag[M_T_LINKS_F] );
	}
	else{
	  cleanup_gather( mtag[M_T_LINKS_F]);
	}
	
      }  
      
  }  
  
   
  g_vecfloatsum(wils_loop1,nxh*nth*NWLOOP1);
  
  for(t=0;t<nth;t++) for(r=0;r<nxh;r++)
    if(this_node==0)printf ("W_LOOP1_%d  %d  %d  %e\n", tot_smear, r, t, 
		 (double)wils_loop1[ ( W_LOOP1 )*nxh*nth + ( t )*nxh +  r  ]  / (double)(9*volume) );
  
  for(t=0;t<nth;t++) for(r=0;r<nxh;r++)
    if(this_node==0)printf ("STAP_SIG_GP1_%d  %d  %d  %e\n", tot_smear, r, t, 
		 (double)wils_loop1[ ( STAP_SIG_GP1 )*nxh*nth + ( t )*nxh +  r  ]  / (double)(9*volume) );
    
  for(t=0;t<nth;t++) for(r=0;r<nxh;r++)
    if(this_node==0)printf ("STAP_PI_U1_%d  %d  %d  %e\n", tot_smear, r, t, 
		 (double)wils_loop1[ ( STAP_PI_U1 )*nxh*nth + ( t )*nxh +  r  ]  / (float)(9*volume) );
    
  for(t=0;t<nth;t++) for(r=0;r<nxh;r++)
    if(this_node==0)printf ("STAP_DELTA_G1_%d  %d  %d  %e\n", tot_smear, r, t, 
		 (double)wils_loop1[ ( STAP_DELTA_G1 )*nxh*nth + ( t )*nxh +  r  ]  / (float)(9*volume) );
  
  free( trans_links ); 
  free( trans_links_f ); 
  free( wils_loop1 );
  free( s_link );
  free( flux_links_f );
  free( shape_links );
  
}  

