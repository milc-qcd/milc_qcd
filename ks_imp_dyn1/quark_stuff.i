# 1 "../generic_ks/quark_stuff.c"
 
 
 


















 
 

# 1 "../generic_ks/generic_ks_includes.h" 1
 






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









# 8 "../generic_ks/generic_ks_includes.h" 2

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







# 9 "../generic_ks/generic_ks_includes.h" 2

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






# 10 "../generic_ks/generic_ks_includes.h" 2

# 1 "../generic_ks/../include/random.h" 1



 

typedef struct {
   
  unsigned long r0,r1,r2,r3,r4,r5,r6;
  unsigned long multiplier,addend,ic_state;
  float scale;
} double_prn;

 

float myrand(double_prn *prn_pt);


# 11 "../generic_ks/generic_ks_includes.h" 2

# 1 "lattice.h" 1


 

 













# 1 "defines.h" 1



 




 








 





# 19 "lattice.h" 2

# 1 "../include/generic_quark_types.h" 1



# 1 "../include/../include/macros.h" 1





 

 
 



 
 





 
 



 
 





 
 

 
 





typedef int field_offset;




 
 










 
 


 



 
 

 
 

 


 

 

 





















 
# 125 "../include/../include/macros.h"



# 137 "../include/../include/macros.h"








# 4 "../include/generic_quark_types.h" 2


 
typedef struct {
  int min;             
  int max;             
  int nrestart;        
  int parity;          
  int start_flag;      
  int nsrc;            
  float resid;         

  float size_r;        
  int converged;       
  field_offset wv1;    
  field_offset wv2;    
  field_offset wv3;    
  field_offset wv4;    
                       
} quark_invert_control;

 

 
 
typedef struct {
  float Kappa;         
  float Clov_c;        
  float U0;            
  field_offset work_f_mn;        
} dirac_clover_param;

 
typedef struct {
  float Kappa;         
} dirac_wilson_param;

 
typedef struct {
  float mass;
} dirac_ks_param;

 
 
 

typedef struct {
  int color;           
  int spin;            
  int type;            
  char descrp[30];     
  int wall_cutoff;     
  int parity;          
  float r0;            
  int x0,y0,z0,t0;      
  int src_pointer ;    

} wilson_quark_source;




# 20 "lattice.h" 2



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







# 38 "../include/io_lat.h" 2



 




# 1 "../include/../include/int32type.h" 1
 



 








# 1 "../include/../include/../include/config.h" 1


 

 


 
 
 

 
 
 



 
 
 

 
 
 
 

 
 
 

 
 

 


 


 


 
   
 



# 14 "../include/../include/int32type.h" 2


 














typedef int int32type;
typedef unsigned int u_int32type;




# 46 "../include/io_lat.h" 2



 
 
 
 









 
 




 

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

 
 

 












# 163 "../include/io_lat.h"

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

 
 
gauge_file *save_lattice( int flag, char *filenamee, float c1, float c2 );
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


# 23 "lattice.h" 2


 

# 1 "../include/su3.h" 1



# 1 "../include/../include/complex.h" 1



 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 

 






 







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

 
								 

								 

								 


								 

								 


								 


								 

								 


								 


								 



								 


								 


								 


								 

								 

								 

								 

								 


                                                                


                                                                 




# 4 "../include/su3.h" 2



 





 
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

void byterevn(int32type w[], int n);

# 27 "lattice.h" 2



 
typedef struct {
     
	 
	short x,y,z,t;
	 
	char parity;
	 
	int index;

	 
	double_prn site_prn;
	 
	int space1;


 
 
 
	 
	su3_matrix link[4];	 

	su3_matrix longlink[4];	 
	su3_matrix fatlink[4];	 







	 
 	anti_hermitmat mom[4];

	 

 	float phase[4];

	 
 	su3_vector phi;		 
 	su3_vector resid;	 
 	su3_vector cg_p;	 
 	su3_vector xxx;		 
 	su3_vector ttt;		 
 	su3_vector g_rand;	 
	 

	



# 90 "lattice.h"


	 
	su3_vector tempvec[4];	 

	su3_vector templongvec[4];	 
        su3_vector templongv1;

	su3_matrix tempmat1,staple;
} site;

 

 







 




extern 	int nx,ny,nz,nt;	 
extern   int volume;		 
extern 	int iseed;		 
extern 	int warms,trajecs,steps,niter,propinterval,nflavors;
extern 	float epsilon;
extern   float beta,mass,u0;
extern 	float rsqmin,rsqprop;
extern 	int startflag;	 

extern 	int saveflag;	 

extern 	char startfile[256 ],savefile[256 ];
extern 	int total_iters;
extern   int phases_in;  
extern   int source_start, source_inc, n_sources;
         

 
 
extern 	int sites_on_node;		 
extern 	int even_sites_on_node;	 
extern 	int odd_sites_on_node;	 
extern 	int number_of_nodes;	 
extern   int this_node;		 

 
extern  int valid_longlinks;
extern  int valid_fatlinks;

extern  gauge_file *startlat_p;

 

extern  double_prn node_prn ;


 

extern  site *lattice;

 
 

extern  char ** gen_pt[16 ];








# 12 "../generic_ks/generic_ks_includes.h" 2




# 1 "../generic_ks/../include/comdefs.h" 1


 






 

 












# 1 "../generic_ks/../include/../include/dirs.h" 1



 
 

 















# 25 "../generic_ks/../include/comdefs.h" 2















 











 






# 84 "../generic_ks/../include/comdefs.h"

struct comlink {
	 
    struct comlink *nextcomlink;
	 
    int othernode;
	 


    int n_even_connected, n_odd_connected;
	 




	 







    int *esitelist, *ositelist;
};

typedef struct comlink comlink;


 
typedef struct {
	 
    int msg_node;
	 
    int msg_size;
	 
    char *msg_buf;
	 



    int msg_id;

} msg_tag;

 
typedef struct {
    int field;	 
    int size;	 
    int index;	 
} msg_request;



 
 

void start_handlers();
void initialize_machine(int argc, char **argv);
void make_nn_gathers();
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
void get_field(char *buf, int size);
char * machine_type();
int mynode();
int numnodes();
void g_sync();
void g_floatsum( float *fpt);
void g_doublesum( double *dpt);
void g_vecdoublesum( double *dpt, int ndoubles);
void g_complexsum( complex *cpt);
void g_dcomplexsum( double_complex *cpt);
void g_veccomplexsum( complex *cpt, int ncomplex);
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
void receive_integer(int *address);

 

 
double dclock();
void time_stamp(char *msg);

void terminate(int status);

 
 
 
 
 


# 16 "../generic_ks/generic_ks_includes.h" 2

# 1 "../generic_ks/../include/generic.h" 1


 







 











 
void ax_gauge();

 
int32type bsd_sum (char *data,int32type total_bytes);

 
float check_unitarity( void );

 
void d_plaquette(double *ss_plaq,double *st_plaq);

 
void gaugefix(int gauge_dir,float relax_boost,int max_gauge_iter,
	      float gauge_fix_tol, field_offset diffmat, field_offset sumvec,
	      int nvector, field_offset vector_offset[], int vector_parity[],
	      int nantiherm, field_offset antiherm_offset[], 
	      int antiherm_parity[] );

 
double imp_gauge_action();
void imp_gauge_force( float eps, field_offset mom_off );
void make_loop_table();
void dsdu_qhb_subl(int dir, int subl);

 
void hvy_pot( field_offset links );

 
void setup_layout( void );
int node_number(int x,int y,int z,int t);
int node_index(int x,int y,int z,int t);
int num_sites(int node);

 
void make_lattice();

 
void make_global_fields();

 
void path_product( int *dir, int length);
void path_prod_subl(int *dir, int length, int subl);

 
void plaquette(float *ss_plaq,float *st_plaq);

 
complex ploop( void );

 
complex ploop_staple(float alpha_fuzz);

 
void rand_gauge(field_offset G);

 
void ranmom( void );

 
void initialize_prn(double_prn *prn_pt, int seed, int index);
float myrand(double_prn *prn_pt);

 
void setup_restrict_fourier( int *key, int *restrict);
void restrict_fourier( 
     field_offset src,	  
     field_offset space,  
     field_offset space2, 
                          
     int size,		  



     int isign);	  

 
void reunitarize( void );
int reunit_su3(su3_matrix *c);

 
void smearing( void );


# 17 "../generic_ks/generic_ks_includes.h" 2


# 1 "../generic_ks/../include/generic_ks.h" 1


 












int congrad( int niter, float rsqmin, int parity, float *rsq );
void copy_latvec(field_offset src, field_offset dest, int parity);
void dslash( field_offset src, field_offset dest, int parity );
void dslash_special( field_offset src, field_offset dest,
    int parity, msg_tag **tag, int start );
void clear_latvec(field_offset v,int parity);

void scalar_mult_latvec(field_offset src, float scalar,
			field_offset dest, int parity);
void scalar_mult_add_latvec(field_offset src1, field_offset src2,
			    float scalar, field_offset dest, int parity);
void scalar2_mult_add_su3_vector(su3_vector *a, float s1, su3_vector *b, 
				 float s2, su3_vector *c);

void scalar2_mult_add_latvec(field_offset src1,float scalar1,
			     field_offset src2,float scalar2,
			     field_offset dest,int parity);
void checkmul();
void phaseset();
void rephase( int flag );

void prefetch_vector( su3_vector * );
void prefetch_matrix( su3_matrix * );

int ks_congrad( field_offset src, field_offset dest, float mass,
     int niter, float rsqmin, int parity, float *rsq );

void cleanup_gathers(msg_tag *tags1[], msg_tag *tags2[]);
void cleanup_dslash_temps();
void dslash_fn( field_offset src, field_offset dest, int parity );
void dslash_fn_alltemp_special(su3_vector *src, su3_vector *dest,
			       int parity, msg_tag **tag, int start );
void dslash_fn_special( field_offset src, field_offset dest,
    int parity, msg_tag **tag, int start );
void dslash_fn_on_temp( su3_vector *src, su3_vector *dest, int parity );
void dslash_fn_on_temp_special(su3_vector *src, su3_vector *dest,
			       int parity, msg_tag **tag, int start );

void dslash_eo( field_offset src, field_offset dest, int parity );
void dslash_eo_special( field_offset src, field_offset dest,
    int parity, msg_tag **tag, int start );

int congrad_ks(             
     field_offset src,        
     field_offset dest,       
     quark_invert_control *qic,  
     void *dmp                
     );

int ks_invert(  
    field_offset src,    

    field_offset dest,   

    int (*invert_func)(field_offset src, field_offset dest,
			quark_invert_control *qic,
			void *dmp),
    quark_invert_control *qic,  
    void *dmp                  
    );

int ks_multicg(	 
    field_offset src,	 
    su3_vector **psim,	 
    float *masses,	 
    int num_masses,	 
    int niter,		 
    float rsqmin,	 
    int parity,		 
    float *final_rsq_ptr	 
    );

 
void f_meas_imp( field_offset phi_off, field_offset xxx_off, float mass );

 
void sym_shift(int dir, field_offset src,field_offset dest) ;
void zeta_shift(int n, int *d, field_offset src, field_offset dest ) ;
void eta_shift(int n, int *d, field_offset src, field_offset dest ) ;


void mult_flavor_vector(int mu, field_offset src, field_offset dest ) ;
void mult_flavor_tensor(int mu, int nu, field_offset src, field_offset dest ) ;
void mult_flavor_pseudovector(int mu, field_offset src, field_offset dest ) ;
void mult_flavor_pseudoscalar(field_offset src, field_offset dest ) ;

void mult_spin_vector(int mu, field_offset src, field_offset dest ) ;
void mult_spin_tensor(int mu, int nu, field_offset src, field_offset dest ) ;
void mult_spin_pseudovector(int mu, field_offset src, field_offset dest ) ;
void mult_spin_pseudoscalar(field_offset src, field_offset dest ) ;

 
void grsource(int parity);

 
void grsource_imp( field_offset dest, float mass, int parity);

 
int mat_invert_cg( field_offset src, field_offset dest, field_offset temp,
		   float mass );
int mat_invert_uml(field_offset src, field_offset dest, field_offset temp,
		   float mass );
void check_invert( field_offset src, field_offset dest, float mass,
		   float tol);
 
int nl_spectrum( float vmass, field_offset tempvec1, field_offset tempvec2,
		 field_offset tempmat1, field_offset tempmat2);

 
void make_path_table();
void eo_fermion_force( float eps, int nflavors, field_offset x_off );
void eo_fermion_force_3f( float eps, int nflav1, field_offset x1_off,
	int nflav2, field_offset x2_off  );
void load_longlinks();
void load_fatlinks();

 
int spectrum();

 
int spectrum2( float vmass, field_offset temp1, field_offset temp2 );

 
int spectrum_mom( float qmass, float amass, field_offset temp, float tol);

 
int spectrum_nd( float mass1, float mass2, float tol );

 
int spectrum_nlpi2( float qmass, float amass, field_offset temp, float tol);


# 19 "../generic_ks/generic_ks_includes.h" 2




# 25 "../generic_ks/quark_stuff.c" 2






     



# 1 "quark_action.h" 1





     




     










    static int path_ind[6 ][7 ] = {
    { 0 , -1 , -1 , -1 , -1 , -1  },	 
    { 0 , 0 , 0 , -1 , -1 , -1  },	 
    { 1 , 0 , 6 , -1 , -1 , -1  },	 
    { 1 , 2 , 0 , 5 , 6 , -1  },	 
    { 1 , 2 , 3 , 0 , 4 , 5 , 6 },	 
    { 1 , 1 , 0 , 6 , 6 , -1  },	 
    };
    static int path_length_in[6 ] = {1,3,3,5,7,5};
    static int quark_action_npaths = 6  ;
    static float path_coeff[6 ] = {
       ( 1.0/8.0)+(6.0/16.0)+(1.0/8.0),         
	     
       (-1.0/24.0),	             
       (-1.0/8.0)*0.5,	             
       ( 1.0/8.0)*0.25*0.5,          
       (-1.0/8.0)*0.125*(1.0/6.0),   
       (-1.0/16 ),                   
    };
    static char quark_action_description[] =
	"\"O(a^2): couplings(pi)=0, Naik term, No O(a^2) errors, tadpole weights\"";

# 35 "../generic_ks/quark_stuff.c" 2

     




void printpath( int *path, int length );
void path_transport( field_offset src, field_offset dest, int parity,
    int *dir, int length );
void path_transport_hwv( field_offset src, field_offset dest, int parity,
    int *dir, int length );

void compute_gen_staple(field_offset staple, int mu, int nu,
                        field_offset link, float coef ) ;





void u_shift_fermion(su3_vector *src, su3_vector *dest, int dir ) ;
void add_force_to_mom(su3_vector *back, su3_vector *forw, int dir, float coef);
void side_link_force(int mu, int nu, float coeff, su3_vector *Path,
		     su3_vector *Path_nu, su3_vector *Path_mu, 
		     su3_vector *Path_numu) ;

void u_shift_hw_fermion(half_wilson_vector *src, 
			half_wilson_vector *dest, int dir ) ;
void add_3f_force_to_mom(half_wilson_vector *back,
			 half_wilson_vector *forw, int dir, float coeff[2]) ;
void side_link_3f_force(int mu, int nu, float coeff[2], 
			half_wilson_vector *Path   , 
			half_wilson_vector *Path_nu, 
			half_wilson_vector *Path_mu, 
			half_wilson_vector *Path_numu) ;


int path_num[6 ];	 

static float act_path_coeff[6 ];  







 

struct {
    int dir[7 ];	 
    int length;		 
    float coeff;	 
    float forwback;	 
}   q_paths[688 ];
int num_q_paths;	 
int num_basic_paths;	 

int is_path_equal( int *path1, int* path2, int length );
int add_basic_path( int *vec, int length, float coeff );

 
void make_path_table() {

    int i,j;

    int k;


     
     
     

    if(this_node==0)printf ("%s\n",quark_action_description);
    num_q_paths = 0;
    num_basic_paths = 0;


     
    if(this_node==0)printf ("path coefficients: npath  path_coeff  multiplicity\n");
    for(j=0;j<quark_action_npaths;j++) {
	float this_coeff;
	this_coeff = path_coeff[j];

	for(k=1;k< path_length_in[j];k++)this_coeff /= u0;

	act_path_coeff[j] = this_coeff ;
	i = add_basic_path( path_ind[j], path_length_in[j],
	    this_coeff );
	if(this_node==0)printf ("                    %d      %e     %d\n",
	    j,this_coeff,i);
    }
}

 

int add_basic_path( int *basic_vec, int length, float coeff ) {

    int perm[8],pp[8],ir[4];
    int j,path_num;
    int vec[7 ];
    int flag;

    path_num = 0;
     



       
      for(perm[0]=0;perm[0]<4;perm[0]++)
      for(perm[1]=0;perm[1]<4;perm[1]++)
      for(perm[2]=0;perm[2]<4;perm[2]++)
      for(perm[3]=0;perm[3]<4;perm[3]++){
	if(perm[0] != perm[1] && perm[0] != perm[2] 
	  && perm[0] != perm[3] && perm[1] != perm[2]
	  && perm[1] != perm[3] && perm[2] != perm[3] ) {
	   
	  for(ir[0]=0;ir[0]<2;ir[0]++)
	  for(ir[1]=0;ir[1]<2;ir[1]++)
	  for(ir[2]=0;ir[2]<2;ir[2]++)
	  for(ir[3]=0;ir[3]<2;ir[3]++){
	    for(j=0;j<4;j++){
	      pp[j]=perm[j];

	      if(ir[j] == 1) pp[j]= (7-( pp[j] )) ;
	      pp[(7-( j )) ]= (7-( pp[j] )) ;
	    }
	     
	    for(j=0;j<length;j++) vec[j]=pp[basic_vec[j]];
	    for(j=length;j< 7 ;j++) vec[j]= -1 ;

            flag=0;
	     
	    for(j=0;j<num_q_paths;j++){
	      flag = is_path_equal( vec, q_paths[j].dir, 7  );
	      if(flag==1)break;
	    }
	    if(flag == 0 ){
	      if(num_q_paths>= 688 ){
		if(this_node==0)printf ("OOPS: MAX_NUM too small\n");
		exit(0);
	      }
	      q_paths[num_q_paths].length=length;
	      for(j=0;j< 7 ;j++) q_paths[num_q_paths].dir[j]=vec[j];
		 
	      if(ir[0]==0){
		q_paths[num_q_paths].coeff =  coeff;
		q_paths[num_q_paths].forwback =  +1;
	      }
	      else{
		q_paths[num_q_paths].coeff = -coeff;
		q_paths[num_q_paths].forwback = -1;
	      }
	      num_q_paths++;
	      path_num++;
	       

	    }

	  }  
        }  
      }  
    num_basic_paths++;
    return(path_num);
}  

 




 

void path_transport( field_offset src, field_offset dest, int parity,
    int *dir, int length ){
    register int i;
    register site *s;
    msg_tag *mtag0;
    int j;
    su3_vector *tmp_src,*tmp_dest,*tmp_work;  
    su3_vector *tmp_pt;  
    int tmp_parity, tmp_otherparity;  

  if( length > 0 ){
    tmp_src = (su3_vector *)malloc( sites_on_node*sizeof(su3_vector) );
    tmp_dest = (su3_vector *)malloc( sites_on_node*sizeof(su3_vector) );
    tmp_work = (su3_vector *)malloc( sites_on_node*sizeof(su3_vector) );

    for( j=length-1; j>=0; j-- ){
	 
	if( j%2==0 ){
	    tmp_parity = parity;
	    switch(tmp_parity){
		case 0x02 : tmp_otherparity= 0x01 ; break;
		case 0x01 : tmp_otherparity= 0x02 ; break;
		case 0x03 : tmp_otherparity= 0x03 ; break;
	    }
	}
	else {  
	    tmp_otherparity = parity;
	    switch(tmp_otherparity){
		case 0x02 : tmp_parity= 0x01 ; break;
		case 0x01 : tmp_parity= 0x02 ; break;
		case 0x03 : tmp_parity= 0x03 ; break;
	    }
	}

	if( j==length-1 ){
	    for(  i =(( tmp_otherparity )== 0x01  ? even_sites_on_node : 0 ),  s = &(lattice[ i ]);  i < ( ( tmp_otherparity )== 0x02  ? even_sites_on_node : sites_on_node);  i ++, s ++) {
	        tmp_src[i] = *(su3_vector *)((char *)(  s  ) + ( src )) ;
	    }
	}

	if( ( dir[j] <= 3 )  ) {
	    mtag0 = start_gather_from_temp( tmp_src, sizeof(su3_vector),
	        dir[j], tmp_parity, gen_pt[0] );
	    wait_gather(mtag0);
	    for(  i =(( tmp_parity )== 0x01  ? even_sites_on_node : 0 ),  s = &(lattice[ i ]);  i < ( ( tmp_parity )== 0x02  ? even_sites_on_node : sites_on_node);  i ++, s ++) {
		mult_su3_mat_vec( &(s->link[dir[j]]),
		    (su3_vector *)(gen_pt[0][i]),
		    &(tmp_dest[i]) );
	    }
	    cleanup_gather(mtag0);
	}

	else{  
	    for(  i =(( tmp_otherparity )== 0x01  ? even_sites_on_node : 0 ),  s = &(lattice[ i ]);  i < ( ( tmp_otherparity )== 0x02  ? even_sites_on_node : sites_on_node);  i ++, s ++) {
		mult_adj_su3_mat_vec( &(s->link[(7-( dir[j] )) ]),
		    &(tmp_src[i]), &(tmp_work[i]) );
	    }
	    mtag0 = start_gather_from_temp( tmp_work, sizeof(su3_vector),
	        dir[j], tmp_parity, gen_pt[0] );
	    wait_gather(mtag0);
	    for(  i =(( tmp_parity )== 0x01  ? even_sites_on_node : 0 ),  s = &(lattice[ i ]);  i < ( ( tmp_parity )== 0x02  ? even_sites_on_node : sites_on_node);  i ++, s ++) {
		 tmp_dest[i] = *(su3_vector *)gen_pt[0][i];
	    }
	    cleanup_gather(mtag0);
	}
	
	 
	tmp_pt=tmp_src; tmp_src=tmp_dest; tmp_dest=tmp_pt;
    }   
     
    for(  i =(( parity )== 0x01  ? even_sites_on_node : 0 ),  s = &(lattice[ i ]);  i < ( ( parity )== 0x02  ? even_sites_on_node : sites_on_node);  i ++, s ++) {
        *(su3_vector *)((char *)(  s  ) + ( dest ))  = tmp_src[i];
    }
    free(tmp_src); free(tmp_dest); free(tmp_work);
  }  
  else if( src != dest ){  
    for(  i =(( parity )== 0x01  ? even_sites_on_node : 0 ),  s = &(lattice[ i ]);  i < ( ( parity )== 0x02  ? even_sites_on_node : sites_on_node);  i ++, s ++) {
        *(su3_vector *)((char *)(  s  ) + ( dest ))  = *(su3_vector *)((char *)(  s  ) + ( src )) ;
    }
  }
}  

 
void path_transport_hwv( field_offset src, field_offset dest, int parity,
    int *dir, int length ){
    register int i;
    register site *s;
    msg_tag *mtag0;
    int j;
    half_wilson_vector *tmp_src,*tmp_dest,*tmp_work;  
    half_wilson_vector *tmp_pt;  
    int tmp_parity, tmp_otherparity;  

  if( length > 0 ){
    tmp_src = (half_wilson_vector *)malloc(
      sites_on_node*sizeof(half_wilson_vector) );
    tmp_dest = (half_wilson_vector *)malloc(
       sites_on_node*sizeof(half_wilson_vector) );
    tmp_work = (half_wilson_vector *)malloc(
       sites_on_node*sizeof(half_wilson_vector) );

    for( j=length-1; j>=0; j-- ){
	 
	if( j%2==0 ){
	    tmp_parity = parity;
	    switch(tmp_parity){
		case 0x02 : tmp_otherparity= 0x01 ; break;
		case 0x01 : tmp_otherparity= 0x02 ; break;
		case 0x03 : tmp_otherparity= 0x03 ; break;
	    }
	}
	else {  
	    tmp_otherparity = parity;
	    switch(tmp_otherparity){
		case 0x02 : tmp_parity= 0x01 ; break;
		case 0x01 : tmp_parity= 0x02 ; break;
		case 0x03 : tmp_parity= 0x03 ; break;
	    }
	}

	if( j==length-1 ){
	    for(  i =(( tmp_otherparity )== 0x01  ? even_sites_on_node : 0 ),  s = &(lattice[ i ]);  i < ( ( tmp_otherparity )== 0x02  ? even_sites_on_node : sites_on_node);  i ++, s ++) {
	        tmp_src[i] = *(half_wilson_vector *)((char *)(  s  ) + ( src )) ;
	    }
	}

	if( ( dir[j] <= 3 )  ) {
	    mtag0 = start_gather_from_temp( tmp_src,
		sizeof(half_wilson_vector), dir[j], tmp_parity, gen_pt[0] );
	    wait_gather(mtag0);
	    for(  i =(( tmp_parity )== 0x01  ? even_sites_on_node : 0 ),  s = &(lattice[ i ]);  i < ( ( tmp_parity )== 0x02  ? even_sites_on_node : sites_on_node);  i ++, s ++) {
		mult_su3_mat_hwvec( &(s->link[dir[j]]),
		    (half_wilson_vector *)(gen_pt[0][i]),
		    &(tmp_dest[i]) );
	    }
	    cleanup_gather(mtag0);
	}

	else{  
	    for(  i =(( tmp_otherparity )== 0x01  ? even_sites_on_node : 0 ),  s = &(lattice[ i ]);  i < ( ( tmp_otherparity )== 0x02  ? even_sites_on_node : sites_on_node);  i ++, s ++) {
		mult_adj_su3_mat_hwvec( &(s->link[(7-( dir[j] )) ]),
		    &(tmp_src[i]), &(tmp_work[i]) );
	    }
	    mtag0 = start_gather_from_temp( tmp_work,
		sizeof(half_wilson_vector), dir[j], tmp_parity, gen_pt[0] );
	    wait_gather(mtag0);
	    for(  i =(( tmp_parity )== 0x01  ? even_sites_on_node : 0 ),  s = &(lattice[ i ]);  i < ( ( tmp_parity )== 0x02  ? even_sites_on_node : sites_on_node);  i ++, s ++) {
		 tmp_dest[i] = *(half_wilson_vector *)gen_pt[0][i];
	    }
	    cleanup_gather(mtag0);
	}
	
	 
	tmp_pt=tmp_src; tmp_src=tmp_dest; tmp_dest=tmp_pt;
    }   
     
    for(  i =(( parity )== 0x01  ? even_sites_on_node : 0 ),  s = &(lattice[ i ]);  i < ( ( parity )== 0x02  ? even_sites_on_node : sites_on_node);  i ++, s ++) {
        *(half_wilson_vector *)((char *)(  s  ) + ( dest ))  = tmp_src[i];
    }
    free(tmp_src); free(tmp_dest); free(tmp_work);
  }  
  else if( src != dest ){  
    for(  i =(( parity )== 0x01  ? even_sites_on_node : 0 ),  s = &(lattice[ i ]);  i < ( ( parity )== 0x02  ? even_sites_on_node : sites_on_node);  i ++, s ++) {
        *(half_wilson_vector *)((char *)(  s  ) + ( dest ))  =
	  *(half_wilson_vector *)((char *)(  s  ) + ( src )) ;
    }
  }
}  


# 412 "../generic_ks/quark_stuff.c"



 

 
 
 
 





 
void load_longlinks() {
  register int i;
  register site *s;
  int ipath,dir;
  int disp[4];
  register su3_matrix *long1;





  if( phases_in != 1){
    if(this_node==0)printf ("BOTCH: load_longlinks needs phases in\n");
    terminate(0);
  }
  for (dir= 0 ; dir<= 3 ; dir++){  
     
    for( i =0, s =lattice; i <sites_on_node; i ++, s ++) {



      long1 = &(s->longlink[dir]);

	clear_su3mat( long1 );
    }

     
    for( ipath=0; ipath<num_q_paths; ipath++ ){   
	 
	for(i= 0 ;i<= 3 ;i++)disp[i]=0;
	for( i=0; i<q_paths[ipath].length; i++){
	  if( ( q_paths[ipath].dir[i] <= 3 )  )
	    disp[        q_paths[ipath].dir[i]  ]++;
	  else
	    disp[(7-( q_paths[ipath].dir[i] ))  ]--;
	}
	for( disp[dir]+=3,i= 0 ; i<= 3 ; i++)if(disp[i]!=0)break;
	if( i<= 3  )continue;   
 



	path_product( q_paths[ipath].dir, q_paths[ipath].length );
	for( i =0, s =lattice; i <sites_on_node; i ++, s ++) {
	  su3_adjoint( &(s->tempmat1), &(s->staple) );



	  long1 = &(s->longlink[dir]);

          scalar_mult_add_su3_matrix( long1,
	    &(s->staple), -q_paths[ipath].coeff, long1 );
		 
	}
    }  

  }  

  valid_longlinks = 1;




}   

 
void load_fatlinks() {
  register int i;
  register site *s;
  int dir;
  register su3_matrix *fat1;

  int  nu,rho,sig ;
  float one_link ;  









  if( phases_in != 1){
    if(this_node==0)printf ("BOTCH: load_fatlinks needs phases in\n");
    terminate(0);
  }

# 557 "../generic_ks/quark_stuff.c"

 



 
  
 one_link = (act_path_coeff[0] - 6.0*act_path_coeff[5]) ; 
 
 for (dir= 0 ; dir<= 3 ; dir++){
   for( i =0, s =lattice; i <sites_on_node; i ++, s ++)   



      fat1 = &(s->fatlink[dir]);

     scalar_mult_su3_matrix(&(s->link[dir]), one_link ,
			    fat1 );
   for(nu= 0 ; nu<= 3 ; nu++) if(nu!=dir)
     {
       compute_gen_staple(((field_offset)(((char *)&(lattice[0].  staple  ))-((char *)&(lattice[0])) )) ,dir,nu,((field_offset)(((char *)&(lattice[0].  link[dir]  ))-((char *)&(lattice[0])) )) ,
			  act_path_coeff[2]);
        
        
       compute_gen_staple(-1 ,dir,nu,((field_offset)(((char *)&(lattice[0].  staple  ))-((char *)&(lattice[0])) )) ,act_path_coeff[5]);
       for(rho= 0 ; rho<= 3 ; rho++) if((rho!=dir)&&(rho!=nu))
	 {
	   compute_gen_staple(((field_offset)(((char *)&(lattice[0].  tempmat1  ))-((char *)&(lattice[0])) )) ,dir,rho,((field_offset)(((char *)&(lattice[0].  staple  ))-((char *)&(lattice[0])) )) ,
			      act_path_coeff[3]);
	   for(sig= 0 ; sig<= 3 ; sig++)
	     if((sig!=dir)&&(sig!=nu)&&(sig!=rho))
	       {
		 compute_gen_staple(-1 ,dir,sig,
				    ((field_offset)(((char *)&(lattice[0].  tempmat1  ))-((char *)&(lattice[0])) )) ,
				    act_path_coeff[4]);
	       }  
	 }  
     }  
 }   


  valid_fatlinks = 1;




}   



 
int is_path_equal( int *path1, int* path2, int length ){
   register int i;
   for(i=0;i<length;i++)if(path1[i]!=path2[i])return(0);
   return(1);
}


 
 

 

# 885 "../generic_ks/quark_stuff.c"

 















void eo_fermion_force( float eps, int nflavors, field_offset x_off ){
   
   
   
   
  register int i ;
  register site *s;
  int mu,nu,rho,sig ;
  int DirectLinks[8] ;
  float ferm_epsilon, coeff;
  float OneLink, Lepage, Naik, FiveSt, ThreeSt, SevenSt ;
  su3_vector *tempvec[8] ;
  su3_vector *temp_x ;






  ferm_epsilon = 2.0*(nflavors/4.0)*eps;
  
   
  OneLink = act_path_coeff[0]*ferm_epsilon ; 
  Naik    = act_path_coeff[1]*ferm_epsilon ;
  ThreeSt = act_path_coeff[2]*ferm_epsilon ;
  FiveSt  = act_path_coeff[3]*ferm_epsilon ;
  SevenSt = act_path_coeff[4]*ferm_epsilon ;
  Lepage  = act_path_coeff[5]*ferm_epsilon ;
   

   
  for(mu=0;mu<8;mu++)
    DirectLinks[mu] = 0 ;

   
  for(mu=0;mu<8;mu++)
    tempvec[mu] = (su3_vector *)malloc( sites_on_node*sizeof(su3_vector) );

   
  temp_x = (su3_vector *)malloc( sites_on_node*sizeof(su3_vector) );
  for( i =0, s =lattice; i <sites_on_node; i ++, s ++)  temp_x[i] = *(su3_vector *)((char *)(  s  ) + ( x_off ))  ;

  for(sig=0;sig<8;sig++)
    {
      for(mu=0;mu<8;mu++)if((mu!=sig)&&(mu!= (7-( sig )) ))
	{
	  u_shift_fermion(temp_x, tempvec[0] , (7-( mu )) );
	  u_shift_fermion(tempvec[0] , tempvec[7] , sig);
	  if(( sig <= 3 ) )
	    {
	       



	      add_force_to_mom(tempvec[7] , tempvec[0] , sig, -ThreeSt) ;
	    }
	  for(nu=0;nu<8;nu++)if((nu!=mu )&&(nu!= (7-( mu  )) )&&
				(nu!=sig)&&(nu!= (7-( sig )) ))
	    {
	      u_shift_fermion(tempvec[0] , tempvec[1] , (7-( nu )) );
	      u_shift_fermion(tempvec[1] , tempvec[6] , sig);
	      if(( sig <= 3 ) )
		{
		   



		  add_force_to_mom(tempvec[6] , tempvec[1] , sig, FiveSt);
		}
	      for(rho=0;rho<8;rho++)if((rho!=mu )&&(rho!= (7-( mu  )) )&&
				       (rho!=nu )&&(rho!= (7-( nu  )) )&&
				       (rho!=sig)&&(rho!= (7-( sig )) ))
		{
		  u_shift_fermion(tempvec[1] , tempvec[2] , (7-( rho )) );
		   
		  u_shift_fermion(tempvec[2] , tempvec[3] ,sig);
		  if(( sig <= 3 ) )
		    {
		       



		      add_force_to_mom(tempvec[3] , tempvec[2] , sig, -SevenSt ) ;
		    }
		   
		  u_shift_fermion(tempvec[3] , tempvec[4] , rho);
		  side_link_force(rho,sig,SevenSt, tempvec[1] , tempvec[3] , tempvec[2] , tempvec[4] );
		   
		  coeff = SevenSt/FiveSt ; 
		  for( i =0, s =lattice; i <sites_on_node; i ++, s ++) 
		    scalar_mult_add_su3_vector(&(tempvec[6] [i]),&(tempvec[4] [i]),coeff,
					       &(tempvec[6] [i]));
		} 
	       
	       
	      u_shift_fermion(tempvec[6] ,tempvec[3] , nu);
	      side_link_force(nu,sig,-FiveSt,tempvec[0] ,tempvec[6] , 
			      tempvec[1] ,tempvec[3] ) ;
	       
	      coeff = FiveSt/ThreeSt ; 
	      for( i =0, s =lattice; i <sites_on_node; i ++, s ++) 
		scalar_mult_add_su3_vector(&(tempvec[7] [i]),&(tempvec[3] [i]),coeff,&(tempvec[7] [i]));
	    } 

	   

	  u_shift_fermion(tempvec[0] , tempvec[1] , (7-( mu )) );
	  u_shift_fermion(tempvec[1] , tempvec[6] , sig);
	  if(( sig <= 3 ) )
	    {
	       



	      add_force_to_mom(tempvec[6] , tempvec[1] , sig, Lepage) ;
	    }
	   
	  u_shift_fermion(tempvec[6] ,tempvec[3] , mu);
	  side_link_force(mu, sig, -Lepage, tempvec[0] , tempvec[6] , tempvec[1] , tempvec[3] ) ;
	   
	  coeff = Lepage/ThreeSt ; 
	  for( i =0, s =lattice; i <sites_on_node; i ++, s ++) 
	    scalar_mult_add_su3_vector(&(tempvec[7] [i]),&(tempvec[3] [i]),coeff,&(tempvec[7] [i]));

	   
	   
	  if(( mu <= 3 ) ) 
	    u_shift_fermion(tempvec[7] ,tempvec[3] , mu );
	   
	  side_link_force(mu, sig, ThreeSt, temp_x, tempvec[7] , tempvec[0] , tempvec[3] );

	   
	   
	  if( (!DirectLinks[mu]) ){
	    DirectLinks[mu]=1 ;
	    if(( mu > 3 ) ) 
	      {
		 

		 
		add_force_to_mom(tempvec[0] , temp_x, (7-( mu )) , OneLink) ;
		 

		 
		u_shift_fermion(temp_x, tempvec[4] , mu);
		 
		 
		add_force_to_mom(tempvec[1] , tempvec[4] , (7-( mu )) , -Naik) ;
		 
		u_shift_fermion(tempvec[1] , tempvec[4] , (7-( mu )) );
		 
		add_force_to_mom(tempvec[4] , temp_x, (7-( mu )) , Naik);
	      }
	    else  
	      {
		u_shift_fermion(temp_x, tempvec[4] , mu);
		 
		 
		add_force_to_mom(tempvec[4] , tempvec[1] , mu, Naik) ;
	      }
	  }
	} 
       
    } 

   
  free(temp_x) ;
  for(mu=0;mu<8;mu++)
    free(tempvec[mu]) ;






}  

























void eo_fermion_force_3f( float eps, int nflav1, field_offset x1_off, 
			  int nflav2, field_offset x2_off){
   
   
   
   
   
  register int i ;
  register site *s;
  int mu,nu,rho,sig ;
  int DirectLinks[8] ;
  float coeff[2],ferm_epsilon ;
  float OneLink[2], Lepage[2], Naik[2], FiveSt[2], ThreeSt[2], SevenSt[2] ;
  float mNaik[2], mLepage[2], mFiveSt[2], mThreeSt[2], mSevenSt[2] ;
  half_wilson_vector *hwvec[8] ;
  half_wilson_vector *temp_x ;






  
   
  ferm_epsilon = 2.0*(nflav1/4.0)*eps;
  OneLink[0] = act_path_coeff[0]*ferm_epsilon ;
  Naik[0]    = act_path_coeff[1]*ferm_epsilon ; mNaik[0]    = -Naik[0] ;
  ThreeSt[0] = act_path_coeff[2]*ferm_epsilon ; mThreeSt[0] = -ThreeSt[0] ;
  FiveSt[0]  = act_path_coeff[3]*ferm_epsilon ; mFiveSt[0]  = -FiveSt[0] ;
  SevenSt[0] = act_path_coeff[4]*ferm_epsilon ; mSevenSt[0] = -SevenSt[0] ;
  Lepage[0]  = act_path_coeff[5]*ferm_epsilon ; mLepage[0]  = -Lepage[0] ;

  ferm_epsilon = 2.0*(nflav2/4.0)*eps;
  OneLink[1] = act_path_coeff[0]*ferm_epsilon ;
  Naik[1]    = act_path_coeff[1]*ferm_epsilon ; mNaik[1]    = -Naik[1] ;
  ThreeSt[1] = act_path_coeff[2]*ferm_epsilon ; mThreeSt[1] = -ThreeSt[1] ;
  FiveSt[1]  = act_path_coeff[3]*ferm_epsilon ; mFiveSt[1]  = -FiveSt[1] ;
  SevenSt[1] = act_path_coeff[4]*ferm_epsilon ; mSevenSt[1] = -SevenSt[1] ;
  Lepage[1]  = act_path_coeff[5]*ferm_epsilon ; mLepage[1]  = -Lepage[1] ;
   

   
  for(mu=0;mu<8;mu++)
    DirectLinks[mu] = 0 ;

   
  for(mu=0;mu<8;mu++)
    hwvec[mu]= 
      (half_wilson_vector *)malloc(sites_on_node*sizeof(half_wilson_vector));

   
  temp_x= 
    (half_wilson_vector *)malloc(sites_on_node*sizeof(half_wilson_vector));
  for( i =0, s =lattice; i <sites_on_node; i ++, s ++) 
    {
      temp_x[i].h[0] = *(su3_vector *)((char *)(  s  ) + ( x1_off ))  ;
      temp_x[i].h[1] = *(su3_vector *)((char *)(  s  ) + ( x2_off ))  ;
    }

  for(sig=0;sig<8;sig++)
    {
      for(mu=0;mu<8;mu++)if((mu!=sig)&&(mu!= (7-( sig )) ))
	{
	  u_shift_hw_fermion(temp_x, hwvec[0] , (7-( mu )) );
	  u_shift_hw_fermion(hwvec[0] , hwvec[7] , sig);
	  if(( sig <= 3 ) )
	    {
	       



	      add_3f_force_to_mom(hwvec[7] , hwvec[0] , sig, mThreeSt) ;
	    }
	  for(nu=0;nu<8;nu++)if((nu!=mu )&&(nu!= (7-( mu  )) )&&
				(nu!=sig)&&(nu!= (7-( sig )) ))
	    {
	      u_shift_hw_fermion(hwvec[0] , hwvec[1] , (7-( nu )) );
	      u_shift_hw_fermion(hwvec[1] , hwvec[6] , sig);
	      if(( sig <= 3 ) )
		{
		   



		  add_3f_force_to_mom(hwvec[6] , hwvec[1] , sig, FiveSt);
		}
	      for(rho=0;rho<8;rho++)if((rho!=mu )&&(rho!= (7-( mu  )) )&&
				       (rho!=nu )&&(rho!= (7-( nu  )) )&&
				       (rho!=sig)&&(rho!= (7-( sig )) ))
		{
		  u_shift_hw_fermion(hwvec[1] , hwvec[2] , (7-( rho )) );
		   
		  u_shift_hw_fermion(hwvec[2] , hwvec[3] ,sig);
		  if(( sig <= 3 ) )
		    {
		       



		      add_3f_force_to_mom(hwvec[3] , hwvec[2] , sig, mSevenSt ) ;
		    }
		   
		  u_shift_hw_fermion(hwvec[3] , hwvec[4] , rho);
		  side_link_3f_force(rho,sig,SevenSt,hwvec[1] ,hwvec[3] ,hwvec[2] ,hwvec[4] );
		   
		  coeff[0] = SevenSt[0]/FiveSt[0] ; 
		  coeff[1] = SevenSt[1]/FiveSt[1] ; 
		  for( i =0, s =lattice; i <sites_on_node; i ++, s ++) 
		    {
		      scalar_mult_add_su3_vector(&(hwvec[6] [i].h[0]),
						 &(hwvec[4] [i].h[0]),coeff[0],
						 &(hwvec[6] [i].h[0]));
		      scalar_mult_add_su3_vector(&(hwvec[6] [i].h[1]),
						 &(hwvec[4] [i].h[1]),coeff[1],
						 &(hwvec[6] [i].h[1]));
		    }
		} 
	       
	       
	      u_shift_hw_fermion(hwvec[6] ,hwvec[3] , nu);
	      side_link_3f_force(nu,sig,mFiveSt,hwvec[0] ,hwvec[6] , 
			      hwvec[1] ,hwvec[3] ) ;
	       
	      coeff[0] = FiveSt[0]/ThreeSt[0] ; 
	      coeff[1] = FiveSt[1]/ThreeSt[1] ; 
	      for( i =0, s =lattice; i <sites_on_node; i ++, s ++) 
		{
		  scalar_mult_add_su3_vector(&(hwvec[7] [i].h[0]),
					     &(hwvec[3] [i].h[0]), coeff[0],
					     &(hwvec[7] [i].h[0]));
		  scalar_mult_add_su3_vector(&(hwvec[7] [i].h[1]),
					     &(hwvec[3] [i].h[1]), coeff[1],
					     &(hwvec[7] [i].h[1]));
		}
	    } 

	   

	  u_shift_hw_fermion(hwvec[0] , hwvec[1] , (7-( mu )) );
	  u_shift_hw_fermion(hwvec[1] , hwvec[6] , sig);
	  if(( sig <= 3 ) )
	    {
	       



	      add_3f_force_to_mom(hwvec[6] , hwvec[1] , sig, Lepage) ;
	    }
	   
	  u_shift_hw_fermion(hwvec[6] ,hwvec[3] , mu);
	  side_link_3f_force(mu, sig, mLepage, hwvec[0] , hwvec[6] , hwvec[1] , hwvec[3] ) ;
	   
	  coeff[0] = Lepage[0]/ThreeSt[0] ; 
	  coeff[1] = Lepage[1]/ThreeSt[1] ; 
	  for( i =0, s =lattice; i <sites_on_node; i ++, s ++) 
	    {
	      scalar_mult_add_su3_vector(&(hwvec[7] [i].h[0]),
					 &(hwvec[3] [i].h[0]),coeff[0],
					 &(hwvec[7] [i].h[0]));
	      scalar_mult_add_su3_vector(&(hwvec[7] [i].h[1]),
					 &(hwvec[3] [i].h[1]),coeff[1],
					 &(hwvec[7] [i].h[1]));
	    }

	   
	   
	  if(( mu <= 3 ) ) 
	    u_shift_hw_fermion(hwvec[7] ,hwvec[3] , mu );
	   
	  side_link_3f_force(mu, sig, ThreeSt, temp_x, hwvec[7] , hwvec[0] , hwvec[3] );

	   
	   
	  if( (!DirectLinks[mu]) ){
	    DirectLinks[mu]=1 ;
	    if(( mu > 3 ) ) 
	      {
		 

		 
		add_3f_force_to_mom(hwvec[0] , temp_x, (7-( mu )) , OneLink) ;
		 

		 
		u_shift_hw_fermion(temp_x, hwvec[4] , mu);
		 
		 
		add_3f_force_to_mom(hwvec[1] , hwvec[4] , (7-( mu )) , mNaik) ;
		 
		u_shift_hw_fermion(hwvec[1] , hwvec[4] , (7-( mu )) );
		 
		add_3f_force_to_mom(hwvec[4] , temp_x, (7-( mu )) , Naik);
	      }
	    else  
	      {
		u_shift_hw_fermion(temp_x, hwvec[4] , mu);
		 
		 
		add_3f_force_to_mom(hwvec[4] , hwvec[1] , mu, Naik) ;
	      }
	  }
	} 
       
    } 

   
  free(temp_x) ;
  for(mu=0;mu<8;mu++)
    free(hwvec[mu]) ;






}  



















void compute_gen_staple(field_offset staple, int mu, int nu, 
			field_offset link, float coef) {
  su3_matrix tmat1,tmat2;
  msg_tag *mtag0,*mtag1;
  su3_matrix *tempmat ;
  register site *s ;
  register int i ;
  register su3_matrix *fat1;

   










   
  mtag0 = start_gather( link, sizeof(su3_matrix), nu, 0x03 , gen_pt[0] );
  mtag1 = start_gather( ((field_offset)(((char *)&(lattice[0].  link[nu]  ))-((char *)&(lattice[0])) )) , sizeof(su3_matrix), mu, 
			0x03 , gen_pt[1] );
  wait_gather(mtag0);
  wait_gather(mtag1);
  
  if(staple!= -1 ){ 
    for( i =0, s =lattice; i <sites_on_node; i ++, s ++) {
      mult_su3_na( (su3_matrix *)gen_pt[0][i],
		   (su3_matrix *)gen_pt[1][i], &tmat1 );
      mult_su3_nn( &(s->link[nu]), &tmat1, (su3_matrix *)((char *)(  s  ) + ( staple ))  );
    }
  }
  else{  
    for( i =0, s =lattice; i <sites_on_node; i ++, s ++) {
      mult_su3_na( (su3_matrix *)gen_pt[0][i],
		   (su3_matrix *)gen_pt[1][i], &tmat1 );
      mult_su3_nn( &(s->link[nu]), &tmat1, &tmat2 );



      fat1 = &(s->fatlink[mu]);

      scalar_mult_add_su3_matrix(fat1, &tmat2, coef,
				 fat1) ;
    }
  }
  cleanup_gather(mtag0);
  cleanup_gather(mtag1);

   
  tempmat = (su3_matrix *)malloc( sites_on_node*sizeof(su3_matrix) );
  mtag0 = start_gather( ((field_offset)(((char *)&(lattice[0].  link[nu]  ))-((char *)&(lattice[0])) )) ,
			sizeof(su3_matrix), mu, 0x03 , gen_pt[0] );
  wait_gather(mtag0);
  for( i =0, s =lattice; i <sites_on_node; i ++, s ++) {
    mult_su3_an( &(s->link[nu]),(su3_matrix *)((char *)(  s  ) + ( link )) , &tmat1 );
    mult_su3_nn( &(tmat1),(su3_matrix *)gen_pt[0][i], &(tempmat[i]) );
  }
  cleanup_gather(mtag0);
  mtag0 = start_gather_from_temp( tempmat, sizeof(su3_matrix),
				  (7-( nu )) , 0x03 , gen_pt[0] );
  wait_gather(mtag0);

  if(staple!= -1 ){ 
    for( i =0, s =lattice; i <sites_on_node; i ++, s ++) {
      add_su3_matrix( (su3_matrix *)((char *)(  s  ) + ( staple )) ,(su3_matrix *)gen_pt[0][i], 
		      (su3_matrix *)((char *)(  s  ) + ( staple ))  );



      fat1 = &(s->fatlink[mu]);

      scalar_mult_add_su3_matrix( fat1,
				 (su3_matrix *)((char *)(  s  ) + ( staple )) , coef, 
				 fat1 );
    }
  }
  else{  
    for( i =0, s =lattice; i <sites_on_node; i ++, s ++) {



      fat1 = &(s->fatlink[mu]);

      scalar_mult_add_su3_matrix( fat1,
				 (su3_matrix *)gen_pt[0][i], coef, 
				 fat1 );
    }
  }

  free(tempmat);
  cleanup_gather(mtag0);
}




 
 
void u_shift_fermion(su3_vector *src, su3_vector *dest, int dir ) {
  su3_vector *tmpvec ; 
  msg_tag *mtag ;
  register site *s ;
  register int i ;
  
  if(( dir <= 3 ) )  
    {
      mtag = start_gather_from_temp(src, sizeof(su3_vector), 
				    dir, 0x03 , gen_pt[0]);
      wait_gather(mtag);
      for( i =0, s =lattice; i <sites_on_node; i ++, s ++) 
	mult_su3_mat_vec(&(s->link[dir]),(su3_vector *)(gen_pt[0][i]),
			 &(dest[i]));
      cleanup_gather(mtag);
    }
  else  
    {
      tmpvec = (su3_vector *)malloc( sites_on_node*sizeof(su3_vector) );
      for( i =0, s =lattice; i <sites_on_node; i ++, s ++) 
	mult_adj_su3_mat_vec(&(s->link[(7-( dir )) ]),&(src[i]), &tmpvec[i]);
      mtag = start_gather_from_temp(tmpvec, sizeof(su3_vector), dir, 
				    0x03 , gen_pt[0]);
      wait_gather(mtag);
       
      for( i =0, s =lattice; i <sites_on_node; i ++, s ++) 
	dest[i] = *(su3_vector *)gen_pt[0][i];
      cleanup_gather(mtag);
      free(tmpvec) ;
    }
}

 

void u_shift_hw_fermion(half_wilson_vector *src, 
			half_wilson_vector *dest, int dir ) {
  half_wilson_vector *tmpvec ; 
  msg_tag *mtag ;
  register site *s ;
  register int i ;
  
  if(( dir <= 3 ) )  
    {
      mtag = start_gather_from_temp(src, sizeof(half_wilson_vector), 
				    dir, 0x03 , gen_pt[0]);
      wait_gather(mtag);
      for( i =0, s =lattice; i <sites_on_node; i ++, s ++) 
	mult_su3_mat_hwvec(&(s->link[dir]),
			   (half_wilson_vector *)(gen_pt[0][i]), &(dest[i]));
      cleanup_gather(mtag);
    }
  else  
    {
      tmpvec = 
	(half_wilson_vector *)malloc(sites_on_node*sizeof(half_wilson_vector));
      for( i =0, s =lattice; i <sites_on_node; i ++, s ++) 
	mult_adj_su3_mat_hwvec(&(s->link[(7-( dir )) ]),&(src[i]), &tmpvec[i]);
      mtag = start_gather_from_temp(tmpvec, sizeof(half_wilson_vector), dir, 
				    0x03 , gen_pt[0]);
      wait_gather(mtag);
       
      for( i =0, s =lattice; i <sites_on_node; i ++, s ++) 
	dest[i] = *(half_wilson_vector *)gen_pt[0][i];
      cleanup_gather(mtag);
      free(tmpvec) ;
    }
}

 
 
void add_force_to_mom(su3_vector *back,su3_vector *forw,int dir,float coeff) {
  register site *s ;
  register int i ;  
  register float tmp_coeff ;

  su3_matrix tmat,tmat2;

  if(( dir > 3 ) )
    {
      dir = (7-( dir ))  ; 
      coeff = -coeff ;
    }
  for( i =0, s =lattice; i <sites_on_node; i ++, s ++) {
    if(s->parity== 0x01 ) 
      tmp_coeff = -coeff ;
    else
      tmp_coeff = coeff ;
    uncompress_anti_hermitian( &(s->mom[dir]), &tmat2 );
    su3_projector(&(back[i]), &(forw[i]), &tmat);
    scalar_mult_add_su3_matrix(&tmat2, &tmat,  tmp_coeff, &tmat2 );
    make_anti_hermitian( &tmat2, &(s->mom[dir]) ); 
  }
}

 
 
void add_3f_force_to_mom(half_wilson_vector *back,
			 half_wilson_vector *forw, int dir, float coeff[2]) {
  register site *s ;
  register int i ;  
  float tmp_coeff[2] ;

  su3_matrix tmat,tmat2;

  if(( dir > 3 ) )
    {
      dir = (7-( dir ))  ; 
      coeff[0] = -coeff[0] ;
      coeff[1] = -coeff[1] ;
    }
  for( i =0, s =lattice; i <sites_on_node; i ++, s ++) {
    if(s->parity== 0x01 )
      {
	tmp_coeff[0] = -coeff[0] ;
	tmp_coeff[1] = -coeff[1] ;
      }
    else
      {
	tmp_coeff[0] = coeff[0] ;
	tmp_coeff[1] = coeff[1] ;
      }
    uncompress_anti_hermitian( &(s->mom[dir]), &tmat2 );
    su3_projector(&(back[i].h[0]), &(forw[i].h[0]), &tmat);
    scalar_mult_add_su3_matrix(&tmat2, &tmat,  tmp_coeff[0], &tmat2 );
    su3_projector(&(back[i].h[1]), &(forw[i].h[1]), &tmat);
    scalar_mult_add_su3_matrix(&tmat2, &tmat,  tmp_coeff[1], &tmat2 );
    make_anti_hermitian( &tmat2, &(s->mom[dir]) ); 
  }
}

 







void side_link_force(int mu, int nu, float coeff, 
		     su3_vector *Path   , su3_vector *Path_nu, 
		     su3_vector *Path_mu, su3_vector *Path_numu) {
  if(( mu <= 3 ) )
    {
       




      if(( nu <= 3 ) )
	add_force_to_mom(Path_numu, Path, mu, coeff ) ;
      else
	add_force_to_mom(Path, Path_numu, (7-( mu )) , -coeff ); 
    }
  else  
    {
       



 
      if(( nu <= 3 ) )
	add_force_to_mom(Path_nu, Path_mu, mu, -coeff) ;  
      else
	add_force_to_mom(Path_mu, Path_nu, (7-( mu )) , coeff) ;
    }
}
 
 

void side_link_3f_force(int mu, int nu, float coeff[2], 
			half_wilson_vector *Path   , 
			half_wilson_vector *Path_nu, 
			half_wilson_vector *Path_mu, 
			half_wilson_vector *Path_numu) {
  float m_coeff[2] ;

  m_coeff[0] = -coeff[0] ;
  m_coeff[1] = -coeff[1] ;

  if(( mu <= 3 ) )
    {
       




      if(( nu <= 3 ) )
	add_3f_force_to_mom(Path_numu, Path, mu, coeff ) ;
      else
	add_3f_force_to_mom(Path,Path_numu,(7-( mu )) ,m_coeff); 
    }
  else  
    {
       



 
      if(( nu <= 3 ) )
	add_3f_force_to_mom(Path_nu, Path_mu, mu, m_coeff) ;  
      else
	add_3f_force_to_mom(Path_mu, Path_nu, (7-( mu )) , coeff) ;
    }
}


 






































































































 







































 
























