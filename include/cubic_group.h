/* Group elements, conjugacy classes and character table for O_h */

/* List of class names */
enum cubic_class {E,C2f,C4,C3,C2e,i,C2fi,C4i,C3i,C2ei,N_CUBIC_CLASS};
enum cubic_irrep {A1p,A2p,Ep,T1p,T2p,A1m,A2m,Em,T1m,T2m,N_CUBIC_IRREP};

char cubic_class_name[N_CUBIC_CLASS][6] = 
{ "E","C2f","C4","C3","C2e","i","C2fi","C4i","C3i","C2ei" };

char cubic_irrep_name[N_CUBIC_IRREP][6] = 
{ "A1+","A2+","E+","T1+","T2+","A1-","A2-","E-","T1-","T2-" };

/* Character table for the classes and irreps as enumerated above */
int cubic_char[N_CUBIC_IRREP][N_CUBIC_CLASS] = 
{
  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
  { 1, 1,-1, 1,-1, 1, 1,-1, 1,-1},
  { 2, 2, 0,-1, 0, 2, 2, 0,-1, 0},
  { 3,-1, 1, 0,-1, 3,-1, 1, 0,-1},
  { 3,-1,-1, 0, 1, 3,-1,-1, 0, 1},
  { 1, 1, 1, 1, 1,-1,-1,-1,-1,-1},
  { 1, 1,-1, 1,-1,-1,-1, 1,-1, 1},
  { 2, 2, 0,-1, 0,-2,-2, 0, 1, 0},
  { 3,-1, 1, 0,-1,-3, 1,-1, 0, 1},
  { 3,-1,-1, 0, 1,-3, 1, 1, 0,-1},
};
  
/* Group elements, expressed as permutations and reflections of x,y,z */

typedef struct {
  int class;
  int p[3];
} cubic_group;

#define N_CUBIC_GROUP 48

cubic_group g[N_CUBIC_GROUP] =
{
  {E   , XUP  , YUP  , ZUP  },
  {C2f , XDOWN, YDOWN, ZUP  },
  {C2f , XDOWN, YUP  , ZDOWN},
  {C2f , XUP  , YDOWN, ZDOWN},
  {C4  , YDOWN, XUP  , ZUP  },
  {C4  , YUP  , XDOWN, ZUP  },
  {C4  , ZUP  , YUP  , XDOWN},
  {C4  , ZDOWN, YUP  , XUP  },
  {C4  , XUP  , ZDOWN, YUP  },
  {C4  , XUP  , ZUP  , YDOWN},
  {C3  , YUP  , ZUP  , XUP  },
  {C3  , ZUP  , XUP  , YUP  },
  {C3  , YUP  , ZDOWN, XDOWN},
  {C3  , ZDOWN, XDOWN, YUP  },
  {C3  , YDOWN, ZUP  , XDOWN},
  {C3  , ZUP  , XDOWN, YDOWN},
  {C3  , YDOWN, ZDOWN, XUP  },
  {C3  , ZDOWN, XUP  , YDOWN},
  {C2e , YUP  , XUP  , ZDOWN},
  {C2e , ZUP  , YDOWN, XUP  },
  {C2e , XDOWN, ZUP  , YUP  },
  {C2e , YDOWN, XDOWN, ZDOWN},
  {C2e , ZDOWN, YDOWN, XDOWN},
  {C2e , XDOWN, ZDOWN, YDOWN},
  {i   , XDOWN, YDOWN, ZDOWN},
  {C2fi, XUP  , YUP  , ZDOWN},
  {C2fi, XUP  , YDOWN, ZUP  },
  {C2fi, XDOWN, YUP  , ZUP  },
  {C4i , YUP  , XDOWN, ZDOWN},
  {C4i , YDOWN, XUP  , ZDOWN},
  {C4i , ZDOWN, YDOWN, XUP  },
  {C4i , ZUP  , YDOWN, XDOWN},
  {C4i , XDOWN, ZUP  , YDOWN},
  {C4i , XDOWN, ZDOWN, YUP  },
  {C3i , YDOWN, ZDOWN, XDOWN},
  {C3i , ZDOWN, XDOWN, YDOWN},
  {C3i , YDOWN, ZUP  , XUP  },
  {C3i , ZUP  , XUP  , YDOWN},
  {C3i , YUP  , ZDOWN, XUP  },
  {C3i , ZDOWN, XUP  , YUP  },
  {C3i , YUP  , ZUP  , XDOWN},
  {C3i , ZUP  , XDOWN, YUP  },
  {C2ei, YDOWN, XDOWN, ZUP  },
  {C2ei, ZDOWN, YUP  , XDOWN},
  {C2ei, XUP  , ZDOWN, YDOWN},
  {C2ei, YUP  , XUP  , ZUP  },
  {C2ei, ZUP  , YUP  , XUP  },
  {C2ei, XUP  , ZUP  , YUP  }
};


