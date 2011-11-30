#ifndef _INFO_H
#define _INFO_H

typedef struct {
  double final_sec;
  double final_flop;
  int status;
  int count1, count2;
} info_t;

#define INFO_ZERO  {0.0,0.0,0,0,0}

#define INFO_HISQ_SVD_COUNTER(info) ((info)->count1)
#define INFO_HISQ_FORCE_FILTER_COUNTER(info) ((info)->count2) 

#endif
