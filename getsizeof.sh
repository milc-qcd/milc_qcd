#! /bin/sh

# We need to define an explicit unsigned 32 bit integer type
# Run this script to find out what base C type to use
# Put the result in include/config.h
# C. DeTar 7/26/01

cat > /tmp/ac_getsizeof.c <<EOF
#include <stdio.h>
int main(){
    printf("short is %d bits\n",8*sizeof(short));
    printf("int is %d bits\n",8*sizeof(int));
    printf("long int is %d bits\n",8*sizeof(long int));
    printf("long long is %d bits\n",8*sizeof(long long));
    printf("float is %d bits\n",8*sizeof(float));
    printf("double is %d bits\n",8*sizeof(double));
    printf("int * is %d bits\n",8*sizeof(int *));
    printf("double * is %d bits\n",8*sizeof(double *));
}
EOF
cc /tmp/ac_getsizeof.c -o /tmp/ac_getsizeof
/tmp/ac_getsizeof
/bin/rm -f /tmp/ac_getsizeof.c /tmp/ac_getsizeof



