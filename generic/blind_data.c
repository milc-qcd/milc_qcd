#ifdef BLIND
#include <openssl/pem.h>
#include <openssl/ssl.h>
#include <openssl/rsa.h>
#include <openssl/evp.h>
#include <openssl/bio.h>
#include <openssl/err.h>
#include <stdio.h>
//#include <complex.h>
#include "../include/complex.h"

#define static_cast(ty,ob) ((ty)(ob))

//========================================

static const int RSAkeylen = 512;

static const int RSApadding = RSA_PKCS1_PADDING;

static unsigned char pub_key[] =
  "-----BEGIN PUBLIC KEY-----\n"\
  "MFwwDQYJKoZIhvcNAQEBBQADSwAwSAJBAJN2UkG1Br8f9/2213mLDUEXvzdNWttT\n"\
  "Jr3Wyp7zMdXPuvs+RMuD7scRiO1d9CNnaSSWJY8/uem8flAB/QfhensCAwEAAQ==\n"\
  "-----END PUBLIC KEY-----\n";
// -------
// 2019-12-20 production version
static unsigned char blind_crypt[] = {
 0x1a, 0x67, 0x84, 0x28, 0x51, 0x9b, 0x13, 0x49,
 0x5b, 0x24, 0x17, 0x5f, 0x05, 0x18, 0x03, 0xa6,
 0x5a, 0x11, 0x47, 0x7d, 0xef, 0x95, 0x0e, 0xea,
 0x57, 0x2a, 0xd3, 0x95, 0x76, 0x83, 0xda, 0x54,
 0x5b, 0xea, 0xf9, 0x4e, 0x8a, 0x8e, 0x58, 0x9b,
 0xfa, 0x26, 0x02, 0x59, 0xd1, 0x50, 0x7c, 0xbd,
 0x5a, 0xe6, 0x7d, 0xb0, 0x93, 0xda, 0xac, 0xc7,
 0x71, 0x05, 0x46, 0xd0, 0x67, 0x2a, 0x34, 0xdb
};
const char* blind_short_tag = "072e62e0";

/*
// 2016 production version!!
static unsigned char blind_crypt[] = {
 0x5d, 0x54, 0x05, 0xf5, 0x67, 0x59, 0x87, 0xf5,
 0x08, 0x56, 0x31, 0xaa, 0xe4, 0xc0, 0x79, 0xf1,
 0x7e, 0x81, 0x63, 0xc8, 0x89, 0x61, 0xbe, 0x70,
 0x66, 0xb1, 0xf8, 0xe8, 0x0b, 0xad, 0x7e, 0x4a,
 0xf1, 0x8b, 0x62, 0x4c, 0xb4, 0x8d, 0x5e, 0x48,
 0xaa, 0xf9, 0x1b, 0x35, 0x0b, 0x1c, 0xa0, 0x9d,
 0x63, 0x58, 0xa7, 0x8a, 0x7e, 0x83, 0xb7, 0xcd,
 0x32, 0xf9, 0x1f, 0x7e, 0x70, 0x57, 0xac, 0xd7
};
const char* blind_short_tag = "e02dbfbc";
*/
/* test version!
static unsigned char blind_crypt[] = {
 0x60, 0x3f, 0x54, 0x46, 0x64, 0x70, 0xa8, 0x6e,
 0xde, 0x5b, 0x63, 0xc0, 0x80, 0xe5, 0x0d, 0x17,
 0x55, 0x47, 0xdb, 0x22, 0x2e, 0xe6, 0xda, 0x42,
 0x0b, 0x61, 0x79, 0x93, 0xe5, 0xa7, 0x76, 0x82,
 0x58, 0x7a, 0xab, 0x7b, 0x2f, 0xca, 0x5d, 0xb7,
 0xe1, 0x97, 0x4b, 0x2d, 0x09, 0xd5, 0xe3, 0x95,
 0xfb, 0x5e, 0xe7, 0x1c, 0x5a, 0x3d, 0x36, 0x3d,
 0xa9, 0x41, 0x0c, 0x2b, 0xed, 0x85, 0xc5, 0xbc
};
const char* blind_short_tag = "f763ec00";
*/

//========================================

static RSA* createRSA(unsigned char* key, int public)
{
  RSA *rsa= NULL;
  BIO *keybio ;
  keybio = BIO_new_mem_buf(key, -1);
  if (keybio==NULL)
    {
      printf( "Failed to create key BIO");
      return 0;
    }
  if(public)
    {
      rsa = PEM_read_bio_RSA_PUBKEY(keybio, &rsa,NULL, NULL);
    }
  else
    {
      rsa = PEM_read_bio_RSAPrivateKey(keybio, &rsa,NULL, NULL);
    }
  if(rsa == NULL)
    {
      printf( "Failed to create RSA");
    }
     return rsa;
}

static int public_decrypt(unsigned char * enc_data,int data_len,unsigned char * key, int padding, unsigned char *decrypted)
{
    RSA * rsa = createRSA(key,1);
    int  result = RSA_public_decrypt(data_len,enc_data,decrypted,rsa,padding);
    return result;
}

static void printLastError(char *msg)
{
  char * err = malloc(130);
  ERR_load_crypto_strings();
  ERR_error_string(ERR_get_error(), err);
  printf("%s ERROR: %s\n",msg, err);
  free(err);
}

static double get_blinding_factor()
{
  unsigned char decrypted[RSAkeylen/8];
  int decrypted_length = public_decrypt(blind_crypt,sizeof(blind_crypt),pub_key,RSApadding,decrypted);
  if(decrypted_length == -1)
    {
      printLastError("decrypt failed");
      exit(0);
    }
  double blind_fact;
  sscanf(static_cast(char*,decrypted),"%lf",&blind_fact);
  return blind_fact;
}

static int do_init = 1;
static double blind_factor;

void blind_vfloat(float v[],int veclen)
{
  int j;
  if(do_init)
    {
      blind_factor = get_blinding_factor();
      do_init = 0;
    }
  for(j=0; j<veclen;++j)
    v[j] *= blind_factor;
}

void blind_vdouble(double v[],int veclen)
{
  int j;
  if(do_init)
    {
      blind_factor = get_blinding_factor();
      do_init = 0;
    }
  for(j=0; j<veclen;++j)
    v[j] *= blind_factor;
}

//void blind_vfcomplex(float complex v[],int veclen) // original
void blind_vfcomplex(complex v[],int veclen) // compiles
{
  int j;
  if(do_init)
    {
      blind_factor = get_blinding_factor();
      do_init = 0;
    }
  for(j=0; j<veclen;++j)
    CMULREAL(v[j],blind_factor,v[j]);
    //v[j] *= blind_factor;
}

//void blind_vdcomplex(double complex v[],int veclen) // original
//void blind_vdcomplex(double_complex v[],int veclen)
void blind_vdcomplex(dcomplex v[],int veclen)
{
  int j;
  if(do_init)
    {
      blind_factor = get_blinding_factor();
      do_init = 0;
    }
  for(j=0; j<veclen;++j)
    CMULREAL(v[j],blind_factor,v[j]);
    //v[j] *= blind_factor;
}
#endif /* BLIND */
