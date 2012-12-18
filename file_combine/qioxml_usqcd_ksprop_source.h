/* TEMPORARY UNTIL ADDED TO QIO! */
#include <qioxml.h>
#include <qio_string.h>
#include <qio_stdint.h>

/*********************************************************************/
/* Top level wrapper for USQCD KS propagator file XML 

   tag           member           description          
   ------------------------------------------------------------
 
*/

typedef struct {
  QIO_TagCharValue usqcdkspropsourcefileinfo_tags;
} QIO_USQCDKSPropSourceFileInfoWrapper;

#define QIO_USQCD_KSPROPSOURCEFILE_INFO_WRAPPER {\
  {"usqcdKSPropSourceFile", "", "" , 0}       \
}

/*******************************************************************/
/* Contents of USQCD kspropsourcefile XML

   tag           member           description          
   ------------------------------------------------------------
   version       version          file XML version
   type          type             file format type string
   info          info             collaboration discretion
*/


typedef struct {
  QIO_TagCharValue version ;
  QIO_TagCharValue type;
  QIO_TagCharValue info;
} QIO_USQCDKSPropSourceFileInfo;


#define QIO_USQCDKSPROPSOURCEFILEFORMATVERSION "1.0"

#define QIO_USQCD_KSPROPSOURCEFILE_INFO_TEMPLATE {\
   {"version", "", "", 0}, \
   {"type"   , "", "", 0}, \
   {"info"   , "", "", 0}  \
}

#define QIO_USQCDKSPROPSOURCEFILETYPESTRING  \
 "USQCD_ColorVectorSource"

#define QIO_USQCDKSPROPSOURCEFILETYPE 0

QIO_USQCDKSPropSourceFileInfo *QIO_create_usqcd_kspropsourcefile_info(char *info);
