/* TEMPORARY UNTIL ADDED TO QIO! */
#include <qioxml.h>
#include <qio_string.h>
#include <qio_stdint.h>

/*********************************************************************/
/* Top level wrapper for USQCD Dirac propagator source file XML 

   tag           member           description          
   ------------------------------------------------------------
 
*/

typedef struct {
  QIO_TagCharValue usqcdpropsourcefileinfo_tags;
} QIO_USQCDPropSourceFileInfoWrapper;

#define QIO_USQCD_PROPSOURCEFILE_INFO_WRAPPER {\
  {"usqcdPropSourceFile", "", "" , 0}       \
}

/*******************************************************************/
/* Contents of USQCD propsourcefile XML

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
} QIO_USQCDPropSourceFileInfo;


#define QIO_USQCDPROPSOURCEFILEFORMATVERSION "1.0"

#define QIO_USQCD_PROPSOURCEFILE_INFO_TEMPLATE {\
   {"version", "", "", 0}, \
   {"type"   , "", "", 0}, \
   {"info"   , "", "", 0}  \
}

#define QIO_USQCDPROPSOURCEFILETYPESTRING  \
 "USQCD_ColorVectorSource"

#define QIO_USQCDPROPSOURCEFILETYPE 0

QIO_USQCDPropSourceFileInfo *QIO_create_usqcd_propsourcefile_info(char *info);
