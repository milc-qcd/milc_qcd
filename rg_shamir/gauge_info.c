/*********************** gauge_info.c *************************/
/* MIMD version 7 */

/* For rg_shamir */

/* Application-dependent routine for writing gauge info file */
/* This file is an ASCII companion to the gauge configuration file
   and contains information about the action used to generate it.
   This information is consistently written in the pattern

       keyword  value

   or

       keyword[n] value1 value2 ... valuen

   where n is an integer.

   To maintain a semblance of consistency, the possible keywords are
   listed in io_lat.h.  Add more as the need arises, but be sure
   to notify the rest of the collaboration.

   */

#include "RG_Shamir_includes.h"
#include <quark_action.h>

/*---------------------------------------------------------------------------*/
/* This routine writes the ASCII info file.  It is called from one of
   the lattice output routines in io_lat4.c.*/

extern char gauge_action_description[128]; /* in gauge_stuff.c */
extern int gauge_action_nloops,gauge_action_nreps;
void write_appl_gauge_info(FILE *fp, gauge_file *gf)
{

  /* Write generic information */
  write_generic_gauge_info(fp, gf);

  /* The rest are optional */

  write_gauge_info_item(fp,"action.description","%s",
			"\"Gauge plus fermion\"",0,0);

  write_gauge_info_item(fp,"gauge.description","%s",
			gauge_action_description,0,0);
  write_gauge_info_item(fp,"gauge.nloops","%d",(char *)&gauge_action_nloops,0,0);
  write_gauge_info_item(fp,"gauge.nreps","%d",(char *)&gauge_action_nreps,0,0);
  write_gauge_info_item(fp,"gauge.tadpole.u0","%f",(char *)&u0,0,0);

  write_gauge_info_item(fp,"quark.description","%s",quark_action_description,0,0);
  write_gauge_info_item(fp,"quark.flavors","%d",(char *)&nflavors,0,0);
  write_gauge_info_item(fp,"quark.mass","%f",(char *)&mass,0,0);
}

#if 0
/*----------------------------------------------------------------------*/
/* This section constructs the QCDML string for a gauge configuration */
#include <gaugeConfigMetadata.h>
#define MAX_STRING 512

/* Fill a gauge configuration structure */
gaugeConfigurationMetadata *createGaugeMD ( )
{
  static gaugeConfigurationMetadata gaugeConfiguration;
  char string[MAX_STRING+1];  /* Temporary string */

  /* initialize the metadata */
  gaugeConfiguration.management.revisions = 0;
  /* We can't know the crcCheckSum until the file is written
     so this field has to be supplied later */
  gaugeConfiguration.management.crcCheckSum = 0;
  gaugeConfigurationHistoryInit (
		   &gaugeConfiguration.management.archiveHistory );
  /* TO DO: Fill in name of collaborator, his/her institution, and date */
  gaugeConfigurationHistoryPushBack (
		   &gaugeConfiguration.management.archiveHistory,
		   0,                           // revision
		   "generate",                  // action
		   "MILC Collaborator", "Some lab", // actor
		   "2003-12-11T10:25:52Z"       // date
		   );

  /* TO DO: Fill in node names somehow */
  gaugeConfiguration.implementation.machine.name =
    "lqcd.fnal.gov::qcd0203+qcd0204+qcd0208+qcd0201";
  /* TO DO: Fill in location of machine somehow */
  gaugeConfiguration.implementation.machine.institution =
    "Fermilab";
  /* TO DO: Fill in machine type */
  gaugeConfiguration.implementation.machine.machineType =
    "cluster::P4E+infiniband";

  /* TO DO: Fill in code name */
  gaugeConfiguration.implementation.code.name = "milc::su3_rmd_symzk1_asqtad";
  /* TO DO: Fill in code version and date */
  gaugeConfiguration.implementation.code.version = "6.22.1-SSE3";
  gaugeConfiguration.implementation.code.date = "2003-11-22";

  parametersInit     ( &gaugeConfiguration.algorithm.parameters );

  /* TO DO: Check on gauge couplings and u0 */
  snprintf(string,MAX_STRING,"%f",epsilon);
  parametersPushBack ( &gaugeConfiguration.algorithm.parameters,
		       "microcanonical_time_step", string );
  snprintf(string,MAX_STRING,"%d",steps);
  parametersPushBack ( &gaugeConfiguration.algorithm.parameters,
		       "steps_per_trajectory", string );
  snprintf(string,MAX_STRING,"%d",niter);
  parametersPushBack ( &gaugeConfiguration.algorithm.parameters,
		       "max_cg_iterations", string );
  snprintf(string,MAX_STRING,"%.2e",sqrt(rsqmin));
  parametersPushBack ( &gaugeConfiguration.algorithm.parameters,
		       "error_per_site", string );
  snprintf(string,MAX_STRING,"%.2e",sqrt(rsqprop));
  parametersPushBack ( &gaugeConfiguration.algorithm.parameters,
		       "error_for_propagator", string );

  gaugeConfiguration.precision = "single";

  /* TO DO: Need encoding for series */
  gaugeConfiguration.markovStep.markovChainURI =
    "http://qcdgrid.fnal.gov/milc/l2896f21b709m0062m031/coulomb";

  /* TO DO: Need series name */
  gaugeConfiguration.markovStep.series = 0;

  /* TO DO: Need sequence number */
  gaugeConfiguration.markovStep.update = 20986;

  gaugeConfiguration.markovStep.avePlaquette = 1.7845913;
  gaugeConfiguration.markovStep.dataLFN = "LQCD001A00BFF32";

  return &gaugeConfiguration;
}

String *createGaugeQCDML ( gaugeConfigurationMetadata *gaugeConfiguration )
{
  char dummy[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>Dummy QCDML</title>";
  size_t size = 4000;
  String* st = new_String ( size, "\n####### String ############\n" );
  outStream* os = new_outStreamString ( st );

  gaugeConfigurationMetadataEncode ( os, "gaugeConfiguration",
				     gaugeConfiguration );

  delete_outStream ( os );
}

void destroyGaugeQCDML(String *st){
  delete_String ( st );
}


#endif

char *create_QCDML(){
  char dummy[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>Dummy QCDML</title>";
  char *qcdml = (char *)malloc(sizeof(dummy)+1);
  strcpy(qcdml, dummy);
  return qcdml;
}

void free_QCDML(char *qcdml){
  if(qcdml != NULL)free(qcdml);
}

