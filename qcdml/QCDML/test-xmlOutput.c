#include <gaugeConfigMetadata.h>

gaugeConfigurationMetadata gaugeConfiguration;

void initMD ( )
{
  // initialize the metadata
  gaugeConfiguration.management.revisions = 0;
  gaugeConfiguration.management.crcCheckSum = 0x126abb95;
  gaugeConfigurationHistoryInit (
		   &gaugeConfiguration.management.archiveHistory );
  gaugeConfigurationHistoryPushBack (
		   &gaugeConfiguration.management.archiveHistory,
		   0,                           // revision
		   "generate",                  // action
		   "John Q. Smith", "Fermilab", // actor
		   "2003-12-11T10:25:52Z"       // date
		   );

  gaugeConfiguration.implementation.machine.name =
    "lqcd.fnal.gov::qcd0203+qcd0204+qcd0208+qcd0201";
  gaugeConfiguration.implementation.machine.institution =
    "Fermilab";
  gaugeConfiguration.implementation.machine.machineType =
    "cluster::P4E+infiniband";

  gaugeConfiguration.implementation.code.name = "milc::su3_rmd_symzk1_asqtad";
  gaugeConfiguration.implementation.code.version = "6.22.1-SSE3";
  gaugeConfiguration.implementation.code.date = "2003-11-22";

  parametersInit     ( &gaugeConfiguration.algorithm.parameters );
  parametersPushBack ( &gaugeConfiguration.algorithm.parameters,
		       "microcanonical_time_step", "0.01" );
  parametersPushBack ( &gaugeConfiguration.algorithm.parameters,
		       "steps_per_trajectory", "10" );
  parametersPushBack ( &gaugeConfiguration.algorithm.parameters,
		       "max_cg_iterations", "900" );
  parametersPushBack ( &gaugeConfiguration.algorithm.parameters,
		       "error_per_site", "0.4e-7" );
  parametersPushBack ( &gaugeConfiguration.algorithm.parameters,
		       "error_for_propagator", "0.2e-7" );

  gaugeConfiguration.precision = "single";

  gaugeConfiguration.markovStep.markovChainURI =
    "http://qcdgrid.fnal.gov/milc/l2896f21b709m0062m031/coulomb";
  gaugeConfiguration.markovStep.series = 0;
  gaugeConfiguration.markovStep.update = 20986;
  gaugeConfiguration.markovStep.avePlaquette = 1.7845913;
  gaugeConfiguration.markovStep.dataLFN = "LQCD001A00BFF32";
}

void testFile ( )
{
  outStream* os = new_outStreamFile ( stdout );

  printf ( "####### file ############\n" );
  gaugeConfigurationMetadataEncode ( os, "gaugeConfiguration",
				     &gaugeConfiguration );

  delete_outStream ( os );
}

void testString ( )
{
  size_t size = 4000;
  String* st = new_String ( size, "\n####### String ############\n" );
  outStream* os = new_outStreamString ( st );

  gaugeConfigurationMetadataEncode ( os, "gaugeConfiguration",
				     &gaugeConfiguration );

  delete_outStream ( os );

  printf ( "%s\n", st -> chars );
  delete_String ( st );
}

int main ( )
{

  initMD ( );

  testFile ( );

  testString ( );

  return 0;
}

