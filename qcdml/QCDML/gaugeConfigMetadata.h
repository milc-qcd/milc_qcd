#include <xmlWriter.h>

// types
typedef struct gaugeConfigurationMetadata_tag       gaugeConfigurationMetadata;
typedef struct gaugeConfigurationManagement_tag     gaugeConfigurationManagement; 
typedef struct Implementation_tag                   Implementation;
typedef struct gaugeConfigurationAlgorithm_tag      gaugeConfigurationAlgorithm;
typedef struct gaugeConfigurationMarkovStep_tag     gaugeConfigurationMarkovStep;
typedef struct gaugeConfigurationHistory_tag        gaugeConfigurationHistory;
typedef struct gaugeConfigurationHistoryEntry_tag   gaugeConfigurationHistoryEntry;
typedef struct participantDescription_tag           participantDescription;
typedef struct machineDescription_tag               machineDescription;
typedef struct codeDescription_tag                  codeDescription;
typedef struct Parameters_tag                       Parameters;
typedef struct Parameter_tag                        Parameter;

// definitions

struct participantDescription_tag
{
  CstringRef name;
  CstringRef institution;
};

struct gaugeConfigurationHistoryEntry_tag
{
  unsigned revision;
  CstringRef revisionAction;
  participantDescription participant;
  CstringRef date;
  gaugeConfigurationHistoryEntry* next;
  gaugeConfigurationHistoryEntry* prev;
};

struct gaugeConfigurationHistory_tag
{
  gaugeConfigurationHistoryEntry* head;
  gaugeConfigurationHistoryEntry* tail;
};

struct gaugeConfigurationManagement_tag
{
  unsigned                       revisions;
  unsigned                       crcCheckSum;
  gaugeConfigurationHistory      archiveHistory;
};

struct machineDescription_tag
{
  CstringRef name;
  CstringRef institution;
  CstringRef machineType;
};

struct codeDescription_tag
{
  CstringRef name;
  CstringRef version;
  CstringRef date;
};

struct Implementation_tag
{
  machineDescription machine;
  codeDescription    code;
};

struct Parameter_tag
{
  CstringRef name;
  CstringRef value;
  Parameter* next;
  Parameter* prev;
};

struct Parameters_tag
{
  Parameter* head;
  Parameter* tail;
};

struct gaugeConfigurationAlgorithm_tag
{
  Parameters parameters;
};

struct gaugeConfigurationMarkovStep_tag
{
  CstringRef markovChainURI;
  unsigned   series;
  unsigned   update;
  float      avePlaquette;
  CstringRef dataLFN;
};

struct gaugeConfigurationMetadata_tag
{
  gaugeConfigurationManagement     management;
  Implementation                   implementation;
  gaugeConfigurationAlgorithm      algorithm;
  CstringRef                       precision;
  gaugeConfigurationMarkovStep     markovStep;
};

// prototypes
void gaugeConfigurationHistoryInit ( gaugeConfigurationHistory* list );

void gaugeConfigurationHistoryPushBack ( gaugeConfigurationHistory* list,
					 unsigned revision,
					 CstringRef revisionAction,
					 CstringRef participantName,
					 CstringRef participantInstitution,
					 CstringRef date
					 );

void parametersInit ( Parameters* list );

void parametersPushBack ( Parameters* list,
			  CstringRef name,
			  CstringRef value
			  );

void gaugeConfigurationMetadataEncode ( outStream* os, char const* name,
					void const* obj );

void gaugeConfigurationManagementEncode ( outStream* os, char const* name,
					  void const* obj );

void ImplementationEncode ( outStream* os, char const* name,
			    void const* obj );

void gaugeConfigurationAlgorithmEncode ( outStream* os, char const* name,
					 void const* obj );

void gaugeConfigurationMarkovStepEncode ( outStream* os, char const* name,
					  void const* obj );

void gaugeConfigurationHistoryEncode ( outStream* os, char const* name,
				       void const* obj );

void gaugeConfigurationHistoryEntryEncode ( outStream* os, char const* name,
					    void const* obj );

void participantDescriptionEncode ( outStream* os, char const* name,
				    void const* obj );

void machineDescriptionEncode ( outStream* os, char const* name,
				void const* obj );

void codeDescriptionEncode ( outStream* os, char const* name,
			     void const* obj );

void ParametersEncode ( outStream* os, char const* name,
			void const* obj );

void ParameterEncode ( outStream* os, char const* name,
		       void const* obj );
