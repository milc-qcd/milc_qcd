#include <gaugeConfigMetadata.h>
#include <stdlib.h>
#include <cast.h>

// doubly linked list objects have data members 'head' and 'tail'
#define listInit(obj) (obj) -> tail = (obj) -> head = NULL

// list member objects have data members 'next' and 'prev'
#define listPushBack(list,obj) \
{ \
  (obj) -> next = NULL; \
  (obj) -> prev = (list) -> tail; \
  (list) -> tail = (obj); \
  if ( (obj) -> prev ) \
    (obj) -> prev -> next = (obj); \
  if ( (list) -> head == NULL ) \
    (list) -> head = (obj); \
}

// list forward iteration
#define foreachForward(iterator,list) \
for ( iterator = list -> head; iterator != NULL; iterator = iterator -> next )


void gaugeConfigurationHistoryInit ( gaugeConfigurationHistory* list )
{
  listInit ( list );
}

void gaugeConfigurationHistoryPushBack ( gaugeConfigurationHistory* list,
					 unsigned revision,
					 CstringRef revisionAction,
					 CstringRef participantName,
					 CstringRef participantInstitution,
					 CstringRef date
					 )
{
  // new
  gaugeConfigurationHistoryEntry* p =
    static_cast ( gaugeConfigurationHistoryEntry*,
		  calloc ( 1, sizeof(gaugeConfigurationHistoryEntry)));

  // set payload
  p -> revision = revision;
  p -> revisionAction = revisionAction;
  p -> participant.name = participantName;
  p -> participant.institution = participantInstitution;
  p -> date = date;

  // insert at back of the list
  listPushBack(list,p);
}

void parametersInit ( Parameters* list )
{
  listInit ( list );
}

void parametersPushBack ( Parameters* list,
			  CstringRef name,
			  CstringRef value
			  )
{
  // new
  Parameter* p = static_cast ( Parameter*,
			       calloc ( 1, sizeof(Parameter) ) );

  // set payload
  p -> name = name;
  p -> value = value;

  // insert at back of the list
  listPushBack(list,p);
}

void gaugeConfigurationMetadataEncode ( outStream* os, char const* name,
					void const* obj )
{
  gaugeConfigurationMetadata const* p =
    static_cast ( gaugeConfigurationMetadata const*, obj );

  os -> print ( os, "<%s>", name );
  gaugeConfigurationManagementEncode ( os, "management", &p -> management );
  ImplementationEncode ( os, "implementation", &p -> implementation );
  gaugeConfigurationAlgorithmEncode ( os, "algorithm", &p -> algorithm );
  CstringRefEncode ( os, "precision", &p -> precision );
  gaugeConfigurationMarkovStepEncode ( os, "markovStep", &p -> markovStep );
  os -> print ( os, "</%s>", name );
}

void gaugeConfigurationManagementEncode ( outStream* os, char const* name,
					  void const* obj )
{
  gaugeConfigurationManagement const* p =
    static_cast (gaugeConfigurationManagement const*, obj );

  os -> print ( os, "<%s>", name );
  unsignedEncode ( os, "revisions", &p -> revisions );
  unsignedEncode ( os, "crcCheckSum", &p -> crcCheckSum );
  gaugeConfigurationHistoryEncode ( os, "archiveHistory", &p -> archiveHistory );
  os -> print ( os, "</%s>", name );
}

void ImplementationEncode ( outStream* os, char const* name,
			    void const* obj )
{
  Implementation const* p =
    static_cast ( Implementation const*, obj );

  os -> print ( os, "<%s>", name );
  machineDescriptionEncode ( os, "machine", &p -> machine );
  codeDescriptionEncode ( os, "code", &p -> code );
  os -> print ( os, "</%s>", name );
}

void gaugeConfigurationAlgorithmEncode ( outStream* os, char const* name,
					 void const* obj )
{
  gaugeConfigurationAlgorithm const* p =
    static_cast ( gaugeConfigurationAlgorithm const*, obj );

  os -> print ( os, "<%s>", name );
  ParametersEncode ( os, "parameters", &p -> parameters );
  os -> print ( os, "</%s>", name );
}

void gaugeConfigurationMarkovStepEncode ( outStream* os, char const* name,
					  void const* obj )
{
  gaugeConfigurationMarkovStep const* p =
    static_cast ( gaugeConfigurationMarkovStep const*, obj );

  os -> print ( os, "<%s>", name );
  CstringRefEncode ( os, "markovChainURI", &p -> markovChainURI );
  unsignedEncode ( os, "series", &p -> series );
  unsignedEncode ( os, "update", &p -> update );
  floatEncode ( os, "avePlaquette", &p -> avePlaquette );
  CstringRefEncode ( os, "dataLFN", &p -> dataLFN );
  os -> print ( os, "</%s>", name );
}

void gaugeConfigurationHistoryEncode ( outStream* os, char const* name,
				       void const* obj )
{
  gaugeConfigurationHistory const* p =
    static_cast ( gaugeConfigurationHistory const*, obj );
  gaugeConfigurationHistoryEntry* elem;

  os -> print ( os, "<%s>", name );
  foreachForward ( elem, p )
    gaugeConfigurationHistoryEntryEncode ( os, "elem", elem );
  os -> print ( os, "</%s>", name );
}

void gaugeConfigurationHistoryEntryEncode ( outStream* os, char const* name,
					    void const* obj )
{
  gaugeConfigurationHistoryEntry const* p =
    static_cast ( gaugeConfigurationHistoryEntry const*, obj );

  os -> print ( os, "<%s>", name );
  unsignedEncode ( os, "revision", &p -> revision );
  CstringRefEncode  ( os, "revisionAction", &p -> revisionAction );
  participantDescriptionEncode ( os, "participant", &p -> participant );
  CstringRefEncode  ( os, "date", &p -> date );
  os -> print ( os, "</%s>", name );
}

void participantDescriptionEncode ( outStream* os, char const* name,
				    void const* obj )
{
  participantDescription const* p =
    static_cast ( participantDescription const*, obj );

  os -> print ( os, "<%s>", name );
  CstringRefEncode ( os, "name",  &p -> name );
  CstringRefEncode ( os, "institution",  &p -> institution );
  os -> print ( os, "</%s>", name );
}

void machineDescriptionEncode ( outStream* os, char const* name,
				void const* obj )
{
  machineDescription const* p =
    static_cast ( machineDescription const*, obj );

  os -> print ( os, "<%s>", name );
  CstringRefEncode ( os, "name",  &p -> name );
  CstringRefEncode ( os, "institution",  &p -> institution );
  CstringRefEncode ( os, "machineType",  &p -> machineType );
  os -> print ( os, "</%s>", name );
}

void codeDescriptionEncode ( outStream* os, char const* name,
			     void const* obj )
{
  codeDescription const* p =
    static_cast ( codeDescription const*, obj );

  os -> print ( os, "<%s>", name );
  CstringRefEncode ( os, "name",  &p -> name );
  CstringRefEncode ( os, "version",  &p -> version );
  CstringRefEncode ( os, "date",  &p -> date );
  os -> print ( os, "</%s>", name );
}

void ParametersEncode ( outStream* os, char const* name,
			void const* obj )
{
  Parameters const* p =
    static_cast ( Parameters const*, obj );
  Parameter* elem;

  os -> print ( os, "<%s>", name );
  foreachForward ( elem, p )
    ParameterEncode ( os, "elem", elem );
  os -> print ( os, "</%s>", name );
}

void ParameterEncode ( outStream* os, char const* name,
		       void const* obj )
{
  Parameter const* p =
    static_cast ( Parameter const*, obj );

  os -> print ( os, "<%s>", name );
  CstringRefEncode ( os, "name",  &p -> name );
  CstringRefEncode ( os, "value", &p -> value );
  os -> print ( os, "</%s>", name );
}
