#ifndef DSLASH_KS_REDEFINE_H
#define DSLASH_KS_REDEFINE_H

/* For Kogut Susskind fermions.
   Normally the generic dslash_site invokes the unimproved KS Dslash.

   Here we map the generic dslash macros to a specific FN or EO
   version for improved actions. */

#ifdef FN
#define dslash_site dslash_fn_site
#define dslash_field dslash_fn_field
#endif
#ifdef EO
#define dslash_site dslash_eo_site
#define dslash_field dslash_eo_field
#endif

#endif
