#include "cmdlib.h"
#include "messages.h"
#include "win32fix.h"
#include "log.h"
#include "hlassert.h"
#include "mathlib.h"
#include "scriplib.h"
#include "winding.h"
#include "threads.h"
#include "bspfile.h"
#include "blockmem.h"
#include "filelib.h"
#ifdef ZHLT_PARAMFILE
#include "cmdlinecfg.h"
#endif
#ifdef RIPENT_PAUSE
#include <conio.h>
#endif

#define DEFAULT_PARSE false
#define DEFAULT_CHART false
#define DEFAULT_INFO true

