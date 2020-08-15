#ifndef WIN32FIX_H__
#define WIN32FIX_H__
#include "cmdlib.h" //--vluzacn

#if _MSC_VER >= 1000
#pragma once
#endif

#include <malloc.h>

/////////////////////////////
#ifdef SYSTEM_WIN32

#define alloca      _alloca

#define strncasecmp _strnicmp
#define strcasecmp  _stricmp

#define snprintf  _snprintf
#if _MSC_VER < 1400 // AdamR: Ignore this definition in Visual Studio 2005 and later
#define vsnprintf _vsnprintf
#endif

#define finite    _finite

#define rotl      _rotl
#define rotr      _rotr

#undef STDCALL
#undef FASTCALL
#undef CDECL

#define STDCALL     __stdcall
#define FASTCALL    __fastcall
#define CDECL       __cdecl

#define INLINE      __inline

#endif
/////////////////////////////

/////////////////////////////
#ifdef SYSTEM_POSIX
#define _MAX_PATH  4096
#define _MAX_DRIVE 4096
#define _MAX_DIR   4096
#define _MAX_FNAME 4096
#define _MAX_EXT   4096

#define STDCALL
#define FASTCALL
#define CDECL

#define INLINE inline

#endif
/////////////////////////////

#endif ///////////////////////////// WIN32FIX_H__
