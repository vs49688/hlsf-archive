#include <intrin.h>

#include "extdll.h"
#include "util.h"
#include "cbase.h"
#include "monsters.h"
#include "weapons.h"
#include "nodes.h"
#include "player.h"

#include "usercmd.h"
#include "entity_state.h"
#include "demo_api.h"
#include "pm_defs.h"
#include "event_api.h"
#include "r_efx.h"

#include "VSSound.h"

extern "C"
{
	__declspec(dllexport) int VSys_PrecacheSound(const char* pszFile);
	__declspec(dllexport) int VSys_SetupChannel( int iChannel, const char *pszSong, bool bPaused, unsigned int mode );
}

/*
void VSys_InitServerHook()
{
	if(g_pServer->GetModule())
		g_pServer->InitHook();
}

void VSys_SetServerInstance(HINSTANCE hServerInstance)
{
	if(!hServerInstance)
		return;
	g_pServer->SetServerInstance(hServerInstance);
}
*/
// Global interface to emulate precaching
int VSys_PrecacheSound(const char* pszFile)
{
	if(strcmp(pszFile, "") == 0)
		return -1;

	int iResult = gMP3.CheckCached(pszFile);

	if(!((iResult >= 0) && (iResult < VSYS_MAX_SOUNDS)))
	{
		iResult = gMP3.CacheSound(pszFile);

		if(!((iResult >= 0) && (iResult < VSYS_MAX_SOUNDS)))
			return iResult;
	}
	return iResult;
}

int VSys_SetupChannel( int iChannel, const char *pszSong, bool bPaused, unsigned int mode )
{
	return gMP3.SetupChannel(iChannel, pszSong, bPaused, mode);
}

//VSys_ClientDeath has been moved to hl_weapons.cpp
