#include "CClientHook.h"
#include <windows.h>

CClientHook gClient;
CClientHook *g_pClient = &gClient;

#define ASSIGNCHECK_LPFN(lpfn, name, type) lpfn = (type)GetProcAddress(hModule, name); if(!lpfn) return 1;

CClientHook::CClientHook()
{
	memset(this, 0, sizeof(CClientHook));
}

int CClientHook::InitHook()
{
	if(hModule)		// Just protection in case we're called more than once
		return 0;

	hModule = GetModuleHandle("client.dll");

	if(!hModule)
		return 1;

	ASSIGNCHECK_LPFN(VSys_PrecacheSound, "VSys_PrecacheSound", I_VSys_CH);
	ASSIGNCHECK_LPFN(VSys_SetupChannel, "VSys_SetupChannel", SETUPCHANNEL);
	//ASSIGNCHECK_LPFN(VSys_ClientDeath, "VSys_ClientDeath", CLIENTDEATH);

	bHookState = true;
	return 0;
}

bool CClientHook::GetHookState()
{
	return bHookState;
}