#include "extdll.h"
#include "util.h"
#include "cbase.h"
#include "weapons.h"
#include "player.h"
#include "skill.h"
#include "items.h"
#include "gamerules.h"

class CFuncSlowMo : public CBaseEntity
{
public:
	void KeyValue(KeyValueData *pkvd);
	void Use(CBaseEntity *pActivator, CBaseEntity *pCaller, USE_TYPE useType, float value);

private:
	BOOL bActive;
	float fFactor;
};

LINK_ENTITY_TO_CLASS(func_slowmo, CFuncSlowMo);

void CFuncSlowMo::KeyValue(KeyValueData *pkvd)
{
	if(FStrEq(pkvd->szKeyName, "factor"))
	{
		fFactor = atof(pkvd->szValue);
		
		if((fFactor > 100.0f) || (fFactor < 0.0f))
			fFactor = 100.0f;

		fFactor = fFactor * 0.0001f;

		pkvd->fHandled = TRUE;
	}
	else 
		CBaseEntity::KeyValue( pkvd );

}

void CFuncSlowMo::Use(CBaseEntity *pActivator, CBaseEntity *pCaller, USE_TYPE useType, float value)
{
	if(bActive)
	{
		SERVER_COMMAND("host_framerate 0\n");
		bActive = FALSE;
	}
	else
	{
		char szCommand[60];
		memset(&szCommand, 0, 60);
		sprintf(szCommand, "host_framerate %f\n", fFactor);

		SERVER_COMMAND(szCommand);

		bActive = TRUE;
	}
	
}

/*
	char *pszCommand = new char[60];
	memset(pszCommand, 0, 60);

	sprintf(pszCommand, "vsys setvol %u %u\n", iChannel, iVolume);

	pClient = g_engfuncs.pfnPEntityOfEntIndex(1);

	if(!pClient)
		return;
	CLIENT_COMMAND(pClient, pszCommand);

	delete[] pszCommand;

	CVSysEntity::Use(pActivator, pCaller, useType, value);
	*/