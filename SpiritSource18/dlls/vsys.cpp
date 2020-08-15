#include "extdll.h"
#include "util.h"
#include "cbase.h"
#include "player.h"
#include "saverestore.h"
#include "gamerules.h"
#include "talkmonster.h"
#include "movewith.h"
#include "locus.h"

#include "../common/vsys/vsys_common.h"
#include "CClientHook.h"

#define FMOD_LOOP_OFF                  0x00000001  /* For non looping sounds. (DEFAULT).  Overrides FMOD_LOOP_NORMAL / FMOD_LOOP_BIDI. */
#define FMOD_LOOP_NORMAL               0x00000002  /* For forward looping sounds. */

class CVSysEntity : public CBaseEntity
{
public:
	void			Spawn( void );
	void			KeyValue( KeyValueData *pkvd );
	void			Use( CBaseEntity *pActivator, CBaseEntity *pCaller, USE_TYPE useType, float value );

	virtual int		Save( CSave &save );
	virtual int		Restore( CRestore &restore );
	static	TYPEDESCRIPTION m_SaveData[];

protected:
	int				iChannel;
	BOOL			m_bPlaying;
	edict_t			*pClient;
};


TYPEDESCRIPTION	CVSysEntity::m_SaveData[] = 
{
	DEFINE_ARRAY( CVSysEntity, iChannel, FIELD_INTEGER, sizeof(int) ),
};

IMPLEMENT_SAVERESTORE( CVSysEntity, CBaseEntity );

void CVSysEntity::Spawn( void )
{
    pev->solid = SOLID_NOT;
    pev->movetype = MOVETYPE_NONE;
    m_bPlaying = FALSE;
}

void CVSysEntity::KeyValue( KeyValueData *pkvd )
{
	if (FStrEq(pkvd->szKeyName, "Channel"))
	{
		if((iChannel > VSYS_MAX_CHANNELS) && ((iChannel < 0) && (iChannel != VSYS_CHANNELID_ALL)))
			iChannel = 0;

		iChannel = atoi(pkvd->szValue);
		pkvd->fHandled = TRUE;
	}
	else
		CBaseEntity::KeyValue( pkvd );
}

void CVSysEntity::Use( CBaseEntity *pActivator, CBaseEntity *pCaller, USE_TYPE useType, float value )
{
	CBaseEntity::Use(pActivator, pCaller, useType, value);
}



class CVSysVol : public CVSysEntity
{
public:
	void			KeyValue( KeyValueData *pkvd );
	void			Use( CBaseEntity *pActivator, CBaseEntity *pCaller, USE_TYPE useType, float value );

private:
	unsigned int	iVolume;
};

LINK_ENTITY_TO_CLASS(vsys_changevol, CVSysVol);

void CVSysVol::KeyValue( KeyValueData *pkvd )
{
	if(FStrEq(pkvd->szKeyName, "Volume"))
	{
		iVolume = atoi(pkvd->szValue);
		
		if(iVolume > 100)
			iVolume = 100;

		pkvd->fHandled = TRUE;
	}
	else 
		CVSysVol::KeyValue( pkvd );

}

void CVSysVol::Use( CBaseEntity *pActivator, CBaseEntity *pCaller, USE_TYPE useType, float value )
{

	char *pszCommand = new char[60];
	memset(pszCommand, 0, 60);

	sprintf(pszCommand, "vsys setvol %u %u\n", iChannel, iVolume);

	pClient = g_engfuncs.pfnPEntityOfEntIndex(1);

	if(!pClient)
		return;
	CLIENT_COMMAND(pClient, pszCommand);

	delete[] pszCommand;

	CVSysEntity::Use(pActivator, pCaller, useType, value);
}


class CVSysToggle : public CVSysEntity
{
public:
	void			Use( CBaseEntity *pActivator, CBaseEntity *pCaller, USE_TYPE useType, float value );

};

LINK_ENTITY_TO_CLASS(vsys_toggle, CVSysToggle);

void CVSysToggle::Use( CBaseEntity *pActivator, CBaseEntity *pCaller, USE_TYPE useType, float value )
{

	char *pszCommand = new char[60];
	memset(pszCommand, 0, 60);

	sprintf(pszCommand, "vsys toggle %u\n", iChannel);

	pClient = g_engfuncs.pfnPEntityOfEntIndex(1);

	if(!pClient)
		return;
	CLIENT_COMMAND(pClient, pszCommand);

	delete[] pszCommand;

	CVSysEntity::Use(pActivator, pCaller, useType, value);
}

class CVSysPause : public CVSysEntity
{
	void			Use( CBaseEntity *pActivator, CBaseEntity *pCaller, USE_TYPE useType, float value );
};

LINK_ENTITY_TO_CLASS(vsys_pause, CVSysPause);

void CVSysPause::Use( CBaseEntity *pActivator, CBaseEntity *pCaller, USE_TYPE useType, float value )
{
	char *pszCommand = new char[60];
	memset(pszCommand, 0, 60);

	sprintf(pszCommand, "vsys pause %u\n", iChannel);
	CLIENT_COMMAND(pClient, pszCommand);

	delete[] pszCommand;

	CVSysEntity::Use(pActivator, pCaller, useType, value);
}

class CVSysResume : public CVSysEntity
{
	void			Use( CBaseEntity *pActivator, CBaseEntity *pCaller, USE_TYPE useType, float value );
};

LINK_ENTITY_TO_CLASS(vsys_resume, CVSysResume);

void CVSysResume::Use( CBaseEntity *pActivator, CBaseEntity *pCaller, USE_TYPE useType, float value )
{
	char *pszCommand = new char[60];
	memset(pszCommand, 0, 60);

	sprintf(pszCommand, "vsys resume %u\n", iChannel);

	pClient = g_engfuncs.pfnPEntityOfEntIndex(1);

	if(!pClient)
		return;
	CLIENT_COMMAND(pClient, pszCommand);

	delete[] pszCommand;

	CVSysEntity::Use(pActivator, pCaller, useType, value);
}

class CVSysObject : public CVSysEntity
{
public:
	void			Spawn( void );
	void			KeyValue( KeyValueData *pkvd );
	void			Use( CBaseEntity *pActivator, CBaseEntity *pCaller, USE_TYPE useType, float value );
	void			Precache( void );

private:
	BOOL			bLooped;
	BOOL			bPaused;
	BOOL			bState;
	char			szFile[MAX_PATH];
};

LINK_ENTITY_TO_CLASS(vsys_object, CVSysObject);

void CVSysObject::Spawn( void )
{
	Precache();
	
	unsigned int mode;

	if(bLooped)
		mode = FMOD_LOOP_NORMAL;
	else
		mode = FMOD_LOOP_OFF;

	bState = bPaused;

	if(g_pClient->GetHookState())
		g_pClient->VSys_SetupChannel(iChannel, (const char*)&szFile, bPaused, mode);
	//int VSys_SetupChannel( int iChannel, const char *pszSong, bool bPaused, unsigned int mode )
	CVSysEntity::Spawn();
}
void CVSysObject::Precache( void )
{
	if(g_pClient->GetHookState())
		g_pClient->VSys_PrecacheSound((const char*)&szFile);
	/*char szCommand[60];
	memset(&szCommand, 0, 60);

	sprintf(szCommand, "vsys cache \"%s\"\n", (char*)&szFile);

	pClient = g_engfuncs.pfnPEntityOfEntIndex(1);

	if(!pClient)
		return;
	CLIENT_COMMAND(pClient, szCommand);*/
}

void CVSysObject::KeyValue( KeyValueData *pkvd )
{
	if(FStrEq(pkvd->szKeyName, "paused"))
	{
		if(atoi(pkvd->szValue) == TRUE)
			bPaused = TRUE;
		else
			bPaused = FALSE;

		pkvd->fHandled = TRUE;
	}
	else if(FStrEq(pkvd->szKeyName, "looped"))
	{
		if(atoi(pkvd->szValue) == TRUE)
			bLooped = TRUE;
		else
			bLooped = FALSE;

		pkvd->fHandled = TRUE;
	}
	else if(FStrEq(pkvd->szKeyName, "file"))
	{
		memset(&szFile, 0, MAX_PATH);
		memcpy(&szFile, pkvd->szValue, MAX_PATH);

		pkvd->fHandled = TRUE;

	}
	else 
		CVSysObject::KeyValue( pkvd );

}

void CVSysObject::Use( CBaseEntity *pActivator, CBaseEntity *pCaller, USE_TYPE useType, float value )
{
	char szCommand[60];
	memset(&szCommand, 0, 60);

	if(bState)
		sprintf(szCommand, "vsys play %u\n", iChannel);
	else
		sprintf(szCommand, "vsys stop %u\n", iChannel);

	pClient = g_engfuncs.pfnPEntityOfEntIndex(1);

	if(!pClient)
		return;

	CLIENT_COMMAND(pClient, szCommand);

	bState = !bState;
		

	//char *pszCommand = new char[MAX_PATH];
	//memset(pszCommand, 0, MAX_PATH);

	//sprintf(pszCommand, "vsys setup %u %u %u \"%s\"\n", iChannel, bPaused, bLooped, (char*)&szFile);



	//delete[] pszCommand;

	CVSysEntity::Use(pActivator, pCaller, useType, value);
}

class CVSysStop : public CVSysEntity
{
	void			Use( CBaseEntity *pActivator, CBaseEntity *pCaller, USE_TYPE useType, float value );
};

LINK_ENTITY_TO_CLASS(vsys_stop, CVSysStop);

void CVSysStop::Use( CBaseEntity *pActivator, CBaseEntity *pCaller, USE_TYPE useType, float value )
{
	char *pszCommand = new char[60];
	memset(pszCommand, 0, 60);

	sprintf(pszCommand, "vsys stop %u\n", iChannel);

	pClient = g_engfuncs.pfnPEntityOfEntIndex(1);

	if(!pClient)
		return;
	CLIENT_COMMAND(pClient, pszCommand);

	delete[] pszCommand;

	CVSysEntity::Use(pActivator, pCaller, useType, value);
}


class CVSysPlay : public CVSysEntity
{
	void			Use( CBaseEntity *pActivator, CBaseEntity *pCaller, USE_TYPE useType, float value );
};

LINK_ENTITY_TO_CLASS(vsys_play, CVSysPlay);

void CVSysPlay::Use( CBaseEntity *pActivator, CBaseEntity *pCaller, USE_TYPE useType, float value )
{
	char *pszCommand = new char[60];
	memset(pszCommand, 0, 60);

	sprintf(pszCommand, "vsys play %u\n", iChannel);

	pClient = g_engfuncs.pfnPEntityOfEntIndex(1);

	if(!pClient)
		return;
	CLIENT_COMMAND(pClient, pszCommand);

	delete[] pszCommand;

	CVSysEntity::Use(pActivator, pCaller, useType, value);
}