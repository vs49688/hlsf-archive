#include "extdll.h"
#include "util.h"
#include "cbase.h"
#include "player.h"
#include "saverestore.h"
#include "gamerules.h"
#include "talkmonster.h"
#include "movewith.h" //LRC
#include "locus.h" //LRC

LINK_ENTITY_TO_CLASS(env_deadplayer, CDeadPlayer);

char *CDeadPlayer::szDeathAnimations[] =
{
	"lying_on_back",
	"lying_on_stomach"
};

void CDeadPlayer::KeyValue( KeyValueData *pkvd )
{
	if (FStrEq(pkvd->szKeyName, "type"))
	{
		m_iType = atoi(pkvd->szValue);
	}
	else if (FStrEq(pkvd->szKeyName, "pose"))
	{
		m_iPose = atoi(pkvd->szValue);
	}
	else
	{
		CBaseEntity::KeyValue( pkvd );
	}
}

void CDeadPlayer::Precache()
{
	PRECACHE_MODEL("models/eli_dead_sci.mdl");
	PRECACHE_MODEL("models/eli_dead_hev.mdl");
}

void CDeadPlayer :: Spawn( )
{
	Precache();

	if (m_iType)
		SET_MODEL(ENT(pev), "models/eli_dead_hev.mdl"); //LRC
	else
		SET_MODEL(ENT(pev), "models/eli_dead_sci.mdl");
	
	pev->effects		= 0;
	pev->sequence		= 0;
	// Corpses have less health
	pev->health			= 100;//gSkillData.scientistHealth;
	
	m_bloodColor = BLOOD_COLOR_RED;

	pev->body = 0;

	if((m_iPose >=0) && (m_iPose <= 1))
		m_iSequence = m_iPose;
	else
		m_iSequence = RANDOM_LONG(0, 1);

	pev->sequence = LookupSequence( szDeathAnimations[m_iSequence] );

	if (pev->sequence == -1)
	{
		ALERT ( at_debug, "Dead player with bad pose\n" );
	}

	//	pev->skin += 2; // use bloody skin -- UNDONE: Turn this back on when we have a bloody skin again!
	MonsterInitDead();
}

CDeadPlayer * CDeadPlayer::Create( int iBodyType, const Vector &vecOrigin, const Vector &vecAngles, edict_t *pentOwner )
{
	edict_t	*pent;
	CDeadPlayer *pEntity;

	pent = CREATE_NAMED_ENTITY( MAKE_STRING( "env_deadplayer" ));
	if ( FNullEnt( pent ) )
	{
		ALERT ( at_debug, "NULL Ent in Create!\n" );
		return NULL;
	}
	pEntity = (CDeadPlayer*)Instance( pent );
	pEntity->pev->owner = pentOwner;
	pEntity->pev->origin = vecOrigin;
	pEntity->pev->angles = vecAngles;
	pEntity->m_iType = iBodyType;
	pEntity->m_iPose = -1;
	DispatchSpawn( pEntity->edict() );
	return pEntity;
}
