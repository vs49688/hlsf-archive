#include "extdll.h"
#include "util.h"
#include "cbase.h"
#include "weapons.h"
#include "player.h"
#include "skill.h"
#include "items.h"
#include "gamerules.h"

#define INPUT_NORMAL	0
#define INPUT_FREE		1

extern int gmsgSetInputType;

class CInputSetType : public CBaseEntity
{
public:
	void			Spawn( void );
	void			KeyValue( KeyValueData *pkvd );
	void			Use( CBaseEntity *pActivator, CBaseEntity *pCaller, USE_TYPE useType, float value );

private:
	int				iInputType;
	int				iMovement;
};

LINK_ENTITY_TO_CLASS(input_settype, CInputSetType);

void CInputSetType::Spawn()
{
}

void CInputSetType::KeyValue(KeyValueData *pkvd)
{
	if(FStrEq(pkvd->szKeyName, "input"))
	{
		iInputType = atoi(pkvd->szValue);

		pkvd->fHandled = TRUE;
	}
	else 
		CBaseEntity::KeyValue( pkvd );
}

void CInputSetType::Use(CBaseEntity *pActivator, CBaseEntity *pCaller, USE_TYPE useType, float value)
{
	pActivator->pev->movetype = MOVETYPE_FLY;
	MESSAGE_BEGIN( MSG_ONE, gmsgSetInputType, NULL, pActivator->pev );
		WRITE_LONG((long)iInputType);
	MESSAGE_END();
}