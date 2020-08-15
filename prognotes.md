# Programmer's Notes

## HUD Colour:

The variable that contains the HUD colour is located in hud.cpp.
It is a hexadecimal value of the form `0x00RRGGBB`.

### hud.cpp

```c
m_iHUDColor = 0x00FFFFFF;
```
### Standard Values:

* Half-Life Orange:
  - 255, 160, 0
  - 0xFFA000
* Blue-Shift Blue:
  - 51, 102, 255
  - 0x3366FF
* Static-Friction White:
  - 255, 255, 255
  - 0xFFFFFF

## VSys Notes

By default, VSys stops all channels on every map load.
As this occurs in the client, a direct call to the `ParseString()` function, with `stop -1` as the argument, is the easiest and most direct way. However, this is also the slowest. Calling the `StopAll()` function is the fastest.

### hl_weapons.cpp

```c++
void CBasePlayer::Spawn( void )
{
	if (m_pActiveItem&&m_pNextItem)
	{
		m_pActiveItem = m_pNextItem;
		m_pActiveItem->Deploy();
		m_pActiveItem->UpdateItemInfo();
		m_pNextItem = NULL;
	}
	g_irunninggausspred = false;
	gMP3.StopAll();
	gMP3.DeathStop();
}
```

All the code relating to the new death scene can be found in `CBasePlayer::StartDeath()`.

### player.cpp:

```c++
void CBasePlayer::StartDeath()
{
	// Nasty way of hiding the player model
	pev->renderamt = 0;
	pev->rendermode = 2; // Texture Rendering Mode

	m_iHideHUD = (HIDEHUD_WEAPONS | HIDEHUD_FLASHLIGHT | HIDEHUD_HEALTH | HIDEHUD_CUSTOMCROSSHAIR | HIDEHUD_ARMOR);

	MESSAGE_BEGIN( MSG_ONE, gmsgStartDeath, NULL, edict() );
		WRITE_BYTE(0);
	MESSAGE_END();

	UTIL_ScreenFade( this, Vector(DEATH_FADE_RED, DEATH_FADE_GREEN, DEATH_FADE_BLUE), DEATH_FADE_DURATION, NULL, DEATH_FADE_RENDER_AMT, DEATH_FADE_FLAGS );

	MESSAGE_BEGIN( MSG_ONE, gmsgShowMissionFailed, NULL, edict() );
		WRITE_BYTE( 0 );
	MESSAGE_END();

	pDeadPlayer = CDeadPlayer::Create(m_iHasSuperSuit, pev->origin, Vector(0, pev->angles.y, pev->angles.x), edict());
}
```

## Server/Client Hooks

Both the server and client hooks are established server-side, as that’s the only time they’re needed.

### h_export.cpp

```c++
void DLLEXPORT GiveFnptrsToDll(	enginefuncs_t* pengfuncsFromEngine, globalvars_t *pGlobals )
{
	// This was the only place I could find to put this.
	g_pClient->InitHook();

	memcpy(&g_engfuncs, pengfuncsFromEngine, sizeof(enginefuncs_t));
	gpGlobals = pGlobals;
}
```

## Server/Client Messages
* Message Names CANNOT be longer than 10 characters…