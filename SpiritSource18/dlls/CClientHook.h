#ifndef _CCLIENTHOOK_H_
#define _CCLIENTHOOK_H_

#include <windows.h>

//typedef CVSSound*(*VGPTR)(void);
typedef int(*I_VSys_I)(int);			// int foo(int)
typedef void(*V_VSys_I)(int);			// void foo(int)
typedef void(*V_VSys_II)(int,int);		// void foo(int, int)
typedef int(*I_VSys_CH)(const char*);	// int foo(const char*)
typedef void(*V_VSys_CH)(const char*);	// int foo(const char*)
typedef void(*V_VSys_V)(void);			// void foo(void)
typedef void(*SETINSTANCE)(HINSTANCE);	// void foo(HINSTANCE)

typedef int(*SETUPCHANNEL)(int,const char*,bool,unsigned int);
//typedef void(*CLIENTDEATH)(void*);


class CClientHook
{
public:
	CClientHook();
	int				InitHook();
	bool			GetHookState();

	I_VSys_CH		VSys_PrecacheSound;	
	SETUPCHANNEL	VSys_SetupChannel;
	//CLIENTDEATH		VSys_ClientDeath;
	
	//SETINSTANCE	VSys_SetServerInstance;
	//V_VSys_V	VSys_InitServerHook;

private:

	bool			bHookState;	
	HMODULE			hModule;

	//VGPTR		lpfnGetVSysPtr;
	//I_VSys_I	lpfnSetupChannel;
	//I_VSys_I	lpfnStop;
	//I_VSys_I	lpfnPlay;
	//V_VSys_I	lpfnPause;
	//V_VSys_I	lpfnResume;
	//V_VSys_I	lpfnToggleChannel;
	//V_VSys_II	lpfnSetVolume;
	//I_VSys_CH	lpfnCheckCached;
	//I_VSys_CH	lpfnCacheSound;
	//V_VSys_I	lpfnUncacheSound;
	//V_VSys_CH	lpfnUncacheSound;
	//V_VSys_V	lpfnUncacheAllSounds;

};

extern CClientHook *g_pClient;
#endif // _CCLIENTHOOK_H_