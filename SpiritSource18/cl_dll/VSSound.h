#ifndef _VSSOUND_H_
#define _VSSOUND_H_

#include "../common/fmod/fmod.hpp"
#include "../common/fmod/fmod_errors.h"
#pragma comment(lib, "../common/fmod/fmodex_vc.lib")

#include "windows.h"
#include "../common/vsys/vsys_common.h"

typedef struct
{
	char *pszName;
	void *pBuffer;
	FMOD_CREATESOUNDEXINFO exinfo;
} CFile;

#ifdef VALVE_DLL
#define VSYS_API __declspec(dllimport)
#elif CLIENT_DLL
#define VSYS_API __declspec(dllexport)
#else
#define VSYS_API
#endif

class VSYS_API CVSSound
{
public:
	CVSSound();

	// System Init/Destroy functions
	int					Initialize();
	int					Shutdown();

	// Channel Init/Destroy functions
	int					SetupChannel(int iChannel, const char *pszSong, bool bPaused, FMOD_MODE mode);
	int					ReleaseChannel(int iChannel);

	// For compatibility with AJH's MP3 player.
	int					PlayMP3(const char *pszSong);
	int					StopMP3();

	// Play controls
	int					Stop(int iChannel);
	int					Play(int iChannel);
	void				Pause(int iChannel);
	void				Resume(int iChannel);
	void				SetVolume(int iChannel, int iVolume);
	void				ToggleChannel(int iChannel);
	void				SetSpeed(int iChannel);

	void				StopAll();
	void				PlayAll();
	void				PauseAll();
	void				ResumeAll();
	void				ToggleAll();
	
	// String Parser
	void				ParseString(const char *pszString);

	// Special Channel Functions
	void				DeathInit();
	void				DeathStop();
	void				DeathPlay();
	void				DeathRelease();

	// Caching Functions
	int					CheckCached(const char *pszFile);
	int					CacheSound(const char *pszFile);
	void				UncacheSound(int iChannel);
	void				UncacheSound(const char *pszFile);
	void				UncacheAllSounds();

private:


	// Miscellaneous Functions
	bool				ERRCHECK(FMOD_RESULT result);
	void				_CommandLineToArgvA(char* CmdLine, int* _argc, char ***argv_buf);

	// Variables
	bool				bInit;
	bool				bIsDynamic;
	FMOD::System		*system;

	CFile				List[VSYS_MAX_SOUNDS];
	FMOD::Sound			*pSound[VSYS_MAX_CHANNELS];
	FMOD::Channel		*pChannel[VSYS_MAX_CHANNELS];

	FMOD::Channel		*pDeathChannel;
	FMOD::Sound			*pDeathSound;

	FMOD::Channel		*pDynChannel;
	FMOD::Sound			*pDynSound[4];

	FMOD_RESULT			result;
	int					key;
	unsigned int		version;
};

extern CVSSound gMP3;
#endif
