#include "hud.h"
#include "cl_util.h"
#include "interface.h"
#include "VSSound.h"


const char * pszCDAList[] = 
{
	0, 
	0, 
	"../valve/media/Half-Life01.mp3",	// Adrenaline Horror
	"../valve/media/Prospero01.mp3",	// Vague Voices
	"../valve/media/Half-Life12.mp3",	// Klaxon Beat
	"../valve/media/Half-Life07.mp3",	// Space Ocean
	"../valve/media/Half-Life10.mp3",	// Cavern Ambience
	"../valve/media/Suspense01.mp3",	// Apprehensive Short
	"../valve/media/Suspense03.mp3",	// Bass String Short
	"../valve/media/Half-Life09.mp3",	// Hurricane Strings
	"../valve/media/Half-Life02.mp3",	// Diabolical Adrenaline Guitar
	"../valve/media/Half-Life13.mp3",	// Valve Theme [Extended]
	"../valve/media/Half-Life04.mp3",	// Nepal Monastery
	"../valve/media/Half-Life15.mp3",	// Alien Shock
	"../valve/media/Half-Life14.mp3",	// Sirens In The Distance
	"../valve/media/Half-Life16.mp3",	// Nuclear Mission Jam
	"../valve/media/Suspense02.mp3",	// Scared Confusion Short
	"../valve/media/Half-Life03.mp3",	// Drums and Riffs
	"../valve/media/Half-Life08.mp3",	// Hard Technology Rock
	"../valve/media/Prospero02.mp3",	// Steam in the Pipes
	"../valve/media/Half-Life05.mp3",	// Electric Guitar Ambience
	"../valve/media/Prospero04.mp3",	// Dimensionless Deepness
	"../valve/media/Half-Life11.mp3",	// Military Precision
	"../valve/media/Half-Life06.mp3",	// Jungle Drums
	"../valve/media/Prospero03.mp3",	// Traveling Through Limbo
	"../valve/media/Half-Life17.mp3",	// Credits / Closing Theme
	"../valve/media/Prospero05.mp3",	// Threatening Short
	"../valve/media/Suspense05.mp3",	// Dark Piano Short
	"../valve/media/Suspense03.mp3"		// Sharp Fear Short
};

CVSSound::CVSSound()
{
	memset(this, 0, sizeof(CVSSound));
}

int CVSSound::Initialize()
{
	if(bInit)
		return 0;

    result = FMOD::System_Create(&system);
    ERRCHECK(result);

    result = system->getVersion(&version);
    ERRCHECK(result);

    if (version < FMOD_VERSION)
    {
		gEngfuncs.Con_Printf("Fatal Error: You are using an old version of FMOD %08x.  This program requires %08x\n", version, FMOD_VERSION);
		bInit = false;
		return 1;
	}

    result = system->init(VSYS_MAX_CHANNELS, FMOD_INIT_NORMAL, 0);
    ERRCHECK(result);

	gEngfuncs.Con_Printf("FMODEx Version %08x loaded succesfully!\n", version);
	bInit = true;

	DeathInit();

	return 1;
}

int CVSSound::Shutdown()
{
	if(!bInit)
		return 1;

	bool bTemp;

	for(int i = 0; i < VSYS_MAX_CHANNELS; i++)
	{
		if(pChannel[i])
		{
			pChannel[i]->isPlaying(&bTemp);
			if(bTemp)
				pChannel[i]->stop();
		}

		if(pSound[i])
		{
			pSound[i]->release();
			pSound[i] = NULL;
		}
	}

	DeathRelease();

    result = system->close();
    ERRCHECK(result);
    result = system->release();
    ERRCHECK(result);

	return 0;
}

int CVSSound::StopMP3( void )
{
	bool bTemp;
	pChannel[0]->isPlaying(&bTemp);
	if(bTemp)
		pChannel[0]->stop();

	pSound[0]->release();
	return 1;
}

int CVSSound::PlayMP3(const char *pszSong)
{
	bool bTemp;
	pChannel[0]->isPlaying(&bTemp);
	if(bTemp)
		ReleaseChannel(0);

	SetupChannel(0, pszSong, false, FMOD_LOOP_NORMAL);

	return 1;
}

int CVSSound::Play( int iChannel )
{
	Stop(iChannel);

	system->playSound(FMOD_CHANNEL_FREE, pSound[iChannel], false, &pChannel[iChannel]);

	return 1;
}

int CVSSound::CheckCached( const char *pszFile )
{
	if(!pszFile)
		return false;

	for(int i = 0; i < VSYS_MAX_SOUNDS; i++)
	{
		if(List[i].pszName)
		{
			if(strcmp(pszFile, List[i].pszName) == 0)
				return i;
		}
	}

	return -1;
}

int CVSSound::CacheSound( const char *pszFile )
{
	if(!pszFile)
		return -1;

	if(pszFile[0] == '\0')
		return -1;

	int iNameLen = strlen(pszFile)+1;

	if(iNameLen == 0)
		return -1;

	int iPos = -1;								// Set our position to INVALID

	for(int i = 0; i < VSYS_MAX_SOUNDS; i++)	// Loop through all the loaded sounds
	{
		if(!List[i].pszName)					// If our name has not been taken
		{
			iPos = i;							// Set the current position
			break;								// Break;
		}
	}

	char szFullPath[60];

	memset(&szFullPath, 0, 60);

	sprintf(szFullPath, "%s/%s", gEngfuncs.pfnGetGameDirectory(), pszFile);

	int iFileLen = 0;
	FILE *pFile = fopen(szFullPath, "rb");

	if(!pFile)
		return -1;

	fseek(pFile, 0, SEEK_END);
	iFileLen = ftell(pFile);
	fseek(pFile, 0, SEEK_SET);

	if(!iFileLen)
		return -1;

	List[iPos].pBuffer = malloc(iFileLen);
	fread(List[iPos].pBuffer, iFileLen, 1, pFile);

	List[iPos].pszName = new char[iNameLen];
	memset(List[iPos].pszName, 0, iNameLen);
	memcpy(List[iPos].pszName, pszFile, iNameLen-1);

	memset(&List[iPos].exinfo, 0, sizeof(List[iPos].exinfo));
	List[iPos].exinfo.cbsize = sizeof(List[iPos].exinfo);
	List[iPos].exinfo.length = iFileLen;

	return iPos;
		/*char szError[128];
		memset(&szError, 0, 128);
		sprintf((char*)&szError, "FMODEx ERROR: (%d) %s", result, FMOD_ErrorString(result));
		ALERT(at_error,(char*)szError);	
		*/
}

void CVSSound::UncacheSound( int iChannel )
{
	if(List[iChannel].pszName)
	{
		delete[] List[iChannel].pszName;
		List[iChannel].pszName = NULL;
	}

	if(List[iChannel].pBuffer)
	{
		free(List[iChannel].pBuffer);
		List[iChannel].pBuffer = NULL;
	}

	memset(&List[iChannel], 0, sizeof(List[iChannel]));
}

void CVSSound::UncacheAllSounds( void )
{
	for(int i = 0; i < VSYS_MAX_SOUNDS; i++)
	{
		UncacheSound(i);
	}
}

void CVSSound::UncacheSound( const char *pszFile )
{
	if(!pszFile)
		return;

	for(int i = 0; i < VSYS_MAX_SOUNDS; i++)
	{
		if(strcmp(pszFile, List[i].pszName) == 0)
			UncacheSound(i);
	}
}

int CVSSound::SetupChannel( int iChannel, const char *pszSong, bool bPaused, FMOD_MODE mode )
{
	if(pChannel[iChannel])
		ReleaseChannel(iChannel);

	int iPos = CheckCached(pszSong);

	if(!((iPos >= 0) && (iPos < VSYS_MAX_SOUNDS)))
	{
		CONPRINT("VSys: WARNING: File %s was not precached.\n", pszSong);
		iPos = CacheSound(pszSong);

		if(!((iPos >= 0) && (iPos < VSYS_MAX_SOUNDS)))
			return 1;
	}


	system->createSound((const char*)List[iPos].pBuffer, FMOD_HARDWARE | FMOD_OPENMEMORY, &List[iPos].exinfo, &pSound[iChannel]);
	pSound[iChannel]->setMode(mode);
	system->playSound(FMOD_CHANNEL_FREE, pSound[iChannel], bPaused, &pChannel[iChannel]);

	UncacheSound(iPos);

	return 0;
}


int CVSSound::ReleaseChannel(int iChannel)
{
	if((iChannel < VSYS_MAX_CHANNELS) && (iChannel >= 0))
	{
		pChannel[iChannel]->stop();
		pSound[iChannel]->release();

		pChannel[iChannel] = NULL;
		pSound[iChannel] = NULL;

	}
	else if(iChannel == -1)
	{
		for(int i = 0; i < VSYS_MAX_CHANNELS; i++)
		{
			pChannel[i]->stop();
			pSound[i]->release();

			pChannel[i] = NULL;
			pSound[i] = NULL;

		}

		UncacheAllSounds();
	}
	return 0;
}

int CVSSound::Stop(int iChannel)
{
	bool bTemp;
	pChannel[iChannel]->isPlaying(&bTemp);
	
	if(!bTemp)
		return 0;

	pChannel[iChannel]->stop();

	return 0;
}

bool CVSSound::ERRCHECK(FMOD_RESULT result)
{
    if (result != FMOD_OK)
    {
		//char szError[128];
		//memset(&szError, 0, 128);
		//sprintf((char*)&szError, "FMODEx ERROR: (%d) %s", result, FMOD_ErrorString(result));
		//ALERT(at_error,(char*)szError);

        return FALSE;
    }
	return TRUE;
}

void CVSSound::ToggleChannel(int iChannel)
{
	bool bTemp;
	pChannel[iChannel]->getPaused(&bTemp);

	if(bTemp)
		Resume(iChannel);
	else
		Pause(iChannel);
}

void CVSSound::Pause(int iChannel)
{
	if(pChannel[iChannel])
		pChannel[iChannel]->setPaused(true);
}

void CVSSound::Resume(int iChannel)
{
	if(pChannel[iChannel])
		pChannel[iChannel]->setPaused(false);
}

void CVSSound::SetVolume(int iChannel, int iVolume)
{
	float fVol = iVolume/100;
	if(!pChannel[iChannel])
		pChannel[iChannel]->setVolume(fVol);
}

void CVSSound::ParseString(const char *pszString)
{
	if(!pszString)
		return;

	int argc;
	char **argv;
	int iChannel;

	_CommandLineToArgvA((char*)pszString, &argc, &argv);

	if(strcmp(argv[0], "toggle") == 0)
	{
		if(argc < 2)
			return;
		iChannel = atoi(argv[1]);

		if((iChannel < VSYS_MAX_CHANNELS) && (iChannel >= 0))
		{
			ToggleChannel(iChannel);
		}
		else if(iChannel == -1)
		{
			ToggleAll();
		}
	}
	else if(strcmp(argv[0], "pause") == 0)
	{
		if(argc < 2)
			return;

		iChannel = atoi(argv[1]);

		if((iChannel < VSYS_MAX_CHANNELS) && (iChannel >= 0))
		{
			Pause(iChannel);
		}
		else if(iChannel == -1)
		{
			PauseAll();
		}
	}
	else if(strcmp(argv[0], "resume") == 0)
	{
		if(argc < 2)
			return;

		iChannel = atoi(argv[1]);

		if((iChannel < VSYS_MAX_CHANNELS) && (iChannel >= 0))
		{
			Resume(iChannel);
		}
		else if(iChannel == -1)
		{
			ResumeAll();
		}
	}
	else if(strcmp(argv[0], "play") == 0)
	{
		if(argc < 2)
			return;

		iChannel = atoi(argv[1]);
		if((iChannel < VSYS_MAX_CHANNELS) && (iChannel >= 0))
		{
			Play(iChannel);
		}
		else if(iChannel == -1)
		{
			PlayAll();
		}
	}
	else if(strcmp(argv[0], "stop") == 0)
	{
		if(argc < 2)
			return;

		iChannel = atoi(argv[1]);
		if((iChannel < VSYS_MAX_CHANNELS) && (iChannel >= 0))
		{
			Stop(iChannel);
		}
		else if(iChannel = -1)
		{
			StopAll();
		}
	}
	else if(strcmp(argv[0], "setup") == 0)
	{
		if(argc < 5)
			return;

		iChannel = atoi(argv[1]);
		
		bool bPaused = false;

		int iTemp = atoi(argv[2]);

		if(iTemp == 1)
			bPaused = true;

		bool bLooped = false;

		iTemp = atoi(argv[3]);
		
		if(iTemp == 1)
			bLooped = true;

		FMOD_MODE mode;
		if(bLooped)
			mode = FMOD_LOOP_NORMAL;
		else
			mode = FMOD_LOOP_OFF;

		if(strcmp(argv[4], "") == 0)
			return;

		if(!((iChannel < VSYS_MAX_CHANNELS) && (iChannel >= 0)))
			return;

		SetupChannel(iChannel, argv[4], bPaused, mode);
	}
	else if(strcmp(argv[0], "release") == 0)
	{
		if(argc < 2)
			return;

		iChannel = atoi(argv[1]);

		ReleaseChannel(iChannel);
	}
	else if(strcmp(argv[0], "setvol") == 0)
	{
		if(argc < 3)
			return;

		iChannel = atoi(argv[1]);
		int iVolume = atoi(argv[2]);

		if(!(iVolume >= 0))
			return;

		if((iChannel < VSYS_MAX_CHANNELS) && (iChannel >= 0))
		{
			SetVolume(iChannel, iVolume);
		}
		else if(iChannel == -1)
		{
			for(int i = 0; i < VSYS_MAX_CHANNELS; i++)
			{
				if(pChannel[i])
					SetVolume(i, iVolume);
			}
		}
	}
	else if(strcmp(argv[0], "cache") == 0)
	{
		if(argc < 2)
			return;

		CacheSound(argv[1]);
	}
	else if(strcmp(argv[0], "playdeath") == 0)
	{
		DeathPlay();
	}
	else if(strcmp(argv[0], "stopdeath") == 0)
	{
		DeathStop();
	}
	else if(strcmp(argv[0], "playcda") == 0)
	{
		if(argc < 3)
			return;

		int iTrack = atoi(argv[1]);

		if (((iTrack < -1) || (iTrack > 30)))
			return;
		

		int iTemp = atoi(argv[2]);

		bool bLooped = false;
		
		if(iTemp == 1)
			bLooped = true;

		FMOD_MODE mode;
		if(bLooped)
			mode = FMOD_LOOP_NORMAL;
		else
			mode = FMOD_LOOP_OFF;

		SetupChannel(0, pszCDAList[iTrack], false, mode);
	}
	else if(strcmp(argv[0], "stopcda") == 0)
	{
		Stop(0);
	}
}

/*
_CommandLineToArgvA(char* pszCmdLine, int * _argc, char ***argv_buf)

Converts a command line to a argc<->argv pair.

pszCmdLine is a pointer to the command line character array.
_argc is a pointer to an integer that will contain the number of arguments
argv_buf is a pointer to a char** which will contain the argumemt array

A big thanks to Alexander A. Telyatnikov or "Alter" for the majority of this code.
http://alter.org.ua
*/
void CVSSound::_CommandLineToArgvA(char* CmdLine, int* _argc, char ***argv_buf)
{
    char** argv;
    char*  _argv;
    unsigned long   len;
    unsigned long   argc;
    char   a;
    unsigned long   i, j;

    unsigned char  in_QM;
    unsigned char  in_TEXT;
    unsigned char  in_SPACE;

    len = strlen(CmdLine);
    i = ((len+2)/2)*sizeof(void*) + sizeof(void*);
	int iLen = i + (len+2)*sizeof(char);

	*argv_buf = (char**)malloc(iLen);

	//Believe me, I tried EVERY other way, but this was the only one that worked
	argv = *argv_buf;
    //argv = (char**)malloc(iLen);

    _argv = (char*)(((unsigned char*)argv)+i);

    argc = 0;
    argv[argc] = _argv;
    in_QM = 0;
    in_TEXT = 0;
    in_SPACE = 1;
    i = 0;
    j = 0;

    while( a = CmdLine[i] )
	{
        if(in_QM)
		{
            if(a == '\"')
			{
                in_QM = 0;
            }
			else
			{
                _argv[j] = a;
                j++;
            }
        }
		else
		{
            switch(a)
			{
            case '\"':
                in_QM = 1;
                in_TEXT = 1;
                if(in_SPACE)
				{
                    argv[argc] = _argv+j;
                    argc++;
                }
                in_SPACE = 0;
                break;
            case ' ':
            case '\t':
            case '\n':
            case '\r':
                if(in_TEXT)
				{
                    _argv[j] = '\0';
                    j++;
                }
                in_TEXT = 0;
                in_SPACE = 1;
                break;
            default:
                in_TEXT = 1;
                if(in_SPACE)
				{
                    argv[argc] = _argv+j;
                    argc++;
                }
                _argv[j] = a;
                j++;
                in_SPACE = 0;
                break;
            }
        }
        i++;
    }
    _argv[j] = '\0';
    argv[argc] = NULL;

    (*_argc) = argc;

	//*argv_buf = (char**)malloc(iLen);
	//memcpy(*argv_buf, argv, iLen);
	//memcpy(*argv_buf, argv, iLen);

    //return argv;
}

void CVSSound::DeathInit()
{
	FMOD_CREATESOUNDEXINFO exinfo;

	void *pBuffer	= NULL;
	FILE *pFile		= NULL;
	int iLength		= 0;

	char szPath[60];
	memset(&szPath, 0, 60);
	sprintf(szPath, "%s/%s", gEngfuncs.pfnGetGameDirectory(), VSYS_DEATH_TUNE);

	pFile = fopen(szPath, "rb");

	if(!pFile)
		return;

	fseek(pFile, 0, SEEK_END);
	iLength = ftell(pFile);
	fseek(pFile, 0, SEEK_SET);

	if(!iLength)
		return;

	pBuffer = malloc(iLength);
	fread(pBuffer, iLength, 1, pFile);


	memset(&exinfo, 0, sizeof(exinfo));
	exinfo.cbsize = sizeof(exinfo);
	exinfo.length = iLength;

	system->createSound((const char*)pBuffer, FMOD_HARDWARE | FMOD_OPENMEMORY, &exinfo, &pDeathSound);
	pDeathSound->setMode(FMOD_LOOP_NORMAL);
	system->playSound(FMOD_CHANNEL_FREE, pDeathSound, true, &pDeathChannel);
	pDeathChannel->setVolume(2.0);

}

void CVSSound::DeathStop()
{
	bool bTemp;
	pDeathChannel->isPlaying(&bTemp);
	
	if(!bTemp)
		return;

	pDeathChannel->stop();
}

void CVSSound::DeathPlay()
{
	//if(g_pServer->VSys_StartDeathFade)
	//	g_pServer->VSys_StartDeathFade();

	DeathStop();

	system->playSound(FMOD_CHANNEL_FREE, pDeathSound, false, &pDeathChannel);
}

void CVSSound::DeathRelease()
{
	bool bTemp;

	if(pDeathChannel)
	{
		pDeathChannel->isPlaying(&bTemp);
		if(bTemp)
			pDeathChannel->stop();
	}

	if(pDeathSound)
	{
		pDeathSound->release();
		pDeathSound = NULL;
	}
}

void CVSSound::StopAll()
{
	for(int i = 0; i < VSYS_MAX_CHANNELS; i++)
	{
		if(pChannel[i])
			Stop(i);
	}
}

void CVSSound::PauseAll()
{
	for(int i = 0; i < VSYS_MAX_CHANNELS; i++)
	{
		if(pChannel[i])
			Pause(i);
	}
}

void CVSSound::ResumeAll()
{
	for(int i = 0; i < VSYS_MAX_CHANNELS; i++)
	{
		if(pChannel[i])
			Resume(i);
	}
}

void CVSSound::PlayAll()
{
	for(int i = 0; i < VSYS_MAX_CHANNELS; i++)
	{
		if(pChannel[i])
			Play(i);
	}
}

void CVSSound::ToggleAll()
{
	for(int i = 0; i < VSYS_MAX_CHANNELS; i++)
	{
		if(pChannel[i])
			ToggleChannel(i);
	}
}

void CVSSound::SetSpeed(int iChannel)
{
	/*pChannel[iChannel]->addDSP(

	pChannel[iChannel]->stop();

	return 0;*/
}