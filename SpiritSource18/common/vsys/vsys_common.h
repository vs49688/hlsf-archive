#ifndef _VSYS_COMMON_H_
#define _VSYS_COMMON_H_

#define VSYS_MAX_CHANNELS			32
#define VSYS_MAX_SOUNDS				32

#define VSYS_NUM_SPECIAL_CHANNELS	2

#define VSYS_CHANNELID_NEXT			-2
#define VSYS_CHANNELID_ALL			-1
#define VSYS_CHANNELID_MUSIC		0
#define VSYS_CHANNELID_DYNAMIC		1

#define VSYS_DEATH_TUNE			"sound/music/death.flac"

#define CONPRINT (gEngfuncs.Con_Printf)

#endif // _VSYS_COMMON_H_