#ifndef _CFREEMOVEMENT_H_
#define _CFREEMOVEMENT_H_

#include "CBaseMovement.h"

class CFreeMovement : public CBaseMovement
{
	void CL_AdjustAngles(float frametime, float *viewangles);
	void CL_CreateMove(float frametime, struct usercmd_s *cmd, int active);
	void IN_Move(float frametime, usercmd_t *cmd, float *deltaAngles);
	void IN_MouseMove(float frametime, usercmd_t *cmd, float *deltaAngles);
	void IN_JoyMove(float frametime, usercmd_t *cmd, float *deltaAngles);
};
#endif // _CFREEMOVEMENT_H_