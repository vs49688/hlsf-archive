#ifndef _CBASEMOVEMENT_H_
#define _CBASEMOVEMENT_H_

#define JOY_ABSOLUTE_AXIS	0x00000000		// control like a joystick
#define JOY_RELATIVE_AXIS	0x00000010		// control like a mouse, spinner, trackball
#define	JOY_MAX_AXES		6				// X, Y, Z, R, U, V
#define JOY_AXIS_X			0
#define JOY_AXIS_Y			1
#define JOY_AXIS_Z			2
#define JOY_AXIS_R			3
#define JOY_AXIS_U			4
#define JOY_AXIS_V			5

enum _ControlList
{
	AxisNada = 0,
	AxisForward,
	AxisLook,
	AxisSide,
	AxisTurn
};

enum _INPUTTYPES
{
	INPUT_NORMAL = 0,
	INPUT_FREE
};

#include "util_vector.h"
#include "usercmd.h"

class CBaseMovement
{
public:
	virtual void CL_CreateMove(float frametime, struct usercmd_s *cmd, int active);
	virtual void CL_AdjustAngles(float frametime, float *viewangles);
	virtual void IN_Move(float frametime, usercmd_t *cmd);
	virtual void IN_MouseMove(float frametime, usercmd_t *cmd);
	virtual void IN_JoyMove(float frametime, usercmd_t *cmd);

};

#endif // _CBASEMOVEMENT_H_