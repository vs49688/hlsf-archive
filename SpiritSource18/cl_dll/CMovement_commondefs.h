#ifndef _CMOVEMENT_COMMONDEFS_H_
#define _CMOVEMENT_COMMONDEFS_H_


// External Functions
extern float CL_KeyState (kbutton_t *key);
extern int CL_ButtonBits( int );
extern void IN_ResetMouse( void );
extern float V_ClampYaw(float oldYaw, float proposedYaw);
extern float V_ClampPitch(float oldPitch, float proposedPitch);
extern Vector V_LimitClampSpeed(Vector& prev, Vector& next, float frametime);
extern void Joy_AdvancedUpdate_f (void);
extern int IN_ReadJoystick (void);
extern int GetMouseActive();

extern "C" 
{
	void DLLEXPORT IN_ActivateMouse( void );
	void DLLEXPORT IN_DeactivateMouse( void );
	void DLLEXPORT IN_MouseEvent (int mstate);
	void DLLEXPORT IN_Accumulate (void);
	void DLLEXPORT IN_ClearStates (void);
	
	float anglemod( float a );
	void QuaternionRotate(float anglesIn[3], float deltaAngles[3], float anglesOut[3]);

}


// External Variables
extern int				g_weaponselect;
extern int				in_impulse;
extern int				in_cancel;
extern int				g_iAlive;
extern int				iMouseInUse;
extern int				g_iVisibleMouse;
extern int				mouse_x, mouse_y, old_mouse_x, old_mouse_y, mx_accum, my_accum;
extern int				joy_avail, joy_advancedinit, joy_haspov;
extern int				joy_id;
extern int *			pMouseActive;

extern DWORD			joy_oldbuttonstate, joy_oldpovstate;
extern DWORD			joy_flags;
extern DWORD			joy_numbuttons;
extern DWORD			dwAxisMap[ JOY_MAX_AXES ];
extern DWORD			dwControlMap[ JOY_MAX_AXES ];

extern PDWORD			pdwRawValue[ JOY_MAX_AXES ];


extern POINT			current_pos;

extern cl_enginefunc_t	gEngfuncs;

extern kbutton_t	in_mlook;
extern kbutton_t	in_klook;
extern kbutton_t	in_jlook;
extern kbutton_t	in_left;
extern kbutton_t	in_right;
extern kbutton_t	in_forward;
extern kbutton_t	in_back;
extern kbutton_t	in_lookup;
extern kbutton_t	in_lookdown;
extern kbutton_t	in_moveleft;
extern kbutton_t	in_moveright;
extern kbutton_t	in_strafe;
extern kbutton_t	in_speed;
extern kbutton_t	in_use;
extern kbutton_t	in_jump;
extern kbutton_t	in_attack;
extern kbutton_t	in_attack2;
extern kbutton_t	in_up;
extern kbutton_t	in_down;
extern kbutton_t	in_duck;
extern kbutton_t	in_reload;
extern kbutton_t	in_alt1;
extern kbutton_t	in_score;
extern kbutton_t	in_break;
extern kbutton_t	in_graph;  // Display the netgraph
extern kbutton_t	in_customhud;	//AJH custom hud
extern kbutton_t	in_briefing;	//AJH show map briefing

extern cvar_t	*m_pitch;
extern cvar_t	*m_yaw;
extern cvar_t	*m_forward;
extern cvar_t	*m_side;

extern cvar_t	*lookstrafe;
extern cvar_t	*lookspring;
extern cvar_t	*cl_pitchup;
extern cvar_t	*cl_pitchdown;
extern cvar_t	*cl_upspeed;
extern cvar_t	*cl_forwardspeed;
extern cvar_t	*cl_backspeed;
extern cvar_t	*cl_sidespeed;
extern cvar_t	*cl_movespeedkey;
extern cvar_t	*cl_yawspeed;
extern cvar_t	*cl_pitchspeed;
extern cvar_t	*cl_anglespeedkey;
extern cvar_t	*cl_vsmoothing;

extern cvar_t	*in_joystick;
extern cvar_t	*joy_name;
extern cvar_t	*joy_advanced;
extern cvar_t	*joy_advaxisx;
extern cvar_t	*joy_advaxisy;
extern cvar_t	*joy_advaxisz;
extern cvar_t	*joy_advaxisr;
extern cvar_t	*joy_advaxisu;
extern cvar_t	*joy_advaxisv;
extern cvar_t	*joy_forwardthreshold;
extern cvar_t	*joy_sidethreshold;
extern cvar_t	*joy_pitchthreshold;
extern cvar_t	*joy_yawthreshold;
extern cvar_t	*joy_forwardsensitivity;
extern cvar_t	*joy_sidesensitivity;
extern cvar_t	*joy_pitchsensitivity;
extern cvar_t	*joy_yawsensitivity;
extern cvar_t	*joy_wwhack1;
extern cvar_t	*joy_wwhack2;

extern cvar_t		*m_filter;
extern cvar_t		*sensitivity;


#endif // _CMOVEMENT_COMMONDEFS_H_