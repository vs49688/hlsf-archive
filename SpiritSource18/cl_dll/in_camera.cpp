//========= Copyright © 1996-2002, Valve LLC, All rights reserved. ============
// New clipping style camera - original idea by Xwider
// Purpose: 
//
// $NoKeywords: $
//=============================================================================

#include "hud.h"
#include "cl_util.h"
#include "camera.h"
#include "in_defs.h"
#include "pmtrace.h"
#include "event_api.h"
#include "pm_defs.h" 

extern "C" 
{
	void DLLEXPORT CAM_Think( void );
	int DLLEXPORT CL_IsThirdPerson( void );
	void DLLEXPORT CL_CameraOffset( float *ofs );
}

int iMouseInUse = 0;

#define CAM_MIN_DIST 4.0
#define CAM_MAX_DIST 48.0
extern Vector v_angles;
void DLLEXPORT CAM_Think( void )
{
    /*if( !cam_thirdperson )
        return;

    float maxDist = cam_idealdist->value;
    float DistFromWall;
    pmtrace_t tr;
    vec3_t origin, camForward;

    if (maxDist > CAM_MAX_DIST)
        maxDist = CAM_MAX_DIST; //cap
   
    DistFromWall = maxDist;
    cam_ofs[ 0 ] = v_angles[0];
    cam_ofs[ 1 ] = v_angles[1];
    cam_ofs[ 2 ] = v_angles[2];

    // Test camera position
//get the local player pointer
    cl_entity_t *player = gEngfuncs.GetLocalPlayer();      

	//static int frames = 0;      //static frame counter
	//if (frames < 5) frames++;    //cap at 5
	if (player && (gEngfuncs.GetMaxClients() > 0) && (gViewPort && gEngfuncs.pEventAPI) && frames >= 5)
    {
        origin = player->origin;
        AngleVectors( cam_ofs, camForward, NULL, NULL );        //get the forward vector
        gEngfuncs.pEventAPI->EV_SetTraceHull(2);         //use duck hull for traces
        gEngfuncs.pEventAPI->EV_SetSolidPlayers(player->index - 1);
        gEngfuncs.pEventAPI->EV_PlayerTrace( origin, origin - (camForward * (maxDist + 8)), PM_STUDIO_BOX, -1, &tr );    //trace to maxDist +8 (leway for the bounding boxes)
       
        if ( tr.fraction < 1.0 )         //if a wall or object was hit..
        {
            DistFromWall = maxDist * tr.fraction;   //calculate the final distance
            if( DistFromWall < CAM_MIN_DIST )    //if it's too close
                DistFromWall = CAM_MIN_DIST;    //push it to min dist
            player->curstate.renderamt = 255 * tr.fraction;      //calculate the players render amount based on distance from camera
            player->curstate.rendermode = kRenderTransColor; //use transcolor for transparency
        }
        cam_ofs[ 2 ] = DistFromWall;    //pass  the distance off for the camera; if you want to keep the roll (z axis) then make DistFromWall a global or something (like v_angles) and read it in view.cpp instead of cam_ofs[2])
    }*/
}
void CAM_Init( void )
{
}

int DLLEXPORT CL_IsThirdPerson( void )
{
	return (gHUD.m_iCameraMode ? 1 : 0) || (g_iUser1 && (g_iUser2 == gEngfuncs.GetLocalPlayer()->index) );
}

void DLLEXPORT CL_CameraOffset( float *ofs )
{
	VectorCopy( vec3_origin, ofs );
}
