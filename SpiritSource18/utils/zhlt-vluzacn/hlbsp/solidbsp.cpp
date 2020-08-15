#pragma warning(disable: 4018) // '<' : signed/unsigned mismatch

#include "bsp5.h"

//  FaceSide
//  ChooseMidPlaneFromList
//  ChoosePlaneFromList
//  SelectPartition

//  CalcSurfaceInfo
//  DivideSurface
//  SplitNodeSurfaces
//  RankForContents
//  ContentsForRank

//  FreeLeafSurfs
//  LinkLeafFaces
//  MakeNodePortal
//  SplitNodePortals
//  CalcNodeBounds
//  CopyFacesToNode
//  BuildBspTree_r
//  SolidBSP

//  Each node or leaf will have a set of portals that completely enclose
//  the volume of the node and pass into an adjacent node.

int             g_maxnode_size = DEFAULT_MAXNODE_SIZE;

static bool g_reportProgress = false;
static int  g_numProcessed = 0;
static int  g_numReported = 0;

static void ResetStatus(bool report_progress)
{
	g_reportProgress = report_progress;
	g_numProcessed = g_numReported = 0;
}

static void UpdateStatus(void)
{
	if(g_reportProgress)
	{
		++g_numProcessed;
		if((g_numProcessed / 500) > g_numReported)
		{
			g_numReported = (g_numProcessed / 500);
			Log("%d...",g_numProcessed);
		}
	}
}	

// =====================================================================================
//  FaceSide
//      For BSP hueristic
// =====================================================================================
static int      FaceSide(face_t* in, const dplane_t* const split
#ifdef HLBSP_AVOIDEPSILONSPLIT
						 , double *epsilonsplit = NULL
#endif
						 )
{
#ifdef HLBSP_AVOIDEPSILONSPLIT
	const vec_t		epsilonmin = 0.002, epsilonmax = 0.2;
	vec_t			d_front, d_back;
#else
    int             frontcount, backcount;
#endif
    vec_t           dot;
    int             i;
    vec_t*          p;

#ifdef HLBSP_AVOIDEPSILONSPLIT
	d_front = d_back = 0;
#else
    frontcount = backcount = 0;
#endif

    // axial planes are fast
    if (split->type <= last_axial)
    {
        vec_t           splitGtEp = split->dist + ON_EPSILON;   // Invariant moved out of loop
        vec_t           splitLtEp = split->dist - ON_EPSILON;   // Invariant moved out of loop

        for (i = 0, p = in->pts[0] + split->type; i < in->numpoints; i++, p += 3)
        {
#ifdef HLBSP_AVOIDEPSILONSPLIT
			dot = *p - split->dist;
			if (dot > d_front)
				d_front = dot;
			if (dot < d_back)
				d_back = dot;
#else
            if (*p > splitGtEp)
            {
                if (backcount)
                {
                    return SIDE_ON;
                }
                frontcount = 1;
            }
            else if (*p < splitLtEp)
            {
                if (frontcount)
                {
                    return SIDE_ON;
                }
                backcount = 1;
            }
#endif
        }
    }
    else
    {
        // sloping planes take longer
        for (i = 0, p = in->pts[0]; i < in->numpoints; i++, p += 3)
        {
            dot = DotProduct(p, split->normal);
            dot -= split->dist;
#ifdef HLBSP_AVOIDEPSILONSPLIT
			if (dot > d_front)
				d_front = dot;
			if (dot < d_back)
				d_back = dot;
#else
            if (dot > ON_EPSILON)
            {
                if (backcount)
                {
                    return SIDE_ON;
                }
                frontcount = 1;
            }
            else if (dot < -ON_EPSILON)
            {
                if (frontcount)
                {
                    return SIDE_ON;
                }
                backcount = 1;
            }
#endif
        }
    }
#ifdef HLBSP_AVOIDEPSILONSPLIT
	if (d_front <= ON_EPSILON)
	{
		if (d_front > epsilonmin || d_back > -epsilonmax)
		{
			if (epsilonsplit)
				(*epsilonsplit)++;
		}
		return SIDE_BACK;
	}
	if (d_back >= -ON_EPSILON)
	{
		if (d_back < -epsilonmin || d_front < epsilonmax)
		{
			if (epsilonsplit)
				(*epsilonsplit)++;
		}
		return SIDE_FRONT;
	}
	if (d_front < epsilonmax || d_back > -epsilonmax)
	{
		if (epsilonsplit)
			(*epsilonsplit)++;
	}
	return SIDE_ON;
#else

    if (!frontcount)
    {
        return SIDE_BACK;
    }
    if (!backcount)
    {
        return SIDE_FRONT;
    }

    return SIDE_ON;
#endif
}

// =====================================================================================
//  ChooseMidPlaneFromList
//      When there are a huge number of planes, just choose one closest
//      to the middle.
// =====================================================================================
static surface_t* ChooseMidPlaneFromList(surface_t* surfaces, const vec3_t mins, const vec3_t maxs
#ifdef ZHLT_DETAILBRUSH
										 , int detaillevel
#endif
										 )
{
    int             j, l;
    surface_t*      p;
    surface_t*      bestsurface;
    vec_t           bestvalue;
    vec_t           value;
    vec_t           dist;
    dplane_t*       plane;

    //
    // pick the plane that splits the least
    //
#ifdef ZHLT_LARGERANGE
	bestvalue = 6.0f * BOGUS_RANGE * BOGUS_RANGE;
#else
    bestvalue = 6 * 8192 * 8192;
#endif
    bestsurface = NULL;

    for (p = surfaces; p; p = p->next)
    {
        if (p->onnode)
        {
            continue;
        }
#ifdef ZHLT_DETAILBRUSH
		if (p->detaillevel != detaillevel)
		{
			continue;
		}
#endif

        plane = &g_dplanes[p->planenum];

        // check for axis aligned surfaces
        l = plane->type;
        if (l > last_axial)
        {
            continue;
        }

        //
        // calculate the split metric along axis l, smaller values are better
        //
        value = 0;

        dist = plane->dist * plane->normal[l];
#ifdef HLBSP_MAXNODESIZE_SKYBOX
		if (maxs[l] - dist < ON_EPSILON || dist - mins[l] < ON_EPSILON)
			continue;
#endif
#ifdef HLBSP_ChooseMidPlane_FIX
		if (maxs[l] - dist < g_maxnode_size/2.0 - ON_EPSILON || dist - mins[l] < g_maxnode_size/2.0 - ON_EPSILON)
			continue;
#endif
        for (j = 0; j < 3; j++)
        {
            if (j == l)
            {
                value += (maxs[l] - dist) * (maxs[l] - dist);
                value += (dist - mins[l]) * (dist - mins[l]);
            }
            else
            {
                value += 2 * (maxs[j] - mins[j]) * (maxs[j] - mins[j]);
            }
        }

        if (value > bestvalue)
        {
            continue;
        }

        //
        // currently the best!
        //
        bestvalue = value;
        bestsurface = p;
    }

    if (!bestsurface)
    {
#ifdef HLBSP_ChooseMidPlane_FIX
		return NULL;
#else
        for (p = surfaces; p; p = p->next)
        {
            if (!p->onnode)
            {
#ifdef ZHLT_DETAILBRUSH
				if (p->detaillevel != detaillevel)
				{
					continue;
				}
#endif
                return p;                                  // first valid surface
            }
        }
        Error("ChooseMidPlaneFromList: no valid planes");
#endif
    }

    return bestsurface;
}

// =====================================================================================
//  ChoosePlaneFromList
//      Choose the plane that splits the least faces
// =====================================================================================
#ifdef HLBSP_ChoosePlane_VL
static surface_t* ChoosePlaneFromList(surface_t* surfaces, const vec3_t mins, const vec3_t maxs
#ifdef ZHLT_DETAILBRUSH
									  // mins and maxs are invalid when detaillevel > 0
									  , int detaillevel
#endif
									  )
{
	surface_t*      p;
	surface_t*      p2;
	surface_t*      bestsurface;
	vec_t           bestvalue;
	vec_t           value;
	dplane_t*       plane;
	face_t*         f;

	double			avesplit;
	double			planecount;
	{
		planecount = 0;
		double totalsplit = 0;
		for (p = surfaces; p; p = p->next)
		{
			if (p->onnode)
			{
				continue;
			}
#ifdef ZHLT_DETAILBRUSH
			if (p->detaillevel != detaillevel)
			{
				continue;
			}
#endif
			planecount++;
			plane = &g_dplanes[p->planenum];
			for (p2 = surfaces; p2; p2 = p2->next)
			{
				if (p2->onnode || p2 == p)
				{
					continue;
				}
				for (f = p2->faces; f; f = f->next)
				{
#ifdef HLCSG_HLBSP_SOLIDHINT
					if (f->facestyle == face_discardable)
					{
						continue;
					}
#endif
					if (FaceSide (f, plane) == SIDE_ON)
					{
						totalsplit++;
					}
				}
			}
		}
		avesplit = (double)totalsplit / (double)planecount;
	}
	//
	// pick the plane that splits the least
	//
	bestvalue = 9e30;
	bestsurface = NULL;

	for (p = surfaces; p; p = p->next)
	{
		if (p->onnode)
		{
			continue;
		}
#ifdef ZHLT_DETAILBRUSH
		if (p->detaillevel != detaillevel)
		{
			continue;
		}
#endif

		double crosscount = 0;
		double frontcount = 0;
		double backcount = 0;
		double coplanarcount = 0;
		double totalcount = 0;
#ifdef HLBSP_AVOIDEPSILONSPLIT
		double epsilonsplit = 0;
#endif

		plane = &g_dplanes[p->planenum];

		for (p2 = surfaces; p2; p2 = p2->next)
		{
			if (p2->onnode)
			{
				for (f = p2->faces; f; f = f->next)
				{
#ifdef HLCSG_HLBSP_SOLIDHINT
					if (f->facestyle == face_discardable)
					{
						continue;
					}
#endif
					totalcount++;
				}
			}
			else if (p2 == p)
			{
				for (f = p2->faces; f; f = f->next)
				{
#ifdef HLCSG_HLBSP_SOLIDHINT
					if (f->facestyle == face_discardable)
					{
						continue;
					}
#endif
					totalcount++;
					coplanarcount++;
				}
			}
			else
			{
				for (f = p2->faces; f; f = f->next)
				{
#ifdef HLCSG_HLBSP_SOLIDHINT
					if (f->facestyle == face_discardable)
					{
#ifdef HLBSP_AVOIDEPSILONSPLIT
						FaceSide (f, plane, &epsilonsplit);
#endif
						continue;
					}
#endif
					totalcount++;
					switch (FaceSide(f, plane
#ifdef HLBSP_AVOIDEPSILONSPLIT
						, &epsilonsplit
#endif
						))
					{
					case SIDE_FRONT:
						frontcount++;
						break;
					case SIDE_BACK:
						backcount++;
						break;
					case SIDE_ON:
						crosscount++;
						break;
					}
				}
				for (f = p->faces; f; f = f->next)
				{
#ifdef HLCSG_HLBSP_SOLIDHINT
					if (f->facestyle == face_discardable)
					{
						continue;
					}
#endif
				}
			}
		}

		value = crosscount - sqrt (coplanarcount); // Not optimized. --vluzacn
#ifdef HLCSG_HLBSP_SOLIDHINT
		if (coplanarcount == 0)
		{
			crosscount += 1;
		}
#endif
#ifdef HLBSP_BALANCE
		// This is the most efficient code among what I have ever tested:
		// (1) BSP file is small, despite possibility of slowing down vis and rad (but still faster than the original non BSP balancing method).
		// (2) Factors need not adjust across various maps.
		double frac = (coplanarcount / 2 + crosscount / 2 + frontcount) / (coplanarcount + frontcount + backcount + crosscount);
		double ent = (0.0001 < frac && frac < 0.9999)? (- frac * log (frac) / log (2.0) - (1 - frac) * log (1 - frac) / log (2.0)): 0.0; // the formula tends to 0 when frac=0,1
		value += crosscount * avesplit * (1 - ent);
#endif
#ifdef HLBSP_AVOIDEPSILONSPLIT
		value += epsilonsplit * 10000;
#endif

		if (value < bestvalue)
		{
			//
			// currently the best!
			//
			bestvalue = value;
			bestsurface = p;
		}
	}

	if (!bestsurface)
		Error("ChoosePlaneFromList: no valid planes");
	return bestsurface;
}
#else
static surface_t* ChoosePlaneFromList(surface_t* surfaces, const vec3_t mins, const vec3_t maxs
#ifdef ZHLT_DETAILBRUSH
									  // mins and maxs are invalid when detaillevel > 0
									  , int detaillevel
#endif
									  )
{
    int             j;
    int             k;
    int             l;
    surface_t*      p;
    surface_t*      p2;
    surface_t*      bestsurface;
    vec_t           bestvalue;
    vec_t           bestdistribution;
    vec_t           value;
    vec_t           dist;
    dplane_t*       plane;
    face_t*         f;

    //
    // pick the plane that splits the least
    //
#define UNDESIREABLE_HINT_FACTOR 10000
#define WORST_VALUE 100000000
    bestvalue = WORST_VALUE;
    bestsurface = NULL;
    bestdistribution = 9e30;

    for (p = surfaces; p; p = p->next)
    {
        if (p->onnode)
        {
            continue;
        }
#ifdef ZHLT_DETAILBRUSH
		if (p->detaillevel != detaillevel)
		{
			continue;
		}
#endif

#ifdef ZHLT_DETAIL
        if (g_bDetailBrushes)
        {
            // AJM: cycle though all faces, and make sure none of them are detail
            // if any of them are, this surface isnt to cause a bsp split
            for (face_t* f = p->faces; f; f = f->next)
            {
                if (f->contents == CONTENTS_DETAIL)
                {
                    //Log("ChoosePlaneFromList::got a detial surface, skipping...\n");
                    continue; // wrong. --vluzacn
                }
            }
        }
#endif
        
        plane = &g_dplanes[p->planenum];
        k = 0;

        for (p2 = surfaces; p2; p2 = p2->next)
        {
            if (p2 == p)
            {
                continue;
            }
            if (p2->onnode)
            {
                continue;
            }

            for (f = p2->faces; f; f = f->next)
            {
                // Give this face (a hint brush fragment) a large 'undesireable' value, only split when we have to)
                if (f->facestyle == face_hint)
                {
                    k += UNDESIREABLE_HINT_FACTOR;
                    hlassert(k < WORST_VALUE);
                    if (k >= WORST_VALUE)
                    {
                        Warning("::ChoosePlaneFromList() surface fragmentation undesireability exceeded WORST_VALUE");
                        k = WORST_VALUE - 1;
                    }
                }
                if (FaceSide(f, plane) == SIDE_ON)
                {
                    k++;
                    if (k >= bestvalue)
                    {
                        break;
                    }
                }

            }
            if (k > bestvalue)
            {
                break;
            }
        }

        if (k > bestvalue)
        {
            continue;
        }

        // if equal numbers, axial planes win, then decide on spatial subdivision

        if (k < bestvalue || (k == bestvalue && (plane->type <= last_axial)))
        {
            // check for axis aligned surfaces
            l = plane->type;

            if (l <= last_axial)
            {                                              // axial aligned                                                
                //
                // calculate the split metric along axis l
                //
                value = 0;

                for (j = 0; j < 3; j++)
                {
                    if (j == l)
                    {
                        dist = plane->dist * plane->normal[l];
                        value += (maxs[l] - dist) * (maxs[l] - dist);
                        value += (dist - mins[l]) * (dist - mins[l]);
                    }
                    else
                    {
                        value += 2 * (maxs[j] - mins[j]) * (maxs[j] - mins[j]);
                    }
                }

                if (value > bestdistribution && k == bestvalue)
                {
                    continue;
                }
                bestdistribution = value;
            }
            //
            // currently the best!
            //
            bestvalue = k;
            bestsurface = p;
        }
    }

    return bestsurface;
}
#endif

// =====================================================================================
//  SelectPartition
//      Selects a surface from a linked list of surfaces to split the group on
//      returns NULL if the surface list can not be divided any more (a leaf)
// =====================================================================================
#ifdef ZHLT_DETAILBRUSH
int CalcSplitDetaillevel (const node_t *node)
{
	int bestdetaillevel = -1;
	surface_t *s;
	face_t *f;
	for (s = node->surfaces; s; s = s->next)
	{
		if (s->onnode)
		{
			continue;
		}
#ifdef HLCSG_HLBSP_SOLIDHINT
		for (f = s->faces; f; f = f->next)
		{
			if (f->facestyle == face_discardable)
			{
				continue;
			}
			if (bestdetaillevel == -1 || f->detaillevel < bestdetaillevel)
			{
				bestdetaillevel = f->detaillevel;
			}
		}
#else
		if (bestdetaillevel == -1 || s->detaillevel < bestdetaillevel)
		{
			bestdetaillevel = s->detaillevel;
		}
#endif
	}
	return bestdetaillevel;
}
#endif
static surface_t* SelectPartition(surface_t* surfaces, const node_t* const node, const bool usemidsplit
#ifdef ZHLT_DETAILBRUSH
								  , int splitdetaillevel
#endif
#ifdef HLBSP_MAXNODESIZE_SKYBOX
								  , vec3_t validmins, vec3_t validmaxs
#endif
								  )
{
#ifdef ZHLT_DETAILBRUSH
	if (splitdetaillevel == -1)
	{
		return NULL;
	}
	// now we MUST choose a surface of this detail level
#else
    int             i;
    surface_t*      p;
    surface_t*      bestsurface;

    //
    // count surface choices
    //
    i = 0;
    bestsurface = NULL;
    for (p = surfaces; p; p = p->next)
    {
        if (!p->onnode)
        {
#ifdef ZHLT_DETAIL
            if (g_bDetailBrushes)
            {
                // AJM: cycle though all faces, and make sure none of them are detail
                // if any of them are, this surface isnt to cause a bsp split
                for (face_t* f = p->faces; f; f = f->next)
                {
                    if (f->contents == CONTENTS_DETAIL)
                    {
                        //Log("SelectPartition::got a detial surface, skipping...\n");
                        continue;
                    }
                }
            }
#endif
#ifdef HLCSG_HLBSP_SOLIDHINT
			for (face_t *f = p->faces; f; f = f->next)
			{
#ifdef ZHLT_DETAILBRUSH
				if (f->detaillevel != splitdetaillevel)
				{
					continue;
				}
#endif
				if (f->facestyle != face_discardable)
				{
					break;
				}
			}
			if (!f)
			{
				continue; // this surface is discardable
				// if all surfaces are discardable, this will become a leaf node
			}
#endif
#ifdef ZHLT_DETAILBRUSH
			if (p->detaillevel != splitdetaillevel)
			{
				continue;
			}
#endif
            i++;
            bestsurface = p;
        }
    }

    if (i == 0)
    {
        return NULL;                                       // this is a leafnode
    }

#ifndef HLCSG_HLBSP_SOLIDHINT // although there is only one undiscardable surface, maybe discardable surfaces should split first.
    if (i == 1)
    {
        return bestsurface;                                // this is a final split
    }
#endif
#endif

#ifdef HLBSP_ChooseMidPlane_FIX
	if (usemidsplit)
	{
		surface_t *s = ChooseMidPlaneFromList(surfaces, 
#ifdef HLBSP_MAXNODESIZE_SKYBOX
			validmins, validmaxs
#else
			node->mins, node->maxs
#endif
#ifdef ZHLT_DETAILBRUSH
			, splitdetaillevel
#endif
			);
		if (s != NULL)
			return s;
	}
	return ChoosePlaneFromList(surfaces, node->mins, node->maxs
#ifdef ZHLT_DETAILBRUSH
		, splitdetaillevel
#endif
		);
#else
    if (usemidsplit)
    {
        // do fast way for clipping hull
        return ChooseMidPlaneFromList(surfaces, 
#ifdef HLBSP_MAXNODESIZE_SKYBOX
			validmins, validmaxs
#else
			node->mins, node->maxs
#endif
#ifdef ZHLT_DETAILBRUSH
			, splitdetaillevel
#endif
			);
    }
    else
    {
        // do slow way to save poly splits for drawing hull
        return ChoosePlaneFromList(surfaces, node->mins, node->maxs
#ifdef ZHLT_DETAILBRUSH
			, splitdetaillevel
#endif
			);
    }
#endif
}

// =====================================================================================
//  CalcSurfaceInfo
//      Calculates the bounding box
// =====================================================================================
static void     CalcSurfaceInfo(surface_t* surf)
{
    int             i;
    int             j;
    face_t*         f;

    hlassume(surf->faces != NULL, assume_ValidPointer);    // "CalcSurfaceInfo() surface without a face"

    //
    // calculate a bounding box
    //
    for (i = 0; i < 3; i++)
    {
        surf->mins[i] = 99999;
        surf->maxs[i] = -99999;
    }

#ifdef ZHLT_DETAILBRUSH
	surf->detaillevel = -1;
#endif
    for (f = surf->faces; f; f = f->next)
    {
        if (f->contents >= 0)
        {
            Error("Bad contents");
        }
        for (i = 0; i < f->numpoints; i++)
        {
            for (j = 0; j < 3; j++)
            {
                if (f->pts[i][j] < surf->mins[j])
                {
                    surf->mins[j] = f->pts[i][j];
                }
                if (f->pts[i][j] > surf->maxs[j])
                {
                    surf->maxs[j] = f->pts[i][j];
                }
            }
        }
#ifdef ZHLT_DETAILBRUSH
		if (surf->detaillevel == -1 || f->detaillevel < surf->detaillevel)
		{
			surf->detaillevel = f->detaillevel;
		}
#endif
    }
}
#ifdef ZHLT_DETAILBRUSH
#ifdef HLCSG_HLBSP_SOLIDHINT
void FixDetaillevelForDiscardable (node_t *node, int detaillevel)
{
	// when we move on to the next detaillevel, some discardable faces of previous detail level remain not on node (because they are discardable). remove them now
	surface_t *s, **psnext;
	face_t *f, **pfnext;
	for (psnext = &node->surfaces; s = *psnext; )
	{
		if (s->onnode)
		{
			psnext = &s->next;
			continue;
		}
		hlassume (s->faces, assume_ValidPointer);
		for (pfnext = &s->faces; f = *pfnext; )
		{
			if (detaillevel == -1 || f->detaillevel < detaillevel)
			{
				*pfnext = f->next;
				FreeFace (f);
			}
			else
			{
				pfnext = &f->next;
			}
		}
		if (!s->faces)
		{
			*psnext = s->next;
			FreeSurface (s);
		}
		else
		{
			psnext = &s->next;
			CalcSurfaceInfo (s);
			hlassume (!(detaillevel == -1 || s->detaillevel < detaillevel), assume_first);
		}
	}
}
#endif
#endif

// =====================================================================================
//  DivideSurface
// =====================================================================================
static void     DivideSurface(surface_t* in, const dplane_t* const split, surface_t** front, surface_t** back)
{
    face_t*         facet;
    face_t*         next;
    face_t*         frontlist;
    face_t*         backlist;
    face_t*         frontfrag;
    face_t*         backfrag;
    surface_t*      news;
    dplane_t*       inplane;

    inplane = &g_dplanes[in->planenum];

    // parallel case is easy

    if (inplane->normal[0] == split->normal[0]
     && inplane->normal[1] == split->normal[1]
     && inplane->normal[2] == split->normal[2])
    {
        if (inplane->dist > split->dist)
        {
            *front = in;
            *back = NULL;
        }
        else if (inplane->dist < split->dist)
        {
            *front = NULL;
            *back = in;
        }
        else
        {                                                  // split the surface into front and back
            frontlist = NULL;
            backlist = NULL;
            for (facet = in->faces; facet; facet = next)
            {
                next = facet->next;
                if (facet->planenum & 1)
                {
                    facet->next = backlist;
                    backlist = facet;
                }
                else
                {
                    facet->next = frontlist;
                    frontlist = facet;
                }
            }
            goto makesurfs;
        }
        return;
    }

    // do a real split.  may still end up entirely on one side
    // OPTIMIZE: use bounding box for fast test
    frontlist = NULL;
    backlist = NULL;

    for (facet = in->faces; facet; facet = next)
    {
        next = facet->next;
        SplitFace(facet, split, &frontfrag, &backfrag);
        if (frontfrag)
        {
            frontfrag->next = frontlist;
            frontlist = frontfrag;
        }
        if (backfrag)
        {
            backfrag->next = backlist;
            backlist = backfrag;
        }
    }

    // if nothing actually got split, just move the in plane
makesurfs:
#ifdef HLBSP_REMOVECOLINEARPOINTS
	if (frontlist == NULL && backlist == NULL)
	{
		*front = NULL;
		*back = NULL;
		return;
	}
#endif
    if (frontlist == NULL)
    {
        *front = NULL;
        *back = in;
        in->faces = backlist;
        return;
    }

    if (backlist == NULL)
    {
        *front = in;
        *back = NULL;
        in->faces = frontlist;
        return;
    }

    // stuff got split, so allocate one new surface and reuse in
    news = AllocSurface();
    *news = *in;
    news->faces = backlist;
    *back = news;

    in->faces = frontlist;
    *front = in;

    // recalc bboxes and flags
    CalcSurfaceInfo(news);
    CalcSurfaceInfo(in);
}

// =====================================================================================
//  SplitNodeSurfaces
// =====================================================================================
static void     SplitNodeSurfaces(surface_t* surfaces, const node_t* const node)
{
    surface_t*      p;
    surface_t*      next;
    surface_t*      frontlist;
    surface_t*      backlist;
    surface_t*      frontfrag;
    surface_t*      backfrag;
    dplane_t*       splitplane;

    splitplane = &g_dplanes[node->planenum];

    frontlist = NULL;
    backlist = NULL;

    for (p = surfaces; p; p = next)
    {
        next = p->next;
        DivideSurface(p, splitplane, &frontfrag, &backfrag);

        if (frontfrag)
        {
            if (!frontfrag->faces)
            {
                Error("surface with no faces");
            }
            frontfrag->next = frontlist;
            frontlist = frontfrag;
        }
        if (backfrag)
        {
            if (!backfrag->faces)
            {
                Error("surface with no faces");
            }
            backfrag->next = backlist;
            backlist = backfrag;
        }
    }

    node->children[0]->surfaces = frontlist;
    node->children[1]->surfaces = backlist;
}
#ifdef ZHLT_DETAILBRUSH
static void SplitNodeBrushes (brush_t *brushes, const node_t *node)
{
	brush_t *frontlist, *frontfrag;
	brush_t *backlist, *backfrag;
	brush_t *b, *next;
	const dplane_t *splitplane;
	frontlist = NULL;
	backlist = NULL;
	splitplane = &g_dplanes[node->planenum];
	for (b = brushes; b; b = next)
	{
		next = b->next;
		SplitBrush (b, splitplane, &frontfrag, &backfrag);
		if (frontfrag)
		{
			frontfrag->next = frontlist;
			frontlist = frontfrag;
		}
		if (backfrag)
		{
			backfrag->next = backlist;
			backlist = backfrag;
		}
	}
	node->children[0]->detailbrushes = frontlist;
	node->children[1]->detailbrushes = backlist;
}
#endif

// =====================================================================================
//  RankForContents
// =====================================================================================
static int      RankForContents(const int contents)
{
    //Log("SolidBSP::RankForContents - contents type is %i ",contents);
    switch (contents)
    {
#ifdef ZHLT_NULLTEX    // AJM
#ifndef HLCSG_HLBSP_CONTENTSNULL_FIX
    case CONTENTS_NULL:
        //Log("(null)\n");
        //return 13;
        return -2;
#endif
#endif

    case CONTENTS_EMPTY:
        //Log("(empty)\n");
        return 0;
    case CONTENTS_WATER:
        //Log("(water)\n");
        return 1;
    case CONTENTS_TRANSLUCENT:
        //Log("(traslucent)\n");
        return 2;
    case CONTENTS_CURRENT_0:
        //Log("(current_0)\n");
        return 3;
    case CONTENTS_CURRENT_90:
        //Log("(current_90)\n");
        return 4;
    case CONTENTS_CURRENT_180:
        //Log("(current_180)\n");
        return 5;
    case CONTENTS_CURRENT_270:
        //Log("(current_270)\n");
        return 6;
    case CONTENTS_CURRENT_UP:
        //Log("(current_up)\n");
        return 7;
    case CONTENTS_CURRENT_DOWN:
        //Log("(current_down)\n");
        return 8;
    case CONTENTS_SLIME:
        //Log("(slime)\n");
        return 9;
    case CONTENTS_LAVA:
        //Log("(lava)\n");
        return 10;
    case CONTENTS_SKY:
        //Log("(sky)\n");
        return 11;
    case CONTENTS_SOLID:
        //Log("(solid)\n");
        return 12;

#ifdef ZHLT_DETAIL
    case CONTENTS_DETAIL:
        return 13;
        //Log("(detail)\n");
#endif

    default:
        hlassert(false);
        Error("RankForContents: bad contents %i", contents);
    }
    return -1;
}

// =====================================================================================
//  ContentsForRank
// =====================================================================================
static int      ContentsForRank(const int rank)
{
    switch (rank)
    {
#ifdef ZHLT_NULLTEX // AJM
#ifndef HLCSG_HLBSP_CONTENTSNULL_FIX
    case -2:
        return CONTENTS_NULL;        // has at leat one face with null
#endif
#endif

    case -1:
#ifdef HLCSG_HLBSP_ALLOWEMPTYENTITY
		return CONTENTS_EMPTY;
#else
        return CONTENTS_SOLID;                             // no faces at all
#endif
    case 0:
        return CONTENTS_EMPTY;
    case 1:
        return CONTENTS_WATER;
    case 2:
        return CONTENTS_TRANSLUCENT;
    case 3:
        return CONTENTS_CURRENT_0;
    case 4:
        return CONTENTS_CURRENT_90;
    case 5:
        return CONTENTS_CURRENT_180;
    case 6:
        return CONTENTS_CURRENT_270;
    case 7:
        return CONTENTS_CURRENT_UP;
    case 8:
        return CONTENTS_CURRENT_DOWN;
    case 9:
        return CONTENTS_SLIME;
    case 10:
        return CONTENTS_LAVA;
    case 11:
        return CONTENTS_SKY;
    case 12:
        return CONTENTS_SOLID;

#ifdef ZHLT_DETAIL // AJM
    case 13:
        return CONTENTS_DETAIL;
#endif

    default:
        hlassert(false);
        Error("ContentsForRank: bad rank %i", rank);
    }
    return -1;
}

// =====================================================================================
//  FreeLeafSurfs
// =====================================================================================
static void     FreeLeafSurfs(node_t* leaf)
{
    surface_t*      surf;
    surface_t*      snext;
    face_t*         f;
    face_t*         fnext;

    for (surf = leaf->surfaces; surf; surf = snext)
    {
        snext = surf->next;
        for (f = surf->faces; f; f = fnext)
        {
            fnext = f->next;
            FreeFace(f);
        }
        FreeSurface(surf);
    }

    leaf->surfaces = NULL;
}
#ifdef ZHLT_DETAILBRUSH
static void FreeLeafBrushes (node_t *leaf)
{
	brush_t *b, *next;
	for (b = leaf->detailbrushes; b; b = next)
	{
		next = b->next;
		FreeBrush (b);
	}
	leaf->detailbrushes = NULL;
}
#endif

// =====================================================================================
//  LinkLeafFaces
//      Determines the contents of the leaf and creates the final list of original faces 
//      that have some fragment inside this leaf
// =====================================================================================
#ifdef HLBSP_MAX_LEAF_FACES
#define	MAX_LEAF_FACES	16384
#else
#define	MAX_LEAF_FACES	1024
#endif

#ifdef HLBSP_WARNMIXEDCONTENTS
const char*     ContentsToString(int contents)
{
    switch (contents)
    {
    case CONTENTS_EMPTY:
        return "EMPTY";
    case CONTENTS_SOLID:
        return "SOLID";
    case CONTENTS_WATER:
        return "WATER";
    case CONTENTS_SLIME:
        return "SLIME";
    case CONTENTS_LAVA:
        return "LAVA";
    case CONTENTS_SKY:
        return "SKY";
    case CONTENTS_CURRENT_0:
        return "CURRENT_0";
    case CONTENTS_CURRENT_90:
        return "CURRENT_90";
    case CONTENTS_CURRENT_180:
        return "CURRENT_180";
    case CONTENTS_CURRENT_270:
        return "CURRENT_270";
    case CONTENTS_CURRENT_UP:
        return "CURRENT_UP";
    case CONTENTS_CURRENT_DOWN:
        return "CURRENT_DOWN";
    case CONTENTS_TRANSLUCENT:
        return "TRANSLUCENT";
    default:
        return "UNKNOWN";
    }
}
#endif
static void     LinkLeafFaces(surface_t* planelist, node_t* leafnode)
{
    face_t*         f;
    surface_t*      surf;
    int             rank, r;
#ifndef ZHLT_DETAILBRUSH
    int             nummarkfaces;
    face_t*         markfaces[MAX_LEAF_FACES];

    leafnode->faces = NULL;
    leafnode->planenum = -1;
#endif

    rank = -1;
    for (surf = planelist; surf; surf = surf->next)
    {
#ifdef ZHLT_DETAILBRUSH
		if (!surf->onnode)
		{
			continue;
		}
#endif
        for (f = surf->faces; f; f = f->next)
        {
            if ((f->contents == CONTENTS_HINT))
            {
                f->contents = CONTENTS_EMPTY;
            }
#ifdef ZHLT_DETAILBRUSH
			if (f->detaillevel)
			{
				continue;
			}
#endif
            r = RankForContents(f->contents);
            if (r > rank)
            {
                rank = r;
            }
        }
    }
#ifdef HLBSP_WARNMIXEDCONTENTS
	for (surf = planelist; surf; surf = surf->next)
	{
#ifdef ZHLT_DETAILBRUSH
		if (!surf->onnode)
		{
			continue;
		}
#endif
		for (f = surf->faces; f; f = f->next)
		{
#ifdef ZHLT_DETAILBRUSH
			if (f->detaillevel)
			{
				continue;
			}
#endif
			r = RankForContents(f->contents);
			if (r != rank)
				break;
		}
		if (f)
			break;
	}
	if (surf)
	{
		entity_t *ent = EntityForModel (g_nummodels - 1);
		if (g_nummodels - 1 != 0 && ent == &g_entities[0])
		{
			ent = NULL;
		}
		Warning ("Ambiguous leafnode content ( %s and %s ) at (%.0f,%.0f,%.0f)-(%.0f,%.0f,%.0f) in hull %d of model %d (entity: classname \"%s\", origin \"%s\", targetname \"%s\")", 
			ContentsToString (ContentsForRank(r)), ContentsToString (ContentsForRank(rank)), 
			leafnode->mins[0], leafnode->mins[1], leafnode->mins[2], leafnode->maxs[0], leafnode->maxs[1], leafnode->maxs[2], 
			g_hullnum, g_nummodels - 1, 
			(ent? ValueForKey (ent, "classname"): "unknown"), 
			(ent? ValueForKey (ent, "origin"): "unknown"), 
			(ent? ValueForKey (ent, "targetname"): "unknown"));
		for (surface_t *surf2 = planelist; surf2; surf2 = surf2->next)
		{
			for (face_t *f2 = surf2->faces; f2; f2 = f2->next)
			{
				Developer (DEVELOPER_LEVEL_SPAM, "content = %d plane = %d normal = (%g,%g,%g)\n", f2->contents, f2->planenum, 
					g_dplanes[f2->planenum].normal[0], g_dplanes[f2->planenum].normal[1], g_dplanes[f2->planenum].normal[2]);
				for (int i = 0; i < f2->numpoints; i++)
				{
					Developer (DEVELOPER_LEVEL_SPAM, "(%g,%g,%g)\n", f2->pts[i][0], f2->pts[i][1], f2->pts[i][2]);
				}
			}
		}
	}
#endif

    leafnode->contents = ContentsForRank(rank);

#ifndef ZHLT_DETAILBRUSH
    if (leafnode->contents != CONTENTS_SOLID)
    {
        nummarkfaces = 0;
        for (surf = leafnode->surfaces; surf; surf = surf->next)
        {
            for (f = surf->faces; f; f = f->next)
            {
#ifdef HLCSG_HLBSP_SOLIDHINT
				if (f->original == NULL)
				{
					continue;
				}
#endif
                hlassume(nummarkfaces < MAX_LEAF_FACES, assume_MAX_LEAF_FACES);

                markfaces[nummarkfaces++] = f->original;
            }
        }

        markfaces[nummarkfaces] = NULL;                    // end marker
        nummarkfaces++;

        leafnode->markfaces = (face_t**)malloc(nummarkfaces * sizeof(*leafnode->markfaces));
        memcpy(leafnode->markfaces, markfaces, nummarkfaces * sizeof(*leafnode->markfaces));
    }

    FreeLeafSurfs(leafnode);
    leafnode->surfaces = NULL;
#endif
}
#ifdef ZHLT_DETAILBRUSH
static void MakeLeaf (node_t *leafnode)
{
	int				nummarkfaces;
	face_t			*markfaces[MAX_LEAF_FACES + 1];
	surface_t		*surf;
	face_t			*f;	

    leafnode->planenum = -1;

	leafnode->iscontentsdetail = leafnode->detailbrushes != NULL;
	FreeLeafBrushes (leafnode);
	leafnode->detailbrushes = NULL;

	if (!(leafnode->isportalleaf && leafnode->contents == CONTENTS_SOLID))
	{
		nummarkfaces = 0;
		for (surf = leafnode->surfaces; surf; surf = surf->next)
		{
			if (!surf->onnode)
			{
				continue;
			}
#ifdef HLCSG_HLBSP_SOLIDHINT
			if (!surf->onnode)
			{
				continue;
			}
#endif
			for (f = surf->faces; f; f = f->next)
			{
				if (f->original == NULL)
				{ // because it is not on node or its content is solid
					continue;
				}
#ifdef HLCSG_HLBSP_SOLIDHINT
				if (f->original == NULL)
				{
					continue;
				}
#endif
				hlassume(nummarkfaces < MAX_LEAF_FACES, assume_MAX_LEAF_FACES);

				markfaces[nummarkfaces++] = f->original;
			}
		}
		markfaces[nummarkfaces] = NULL;                    // end marker
		nummarkfaces++;

		leafnode->markfaces = (face_t**)malloc(nummarkfaces * sizeof(*leafnode->markfaces));
		memcpy(leafnode->markfaces, markfaces, nummarkfaces * sizeof(*leafnode->markfaces));
	}

	FreeLeafSurfs(leafnode);
	leafnode->surfaces = NULL;

}
#endif

// =====================================================================================
//  MakeNodePortal
//      Create the new portal by taking the full plane winding for the cutting plane and 
//      clipping it by all of the planes from the other portals.
//      Each portal tracks the node that created it, so unused nodes can be removed later.
// =====================================================================================
static void     MakeNodePortal(node_t* node)
{
    portal_t*       new_portal;
    portal_t*       p;
    dplane_t*       plane;
    dplane_t        clipplane;
    Winding *       w;
    int             side = 0;

    plane = &g_dplanes[node->planenum];
    w = new Winding(*plane);

    new_portal = AllocPortal();
    new_portal->plane = *plane;
    new_portal->onnode = node;

    for (p = node->portals; p; p = p->next[side])
    {
        clipplane = p->plane;
        if (p->nodes[0] == node)
        {
            side = 0;
        }
        else if (p->nodes[1] == node)
        {
            clipplane.dist = -clipplane.dist;
            VectorSubtract(vec3_origin, clipplane.normal, clipplane.normal);
            side = 1;
        }
        else
        {
            Error("MakeNodePortal: mislinked portal");
        }

        w->Clip(clipplane, true);
#ifdef ZHLT_WINDING_RemoveColinearPoints_VL
		if (w->m_NumPoints == 0)
#else
        if (!w)
#endif
        {
#ifdef ZHLT_WINDING_RemoveColinearPoints_VL
			Developer (DEVELOPER_LEVEL_WARNING, 
#else
            Warning(
#endif
				"MakeNodePortal:new portal was clipped away from node@(%.0f,%.0f,%.0f)-(%.0f,%.0f,%.0f)",
                    node->mins[0], node->mins[1], node->mins[2], node->maxs[0], node->maxs[1], node->maxs[2]);
            FreePortal(new_portal);
            return;
        }
    }

    new_portal->winding = w;
    AddPortalToNodes(new_portal, node->children[0], node->children[1]);
}

// =====================================================================================
//  SplitNodePortals
//      Move or split the portals that bound node so that the node's children have portals instead of node.
// =====================================================================================
static void     SplitNodePortals(node_t *node)
{
    portal_t*       p;
    portal_t*       next_portal;
    portal_t*       new_portal;
    node_t*         f;
    node_t*         b;
    node_t*         other_node;
    int             side = 0;
    dplane_t*       plane;
    Winding*        frontwinding;
    Winding*        backwinding;

    plane = &g_dplanes[node->planenum];
    f = node->children[0];
    b = node->children[1];

    for (p = node->portals; p; p = next_portal)
    {
        if (p->nodes[0] == node)
        {
            side = 0;
        }
        else if (p->nodes[1] == node)
        {
            side = 1;
        }
        else
        {
            Error("SplitNodePortals: mislinked portal");
        }
        next_portal = p->next[side];

        other_node = p->nodes[!side];
        RemovePortalFromNode(p, p->nodes[0]);
        RemovePortalFromNode(p, p->nodes[1]);

        // cut the portal into two portals, one on each side of the cut plane
        p->winding->Divide(*plane, &frontwinding, &backwinding);

#ifdef ZHLT_WINDING_RemoveColinearPoints_VL
		if (!frontwinding && !backwinding)
		{
			continue;
		}
#endif
        if (!frontwinding)
        {
            if (side == 0)
            {
                AddPortalToNodes(p, b, other_node);
            }
            else
            {
                AddPortalToNodes(p, other_node, b);
            }
            continue;
        }
        if (!backwinding)
        {
            if (side == 0)
            {
                AddPortalToNodes(p, f, other_node);
            }
            else
            {
                AddPortalToNodes(p, other_node, f);
            }
            continue;
        }

        // the winding is split
        new_portal = AllocPortal();
        *new_portal = *p;
        new_portal->winding = backwinding;
        delete p->winding;
        p->winding = frontwinding;

        if (side == 0)
        {
            AddPortalToNodes(p, f, other_node);
            AddPortalToNodes(new_portal, b, other_node);
        }
        else
        {
            AddPortalToNodes(p, other_node, f);
            AddPortalToNodes(new_portal, other_node, b);
        }
    }

    node->portals = NULL;
}

// =====================================================================================
//  CalcNodeBounds
//      Determines the boundaries of a node by minmaxing all the portal points, whcih 
//      completely enclose the node.
//      Returns true if the node should be midsplit.(very large)
// =====================================================================================
static bool     CalcNodeBounds(node_t* node
#ifdef HLBSP_MAXNODESIZE_SKYBOX
							   , vec3_t validmins, vec3_t validmaxs
#endif
							   )
{
    int             i;
    int             j;
    vec_t           v;
    portal_t*       p;
    portal_t*       next_portal;
    int             side = 0;

#ifdef ZHLT_DETAILBRUSH
	if (node->isdetail)
	{
		return false;
	}
#endif
#ifdef ZHLT_LARGERANGE
    node->mins[0] = node->mins[1] = node->mins[2] = BOGUS_RANGE;
    node->maxs[0] = node->maxs[1] = node->maxs[2] = -BOGUS_RANGE;
#else
    node->mins[0] = node->mins[1] = node->mins[2] = 9999;
    node->maxs[0] = node->maxs[1] = node->maxs[2] = -9999;
#endif

    for (p = node->portals; p; p = next_portal)
    {
        if (p->nodes[0] == node)
        {
            side = 0;
        }
        else if (p->nodes[1] == node)
        {
            side = 1;
        }
        else
        {
            Error("CalcNodeBounds: mislinked portal");
        }
        next_portal = p->next[side];

        for (i = 0; i < p->winding->m_NumPoints; i++)
        {
            for (j = 0; j < 3; j++)
            {
                v = p->winding->m_Points[i][j];
                if (v < node->mins[j])
                {
                    node->mins[j] = v;
                }
                if (v > node->maxs[j])
                {
                    node->maxs[j] = v;
                }
            }
        }
    }

#ifdef ZHLT_DETAILBRUSH
	if (node->isportalleaf)
	{
		return false;
	}
#endif
#ifdef HLBSP_MAXNODESIZE_SKYBOX
	for (i = 0; i < 3; i++)
	{
		validmins[i] = max (node->mins[i], -(ENGINE_ENTITY_RANGE + g_maxnode_size));
		validmaxs[i] = min (node->maxs[i], ENGINE_ENTITY_RANGE + g_maxnode_size);
	}
	for (i = 0; i < 3; i++)
	{
		if (validmaxs[i] - validmins[i] <= ON_EPSILON)
		{
			return false;
		}
	}
	for (i = 0; i < 3; i++)
	{
		if (validmaxs[i] - validmins[i] > g_maxnode_size + ON_EPSILON)
		{
			return true;
		}
	}
#else
    for (i = 0; i < 3; i++)
    {
        if (node->maxs[i] - node->mins[i] > g_maxnode_size)
        {
            return true;
        }
    }
#endif
    return false;
}

// =====================================================================================
//  CopyFacesToNode
//      Do a final merge attempt, then subdivide the faces to surface cache size if needed.
//      These are final faces that will be drawable in the game.
//      Copies of these faces are further chopped up into the leafs, but they will reference these originals.
// =====================================================================================
static void     CopyFacesToNode(node_t* node, surface_t* surf)
{
    face_t**        prevptr;
    face_t*         f;
    face_t*         newf;

    // merge as much as possible
    MergePlaneFaces(surf);

    // subdivide large faces
    prevptr = &surf->faces;
    while (1)
    {
        f = *prevptr;
        if (!f)
        {
            break;
        }
        SubdivideFace(f, prevptr);
        f = *prevptr;
        prevptr = &f->next;
    }

    // copy the faces to the node, and consider them the originals
    node->surfaces = NULL;
    node->faces = NULL;
    for (f = surf->faces; f; f = f->next)
    {
#ifdef HLCSG_HLBSP_SOLIDHINT
		if (f->facestyle == face_discardable)
		{
			continue;
		}
#endif
        if (f->contents != CONTENTS_SOLID)
        {
            newf = AllocFace();
            *newf = *f;
            f->original = newf;
            newf->next = node->faces;
            node->faces = newf;
        }
    }
}

// =====================================================================================
//  BuildBspTree_r
// =====================================================================================
static void     BuildBspTree_r(node_t* node)
{
    surface_t*      split;
    bool            midsplit;
    surface_t*      allsurfs;
#ifdef HLBSP_MAXNODESIZE_SKYBOX
	vec3_t			validmins, validmaxs;
#endif

    midsplit = CalcNodeBounds(node
#ifdef HLBSP_MAXNODESIZE_SKYBOX
		, validmins, validmaxs
#endif
		);

#ifdef ZHLT_DETAILBRUSH
	int splitdetaillevel = CalcSplitDetaillevel (node);
#ifdef HLCSG_HLBSP_SOLIDHINT
	FixDetaillevelForDiscardable (node, splitdetaillevel);
#endif
#endif
    split = SelectPartition(node->surfaces, node, midsplit
#ifdef ZHLT_DETAILBRUSH
		, splitdetaillevel
#endif
#ifdef HLBSP_MAXNODESIZE_SKYBOX
		, validmins, validmaxs
#endif
		);
#ifdef ZHLT_DETAILBRUSH
	if (!node->isdetail && (!split || split->detaillevel > 0))
	{
		node->isportalleaf = true;
		LinkLeafFaces (node->surfaces, node); // set contents
		if (node->contents == CONTENTS_SOLID)
		{
			split = NULL;
		}
	}
	else
	{
		node->isportalleaf = false;
	}
#endif
    if (!split)
    {                                                      // this is a leaf node
#ifdef ZHLT_DETAILBRUSH
		MakeLeaf (node);
#else
        node->planenum = PLANENUM_LEAF;
        LinkLeafFaces(node->surfaces, node);
#endif
        return;
    }

    // these are final polygons
    split->onnode = node;                                  // can't use again
    allsurfs = node->surfaces;
    node->planenum = split->planenum;
    node->faces = NULL;
    CopyFacesToNode(node, split);

    node->children[0] = AllocNode();
    node->children[1] = AllocNode();
#ifdef ZHLT_DETAILBRUSH
	node->children[0]->isdetail = split->detaillevel > 0;
	node->children[1]->isdetail = split->detaillevel > 0;
#endif

    // split all the polysurfaces into front and back lists
    SplitNodeSurfaces(allsurfs, node);
#ifdef ZHLT_DETAILBRUSH
	SplitNodeBrushes (node->detailbrushes, node);
#endif

#ifdef ZHLT_DETAILBRUSH
	if (!split->detaillevel)
	{
		MakeNodePortal (node);
		SplitNodePortals (node);
	}
#else
    // create the portal that seperates the two children
    MakeNodePortal(node);

    // carve the portals on the boundaries of the node
    SplitNodePortals(node);
#endif

    // recursively do the children
    BuildBspTree_r(node->children[0]);
    BuildBspTree_r(node->children[1]);
	UpdateStatus();
}

// =====================================================================================
//  SolidBSP
//      Takes a chain of surfaces plus a split type, and returns a bsp tree with faces 
//      off the nodes.
//      The original surface chain will be completely freed.
// =====================================================================================
node_t*         SolidBSP(const surfchain_t* const surfhead, 
#ifdef ZHLT_DETAILBRUSH
						 brush_t *detailbrushes, 
#endif
						 bool report_progress)
{
    node_t*         headnode;

	ResetStatus(report_progress);
	double start_time = I_FloatTime();
	if(report_progress)
	{
		Log("SolidBSP [hull %d] ",g_hullnum);
	}
	else
	{
	    Verbose("----- SolidBSP -----\n");
	}

    headnode = AllocNode();
    headnode->surfaces = surfhead->surfaces;
#ifdef ZHLT_DETAILBRUSH
	headnode->detailbrushes = detailbrushes;
	headnode->isdetail = false;
#endif

#ifndef HLCSG_HLBSP_ALLOWEMPTYENTITY
    if (!surfhead->surfaces)
    {
        // nothing at all to build
#ifdef ZHLT_DETAILBRUSH
		headnode->isportalleaf = true;
		headnode->iscontentsdetail = (headnode->detailbrushes != NULL);
		FreeLeafBrushes (headnode->detailbrushes);
#endif
        headnode->planenum = -1;
        headnode->contents = CONTENTS_EMPTY;
        return headnode;
    }
#endif

    // generate six portals that enclose the entire world
    MakeHeadnodePortals(headnode, surfhead->mins, surfhead->maxs);

    // recursively partition everything
    BuildBspTree_r(headnode);

	double end_time = I_FloatTime();
	if(report_progress)
	{
		Log("%d (%.2f seconds)\n",++g_numProcessed,(end_time - start_time));
	}

    return headnode;
}
