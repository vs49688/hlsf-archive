#include "qrad.h"

edgeshare_t     g_edgeshare[MAX_MAP_EDGES];
vec3_t          g_face_centroids[MAX_MAP_EDGES];
bool            g_sky_lighting_fix = DEFAULT_SKY_LIGHTING_FIX;
#ifdef HLRAD_SMOOTH_TEXNORMAL
vec3_t          g_face_texnormals[MAX_MAP_FACES];
#endif

#define TEXTURE_STEP   16.0

#ifdef HLRAD_SMOOTH_TEXNORMAL
bool GetIntertexnormal (int facenum1, int facenum2, vec_t *out)
{
	vec3_t normal;
	const dplane_t *p1 = getPlaneFromFaceNumber (facenum1);
	const dplane_t *p2 = getPlaneFromFaceNumber (facenum2);
	VectorAdd (g_face_texnormals[facenum1], g_face_texnormals[facenum2], normal);
	if (!VectorNormalize (normal)
		|| DotProduct (normal, p1->normal) <= NORMAL_EPSILON
		|| DotProduct (normal, p2->normal) <= NORMAL_EPSILON
		)
	{
		return false;
	}
	if (out)
	{
		VectorCopy (normal, out);
	}
	return true;
}
#endif
// =====================================================================================
//  PairEdges
// =====================================================================================
#ifdef HLRAD_SMOOTH_FACELIST
typedef struct
{
	int numclipplanes;
	dplane_t *clipplanes;
}
intersecttest_t;
bool TestFaceIntersect (intersecttest_t *t, int facenum)
{
	dface_t *f2 = &g_dfaces[facenum];
	Winding *w = new Winding (*f2);
	int k;
	for (k = 0; k < w->m_NumPoints; k++)
	{
		VectorAdd (w->m_Points[k], g_face_offset[facenum], w->m_Points[k]);
	}
	for (k = 0; k < t->numclipplanes; k++)
	{
		if (!w->Clip (t->clipplanes[k], false
#ifdef ZHLT_WINDING_EPSILON
			, ON_EPSILON*4
#endif
			))
		{
			break;
		}
	}
	bool intersect = w->m_NumPoints > 0;
	delete w;
	return intersect;
}
intersecttest_t *CreateIntersectTest (const dplane_t *p, int facenum)
{
	dface_t *f = &g_dfaces[facenum];
	intersecttest_t *t;
	t = (intersecttest_t *)malloc (sizeof (intersecttest_t));
	hlassume (t != NULL, assume_NoMemory);
	t->clipplanes = (dplane_t *)malloc (f->numedges * sizeof (dplane_t));
	hlassume (t->clipplanes != NULL, assume_NoMemory);
	t->numclipplanes = 0;
	int j;
	for (j = 0; j < f->numedges; j++)
	{
		// should we use winding instead?
		int edgenum = g_dsurfedges[f->firstedge + j];
		{
			vec3_t v0, v1;
			vec3_t dir, normal;
			if (edgenum < 0)
			{
				VectorCopy (g_dvertexes[g_dedges[-edgenum].v[1]].point, v0);
				VectorCopy (g_dvertexes[g_dedges[-edgenum].v[0]].point, v1);
			}
			else
			{
				VectorCopy (g_dvertexes[g_dedges[edgenum].v[0]].point, v0);
				VectorCopy (g_dvertexes[g_dedges[edgenum].v[1]].point, v1);
			}
			VectorAdd (v0, g_face_offset[facenum], v0);
			VectorAdd (v1, g_face_offset[facenum], v1);
			VectorSubtract (v1, v0, dir);
			CrossProduct (dir, p->normal, normal); // facing inward
			if (!VectorNormalize (normal))
			{
				continue;
			}
			VectorCopy (normal, t->clipplanes[t->numclipplanes].normal);
			t->clipplanes[t->numclipplanes].dist = DotProduct (v0, normal);
			t->numclipplanes++;
		}
	}
	return t;
}
void FreeIntersectTest (intersecttest_t *t)
{
	free (t->clipplanes);
	free (t);
}
#endif
#ifdef HLRAD_GetPhongNormal_VL
void AddFaceForVertexNormal_printerror (const int edgeabs, const int edgeend, dface_t *const f)
{
	if (DEVELOPER_LEVEL_WARNING <= g_developer)
	{
		int i, e;
		Log ("AddFaceForVertexNormal - bad face:\n");
		Log (" edgeabs=%d edgeend=%d\n", edgeabs, edgeend);
		for (i = 0; i < f->numedges; i++)
		{
			e = g_dsurfedges[f->firstedge + i];
			edgeshare_t *es = &g_edgeshare[abs(e)];
			int v0 = g_dedges[abs(e)].v[0], v1 = g_dedges[abs(e)].v[1];
			Log (" e=%d v0=%d(%f,%f,%f) v1=%d(%f,%f,%f) share0=%d share1=%d\n", e,
				v0, g_dvertexes[v0].point[0], g_dvertexes[v0].point[1], g_dvertexes[v0].point[2],
				v1, g_dvertexes[v1].point[0], g_dvertexes[v1].point[1], g_dvertexes[v1].point[2],
				(es->faces[0]==NULL? -1: es->faces[0]-g_dfaces), (es->faces[1]==NULL? -1: es->faces[1]-g_dfaces));
		}
	}
}
int AddFaceForVertexNormal (const int edgeabs, int &edgeabsnext, const int edgeend, int &edgeendnext, dface_t *const f, dface_t *&fnext, vec_t &angle, vec3_t &normal)
// Must guarantee these faces will form a loop or a chain, otherwise will result in endless loop.
//
//   e[end]/enext[endnext]
//  *
//  |\
//  |a\ fnext
//  |  \
//  | f \
//  |    \
//  e   enext
//
{
	VectorCopy(getPlaneFromFace(f)->normal, normal);
	int vnum = g_dedges[edgeabs].v[edgeend];
	int iedge, iedgenext, edge, edgenext;
	int i, e, count1, count2;
	vec_t dot;
	for (count1 = count2 = 0, i = 0; i < f->numedges; i++)
	{
		e = g_dsurfedges[f->firstedge + i];
		if (g_dedges[abs(e)].v[0] == g_dedges[abs(e)].v[1])
			continue;
		if (abs(e) == edgeabs)
		{
			iedge = i;
			edge = e;
			count1 ++;
		}
		else if (g_dedges[abs(e)].v[0] == vnum || g_dedges[abs(e)].v[1] == vnum)
		{
			iedgenext = i;
			edgenext = e;
			count2 ++;
		}
	}
	if (count1 != 1 || count2 != 1)
	{
		AddFaceForVertexNormal_printerror (edgeabs, edgeend, f);
		return -1;
	}
	int vnum11, vnum12, vnum21, vnum22;
	vec3_t vec1, vec2;
	vnum11 = g_dedges[abs(edge)].v[edge>0?0:1];
	vnum12 = g_dedges[abs(edge)].v[edge>0?1:0];
	vnum21 = g_dedges[abs(edgenext)].v[edgenext>0?0:1];
	vnum22 = g_dedges[abs(edgenext)].v[edgenext>0?1:0];
	if (vnum == vnum12 && vnum == vnum21 && vnum != vnum11 && vnum != vnum22)
	{
		VectorSubtract(g_dvertexes[vnum11].point, g_dvertexes[vnum].point, vec1);
		VectorSubtract(g_dvertexes[vnum22].point, g_dvertexes[vnum].point, vec2);
		edgeabsnext = abs(edgenext);
		edgeendnext = edgenext>0?0:1;
	}
	else if (vnum == vnum11 && vnum == vnum22 && vnum != vnum12 && vnum != vnum21)
	{
		VectorSubtract(g_dvertexes[vnum12].point, g_dvertexes[vnum].point, vec1);
		VectorSubtract(g_dvertexes[vnum21].point, g_dvertexes[vnum].point, vec2);
		edgeabsnext = abs(edgenext);
		edgeendnext = edgenext>0?1:0;
	}
	else
	{
		AddFaceForVertexNormal_printerror (edgeabs, edgeend, f);
		return -1;
	}
	VectorNormalize(vec1);
	VectorNormalize(vec2);
	dot = DotProduct(vec1,vec2);
	dot = dot>1? 1: dot<-1? -1: dot;
	angle = acos(dot);
	edgeshare_t *es = &g_edgeshare[edgeabsnext];
	if (!(es->faces[0] && es->faces[1]))
		return 1;
	if (es->faces[0] == f && es->faces[1] != f)
		fnext = es->faces[1];
	else if (es->faces[1] == f && es->faces[0] != f)
		fnext = es->faces[0];
	else
	{
		AddFaceForVertexNormal_printerror (edgeabs, edgeend, f);
		return -1;
	}
	return 0;
}
#endif
void            PairEdges()
{
    int             i, j, k;
    dface_t*        f;
    edgeshare_t*    e;

    memset(&g_edgeshare, 0, sizeof(g_edgeshare));

    f = g_dfaces;
    for (i = 0; i < g_numfaces; i++, f++)
    {
#ifdef HLRAD_SMOOTH_TEXNORMAL
		{
			const dplane_t *fp = getPlaneFromFace (f);
			vec3_t texnormal;
			const texinfo_t *tex = &g_texinfo[f->texinfo];
			CrossProduct (tex->vecs[1], tex->vecs[0], texnormal);
			VectorNormalize (texnormal);
			if (DotProduct (texnormal, fp->normal) < 0)
			{
				VectorSubtract (vec3_origin, texnormal, texnormal);
			}
			VectorCopy (texnormal, g_face_texnormals[i]);
		}
#endif
#ifdef HLRAD_EDGESHARE_NOSPECIAL
		if (g_texinfo[f->texinfo].flags & TEX_SPECIAL)
		{
			continue;
		}
#endif
        for (j = 0; j < f->numedges; j++)
        {
            k = g_dsurfedges[f->firstedge + j];
            if (k < 0)
            {
                e = &g_edgeshare[-k];

                hlassert(e->faces[1] == NULL);
                e->faces[1] = f;
            }
            else
            {
                e = &g_edgeshare[k];

                hlassert(e->faces[0] == NULL);
                e->faces[0] = f;
            }

            if (e->faces[0] && e->faces[1]) {
				// determine if coplanar
				if (e->faces[0]->planenum == e->faces[1]->planenum
#ifdef HLRAD_PairEdges_FACESIDE_FIX
					&& e->faces[0]->side == e->faces[1]->side
#endif
					) {
						e->coplanar = true;
				} else {
                    // see if they fall into a "smoothing group" based on angle of the normals
                    vec3_t          normals[2];

                    VectorCopy(getPlaneFromFace(e->faces[0])->normal, normals[0]);
                    VectorCopy(getPlaneFromFace(e->faces[1])->normal, normals[1]);

                    e->cos_normals_angle = DotProduct(normals[0], normals[1]);

#ifdef HLRAD_CUSTOMSMOOTH
					vec_t smoothvalue;
					int m0 = g_texinfo[e->faces[0]->texinfo].miptex;
					int m1 = g_texinfo[e->faces[1]->texinfo].miptex;
					smoothvalue = max (smoothvalues[m0], smoothvalues[m1]);
					if (m0 != m1)
					{
						smoothvalue = max (smoothvalue, g_smoothing_threshold_2);
					}
					if (smoothvalue >= 1.0 - NORMAL_EPSILON)
					{
						smoothvalue = 2.0;
					}
#endif
                    if (e->cos_normals_angle > (1.0 - NORMAL_EPSILON))
                    {
                        e->coplanar = true;
                    }
#ifndef HLRAD_CUSTOMSMOOTH
                    else if (g_smoothing_threshold > 0.0)
#endif
                    {
#ifdef HLRAD_CUSTOMSMOOTH
                        if (e->cos_normals_angle >= max (smoothvalue - NORMAL_EPSILON, NORMAL_EPSILON))
#else
                        if (e->cos_normals_angle >= g_smoothing_threshold)
#endif
                        {
                            VectorAdd(normals[0], normals[1], e->interface_normal);
                            VectorNormalize(e->interface_normal);
                        }
                    }
                }
#ifdef HLRAD_TRANSLUCENT
				if (!VectorCompare (g_translucenttextures[g_texinfo[e->faces[0]->texinfo].miptex], g_translucenttextures[g_texinfo[e->faces[1]->texinfo].miptex]))
				{
					e->coplanar = false;
					VectorClear (e->interface_normal);
				}
#endif
#ifdef HLRAD_GetPhongNormal_VL
				if (e->coplanar)
					VectorCopy(getPlaneFromFace(e->faces[0])->normal, e->interface_normal);
				if (!VectorCompare(e->interface_normal, vec3_origin))
					e->smooth = true;
#ifdef HLRAD_SMOOTH_TEXNORMAL
				if (!GetIntertexnormal (e->faces[0] - g_dfaces, e->faces[1] - g_dfaces))
				{
					e->coplanar = false;
					VectorClear (e->interface_normal);
					e->smooth = false;
				}
#endif
#endif
            }
        }
    }
#ifdef HLRAD_GetPhongNormal_VL
	{
		int edgeabs, edgeabsnext;
		int edgeend, edgeendnext;
		int d;
		dface_t *f, *fcurrent, *fnext;
		vec_t angle, angles;
		vec3_t normal, normals;
		vec3_t edgenormal;
		int r, count;
		for (edgeabs = 0; edgeabs < MAX_MAP_EDGES; edgeabs++)
		{
			e = &g_edgeshare[edgeabs];
			if (!e->smooth)
				continue;
			VectorCopy(e->interface_normal, edgenormal);
			if (g_dedges[edgeabs].v[0] == g_dedges[edgeabs].v[1])
			{
				vec3_t errorpos;
				VectorCopy (g_dvertexes[g_dedges[edgeabs].v[0]].point, errorpos);
				VectorAdd (errorpos, g_face_offset[e->faces[0] - g_dfaces], errorpos);
				Warning ("PairEdges: invalid edge at (%f,%f,%f)", errorpos[0], errorpos[1], errorpos[2]);
				VectorCopy(edgenormal, e->vertex_normal[0]);
				VectorCopy(edgenormal, e->vertex_normal[1]);
			}
			else
			{
				const dplane_t *p0 = getPlaneFromFace (e->faces[0]);
				const dplane_t *p1 = getPlaneFromFace (e->faces[1]);
#ifdef HLRAD_SMOOTH_FACELIST
				intersecttest_t *test0 = CreateIntersectTest (p0, e->faces[0] - g_dfaces);
				intersecttest_t *test1 = CreateIntersectTest (p1, e->faces[1] - g_dfaces);
#endif
				for (edgeend = 0; edgeend < 2; edgeend++)
				{
					vec3_t errorpos;
					VectorCopy (g_dvertexes[g_dedges[edgeabs].v[edgeend]].point, errorpos);
					VectorAdd (errorpos, g_face_offset[e->faces[0] - g_dfaces], errorpos);
					angles = 0;
					VectorClear (normals);

					for (d = 0; d < 2; d++)
					{
						f = e->faces[d];
						count = 0, fnext = f, edgeabsnext = edgeabs, edgeendnext = edgeend;
						while (1)
						{
							fcurrent = fnext;
							r = AddFaceForVertexNormal (edgeabsnext, edgeabsnext, edgeendnext, edgeendnext, fcurrent, fnext, angle, normal);
							count++;
							if (r == -1)
							{
								Warning ("PairEdges: face edges mislink at (%f,%f,%f)", errorpos[0], errorpos[1], errorpos[2]);
								break;
							}
							if (count >= 100)
							{
								Warning ("PairEdges: faces mislink at (%f,%f,%f)", errorpos[0], errorpos[1], errorpos[2]);
								break;
							}
							if (DotProduct (normal, p0->normal) <= NORMAL_EPSILON || DotProduct(normal, p1->normal) <= NORMAL_EPSILON)
								break;
	#ifdef HLRAD_CUSTOMSMOOTH
							vec_t smoothvalue;
							int m0 = g_texinfo[f->texinfo].miptex;
							int m1 = g_texinfo[fcurrent->texinfo].miptex;
							smoothvalue = max (smoothvalues[m0], smoothvalues[m1]);
							if (m0 != m1)
							{
								smoothvalue = max (smoothvalue, g_smoothing_threshold_2);
							}
							if (smoothvalue >= 1.0 - NORMAL_EPSILON)
							{
								smoothvalue = 2.0;
							}
							if (DotProduct (edgenormal, normal) < max (smoothvalue - NORMAL_EPSILON, NORMAL_EPSILON))
	#else
							if (DotProduct (edgenormal, normal) + NORMAL_EPSILON < g_smoothing_threshold)
	#endif
								break;
	#ifdef HLRAD_SMOOTH_TEXNORMAL
							if (!GetIntertexnormal (fcurrent - g_dfaces, e->faces[0] - g_dfaces) || !GetIntertexnormal (fcurrent - g_dfaces, e->faces[1] - g_dfaces))
								break;
	#endif
	#ifdef HLRAD_SMOOTH_FACELIST
							if (fcurrent != e->faces[0] && fcurrent != e->faces[1] &&
								(TestFaceIntersect (test0, fcurrent - g_dfaces) || TestFaceIntersect (test1, fcurrent - g_dfaces)))
							{
								Developer (DEVELOPER_LEVEL_WARNING, "Overlapping faces around corner (%f,%f,%f)\n", errorpos[0], errorpos[1], errorpos[2]);
								break;
							}
	#endif
							angles += angle;
							VectorMA(normals, angle, normal, normals);
	#ifdef HLRAD_SMOOTH_FACELIST
							{
								bool in = false;
								if (fcurrent == e->faces[0] || fcurrent == e->faces[1])
								{
									in = true;
								}
								for (facelist_t *l = e->vertex_facelist[edgeend]; l; l = l->next)
								{
									if (fcurrent == l->face)
									{
										in = true;
									}
								}
								if (!in)
								{
									facelist_t *l = (facelist_t *)malloc (sizeof (facelist_t));
									hlassume (l != NULL, assume_NoMemory);
									l->face = fcurrent;
									l->next = e->vertex_facelist[edgeend];
									e->vertex_facelist[edgeend] = l;
								}
							}
	#endif
							if (r != 0 || fnext == f)
								break;
						}
					}

					if (angles < NORMAL_EPSILON)
					{
						VectorCopy(edgenormal, e->vertex_normal[edgeend]);
						Warning ("PairEdges: no valid faces at (%f,%f,%f)", errorpos[0], errorpos[1], errorpos[2]);
					}
					else
					{
						VectorNormalize(normals);
						VectorCopy(normals, e->vertex_normal[edgeend]);
					}
				}
#ifdef HLRAD_SMOOTH_FACELIST
				FreeIntersectTest (test0);
				FreeIntersectTest (test1);
#endif
			}
			if (e->coplanar)
			{
				if (!VectorCompare (e->vertex_normal[0], e->interface_normal) || !VectorCompare (e->vertex_normal[1], e->interface_normal))
					e->coplanar = false;
			}
		}
	}
#endif
}

#define	MAX_SINGLEMAP	(18*18*4)

typedef struct
{
    vec_t*          light;
    vec_t           facedist;
    vec3_t          facenormal;
#ifdef HLRAD_TRANSLUCENT
	bool			translucent_b;
	vec3_t			translucent_v;
#endif

    int             numsurfpt;
    vec3_t          surfpt[MAX_SINGLEMAP];
#ifdef HLRAD_AddSampleToPatch_PRECISE
	vec3_t			surfpt_original[MAX_SINGLEMAP];
#endif

    vec3_t          texorg;
    vec3_t          worldtotex[2];                         // s = (world - texorg) . worldtotex[0]
    vec3_t          textoworld[2];                         // world = texorg + s * textoworld[0]
#ifdef HLRAD_CalcPoints_NEW
	vec3_t			texnormal;
#endif

    vec_t           exactmins[2], exactmaxs[2];

    int             texmins[2], texsize[2];
    int             lightstyles[256];
    int             surfnum;
    dface_t*        face;
}
lightinfo_t;
#ifdef HLRAD_MDL_LIGHT_HACK
#ifndef HLRAD_MDL_LIGHT_HACK_NEW
typedef struct
{
	vec3_t			texorg;
	vec3_t			offset;
	vec3_t			textoworld[2];
	vec3_t			worldtotex[2];
	int				texmins[2], texsize[2];
}
facesampleinfo_t;
static facesampleinfo_t facesampleinfo[MAX_MAP_FACES];
#endif
#endif

// =====================================================================================
//  TextureNameFromFace
// =====================================================================================
static const char* TextureNameFromFace(const dface_t* const f)
{
    texinfo_t*      tx;
    miptex_t*       mt;
    int             ofs;

    //
    // check for light emited by texture
    //
    tx = &g_texinfo[f->texinfo];

    ofs = ((dmiptexlump_t*)g_dtexdata)->dataofs[tx->miptex];
    mt = (miptex_t*)((byte*) g_dtexdata + ofs);

	return mt->name;
}

// =====================================================================================
//  CalcFaceExtents
//      Fills in s->texmins[] and s->texsize[]
//      also sets exactmins[] and exactmaxs[]
// =====================================================================================
#ifdef HLRAD_MAXEXTENT
bool g_warnedextent = false;
#endif
static void     CalcFaceExtents(lightinfo_t* l)
{
    const int       facenum = l->surfnum;
    dface_t*        s;
    vec_t           mins[2], maxs[2], val;
    int             i, j, e;
    dvertex_t*      v;
    texinfo_t*      tex;

    s = l->face;

    mins[0] = mins[1] = 999999;
    maxs[0] = maxs[1] = -99999; // a little small, but same with HL. --vluzacn

    tex = &g_texinfo[s->texinfo];

    for (i = 0; i < s->numedges; i++)
    {
        e = g_dsurfedges[s->firstedge + i];
        if (e >= 0)
        {
            v = g_dvertexes + g_dedges[e].v[0];
        }
        else
        {
            v = g_dvertexes + g_dedges[-e].v[1];
        }

        for (j = 0; j < 2; j++)
        {
            val = v->point[0] * tex->vecs[j][0] +
                v->point[1] * tex->vecs[j][1] + v->point[2] * tex->vecs[j][2] + tex->vecs[j][3];
            if (val < mins[j])
            {
                mins[j] = val;
            }
            if (val > maxs[j])
            {
                maxs[j] = val;
            }
        }
    }

    for (i = 0; i < 2; i++)
    {
        l->exactmins[i] = mins[i];
        l->exactmaxs[i] = maxs[i];

        mins[i] = floor(mins[i] / 16.0);
        maxs[i] = ceil(maxs[i] / 16.0);

		l->texmins[i] = mins[i];
        l->texsize[i] = maxs[i] - mins[i];
	}

	if (!(tex->flags & TEX_SPECIAL))
	{
#ifdef HLRAD_MAXEXTENT
		if (l->texsize[0] < 0 || l->texsize[1] < 0 || l->texsize[0] > 17 || l->texsize[1] > 17)
#else
		if ((l->texsize[0] > 16) || (l->texsize[1] > 16))
#endif
		{
			ThreadLock();
			PrintOnce("\nfor Face %d (texture %s) at ", s - g_dfaces, TextureNameFromFace(s));

			for (i = 0; i < s->numedges; i++)
			{
				e = g_dsurfedges[s->firstedge + i];
				if (e >= 0)
                {
					v = g_dvertexes + g_dedges[e].v[0];
                }
				else
                {
					v = g_dvertexes + g_dedges[-e].v[1];
                }
                VectorAdd(v->point, g_face_offset[facenum], v->point);
				Log("(%4.3f %4.3f %4.3f) ", v->point[0], v->point[1], v->point[2]);
			}
			Log("\n");

			Error( "Bad surface extents (%d x %d)\nCheck the file ZHLTProblems.html for a detailed explanation of this problem", l->texsize[0], l->texsize[1]);
		}
#ifdef HLRAD_MAXEXTENT
		if (l->texsize[0] > 16 || l->texsize[1] > 16)
		{
			if (!g_warnedextent)
			{ // only warn once
				ThreadLock ();
				if (!g_warnedextent)
				{
					g_warnedextent = true;
					Warning ("\nfor Face %d (texture %s) at ", s - g_dfaces, TextureNameFromFace(s));
					for (i = 0; i < s->numedges; i++)
					{
						e = g_dsurfedges[s->firstedge + i];
						if (e >= 0)
						{
							v = g_dvertexes + g_dedges[e].v[0];
						}
						else
						{
							v = g_dvertexes + g_dedges[-e].v[1];
						}
						VectorAdd(v->point, g_face_offset[facenum], v->point);
						Log("(%4.3f %4.3f %4.3f) ", v->point[0], v->point[1], v->point[2]);
					}
					Log("\n");
					Warning ("Surface extents (%d x %d) exceeded (%d x %d)\nThis map will not work in 'Software' video mode.", l->texsize[0], l->texsize[1], 16, 16);
				}
				ThreadUnlock ();
			}
		}
#endif
	}
}

// =====================================================================================
//  CalcFaceVectors
//      Fills in texorg, worldtotex. and textoworld
// =====================================================================================
static void     CalcFaceVectors(lightinfo_t* l)
{
    texinfo_t*      tex;
    int             i, j;
    vec3_t          texnormal;
    vec_t           distscale;
    vec_t           dist, len;

    tex = &g_texinfo[l->face->texinfo];

    // convert from float to double
    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 3; j++)
        {
            l->worldtotex[i][j] = tex->vecs[i][j];
        }
    }

    // calculate a normal to the texture axis.  points can be moved along this
    // without changing their S/T
    CrossProduct(tex->vecs[1], tex->vecs[0], texnormal);
    VectorNormalize(texnormal);

    // flip it towards plane normal
    distscale = DotProduct(texnormal, l->facenormal);
    if (distscale == 0.0)
    {
        const unsigned facenum = l->face - g_dfaces;
    
        ThreadLock();
        Log("Malformed face (%d) normal @ \n", facenum);
        Winding* w = new Winding(*l->face);
        {
            const unsigned numpoints = w->m_NumPoints;
            unsigned x;
            for (x=0; x<numpoints; x++)
            {
                VectorAdd(w->m_Points[x], g_face_offset[facenum], w->m_Points[x]);
            }
        }
        w->Print();
        delete w;
        ThreadUnlock();

        hlassume(false, assume_MalformedTextureFace);
    }

    if (distscale < 0)
    {
        distscale = -distscale;
        VectorSubtract(vec3_origin, texnormal, texnormal);
    }

    // distscale is the ratio of the distance along the texture normal to
    // the distance along the plane normal
    distscale = 1.0 / distscale;

    for (i = 0; i < 2; i++)
    {
        len = (float)VectorLength(l->worldtotex[i]);
        dist = DotProduct(l->worldtotex[i], l->facenormal);
        dist *= distscale;
        VectorMA(l->worldtotex[i], -dist, texnormal, l->textoworld[i]);
        VectorScale(l->textoworld[i], (1 / len) * (1 / len), l->textoworld[i]);
    }

    // calculate texorg on the texture plane
    for (i = 0; i < 3; i++)
    {
        l->texorg[i] = -tex->vecs[0][3] * l->textoworld[0][i] - tex->vecs[1][3] * l->textoworld[1][i];
    }

    // project back to the face plane
    dist = DotProduct(l->texorg, l->facenormal) - l->facedist - DEFAULT_HUNT_OFFSET;
    dist *= distscale;
    VectorMA(l->texorg, -dist, texnormal, l->texorg);
#ifdef HLRAD_CalcPoints_NEW
	VectorCopy (texnormal, l->texnormal);
#endif

}

// =====================================================================================
//  SetSurfFromST
// =====================================================================================
static void     SetSurfFromST(const lightinfo_t* const l, vec_t* surf, const vec_t s, const vec_t t)
{
    const int       facenum = l->surfnum;
    int             j;

    for (j = 0; j < 3; j++)
    {
        surf[j] = l->texorg[j] + l->textoworld[0][j] * s + l->textoworld[1][j] * t;
    }

    // Adjust for origin-based models
    VectorAdd(surf, g_face_offset[facenum], surf);
}

#ifndef HLRAD_CalcPoints_NEW
// =====================================================================================
//  FindSurfaceMidpoint
// =====================================================================================
static dleaf_t* FindSurfaceMidpoint(const lightinfo_t* const l, vec_t* midpoint)
{
    int             s, t;
    int             w, h;
    vec_t           starts, startt;
    vec_t           us, ut;

    vec3_t          broken_midpoint;
    vec3_t          surface_midpoint;
    int             inside_point_count;

    dleaf_t*        last_valid_leaf = NULL;
    dleaf_t*        leaf_mid;

    const int       facenum = l->surfnum;
    const dface_t*  f = g_dfaces + facenum;
    const dplane_t* p = getPlaneFromFace(f);

    const vec_t*    face_delta = g_face_offset[facenum];

    h = l->texsize[1] + 1;
    w = l->texsize[0] + 1;
    starts = (float)l->texmins[0] * 16;
    startt = (float)l->texmins[1] * 16;

    // General case
    inside_point_count = 0;
    VectorClear(surface_midpoint);
    for (t = 0; t < h; t++)
    {
        for (s = 0; s < w; s++)
        {
            us = starts + s * TEXTURE_STEP;
            ut = startt + t * TEXTURE_STEP;

            SetSurfFromST(l, midpoint, us, ut);
            if ((leaf_mid = PointInLeaf(midpoint)) != g_dleafs)
            {
                if ((leaf_mid->contents != CONTENTS_SKY) && (leaf_mid->contents != CONTENTS_SOLID))
                {
                    last_valid_leaf = leaf_mid;
                    inside_point_count++;
                    VectorAdd(surface_midpoint, midpoint, surface_midpoint);
                }
            }
        }
    }

    if (inside_point_count > 1)
    {
        vec_t           tmp = 1.0 / inside_point_count;

        VectorScale(surface_midpoint, tmp, midpoint);

        //Verbose("Trying general at (%4.3f %4.3f %4.3f) %d\n", surface_midpoint[0], surface_midpoint[1], surface_midpoint[2], inside_point_count);
        if (
            (leaf_mid =
             HuntForWorld(midpoint, face_delta, p, DEFAULT_HUNT_SIZE, DEFAULT_HUNT_SCALE, DEFAULT_HUNT_OFFSET)))
        {
            //Verbose("general method succeeded at (%4.3f %4.3f %4.3f)\n", midpoint[0], midpoint[1], midpoint[2]);
            return leaf_mid;
        }
        //Verbose("Tried general , failed at (%4.3f %4.3f %4.3f)\n", midpoint[0], midpoint[1], midpoint[2]);
    }
    else if (inside_point_count == 1)
    {
        //Verbose("Returning single point from general\n");
        VectorCopy(surface_midpoint, midpoint);
        return last_valid_leaf;
    }
    else
    {
        //Verbose("general failed (no points)\n");
    }

    // Try harder
    inside_point_count = 0;
    VectorClear(surface_midpoint);
    for (t = 0; t < h; t++)
    {
        for (s = 0; s < w; s++)
        {
            us = starts + s * TEXTURE_STEP;
            ut = startt + t * TEXTURE_STEP;

            SetSurfFromST(l, midpoint, us, ut);
            leaf_mid =
                HuntForWorld(midpoint, face_delta, p, DEFAULT_HUNT_SIZE, DEFAULT_HUNT_SCALE, DEFAULT_HUNT_OFFSET);
            if (leaf_mid != g_dleafs)
            {
                last_valid_leaf = leaf_mid;
                inside_point_count++;
                VectorAdd(surface_midpoint, midpoint, surface_midpoint);
            }
        }
    }

    if (inside_point_count > 1)
    {
        vec_t           tmp = 1.0 / inside_point_count;

        VectorScale(surface_midpoint, tmp, midpoint);

        if (
            (leaf_mid =
             HuntForWorld(midpoint, face_delta, p, DEFAULT_HUNT_SIZE, DEFAULT_HUNT_SCALE, DEFAULT_HUNT_OFFSET)))
        {
            //Verbose("best method succeeded at (%4.3f %4.3f %4.3f)\n", midpoint[0], midpoint[1], midpoint[2]);
            return leaf_mid;
        }
        //Verbose("Tried best, failed at (%4.3f %4.3f %4.3f)\n", midpoint[0], midpoint[1], midpoint[2]);
    }
    else if (inside_point_count == 1)
    {
        //Verbose("Returning single point from best\n");
        VectorCopy(surface_midpoint, midpoint);
        return last_valid_leaf;
    }
    else
    {
        //Verbose("best failed (no points)\n");
    }

    // Original broken code
    {
        vec_t           mids = (l->exactmaxs[0] + l->exactmins[0]) / 2;
        vec_t           midt = (l->exactmaxs[1] + l->exactmins[1]) / 2;

        SetSurfFromST(l, midpoint, mids, midt);

        if ((leaf_mid = PointInLeaf(midpoint)) != g_dleafs)
        {
            if ((leaf_mid->contents != CONTENTS_SKY) && (leaf_mid->contents != CONTENTS_SOLID))
            {
                return leaf_mid;
            }
        }

        VectorCopy(midpoint, broken_midpoint);
        //Verbose("Tried original method, failed at (%4.3f %4.3f %4.3f)\n", midpoint[0], midpoint[1], midpoint[2]);
    }

    VectorCopy(broken_midpoint, midpoint);
    return HuntForWorld(midpoint, face_delta, p, DEFAULT_HUNT_SIZE, DEFAULT_HUNT_SCALE, DEFAULT_HUNT_OFFSET);
}

// =====================================================================================
//  SimpleNudge
//      Return vec_t in point only valid when function returns true
//      Use negative scales to push away from center instead
// =====================================================================================
static bool     SimpleNudge(vec_t* const point, const lightinfo_t* const l, vec_t* const s, vec_t* const t, const vec_t delta)
{
    const int       facenum = l->surfnum;
    const dface_t*  f = g_dfaces + facenum;
    const dplane_t* p = getPlaneFromFace(f);
    const vec_t*    face_delta = g_face_offset[facenum];
    const int       h = l->texsize[1] + 1;
    const int       w = l->texsize[0] + 1;
    const vec_t     half_w = (vec_t)(w - 1) / 2.0;
    const vec_t     half_h = (vec_t)(h - 1) / 2.0;
    const vec_t     s_vec = *s;
    const vec_t     t_vec = *t;
    vec_t           s1;
    vec_t           t1;

    if (s_vec > half_w)
    {
        s1 = s_vec - delta;
    }
    else
    {
        s1 = s_vec + delta;
    }

    SetSurfFromST(l, point, s1, t_vec);
    if (HuntForWorld(point, face_delta, p, DEFAULT_HUNT_SIZE, DEFAULT_HUNT_SCALE, DEFAULT_HUNT_OFFSET))
    {
        *s = s1;
        return true;
    }

    if (t_vec > half_h)
    {
        t1 = t_vec - delta;
    }
    else
    {
        t1 = t_vec + delta;
    }

    SetSurfFromST(l, point, s_vec, t1);
    if (HuntForWorld(point, face_delta, p, DEFAULT_HUNT_SIZE, DEFAULT_HUNT_SCALE, DEFAULT_HUNT_OFFSET))
    {
        *t = t1;
        return true;
    }

    return false;
}
#endif

typedef enum
{
    LightOutside,                                          // Not lit
    LightShifted,                                          // used HuntForWorld on 100% dark face
    LightShiftedInside,                                    // moved to neighbhor on 2nd cleanup pass
    LightNormal,                                           // Normally lit with no movement
    LightPulledInside,                                     // Pulled inside by bleed code adjustments
    LightSimpleNudge,                                      // A simple nudge 1/3 or 2/3 towards center along S or T axist
#ifndef HLRAD_NUDGE_VL
    LightSimpleNudgeEmbedded                               // A nudge either 1 full unit in each of S and T axis, or 1/3 or 2/3 AWAY from center
#endif
}
light_flag_t;

// =====================================================================================
//  CalcPoints
//      For each texture aligned grid point, back project onto the plane
//      to get the world xyz value of the sample point
// =====================================================================================
#ifdef HLRAD_CalcPoints_NEW
static int		PointInFace(const lightinfo_t *l, const vec_t* point)
{
	int facenum = l->surfnum;
	const dface_t* f = &g_dfaces[facenum];
	Winding *w;
	dplane_t plane;
	VectorCopy (l->texnormal, plane.normal);
	const dplane_t *p = &plane;
	vec3_t new_point;
	VectorSubtract (point, g_face_offset[facenum], new_point);
	w = new Winding (*f);
	if (point_in_winding (*w, *p, new_point))
	{
		delete w;
		return facenum;
	}
	delete w;

	int j;
	for (j = 0; j < f->numedges; j++)
	{
		int e;
		edgeshare_t *es;
		dface_t* f2;
		e = g_dsurfedges[f->firstedge + j];
		es = &g_edgeshare[abs(e)];
		if (!es->smooth)
			continue;
		f2 = es->faces[!(e<0)];
		const dplane_t *p2 = getPlaneFromFace (f2);
		if (DotProduct (p->normal, p2->normal) < NORMAL_EPSILON)
			continue;
		w = new Winding (*f2);
		if (point_in_winding (*w, *p, new_point))
		{
			delete w;
			return f2 - g_dfaces;
		}
		delete w;
	}
#ifdef HLRAD_SMOOTH_FACELIST
	for (j = 0; j < f->numedges; j++)
	{
		int e;
		edgeshare_t *es;
		dface_t* f2;
		e = g_dsurfedges[f->firstedge + j];
		es = &g_edgeshare[abs(e)];
		if (!es->smooth)
			continue;
		for (int edgeend = 0; edgeend < 2; edgeend++)
		{
			for (facelist_t *l = es->vertex_facelist[edgeend]; l; l = l->next)
			{
				f2 = l->face;
				const dplane_t *p2 = getPlaneFromFace (f2);
				if (DotProduct (p->normal, p2->normal) < NORMAL_EPSILON)
					continue;
				w = new Winding (*f2);
				if (point_in_winding (*w, *p, new_point))
				{
					delete w;
					return f2 - g_dfaces;
				}
				delete w;
			}
		}
	}
#endif
	return facenum;
}
static void		SetSTFromSurf(const lightinfo_t* const l, const vec_t* surf, vec_t& s, vec_t& t)
{
    const int       facenum = l->surfnum;
    int             j;

	s = t = 0;
	for (j = 0; j < 3; j++)
	{
		s += (surf[j] - g_face_offset[facenum][j] - l->texorg[j]) * l->worldtotex[0][j];
		t += (surf[j] - g_face_offset[facenum][j] - l->texorg[j]) * l->worldtotex[1][j];
	}
}
static light_flag_t SetSampleFromST(vec_t* const point, const lightinfo_t* const l, const vec_t original_s, const vec_t original_t, eModelLightmodes lightmode)
{
	light_flag_t	LuxelFlag = LightOutside;
	int				huntsize = 3;
	vec_t			huntscale = 0.2;
	vec_t			width = DEFAULT_EDGE_WIDTH;

	int				facenum = l->surfnum;
	const vec_t*	face_delta = g_face_offset[facenum];

	const dface_t*	f = &g_dfaces[facenum];
	const dplane_t*	p = getPlaneFromFace (f);
	Winding			*wd = new Winding (*f);
	{
		int				j;
		for (j = 0; j < wd->m_NumPoints; j++)
		{
			VectorAdd (wd->m_Points[j], face_delta, wd->m_Points[j]);
		}
	}

	const vec_t*	face_centroid = g_face_centroids[facenum];
	vec_t			mids, midt;
	SetSTFromSurf (l, face_centroid, mids, midt);

	vec3_t			surf_original;
	dleaf_t*		leaf_original;
	SetSurfFromST (l, surf_original, original_s, original_t);
	leaf_original = HuntForWorld (surf_original, face_delta, p, 1, 0.0, DEFAULT_HUNT_OFFSET);

	int				facenum_tosnap = PointInFace (l, surf_original);
	const dface_t*	f_tosnap = &g_dfaces[facenum_tosnap];
	const dplane_t*	p_tosnap = getPlaneFromFace (f_tosnap);
#ifdef HLRAD_SMOOTH_TEXNORMAL
	vec3_t			snapdir;
	if (!GetIntertexnormal (facenum, facenum_tosnap, snapdir))
	{
		facenum_tosnap = facenum;
		f_tosnap = f;
		p_tosnap = p;
	}
#endif

	vec3_t			surf_direct;
	dleaf_t*		leaf_direct;
	VectorCopy (surf_original, surf_direct);
	{
		vec_t dist;
		vec_t scale;
#ifdef HLRAD_SMOOTH_TEXNORMAL
		scale = DotProduct (snapdir, p_tosnap->normal);
#else
		scale = DotProduct (l->texnormal, p_tosnap->normal);
#endif
		dist = DotProduct (surf_direct, p_tosnap->normal) - DotProduct (face_delta, p_tosnap->normal) - p_tosnap->dist - DEFAULT_HUNT_OFFSET;
#ifdef HLRAD_SMOOTH_TEXNORMAL
		VectorMA (surf_direct, - dist / scale, snapdir, surf_direct);
#else
		VectorMA (surf_direct, - dist / scale, l->texnormal, surf_direct);
#endif
	}
	leaf_direct = HuntForWorld (surf_direct, face_delta, p_tosnap, huntsize, huntscale, DEFAULT_HUNT_OFFSET);

	if (LuxelFlag == LightOutside)
	{
		if (leaf_direct && point_in_winding_noedge (*wd, *p, surf_direct, width))
		{
			LuxelFlag = LightNormal;
			VectorCopy (surf_direct, point);
		}
	}

	if (LuxelFlag == LightOutside)
	{
		bool	blocked_direct;
		bool	blocked_inwinding;
		bool	blocked_inwinding_noedge;
		vec3_t	surf_inwinding;
		vec3_t	surf_inwinding_noedge;
		dleaf_t*leaf_inwinding;
		dleaf_t*leaf_inwinding_noedge;
#ifdef HLRAD_HULLU
		vec3_t transparency = { 1.0, 1.0, 1.0 };
#endif
#ifdef HLRAD_OPAQUE_STYLE
		int opaquestyle;
#endif
		{
			blocked_direct = (leaf_direct == NULL);
			if (!point_in_winding (*wd, *p, surf_original))
			{
				VectorCopy (surf_original, surf_inwinding);
				snap_to_winding (*wd, *p, surf_inwinding);
				leaf_inwinding = HuntForWorld (surf_inwinding, face_delta, p, huntsize, huntscale, DEFAULT_HUNT_OFFSET);
				if ( blocked_direct
					|| !leaf_inwinding
					|| TestLine (surf_direct, surf_inwinding) != CONTENTS_EMPTY
					|| TestSegmentAgainstOpaqueList (surf_direct, surf_inwinding
	#ifdef HLRAD_HULLU
						, transparency
	#endif
	#ifdef HLRAD_OPAQUE_STYLE
						, opaquestyle
	#endif
						) == true
	#ifdef HLRAD_OPAQUE_STYLE
					|| opaquestyle != -1
	#endif
	#ifdef HLRAD_TRANSLUCENT
					|| l->translucent_b && facenum_tosnap == facenum
	#endif
					)
				{
					blocked_direct = true;
				}
			}
			else
			{
				VectorCopy (surf_original, surf_inwinding);
				leaf_inwinding = leaf_original;
			}
			blocked_inwinding = (leaf_inwinding == NULL);
			if (!point_in_winding_noedge (*wd, *p, surf_inwinding, width))
			{
				VectorCopy (surf_inwinding, surf_inwinding_noedge);
				snap_to_winding_noedge (*wd, *p, surf_inwinding_noedge, width);
				leaf_inwinding_noedge = HuntForWorld (surf_inwinding_noedge, face_delta, p, huntsize, huntscale, DEFAULT_HUNT_OFFSET);
				if ( blocked_inwinding
					|| !leaf_inwinding_noedge
					|| TestLine (surf_inwinding, surf_inwinding_noedge) != CONTENTS_EMPTY
					|| TestSegmentAgainstOpaqueList (surf_inwinding, surf_inwinding_noedge
	#ifdef HLRAD_HULLU
						, transparency
	#endif
	#ifdef HLRAD_OPAQUE_STYLE
						, opaquestyle
	#endif
						) == true
	#ifdef HLRAD_OPAQUE_STYLE
					|| opaquestyle != -1
	#endif
					)
				{
					blocked_inwinding = true;
				}
			}
			else
			{
				VectorCopy (surf_inwinding, surf_inwinding_noedge);
				leaf_inwinding_noedge = leaf_inwinding;
			}
			blocked_inwinding_noedge = (leaf_inwinding_noedge == NULL);
			if (blocked_inwinding_noedge == true)
			{
				blocked_inwinding = true;
			}
			if (blocked_inwinding == true)
			{
				blocked_direct = true;
			}
		}
		if (!blocked_direct)
		{
			LuxelFlag = LightNormal;
			VectorCopy (surf_direct, point);
		}
		else if (!blocked_inwinding)
		{
			LuxelFlag = LightPulledInside;
			VectorCopy (surf_inwinding, point);
		}
		else if (!blocked_inwinding_noedge)
		{
			LuxelFlag = LightPulledInside;
			VectorCopy (surf_inwinding_noedge, point);
		}
	}

	if (LuxelFlag == LightOutside)
	{
		const int numnudges = 13;
		vec_t nudgelist[numnudges][2] = {{0,0},{0.6,0},{0,0.6},{-0.6,0},{0,-0.6},{1.1,1.1},{1.1,-1.1},{-1.1,1.1},{-1.1,-1.1},{1.6,0},{0,1.6},{-1.6,0},{0,-1.6}};
		vec_t nudgescale_s, nudgescale_t;
		nudgescale_s = original_s <= mids? TEXTURE_STEP: -TEXTURE_STEP;
		nudgescale_t = original_t <= midt? TEXTURE_STEP: -TEXTURE_STEP;
		int i;
		for (i = 0; i < numnudges; i++)
		{
			vec_t s1 = original_s + nudgelist[i][0] * nudgescale_s;
			vec_t t1 = original_t + nudgelist[i][1] * nudgescale_t;
			vec3_t surf;
			SetSurfFromST(l, surf, s1, t1);
			if (point_in_winding (*wd, *p, surf) && HuntForWorld (surf, face_delta, p, huntsize, huntscale, DEFAULT_HUNT_OFFSET) && point_in_winding (*wd, *p, surf))
			{
				LuxelFlag = LightSimpleNudge;
				VectorCopy (surf, point);
				break;
			}
		}
	}

	if (LuxelFlag == LightOutside)
	{
		VectorCopy (surf_original, point);
	}

	delete wd;
	return LuxelFlag;
}
static void		CalcPoints(lightinfo_t* l)
{
	const int       facenum = l->surfnum;
	const dface_t*  f = g_dfaces + facenum;
	const dplane_t* p = getPlaneFromFace	(f);
	const vec_t*    face_delta = g_face_offset[facenum];
	const eModelLightmodes lightmode = g_face_lightmode[facenum];
	const int       h = l->texsize[1] + 1;
	const int       w = l->texsize[0] + 1;
	const vec_t     starts = (l->texmins[0] * 16);
	const vec_t     startt = (l->texmins[1] * 16);
	light_flag_t    LuxelFlags[MAX_SINGLEMAP];
	light_flag_t*   pLuxelFlags;
	vec_t           us, ut;
	vec_t*          surf;
	int             s, t;
	l->numsurfpt = w * h;
	for (t = 0; t < h; t++)
	{
		for (s = 0; s < w; s++)
		{
			surf = l->surfpt[s+w*t];
			pLuxelFlags = &LuxelFlags[s+w*t];
			us = starts + s * TEXTURE_STEP;
			ut = startt + t * TEXTURE_STEP;
			*pLuxelFlags = SetSampleFromST (surf, l, us, ut, lightmode);
		}
	}
#if 1
    {
		int i, n;
		int s_other, t_other;
		light_flag_t* pLuxelFlags_other;
		vec_t* surf_other;
		bool adjusted;
		for (i = 0; i < h + w; i++)
		{
			adjusted = false;
			for (t = 0; t < h; t++)
			{
				for (s = 0; s < w; s++)
				{
					surf = l->surfpt[s+w*t];
					pLuxelFlags = &LuxelFlags[s+w*t];
					if (*pLuxelFlags != LightOutside)
						continue;
					for (n = 0; n < 4; n++)
					{
						switch (n)
						{
						case 0: s_other = s + 1; t_other = t; break;
						case 1: s_other = s - 1; t_other = t; break;
						case 2: s_other = s; t_other = t + 1; break;
						case 3: s_other = s; t_other = t - 1; break;
						}
						if (t_other < 0 || t_other >= h || s_other < 0 || s_other >= w)
							continue;
						surf_other = l->surfpt[s_other+w*t_other];
						pLuxelFlags_other = &LuxelFlags[s_other+w*t_other];
						if (*pLuxelFlags_other != LightOutside)
						{
							*pLuxelFlags = LightShiftedInside;
							VectorCopy (surf_other, surf);
							adjusted = true;
							break;
						}
					}
				}
			}
			if (!adjusted)
				break;
		}
	}
#endif
}
#else
static void     CalcPoints(lightinfo_t* l)
{
    const int       facenum = l->surfnum;
    const dface_t*  f = g_dfaces + facenum;
    const dplane_t* p = getPlaneFromFace	(f);

    const vec_t*    face_delta = g_face_offset[facenum];
    const eModelLightmodes lightmode = g_face_lightmode[facenum];

#ifdef HLRAD_NUDGE_VL
	vec_t mids, midt;
	{
		// use winding center instead
		vec3_t surf;
		VectorSubtract (g_face_centroids[facenum], g_face_offset[facenum], surf);
		VectorSubtract (surf, l->texorg, surf);
		mids = DotProduct (surf, l->worldtotex[0]);
		midt = DotProduct (surf, l->worldtotex[1]);
	}
#else
    const vec_t     mids = (l->exactmaxs[0] + l->exactmins[0]) / 2;
    const vec_t     midt = (l->exactmaxs[1] + l->exactmins[1]) / 2;
#endif

    const int       h = l->texsize[1] + 1;
    const int       w = l->texsize[0] + 1;

    const vec_t     starts = (l->texmins[0] * 16);
    const vec_t     startt = (l->texmins[1] * 16);

    light_flag_t    LuxelFlags[MAX_SINGLEMAP];
    light_flag_t*   pLuxelFlags;
    vec_t           us, ut;
    vec_t*          surf;
    vec3_t          surface_midpoint;
    dleaf_t*        leaf_mid;
    dleaf_t*        leaf_surf;
    int             s, t;
    int             i;

    l->numsurfpt = w * h;

    memset(LuxelFlags, 0, sizeof(LuxelFlags));

    leaf_mid = FindSurfaceMidpoint(l, surface_midpoint);
#if 0
    if (!leaf_mid)
    {
        Developer(DEVELOPER_LEVEL_FLUFF, "CalcPoints [face %d] (%4.3f %4.3f %4.3f) midpoint outside world\n",
                  facenum, surface_midpoint[0], surface_midpoint[1], surface_midpoint[2]);
    }
    else
    {
        Developer(DEVELOPER_LEVEL_FLUFF, "FindSurfaceMidpoint [face %d] @ (%4.3f %4.3f %4.3f)\n",
                  facenum, surface_midpoint[0], surface_midpoint[1], surface_midpoint[2]);
    }
#endif

    // First pass, light normally, and pull any faces toward the center for bleed adjustment

    surf = l->surfpt[0];
    pLuxelFlags = LuxelFlags;
    for (t = 0; t < h; t++)
    {
        for (s = 0; s < w; s++, surf += 3, pLuxelFlags++)
        {
            vec_t           original_s = us = starts + s * TEXTURE_STEP;
            vec_t           original_t = ut = startt + t * TEXTURE_STEP;

#ifdef HLRAD_AddSampleToPatch_PRECISE
			SetSurfFromST (l, l->surfpt_original[s+w*t], original_s, original_t);
#endif
            SetSurfFromST(l, surf, us, ut);
            leaf_surf = HuntForWorld(surf, face_delta, p, DEFAULT_HUNT_SIZE, DEFAULT_HUNT_SCALE, DEFAULT_HUNT_OFFSET);

            if (!leaf_surf)
            {
                // At first try a 1/3 and 2/3 distance to nearest in each S and T axis towards the face midpoint
                if (SimpleNudge(surf, l, &us, &ut, TEXTURE_STEP * (1.0 / 3.0)))
                {
                    *pLuxelFlags = LightSimpleNudge;
                }
                else if (SimpleNudge(surf, l, &us, &ut, -TEXTURE_STEP * (1.0 / 3.0)))
                {
                    *pLuxelFlags = LightSimpleNudge;
                }
                else if (SimpleNudge(surf, l, &us, &ut, TEXTURE_STEP * (2.0 / 3.0)))
                {
                    *pLuxelFlags = LightSimpleNudge;
                }
                else if (SimpleNudge(surf, l, &us, &ut, -TEXTURE_STEP * (2.0 / 3.0)))
                {
                    *pLuxelFlags = LightSimpleNudge;
                }
#ifdef HLRAD_NUDGE_VL
                else if (SimpleNudge(surf, l, &us, &ut, TEXTURE_STEP))
                {
                    *pLuxelFlags = LightSimpleNudge;
                }
                else if (SimpleNudge(surf, l, &us, &ut, -TEXTURE_STEP))
                {
                    *pLuxelFlags = LightSimpleNudge;
                }
#else
                // Next, if this is a model flagged with the 'Embedded' mode, try away from the facemid too
                else if (lightmode & eModelLightmodeEmbedded)
                {
                    SetSurfFromST(l, surf, us, ut);
                    if (SimpleNudge(surf, l, &us, &ut, TEXTURE_STEP))
                    {
                        *pLuxelFlags = LightSimpleNudgeEmbedded;
                        continue;
                    }
                    if (SimpleNudge(surf, l, &us, &ut, -TEXTURE_STEP))
                    {
                        *pLuxelFlags = LightSimpleNudgeEmbedded;
                        continue;
                    }

                    SetSurfFromST(l, surf, original_s, original_t);
                    *pLuxelFlags = LightOutside;
                    continue;
                }
#endif
            }

#ifndef HLRAD_NUDGE_VL
            if (!(lightmode & eModelLightmodeEmbedded))
#endif
            {
#ifdef HLRAD_NUDGE_VL
				// HLRAD_NUDGE_VL: only pull when light is blocked AND point is outside face.
				vec3_t			surf_nopull;
				vec_t			us_nopull = us, ut_nopull = ut;
				Winding			*wd = new Winding (*f);
				int				j;
				for (j = 0; j < wd->m_NumPoints; j++)
				{
					VectorAdd (wd->m_Points[j], face_delta, wd->m_Points[j]);
				}
#endif
#ifdef HLRAD_SNAPTOWINDING
				bool nudge_succeeded = false;
				SetSurfFromST(l, surf, us, ut);
				leaf_surf = HuntForWorld(surf, face_delta, p, DEFAULT_HUNT_SIZE, DEFAULT_HUNT_SCALE, DEFAULT_HUNT_OFFSET);
				if (leaf_surf && point_in_winding_noedge (*wd, *p, surf, 1.0))
				{
					*pLuxelFlags = LightNormal;
					nudge_succeeded = true;
				}
				else
				{
					SetSurfFromST(l, surf, us, ut);
					snap_to_winding (*wd, *p, surf);
					if (lightmode & eModelLightmodeConcave)
					{
						VectorScale (surf, 0.99, surf);
						VectorMA (surf, 0.01, g_face_centroids[facenum], surf);
					}
					leaf_surf = HuntForWorld(surf, face_delta, p, DEFAULT_HUNT_SIZE, DEFAULT_HUNT_SCALE, DEFAULT_HUNT_OFFSET);
					if (leaf_surf)
					{
						*pLuxelFlags = LightPulledInside;
						nudge_succeeded = true;
					}
				}

#else
				// Pull the sample points towards the facemid if visibility is blocked
				// and the facemid is inside the world
	#ifdef HLRAD_NUDGE_SMALLSTEP
				int             nudge_divisor = 4 * max(max(w, h), 4);
	#else
				int             nudge_divisor = max(max(w, h), 4);
	#endif
				int             max_nudge = nudge_divisor + 1;
				bool            nudge_succeeded = false;

				vec_t           nudge_s = (mids - us) / (vec_t)nudge_divisor;
				vec_t           nudge_t = (midt - ut) / (vec_t)nudge_divisor;

				// if a line can be traced from surf to facemid, the point is good
				for (i = 0; i < max_nudge; i++)
				{
					// Make sure we are "in the world"(Not the zero leaf)
	#ifndef HLRAD_NUDGE_VL
					if (leaf_mid)
					{
	#endif
						SetSurfFromST(l, surf, us, ut);
						leaf_surf =
							HuntForWorld(surf, face_delta, p, DEFAULT_HUNT_SIZE, DEFAULT_HUNT_SCALE,
										 DEFAULT_HUNT_OFFSET);

						if (leaf_surf)
						{
	#ifdef HLRAD_NUDGE_VL
							if (point_in_winding_noedge (*wd, *p, surf, 1.0))
							{
	#else
							if (TestLine(surface_midpoint, surf) == CONTENTS_EMPTY)
							{
								if (lightmode & eModelLightmodeConcave)
								{
		#ifdef HLRAD_HULLU
									vec3_t transparency = { 1.0, 1.0, 1.0 };
		#endif
		#ifdef HLRAD_OPAQUE_STYLE
									int opaquestyle;
		#endif
									if (TestSegmentAgainstOpaqueList(surface_midpoint, surf
		#ifdef HLRAD_HULLU
										, transparency
		#endif
		#ifdef HLRAD_OPAQUE_STYLE
										, opaquestyle
		#endif
										)
		#ifdef HLRAD_OPAQUE_STYLE
										|| opaquestyle != -1
		#endif
										)
									{
										Log("SDF::4\n");
										us += nudge_s;
										ut += nudge_t;
										continue;   // Try nudge again, we hit an opaque face
									}
								}
	#endif
								if (i)
								{
									*pLuxelFlags = LightPulledInside;
								}
								else
								{
									*pLuxelFlags = LightNormal;
								}
								nudge_succeeded = true;
								break;
							}
						}
	#ifndef HLRAD_NUDGE_VL
					}
					else
					{
						leaf_surf = PointInLeaf(surf);
						if (leaf_surf != g_dleafs)
						{
							if ((leaf_surf->contents != CONTENTS_SKY) && (leaf_surf->contents != CONTENTS_SOLID))
							{
								*pLuxelFlags = LightNormal;
								nudge_succeeded = true;
								break;
							}
						}
					}
	#endif

					us += nudge_s;
					ut += nudge_t;
				}
#endif /*HLRAD_SNAPTOWINDING*/

                if (!nudge_succeeded)
                {
                    SetSurfFromST(l, surf, original_s, original_t);
                    *pLuxelFlags = LightOutside;
                }
#ifdef HLRAD_NUDGE_VL
				delete wd;
				if (*pLuxelFlags == LightPulledInside)
				{
					SetSurfFromST(l, surf_nopull, us_nopull, ut_nopull);
					leaf_surf =
						HuntForWorld(surf_nopull, face_delta, p, DEFAULT_HUNT_SIZE, DEFAULT_HUNT_SCALE,
							DEFAULT_HUNT_OFFSET);
					if (leaf_surf)
					{
						if (TestLine(surf, surf_nopull) == CONTENTS_EMPTY)
						{
#ifdef HLRAD_HULLU
							vec3_t transparency = { 1.0, 1.0, 1.0 };
#endif
#ifdef HLRAD_OPAQUE_STYLE
							int opaquestyle;
#endif
							if (!TestSegmentAgainstOpaqueList(surf, surf_nopull
#ifdef HLRAD_HULLU
								, transparency
#endif
#ifdef HLRAD_OPAQUE_STYLE
								, opaquestyle
#endif
								)
#ifdef HLRAD_OPAQUE_STYLE
								&& opaquestyle == -1
#endif
								)
							{
								*pLuxelFlags = LightNormal;
								VectorCopy (surf_nopull, surf);
							}
						}
					}
				}
#endif
            }
        }
    }

    // 2nd Pass, find units that are not lit and try to move them one half or unit worth 
    // in each direction and see if that is lit.
    // This handles 1 x N lightmaps which are all dark everywhere and have no frame of refernece
    // for a good center or directly lit areas
    surf = l->surfpt[0];
    pLuxelFlags = LuxelFlags;
#if 0
    Developer(DEVELOPER_LEVEL_SPAM,
              "w (%d) h (%d) dim (%d) leafmid (%4.3f %4.3f %4.3f) plane normal (%4.3f) (%4.3f) (%4.3f) dist (%f)\n", w,
              h, w * h, surface_midpoint[0], surface_midpoint[1], surface_midpoint[2], p->normal[0], p->normal[1],
              p->normal[2], p->dist);
#endif
    {
        int             total_dark = 0;
        int             total_adjusted = 0;

        for (t = 0; t < h; t++)
        {
            for (s = 0; s < w; s++, surf += 3, pLuxelFlags++)
            {
                if (!*pLuxelFlags)
                {
#if 0
                    Developer(DEVELOPER_LEVEL_FLUFF, "Dark (%d %d) (%4.3f %4.3f %4.3f)\n",
                              s, t, surf[0], surf[1], surf[2]);
#endif
                    total_dark++;
                    if (HuntForWorld(surf, face_delta, p, DEFAULT_HUNT_SIZE, DEFAULT_HUNT_SCALE, DEFAULT_HUNT_OFFSET))
                    {
#if 0
                        Developer(DEVELOPER_LEVEL_FLUFF, "Shifted %d %d to (%4.3f %4.3f %4.3f)\n", s, t, surf[0],
                                  surf[1], surf[2]);
#endif
                        *pLuxelFlags = LightShifted;
                        total_adjusted++;
                    }
                    else if (HuntForWorld(surf, face_delta, p, 101, 0.5, DEFAULT_HUNT_OFFSET))
                    {
#if 0
                        Developer(DEVELOPER_LEVEL_FLUFF, "Shifted %d %d to (%4.3f %4.3f %4.3f)\n", s, t, surf[0],
                                  surf[1], surf[2]);
#endif
                        *pLuxelFlags = LightShifted;
                        total_adjusted++;
                    }
                }
            }
        }
#if 0
        if (total_dark)
        {
            Developer(DEVELOPER_LEVEL_FLUFF, "Pass 2 : %d dark, %d corrected\n", total_dark, total_adjusted);
        }
#endif
    }

    // 3rd Pass, find units that are not lit and move them towards neighbhors who are
    // Currently finds the first lit neighbhor and uses its data
    surf = l->surfpt[0];
    pLuxelFlags = LuxelFlags;
    {
        int             total_dark = 0;
        int             total_adjusted = 0;

        for (t = 0; t < h; t++)
        {
            for (s = 0; s < w; s++, surf += 3, pLuxelFlags++)
            {
                if (!*pLuxelFlags)
                {
                    int             x_min = max(0, s - 1);
                    int             x_max = min(w, s + 1);
                    int             y_min = max(0, t - 1);
                    int             y_max = min(t, t + 1);

                    int             x, y;

#if 0
                    Developer(DEVELOPER_LEVEL_FLUFF, "Point outside (%d %d) (%4.3f %4.3f %4.3f)\n",
                              s, t, surf[0], surf[1], surf[2]);
#endif

                    total_dark++;

                    for (x = x_min; x < x_max; x++)
                    {
                        for (y = y_min; y < y_max; y++)
                        {
                            if (*pLuxelFlags >= LightNormal)
                            {
                                dleaf_t*        leaf;
                                vec_t*          other_surf = l->surfpt[0];

                                other_surf += ((y * w) + x) * 3;

                                leaf = PointInLeaf(other_surf);
                                if ((leaf->contents != CONTENTS_SKY && leaf->contents != CONTENTS_SOLID))
                                {
                                    *pLuxelFlags = LightShiftedInside;
#if 0
                                    Developer(DEVELOPER_LEVEL_MESSAGE,
                                              "Nudged (%d %d) (%4.3f %4.3f %4.3f) to (%d %d) (%4.3f %4.3f %4.3f) \n",
                                              s, t, surf[0], surf[1], surf[2], x, y, other_surf[0], other_surf[1],
                                              other_surf[2]);
#endif
                                    VectorCopy(other_surf, surf);
                                    total_adjusted++;
                                    goto found_it;
                                }
                            }
                        }
                    }
                }
              found_it:;
            }
        }
#if 0
        if (total_dark)
        {
            Developer(DEVELOPER_LEVEL_FLUFF, "Pass 2 : %d dark, %d corrected\n", total_dark, total_adjusted);
        }
#endif
    }
}
#endif

//==============================================================

typedef struct
{
    vec3_t          pos;
    vec3_t          light;
}
sample_t;

typedef struct
{
    int             numsamples;
    sample_t*       samples[MAXLIGHTMAPS];
}
facelight_t;

static directlight_t* directlights[MAX_MAP_LEAFS];
static facelight_t facelight[MAX_MAP_FACES];
static int      numdlights;

#ifndef HLRAD_REFLECTIVITY
#define	DIRECT_SCALE	0.1f
#endif

// =====================================================================================
//  CreateDirectLights
// =====================================================================================
void            CreateDirectLights()
{
    unsigned        i;
    patch_t*        p;
    directlight_t*  dl;
    dleaf_t*        leaf;
    int             leafnum;
    entity_t*       e;
    entity_t*       e2;
    const char*     name;
    const char*     target;
    float           angle;
    vec3_t          dest;

    // AJM: coplaner lighting
    vec3_t          temp_normal;

    numdlights = 0;
#ifdef HLRAD_STYLEREPORT
	int styleused[ALLSTYLES];
	memset (styleused, 0, ALLSTYLES * sizeof(styleused[0]));
	styleused[0] = true;
	int numstyles = 1;
#endif

    //
    // surfaces
    //
    for (i = 0, p = g_patches; i < g_num_patches; i++, p++)
    {
#ifdef ZHLT_TEXLIGHT
#ifdef HLRAD_STYLEREPORT
		if (p->emitstyle >= 0 && p->emitstyle < ALLSTYLES)
		{
			if (styleused[p->emitstyle] == false)
			{
				styleused[p->emitstyle] = true;
				numstyles++;
			}
		}
#endif
        if (
	#ifdef HLRAD_REFLECTIVITY
			DotProduct (p->baselight, p->texturereflectivity) / 3
	#else
			VectorAvg(p->baselight)
	#endif
	#ifdef HLRAD_TEXLIGHTTHRESHOLD_FIX
			> 0.0
	#else
			>= g_dlight_threshold
	#endif
	#ifdef HLRAD_CUSTOMTEXLIGHT
			&& !(g_face_texlights[p->faceNumber]
				&& *ValueForKey (g_face_texlights[p->faceNumber], "_scale")
				&& FloatForKey (g_face_texlights[p->faceNumber], "_scale") <= 0)
	#endif
			) //LRC
#else
        if (
	#ifdef HLRAD_REFLECTIVITY
			DotProduct (p->totallight, p->texturereflectivity) / 3
	#else
			VectorAvg(p->totallight)
	#endif
			>= g_dlight_threshold
			)
#endif
        {
            numdlights++;
            dl = (directlight_t*)calloc(1, sizeof(directlight_t));
#ifdef HLRAD_HLASSUMENOMEMORY
			hlassume (dl != NULL, assume_NoMemory);
#endif

            VectorCopy(p->origin, dl->origin);

            leaf = PointInLeaf(dl->origin);
            leafnum = leaf - g_dleafs;

            dl->next = directlights[leafnum];
            directlights[leafnum] = dl;
#ifdef ZHLT_TEXLIGHT
            dl->style = p->emitstyle; //LRC
#endif
#ifdef HLRAD_GatherPatchLight
			dl->topatch = false;
	#ifdef HLRAD_TEXLIGHTTHRESHOLD_FIX
			if (!p->emitmode)
			{
				dl->topatch = true;
			}
	#endif
#endif
#ifdef HLRAD_CUSTOMTEXLIGHT
			dl->stopdot = 0.0;
			dl->stopdot2 = 0.0;
			if (g_face_texlights[p->faceNumber])
			{
				if (*ValueForKey (g_face_texlights[p->faceNumber], "_cone"))
				{
					dl->stopdot = FloatForKey (g_face_texlights[p->faceNumber], "_cone");
					dl->stopdot = dl->stopdot >= 90? 0: (float)cos (dl->stopdot / 180 * Q_PI);
				}
				if (*ValueForKey (g_face_texlights[p->faceNumber], "_cone2"))
				{
					dl->stopdot2 = FloatForKey (g_face_texlights[p->faceNumber], "_cone2");
					dl->stopdot2 = dl->stopdot2 >= 90? 0: (float)cos (dl->stopdot2 / 180 * Q_PI);
				}
				if (dl->stopdot2 > dl->stopdot)
					dl->stopdot2 = dl->stopdot;
			}
#endif

            dl->type = emit_surface;
            VectorCopy(getPlaneFromFaceNumber(p->faceNumber)->normal, dl->normal);
#ifdef ZHLT_TEXLIGHT
            VectorCopy(p->baselight, dl->intensity); //LRC
#else
            VectorCopy(p->totallight, dl->intensity);
#endif
#ifdef HLRAD_CUSTOMTEXLIGHT
			if (g_face_texlights[p->faceNumber])
			{
				if (*ValueForKey (g_face_texlights[p->faceNumber], "_scale"))
				{
					vec_t scale = FloatForKey (g_face_texlights[p->faceNumber], "_scale");
					VectorScale (dl->intensity, scale, dl->intensity);
				}
			}
#endif
            VectorScale(dl->intensity, p->area, dl->intensity);
#ifdef HLRAD_REFLECTIVITY
			VectorScale (dl->intensity, 1.0 / Q_PI, dl->intensity);
			VectorMultiply (dl->intensity, p->texturereflectivity, dl->intensity);
#else
            VectorScale(dl->intensity, DIRECT_SCALE, dl->intensity);
#endif
        
#ifdef HLRAD_WATERBACKFACE_FIX
			dface_t *f = &g_dfaces[p->faceNumber];
			if (g_face_entity[p->faceNumber] - g_entities != 0 && !strncasecmp (GetTextureByNumber (f->texinfo), "!", 1))
			{
				directlight_t *dl2;
				numdlights++;
				dl2 = (directlight_t *)calloc (1, sizeof (directlight_t));
				hlassume (dl2 != NULL, assume_NoMemory);
				*dl2 = *dl;
				VectorMA (dl->origin, -2, dl->normal, dl2->origin);
				VectorSubtract (vec3_origin, dl->normal, dl2->normal);
				leaf = PointInLeaf (dl2->origin);
				leafnum = leaf - g_dleafs;
				dl2->next = directlights[leafnum];
				directlights[leafnum] = dl2;
			}
#endif
#ifndef HLRAD_CUSTOMTEXLIGHT // no softlight hack
        	// --------------------------------------------------------------
	        // Changes by Adam Foster - afoster@compsoc.man.ac.uk
	        // mazemaster's l33t backwards lighting (I still haven't a clue
	        // what it's supposed to be for) :-)
#ifdef HLRAD_WHOME

	        if (g_softlight_hack[0] || g_softlight_hack[1] || g_softlight_hack[2]) 
            {
		        numdlights++;
		        dl = (directlight_t *) calloc(1, sizeof(directlight_t));
#ifdef HLRAD_HLASSUMENOMEMORY
				hlassume (dl != NULL, assume_NoMemory);
#endif

		        VectorCopy(p->origin, dl->origin);

		        leaf = PointInLeaf(dl->origin);
		        leafnum = leaf - g_dleafs;

		        dl->next = directlights[leafnum];
		        directlights[leafnum] = dl;

#ifdef HLRAD_GatherPatchLight
				dl->topatch = false;
	#ifdef HLRAD_TEXLIGHTTHRESHOLD_FIX
				if (!p->emitmode)
				{
					dl->topatch = true;
				}
	#endif
#endif
		        dl->type = emit_surface;
		        VectorCopy(getPlaneFromFaceNumber(p->faceNumber)->normal, dl->normal);
		        VectorScale(dl->normal, g_softlight_hack_distance, temp_normal);
		        VectorAdd(dl->origin, temp_normal, dl->origin);
		        VectorScale(dl->normal, -1, dl->normal);

#ifdef ZHLT_TEXLIGHT
                VectorCopy(p->baselight, dl->intensity); //LRC
#else
		        VectorCopy(p->totallight, dl->intensity);
#endif
		        VectorScale(dl->intensity, p->area, dl->intensity);
#ifdef HLRAD_REFLECTIVITY
				VectorScale (dl->intensity, 1.0 / Q_PI, dl->intensity);
				VectorMultiply (dl->intensity, p->texturereflectivity, dl->intensity);
#else
		        VectorScale(dl->intensity, DIRECT_SCALE, dl->intensity);
#endif

		        dl->intensity[0] *= g_softlight_hack[0];
		        dl->intensity[1] *= g_softlight_hack[1];
		        dl->intensity[2] *= g_softlight_hack[2];
	        }

#endif
	        // --------------------------------------------------------------
#endif
        }

#ifdef ZHLT_TEXLIGHT
        //LRC        VectorClear(p->totallight[0]);
#else
        VectorClear(p->totallight);
#endif
    }

    //
    // entities
    //
    for (i = 0; i < (unsigned)g_numentities; i++)
    {
        const char*     pLight;
        double          r, g, b, scaler;
        float           l1;
        int             argCnt;

        e = &g_entities[i];
        name = ValueForKey(e, "classname");
        if (strncmp(name, "light", 5))
            continue;
#ifdef HLRAD_STYLE_CORING
		{
			int style = IntForKey (e, "style");
	#ifdef ZHLT_TEXLIGHT
			if (style < 0)
			{
				style = -style;
			}
	#endif
			style = (unsigned char)style;
			if (style > 0 && style < ALLSTYLES && *ValueForKey (e, "zhlt_stylecoring"))
			{
				g_corings[style] = FloatForKey (e, "zhlt_stylecoring");
			}
		}
#endif
#ifdef HLRAD_OPAQUE_STYLE
		if (!strcmp (ValueForKey (e, "classname"), "light_shadow"))
		{
#ifdef HLRAD_STYLEREPORT
			int style = IntForKey (e, "style");
	#ifdef ZHLT_TEXLIGHT
			if (style < 0)
			{
				style = -style;
			}
	#endif
			style = (unsigned char)style;
			if (style >= 0 && style < ALLSTYLES)
			{
				if (styleused[style] == false)
				{
					styleused[style] = true;
					numstyles++;
				}
			}
#endif
			continue;
		}
#endif
#ifdef HLRAD_CUSTOMTEXLIGHT
		if (!strcmp (ValueForKey (e, "classname"), "light_surface"))
		{
			continue;
		}
#endif

        numdlights++;
        dl = (directlight_t*)calloc(1, sizeof(directlight_t));
#ifdef HLRAD_HLASSUMENOMEMORY
		hlassume (dl != NULL, assume_NoMemory);
#endif

        GetVectorForKey(e, "origin", dl->origin);

        leaf = PointInLeaf(dl->origin);
        leafnum = leaf - g_dleafs;

        dl->next = directlights[leafnum];
        directlights[leafnum] = dl;

        dl->style = IntForKey(e, "style");
#ifdef ZHLT_TEXLIGHT
        if (dl->style < 0) 
            dl->style = -dl->style; //LRC
#endif
#ifdef HLRAD_STYLE_CORING
		dl->style = (unsigned char)dl->style;
		if (dl->style >= ALLSTYLES)
		{
			Error ("invalid light style: style (%d) >= ALLSTYLES (%d)", dl->style, ALLSTYLES);
		}
#endif
#ifdef HLRAD_STYLEREPORT
		if (dl->style >= 0 && dl->style < ALLSTYLES)
		{
			if (styleused[dl->style] == false)
			{
				styleused[dl->style] = true;
				numstyles++;
			}
		}
#endif
#ifdef HLRAD_GatherPatchLight
		dl->topatch = false;
		if (IntForKey (e, "_fast") == 1)
		{
			dl->topatch = true;
		}
#endif
        pLight = ValueForKey(e, "_light");
        // scanf into doubles, then assign, so it is vec_t size independent
        r = g = b = scaler = 0;
        argCnt = sscanf(pLight, "%lf %lf %lf %lf", &r, &g, &b, &scaler);
        dl->intensity[0] = (float)r;
        if (argCnt == 1)
        {
            // The R,G,B values are all equal.
            dl->intensity[1] = dl->intensity[2] = (float)r;
        }
        else if (argCnt == 3 || argCnt == 4)
        {
            // Save the other two G,B values.
            dl->intensity[1] = (float)g;
            dl->intensity[2] = (float)b;

            // Did we also get an "intensity" scaler value too?
            if (argCnt == 4)
            {
                // Scale the normalized 0-255 R,G,B values by the intensity scaler
                dl->intensity[0] = dl->intensity[0] / 255 * (float)scaler;
                dl->intensity[1] = dl->intensity[1] / 255 * (float)scaler;
                dl->intensity[2] = dl->intensity[2] / 255 * (float)scaler;
            }
        }
        else
        {
            Log("light at (%f,%f,%f) has bad or missing '_light' value : '%s'\n",
                dl->origin[0], dl->origin[1], dl->origin[2], pLight);
            continue;
        }

        dl->fade = FloatForKey(e, "_fade");
        if (dl->fade == 0.0)
        {
            dl->fade = g_fade;
        }

#ifndef HLRAD_ARG_MISC
        dl->falloff = IntForKey(e, "_falloff");
        if (dl->falloff == 0)
        {
            dl->falloff = g_falloff;
        }
#endif

        target = ValueForKey(e, "target");

        if (!strcmp(name, "light_spot") || !strcmp(name, "light_environment") || target[0])
        {
            if (!VectorAvg(dl->intensity))
            {
#ifndef HLRAD_ALLOWZEROBRIGHTNESS
                VectorFill(dl->intensity, 500);
#endif
            }
            dl->type = emit_spotlight;
            dl->stopdot = FloatForKey(e, "_cone");
            if (!dl->stopdot)
            {
                dl->stopdot = 10;
            }
            dl->stopdot2 = FloatForKey(e, "_cone2");
            if (!dl->stopdot2)
            {
                dl->stopdot2 = dl->stopdot;
            }
            if (dl->stopdot2 < dl->stopdot)
            {
                dl->stopdot2 = dl->stopdot;
            }
            dl->stopdot2 = (float)cos(dl->stopdot2 / 180 * Q_PI);
            dl->stopdot = (float)cos(dl->stopdot / 180 * Q_PI);

            if (!FindTargetEntity(target)) //--vluzacn
            {
                Warning("light at (%i %i %i) has missing target",
                        (int)dl->origin[0], (int)dl->origin[1], (int)dl->origin[2]);
				target = "";
            }
            if (target[0])
            {                                              // point towards target
                e2 = FindTargetEntity(target);
                if (!e2)
                {
                    Warning("light at (%i %i %i) has missing target",
                            (int)dl->origin[0], (int)dl->origin[1], (int)dl->origin[2]);
                }
                else
                {
                    GetVectorForKey(e2, "origin", dest);
                    VectorSubtract(dest, dl->origin, dl->normal);
                    VectorNormalize(dl->normal);
                }
            }
            else
            {                                              // point down angle
                vec3_t          vAngles;

                GetVectorForKey(e, "angles", vAngles);

                angle = (float)FloatForKey(e, "angle");
                if (angle == ANGLE_UP)
                {
                    dl->normal[0] = dl->normal[1] = 0;
                    dl->normal[2] = 1;
                }
                else if (angle == ANGLE_DOWN)
                {
                    dl->normal[0] = dl->normal[1] = 0;
                    dl->normal[2] = -1;
                }
                else
                {
                    // if we don't have a specific "angle" use the "angles" YAW
                    if (!angle)
                    {
                        angle = vAngles[1];
                    }

                    dl->normal[2] = 0;
                    dl->normal[0] = (float)cos(angle / 180 * Q_PI);
                    dl->normal[1] = (float)sin(angle / 180 * Q_PI);
                }

                angle = FloatForKey(e, "pitch");
                if (!angle)
                {
                    // if we don't have a specific "pitch" use the "angles" PITCH
                    angle = vAngles[0];
                }

                dl->normal[2] = (float)sin(angle / 180 * Q_PI);
                dl->normal[0] *= (float)cos(angle / 180 * Q_PI);
                dl->normal[1] *= (float)cos(angle / 180 * Q_PI);
            }

            if (FloatForKey(e, "_sky") || !strcmp(name, "light_environment"))
            {
				// -----------------------------------------------------------------------------------
				// Changes by Adam Foster - afoster@compsoc.man.ac.uk
				// diffuse lighting hack - most of the following code nicked from earlier
				// need to get diffuse intensity from new _diffuse_light key
				//
				// What does _sky do for spotlights, anyway?
				// -----------------------------------------------------------------------------------
#ifdef HLRAD_WHOME
				pLight = ValueForKey(e, "_diffuse_light");
        		r = g = b = scaler = 0;
        		argCnt = sscanf(pLight, "%lf %lf %lf %lf", &r, &g, &b, &scaler);
        		dl->diffuse_intensity[0] = (float)r;
        		if (argCnt == 1)
        		{
            		// The R,G,B values are all equal.
            		dl->diffuse_intensity[1] = dl->diffuse_intensity[2] = (float)r;
        		}
        		else if (argCnt == 3 || argCnt == 4)
        		{
            		// Save the other two G,B values.
            		dl->diffuse_intensity[1] = (float)g;
            		dl->diffuse_intensity[2] = (float)b;

            		// Did we also get an "intensity" scaler value too?
      		    	if (argCnt == 4)
     	       		{
                		// Scale the normalized 0-255 R,G,B values by the intensity scaler
                		dl->diffuse_intensity[0] = dl->diffuse_intensity[0] / 255 * (float)scaler;
                		dl->diffuse_intensity[1] = dl->diffuse_intensity[1] / 255 * (float)scaler;
                		dl->diffuse_intensity[2] = dl->diffuse_intensity[2] / 255 * (float)scaler;
            		}
        		}
				else
        		{
					// backwards compatibility with maps without _diffuse_light

					dl->diffuse_intensity[0] = dl->intensity[0];
					dl->diffuse_intensity[1] = dl->intensity[1];
					dl->diffuse_intensity[2] = dl->intensity[2];
        		}
#endif
				// -----------------------------------------------------------------------------------

                dl->type = emit_skylight;
                dl->stopdot2 = FloatForKey(e, "_sky");     // hack stopdot2 to a sky key number
            }
        }
        else
        {
            if (!VectorAvg(dl->intensity))
			{
#ifndef HLRAD_ALLOWZEROBRIGHTNESS
                VectorFill(dl->intensity, 300);
#endif
			}
            dl->type = emit_point;
        }

        if (dl->type != emit_skylight)
        {
			//why? --vluzacn
            l1 = max(dl->intensity[0], max(dl->intensity[1], dl->intensity[2]));
            l1 = l1 * l1 / 10;

            dl->intensity[0] *= l1;
            dl->intensity[1] *= l1;
            dl->intensity[2] *= l1;
        }
    }

    hlassume(numdlights, assume_NoLights);
    Log("%i direct lights\n", numdlights);
#ifdef HLRAD_STYLEREPORT
	Log("%i light styles\n", numstyles);
#endif
}

// =====================================================================================
//  DeleteDirectLights
// =====================================================================================
void            DeleteDirectLights()
{
    int             l;
    directlight_t*  dl;

#ifdef HLRAD_VIS_FIX
	for (l = 0; l < 1 + g_dmodels[0].visleafs; l++)
#else
    for (l = 0; l < g_numleafs; l++)
#endif
    {
        dl = directlights[l];
        while (dl)
        {
            directlights[l] = dl->next;
            free(dl);
            dl = directlights[l];
        }
    }

    // AJM: todo: strip light entities out at this point
}

// =====================================================================================
//  GatherSampleLight
// =====================================================================================
#ifdef HLRAD_SOFTSKY
#define NUMVERTEXNORMALS_DEFAULT 162
double          r_avertexnormals_default[NUMVERTEXNORMALS_DEFAULT][3] = {
	#include "anorms.h"
};
int NUMVERTEXNORMALS = NUMVERTEXNORMALS_DEFAULT;
double (*r_avertexnormals)[3] = r_avertexnormals_default;
#else
#define NUMVERTEXNORMALS	162
double          r_avertexnormals[NUMVERTEXNORMALS][3] = {
//#include "../common/anorms.h"
	#include "anorms.h" //--vluzacn
};
#endif
#ifdef HLRAD_SOFTSKY
#define SKYLEVEL 7
typedef double point_t[3];
typedef int triangle_t[3];
void BuildDiffuseNormals (int n)
{
	assume (n==0, assume_first);
	int i, j, k, l, m;
	int numpoints = 6;
	point_t *points = (point_t *)malloc (((1 << (2 * SKYLEVEL)) + 2) * sizeof (point_t));
	points[0][0] = 1, points[0][1] = 0, points[0][2] = 0;
	points[1][0] = -1,points[1][1] = 0, points[1][2] = 0;
	points[2][0] = 0, points[2][1] = 1, points[2][2] = 0;
	points[3][0] = 0, points[3][1] = -1,points[3][2] = 0;
	points[4][0] = 0, points[4][1] = 0, points[4][2] = 1;
	points[5][0] = 0, points[5][1] = 0, points[5][2] = -1;
	int numtriangles = 8;
	triangle_t *triangles = (triangle_t *)malloc (((1 << (2 * SKYLEVEL)) * 2) * sizeof (triangle_t));
	triangles[0][0] = 0, triangles[0][1] = 2, triangles[0][2] = 4;
	triangles[1][0] = 0, triangles[1][1] = 2, triangles[1][2] = 5;
	triangles[2][0] = 0, triangles[2][1] = 3, triangles[2][2] = 4;
	triangles[3][0] = 0, triangles[3][1] = 3, triangles[3][2] = 5;
	triangles[4][0] = 1, triangles[4][1] = 2, triangles[4][2] = 4;
	triangles[5][0] = 1, triangles[5][1] = 2, triangles[5][2] = 5;
	triangles[6][0] = 1, triangles[6][1] = 3, triangles[6][2] = 4;
	triangles[7][0] = 1, triangles[7][1] = 3, triangles[7][2] = 5;
	for (i = 1; i < SKYLEVEL; i++)
	{
		for (j = 0; j < (1 << (2 * i)) * 2; j++)
		{
			int e[3];
			for (k = 0; k < 3; k++)
			{
				point_t mid;
				double len;
				VectorAdd (points[triangles[j][k]], points[triangles[j][(k+1)%3]], mid);
				len = sqrt (mid[0]*mid[0]+mid[1]*mid[1]+mid[2]*mid[2]);
				VectorScale (mid, 1/len, mid);
				for (l = 0; l < numpoints; l++)
				{
					for (m = 0; m < 3; m++)
					{
						if (fabs (points[l][m] - mid[m]) > 0.00001)
							break;
					}
					if (m == 3)
						break;
				}
				e[k] = l;
				if (l == numpoints)
				{
					hlassume (numpoints < (1 << (2 * SKYLEVEL)) + 2, assume_first);
					VectorCopy (mid, points[l]);
					numpoints++;
				}
			}
			for (k = 0; k < 3; k++)
			{
				hlassume (numtriangles < (1 << (2 * SKYLEVEL)) * 2, assume_first);
				triangles[numtriangles][0] = e[k], triangles[numtriangles][1] = triangles[j][(k+1)%3], triangles[numtriangles][2] = e[(k+1)%3];
				numtriangles++;
			}
			VectorCopy (e, triangles[j]);
		}
	}
	hlassume (numpoints == (1 << (2 * SKYLEVEL)) + 2, assume_first);
	hlassume (numtriangles == (1 << (2 * SKYLEVEL)) * 2, assume_first);
	NUMVERTEXNORMALS = (1 << (2 * SKYLEVEL)) + 2;
	r_avertexnormals = points;
	free (triangles);
}
#endif
static void     GatherSampleLight(const vec3_t pos, const byte* const pvs, const vec3_t normal, vec3_t* sample, byte* styles
#ifdef HLRAD_GatherPatchLight
								  , int step
#endif
								  )
{
    int             i;
    directlight_t*  l;
    vec3_t          add;
    vec3_t          delta;
    float           dot, dot2;
    float           dist;
    float           ratio;
#ifdef HLRAD_OPACITY // AJM
    float           l_opacity;
#endif
    int             style_index;
#ifdef HLRAD_GatherPatchLight
	int				step_match;
#endif
#ifdef HLRAD_MULTISKYLIGHT
	bool sky_used = false;
	bool sky_diffuse_used = false;
#ifdef HLRAD_OPAQUE_STYLE
	vec3_t sky_diffuse_total[ALLSTYLES+1];
#else
	vec3_t sky_diffuse_total;
#endif
#else
    directlight_t*  sky_used = NULL;
#endif
#ifdef HLRAD_STYLE_CORING
	vec3_t			adds[ALLSTYLES];
	int				style;
	memset (adds, 0, ALLSTYLES * sizeof(vec3_t));
#endif

#ifdef HLRAD_SKYFIX_FIX
#ifdef HLRAD_VIS_FIX
    for (i = 0; i < 1 + g_dmodels[0].visleafs; i++)
#else
    for (i = 0; i < g_numleafs; i++)
#endif
#else
#ifdef HLRAD_VIS_FIX
    for (i = 1; i < 1 + g_dmodels[0].visleafs; i++)
#else
    for (i = 1; i < g_numleafs; i++)
#endif
#endif
    {
        l = directlights[i];
#ifdef HLRAD_SKYFIX_FIX
		{
            for (; l; l = l->next)
            {
                if (((l->type == emit_skylight) && (g_sky_lighting_fix)) || i != 0 && (pvs[(i - 1) >> 3] & (1 << ((i - 1) & 7))))
                {
#else
        if (l)
        {
            if (((l->type == emit_skylight) && (g_sky_lighting_fix)) || (pvs[(i - 1) >> 3] & (1 << ((i - 1) & 7))))
            {
                for (; l; l = l->next)
                {
#endif
#ifdef HLRAD_OPAQUE_STYLE
					int opaquestyle = -1;
#endif
                    // skylights work fundamentally differently than normal lights
                    if (l->type == emit_skylight)
                    {
#ifdef HLRAD_MULTISKYLIGHT
						if (sky_used && !g_sky_lighting_fix)
							continue;
						sky_used = true;
						VectorClear (add);
						do
						{
#ifdef HLRAD_GatherPatchLight
							step_match = (int)l->topatch;
							if (step != step_match)
								continue;
#endif
#ifdef HLRAD_ALLOWZEROBRIGHTNESS
							if (VectorCompare (l->intensity, vec3_origin))
								continue;
#endif
							// make sure the angle is okay
							dot = -DotProduct(normal, l->normal);
							if (dot <= NORMAL_EPSILON) //ON_EPSILON / 10 //--vluzacn
							{
								continue;
							}

							// search back to see if we can hit a sky brush
#ifdef ZHLT_LARGERANGE
							VectorScale(l->normal, -BOGUS_RANGE, delta);
#else
							VectorScale(l->normal, -10000, delta);
#endif
							VectorAdd(pos, delta, delta);
	#ifdef HLRAD_OPAQUEINSKY_FIX
							vec3_t skyhit;
							VectorCopy (delta, skyhit);
	#endif
							if (TestLine(pos, delta
	#ifdef HLRAD_OPAQUEINSKY_FIX
								, skyhit
	#endif
								) != CONTENTS_SKY)
							{
								continue;                      // occluded
							}

	#ifdef HLRAD_HULLU
							vec3_t transparency = {1.0,1.0,1.0};
	#endif
							if (TestSegmentAgainstOpaqueList(pos, 
	#ifdef HLRAD_OPAQUEINSKY_FIX
								skyhit
	#else
								delta
	#endif
	#ifdef HLRAD_HULLU
								, transparency
	#endif
	#ifdef HLRAD_OPAQUE_STYLE
								, opaquestyle
	#endif
								))
							{
								continue;
							}

							VectorScale(l->intensity, dot, add);
	#ifdef HLRAD_HULLU
							VectorMultiply(add, transparency, add);
	#endif
						}
						while (0);
						do
						{
#ifdef HLRAD_GatherPatchLight
							step_match = 0;
	#ifdef HLRAD_SOFTSKY
							if (g_softsky)
								step_match = 1;
	#endif
							if (step != step_match)
								continue;
#endif
							vec3_t          sky_intensity;
					#ifdef HLRAD_WHOME
							VectorScale(l->diffuse_intensity, g_indirect_sun / (NUMVERTEXNORMALS * 2), sky_intensity);
					#else
							VectorScale(l->intensity, g_indirect_sun / (NUMVERTEXNORMALS * 2), sky_intensity);
					#endif
							if (g_indirect_sun <= 0.0 || VectorMaximum (sky_intensity) <= 0.0)
								continue;

							if (sky_diffuse_used == false)
							{
								sky_diffuse_used = true;
								vec3_t			add;
								int             j;
								vec3_t          delta;
						#ifdef HLRAD_OPAQUE_STYLE
								memset (sky_diffuse_total, 0, (ALLSTYLES+1) * sizeof (vec3_t));
						#else
								vec3_t          total;
								VectorClear (total);
						#endif
								for (j = 0; j < NUMVERTEXNORMALS; j++)
								{
						#ifdef HLRAD_OPAQUE_STYLE
									int				opaquestyle = -1;
						#endif
									// make sure the angle is okay
									dot = -DotProduct(normal, r_avertexnormals[j]);
									if (dot <= ON_EPSILON / 10)
									{
										continue;
									}

									// search back to see if we can hit a sky brush
#ifdef ZHLT_LARGERANGE
									VectorScale(r_avertexnormals[j], -BOGUS_RANGE, delta);
#else
									VectorScale(r_avertexnormals[j], -10000, delta);
#endif
									VectorAdd(pos, delta, delta);
						#ifdef HLRAD_OPAQUEINSKY_FIX
									vec3_t skyhit;
									VectorCopy (delta, skyhit);
						#endif
									if (TestLine(pos, delta
						#ifdef HLRAD_OPAQUEINSKY_FIX
										, skyhit
						#endif
										) != CONTENTS_SKY)
									{
										continue;                                  // occluded
									}

						#ifdef HLRAD_OPAQUE_DIFFUSE_FIX
						#ifdef HLRAD_HULLU
									vec3_t transparency = {1.0,1.0,1.0};
						#endif
									if (TestSegmentAgainstOpaqueList(pos, 
						#ifdef HLRAD_OPAQUEINSKY_FIX
										skyhit
						#else
										delta
						#endif
						#ifdef HLRAD_HULLU
										, transparency
						#endif
						#ifdef HLRAD_OPAQUE_STYLE
										, opaquestyle
						#endif
										))
									{
										continue;
									}
						#endif /*HLRAD_OPAQUE_DIFFUSE_FIX*/
									//VectorScale(sky_intensity, dot, add);
									VectorFill (add, dot);
						#ifdef HLRAD_OPAQUE_DIFFUSE_FIX
						#ifdef HLRAD_HULLU
												VectorMultiply(add, transparency, add);
						#endif
						#endif /*HLRAD_OPAQUE_DIFFUSE_FIX*/
						#ifdef HLRAD_OPAQUE_STYLE
									VectorAdd (sky_diffuse_total[opaquestyle+1], add, sky_diffuse_total[opaquestyle+1]);
						#else
									VectorAdd(total, add, total);
						#endif
								}
						#ifndef HLRAD_OPAQUE_STYLE
								VectorCopy (total, sky_diffuse_total);
						#endif
							} /*if (sky_diffuse_used == false)*/

	#ifdef HLRAD_OPAQUE_STYLE
							{
								int addstyle;
								int opaquestyle;
								vec3_t diffuseadd;
								for (opaquestyle = -1; opaquestyle < ALLSTYLES; opaquestyle++)
								{
									addstyle = l->style;
									if (opaquestyle != -1)
									{
										if (addstyle == 0 || addstyle == opaquestyle)
											addstyle = opaquestyle;
										else
											continue;
									}
									VectorMultiply (sky_intensity, sky_diffuse_total[opaquestyle+1], diffuseadd);
									VectorAdd (adds[addstyle], diffuseadd, adds[addstyle]);
								}
							}
	#else
							VectorMultiply (sky_intensity, sky_diffuse_total, sky_intensity);
							VectorAdd (add, sky_intensity, add);
	#endif
						}
						while (0);
						if (VectorMaximum (add) <= 0.0)
							continue;

#else /*HLRAD_MULTISKYLIGHT*/
                        // only allow one of each sky type to hit any given point
                        if (sky_used)
                        {
                            continue;
                        }
                        sky_used = l;

                        // make sure the angle is okay
                        dot = -DotProduct(normal, l->normal);
                        if (dot <= NORMAL_EPSILON) //ON_EPSILON / 10 //--vluzacn
                        {
                            continue;
                        }

                        // search back to see if we can hit a sky brush
#ifdef ZHLT_LARGERANGE
                        VectorScale(l->normal, -BOGUS_RANGE, delta);
#else
                        VectorScale(l->normal, -10000, delta);
#endif
                        VectorAdd(pos, delta, delta);
#ifdef HLRAD_OPAQUEINSKY_FIX
						vec3_t skyhit;
						VectorCopy (delta, skyhit);
#endif
                        if (TestLine(pos, delta
#ifdef HLRAD_OPAQUEINSKY_FIX
							, skyhit
#endif
							) != CONTENTS_SKY)
                        {
                            continue;                      // occluded
                        }

#ifdef HLRAD_HULLU
			            vec3_t transparency = {1.0,1.0,1.0};
#endif
                        if (TestSegmentAgainstOpaqueList(pos, 
#ifdef HLRAD_OPAQUEINSKY_FIX
							skyhit
#else
							delta
#endif
#ifdef HLRAD_HULLU
							, transparency
#endif
							))
                        {
                            continue;
                        }

                        VectorScale(l->intensity, dot, add);
#ifdef HLRAD_HULLU
                        VectorMultiply(add, transparency, add);
#endif

#endif /*HLRAD_MULTISKYLIGHT*/
                    }
                    else
                    {
#ifdef HLRAD_GatherPatchLight
						step_match = (int)l->topatch;
						if (step != step_match)
							continue;
#endif
#ifdef HLRAD_ALLOWZEROBRIGHTNESS
						if (VectorCompare (l->intensity, vec3_origin))
							continue;
#endif
                        float           denominator;

                        VectorSubtract(l->origin, pos, delta);
                        dist = VectorNormalize(delta);
                        dot = DotProduct(delta, normal);
                        //                        if (dot <= 0.0)
                        //                            continue;
                        if (dot <= NORMAL_EPSILON) //ON_EPSILON / 10 //--vluzacn
                        {
                            continue;                      // behind sample surface
                        }

                        if (dist < 1.0)
                        {
                            dist = 1.0;
                        }

                        // Variable power falloff (1 = inverse linear, 2 = inverse square
                        denominator = dist * l->fade;
#ifndef HLRAD_ARG_MISC
                        if (l->falloff == 2)
#endif
                        {
                            denominator *= dist;
                        }

                        switch (l->type)
                        {
                        case emit_point:
                        {
                            // Variable power falloff (1 = inverse linear, 2 = inverse square
                            vec_t           denominator = dist * l->fade;

#ifndef HLRAD_ARG_MISC
                            if (l->falloff == 2)
#endif
                            {
                                denominator *= dist;
                            }
                            ratio = dot / denominator;
                            VectorScale(l->intensity, ratio, add);
                            break;
                        }

                        case emit_surface:
                        {
                            dot2 = -DotProduct(delta, l->normal);
#ifdef HLRAD_CUSTOMTEXLIGHT
							if (dot2 <= l->stopdot2 + NORMAL_EPSILON)
							{
								continue;
							}
							ratio = dot * dot2 / (dist * dist);
							if (dot2 <= l->stopdot)
							{
								ratio *= (dot2 - l->stopdot2) / (l->stopdot - l->stopdot2);
							}
#else
                            if (dot2 <= NORMAL_EPSILON) //ON_EPSILON / 10 //--vluzacn
                            {
                                continue;                  // behind light surface
                            }

                            // Variable power falloff (1 = inverse linear, 2 = inverse square
                            vec_t           denominator = dist * g_fade;
#ifndef HLRAD_ARG_MISC
                            if (g_falloff == 2)
#endif
                            {
                                denominator *= dist;
                            }
                            ratio = dot * dot2 / denominator;
#endif

                            VectorScale(l->intensity, ratio, add);
                            break;
                        }

                        case emit_spotlight:
                        {
                            dot2 = -DotProduct(delta, l->normal);
                            if (dot2 <= l->stopdot2)
                            {
                                continue;                  // outside light cone
                            }

                            // Variable power falloff (1 = inverse linear, 2 = inverse square
                            vec_t           denominator = dist * l->fade;
#ifndef HLRAD_ARG_MISC
                            if (l->falloff == 2)
#endif
                            {
                                denominator *= dist;
                            }
                            ratio = dot * dot2 / denominator;

                            if (dot2 <= l->stopdot)
                            {
                                ratio *= (dot2 - l->stopdot2) / (l->stopdot - l->stopdot2);
                            }
                            VectorScale(l->intensity, ratio, add);
                            break;
                        }

                        default:
                        {
                            hlassume(false, assume_BadLightType);
                            break;
                        }
                        }
                    }

#ifndef HLRAD_STYLE_CORING
                    if (VectorMaximum(add) > (l->style ? g_coring : 0))
#endif
                    {
#ifdef HLRAD_HULLU
                 	    vec3_t transparency = {1.0,1.0,1.0};
#endif 

                        if (l->type != emit_skylight && TestLine(pos, l->origin) != CONTENTS_EMPTY)
                        {
                            continue;                      // occluded
                        }

                        if (l->type != emit_skylight)
                        {                                  // Don't test from light_environment entities to face, the special sky code occludes correctly
                            if (TestSegmentAgainstOpaqueList(pos, l->origin
#ifdef HLRAD_HULLU
								, transparency
#endif
#ifdef HLRAD_OPAQUE_STYLE
								, opaquestyle
#endif
								))
                            {
                                continue;
                            }
                        }

#ifdef HLRAD_OPACITY
                        //VectorScale(add, l_opacity, add);
#endif

#ifdef HLRAD_STYLE_CORING
#ifdef HLRAD_HULLU
                        VectorMultiply(add,transparency,add);
#endif
#ifdef HLRAD_OPAQUE_STYLE
						{
							int addstyle = l->style;
							if (opaquestyle != -1)
							{
								if (addstyle == 0 || addstyle == opaquestyle)
									addstyle = opaquestyle;
								else
									continue;
							}
							VectorAdd (adds[addstyle], add, adds[addstyle]);
						}
#else
						VectorAdd (adds[l->style], add, adds[l->style]);
#endif
#else
                        for (style_index = 0; style_index < MAXLIGHTMAPS; style_index++)
                        {
                            if (styles[style_index] == l->style || styles[style_index] == 255)
                            {
                                break;
                            }
                        }

                        if (style_index == MAXLIGHTMAPS)
                        {
#ifdef HLRAD_READABLE_EXCEEDSTYLEWARNING
							if (++stylewarningcount >= stylewarningnext)
							{
								stylewarningnext = stylewarningcount * 2;
								Warning("Too many direct light styles on a face(%f,%f,%f)", pos[0], pos[1], pos[2]);
								Warning(" total %d warnings for too many styles", stylewarningcount);
							}
#else
                            Warning("Too many direct light styles on a face(%f,%f,%f)", pos[0], pos[1], pos[2]);
#endif
                            continue;
                        }

                        if (styles[style_index] == 255)
                        {
                            styles[style_index] = l->style;
                        }

#ifdef HLRAD_HULLU
                        VectorMultiply(add,transparency,add);
#endif
                        VectorAdd(sample[style_index], add, sample[style_index]);
#endif
                    }
                }
            }
        }
    }

#ifndef HLRAD_MULTISKYLIGHT
    if (sky_used && g_indirect_sun != 0.0)
    {
        vec3_t          total;
        int             j;
		vec3_t          sky_intensity;

		// -----------------------------------------------------------------------------------
		// Changes by Adam Foster - afoster@compsoc.man.ac.uk
		// Instead of using intensity from sky_used->intensity, get it from the new sky_used->diffuse_intensity
#ifdef HLRAD_WHOME
		VectorScale(sky_used->diffuse_intensity, g_indirect_sun / (NUMVERTEXNORMALS * 2), sky_intensity);
#else
        VectorScale(sky_used->intensity, g_indirect_sun / (NUMVERTEXNORMALS * 2), sky_intensity);
#endif
		// That should be it. Who knows - it might actually work!
        // AJM: It DOES actually work. Havent you ever heard of beta testing....
		// -----------------------------------------------------------------------------------

        total[0] = total[1] = total[2] = 0.0;
        for (j = 0; j < NUMVERTEXNORMALS; j++)
        {
            // make sure the angle is okay
            dot = -DotProduct(normal, r_avertexnormals[j]);
            if (dot <= NORMAL_EPSILON) //ON_EPSILON / 10 //--vluzacn
            {
                continue;
            }

            // search back to see if we can hit a sky brush
#ifdef ZHLT_LARGERANGE
            VectorScale(r_avertexnormals[j], -BOGUS_RANGE, delta);
#else
            VectorScale(r_avertexnormals[j], -10000, delta);
#endif
            VectorAdd(pos, delta, delta);
#ifdef HLRAD_OPAQUEINSKY_FIX
			vec3_t skyhit;
			VectorCopy (delta, skyhit);
#endif
            if (TestLine(pos, delta
#ifdef HLRAD_OPAQUEINSKY_FIX
				, skyhit
#endif
				) != CONTENTS_SKY)
            {
                continue;                                  // occluded
            }

#ifdef HLRAD_OPAQUE_DIFFUSE_FIX
#ifdef HLRAD_HULLU
			vec3_t transparency = {1.0,1.0,1.0};
#endif
			if (TestSegmentAgainstOpaqueList(pos, 
#ifdef HLRAD_OPAQUEINSKY_FIX
				skyhit
#else
				delta
#endif
#ifdef HLRAD_HULLU
				, transparency
#endif
				))
			{
				continue;
			}
#endif /*HLRAD_OPAQUE_DIFFUSE_FIX*/
            VectorScale(sky_intensity, dot, add);
#ifdef HLRAD_OPAQUE_DIFFUSE_FIX
#ifdef HLRAD_HULLU
                        VectorMultiply(add, transparency, add);
#endif
#endif /*HLRAD_OPAQUE_DIFFUSE_FIX*/
            VectorAdd(total, add, total);
        }
        if (VectorMaximum(total) > 0)
        {
#ifdef HLRAD_STYLE_CORING
			VectorAdd (adds[sky_used->style], total, adds[sky_used->style]);
#else
            for (style_index = 0; style_index < MAXLIGHTMAPS; style_index++)
            {
                if (styles[style_index] == sky_used->style || styles[style_index] == 255)
                {
                    break;
                }
            }

            if (style_index == MAXLIGHTMAPS)
            {
#ifdef HLRAD_READABLE_EXCEEDSTYLEWARNING
				if (++stylewarningcount >= stylewarningnext)
				{
					stylewarningnext = stylewarningcount * 2;
					Warning("Too many direct light styles on a face(%f,%f,%f)\n", pos[0], pos[1], pos[2]);
					Warning(" total %d warnings for too many styles", stylewarningcount);
				}
#else
                Warning("Too many direct light styles on a face(%f,%f,%f)\n", pos[0], pos[1], pos[2]);
#endif
                return;
            }

            if (styles[style_index] == 255)
            {
                styles[style_index] = sky_used->style;
            }

            VectorAdd(sample[style_index], total, sample[style_index]);
#endif
        }
    }
#endif /*HLRAD_MULTISKYLIGHT*/
#ifdef HLRAD_STYLE_CORING
	for (style = 0; style < ALLSTYLES; ++style)
	{
#ifdef HLRAD_AUTOCORING
		if (VectorMaximum(adds[style]) > g_corings[style] * 0.5)
#else
		if (VectorMaximum(adds[style]) > g_corings[style])
#endif
		{
	#ifdef HLRAD_AUTOCORING
			for (style_index = 0; style_index < ALLSTYLES; style_index++)
	#else
			for (style_index = 0; style_index < MAXLIGHTMAPS; style_index++)
	#endif
			{
				if (styles[style_index] == style || styles[style_index] == 255)
				{
					break;
				}
			}

	#ifdef HLRAD_AUTOCORING
			if (style_index == ALLSTYLES)
	#else
			if (style_index == MAXLIGHTMAPS)
	#endif
			{
#ifdef HLRAD_READABLE_EXCEEDSTYLEWARNING
				if (++stylewarningcount >= stylewarningnext)
				{
					stylewarningnext = stylewarningcount * 2;
					Warning("Too many direct light styles on a face(%f,%f,%f)", pos[0], pos[1], pos[2]);
					Warning(" total %d warnings for too many styles", stylewarningcount);
				}
#else
				Warning("Too many direct light styles on a face(%f,%f,%f)", pos[0], pos[1], pos[2]);
#endif
				return;
			}

			if (styles[style_index] == 255)
			{
				styles[style_index] = style;
			}

			VectorAdd(sample[style_index], adds[style], sample[style_index]);
		}
	}
#endif
}

// =====================================================================================
//  AddSampleToPatch
//      Take the sample's collected light and add it back into the apropriate patch for the radiosity pass.
// =====================================================================================
#ifdef ZHLT_TEXLIGHT
static void     AddSampleToPatch(const sample_t* const s, const int facenum, int style
	#ifdef HLRAD_AddSampleToPatch_PRECISE
								 , const vec3_t pos
	#endif
								 ) //LRC
#else
static void     AddSampleToPatch(const sample_t* const s, const int facenum
	#ifdef HLRAD_AddSampleToPatch_PRECISE
								 , const vec3_t pos
	#endif
								 )
#endif
{
    patch_t*        patch;
    BoundingBox     bounds;
    int             i;

#ifndef HLRAD_GatherPatchLight
    if (g_numbounce == 0)
    {
        return;
    }
#endif

    for (patch = g_face_patches[facenum]; patch; patch = patch->next)
    {
#ifdef HLRAD_AddSampleToPatch_PRECISE
		if (!point_in_winding (*(patch->winding), *getPlaneFromFaceNumber (patch->faceNumber), pos))
			goto nextpatch;
#else
        // see if the point is in this patch (roughly)
        patch->winding->getBounds(bounds);
        for (i = 0; i < 3; i++)
        {
            if (bounds.m_Mins[i] > s->pos[i] + 16)
            {
                goto nextpatch;
            }
            if (bounds.m_Maxs[i] < s->pos[i] - 16)
            {
                goto nextpatch;
            }
        }
#endif
#ifdef HLRAD_AUTOCORING
		if (style == 0)
		{
			patch->samples++;
		}
#endif

        // add the sample to the patch
#ifdef ZHLT_TEXLIGHT
        //LRC:
	#ifdef HLRAD_AUTOCORING
		for (i = 0; i < ALLSTYLES && patch->totalstyle_all[i] != 255; i++)
		{
			if (patch->totalstyle_all[i] == style)
				break;
		}
		if (i == ALLSTYLES)
	#else
		for (i = 0; i < MAXLIGHTMAPS && patch->totalstyle[i] != 255; i++)
		{
			if (patch->totalstyle[i] == style)
				break;
		}
		if (i == MAXLIGHTMAPS)
	#endif
		{
#ifdef HLRAD_READABLE_EXCEEDSTYLEWARNING
			if (++stylewarningcount >= stylewarningnext)
			{
				stylewarningnext = stylewarningcount * 2;
				Warning("Too many direct light styles on a face(?,?,?)\n");
				Warning(" total %d warnings for too many styles", stylewarningcount);
			}
#else
			Warning("Too many direct light styles on a face(?,?,?)\n");
#endif
		}
		else
		{
	#ifdef HLRAD_AUTOCORING
			if (patch->totalstyle_all[i] == 255)
			{
				patch->totalstyle_all[i] = style;
			}
			VectorAdd(patch->samplelight_all[i], s->light, patch->samplelight_all[i]);
	#else
			if (patch->totalstyle[i] == 255)
			{
				patch->totalstyle[i] = style;
			}

	        patch->samples[i]++;
			VectorAdd(patch->samplelight[i], s->light, patch->samplelight[i]);
	#endif
		}
        //LRC (ends)
#else
        patch->samples++;
        VectorAdd(patch->samplelight, s->light, patch->samplelight);
#endif
        //return;

      nextpatch:;
    }

    // don't worry if some samples don't find a patch
}

// =====================================================================================
//  GetPhongNormal
// =====================================================================================
void            GetPhongNormal(int facenum, vec3_t spot, vec3_t phongnormal)
{
    int             j;
#ifdef HLRAD_GetPhongNormal_VL
	int				s; // split every edge into two parts
#endif
    const dface_t*  f = g_dfaces + facenum;
    const dplane_t* p = getPlaneFromFace(f);
    vec3_t          facenormal;

    VectorCopy(p->normal, facenormal);
    VectorCopy(facenormal, phongnormal);

#ifndef HLRAD_CUSTOMSMOOTH
    if (g_smoothing_threshold > 0.0)
#endif
    {
        // Calculate modified point normal for surface
        // Use the edge normals iff they are defined.  Bend the surface towards the edge normal(s)
        // Crude first attempt: find nearest edge normal and do a simple interpolation with facenormal.
        // Second attempt: find edge points+center that bound the point and do a three-point triangulation(baricentric)
        // Better third attempt: generate the point normals for all vertices and do baricentric triangulation.

        for (j = 0; j < f->numedges; j++)
        {
            vec3_t          p1;
            vec3_t          p2;
            vec3_t          v1;
            vec3_t          v2;
            vec3_t          vspot;
            unsigned        prev_edge;
            unsigned        next_edge;
            int             e;
            int             e1;
            int             e2;
            edgeshare_t*    es;
            edgeshare_t*    es1;
            edgeshare_t*    es2;
            float           a1;
            float           a2;
            float           aa;
            float           bb;
            float           ab;

            if (j)
            {
#ifdef HLRAD_NEGATIVEDIVIDEND_MISCFIX
                prev_edge = f->firstedge + ((j + f->numedges - 1) % f->numedges);
#else
                prev_edge = f->firstedge + ((j - 1) % f->numedges);
#endif
            }
            else
            {
                prev_edge = f->firstedge + f->numedges - 1;
            }

            if ((j + 1) != f->numedges)
            {
                next_edge = f->firstedge + ((j + 1) % f->numedges);
            }
            else
            {
                next_edge = f->firstedge;
            }

            e = g_dsurfedges[f->firstedge + j];
            e1 = g_dsurfedges[prev_edge];
            e2 = g_dsurfedges[next_edge];

            es = &g_edgeshare[abs(e)];
            es1 = &g_edgeshare[abs(e1)];
            es2 = &g_edgeshare[abs(e2)];

#ifdef HLRAD_GetPhongNormal_VL
			if ((es->coplanar || !es->smooth) && (es1->coplanar || !es1->smooth) && (es2->coplanar || !es2->smooth))
#else
            if (
                (es->coplanar && es1->coplanar && es2->coplanar)
                ||
                (VectorCompare(es->interface_normal, vec3_origin) &&
                 VectorCompare(es1->interface_normal, vec3_origin) &&
                 VectorCompare(es2->interface_normal, vec3_origin)))
#endif
            {
                continue;
            }

            if (e > 0)
            {
                VectorCopy(g_dvertexes[g_dedges[e].v[0]].point, p1);
                VectorCopy(g_dvertexes[g_dedges[e].v[1]].point, p2);
            }
            else
            {
                VectorCopy(g_dvertexes[g_dedges[-e].v[1]].point, p1);
                VectorCopy(g_dvertexes[g_dedges[-e].v[0]].point, p2);
            }

            // Adjust for origin-based models
            VectorAdd(p1, g_face_offset[facenum], p1);
            VectorAdd(p2, g_face_offset[facenum], p2);
#ifdef HLRAD_GetPhongNormal_VL
		for (s = 0; s < 2; s++)
		{
			vec3_t s1, s2;
			if (s == 0)
			{
				VectorCopy(p1, s1);
			}
			else
			{
				VectorCopy(p2, s1);
			}

			VectorAdd(p1,p2,s2); // edge center
			VectorScale(s2,0.5,s2);

            VectorSubtract(s1, g_face_centroids[facenum], v1);
            VectorSubtract(s2, g_face_centroids[facenum], v2);
#else

            // Build vectors from the middle of the face to the edge vertexes and the sample pos.
            VectorSubtract(p1, g_face_centroids[facenum], v1);
            VectorSubtract(p2, g_face_centroids[facenum], v2);
#endif
            VectorSubtract(spot, g_face_centroids[facenum], vspot);

            aa = DotProduct(v1, v1);
            bb = DotProduct(v2, v2);
            ab = DotProduct(v1, v2);
            a1 = (bb * DotProduct(v1, vspot) - ab * DotProduct(vspot, v2)) / (aa * bb - ab * ab);
            a2 = (DotProduct(vspot, v2) - a1 * ab) / bb;

            // Test center to sample vector for inclusion between center to vertex vectors (Use dot product of vectors)
#ifdef HLRAD_GetPhongNormal_VL
            if (a1 >= -ON_EPSILON && a2 >= -ON_EPSILON)
#else
            if (a1 >= 0.0 && a2 >= 0.0)
#endif
            {
                // calculate distance from edge to pos
                vec3_t          n1, n2;
                vec3_t          temp;

#ifdef HLRAD_GetPhongNormal_VL
				/*
				if (s == 0)
				{VectorCopy(es1->interface_normal, n1);}
				else
				{VectorCopy(es2->interface_normal, n1);}
				if (VectorCompare(n1, vec3_origin))
				{VectorCopy(facenormal, n1);}
				if (VectorCompare(es->interface_normal, vec3_origin))
				{VectorAdd(n1, facenormal, n1);}
				else
				{VectorAdd(n1, es->interface_normal, n1);}
				VectorSubtract(n1, facenormal, n1);
				VectorNormalize(n1);
				*/
				if (es->smooth)
					if (s == 0)
					{VectorCopy(es->vertex_normal[e>0?0:1], n1);}
					else
					{VectorCopy(es->vertex_normal[e>0?1:0], n1);}
				else if (s == 0 && es1->smooth)
				{VectorCopy(es1->vertex_normal[e1>0?1:0], n1);}
				else if (s == 1 && es2->smooth)
				{VectorCopy(es2->vertex_normal[e2>0?0:1], n1);}
				else
				{VectorCopy(facenormal, n1);}

				if (es->smooth)
				{VectorCopy(es->interface_normal, n2);}
				else
				{VectorCopy(facenormal, n2);}
#else
                VectorAdd(es->interface_normal, es1->interface_normal, n1)

                if (VectorCompare(n1, vec3_origin))
                {
                    VectorCopy(facenormal, n1);
                }
                VectorNormalize(n1);

                VectorAdd(es->interface_normal, es2->interface_normal, n2);

                if (VectorCompare(n2, vec3_origin))
                {
                    VectorCopy(facenormal, n2);
                }
                VectorNormalize(n2);
#endif

                // Interpolate between the center and edge normals based on sample position
                VectorScale(facenormal, 1.0 - a1 - a2, phongnormal);
                VectorScale(n1, a1, temp);
                VectorAdd(phongnormal, temp, phongnormal);
                VectorScale(n2, a2, temp);
                VectorAdd(phongnormal, temp, phongnormal);
                VectorNormalize(phongnormal);
                break;
            }
#ifdef HLRAD_GetPhongNormal_VL
		} // s=0,1
#endif
        }
    }
}

const vec3_t    s_circuscolors[] = {
    {100000.0,  100000.0,   100000.0},                              // white
    {100000.0,  0.0,        0.0     },                              // red
    {0.0,       100000.0,   0.0     },                              // green
    {0.0,       0.0,        100000.0},                              // blue
    {0.0,       100000.0,   100000.0},                              // cyan
    {100000.0,  0.0,        100000.0},                              // magenta
    {100000.0,  100000.0,   0.0     }                               // yellow
};

// =====================================================================================
//  BuildFacelights
// =====================================================================================
void            BuildFacelights(const int facenum)
{
    dface_t*        f;
#ifdef HLRAD_AUTOCORING
	unsigned char	f_styles[ALLSTYLES];
	sample_t		*fl_samples[ALLSTYLES];
	vec3_t			sampled[ALLSTYLES];
#else
    vec3_t          sampled[MAXLIGHTMAPS];
#endif
    lightinfo_t     l;
    int             i;
    int             j;
    int             k;
    sample_t*       s;
    vec_t*          spot;
    patch_t*        patch;
    const dplane_t* plane;
    byte            pvs[(MAX_MAP_LEAFS + 7) / 8];
    int             thisoffset = -1, lastoffset = -1;
    int             lightmapwidth;
    int             lightmapheight;
    int             size;
#ifdef HLRAD_TRANSLUCENT
	vec3_t			spot2, normal2;
	vec3_t			delta;
	byte			pvs2[(MAX_MAP_LEAFS + 7) / 8];
	int				thisoffset2 = -1, lastoffset2 = -1;
#endif

#ifndef HLRAD_TRANCPARENCYLOSS_FIX
#ifdef HLRAD_HULLU
    bool            b_transparency_loss = false;
    vec_t           light_left_for_facelight = 1.0;
#endif
#endif

    f = &g_dfaces[facenum];

    //
    // some surfaces don't need lightmaps
    //
    f->lightofs = -1;
#ifdef HLRAD_AUTOCORING
    for (j = 0; j < ALLSTYLES; j++)
    {
        f_styles[j] = 255;
    }
#else
    for (j = 0; j < MAXLIGHTMAPS; j++)
    {
        f->styles[j] = 255;
    }
#endif

    if (g_texinfo[f->texinfo].flags & TEX_SPECIAL)
    {
#ifdef HLRAD_AUTOCORING
		for (j = 0; j < MAXLIGHTMAPS; j++)
		{
			f->styles[j] = 255;
		}
#endif
        return;                                            // non-lit texture
    }

#ifdef HLRAD_AUTOCORING
	f_styles[0] = 0;
#else
    f->styles[0] = 0;                                      // Everyone gets the style zero map.
#endif
#ifdef HLRAD_STYLE_CORING
#ifdef ZHLT_TEXLIGHT
	if (g_face_patches[facenum] && g_face_patches[facenum]->emitstyle)
	{
#ifdef HLRAD_AUTOCORING
		f_styles[1] = g_face_patches[facenum]->emitstyle;
#else
		f->styles[1] = g_face_patches[facenum]->emitstyle;
#endif
	}
#endif
#endif

    memset(&l, 0, sizeof(l));

    l.surfnum = facenum;
    l.face = f;

#ifdef HLRAD_TRANSLUCENT
	VectorCopy (g_translucenttextures[g_texinfo[f->texinfo].miptex], l.translucent_v);
	l.translucent_b = !VectorCompare (l.translucent_v, vec3_origin);
#endif
#ifndef HLRAD_TRANCPARENCYLOSS_FIX
    //
    // get transparency loss (part of light go through transparency faces.. reduce facelight on these)
    //
#ifndef HLRAD_OPAQUE_NODE
#ifdef HLRAD_HULLU
    for(unsigned int m = 0; m < g_opaque_face_count; m++)
    {
        opaqueList_t* opaque = &g_opaque_face_list[m];
        if(opaque->facenum == facenum && opaque->transparency)
        {
            vec_t transparency = VectorAvg (opaque->transparency_scale); //vec_t transparency = opaque->transparency; //--vluzacn
            
            b_transparency_loss = true;
            
            light_left_for_facelight = 1.0 - transparency;
            if( light_left_for_facelight < 0.0 ) light_left_for_facelight = 0.0;
            if( light_left_for_facelight > 1.0 ) light_left_for_facelight = 1.0;
            
            break;
        }
    }
#endif
#endif
#endif

    //
    // rotate plane
    //
    plane = getPlaneFromFace(f);
    VectorCopy(plane->normal, l.facenormal);
    l.facedist = plane->dist;

    CalcFaceVectors(&l);
    CalcFaceExtents(&l);
    CalcPoints(&l);
#ifdef HLRAD_MDL_LIGHT_HACK
#ifndef HLRAD_MDL_LIGHT_HACK_NEW
	VectorCopy (g_face_offset[facenum], facesampleinfo[facenum].offset);
	for (i=0; i<2; ++i)
	{
		facesampleinfo[facenum].texmins[i] = l.texmins[i];
		facesampleinfo[facenum].texsize[i] = l.texsize[i];
		VectorCopy (l.textoworld[i], facesampleinfo[facenum].textoworld[i]);
		VectorCopy (l.worldtotex[i], facesampleinfo[facenum].worldtotex[i]);
	}
	VectorCopy (l.texorg, facesampleinfo[facenum].texorg);
#endif
#endif

    lightmapwidth = l.texsize[0] + 1;
    lightmapheight = l.texsize[1] + 1;

    size = lightmapwidth * lightmapheight;
    hlassume(size <= MAX_SINGLEMAP, assume_MAX_SINGLEMAP);

    facelight[facenum].numsamples = l.numsurfpt;

#ifdef HLRAD_AUTOCORING
	for (k = 0; k < ALLSTYLES; k++)
	{
		fl_samples[k] = (sample_t *)calloc (l.numsurfpt, sizeof(sample_t));
		hlassume (fl_samples[k] != NULL, assume_NoMemory);
	}
#else
    for (k = 0; k < MAXLIGHTMAPS; k++)
    {
        facelight[facenum].samples[k] = (sample_t*)calloc(l.numsurfpt, sizeof(sample_t));
#ifdef HLRAD_HLASSUMENOMEMORY
		hlassume (facelight[facenum].samples[k] != NULL, assume_NoMemory);
#endif
    }
#endif
#ifdef HLRAD_AUTOCORING
	for (patch = g_face_patches[facenum]; patch; patch = patch->next)
	{
		patch->totalstyle_all = (unsigned char *)malloc (ALLSTYLES * sizeof (unsigned char));
		hlassume (patch->totalstyle_all != NULL, assume_NoMemory);
		patch->samplelight_all = (vec3_t *)malloc (ALLSTYLES * sizeof (vec3_t));
		hlassume (patch->samplelight_all != NULL, assume_NoMemory);
		patch->totallight_all = (vec3_t *)malloc (ALLSTYLES * sizeof (vec3_t));
		hlassume (patch->totallight_all != NULL, assume_NoMemory);
		patch->directlight_all = (vec3_t *)malloc (ALLSTYLES * sizeof (vec3_t));
		hlassume (patch->directlight_all != NULL, assume_NoMemory);
		for (j = 0; j < ALLSTYLES; j++)
		{
			patch->totalstyle_all[j] = 255;
			VectorClear (patch->samplelight_all[j]);
			VectorClear (patch->totallight_all[j]);
			VectorClear (patch->directlight_all[j]);
		}
		patch->totalstyle_all[0] = 0;
	}
#endif

    spot = l.surfpt[0];
    for (i = 0; i < l.numsurfpt; i++, spot += 3)
    {
        vec3_t          pointnormal = { 0, 0, 0 };

#ifdef HLRAD_LERP_TEXNORMAL
		vec3_t spot_original;
		{
			vec_t s_vec = l.texmins[0] * 16 + (i % lightmapwidth) * TEXTURE_STEP;
			vec_t t_vec = l.texmins[1] * 16 + (i / lightmapwidth) * TEXTURE_STEP;
			SetSurfFromST (&l, spot_original, s_vec, t_vec);
			{
				// adjust sample's offset to 0
				vec_t scale;
				scale = DotProduct (l.texnormal, l.facenormal);
				VectorMA (spot_original, - DEFAULT_HUNT_OFFSET / scale, l.texnormal, spot_original);
			}
		}
#endif
#ifdef HLRAD_AUTOCORING
        for (k = 0; k < ALLSTYLES; k++)
        {
#ifdef HLRAD_LERP_TEXNORMAL
            VectorCopy(spot_original, fl_samples[k][i].pos);
#else
            VectorCopy(spot, fl_samples[k][i].pos);
#endif
        }
#else
        for (k = 0; k < MAXLIGHTMAPS; k++)
        {
#ifdef HLRAD_LERP_TEXNORMAL
            VectorCopy(spot_original, facelight[facenum].samples[k][i].pos);
#else
            VectorCopy(spot, facelight[facenum].samples[k][i].pos);
#endif
        }
#endif

        // get the PVS for the pos to limit the number of checks
        if (!g_visdatasize)
        {
            memset(pvs, 255, (g_numleafs + 7) / 8);
            lastoffset = -1;
        }
        else
        {
            dleaf_t*        leaf = PointInLeaf(spot);

            thisoffset = leaf->visofs;
            if (i == 0 || thisoffset != lastoffset)
            {
	#ifdef HLRAD_VIS_FIX
				if (thisoffset == -1)
				{
					memset (pvs, 0, (g_numleafs + 7) / 8);
				}
				else
				{
					DecompressVis(&g_dvisdata[leaf->visofs], pvs, sizeof(pvs));
				}
	#else
                hlassert(thisoffset != -1);
                DecompressVis(&g_dvisdata[leaf->visofs], pvs, sizeof(pvs));
	#endif
            }
            lastoffset = thisoffset;
        }
#ifdef HLRAD_TRANSLUCENT
		if (l.translucent_b)
		{
			VectorSubtract (g_face_centroids[facenum], spot, delta);
			VectorNormalize (delta);
			VectorMA (spot, 0.2, delta, spot2);
			VectorMA (spot2, -(g_translucentdepth + 2*DEFAULT_HUNT_OFFSET), l.facenormal, spot2);
			if (!g_visdatasize)
			{
				memset(pvs2, 255, (g_numleafs + 7) / 8);
				lastoffset2 = -1;
			}
			else
			{
				dleaf_t*        leaf2 = PointInLeaf(spot2);

				thisoffset2 = leaf2->visofs;
				if (i == 0 || thisoffset2 != lastoffset2)
				{
	#ifdef HLRAD_VIS_FIX
					if (thisoffset2 == -1)
					{
						memset(pvs2, 0, (g_numleafs + 7) / 8);
					}
					else
					{
						DecompressVis(&g_dvisdata[leaf2->visofs], pvs2, sizeof(pvs2));
					}
	#else
					hlassert(thisoffset2 != -1);
					DecompressVis(&g_dvisdata[leaf2->visofs], pvs2, sizeof(pvs2));
	#endif
				}
				lastoffset2 = thisoffset2;
			}
		}
#endif

        memset(sampled, 0, sizeof(sampled));

        // If we are doing "extra" samples, oversample the direct light around the point.
        if (g_extra)
        {
#ifdef HLRAD_WEIGHT_FIX
            int             weighting[3][3] = { {1, 1, 1}, {1, 1, 1}, {1, 1, 1} }; // because we are using 1/3 dist not 1/2
#else
            int             weighting[3][3] = { {5, 9, 5}, {9, 16, 9}, {5, 9, 5} };
#endif
            vec3_t          pos;
            int             s, t, subsamples = 0;

            for (t = -1; t <= 1; t++)
            {
                for (s = -1; s <= 1; s++)
                {
#ifndef HLRAD_CalcPoints_NEW
                    int             subsample = i + t * lightmapwidth + s;
                    int             sample_s = i % lightmapwidth;
                    int sample_t = i / lightmapwidth;

                    if ((0 <= s + sample_s) && (s + sample_s < lightmapwidth)
                        && (0 <= t + sample_t)&&(t + sample_t <lightmapheight))
#endif
                    {
#ifdef HLRAD_AUTOCORING
                        vec3_t          subsampled[ALLSTYLES];

                        for (j = 0; j < ALLSTYLES; j++)
#else
                        vec3_t          subsampled[MAXLIGHTMAPS];

                        for (j = 0; j < MAXLIGHTMAPS; j++)
#endif
                        {
                            VectorClear(subsampled[j]);
                        }

#ifdef HLRAD_CalcPoints_NEW
						vec_t s_vec = l.texmins[0] * 16 + (i % lightmapwidth + s * 1.0/3.0) * TEXTURE_STEP;
						vec_t t_vec = l.texmins[1] * 16 + (i / lightmapwidth + t * 1.0/3.0) * TEXTURE_STEP;
						if (SetSampleFromST (pos, &l, s_vec, t_vec, g_face_lightmode[facenum]) == LightOutside)
						{
							VectorCopy(l.surfpt[i], pos);
						}
#else
                        // Calculate the point one third of the way toward the "subsample point"
                        VectorCopy(l.surfpt[i], pos);
                        VectorAdd(pos, l.surfpt[i], pos);
                        VectorAdd(pos, l.surfpt[subsample], pos);
                        VectorScale(pos, 1.0 / 3.0, pos);
#endif

#ifdef HLRAD_PHONG_FROMORIGINAL
						// this will generate smoother light for cylinders partially embedded in solid,
						vec3_t pos_original;
						SetSurfFromST (&l, pos_original, s_vec, t_vec);
						{
							// adjust sample's offset to 0
							vec_t scale;
							scale = DotProduct (l.texnormal, l.facenormal);
							VectorMA (pos_original, - DEFAULT_HUNT_OFFSET / scale, l.texnormal, pos_original);
						}
                        GetPhongNormal(facenum, pos_original, pointnormal);
#else
                        GetPhongNormal(facenum, pos, pointnormal);
#endif
                        GatherSampleLight(pos, pvs, pointnormal, subsampled, 
#ifdef HLRAD_AUTOCORING
							f_styles
#else
							f->styles
#endif
#ifdef HLRAD_GatherPatchLight
							, 0
#endif
							);
#ifdef HLRAD_TRANSLUCENT
						if (l.translucent_b)
						{
#ifdef HLRAD_AUTOCORING
							vec3_t subsampled2[ALLSTYLES];
							for (j = 0; j < ALLSTYLES; j++)
#else
							vec3_t subsampled2[MAXLIGHTMAPS];
							for (j = 0; j < MAXLIGHTMAPS; j++)
#endif
							{
								VectorFill(subsampled2[j], 0);
							}
							VectorSubtract (g_face_centroids[facenum], pos, delta);
							VectorNormalize (delta);
							VectorMA (pos, 0.2, delta, spot2);
							VectorMA (spot2, -(g_translucentdepth + 2*DEFAULT_HUNT_OFFSET), l.facenormal, spot2);
							VectorSubtract (vec3_origin, pointnormal, normal2);
							GatherSampleLight(spot2, pvs2, normal2, subsampled2, 
#ifdef HLRAD_AUTOCORING
								f_styles
#else
								f->styles
#endif
	#ifdef HLRAD_GatherPatchLight
								, 0
	#endif
								);
#ifdef HLRAD_AUTOCORING
							for (j = 0; j < ALLSTYLES && f_styles[j] != 255; j++)
#else
							for (j = 0; j < MAXLIGHTMAPS && (f->styles[j] != 255); j++)
#endif
							{
								for (int x = 0; x < 3; x++)
								{
									subsampled[j][x] = (1.0 - l.translucent_v[x]) * subsampled[j][x] + l.translucent_v[x] * subsampled2[j][x];
								}
							}
						}
#endif
#ifdef HLRAD_AUTOCORING
						for (j = 0; j < ALLSTYLES && f_styles[j] != 255; j++)
#else
                        for (j = 0; j < MAXLIGHTMAPS && (f->styles[j] != 255); j++)
#endif
                        {
                            VectorScale(subsampled[j], weighting[s + 1][t + 1], subsampled[j]);
                            VectorAdd(sampled[j], subsampled[j], sampled[j]);
                        }
                        subsamples += weighting[s + 1][t + 1];
                    }
                }
            }
#ifdef HLRAD_AUTOCORING
			for (j = 0; j < ALLSTYLES && f_styles[j] != 255; j++)
#else
            for (j = 0; j < MAXLIGHTMAPS && (f->styles[j] != 255); j++)
#endif
            {
                VectorScale(sampled[j], 1.0 / subsamples, sampled[j]);
            }
        }
        else
        {
            GetPhongNormal(facenum, spot, pointnormal);
            GatherSampleLight(spot, pvs, pointnormal, sampled, 
#ifdef HLRAD_AUTOCORING
				f_styles
#else
				f->styles
#endif
#ifdef HLRAD_GatherPatchLight
				, 0
#endif
				);
#ifdef HLRAD_TRANSLUCENT
			if (l.translucent_b)
			{
#ifdef HLRAD_AUTOCORING
				vec3_t sampled2[ALLSTYLES];
				for (j = 0; j < ALLSTYLES; j++)
#else
				vec3_t sampled2[MAXLIGHTMAPS];
				for (j = 0; j < MAXLIGHTMAPS; j++)
#endif
				{
					VectorFill(sampled2[j], 0);
				}
				VectorSubtract (vec3_origin, pointnormal, normal2);
				GatherSampleLight(spot2, pvs2, normal2, sampled2, 
#ifdef HLRAD_AUTOCORING
					f_styles
#else
					f->styles
#endif
	#ifdef HLRAD_GatherPatchLight
					, 0
	#endif
					);
#ifdef HLRAD_AUTOCORING
				for (j = 0; j < ALLSTYLES && f_styles[j] != 255; j++)
#else
				for (j = 0; j < MAXLIGHTMAPS && (f->styles[j] != 255); j++)
#endif
				{
					for (int x = 0; x < 3; x++)
					{
						sampled[j][x] = (1.0 - l.translucent_v[x]) * sampled[j][x] + l.translucent_v[x] * sampled2[j][x];
					}
				}
			}
#endif
        }

#ifdef HLRAD_AUTOCORING
		for (j = 0; j < ALLSTYLES && f_styles[j] != 255; j++)
#else
        for (j = 0; j < MAXLIGHTMAPS && (f->styles[j] != 255); j++)
#endif
        {
#ifdef HLRAD_AUTOCORING
			VectorCopy (sampled[j], fl_samples[j][i].light);
#else
            VectorCopy(sampled[j], facelight[facenum].samples[j][i].light);
#endif

#ifndef HLRAD_TRANCPARENCYLOSS_FIX
#ifdef HLRAD_HULLU            
            if(b_transparency_loss)
            {
#ifdef HLRAD_AUTOCORING
				VectorScale (fl_samples[j][i].light, light_left_for_facelight, fl_samples[j][i].light);
#else
            	 VectorScale(facelight[facenum].samples[j][i].light, light_left_for_facelight, facelight[facenum].samples[j][i].light);
#endif
            }
#endif
#endif

#ifdef ZHLT_TEXLIGHT
	#ifdef HLRAD_AUTOCORING
			AddSampleToPatch (&fl_samples[j][i], facenum, f_styles[j]
	#else
            AddSampleToPatch(&facelight[facenum].samples[j][i], facenum, f->styles[j]
	#endif
	#ifdef HLRAD_AddSampleToPatch_PRECISE
				, l.surfpt_original[i]
	#endif
				); //LRC
#else
            if (f->styles[j] == 0)
            {
                AddSampleToPatch(&facelight[facenum].samples[j][i], facenum
	#ifdef HLRAD_AddSampleToPatch_PRECISE
					, l.surfpt_original[i]
	#endif
					);
            }
#endif
        }
    } // end of i loop

    // average up the direct light on each patch for radiosity
#ifndef HLRAD_GatherPatchLight
    if (g_numbounce > 0)
#endif
    {
        for (patch = g_face_patches[facenum]; patch; patch = patch->next)
        {
#ifdef ZHLT_TEXLIGHT
            //LRC:
			unsigned istyle;
	#ifdef HLRAD_AUTOCORING
			if (patch->samples)
			{
				for (istyle = 0; istyle < ALLSTYLES && patch->totalstyle_all[istyle] != 255; istyle++)
				{
					vec3_t v;
					VectorScale (patch->samplelight_all[istyle], 1.0f / patch->samples, v);
					VectorAdd (patch->directlight_all[istyle], v, patch->directlight_all[istyle]);
				}
			}
	#else
			for (istyle = 0; istyle < MAXLIGHTMAPS && patch->totalstyle[istyle] != 255; istyle++)
			{
				if (patch->samples[istyle])
		        {
		            vec3_t          v;                         // BUGBUG: Use a weighted average instead?

					VectorScale(patch->samplelight[istyle], (1.0f / patch->samples[istyle]), v);
					VectorAdd(patch->totallight[istyle], v, patch->totallight[istyle]);
	                VectorAdd(patch->directlight[istyle], v, patch->directlight[istyle]);
				}
			}
	#endif
            //LRC (ends)
#else
            if (patch->samples)
            {
                vec3_t          v;                         // BUGBUG: Use a weighted average instead?

                VectorScale(patch->samplelight, (1.0f / patch->samples), v);
                VectorAdd(patch->totallight, v, patch->totallight);
                VectorAdd(patch->directlight, v, patch->directlight);
            }
#endif
        }
    }
#ifdef HLRAD_GatherPatchLight
	for (patch = g_face_patches[facenum]; patch; patch = patch->next)
	{
		// get the PVS for the pos to limit the number of checks
		if (!g_visdatasize)
		{
			memset(pvs, 255, (g_numleafs + 7) / 8);
			lastoffset = -1;
		}
		else
		{
			dleaf_t*        leaf = PointInLeaf(patch->origin);

			thisoffset = leaf->visofs;
			if (l.numsurfpt == 0 || thisoffset != lastoffset)
			{
	#ifdef HLRAD_VIS_FIX
				if (thisoffset == -1)
				{
					memset(pvs, 0, (g_numleafs + 7) / 8);
				}
				else
				{
					DecompressVis(&g_dvisdata[leaf->visofs], pvs, sizeof(pvs));
				}
	#else
				hlassert(thisoffset != -1);
				DecompressVis(&g_dvisdata[leaf->visofs], pvs, sizeof(pvs));
	#endif
			}
			lastoffset = thisoffset;
		}
#ifdef HLRAD_TRANSLUCENT
		if (l.translucent_b)
		{
			if (!g_visdatasize)
			{
				memset(pvs2, 255, (g_numleafs + 7) / 8);
				lastoffset2 = -1;
			}
			else
			{
				VectorMA (patch->origin, -(g_translucentdepth+2*PATCH_HUNT_OFFSET), l.facenormal, spot2);
				dleaf_t*        leaf2 = PointInLeaf(spot2);

				thisoffset2 = leaf2->visofs;
				if (l.numsurfpt == 0 || thisoffset2 != lastoffset2)
				{
	#ifdef HLRAD_VIS_FIX
					if (thisoffset2 == -1)
					{
						memset(pvs2, 0, (g_numleafs + 7) / 8);
					}
					else
					{
						DecompressVis(&g_dvisdata[leaf2->visofs], pvs2, sizeof(pvs2));
					}
	#else
					hlassert(thisoffset2 != -1);
					DecompressVis(&g_dvisdata[leaf2->visofs], pvs2, sizeof(pvs2));
	#endif
				}
				lastoffset2 = thisoffset2;
			}
	#ifdef HLRAD_AUTOCORING
			vec3_t frontsampled[ALLSTYLES], backsampled[ALLSTYLES];
			for (j = 0; j < ALLSTYLES; j++)
	#else
			vec3_t frontsampled[MAXLIGHTMAPS], backsampled[MAXLIGHTMAPS];
			for (j = 0; j < MAXLIGHTMAPS; j++)
	#endif
			{
				VectorClear (frontsampled[j]);
				VectorClear (backsampled[j]);
			}
			VectorSubtract (vec3_origin, l.facenormal, normal2);
			GatherSampleLight (patch->origin, pvs, l.facenormal, frontsampled, 
	#ifdef HLRAD_AUTOCORING
				patch->totalstyle_all
	#else
				patch->totalstyle
	#endif
				, 1);
			GatherSampleLight (spot2, pvs2, normal2, backsampled, 
	#ifdef HLRAD_AUTOCORING
				patch->totalstyle_all
	#else
				patch->totalstyle
	#endif
				, 1);
	#ifdef HLRAD_AUTOCORING
			for (j = 0; j < ALLSTYLES && patch->totalstyle_all[j] != 255; j++)
	#else
			for (j = 0; j < MAXLIGHTMAPS && (patch->totalstyle[j] != 255); j++)
	#endif
			{
				for (int x = 0; x < 3; x++)
				{
	#ifdef HLRAD_AUTOCORING
					patch->totallight_all[j][x] += (1.0 - l.translucent_v[x]) * frontsampled[j][x] + l.translucent_v[x] * backsampled[j][x];
	#else
					patch->totallight[j][x] += (1.0 - l.translucent_v[x]) * frontsampled[j][x] + l.translucent_v[x] * backsampled[j][x];
	#endif
				}
			}
		}
		else
		{
			GatherSampleLight (patch->origin, pvs, l.facenormal, 
	#ifdef HLRAD_AUTOCORING
				patch->totallight_all, patch->totalstyle_all
	#else
				patch->totallight, patch->totalstyle
	#endif
				, 1);
		}
#else
		GatherSampleLight (patch->origin, pvs, l.facenormal, 
	#ifdef HLRAD_AUTOCORING
			patch->totallight_all, patch->totalstyle_all
	#else
			patch->totallight, patch->totalstyle
	#endif
			, 1);
#endif
	}
#endif

    // add an ambient term if desired
    if (g_ambient[0] || g_ambient[1] || g_ambient[2])
    {
#ifdef HLRAD_AUTOCORING
		for (j = 0; j < ALLSTYLES && f_styles[j] != 255; j++)
		{
			if (f_styles[j] == 0)
			{
				s = fl_samples[j];
#else
        for (j = 0; j < MAXLIGHTMAPS && f->styles[j] != 255; j++)
        {
            if (f->styles[j] == 0)
            {
                s = facelight[facenum].samples[j];
#endif
                for (i = 0; i < l.numsurfpt; i++, s++)
                {
                    VectorAdd(s->light, g_ambient, s->light);
                }
                break;
            }
        }

    }

    // add circus lighting for finding black lightmaps
    if (g_circus)
    {
#ifdef HLRAD_AUTOCORING
		for (j = 0; j < ALLSTYLES && f_styles[j] != 255; j++)
		{
			if (f_styles[j] == 0)
			{
#else
        for (j = 0; j < MAXLIGHTMAPS && f->styles[j] != 255; j++)
        {
            if (f->styles[j] == 0)
            {
#endif
                int             amt = 7;

#ifdef HLRAD_AUTOCORING
				s = fl_samples[j];
#else
                s = facelight[facenum].samples[j];
#endif

                while ((l.numsurfpt % amt) == 0)
                {
                    amt--;
                }
                if (amt < 2)
                {
                    amt = 7;
                }

                for (i = 0; i < l.numsurfpt; i++, s++)
                {
                    if ((s->light[0] == 0) && (s->light[1] == 0) && (s->light[2] == 0))
                    {
                        VectorAdd(s->light, s_circuscolors[i % amt], s->light);
                    }
                }
                break;
            }
        }
    }

    // light from dlight_threshold and above is sent out, but the
    // texture itself should still be full bright

    // if( VectorAvg( face_patches[facenum]->baselight ) >= dlight_threshold)       // Now all lighted surfaces glow
    {
#ifdef ZHLT_TEXLIGHT
        //LRC:
		if (g_face_patches[facenum])
		{
	#ifdef HLRAD_AUTOCORING
			for (j = 0; j < ALLSTYLES && f_styles[j] != 255; j++)
			{
				if (f_styles[j] == g_face_patches[facenum]->emitstyle)
				{
					break;
				}
			}
			if (j == ALLSTYLES)
	#else
			for (j = 0; j < MAXLIGHTMAPS && f->styles[j] != 255; j++)
			{
                if (f->styles[j] == g_face_patches[facenum]->emitstyle) //LRC
				{
					break;
				}
			}

			if (j == MAXLIGHTMAPS)
	#endif
			{
#ifdef HLRAD_READABLE_EXCEEDSTYLEWARNING
				if (++stylewarningcount >= stylewarningnext)
				{
					stylewarningnext = stylewarningcount * 2;
					Warning("Too many direct light styles on a face(?,?,?)");
					Warning(" total %d warnings for too many styles", stylewarningcount);
				}
#else
				Warning("Too many direct light styles on a face(?,?,?)");
#endif
			}
			else
			{
	#ifdef HLRAD_AUTOCORING
				if (f_styles[j] == 255)
				{
					f_styles[j] = g_face_patches[facenum]->emitstyle;
				}

				s = fl_samples[j];
	#else
				if (f->styles[j] == 255)
				{
					f->styles[j] = g_face_patches[facenum]->emitstyle;
				}

				s = facelight[facenum].samples[j];
	#endif
				for (i = 0; i < l.numsurfpt; i++, s++)
				{
					VectorAdd(s->light, g_face_patches[facenum]->baselight, s->light);
				}
			}
		}
        //LRC (ends)
#else
        for (j = 0; j < MAXLIGHTMAPS && f->styles[j] != 255; j++)
        {
            if (f->styles[j] == 0)
            {
                if (g_face_patches[facenum])
                {
                    s = facelight[facenum].samples[j];
                    for (i = 0; i < l.numsurfpt; i++, s++)
                    {
                        VectorAdd(s->light, g_face_patches[facenum]->baselight, s->light);
                    }
                    break;
                }
            }
        }
#endif
    }
#ifdef HLRAD_AUTOCORING
	{
		facelight_t *fl = &facelight[facenum];
		vec_t maxlights[ALLSTYLES];
		for (j = 0; j < ALLSTYLES && f_styles[j] != 255; j++)
		{
			maxlights[j] = 0;
			for (i = 0; i < fl->numsamples; i++)
			{
				vec_t b = VectorMaximum (fl_samples[j][i].light);
				maxlights[j] = max (maxlights[j], b);
			}
			if (maxlights[j] < g_corings[f_styles[j]] * 0.5)
			{
				maxlights[j] = 0;
			}
		}
		for (k = 0; k < MAXLIGHTMAPS; k++)
		{
			int bestindex = -1;
			if (k == 0)
			{
				bestindex = 0;
			}
			else
			{
				vec_t bestmaxlight = 0;
				for (j = 1; j < ALLSTYLES && f_styles[j] != 255; j++)
				{
					if (maxlights[j] > bestmaxlight + NORMAL_EPSILON)
					{
						bestmaxlight = maxlights[j];
						bestindex = j;
					}
				}
			}
			if (bestindex != -1)
			{
				maxlights[bestindex] = 0;
				f->styles[k] = f_styles[bestindex];
				fl->samples[k] = (sample_t *)malloc (fl->numsamples * sizeof (sample_t));
				hlassume (fl->samples[k] != NULL, assume_NoMemory);
				memcpy (fl->samples[k], fl_samples[bestindex], fl->numsamples * sizeof (sample_t));
			}
			else
			{
				f->styles[k] = 255;
				fl->samples[k] = NULL;
			}
		}
		for (j = 0; j < ALLSTYLES; j++)
		{
			free (fl_samples[j]);
		}
	}
	for (patch = g_face_patches[facenum]; patch; patch = patch->next)
	{
		vec_t maxlights[ALLSTYLES];
		for (j = 0; j < ALLSTYLES && patch->totalstyle_all[j] != 255; j++)
		{
			maxlights[j] = VectorMaximum (patch->totallight_all[j]);
		}
		for (k = 0; k < MAXLIGHTMAPS; k++)
		{
			int bestindex = -1;
			if (k == 0)
			{
				bestindex = 0;
			}
			else
			{
				vec_t bestmaxlight = 0;
				for (j = 1; j < ALLSTYLES && patch->totalstyle_all[j] != 255; j++)
				{
					if (maxlights[j] > bestmaxlight + NORMAL_EPSILON)
					{
						bestmaxlight = maxlights[j];
						bestindex = j;
					}
				}
			}
			if (bestindex != -1)
			{
				maxlights[bestindex] = 0;
				patch->totalstyle[k] = patch->totalstyle_all[bestindex];
				VectorCopy (patch->totallight_all[bestindex], patch->totallight[k]);
			}
			else
			{
				patch->totalstyle[k] = 255;
			}
		}
		for (j = 0; j < ALLSTYLES && patch->totalstyle_all[j] != 255; j++)
		{
			maxlights[j] = VectorMaximum (patch->directlight_all[j]);
		}
		for (k = 0; k < MAXLIGHTMAPS; k++)
		{
			int bestindex = -1;
			if (k == 0)
			{
				bestindex = 0;
			}
			else
			{
				vec_t bestmaxlight = 0;
				for (j = 1; j < ALLSTYLES && patch->totalstyle_all[j] != 255; j++)
				{
					if (maxlights[j] > bestmaxlight + NORMAL_EPSILON)
					{
						bestmaxlight = maxlights[j];
						bestindex = j;
					}
				}
			}
			if (bestindex != -1)
			{
				maxlights[bestindex] = 0;
				patch->directstyle[k] = patch->totalstyle_all[bestindex];
				VectorCopy (patch->directlight_all[bestindex], patch->directlight[k]);
			}
			else
			{
				patch->directstyle[k] = 255;
			}
		}
		free (patch->totalstyle_all);
		patch->totalstyle_all = NULL;
		free (patch->samplelight_all);
		patch->samplelight_all = NULL;
		free (patch->totallight_all);
		patch->totallight_all = NULL;
		free (patch->directlight_all);
		patch->directlight_all = NULL;
	}
#endif
}

// =====================================================================================
//  PrecompLightmapOffsets
// =====================================================================================
void            PrecompLightmapOffsets()
{
    int             facenum;
    dface_t*        f;
    facelight_t*    fl;
    int             lightstyles;

#ifdef ZHLT_TEXLIGHT
    int             i; //LRC
	patch_t*        patch; //LRC
#endif

    g_lightdatasize = 0;

    for (facenum = 0; facenum < g_numfaces; facenum++)
    {
        f = &g_dfaces[facenum];
        fl = &facelight[facenum];

        if (g_texinfo[f->texinfo].flags & TEX_SPECIAL)
        {
            continue;                                      // non-lit texture
        }

#ifdef HLRAD_ENTSTRIPRAD
		if (IntForKey (g_face_entity[facenum], "zhlt_striprad"))
		{
			continue;
		}
#endif

#ifdef HLRAD_AUTOCORING
		{
			int i, j, k;
			vec_t maxlights[ALLSTYLES];
			{
				vec3_t maxlights1[ALLSTYLES];
				vec3_t maxlights2[ALLSTYLES];
				for (j = 0; j < ALLSTYLES; j++)
				{
					VectorClear (maxlights1[j]);
					VectorClear (maxlights2[j]);
				}
				for (k = 0; k < MAXLIGHTMAPS && f->styles[k] != 255; k++)
				{
					for (i = 0; i < fl->numsamples; i++)
					{
						VectorCompareMaximum (maxlights1[f->styles[k]], fl->samples[k][i].light, maxlights1[f->styles[k]]);
					}
				}
				for (patch = g_face_patches[facenum]; patch; patch = patch->next)
				{
					for (k = 0; k < MAXLIGHTMAPS && patch->totalstyle[k] != 255; k++)
					{
						VectorCompareMaximum (maxlights2[patch->totalstyle[k]], patch->totallight[k], maxlights2[patch->totalstyle[k]]);
					}
				}
				for (j = 0; j < ALLSTYLES; j++)
				{
					vec3_t v;
					VectorAdd (maxlights1[j], maxlights2[j], v);
					maxlights[j] = VectorMaximum (v);
					if (maxlights[j] < g_corings[j])
					{
						maxlights[j] = 0;
					}
				}
			}
			unsigned char oldstyles[MAXLIGHTMAPS];
			sample_t *oldsamples[MAXLIGHTMAPS];
			for (k = 0; k < MAXLIGHTMAPS; k++)
			{
				oldstyles[k] = f->styles[k];
				oldsamples[k] = fl->samples[k];
			}
			for (k = 0; k < MAXLIGHTMAPS; k++)
			{
				unsigned char beststyle = 255;
				if (k == 0)
				{
					beststyle = 0;
				}
				else
				{
					vec_t bestmaxlight = 0;
					for (j = 1; j < ALLSTYLES; j++)
					{
						if (maxlights[j] > bestmaxlight + NORMAL_EPSILON)
						{
							bestmaxlight = maxlights[j];
							beststyle = j;
						}
					}
				}
				if (beststyle != 255)
				{
					maxlights[beststyle] = 0;
					f->styles[k] = beststyle;
					fl->samples[k] = (sample_t *)malloc (fl->numsamples * sizeof (sample_t));
					hlassume (fl->samples[k] != NULL, assume_NoMemory);
					for (i = 0; i < MAXLIGHTMAPS && oldstyles[i] != 255; i++)
					{
						if (oldstyles[i] == f->styles[k])
						{
							break;
						}
					}
					if (i < MAXLIGHTMAPS && oldstyles[i] != 255)
					{
						memcpy (fl->samples[k], oldsamples[i], fl->numsamples * sizeof (sample_t));
					}
					else
					{
						memcpy (fl->samples[k], oldsamples[0], fl->numsamples * sizeof (sample_t));
						for (j = 0; j < fl->numsamples; j++)
						{
							VectorClear (fl->samples[k][j].light);
						}
					}
				}
				else
				{
					f->styles[k] = 255;
					fl->samples[k] = NULL;
				}
			}
			for (k = 0; k < MAXLIGHTMAPS && oldstyles[k] != 255; k++)
			{
				free (oldsamples[k]);
			}
		}
#else
#ifdef ZHLT_TEXLIGHT
        		//LRC - find all the patch lightstyles, and add them to the ones used by this face
#ifdef HLRAD_STYLE_CORING
		for (patch = g_face_patches[facenum]; patch; patch = patch->next)
#else
		patch = g_face_patches[facenum];
		if (patch)
#endif
		{
			for (i = 0; i < MAXLIGHTMAPS && patch->totalstyle[i] != 255; i++)
			{
				for (lightstyles = 0; lightstyles < MAXLIGHTMAPS && f->styles[lightstyles] != 255; lightstyles++)
				{
					if (f->styles[lightstyles] == patch->totalstyle[i])
						break;
				}
				if (lightstyles == MAXLIGHTMAPS)
				{
#ifdef HLRAD_READABLE_EXCEEDSTYLEWARNING
					if (++stylewarningcount >= stylewarningnext)
					{
						stylewarningnext = stylewarningcount * 2;
						Warning("Too many direct light styles on a face(?,?,?)\n");
						Warning(" total %d warnings for too many styles", stylewarningcount);
					}
#else
					Warning("Too many direct light styles on a face(?,?,?)\n");
#endif
				}
				else if (f->styles[lightstyles] == 255)
				{
					f->styles[lightstyles] = patch->totalstyle[i];
//					Log("Face acquires new lightstyle %d at offset %d\n", f->styles[lightstyles], lightstyles);
				}
			}
		}
		//LRC (ends)
#endif
#endif

        for (lightstyles = 0; lightstyles < MAXLIGHTMAPS; lightstyles++)
        {
            if (f->styles[lightstyles] == 255)
            {
                break;
            }
        }

        if (!lightstyles)
        {
            continue;
        }

        f->lightofs = g_lightdatasize;
        g_lightdatasize += fl->numsamples * 3 * lightstyles;
		hlassume (g_lightdatasize <= g_max_map_lightdata, assume_MAX_MAP_LIGHTING); //lightdata

    }
}
#ifdef HLRAD_REDUCELIGHTMAP
void ReduceLightmap ()
{
	byte *oldlightdata = (byte *)malloc (g_lightdatasize);
	hlassume (oldlightdata != NULL, assume_NoMemory);
	memcpy (oldlightdata, g_dlightdata, g_lightdatasize);
	g_lightdatasize = 0;

	int facenum;
	for (facenum = 0; facenum < g_numfaces; facenum++)
	{
		dface_t *f = &g_dfaces[facenum];
		facelight_t *fl = &facelight[facenum];
		if (g_texinfo[f->texinfo].flags & TEX_SPECIAL)
		{
			continue;                                      // non-lit texture
		}
#ifdef HLRAD_ENTSTRIPRAD
		if (IntForKey (g_face_entity[facenum], "zhlt_striprad"))
		{
			continue;
		}
#endif
		if (f->lightofs == -1)
		{
			continue;
		}

		int i, j, k;
		int oldofs;
		unsigned char oldstyles[MAXLIGHTMAPS];
		oldofs = f->lightofs;
		f->lightofs = g_lightdatasize;
		for (k = 0; k < MAXLIGHTMAPS; k++)
		{
			oldstyles[k] = f->styles[k];
			f->styles[k] = 255;
		}
		int numstyles = 0;
		for (k = 0; k < MAXLIGHTMAPS && oldstyles[k] != 255; k++)
		{
			unsigned char maxb = 0;
			for (i = 0; i < fl->numsamples; i++)
			{
				unsigned char *v = &oldlightdata[oldofs + fl->numsamples * 3 * k + i * 3];
				maxb = max (maxb, VectorMaximum (v));
			}
			if (maxb <= 1) // very dark
			{
				continue;
			}
			f->styles[numstyles] = oldstyles[k];
			hlassume (g_lightdatasize + fl->numsamples * 3 * (numstyles + 1) <= g_max_map_lightdata, assume_MAX_MAP_LIGHTING);
			memcpy (&g_dlightdata[f->lightofs + fl->numsamples * 3 * numstyles], &oldlightdata[oldofs + fl->numsamples * 3 * k], fl->numsamples * 3);
			numstyles++;
		}
		g_lightdatasize += fl->numsamples * 3 * numstyles;
	}
	free (oldlightdata);
}
#endif

#ifdef HLRAD_MDL_LIGHT_HACK

// Change the sample light right under a mdl file entity's origin.
// Use this when "mdl" in shadow has incorrect brightness.

const int MLH_MAXFACECOUNT = 16;
const int MLH_MAXSAMPLECOUNT = 4;
const vec_t MLH_LEFT = 0;
const vec_t MLH_RIGHT = 1;

typedef struct
{
	vec3_t origin;
	vec3_t floor;
	struct
	{
		int num;
		struct
		{
			bool exist;
			int seq;
		}
		style[ALLSTYLES];
		struct
		{
			int num;
			vec3_t pos;
			unsigned char* (style[ALLSTYLES]);
		}
		sample[MLH_MAXSAMPLECOUNT];
		int samplecount;
	}
	face[MLH_MAXFACECOUNT];
	int facecount;
} mdllight_t;

#ifdef HLRAD_MDL_LIGHT_HACK_NEW
int MLH_AddFace (mdllight_t *ml, int facenum)
{
	dface_t *f = &g_dfaces[facenum];
	int i, j;
	for (i = 0; i < ml->facecount; i++)
	{
		if (ml->face[i].num == facenum)
		{
			return -1;
		}
	}
	if (ml->facecount >= MLH_MAXFACECOUNT)
	{
		return -1;
	}
	i = ml->facecount;
	ml->facecount++;
	ml->face[i].num = facenum;
	ml->face[i].samplecount = 0;
	for (j = 0; j < ALLSTYLES; j++)
	{
		ml->face[i].style[j].exist = false;
	}
	for (j = 0; j < MAXLIGHTMAPS && f->styles[j] != 255; j++)
	{
		ml->face[i].style[f->styles[j]].exist = true;
		ml->face[i].style[f->styles[j]].seq = j;
	}
	return i;
}
void MLH_AddSample (mdllight_t *ml, int facenum, int w, int h, int s, int t, const vec3_t pos)
{
	dface_t *f = &g_dfaces[facenum];
	int i, j;
	int r = MLH_AddFace (ml, facenum);
	if (r == -1)
	{
		return;
	}
	int size = w * h;
	int num = s + w * t;
	for (i = 0; i < ml->face[r].samplecount; i++)
	{
		if (ml->face[r].sample[i].num == num)
		{
			return;
		}
	}
	if (ml->face[r].samplecount >= MLH_MAXSAMPLECOUNT)
	{
		return;
	}
	i = ml->face[r].samplecount;
	ml->face[r].samplecount++;
	ml->face[r].sample[i].num = num;
	VectorCopy (pos, ml->face[r].sample[i].pos);
	for (j = 0; j < ALLSTYLES; j++)
	{
		if (ml->face[r].style[j].exist)
		{
			ml->face[r].sample[i].style[j] = &g_dlightdata[f->lightofs + (num + size * ml->face[r].style[j].seq) * 3];
		}
	}
}
void MLH_CalcExtents (const dface_t *f, int *texturemins, int *extents)
{
	float mins[2], maxs[2];
	int bmins[2], bmaxs[2];
	texinfo_t *tex;
	tex = &g_texinfo[f->texinfo];
	mins[0] = mins[1] = 999999;
	maxs[0] = maxs[1] = -99999;
	int i;
	for (i = 0; i < f->numedges; i++)
	{
		int e;
		dvertex_t *v;
		int j;
		e = g_dsurfedges[f->firstedge + i];
		if (e >= 0)
		{
			v = &g_dvertexes[g_dedges[e].v[0]];
		}
		else
		{
			v = &g_dvertexes[g_dedges[-e].v[1]];
		}
		for (j = 0; j < 2; j++)
		{
			float val = v->point[0] * tex->vecs[j][0] + v->point[1] * tex->vecs[j][1]
				+ v->point[2] * tex->vecs[j][2] + tex->vecs[j][3];
			if (val < mins[j])
			{
				mins[j] = val;
			}
			if (val > maxs[j])
			{
				maxs[j] = val;
			}
		}
	}
	for (i = 0; i < 2; i++)
	{
		bmins[i] = floor (mins[i] / 16);
		bmaxs[i] = ceil (maxs[i] / 16);
		texturemins[i] = bmins[i] * 16;
		extents[i] = (bmaxs[i] - bmins[i]) * 16;
	}
}
void MLH_GetSamples_r (mdllight_t *ml, int nodenum, const float *start, const float *end)
{
	if (nodenum < 0)
		return;
	dnode_t *node = &g_dnodes[nodenum];
	dplane_t *plane;
	float front, back, frac;
	float mid[3];
	int side;
	plane = &g_dplanes[node->planenum];
	front = DotProduct (start, plane->normal) - plane->dist;
	back = DotProduct (end, plane->normal) - plane->dist;
	side = front < 0;
	if ((back < 0) == side)
	{
		MLH_GetSamples_r (ml, node->children[side], start, end);
		return;
	}
	frac = front / (front - back);
	mid[0] = start[0] + (end[0] - start[0]) * frac;
	mid[1] = start[1] + (end[1] - start[1]) * frac;
	mid[2] = start[2] + (end[2] - start[2]) * frac;
	MLH_GetSamples_r (ml, node->children[side], start, mid);
	if (ml->facecount > 0)
	{
		return;
	}
	{
		int i;
		for (i = 0; i < node->numfaces; i++)
		{
			dface_t *f = &g_dfaces[node->firstface + i];
			texinfo_t *tex = &g_texinfo[f->texinfo];
			const char *texname = GetTextureByNumber (f->texinfo);
			if (!strncmp (texname, "sky", 3))
			{
				continue;
			}
			if (f->lightofs == -1)
			{
				continue;
			}
			int s = DotProduct (mid, tex->vecs[0]) + tex->vecs[0][3];
			int t = DotProduct (mid, tex->vecs[1]) + tex->vecs[1][3];
			int texturemins[2], extents[2];
			MLH_CalcExtents (f, texturemins, extents);
			if (s < texturemins[0] || t < texturemins[1])
			{
				continue;
			}
			int ds = s - texturemins[0];
			int dt = t - texturemins[1];
			if (ds > extents[0] || dt > extents[1])
			{
				continue;
			}
			ds >>= 4;
			dt >>= 4;
			MLH_AddSample (ml, node->firstface + i, extents[0] / 16 + 1, extents[1] / 16 + 1, ds, dt, mid);
			break;
		}
	}
	if (ml->facecount > 0)
	{
		VectorCopy (mid, ml->floor);
		return;
	}
	MLH_GetSamples_r (ml, node->children[!side], mid, end);
}
void MLH_mdllightCreate (mdllight_t *ml)
{
	// code from Quake
	float p[3];
	float end[3];
	ml->facecount = 0;
	VectorCopy (ml->origin, ml->floor);
	VectorCopy (ml->origin, p);
	VectorCopy (ml->origin, end);
	end[2] -= 2048;
	MLH_GetSamples_r (ml, 0, p, end);
}
#else
void MLH_mdllightCreate (mdllight_t *ml)
{
	int i, j, k;
	vec_t height, minheight = BOGUS_RANGE;
	ml->facecount = 0;
	for (i = 0; i < g_numfaces; ++i)
	{
		if (stricmp (ValueForKey (g_face_entity[i], "classname"), "worldspawn"))
			continue;
		const dface_t *f = &g_dfaces[i];
		const dplane_t *p = getPlaneFromFace(f);
		Winding *w=new Winding (*f);
		for (j = 0;j < w->m_NumPoints; j++)
		{
			VectorAdd(w->m_Points[j], g_face_offset[i], w->m_Points[j]);
		}
		vec3_t delta , sect;
		VectorCopy (ml->origin, delta);
		delta[2] -= BOGUS_RANGE;
		if (intersect_linesegment_plane(p, ml->origin, delta, sect) && point_in_winding (*w, *p, sect))
		{
			height = ml->origin[2] - sect[2];
			if (height >= 0 && height <= minheight)
				minheight = height;
		}
		delete w;
	}
	VectorCopy (ml->origin, ml->floor);
	ml->floor[2] -= minheight;
	for (i = 0; i < g_numfaces; ++i)
	{
		if (stricmp (ValueForKey (g_face_entity[i], "classname"), "worldspawn"))
			continue;
		const dface_t *f = &g_dfaces[i];
		const dplane_t *p = getPlaneFromFace(f);
		Winding *w=new Winding (*f);
		if (g_texinfo[f->texinfo].flags & TEX_SPECIAL)
		{
			continue;                                            // non-lit texture
		}
		for (j = 0;j < w->m_NumPoints; j++)
		{
			VectorAdd(w->m_Points[j], g_face_offset[i], w->m_Points[j]);
		}
		vec3_t delta , sect;
		VectorCopy (ml->origin, delta);
		delta[2] -= BOGUS_RANGE;
		if (intersect_linesegment_plane(p, ml->origin, delta, sect) && VectorCompare (sect, ml->floor))
		{
			bool inlightmap = false;
			{
				vec3_t v;
				facesampleinfo_t *info = &facesampleinfo[i];
				int w = info->texsize[0] + 1;
				int h = info->texsize[1] + 1;
				vec_t vs, vt;
				int s1, s2, t1, t2, s, t;
				VectorCopy (ml->floor, v);
				VectorSubtract (v, info->offset, v);
				VectorSubtract (v, info->texorg, v);
				vs = DotProduct (v, info->worldtotex[0]);
				vt = DotProduct (v, info->worldtotex[1]);
				s1 = (int)floor((vs-MLH_LEFT)/16) - info->texmins[0];
				s2 = (int)floor((vs+MLH_RIGHT)/16) - info->texmins[0];
				t1 = (int)floor((vt-MLH_LEFT)/16) - info->texmins[1];
				t2 = (int)floor((vt+MLH_RIGHT)/16) - info->texmins[1];
				for (s=s1; s<=s2; ++s)
					for (t=t1; t<=t2; ++t)
						if (s>=0 && s<w && t>=0 && t<h)
							inlightmap = true;
			}
			if (inlightmap && ml->facecount < MLH_MAXFACECOUNT)
			{
				ml->face[ml->facecount].num = i;
				ml->facecount++;
			}
		}
		delete w;
	}
	for (i = 0; i < ml->facecount; ++i)
	{
		const dface_t *f = &g_dfaces[ml->face[i].num];
		for (j = 0; j < ALLSTYLES; ++j)
			ml->face[i].style[j].exist = false;
		for (j = 0; j < MAXLIGHTMAPS && f->styles[j] != 255; ++j)
		{
			ml->face[i].style[f->styles[j]].exist = true;
			ml->face[i].style[f->styles[j]].seq = j;
		}
		ml->face[i].samplecount = 0;
		if (j == 0)
			continue;

	    const facelight_t *fl=&facelight[ml->face[i].num];
		{
			vec3_t v;
			facesampleinfo_t *info = &facesampleinfo[ml->face[i].num];
			int w = info->texsize[0] + 1;
			int h = info->texsize[1] + 1;
			vec_t vs, vt;
			int s1, s2, t1, t2, s, t;
			VectorCopy (ml->floor, v);
			VectorSubtract (v, info->offset, v);
			VectorSubtract (v, info->texorg, v);
			vs = DotProduct (v, info->worldtotex[0]);
			vt = DotProduct (v, info->worldtotex[1]);
			s1 = (int)floor((vs-MLH_LEFT)/16) - info->texmins[0];
			s2 = (int)floor((vs+MLH_RIGHT)/16) - info->texmins[0];
			t1 = (int)floor((vt-MLH_LEFT)/16) - info->texmins[1];
			t2 = (int)floor((vt+MLH_RIGHT)/16) - info->texmins[1];
			for (s=s1; s<=s2; ++s)
				for (t=t1; t<=t2; ++t)
					if (s>=0 && s<w && t>=0 && t<h)
						if (ml->face[i].samplecount < MLH_MAXSAMPLECOUNT)
						{
							ml->face[i].sample[ml->face[i].samplecount].num = s + t * w;
							VectorAdd (info->offset, info->texorg, v);
							vs = 16.0 * (s + info->texmins[0]);
							vt = 16.0 * (t + info->texmins[1]);
							VectorMA (v, vs, info->textoworld[0], v);
							VectorMA (v, vt, info->textoworld[1], v);
							VectorCopy (v, ml->face[i].sample[ml->face[i].samplecount].pos);
							ml->face[i].samplecount++;
						}
		}

		for (j = 0; j < ml->face[i].samplecount; ++j)
		{
			for (k = 0; k < ALLSTYLES; ++k)
				if (ml->face[i].style[k].exist)
				{
					ml->face[i].sample[j].style[k] = 
						&g_dlightdata[f->lightofs + ml->face[i].style[k].seq * fl->numsamples * 3 + ml->face[i].sample[j].num * 3];
				}
		}
	}
}
#endif

int MLH_CopyLight (const vec3_t from, const vec3_t to)
{
	int i, j, k, count = 0;
	mdllight_t mlfrom, mlto;
	VectorCopy (from, mlfrom.origin);
	VectorCopy (to, mlto.origin);
	MLH_mdllightCreate (&mlfrom);
	MLH_mdllightCreate (&mlto);
	if (mlfrom.facecount == 0 || mlfrom.face[0].samplecount == 0)
		return -1;
	for (i = 0; i < mlto.facecount; ++i)
		for (j = 0; j < mlto.face[i].samplecount; ++j, ++count)
			for (k = 0; k < ALLSTYLES; ++k)
				if (mlto.face[i].style[k].exist && mlfrom.face[0].style[k].exist)
				{
					VectorCopy (mlfrom.face[0].sample[0].style[k],mlto.face[i].sample[j].style[k]);
					Developer (DEVELOPER_LEVEL_SPAM, "Mdl Light Hack: face (%d) sample (%d) style (%d) position (%f,%f,%f)\n",
						mlto.face[i].num, mlto.face[i].sample[j].num, k, 
						mlto.face[i].sample[j].pos[0], mlto.face[i].sample[j].pos[1], mlto.face[i].sample[j].pos[2]);
				}
	Developer (DEVELOPER_LEVEL_MESSAGE, "Mdl Light Hack: %d sample light copied from (%f,%f,%f) to (%f,%f,%f)\n", 
		count, mlfrom.floor[0], mlfrom.floor[1], mlfrom.floor[2], mlto.floor[0], mlto.floor[1], mlto.floor[2]);
	return count;
}

void MdlLightHack ()
{
	int ient;
	entity_t *ent1, *ent2;
	vec3_t origin1, origin2;
	const char *target;
#ifndef HLRAD_MDL_LIGHT_HACK_NEW
    double start, end;
#endif
	int used = 0, countent = 0, countsample = 0, r;
#ifndef HLRAD_MDL_LIGHT_HACK_NEW
    start = I_FloatTime();
#endif
	for (ient = 0; ient < g_numentities; ++ient)
	{
		ent1 = &g_entities[ient];
		target = ValueForKey (ent1, "zhlt_copylight");
		if (!strcmp (target, ""))
			continue;
		used = 1;
		ent2 = FindTargetEntity (target);
		if (ent2 == NULL)
		{
			Warning ("target entity '%s' not found", target);
			continue;
		}
		GetVectorForKey (ent1, "origin", origin1);
		GetVectorForKey (ent2, "origin", origin2);
		r = MLH_CopyLight (origin2, origin1);
		if (r < 0)
			Warning ("can not copy light from (%f,%f,%f)", origin2[0], origin2[1], origin2[2]);
		else
		{
			countent += 1;
			countsample += r;
		}
	}
#ifndef HLRAD_MDL_LIGHT_HACK_NEW
    end = I_FloatTime();
#endif
	if (used)
#ifdef HLRAD_MDL_LIGHT_HACK_NEW
		Log ("Adjust mdl light: modified %d samples for %d entities\n", countsample, countent);
#else
		Log("Mdl Light Hack: %d entities %d samples (%.2f seconds)\n", countent, countsample, end - start);
#endif
}
#endif /*HLRAD_MDL_LIGHT_HACK*/

// =====================================================================================
//  FinalLightFace
//      Add the indirect lighting on top of the direct lighting and save into final map format
// =====================================================================================
void            FinalLightFace(const int facenum)
{
#ifdef HLRAD_DEBUG_DRAWPOINTS
	if (facenum == 0 && g_drawsample)
	{
		char name[_MAX_PATH+20];
		sprintf (name, "%s_sample.pts", g_Mapname);
		Log ("Writing '%s' ...\n", name);
		FILE *f;
		f = fopen(name, "w");
		if (f)
		{
			const int pos_count = 15;
			const vec3_t pos[pos_count] = {{0,0,0},{1,0,0},{0,1,0},{-1,0,0},{0,-1,0},{1,0,0},{0,0,1},{-1,0,0},{0,0,-1},{0,-1,0},{0,0,1},{0,1,0},{0,0,-1},{1,0,0},{0,0,0}};
			int i, j, k;
			vec3_t v, dist;
			for (i = 0; i < g_numfaces; ++i)
			{
				const facelight_t *fl=&facelight[i];
				for (j = 0; j < fl->numsamples; ++j)
				{
					VectorCopy (fl->samples[0][j].pos, v);
					VectorSubtract (v, g_drawsample_origin, dist);
					if (DotProduct (dist, dist) < g_drawsample_radius * g_drawsample_radius)
					{
						for (k = 0; k < pos_count; ++k)
							fprintf (f, "%g %g %g\n", v[0]+pos[k][0], v[1]+pos[k][1], v[2]+pos[k][2]);
					}
				}
			}
			fclose(f);
			Log ("OK.\n");
		}
		else
			Log ("Error.\n");
	}
#endif
    int             i, j, k;
    vec3_t          lb, v;
    facelight_t*    fl;
    sample_t*       samp;
    float           minlight;
    int             lightstyles;
    dface_t*        f;
    lerpTriangulation_t* trian = NULL;
#ifdef HLRAD_FinalLightFace_VL
	vec3_t			*original_basiclight;
	int				(*final_basiclight)[3];
	int				lbi[3];
#endif

    // ------------------------------------------------------------------------
    // Changes by Adam Foster - afoster@compsoc.man.ac.uk
#ifdef HLRAD_WHOME
    float           temp_rand;
#endif
    // ------------------------------------------------------------------------

    f = &g_dfaces[facenum];
    fl = &facelight[facenum];

    if (g_texinfo[f->texinfo].flags & TEX_SPECIAL)
    {
        return;                                            // non-lit texture
    }

#ifdef HLRAD_ENTSTRIPRAD
	if (IntForKey (g_face_entity[facenum], "zhlt_striprad"))
	{
		return;
	}
#endif

    for (lightstyles = 0; lightstyles < MAXLIGHTMAPS; lightstyles++)
    {
        if (f->styles[lightstyles] == 255)
        {
            break;
        }
    }

    if (!lightstyles)
    {
        return;
    }

    //
    // set up the triangulation
    //
#ifndef HLRAD_GatherPatchLight
    if (g_numbounce)
#endif
    {
        trian = CreateTriangulation(facenum);
    }
    //
    // sample the triangulation
    //
    minlight = FloatForKey(g_face_entity[facenum], "_minlight") * 128;

#ifdef HLRAD_FinalLightFace_VL
	original_basiclight = (vec3_t *)calloc (fl->numsamples, sizeof(vec3_t));
	final_basiclight = (int (*)[3])calloc (fl->numsamples, sizeof(int [3]));
	hlassume (original_basiclight != NULL, assume_NoMemory);
	hlassume (final_basiclight != NULL, assume_NoMemory);
#endif
    for (k = 0; k < lightstyles; k++)
    {
        samp = fl->samples[k];
        for (j = 0; j < fl->numsamples; j++, samp++)
        {
            // Should be a VectorCopy, but we scale by 2 to compensate for an earlier lighting flaw
            // Specifically, the directlight contribution was included in the bounced light AND the directlight
            // Since many of the levels were built with this assumption, this "fudge factor" compensates for it.

			// Default direct_scale has been changed from 2 to 1 and default scale has been changed from 1 to 2. --vluzacn
            VectorScale(samp->light, g_direct_scale, lb);

#ifdef ZHLT_TEXLIGHT
#ifndef HLRAD_GatherPatchLight
            if (g_numbounce)//LRC && (k == 0))
#endif
            {
                SampleTriangulation(trian, samp->pos, v, f->styles[k]); //LRC
#else
            if (
#ifndef HLRAD_GatherPatchLight
				g_numbounce &&
#endif
				(k == 0))
            {
                SampleTriangulation(trian, samp->pos, v);
#endif

                if (isPointFinite(v))
                {
#ifdef HLRAD_STYLE_CORING
					VectorAdd (lb, v, v);
					if (VectorMaximum (v) >= g_corings[f->styles[k]])
					{
						VectorCopy (v, lb);
					}
#else
                    VectorAdd(lb, v, lb);
#endif
                }
                else
                {
                    Warning("point (%4.3f %4.3f %4.3f) infinite v (%4.3f %4.3f %4.3f)\n",
                            samp->pos[0], samp->pos[1], samp->pos[2], v[0], v[1], v[2]);
                }


            }
#ifdef HLRAD_FinalLightFace_VL
			// Warning: always assume the first style is 0
			if (k == 0)
			{
				VectorCopy (lb, original_basiclight[j])
			}
			else
			{
				VectorAdd (lb, original_basiclight[j], lb);
			}
#endif
            // ------------------------------------------------------------------------
	        // Changes by Adam Foster - afoster@compsoc.man.ac.uk
	        // colour lightscale:
#ifdef HLRAD_WHOME
	        lb[0] *= g_colour_lightscale[0];
	        lb[1] *= g_colour_lightscale[1];
	        lb[2] *= g_colour_lightscale[2];
#else
            VectorScale(lb, g_lightscale, lb);
#endif
            // ------------------------------------------------------------------------

            // clip from the bottom first
            for (i = 0; i < 3; i++)
            {
                if (lb[i] < minlight)
                {
                    lb[i] = minlight;
                }
            }

#ifndef HLRAD_FinalLightFace_VL
            // clip from the top
            {
                vec_t           max = VectorMaximum(lb);

                if (max > g_maxlight)
                {
                    vec_t           scale = g_maxlight / max;

                    lb[0] *= scale;
                    lb[1] *= scale;
                    lb[2] *= scale;
                }
            }
#endif

	        // ------------------------------------------------------------------------
	        // Changes by Adam Foster - afoster@compsoc.man.ac.uk
#ifdef HLRAD_WHOME

            // AJM: your code is formatted really wierd, and i cant understand a damn thing. 
            //      so i reformatted it into a somewhat readable "normal" fashion. :P

	        if ( g_colour_qgamma[0] != 1.0 ) 
		        lb[0] = (float) pow(lb[0] / 256.0f, g_colour_qgamma[0]) * 256.0f;

	        if ( g_colour_qgamma[1] != 1.0 ) 
		        lb[1] = (float) pow(lb[1] / 256.0f, g_colour_qgamma[1]) * 256.0f;

	        if ( g_colour_qgamma[2] != 1.0 ) 
		        lb[2] = (float) pow(lb[2] / 256.0f, g_colour_qgamma[2]) * 256.0f;

	        // Two different ways of adding noise to the lightmap - colour jitter
	        // (red, green and blue channels are independent), and mono jitter
	        // (monochromatic noise). For simulating dithering, on the cheap. :)

	        // Tends to create seams between adjacent polygons, so not ideal.

	        // Got really weird results when it was set to limit values to 256.0f - it
	        // was as if r, g or b could wrap, going close to zero.

	#ifndef HLRAD_FinalLightFace_VL
	        if (g_colour_jitter_hack[0] || g_colour_jitter_hack[1] || g_colour_jitter_hack[2]) 
            {
		        for (i = 0; i < 3; i++) 
                {
		            lb[i] += g_colour_jitter_hack[i] * ((float)rand() / RAND_MAX - 0.5);
		            if (lb[i] < 0.0f)
                    {
			            lb[i] = 0.0f;
                    }
		            else if (lb[i] > 255.0f)
                    {
			            lb[i] = 255.0f;
                    }
		        }
	        }

	        if (g_jitter_hack[0] || g_jitter_hack[1] || g_jitter_hack[2]) 
            {
		        temp_rand = (float)rand() / RAND_MAX - 0.5;
		        for (i = 0; i < 3; i++) 
                {
		            lb[i] += g_jitter_hack[i] * temp_rand;
		            if (lb[i] < 0.0f)
                    {
			            lb[i] = 0.0f;
                    }
		            else if (lb[i] > 255.0f)
                    {
			            lb[i] = 255.0f;
                    }
		        }
	        }
	#endif
#else
            if (g_qgamma != 1.0) {
	            for (i = 0; i < 3; i++) {
	                lb[i] = (float) pow(lb[i] / 256.0f, g_qgamma) * 256.0f;
	            }
	        }
#endif

#ifdef HLRAD_MINLIGHT
			for (i = 0; i < 3; ++i)
				if (lb[i] < g_minlight)
					lb[i] = g_minlight;
#endif
	        // ------------------------------------------------------------------------
#ifdef HLRAD_FinalLightFace_VL
			for (i = 0; i < 3; ++i)
			{
				lbi[i] = (int) floor (lb[i] + 0.5);
				if (lbi[i] < 0) lbi[i] = 0;
			}
			if (k == 0)
			{
				VectorCopy (lbi, final_basiclight[j]);
			}
			else
			{
				VectorSubtract (lbi, final_basiclight[j], lbi);
				if (VectorMinimum (lbi) < -1)
					Warning ("HLRAD_FinalLightFace_VL: bad internal assumption @(%f,%f,%f)", samp->pos[0], samp->pos[1], samp->pos[2]);
			}
	#ifdef HLRAD_WHOME
			if (k == 0)
			{
				if (g_colour_jitter_hack[0] || g_colour_jitter_hack[1] || g_colour_jitter_hack[2]) 
					for (i = 0; i < 3; i++) 
						lbi[i] += g_colour_jitter_hack[i] * ((float)rand() / RAND_MAX - 0.5);
				if (g_jitter_hack[0] || g_jitter_hack[1] || g_jitter_hack[2]) 
				{
					temp_rand = (float)rand() / RAND_MAX - 0.5;
					for (i = 0; i < 3; i++) 
						lbi[i] += g_jitter_hack[i] * temp_rand;
				}
			}
	#endif
			for (i = 0; i < 3; ++i)
			{
				if (lbi[i] < 0) lbi[i] = 0;
				if (lbi[i] > 255) lbi[i] = 255;
			}
            {
                unsigned char* colors = &g_dlightdata[f->lightofs + k * fl->numsamples * 3 + j * 3];

                colors[0] = (unsigned char)lbi[0];
                colors[1] = (unsigned char)lbi[1];
                colors[2] = (unsigned char)lbi[2];
            }
#else
            {
                unsigned char* colors = &g_dlightdata[f->lightofs + k * fl->numsamples * 3 + j * 3];

                colors[0] = (unsigned char)lb[0];
                colors[1] = (unsigned char)lb[1];
                colors[2] = (unsigned char)lb[2];
            }
#endif
        }
    }
#ifdef HLRAD_FinalLightFace_VL
	free (original_basiclight);
	free (final_basiclight);
#endif

#ifndef HLRAD_GatherPatchLight
    if (g_numbounce)
#endif
    {
        FreeTriangulation(trian);
    }
}


#ifdef ZHLT_TEXLIGHT
//LRC
vec3_t    totallight_default = { 0, 0, 0 };

//LRC - utility for getting the right totallight value from a patch
vec3_t* GetTotalLight(patch_t* patch, int style)
{
	int i;
	for (i = 0; i < MAXLIGHTMAPS && patch->totalstyle[i] != 255; i++)
	{
		if (patch->totalstyle[i] == style)
			return &(patch->totallight[i]);
	}
	return &totallight_default;
}

#endif
