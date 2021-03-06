#include "bsp5.h"

//  WriteClipNodes_r
//  WriteClipNodes
//  WriteDrawLeaf
//  WriteFace
//  WriteDrawNodes_r
//  FreeDrawNodes_r
//  WriteDrawNodes
//  BeginBSPFile
//  FinishBSPFile

#include <map>

typedef std::map<int,int> PlaneMap;
static PlaneMap gPlaneMap;
static int gNumMappedPlanes;
static dplane_t gMappedPlanes[MAX_MAP_PLANES];
extern bool g_noopt;

// =====================================================================================
//  WritePlane
//  hook for plane optimization
// =====================================================================================
static int WritePlane(int planenum)
{
	planenum = planenum & (~1);

	if(g_noopt)
	{
		return planenum;
	}

	PlaneMap::iterator item = gPlaneMap.find(planenum);
	if(item != gPlaneMap.end())
	{
		return item->second;
	}
	//add plane to BSP
	hlassume(gNumMappedPlanes < MAX_MAP_PLANES, assume_MAX_MAP_PLANES);
	gMappedPlanes[gNumMappedPlanes] = g_dplanes[planenum];
	gPlaneMap.insert(PlaneMap::value_type(planenum,gNumMappedPlanes));

	return gNumMappedPlanes++;
}

// =====================================================================================
//  WriteClipNodes_r
// =====================================================================================
static int      WriteClipNodes_r(node_t* node
#ifdef ZHLT_DETAILBRUSH
								 , const node_t *portalleaf
#endif
								 )
{
    int             i, c;
    dclipnode_t*    cn;
    int             num;

#ifdef ZHLT_DETAILBRUSH
	if (node->isportalleaf)
	{
		if (node->contents == CONTENTS_SOLID)
		{
			free (node);
			return CONTENTS_SOLID;
		}
		else
		{
			portalleaf = node;
		}
	}
	if (node->planenum == -1)
	{
		if (node->iscontentsdetail)
		{
			num = CONTENTS_SOLID;
		}
		else
		{
			num = portalleaf->contents;
		}
		free (node->markfaces);
		free (node);
		return num;
	}
#else
    if (node->planenum == -1)
    {
        num = node->contents;
        free(node->markfaces);
        free(node);
        return num;
    }
#endif

    // emit a clipnode
    hlassume(g_numclipnodes < MAX_MAP_CLIPNODES, assume_MAX_MAP_CLIPNODES);

    c = g_numclipnodes;
    cn = &g_dclipnodes[g_numclipnodes];
    g_numclipnodes++;
    if (node->planenum & 1)
    {
        Error("WriteClipNodes_r: odd planenum");
    }
    cn->planenum = WritePlane(node->planenum);
    for (i = 0; i < 2; i++)
    {
        cn->children[i] = WriteClipNodes_r(node->children[i]
#ifdef ZHLT_DETAILBRUSH
			, portalleaf
#endif
			);
    }

    free(node);
    return c;
}

// =====================================================================================
//  WriteClipNodes
//      Called after the clipping hull is completed.  Generates a disk format
//      representation and frees the original memory.
// =====================================================================================
void            WriteClipNodes(node_t* nodes)
{
    WriteClipNodes_r(nodes
#ifdef ZHLT_DETAILBRUSH
		, NULL
#endif
		);
}

// =====================================================================================
//  WriteDrawLeaf
// =====================================================================================
#ifdef ZHLT_DETAILBRUSH
static int		WriteDrawLeaf (node_t *node, const node_t *portalleaf)
#else
static void     WriteDrawLeaf(const node_t* const node)
#endif
{
    face_t**        fp;
    face_t*         f;
    dleaf_t*        leaf_p;
#ifdef ZHLT_DETAILBRUSH
	int				leafnum = g_numleafs;
#endif

    // emit a leaf
#ifdef ZHLT_MAX_MAP_LEAFS
	hlassume (g_numleafs < MAX_MAP_LEAFS, assume_MAX_MAP_LEAFS);
#endif
    leaf_p = &g_dleafs[g_numleafs];
    g_numleafs++;

#ifdef ZHLT_DETAILBRUSH
	leaf_p->contents = portalleaf->contents;
#else
    leaf_p->contents = node->contents;
#endif

    //
    // write bounding box info
    //
#ifdef ZHLT_DETAILBRUSH
	VectorCopy (portalleaf->mins, leaf_p->mins);
	VectorCopy (portalleaf->maxs, leaf_p->maxs);
#else
    VectorCopy(node->mins, leaf_p->mins);
    VectorCopy(node->maxs, leaf_p->maxs);
#endif

    leaf_p->visofs = -1;                                   // no vis info yet

    //
    // write the marksurfaces
    //
    leaf_p->firstmarksurface = g_nummarksurfaces;

    hlassume(node->markfaces != NULL, assume_EmptySolid);

    for (fp = node->markfaces; *fp; fp++)
    {
        // emit a marksurface
        f = *fp;
        do
        {
#ifdef HLBSP_NULLFACEOUTPUT_FIX
			// fix face 0 being seen everywhere
			if (f->outputnumber == -1)
			{
				f = f->original;
				continue;
			}
#endif
            g_dmarksurfaces[g_nummarksurfaces] = f->outputnumber;
            hlassume(g_nummarksurfaces < MAX_MAP_MARKSURFACES, assume_MAX_MAP_MARKSURFACES);
            g_nummarksurfaces++;
            f = f->original;                               // grab tjunction split faces
        }
        while (f);
    }
    free(node->markfaces);

    leaf_p->nummarksurfaces = g_nummarksurfaces - leaf_p->firstmarksurface;
#ifdef ZHLT_DETAILBRUSH
	return leafnum;
#endif
}

// =====================================================================================
//  WriteFace
// =====================================================================================
static void     WriteFace(face_t* f)
{
    dface_t*        df;
    int             i;
    int             e;

    if (    CheckFaceForHint(f)
        ||  CheckFaceForSkip(f)
#ifdef ZHLT_NULLTEX
        ||  CheckFaceForNull(f)  // AJM
#endif
#ifdef HLCSG_HLBSP_SOLIDHINT
		|| CheckFaceForDiscardable (f)
#endif
#ifdef HLCSG_HLBSP_VOIDTEXINFO
		|| f->texturenum == -1
#endif

// =====================================================================================
//Cpt_Andrew - Env_Sky Check
// =====================================================================================
       ||  CheckFaceForEnv_Sky(f)
// =====================================================================================

       )
    {
#ifdef HLBSP_NULLFACEOUTPUT_FIX
		f->outputnumber = -1;
#endif
        return;
    }

    f->outputnumber = g_numfaces;

    df = &g_dfaces[g_numfaces];
    hlassume(g_numfaces < MAX_MAP_FACES, assume_MAX_MAP_FACES);
    g_numfaces++;

	df->planenum = WritePlane(f->planenum);
	df->side = f->planenum & 1;
    df->firstedge = g_numsurfedges;
    df->numedges = f->numpoints;
    df->texinfo = f->texturenum;
    for (i = 0; i < f->numpoints; i++)
    {
#ifdef ZHLT_DETAILBRUSH
		e = f->outputedges[i];
#else
        e = GetEdge(f->pts[i], f->pts[(i + 1) % f->numpoints], f);
#endif
        hlassume(g_numsurfedges < MAX_MAP_SURFEDGES, assume_MAX_MAP_SURFEDGES);
        g_dsurfedges[g_numsurfedges] = e;
        g_numsurfedges++;
    }
#ifdef ZHLT_DETAILBRUSH
	free (f->outputedges);
	f->outputedges = NULL;
#endif
}

// =====================================================================================
//  WriteDrawNodes_r
// =====================================================================================
#ifdef ZHLT_DETAILBRUSH
static int WriteDrawNodes_r (node_t *node, const node_t *portalleaf)
#else
static void     WriteDrawNodes_r(const node_t* const node)
#endif
{
#ifdef ZHLT_DETAILBRUSH
	if (node->isportalleaf)
	{
		if (node->contents == CONTENTS_SOLID)
		{
			return -1;
		}
		else
		{
			portalleaf = node;
			// Warning: make sure parent data are not freed when writing children.
		}
	}
	if (node->planenum == -1)
	{
		if (node->iscontentsdetail)
		{
			free(node->markfaces);
			return -1;
		}
		else
		{
			int leafnum = WriteDrawLeaf (node, portalleaf);
			return -1 - leafnum;
		}
	}
#endif
    dnode_t*        n;
    int             i;
    face_t*         f;
#ifdef ZHLT_DETAILBRUSH
	int nodenum = g_numnodes;
#endif

    // emit a node
    hlassume(g_numnodes < MAX_MAP_NODES, assume_MAX_MAP_NODES);
    n = &g_dnodes[g_numnodes];
    g_numnodes++;

#ifdef ZHLT_DETAILBRUSH
	if (node->isdetail)
	{
		VectorCopy (portalleaf->mins, n->mins);
		VectorCopy (portalleaf->maxs, n->maxs);
	}
	else
	{
		VectorCopy (node->mins, n->mins);
		VectorCopy (node->maxs, n->maxs);
	}
#else
    VectorCopy(node->mins, n->mins);
    VectorCopy(node->maxs, n->maxs);
#endif

    if (node->planenum & 1)
    {
        Error("WriteDrawNodes_r: odd planenum");
    }
    n->planenum = WritePlane(node->planenum);
    n->firstface = g_numfaces;

    for (f = node->faces; f; f = f->next)
    {
        WriteFace(f);
    }

    n->numfaces = g_numfaces - n->firstface;

    //
    // recursively output the other nodes
    //
    for (i = 0; i < 2; i++)
    {
#ifdef ZHLT_DETAILBRUSH
		n->children[i] = WriteDrawNodes_r (node->children[i], portalleaf);
#else
        if (node->children[i]->planenum == -1)
        {
            if (node->children[i]->contents == CONTENTS_SOLID)
            {
                n->children[i] = -1;
            }
            else
            {
                n->children[i] = -(g_numleafs + 1);
                WriteDrawLeaf(node->children[i]);
            }
        }
        else
        {
            n->children[i] = g_numnodes;
            WriteDrawNodes_r(node->children[i]);
        }
#endif
    }
#ifdef ZHLT_DETAILBRUSH
	return nodenum;
#endif
}

// =====================================================================================
//  FreeDrawNodes_r
// =====================================================================================
static void     FreeDrawNodes_r(node_t* node)
{
    int             i;
    face_t*         f;
    face_t*         next;

    for (i = 0; i < 2; i++)
    {
        if (node->children[i]->planenum != -1)
        {
            FreeDrawNodes_r(node->children[i]);
        }
    }

    //
    // free the faces on the node
    //
    for (f = node->faces; f; f = next)
    {
        next = f->next;
        FreeFace(f);
    }

    free(node);
}

// =====================================================================================
//  WriteDrawNodes
//      Called after a drawing hull is completed
//      Frees all nodes and faces
// =====================================================================================
#ifdef ZHLT_DETAILBRUSH
void OutputEdges_face (face_t *f)
{
	if (CheckFaceForHint(f)
		|| CheckFaceForSkip(f)
#ifdef ZHLT_NULLTEX
        || CheckFaceForNull(f)  // AJM
#endif
#ifdef HLCSG_HLBSP_SOLIDHINT
		|| CheckFaceForDiscardable (f)
#endif
#ifdef HLCSG_HLBSP_VOIDTEXINFO
		|| f->texturenum == -1
#endif
		|| CheckFaceForEnv_Sky(f)//Cpt_Andrew - Env_Sky Check
		)
	{
		return;
	}
	f->outputedges = (int *)malloc (f->numpoints * sizeof (int));
	hlassume (f->outputedges != NULL, assume_NoMemory);
	int i;
	for (i = 0; i < f->numpoints; i++)
	{
		int e = GetEdge (f->pts[i], f->pts[(i + 1) % f->numpoints], f);
		f->outputedges[i] = e;
	}
}
int OutputEdges_r (node_t *node, int detaillevel)
{
	int next = -1;
	if (node->planenum == -1)
	{
		return next;
	}
	face_t *f;
	for (f = node->faces; f; f = f->next)
	{
		if (f->detaillevel > detaillevel)
		{
			if (next == -1? true: f->detaillevel < next)
			{
				next = f->detaillevel;
			}
		}
		if (f->detaillevel == detaillevel)
		{
			OutputEdges_face (f);
		}
	}
	int i;
	for (i = 0; i < 2; i++)
	{
		int r = OutputEdges_r (node->children[i], detaillevel);
		if (r == -1? false: next == -1? true: r < next)
		{
			next = r;
		}
	}
	return next;
}
#endif
void            WriteDrawNodes(node_t* headnode)
{
#ifdef ZHLT_DETAILBRUSH
	// higher detail level should not compete for edge pairing with lower detail level.
	int detaillevel, nextdetaillevel;
	for (detaillevel = 0; detaillevel != -1; detaillevel = nextdetaillevel)
	{
		nextdetaillevel = OutputEdges_r (headnode, detaillevel);
	}
	WriteDrawNodes_r (headnode, NULL);
#else
    if (headnode->contents < 0)
    {
        WriteDrawLeaf(headnode);
    }
    else
    {
        WriteDrawNodes_r(headnode);
        FreeDrawNodes_r(headnode);
    }
#endif
}


// =====================================================================================
//  BeginBSPFile
// =====================================================================================
void            BeginBSPFile()
{
    // these values may actually be initialized
    // if the file existed when loaded, so clear them explicitly
	gNumMappedPlanes = 0;
	gPlaneMap.clear();
    g_nummodels = 0;
    g_numfaces = 0;
    g_numnodes = 0;
    g_numclipnodes = 0;
    g_numvertexes = 0;
    g_nummarksurfaces = 0;
    g_numsurfedges = 0;

    // edge 0 is not used, because 0 can't be negated
    g_numedges = 1;

    // leaf 0 is common solid with no faces
    g_numleafs = 1;
    g_dleafs[0].contents = CONTENTS_SOLID;
}

// =====================================================================================
//  FinishBSPFile
// =====================================================================================
void            FinishBSPFile()
{
    Verbose("--- FinishBSPFile ---\n");

#ifdef ZHLT_MAX_MAP_LEAFS
	if (g_dmodels[0].visleafs > MAX_MAP_LEAFS_ENGINE)
	{
		Warning ("Number of visleafs(%d) exceeded MAX_MAP_LEAFS(%d)\nIf you encounter problems when running your map, consider this the most likely cause.\n", g_dmodels[0].visleafs, MAX_MAP_LEAFS_ENGINE);
	}
#endif
	if(!g_noopt)
	{
#ifdef HLBSP_REDUCETEXTURE
		{
			bool *Used = (bool *)calloc (g_numtexinfo, sizeof(bool));
			int Num = 0;
			int *Map = (int *)malloc (g_numtexinfo * sizeof(int));
			int i;
			hlassume (Used != NULL && Map != NULL, assume_NoMemory);
			for (i = 0; i < g_numfaces; i++)
			{
				dface_t *f = &g_dfaces[i];
				if (f->texinfo < 0 || f->texinfo >= g_numtexinfo)
				{
					Warning ("Bad texinfo number %d.\n", f->texinfo);
					goto skipReduceTexinfo;
				}
				Used[f->texinfo] = true;
			}
			for (i = 0; i < g_numtexinfo; i++)
			{
				if (Used[i])
				{
					g_texinfo[Num] = g_texinfo[i];
					Map[i] = Num;
					Num++;
				}
				else
				{
					Map[i] = -1;
				}
			}
			for (i = 0; i < g_numfaces; i++)
			{
				dface_t *f = &g_dfaces[i];
				f->texinfo = Map[f->texinfo];
			}
			Log ("Reduced %d texinfos to %d\n", g_numtexinfo, Num);
			g_numtexinfo = Num;
			skipReduceTexinfo:;
			free (Used);
			free (Map);
		}
		{
			dmiptexlump_t *l = (dmiptexlump_t *)g_dtexdata;
			int &g_nummiptex = l->nummiptex;
			bool *Used = (bool *)calloc (g_nummiptex, sizeof(bool));
			int Num = 0, Size = 0;
			int *Map = (int *)malloc (g_nummiptex * sizeof(int));
			int i;
			hlassume (Used != NULL && Map != NULL, assume_NoMemory);
			int *lumpsizes = (int *)malloc (g_nummiptex * sizeof (int));
			const int newdatasizemax = g_texdatasize - ((byte *)&l->dataofs[g_nummiptex] - (byte *)l);
			byte *newdata = (byte *)malloc (newdatasizemax);
			int newdatasize = 0;
			hlassume (lumpsizes != NULL && newdata != NULL, assume_NoMemory);
			int total = 0;
			for (i = 0; i < g_nummiptex; i++)
			{
				if (l->dataofs[i] == -1)
				{
					lumpsizes[i] = -1;
					continue;
				}
				lumpsizes[i] = g_texdatasize - l->dataofs[i];
				for (int j = 0; j < g_nummiptex; j++)
				{
					int lumpsize = l->dataofs[j] - l->dataofs[i];
					if (l->dataofs[j] == -1 || lumpsize < 0 || lumpsize == 0 && j <= i)
						continue;
					if (lumpsize < lumpsizes[i])
						lumpsizes[i] = lumpsize;
				}
				total += lumpsizes[i];
			}
			if (total != newdatasizemax)
			{
				Warning ("Bad texdata structure.\n");
				goto skipReduceTexdata;
			}
			for (i = 0; i < g_numtexinfo; i++)
			{
				texinfo_t *t = &g_texinfo[i];
				if (t->miptex < 0 || t->miptex >= g_nummiptex)
				{
					Warning ("Bad miptex number %d.\n", t->miptex);
					goto skipReduceTexdata;
				}
				Used[t->miptex] = true;
			}
			for (i = 0; i < g_nummiptex; i++)
			{
				const int MAXWADNAME = 16;
				char name[MAXWADNAME];
				int j, k;
				if (l->dataofs[i] < 0)
					continue;
				if (Used[i] == true)
				{
					miptex_t *m = (miptex_t *)((byte *)l + l->dataofs[i]);
					if (m->name[0] != '+' && m->name[0] != '-')
						continue;
					safe_strncpy (name, m->name, MAXWADNAME);
					if (name[1] == '\0')
						continue;
					for (j = 0; j < 20; j++)
					{
						if (j < 10)
							name[1] = '0' + j;
						else
							name[1] = 'A' + j - 10;
						for (k = 0; k < g_nummiptex; k++)
						{
							if (l->dataofs[k] < 0)
								continue;
							miptex_t *m2 = (miptex_t *)((byte *)l + l->dataofs[k]);
							if (!stricmp (name, m2->name))
								Used[k] = true;
						}
					}
				}
			}
			for (i = 0; i < g_nummiptex; i++)
			{
				if (Used[i])
				{
					Map[i] = Num;
					Num++;
				}
				else
				{
					Map[i] = -1;
				}
			}
			for (i = 0; i < g_numtexinfo; i++)
			{
				texinfo_t *t = &g_texinfo[i];
				t->miptex = Map[t->miptex];
			}
			Size += (byte *)&l->dataofs[Num] - (byte *)l;
			for (i = 0; i < g_nummiptex; i++)
			{
				if (Used[i])
				{
					if (lumpsizes[i] == -1)
					{
						l->dataofs[Map[i]] = -1;
					}
					else
					{
						memcpy ((byte *)newdata + newdatasize, (byte *)l + l->dataofs[i], lumpsizes[i]);
						l->dataofs[Map[i]] = Size;
						newdatasize += lumpsizes[i];
						Size += lumpsizes[i];
					}
				}
			}
			memcpy (&l->dataofs[Num], newdata, newdatasize);
			Log ("Reduced %d texdatas to %d (%d bytes to %d)\n", g_nummiptex, Num, g_texdatasize, Size);
			g_nummiptex = Num;
			g_texdatasize = Size;
			skipReduceTexdata:;
			free (lumpsizes);
			free (newdata);
			free (Used);
			free (Map);
		}
		Log ("Reduced %d planes to %d\n", g_numplanes, gNumMappedPlanes);
#endif
		for(int counter = 0; counter < gNumMappedPlanes; counter++)
		{
			g_dplanes[counter] = gMappedPlanes[counter];
		}
		g_numplanes = gNumMappedPlanes;
	}

	if (g_chart)
    {
        PrintBSPFileSizes();
    }

#ifdef HLCSG_HLBSP_DOUBLEPLANE
#undef dplane_t
#undef g_dplanes
	for (int i = 0; i < g_numplanes; i++)
	{
		plane_t *mp = &g_mapplanes[i];
		dplane_t *dp = &g_dplanes[i];
		VectorCopy (mp->normal, dp->normal);
		dp->dist = mp->dist;
		dp->type = mp->type;
	}
#define dplane_t plane_t
#define g_dplanes g_mapplanes
#endif
    WriteBSPFile(g_bspfilename);
}
