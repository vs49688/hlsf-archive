#include "qrad.h"
#ifdef HLRAD_TEXTURE

#ifdef WORDS_BIGENDIAN
#error WORDS_BIGENDIAN
#endif

int g_numtextures;
radtexture_t *g_textures;

typedef struct waddir_s
{
	struct waddir_s *next;
	char path[_MAX_PATH];
} waddir_t;
waddir_t *g_waddirs = NULL;

void AddWadFolder (const char *path)
{
	waddir_t *waddir;
	waddir = (waddir_t *)malloc (sizeof (waddir_t));
	hlassume (waddir != NULL, assume_NoMemory);
	{
		waddir_t **pos;
		for (pos = &g_waddirs; *pos; pos = &(*pos)->next)
			;
		waddir->next = *pos;
		*pos = waddir;
	}
	safe_snprintf (waddir->path, _MAX_PATH, "%s", path);
}

typedef struct
{
	int filepos;
	int disksize;
	int size;
	char type;
	char compression;
	char pad1, pad2;
	char name[16];
} lumpinfo_t;

typedef struct wadfile_s
{
	struct wadfile_s *next;
	char path[_MAX_PATH];
	FILE *file;
	int numlumps;
	lumpinfo_t *lumpinfos;
} wadfile_t;

wadfile_t *g_wadfiles = NULL;
bool g_wadfiles_opened;

static int CDECL lump_sorter_by_name (const void *lump1, const void *lump2)
{
	lumpinfo_t *plump1 = (lumpinfo_t *)lump1;
	lumpinfo_t *plump2 = (lumpinfo_t *)lump2;
	return strcasecmp (plump1->name, plump2->name);
}

void OpenWadFile (const char *name)
{
	int i;
	wadfile_t *wad;
	wad = (wadfile_t *)malloc (sizeof (wadfile_t));
	hlassume (wad != NULL, assume_NoMemory);
	{
		wadfile_t **pos;
		for (pos = &g_wadfiles; *pos; pos = &(*pos)->next)
			;
		wad->next = *pos;
		*pos = wad;
	}
	waddir_t *dir;
	for (dir = g_waddirs; dir; dir = dir->next)
	{
		safe_snprintf (wad->path, _MAX_PATH, "%s\\%s", dir->path, name);
		if (wad->file = fopen (wad->path, "rb"))
		{
			break;
		}
	}
	if (!dir)
	{
		Fatal (assume_COULD_NOT_LOCATE_WAD, "Could not locate wad file %s", name);
		return;
	}
	Log ("Using Wadfile: %s\n", wad->path);
	struct
	{
		char identification[4];
		int numlumps;
		int infotableofs;
	} wadinfo;
	SafeRead (wad->file, &wadinfo, sizeof(wadinfo));
	wadinfo.numlumps  = LittleLong(wadinfo.numlumps);
	wadinfo.infotableofs = LittleLong(wadinfo.infotableofs);
	if (strncmp (wadinfo.identification, "WAD2", 4) && strncmp (wadinfo.identification, "WAD3", 4))
		Error ("%s isn't a Wadfile!", wad->path);
	wad->numlumps = wadinfo.numlumps;
	wad->lumpinfos = (lumpinfo_t *)malloc (wad->numlumps * sizeof (lumpinfo_t));
	hlassume (wad->lumpinfos != NULL, assume_NoMemory);
    if (fseek (wad->file, wadinfo.infotableofs, SEEK_SET))
		Error ("File read failure",wad->path);
	for (i = 0; i < wad->numlumps; i++)
	{
		SafeRead (wad->file, &wad->lumpinfos[i], sizeof (lumpinfo_t));
		if (!TerminatedString(wad->lumpinfos[i].name, 16))
		{
			wad->lumpinfos[i].name[16 - 1] = 0;
			Warning("Unterminated texture name : wad[%s] texture[%d] name[%s]\n", wad->path, i, wad->lumpinfos[i].name);
		}
		wad->lumpinfos[i].filepos = LittleLong(wad->lumpinfos[i].filepos);
		wad->lumpinfos[i].disksize = LittleLong(wad->lumpinfos[i].disksize);
		wad->lumpinfos[i].size = LittleLong(wad->lumpinfos[i].size);
	}
	qsort (wad->lumpinfos, wad->numlumps, sizeof (lumpinfo_t), lump_sorter_by_name);
}

void OpenWadFiles ()
{
	if (!g_wadfiles_opened)
	{
		g_wadfiles_opened = true;
		Log ("Opening wad files from directories:\n");
		if (!g_waddirs)
		{
			Warning ("No wad directories have been set.");
		}
		else
		{
			waddir_t *dir;
			for (dir = g_waddirs; dir; dir = dir->next)
			{
				Log ("  %s\n", dir->path);
			}
		}
		const char *value = ValueForKey (&g_entities[0], "wad");
		char path[MAX_VAL];
		int i, j;
		for (i = 0, j = 0; i < strlen(value) + 1; i++)
		{
			if (value[i] == ';' || value[i] == '\0')
			{
				path[j] = '\0';
				if (path[0])
				{
					char name[MAX_VAL];
					ExtractFile (path, name);
					DefaultExtension (name, ".wad");
					OpenWadFile (name);
				}
				j = 0;
			}
			else
			{
				path[j] = value[i];
				j++;
			}
		}
		CheckFatal ();
	}
}

void CloseWadFiles ()
{
	if (g_wadfiles_opened)
	{
		g_wadfiles_opened = false;
		wadfile_t *wadfile, *next;
		for (wadfile = g_wadfiles; wadfile; wadfile = next)
		{
			next = wadfile->next;
			free (wadfile->lumpinfos);
			fclose (wadfile->file);
			free (wadfile);
		}
		g_wadfiles = NULL;
	}
}

void DefaultTexture (radtexture_t *tex, const char *name)
{
	int i;
	tex->width = 16;
	tex->height = 16;
	strcpy (tex->name, name);
	tex->name[16 - 1] = '\0';
	tex->canvas = (byte *)malloc (tex->width * tex->height);
	hlassume (tex->canvas != NULL, assume_NoMemory);
	for (i = 0; i < 256; i++)
	{
		VectorFill (tex->palette[i], 0x80);
	}
	for (i = 0; i < tex->width * tex->height; i++)
	{
		tex->canvas[i] = 0x00;
	}
}

void LoadTexture (radtexture_t *tex, const miptex_t *mt, int size)
{
	int i, j;
	const miptex_t *header = mt;
	const byte *data = (const byte *)mt;
	tex->width = header->width;
	tex->height = header->height;
	strcpy (tex->name, header->name);
	tex->name[16 - 1] = '\0';
	if (tex->width <= 0 || tex->height <= 0 ||
		tex->width % (2 * 1 << (MIPLEVELS - 1)) != 0 || tex->height % (2 * (1 << (MIPLEVELS - 1))) != 0)
	{
		Error ("Texture '%s' (%dx%d): dimension is not multiple of %d.", tex->name, tex->width, tex->height, 2 * (1 << (MIPLEVELS - 1)));
	}
	int mipsize;
	for (mipsize = 0, i = 0; i < MIPLEVELS; i++)
	{
		if (mt->offsets[i] != sizeof (miptex_t) + mipsize)
		{
			Error ("Texture '%s' (%dx%d): unexpected miptex offset.", tex->name, tex->width, tex->height);
		}
		mipsize += (tex->width >> i) * (tex->height >> i);
	}
	if (size < sizeof (miptex_t) + mipsize + 2 + 256 * 3)
	{
		Error ("Texture '%s' (%dx%d): no enough data.", tex->width, tex->height);
	}
	if (*(unsigned short *)&data[sizeof (miptex_t) + mipsize] != 256)
	{
		Error ("Texture '%s' (%dx%d): palette size is not 256.", tex->name, tex->width, tex->height);
	}
	tex->canvas = (byte *)malloc (tex->width * tex->height);
	hlassume (tex->canvas != NULL, assume_NoMemory);
	for (i = 0; i < tex->height; i++)
	{
		for (j = 0; j < tex->width; j++)
		{
			tex->canvas[i * tex->width + j] = data[sizeof (miptex_t) + i * tex->width + j];
		}
	}
	for (i = 0; i < 256; i++)
	{
		for (j = 0; j < 3; j++)
		{
			tex->palette[i][j] = data[sizeof (miptex_t) + mipsize + 2 + i * 3 + j];
		}
	}
}

void LoadTextureFromWad (radtexture_t *tex, const miptex_t *header)
{
	tex->width = header->width;
	tex->height = header->height;
	strcpy (tex->name, header->name);
	tex->name[16 - 1] = '\0';
	wadfile_t *wad;
	for (wad = g_wadfiles; wad; wad = wad->next)
	{
		lumpinfo_t temp, *found;
		strcpy (temp.name, tex->name);
		found = (lumpinfo_t *)bsearch (&temp, wad->lumpinfos, wad->numlumps, sizeof (lumpinfo_t), lump_sorter_by_name);
		if (found)
		{
			Developer (DEVELOPER_LEVEL_MESSAGE, "Texture '%s' (%dx%d): found in '%s'.\n", header->name, header->width, header->height, wad->path);
			if (found->type != 67 || found->compression != 0)
				continue;
			if (found->disksize < sizeof (miptex_t))
				Error ("Invalid texture data in '%s'.", wad->path);
			miptex_t *mt = (miptex_t *)malloc (found->disksize);
			hlassume (mt != NULL, assume_NoMemory);
			if (fseek (wad->file, found->filepos, SEEK_SET))
				Error ("Invalid texture data in '%s'.", wad->path);
			SafeRead (wad->file, mt, found->disksize);
			if (strcasecmp (mt->name, tex->name) || mt->width != tex->width || mt->height != tex->height)
			{
				free (mt);
				continue;
			}
			LoadTexture (tex, mt, found->disksize);
			free (mt);
			break;
		}
	}
	if (!wad)
	{
		Warning ("Texture '%s' (%dx%d): texture is not found in wad files.", tex->name, tex->width, tex->height);
		DefaultTexture (tex, tex->name);
		return;
	}
}

void LoadTextures ()
{
	if (!g_notextures)
	{
		Log ("Load Textures:\n");
	}
	g_numtextures = g_texdatasize? ((dmiptexlump_t *)g_dtexdata)->nummiptex: 0;
	g_textures = (radtexture_t *)malloc (g_numtextures * sizeof (radtexture_t));
	hlassume (g_textures != NULL, assume_NoMemory);
	int i;
	for (i = 0; i < g_numtextures; i++)
	{
		int offset = ((dmiptexlump_t *)g_dtexdata)->dataofs[i];
		int size = g_texdatasize - offset;
		radtexture_t *tex = &g_textures[i];
		if (g_notextures)
		{
			DefaultTexture (tex, "DEFAULT");
		}
		else if (offset < 0 || size < sizeof (miptex_t))
		{
			Warning ("Invalid texture data in '%s'.", g_source);
			DefaultTexture (tex, "");
		}
		else
		{
			miptex_t *mt = (miptex_t *)&g_dtexdata[offset];
			if (mt->offsets[0])
			{
				Developer (DEVELOPER_LEVEL_MESSAGE, "Texture '%s' (%dx%d): found in '%s'.\n", mt->name, mt->width, mt->height, g_source);
				LoadTexture (tex, mt, size);
			}
			else
			{
				OpenWadFiles ();
				LoadTextureFromWad (tex, mt);
			}
		}
#ifdef HLRAD_REFLECTIVITY
		{
			vec3_t total;
			VectorClear (total);
			for (int j = 0; j < tex->width * tex->height; j++)
			{
				vec3_t reflectivity;
				if (tex->name[0] == '{' && tex->canvas[j] == 0xFF)
				{
					VectorFill (reflectivity, 0.0);
				}
				else
				{
					VectorScale (tex->palette[tex->canvas[j]], 1.0/255.0, reflectivity);
					for (int k = 0; k < 3; k++)
					{
						reflectivity[k] = pow (reflectivity[k], g_texreflectgamma);
					}
					VectorScale (reflectivity, g_texreflectscale, reflectivity);
				}
				VectorAdd (total, reflectivity, total);
			}
			VectorScale (total, 1.0 / (double)(tex->width * tex->height), total);
			VectorCopy (total, tex->reflectivity);
			Developer (DEVELOPER_LEVEL_MESSAGE, "Texture '%s' (%dx%d): reflectivity is (%f,%f,%f).\n",
				tex->name, tex->width, tex->height, tex->reflectivity[0], tex->reflectivity[1], tex->reflectivity[2]);
		}
#endif
	}
	if (!g_notextures)
	{
		Log ("%i textures referenced\n", g_numtextures);
		CloseWadFiles ();
	}
}

#endif
