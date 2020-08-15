/***
*
*	Copyright (c) 1998, Valve LLC. All rights reserved.
*	
*	This product contains software technology licensed from Id 
*	Software, Inc. ("Id Technology").  Id Technology (c) 1996 Id Software, Inc. 
*	All Rights Reserved.
*
****/

// csg4.c

#include "ripent.h"

typedef enum
{
    hl_undefined = -1,
    hl_export = 0,
    hl_import = 1
}
hl_types;

static hl_types g_mode = hl_undefined;

// g_parse: command line switch (-parse).
// Added by: Ryan Gregg aka Nem
bool g_parse = DEFAULT_PARSE;

bool g_chart = DEFAULT_CHART;

bool g_info = DEFAULT_INFO;

#ifdef RIPENT_PAUSE
bool g_pause = false;
#endif


// ScanForToken()
// Added by: Ryan Gregg aka Nem
// 
// Scans entity data starting  at iIndex for cToken.  Every time a \n char
// is encountered iLine is incremented.  If iToken is not null, the index
// cToken was found at is inserted into it. 
bool ScanForToken(char cToken, int &iIndex, int &iLine, bool bIgnoreWhiteSpace, bool bIgnoreOthers, int *iToken = 0)
{
	for(; iIndex < g_entdatasize; iIndex++)
	{
		// If we found a null char, consider it end of data.
		if(g_dentdata[iIndex] == '\0')
		{
			iIndex = g_entdatasize;
			return false;
		}

		// Count lines (for error message).
		if(g_dentdata[iIndex] == '\n')
		{
			iLine++;
		}

		// Ignore white space, if we are ignoring it.
		if(!bIgnoreWhiteSpace && isspace(g_dentdata[iIndex]))
		{
			continue;
		}

		if(g_dentdata[iIndex] != cToken)
		{
			if(bIgnoreOthers)
				continue;
			else
				return false;
		}

		// Advance the index past the token.
		iIndex++;

		// Return the index of the token if requested.
		if(iToken != 0)
		{
			*iToken = iIndex - 1;
		}

		return true;
	}

	// End of data.
	return false;
}

#include <list>
typedef std::list<char *> CEntityPairList;
typedef std::list<CEntityPairList *> CEntityList;

// ParseEntityData()
// Added by: Ryan Gregg aka Nem
// 
// Pareses and reformats entity data stripping all non essential
// formatting  and using the formatting  options passed through this
// function.  The length is specified because in some cases (i.e. the
// terminator) a null char is desired to be printed.
void ParseEntityData(const char *cTab, int iTabLength, const char *cNewLine, int iNewLineLength, const char *cTerminator, int iTerminatorLength)
{
	CEntityList EntityList;		// Parsed entities.
	CEntityList::iterator i;	// eliminate redefine error. --vluzacn

	int iIndex = 0;				// Current char in g_dentdata.
	int iLine = 0;				// Current line in g_dentdata.

	char cError[256] = "";

	try
	{
		//
		// Parse entity data.
		//

		Log("\nParsing entity data.\n");

		while(true)
		{
			// Parse the start of an entity.
			if(!ScanForToken('{', iIndex, iLine, false, false))
			{
				if(iIndex == g_entdatasize)
				{
					// We read all the entities.
					break;
				}
				else
				{
					sprintf_s(cError, "expected token %s on line %d.", "{", iLine);
					throw cError;
				}
			}

			CEntityPairList *EntityPairList = new CEntityPairList();

			// Parse the rest of the entity.
			while(true)
			{
				// Parse the key and value.
				for(int j = 0; j < 2; j++)
				{
					int iStart;
					// Parse the start of a string.
					if(!ScanForToken('\"', iIndex, iLine, false, false, &iStart))
					{
						sprintf_s(cError, "expected token %s on line %d.", "\"", iLine);
						throw cError;
					}

					int iEnd;
					// Parse the end of a string.
					if(!ScanForToken('\"', iIndex, iLine, true, true, &iEnd))
					{
						sprintf_s(cError, "expected token %s on line %d.", "\"", iLine);
						throw cError;
					}

					// Extract the string.
					int iLength = iEnd - iStart - 1;
					char *cString = new char[iLength + 1];
					memcpy(cString, &g_dentdata[iStart + 1], iLength);
					cString[iLength] = '\0';

					// Save it.
					EntityPairList->push_back(cString);
				}

				// Parse the end of an entity.
				if(!ScanForToken('}', iIndex, iLine, false, false))
				{
					if(g_dentdata[iIndex] == '\"')
					{
						// We arn't done the entity yet.
						continue;
					}
					else
					{
						sprintf_s(cError, "expected token %s on line %d.", "}", iLine);
						throw cError;
					}
				}

				// We read the entity.
				EntityList.push_back(EntityPairList);
				break;
			}
		}

		Log("%d entities parsed.\n", EntityList.size());

		//
		// Calculate new data length.
		//

		int iNewLength = 0;

		for(/*CEntityList::iterator*/ i = EntityList.begin(); i != EntityList.end(); ++i) //--vluzacn
		{
			// Opening brace.
			iNewLength += 1;

			// New line.
			iNewLength += iNewLineLength;

			CEntityPairList *EntityPairList = *i;

			for(CEntityPairList::iterator j = EntityPairList->begin(); j != EntityPairList->end(); ++j)
			{
				// Tab.
				iNewLength += iTabLength;

				// String.
				iNewLength += 1;
				iNewLength += (int)strlen(*j);
				iNewLength += 1;

				// String seperator.
				iNewLength += 1;

				++j;

				// String.
				iNewLength += 1;
				iNewLength += (int)strlen(*j);
				iNewLength += 1;

				// New line.
				iNewLength += iNewLineLength;
			}

			// Closing brace.
			iNewLength += 1;

			// New line.
			iNewLength += iNewLineLength;
		}

		// Terminator.
		iNewLength += iTerminatorLength;

		//
		// Check our parsed data.
		//

		assume(iNewLength != 0, "No entity data.");
		assume(iNewLength < sizeof(g_dentdata), "Entity data size exceedes dentdata limit.");

		//
		// Clear current data.
		//

		g_entdatasize = 0;

		//
		// Fill new data.
		//

		Log("Formating entity data.\n\n");

		for(/*CEntityList::iterator*/ i = EntityList.begin(); i != EntityList.end(); ++i) //--vluzacn
		{
			// Opening brace.
			g_dentdata[g_entdatasize] = '{';
			g_entdatasize += 1;

			// New line.
			memcpy(&g_dentdata[g_entdatasize], cNewLine, iNewLineLength);
			g_entdatasize += iNewLineLength;

			CEntityPairList *EntityPairList = *i;

			for(CEntityPairList::iterator j = EntityPairList->begin(); j != EntityPairList->end(); ++j)
			{
				// Tab.
				memcpy(&g_dentdata[g_entdatasize], cTab, iTabLength);
				g_entdatasize += iTabLength;

				// String.
				g_dentdata[g_entdatasize] = '\"';
				g_entdatasize += 1;
				memcpy(&g_dentdata[g_entdatasize], *j, strlen(*j));
				g_entdatasize += (int)strlen(*j);
				g_dentdata[g_entdatasize] = '\"';
				g_entdatasize += 1;

				// String seperator.
				g_dentdata[g_entdatasize] = ' ';
				g_entdatasize += 1;

				++j;

				// String.
				g_dentdata[g_entdatasize] = '\"';
				g_entdatasize += 1;
				memcpy(&g_dentdata[g_entdatasize], *j, strlen(*j));
				g_entdatasize += (int)strlen(*j);
				g_dentdata[g_entdatasize] = '\"';
				g_entdatasize += 1;

				// New line.
				memcpy(&g_dentdata[g_entdatasize], cNewLine, iNewLineLength);
				g_entdatasize += iNewLineLength;
			}

			// Closing brace.
			g_dentdata[g_entdatasize] = '}';
			g_entdatasize += 1;

			// New line.
			memcpy(&g_dentdata[g_entdatasize], cNewLine, iNewLineLength);
			g_entdatasize += iNewLineLength;
		}

		// Terminator.
		memcpy(&g_dentdata[g_entdatasize], cTerminator, iTerminatorLength);
		g_entdatasize += iTerminatorLength;

		//
		// Delete entity data.
		//

		for(/*CEntityList::iterator*/ i = EntityList.begin(); i != EntityList.end(); ++i) //--vluzacn
		{
			CEntityPairList *EntityPairList = *i;

			for(CEntityPairList::iterator j = EntityPairList->begin(); j != EntityPairList->end(); ++j)
			{
				delete []*j;
			}

			delete EntityPairList;
		}

		//return true;
	}
	catch(...)
	{
		//
		// Delete entity data.
		//

		for(CEntityList::iterator i = EntityList.begin(); i != EntityList.end(); ++i)
		{
			CEntityPairList *EntityPairList = *i;

			for(CEntityPairList::iterator j = EntityPairList->begin(); j != EntityPairList->end(); ++j)
			{
				delete []*j;
			}

			delete EntityPairList;
		}

		// If we threw the error cError wont be null, this is
		// a message, print it.
		if(*cError != '\0')
		{
			Error(cError);
		}
		Error("unknowen exception.");

		//return false;
	}
}

static void     ReadBSP(const char* const name)
{
    char            filename[_MAX_PATH];

    safe_strncpy(filename, name, _MAX_PATH);
    StripExtension(filename);
    DefaultExtension(filename, ".bsp");

    LoadBSPFile(name);
}

static void     WriteBSP(const char* const name)
{
    char            filename[_MAX_PATH];

    safe_strncpy(filename, name, _MAX_PATH);
    StripExtension(filename);
    DefaultExtension(filename, ".bsp");

    WriteBSPFile(filename);
}
#ifdef RIPENT_TEXTURE
#ifdef WORDS_BIGENDIAN
#error
#endif
typedef struct
{
    char            identification[4];                     // should be WAD2/WAD3
    int             numlumps;
    int             infotableofs;
} wadinfo_t;
typedef struct
{
	int			filepos;
	int			disksize;
	int			size;					// uncompressed
	char		type;
	char		compression;
	char		pad1, pad2;
	char		name[16];				// must be null terminated
} lumpinfo_t;
/*int TextureSize(const miptex_t *tex)
{
	int size = 0;
	int w, h;
	size += sizeof(miptex_t);
	w = tex->width, h = tex->height;
	for (int imip = 0; imip < MIPLEVELS; ++imip, w/=2, h/=2)
	{
		size += w * h;
	}
	size += 256 * 3 + 4;
	return size;
}
void WriteEmptyTexture(FILE *outwad, const miptex_t *tex)
{
	miptex_t outtex;
	int start, end;
	int w, h;
	memcpy (&outtex, tex, sizeof(miptex_t));
	start = ftell (outwad);
	fseek (outwad, sizeof(miptex_t), SEEK_CUR);
	w = tex->width, h = tex->height;
	for (int imip = 0; imip < MIPLEVELS; ++imip, w/=2, h/=2)
	{
		void *tmp = calloc (w * h, 1);
		outtex.offsets[imip] = ftell (outwad) - start;
		SafeWrite (outwad, tmp, w * h);
		free (tmp);
	}
	short s = 256;
	SafeWrite (outwad, &s, sizeof(short));
	void *tmp = calloc (256 * 3 + 2, 1); // assume width and height are multiples of 16
	SafeWrite (outwad, tmp, 256 * 3 + 2);
	free (tmp);
	end = ftell (outwad);
	fseek (outwad, start, SEEK_SET);
	SafeWrite (outwad, &outtex, sizeof (miptex_t));
	fseek (outwad, end, SEEK_SET);
}*/
static void		WriteTextures(const char* const name)
{
	char wadfilename[_MAX_PATH];
	FILE *wadfile;
    safe_strncpy(wadfilename, name, _MAX_PATH);
    StripExtension(wadfilename);
    DefaultExtension(wadfilename, ".wad");
    _unlink(wadfilename);
	wadfile = SafeOpenWrite (wadfilename);
	Log("\nWriting %s.\n", wadfilename);

    char texfilename[_MAX_PATH];
	FILE *texfile;
    safe_strncpy(texfilename, name, _MAX_PATH);
    StripExtension(texfilename);
    DefaultExtension(texfilename, ".tex");
    _unlink(texfilename);
	if (!g_parse)
	{
		int dataofs = (int)&((dmiptexlump_t*)NULL)->dataofs[((dmiptexlump_t*)g_dtexdata)->nummiptex];
		int wadofs = sizeof(wadinfo_t);

		wadinfo_t header;
		header.identification[0] = 'W';
		header.identification[1] = 'A';
		header.identification[2] = 'D';
		header.identification[3] = '3';
		header.numlumps = ((dmiptexlump_t*)g_dtexdata)->nummiptex;
		header.infotableofs = g_texdatasize - dataofs + wadofs;
		SafeWrite (wadfile, &header, wadofs);

		SafeWrite (wadfile, (byte *)g_dtexdata + dataofs, g_texdatasize - dataofs);

		lumpinfo_t *info;
		info = (lumpinfo_t *)malloc (((dmiptexlump_t*)g_dtexdata)->nummiptex * sizeof (lumpinfo_t));
		hlassume (info != NULL, assume_NoMemory);
		memset (info, 0, header.numlumps * sizeof(lumpinfo_t));

		for (int i = 0; i < header.numlumps; i++)
		{
			int ofs = ((dmiptexlump_t*)g_dtexdata)->dataofs[i];
			int size = 0;
			if (ofs >= 0)
			{
				size = g_texdatasize - ofs;
				for (int j = 0; j < ((dmiptexlump_t*)g_dtexdata)->nummiptex; ++j)
					if (ofs < ((dmiptexlump_t*)g_dtexdata)->dataofs[j] &&
						ofs + size > ((dmiptexlump_t*)g_dtexdata)->dataofs[j])
						size = ((dmiptexlump_t*)g_dtexdata)->dataofs[j] - ofs;
			}
			info[i].filepos = ofs - dataofs + wadofs;
			info[i].disksize = size;
			info[i].size = size;
			info[i].type = (ofs >= 0 && ((miptex_t*)(g_dtexdata+ofs))->offsets[0] > 0)? 67: 0;
			info[i].compression = 0;
			strcpy (info[i].name, ofs >= 0? ((miptex_t*)(g_dtexdata+ofs))->name: "\rTEXTUREMISSING");
		}
		SafeWrite (wadfile, info, header.numlumps * sizeof(lumpinfo_t));
		free (info);
	}
	else
	{
		texfile = SafeOpenWrite (texfilename);
		Log("\nWriting %s.\n", texfilename);

		wadinfo_t header;
		header.identification[0] = 'W';
		header.identification[1] = 'A';
		header.identification[2] = 'D';
		header.identification[3] = '3';
		header.numlumps = 0;
		
		lumpinfo_t *info;
		info = (lumpinfo_t *)malloc (((dmiptexlump_t*)g_dtexdata)->nummiptex * sizeof (lumpinfo_t)); // might be more than needed
		hlassume (info != NULL, assume_NoMemory);

		fprintf (texfile, "%d\r\n", ((dmiptexlump_t*)g_dtexdata)->nummiptex);
		fseek (wadfile, sizeof(wadinfo_t), SEEK_SET);

		for (int itex = 0; itex < ((dmiptexlump_t*)g_dtexdata)->nummiptex; ++itex)
		{
			int ofs = ((dmiptexlump_t*)g_dtexdata)->dataofs[itex];
			miptex_t *tex = (miptex_t*)(g_dtexdata+ofs);
			if (ofs < 0)
			{
				fprintf (texfile, "[-1]\r\n");
			}
			else
			{
				int size = g_texdatasize - ofs;
				for (int j = 0; j < ((dmiptexlump_t*)g_dtexdata)->nummiptex; ++j)
					if (ofs < ((dmiptexlump_t*)g_dtexdata)->dataofs[j] &&
						ofs + size > ((dmiptexlump_t*)g_dtexdata)->dataofs[j])
						size = ((dmiptexlump_t*)g_dtexdata)->dataofs[j] - ofs;
				bool included = false;
				if (tex->offsets[0] > 0)
					included = true;
				if (included)
				{
					memset (&info[header.numlumps], 0, sizeof (lumpinfo_t));
					info[header.numlumps].filepos = ftell (wadfile);
					SafeWrite (wadfile, tex, size);
					info[header.numlumps].disksize = ftell (wadfile) - info[header.numlumps].filepos;
					info[header.numlumps].size = info[header.numlumps].disksize;
					info[header.numlumps].type = 67;
					info[header.numlumps].compression = 0;
					strcpy (info[header.numlumps].name, tex->name);
					header.numlumps++;
				}
				fprintf (texfile, "[%d]", strlen(tex->name));
				SafeWrite (texfile, tex->name, strlen(tex->name));
				fprintf (texfile, " %d %d\r\n", tex->width, tex->height);
			}
		}
		header.infotableofs = ftell (wadfile);
		SafeWrite (wadfile, info, header.numlumps * sizeof(lumpinfo_t));
		fseek (wadfile, 0, SEEK_SET);
		SafeWrite (wadfile, &header, sizeof(wadinfo_t));

		fclose (texfile);
		free (info);
	}
	fclose (wadfile);
}
inline void skipspace (FILE *f) { fscanf (f, "%*[ \t\r\n]s");}
inline void skipline (FILE *f) { fscanf (f, "%*[^\r\n]s");}
static void		ReadTextures(const char *name)
{
	char wadfilename[_MAX_PATH];
	FILE *wadfile;
    safe_strncpy(wadfilename, name, _MAX_PATH);
    StripExtension(wadfilename);
    DefaultExtension(wadfilename, ".wad");
	wadfile = SafeOpenRead (wadfilename);
	Log("\nReading %s.\n", wadfilename);

    char texfilename[_MAX_PATH];
	FILE *texfile;
    safe_strncpy(texfilename, name, _MAX_PATH);
    StripExtension(texfilename);
    DefaultExtension(texfilename, ".tex");
	if (!g_parse)
	{
		wadinfo_t header;
		int wadofs = sizeof(wadinfo_t);
		SafeRead (wadfile, &header, wadofs);
		((dmiptexlump_t*)g_dtexdata)->nummiptex = header.numlumps;
		int dataofs = (int)&((dmiptexlump_t*)NULL)->dataofs[((dmiptexlump_t*)g_dtexdata)->nummiptex];
		g_texdatasize = header.infotableofs - wadofs + dataofs;

		SafeRead (wadfile, (byte *)g_dtexdata + dataofs, g_texdatasize - dataofs);
		
		lumpinfo_t *info;
		info = (lumpinfo_t *)malloc (header.numlumps * sizeof (lumpinfo_t));
		hlassume (info != NULL, assume_NoMemory);
		SafeRead (wadfile, info, header.numlumps * sizeof(lumpinfo_t));

		for (int i = 0; i < header.numlumps; i++)
		{
			((dmiptexlump_t*)g_dtexdata)->dataofs[i] = info[i].filepos - wadofs + dataofs;
		}

		free (info);
	}
	else
	{
		texfile = SafeOpenRead (texfilename);
		Log("\nReading %s.\n", texfilename);

		wadinfo_t header;
		SafeRead (wadfile, &header, sizeof(wadinfo_t));
		fseek (wadfile, header.infotableofs, SEEK_SET);

		lumpinfo_t *info;
		info = (lumpinfo_t *)malloc (header.numlumps * sizeof (lumpinfo_t));
		hlassume (info != NULL, assume_NoMemory);
		SafeRead (wadfile, info, header.numlumps * sizeof(lumpinfo_t));

		int nummiptex = 0;
		if (skipspace (texfile), fscanf (texfile, "%d", &nummiptex) != 1)
			Error ("File read failure");
		((dmiptexlump_t*)g_dtexdata)->nummiptex = nummiptex;
		g_texdatasize = (byte *)(&((dmiptexlump_t*)g_dtexdata)->dataofs[nummiptex]) - g_dtexdata;

		for (int itex = 0; itex < nummiptex; ++itex)
		{
			int len;
			if (skipspace (texfile), fscanf (texfile, "[%d]", &len) != 1)
				Error ("File read failure");
			if (len < 0)
			{
				((dmiptexlump_t*)g_dtexdata)->dataofs[itex] = -1;
			}
			else
			{
				char name[16];
				if (len > 15)
					Error ("Texture name is too long");
				memset (name, '\0', 16);
				SafeRead (texfile, name, len);
				((dmiptexlump_t*)g_dtexdata)->dataofs[itex] = g_texdatasize;
				miptex_t *tex = (miptex_t*)(g_dtexdata + g_texdatasize);
				int j;
				for (j = 0; j < header.numlumps; ++j)
					if (!stricmp (name, info[j].name))
						break;
				if (j == header.numlumps)
				{
					int w, h;
					if (skipspace (texfile), fscanf (texfile, "%d", &w) != 1)
						Error ("File read failure");
					if (skipspace (texfile), fscanf (texfile, "%d", &h) != 1)
						Error ("File read failure");
					g_texdatasize += sizeof(miptex_t);
					hlassume(g_texdatasize < g_max_map_miptex, assume_MAX_MAP_MIPTEX);
					memset (tex, 0, sizeof(miptex_t));
					strcpy (tex->name, name);
					tex->width = w;
					tex->height = h;
					for (int k = 0; k < MIPLEVELS; k++)
						tex->offsets[k] = 0;
				}
				else
				{
					fseek (wadfile, info[j].filepos, SEEK_SET);
					g_texdatasize += info[j].disksize;
					hlassume(g_texdatasize < g_max_map_miptex, assume_MAX_MAP_MIPTEX);
					SafeRead (wadfile, tex, info[j].disksize);
					strcpy (tex->name, name);
				}
			}
			skipline (texfile);
		}

		fclose (texfile);
		free (info);
	}
	fclose (wadfile);
}
#endif /*RIPENT_TEXTURE*/

static void     WriteEntities(const char* const name)
{
    char filename[_MAX_PATH];

    safe_strncpy(filename, name, _MAX_PATH);
    StripExtension(filename);
    DefaultExtension(filename, ".ent");
    _unlink(filename);

    {
		if(g_parse)  // Added by Nem.
		{
			ParseEntityData("  ", 2, "\r\n", 2, "", 0);
		}

        FILE *f = SafeOpenWrite(filename);
		Log("\nWriting %s.\n", filename);  // Added by Nem.
        SafeWrite(f, g_dentdata, g_entdatasize);
        fclose(f);
    }
}

static void     ReadEntities(const char* const name)
{
    char filename[_MAX_PATH];

    safe_strncpy(filename, name, _MAX_PATH);
    StripExtension(filename);
    DefaultExtension(filename, ".ent");

    {
        FILE *f = SafeOpenRead(filename);
		Log("\nReading %s.\n", filename);  // Added by Nem.

        g_entdatasize = q_filelength(f);

		assume(g_entdatasize != 0, "No entity data.");
        assume(g_entdatasize < sizeof(g_dentdata), "Entity data size exceedes dentdata limit.");

        SafeRead(f, g_dentdata, g_entdatasize);

        fclose(f);

        if (g_dentdata[g_entdatasize-1] != 0)
        {
//            Log("g_dentdata[g_entdatasize-1] = %d\n", g_dentdata[g_entdatasize-1]);

			if(g_parse)  // Added by Nem.
			{
				ParseEntityData("", 0, "\n", 1, "\0", 1);
			}
			else
			{
				if(g_dentdata[g_entdatasize - 1] != '\0')
				{
					g_dentdata[g_entdatasize] = '\0';
					g_entdatasize++;
				}
			}
        }
    }
}

//======================================================================

static void     Usage(void)
{
    //Log("%s " ZHLT_VERSIONSTRING "\n" MODIFICATIONS_STRING "\n", g_Program);
    //Log("  Usage: ripent [-import|-export] [-texdata n] bspname\n");

	// Modified to behave like other tools by Nem.

	Banner();
	Log("\n-= %s Options =-\n\n", g_Program);

#ifdef ZHLT_CONSOLE
	Log("    -console 0      : turn console off\n");
#endif
#ifdef ZHLT_LANGFILE
	Log("    -lang file      : localization file\n");
#endif
#ifdef RIPENT_TEXTURE
	Log("    -export         : Export entity data and texture data\n");
	Log("    -import         : Import entity data and texture data\n\n");
#else
	Log("    -export         : Export entity data\n");
	Log("    -import         : Import entity data\n\n");
#endif

	Log("    -parse          : Parse and format entity data\n\n");

    Log("    -texdata #      : Alter maximum texture memory limit (in kb)\n");
    Log("    -lightdata #    : Alter maximum lighting memory limit (in kb)\n");
	Log("    -chart          : Display bsp statitics\n");
	Log("    -noinfo         : Do not show tool configuration information\n\n");
#ifdef RIPENT_PAUSE
	Log("    -pause          : Pause before exit\n\n");
#endif

	Log("    mapfile         : The mapfile to process\n\n");

    exit(1);
}

#ifdef RIPENT_PAUSE
void pause ()
{
	if (g_pause)
	{
		Log("\nPress any key to continue\n");
		getch ();
	}
}
#endif

// =====================================================================================
//  Settings
// =====================================================================================
static void     Settings()
{
    char*           tmp;

    if (!g_info)
    {
        return; 
    }

    Log("\n-= Current %s Settings =-\n", g_Program);
    Log("Name               |  Setting  |  Default\n" "-------------------|-----------|-------------------------\n");

    // ZHLT Common Settings
    Log("chart               [ %7s ] [ %7s ]\n", g_chart ? "on" : "off", DEFAULT_CHART ? "on" : "off");
    Log("max texture memory  [ %7d ] [ %7d ]\n", g_max_map_miptex, DEFAULT_MAX_MAP_MIPTEX);
	Log("max lighting memory [ %7d ] [ %7d ]\n", g_max_map_lightdata, DEFAULT_MAX_MAP_LIGHTDATA);

    switch (g_mode)
    {
    case hl_import:
    default:
        tmp = "Import";
        break;
    case hl_export:
        tmp = "Export";
        break;
    }

	Log("\n");

    // RipEnt Specific Settings
	Log("mode                [ %7s ] [ %7s ]\n", tmp, "N/A");
    Log("parse               [ %7s ] [ %7s ]\n", g_parse ? "on" : "off", DEFAULT_PARSE ? "on" : "off");

    Log("\n\n");
}

/*
 * ============
 * main
 * ============
 */
int             main(int argc, char** argv)
{
    int             i;
    double          start, end;

    g_Program = "ripent";

#ifdef RIPENT_PAUSE
	atexit (&pause);
#endif
#ifdef ZHLT_PARAMFILE
	int argcold = argc;
	char ** argvold = argv;
	{
		int argc;
		char ** argv;
		ParseParamFile (argcold, argvold, argc, argv);
		{
#endif
#ifdef ZHLT_CONSOLE
	if (InitConsole (argc, argv) < 0)
		Usage();
#endif
    if (argc == 1)
    {
        Usage();
    }

    for (i = 1; i < argc; i++)
    {
        if (!strcasecmp(argv[i], "-import"))
        {
            g_mode = hl_import;
        }
#ifdef ZHLT_CONSOLE
		else if (!strcasecmp(argv[i], "-console"))
		{
			if (i + 1 < argc)
				++i;
			else
				Usage();
		}
#endif
        else if (!strcasecmp(argv[i], "-export"))
        {
            g_mode = hl_export;
        }
		// g_parse: command line switch (-parse).
		// Added by: Ryan Gregg aka Nem
		else if(!strcasecmp(argv[i], "-parse"))
		{
			g_parse = true;
		}
        else if (!strcasecmp(argv[i], "-texdata"))
        {
            if (i + 1 < argc)	//added "1" .--vluzacn
            {
                int             x = atoi(argv[++i]) * 1024;

                //if (x > g_max_map_miptex) //--vluzacn
                {
                    g_max_map_miptex = x;
                }
            }
            else
            {
                Usage();
            }
        }
        else if (!strcasecmp(argv[i], "-lightdata"))
        {
            if (i + 1 < argc)	//added "1" .--vluzacn
            {
                int             x = atoi(argv[++i]) * 1024;

                //if (x > g_max_map_lightdata) //--vluzacn
                {
                    g_max_map_lightdata = x;
                }
            }
            else
            {
                Usage();
            }
        }
        else if (!strcasecmp(argv[i], "-chart"))
        {
            g_chart = true;
        }
        else if (!strcasecmp(argv[i], "-noinfo"))
        {
            g_info = false;
        }
#ifdef RIPENT_PAUSE
		else if (!strcasecmp(argv[i], "-pause"))
		{
			g_pause = true;
		}
#endif
#ifdef ZHLT_LANGFILE
		else if (!strcasecmp (argv[i], "-lang"))
		{
			if (i + 1 < argc)
			{
				LoadLangFile (argv[++i], argv[0]);
			}
			else
			{
				Usage();
			}
		}
#endif
		else if (argv[i][0] == '-') //--vluzacn
		{
			Log("Unknown option: '%s'\n", argv[i]);
			Usage ();
		}
        else
        {
            safe_strncpy(g_Mapname, argv[i], _MAX_PATH);
            StripExtension(g_Mapname);
            DefaultExtension(g_Mapname, ".bsp");
        }
    }

    if (g_mode == hl_undefined)
    {
        Log("Must specify either -import or -export\n"); //--vluzacn
        Usage();
    }

    if (!q_exists(g_Mapname))
    {
        Log("bspfile '%s' does not exist\n", g_Mapname); //--vluzacn
        Usage();
    }

#ifdef ZHLT_PARAMFILE
    LogStart(argcold, argvold);
	{
		int			 i;
		Log("Arguments: ");
		for (i = 1; i < argc; i++)
		{
			if (strchr(argv[i], ' '))
			{
				Log("\"%s\" ", argv[i]);
			}
			else
			{
				Log("%s ", argv[i]);
			}
		}
		Log("\n");
	}
#else
    LogStart(argc, argv);
#endif
	atexit(LogEnd);

	Settings();

    dtexdata_init();
    atexit(dtexdata_free);

    // BEGIN RipEnt
    start = I_FloatTime();

    switch (g_mode)
    {
    case hl_import:
		ReadBSP(g_Mapname);
        ReadEntities(g_Mapname);
#ifdef RIPENT_TEXTURE
		ReadTextures(g_Mapname);
#endif
        WriteBSP(g_Mapname);
        break;
    case hl_export:
		ReadBSP(g_Mapname);
        WriteEntities(g_Mapname);
#ifdef RIPENT_TEXTURE
		WriteTextures(g_Mapname);
#endif
        break;
    }

    if (g_chart)
        PrintBSPFileSizes();

    end = I_FloatTime();
    LogTimeElapsed(end - start);
    // END RipEnt
#ifdef ZHLT_PARAMFILE
		}
	}
#endif

    return 0;
}

// do nothing - we don't have params to fetch
void GetParamsFromEnt(entity_t* mapent) {}
