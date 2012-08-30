#include <stdio.h>
#include <stdlib.h>

#include "groupstuff.h"

/*! \mainpage Reference documentation for IDL-Grouplib

\author Volker Springel \n
        Max-Planck-Institute for Astrophysics \n
        Garching, Germany \n
        volker@mpa-garching.mpg.de \n
\n

This C-library is meant to be called by IDL (or other languages) and
provides an easy-to-use interface for the groupcatalogues generated by
the FOF algorithm of <b>P-GADGET2</b>.
*/



typedef struct
{
  unsigned int slen;
  unsigned int stype;
  char *s;
}
IDL_STRING;


typedef struct
{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  int flag_stellarage;
  int flag_metals;
  int hashtabsize;
  char fill[84];		/* fills to 256 Bytes */
}
io_header;



typedef struct
{
  int len;
  int file;
  int offset;
  int typecount[6];
  double typemass[6];
  float center[3];
  float sfr;
  float mdot;
  float mbh;
}
cat_data;


typedef struct
{
  int ID;
  int rank;
}
idsort_data, idsort_gasdata, idsort_starsdata;

static double BoxSize;

/*  Expected input values:
 *  char*    Output directory
 *  int      snapshot number
 */

int get_total_number_of_groups(int argc, void *argv[])
{
  int Ngroups, Nids, TotNgroups, NTask, Num;
  IDL_STRING *idl_s;
  char *OutputDir;
  char buf[1000];
  FILE *fd;

  idl_s = (IDL_STRING *) argv[0];
  OutputDir = idl_s->s;
  Num = *(int *) argv[1];

  sprintf(buf, "%s/groups_%03d/group_tab_%03d.0", OutputDir, Num, Num);
  if(!(fd = fopen(buf, "r")))
    {
      sprintf(buf, "%s/group_tab_%03d.0", OutputDir, Num);
      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s'\n", buf);
	  return -1;
	}
    }

  fread(&Ngroups, sizeof(int), 1, fd);
  fread(&Nids, sizeof(int), 1, fd);
  fread(&TotNgroups, sizeof(int), 1, fd);
  fread(&NTask, sizeof(int), 1, fd);
  fclose(fd);

  return TotNgroups;
}


int get_minimum_group_len(int argc, void *argv[])
{
  int Ngroups, Nids, TotNgroups, NTask, Num, GroupMinLen;
  IDL_STRING *idl_s;
  char *OutputDir;
  char buf[1000];
  FILE *fd;

  idl_s = (IDL_STRING *) argv[0];
  OutputDir = idl_s->s;
  Num = *(int *) argv[1];

  sprintf(buf, "%s/groups_%03d/group_tab_%03d.0", OutputDir, Num, Num);
  if(!(fd = fopen(buf, "r")))
    {
      sprintf(buf, "%s/group_tab_%03d.0", OutputDir, Num);
      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s'\n", buf);
	  return -1;
	}
    }

  fread(&Ngroups, sizeof(int), 1, fd);
  fread(&Nids, sizeof(int), 1, fd);
  fread(&TotNgroups, sizeof(int), 1, fd);
  fread(&NTask, sizeof(int), 1, fd);
  fseek(fd, sizeof(int) * Ngroups, SEEK_CUR);	/* skip GroupLen */
  fseek(fd, sizeof(int) * Ngroups, SEEK_CUR);	/* skip GroupOffset */
  fread(&GroupMinLen, sizeof(int), 1, fd);
  fclose(fd);

  return GroupMinLen;
}


int get_group_catalogue(int argc, void *argv[])
{
  int Ngroups, Nids, TotNgroups, NTask, Num, count, i, j, filenr;
  int *GroupLen, *GroupFileNr, *GroupNr, *GroupTypeLen;
  double *GroupTypeMass;
  float *GroupCenter, *GroupSfr;
  IDL_STRING *idl_s;
  cat_data *temp;
  char *OutputDir;
  char buf[1000];
  FILE *fd;

  idl_s = (IDL_STRING *) argv[0];
  OutputDir = idl_s->s;
  Num = *(int *) argv[1];

  GroupLen = (int *) argv[2];
  GroupFileNr = (int *) argv[3];
  GroupNr = (int *) argv[4];
  GroupTypeLen = (int *) argv[5];
  GroupTypeMass = (double *) argv[6];
  GroupCenter = (float *) argv[7];
  GroupSfr = (float *) argv[8];

  sprintf(buf, "%s/groups_%03d/group_tab_%03d.0", OutputDir, Num, Num);
  if(!(fd = fopen(buf, "r")))
    {
      sprintf(buf, "%s/group_tab_%03d.0", OutputDir, Num);
      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s'\n", buf);
	  return -1;
	}
    }

  fread(&Ngroups, sizeof(int), 1, fd);
  fread(&Nids, sizeof(int), 1, fd);
  fread(&TotNgroups, sizeof(int), 1, fd);
  fread(&NTask, sizeof(int), 1, fd);
  fclose(fd);


  for(filenr = 0, count = 0; filenr < NTask; filenr++)
    {
      sprintf(buf, "%s/groups_%03d/group_tab_%03d.%d", OutputDir, Num, Num, filenr);
      if(!(fd = fopen(buf, "r")))
	{
	  sprintf(buf, "%s/group_tab_%03d.%d", OutputDir, Num, filenr);
	  if(!(fd = fopen(buf, "r")))
	    {
	      printf("can't open file `%s'\n", buf);
	      return -1;
	    }
	}

      fread(&Ngroups, sizeof(int), 1, fd);
      fread(&Nids, sizeof(int), 1, fd);
      fread(&TotNgroups, sizeof(int), 1, fd);
      fread(&NTask, sizeof(int), 1, fd);

      fread(&GroupLen[count], sizeof(int), Ngroups, fd);

      /* skip offset table */
      fseek(fd, sizeof(int) * Ngroups, SEEK_CUR);

      fread(&GroupTypeLen[6*count], 6*sizeof(int), Ngroups, fd);
      fread(&GroupTypeMass[6*count], 6*sizeof(double), Ngroups, fd);
      fread(&GroupCenter[3*count], 3*sizeof(float), Ngroups, fd);
      fread(&GroupSfr[count], sizeof(float), Ngroups, fd);

      for(i = 0; i < Ngroups; i++)
	{
	  GroupFileNr[i + count] = filenr;
	  GroupNr[i + count] = i;
	}

      count += Ngroups;

      fclose(fd);
    }

  temp = malloc(sizeof(cat_data) * TotNgroups);

  for(i = 0; i < TotNgroups; i++)
    {
      temp[i].len = GroupLen[i];
      temp[i].file = GroupFileNr[i];
      temp[i].offset = GroupNr[i];
      for(j=0; j<6; j++)
	{
	  temp[i].typecount[j] = GroupTypeLen[6*i+j];
	  temp[i].typemass[j] = GroupTypeMass[6*i+j];
	}
      for(j=0; j<3; j++)
	temp[i].center[j] = GroupCenter[3*i+j];
      temp[i].sfr = GroupSfr[i];
    }

  qsort(temp, TotNgroups, sizeof(cat_data), id_sort_groups);

  for(i = 0; i < TotNgroups; i++)
    {
      GroupLen[i] = temp[i].len;
      GroupFileNr[i] = temp[i].file;
      GroupNr[i] = temp[i].offset;
      for(j=0; j<6; j++)
	{
	  GroupTypeLen[6*i+j] = temp[i].typecount[j];
	  GroupTypeMass[6*i+j] = temp[i].typemass[j];
	}

      for(j=0; j<3; j++)
	GroupCenter[3*i+j]= temp[i].center[j];
      GroupSfr[i] = temp[i].sfr;
    }

  free(temp);

  return TotNgroups;
}


int get_group_catalogue_bh(int argc, void *argv[])
{
  int Ngroups, Nids, TotNgroups, NTask, Num, count, i, j, filenr;
  int *GroupLen, *GroupFileNr, *GroupNr, *GroupTypeLen;
  double *GroupTypeMass;
  float *GroupCenter, *GroupSfr, *GroupMdot, *GroupMbh;
  IDL_STRING *idl_s;
  cat_data *temp;
  char *OutputDir;
  char buf[1000];
  FILE *fd;

  idl_s = (IDL_STRING *) argv[0];
  OutputDir = idl_s->s;
  Num = *(int *) argv[1];

  GroupLen = (int *) argv[2];
  GroupFileNr = (int *) argv[3];
  GroupNr = (int *) argv[4];
  GroupTypeLen = (int *) argv[5];
  GroupTypeMass = (double *) argv[6];
  GroupCenter = (float *) argv[7];
  GroupSfr = (float *) argv[8];
  GroupMdot =  (float *) argv[9];
  GroupMbh = (float *) argv[10];


  sprintf(buf, "%s/groups_%03d/group_tab_%03d.0", OutputDir, Num, Num);
  if(!(fd = fopen(buf, "r")))
    {
      sprintf(buf, "%s/group_tab_%03d.0", OutputDir, Num);
      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s'\n", buf);
	  return -1;
	}
    }

  fread(&Ngroups, sizeof(int), 1, fd);
  fread(&Nids, sizeof(int), 1, fd);
  fread(&TotNgroups, sizeof(int), 1, fd);
  fread(&NTask, sizeof(int), 1, fd);
  fclose(fd);


  for(filenr = 0, count = 0; filenr < NTask; filenr++)
    {
      sprintf(buf, "%s/groups_%03d/group_tab_%03d.%d", OutputDir, Num, Num, filenr);
      if(!(fd = fopen(buf, "r")))
	{
	  sprintf(buf, "%s/group_tab_%03d.%d", OutputDir, Num, filenr);
	  if(!(fd = fopen(buf, "r")))
	    {
	      printf("can't open file `%s'\n", buf);
	      return -1;
	    }
	}

      fread(&Ngroups, sizeof(int), 1, fd);
      fread(&Nids, sizeof(int), 1, fd);
      fread(&TotNgroups, sizeof(int), 1, fd);
      fread(&NTask, sizeof(int), 1, fd);

      fread(&GroupLen[count], sizeof(int), Ngroups, fd);

      /* skip offset table */
      fseek(fd, sizeof(int) * Ngroups, SEEK_CUR);

      fread(&GroupTypeLen[6*count], 6*sizeof(int), Ngroups, fd);
      fread(&GroupTypeMass[6*count], 6*sizeof(double), Ngroups, fd);
      fread(&GroupCenter[3*count], 3*sizeof(float), Ngroups, fd);
      fread(&GroupSfr[count], sizeof(float), Ngroups, fd);
      fread(&GroupMdot[count], sizeof(float), Ngroups, fd);
      fread(&GroupMbh[count], sizeof(float), Ngroups, fd);

      for(i = 0; i < Ngroups; i++)
	{
	  GroupFileNr[i + count] = filenr;
	  GroupNr[i + count] = i;
	}

      count += Ngroups;

      fclose(fd);
    }

  temp = malloc(sizeof(cat_data) * TotNgroups);

  for(i = 0; i < TotNgroups; i++)
    {
      temp[i].len = GroupLen[i];
      temp[i].file = GroupFileNr[i];
      temp[i].offset = GroupNr[i];
      for(j=0; j<6; j++)
	{
	  temp[i].typecount[j] = GroupTypeLen[6*i+j];
	  temp[i].typemass[j] = GroupTypeMass[6*i+j];
	}
      for(j=0; j<3; j++)
	temp[i].center[j] = GroupCenter[3*i+j];
      temp[i].sfr = GroupSfr[i];
      temp[i].mdot = GroupMdot[i];
      temp[i].mbh = GroupMbh[i];
    }

  qsort(temp, TotNgroups, sizeof(cat_data), id_sort_groups);

  for(i = 0; i < TotNgroups; i++)
    {
      GroupLen[i] = temp[i].len;
      GroupFileNr[i] = temp[i].file;
      GroupNr[i] = temp[i].offset;
      for(j=0; j<6; j++)
	{
	  GroupTypeLen[6*i+j] = temp[i].typecount[j];
	  GroupTypeMass[6*i+j] = temp[i].typemass[j];
	}

      for(j=0; j<3; j++)
	GroupCenter[3*i+j]= temp[i].center[j];
      GroupSfr[i] = temp[i].sfr;
      GroupMdot[i] = temp[i].mdot;
      GroupMbh[i] = temp[i].mbh;
    }

  free(temp);

  return TotNgroups;
}



int id_sort_groups(const void *a, const void *b)
{
  if(((cat_data *) a)->len > ((cat_data *) b)->len)
    return -1;

  if(((cat_data *) a)->len < ((cat_data *) b)->len)
    return +1;

  return 0;
}

int get_total_particle_count(int argc, void *argv[])
{
  int Num, i, numpart, *Files, dummy;
  IDL_STRING *idl_s;
  char *OutputDir, *Snapbase;
  char buf [1000];
  FILE *fd;
  io_header header;

  idl_s = (IDL_STRING *) argv[0];
  OutputDir = idl_s->s;
  Num = *(int *) argv[1];
  idl_s = (IDL_STRING *) argv[2];
  Snapbase = idl_s->s;

  Files = (int *) argv[3];

  sprintf(buf, "%s/%s_%03d.0", OutputDir, Snapbase, Num);
  if(!(fd = fopen(buf, "r")))
    {
      sprintf(buf, "%s/%s_%03d", OutputDir, Snapbase, Num);
      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s' or '%s.0'\n", buf, buf);
	  return -1;
	}
    }
  fread(&dummy, sizeof(int), 1, fd);
  fread(&header, sizeof(io_header), 1, fd);
  fread(&dummy, sizeof(int), 1, fd);
  fclose(fd);


  *Files = header.num_files;
  
  for(i=0, numpart=0; i<6; i++)
    numpart += header.npartTotal[i];

  printf("Total particles = %d\n", numpart);

  return numpart;
}

int get_gas_particle_count(int argc, void *argv[])
{
  int Num, numpart, *Files, dummy;
  IDL_STRING *idl_s;
  char *OutputDir, *Snapbase;
  char buf [1000];
  FILE *fd;
  io_header header;

  idl_s = (IDL_STRING *) argv[0];
  OutputDir = idl_s->s;
  Num = *(int *) argv[1];
  idl_s = (IDL_STRING *) argv[2];
  Snapbase = idl_s->s;

  Files = (int *) argv[3];

  sprintf(buf, "%s/%s_%03d.0", OutputDir, Snapbase, Num);
  if(!(fd = fopen(buf, "r")))
    {
      sprintf(buf, "%s/%s_%03d", OutputDir, Snapbase, Num);
      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s' or '%s.0'\n", buf, buf);
	  return -1;
	}
    }
  fread(&dummy, sizeof(int), 1, fd);
  fread(&header, sizeof(io_header), 1, fd);
  fread(&dummy, sizeof(int), 1, fd);
  fclose(fd);


  *Files = header.num_files;
  
    numpart=0; 
    numpart += header.npartTotal[0];

  printf("Gas particles = %d\n", numpart);

  return numpart;
}

int get_dm_particle_count(int argc, void *argv[])
{
  int Num, numpart, *Files, dummy;
  IDL_STRING *idl_s;
  char *OutputDir, *Snapbase;
  char buf [1000];
  FILE *fd;
  io_header header;

  idl_s = (IDL_STRING *) argv[0];
  OutputDir = idl_s->s;
  Num = *(int *) argv[1];
  idl_s = (IDL_STRING *) argv[2];
  Snapbase = idl_s->s;

  Files = (int *) argv[3];

  sprintf(buf, "%s/%s_%03d.0", OutputDir, Snapbase, Num);
  if(!(fd = fopen(buf, "r")))
    {
      sprintf(buf, "%s/%s_%03d", OutputDir, Snapbase, Num);
      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s' or '%s.0'\n", buf, buf);
	  return -1;
	}
    }
  fread(&dummy, sizeof(int), 1, fd);
  fread(&header, sizeof(io_header), 1, fd);
  fread(&dummy, sizeof(int), 1, fd);
  fclose(fd);


  *Files = header.num_files;
  
    numpart=0; 
    numpart += header.npartTotal[1]+header.npartTotal[2]+header.npartTotal[3];

  printf("Dm particles = %d\n", numpart);

  return numpart;
}


int get_mass_particle_count(int argc, void *argv[])
{
  int Num, numpart, *Files, dummy, i;
  IDL_STRING *idl_s;
  char *OutputDir, *Snapbase;
  char buf [1000];
  FILE *fd;
  io_header header;

  idl_s = (IDL_STRING *) argv[0];
  OutputDir = idl_s->s;
  Num = *(int *) argv[1];
  idl_s = (IDL_STRING *) argv[2];
  Snapbase = idl_s->s;

  Files = (int *) argv[3];

  sprintf(buf, "%s/%s_%03d.0", OutputDir, Snapbase, Num);
  if(!(fd = fopen(buf, "r")))
    {
      sprintf(buf, "%s/%s_%03d", OutputDir, Snapbase, Num);
      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s' or '%s.0'\n", buf, buf);
	  return -1;
	}
    }
  fread(&dummy, sizeof(int), 1, fd);
  fread(&header, sizeof(io_header), 1, fd);
  fread(&dummy, sizeof(int), 1, fd);
  fclose(fd);


  *Files = header.num_files;
  
  for(i=0, numpart=0; i < 6; i++)
    {
      if( (header.mass[i] == 0) && (header.npart[i] > 0))
	numpart += header.npartTotal[i];
    }

  printf("Particles with variable mass= %d\n", numpart);

  return numpart;
}

int get_stars_particle_count(int argc, void *argv[])
{
  int Num, numpart, *Files, dummy;
  IDL_STRING *idl_s;
  char *OutputDir, *Snapbase;
  char buf [1000];
  FILE *fd;
  io_header header;

  idl_s = (IDL_STRING *) argv[0];
  OutputDir = idl_s->s;
  Num = *(int *) argv[1];
  idl_s = (IDL_STRING *) argv[2];
  Snapbase = idl_s->s;

  Files = (int *) argv[3];

  sprintf(buf, "%s/%s_%03d.0", OutputDir, Snapbase, Num);
  if(!(fd = fopen(buf, "r")))
    {
      sprintf(buf, "%s/%s_%03d", OutputDir, Snapbase, Num);
      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s' or '%s.0'\n", buf, buf);
	  return -1;
	}
    }
  fread(&dummy, sizeof(int), 1, fd);
  fread(&header, sizeof(io_header), 1, fd);
  fread(&dummy, sizeof(int), 1, fd);
  fclose(fd);


  *Files = header.num_files;
  
    numpart=0; 
    numpart += header.npartTotal[4];

  printf("Stars particles = %d\n", numpart);

  return numpart;
}


int get_particle_data(int argc, void *argv[])
{
  int Num, fnr, i, j, Files, NumPart, *Type, *RankID;
  unsigned long *ID;
  int dummy, count, len;
  IDL_STRING *idl_s;
  char *OutputDir, *Snapbase;
  char buf [1000];
  FILE *fd;
  float *Pos;
  float *Vel;
  float *U;
  float *Rho;
  float *Mass;
  float *Ne;
  float *NH;
  float *Hsml;
  float *SFR;
  float *stellarage;
  float *Metallicity;
  int gaslen, gascount, *RankIDgas, NumGas, NumDm;
  int masslen, masscount, *RankIDmass, NumMass;
  int starslen, starscount, *RankIDstars, NumStars;

  idsort_data *iddat;
  idsort_gasdata *iddat_gas;
  idsort_starsdata *iddat_stars;
 
  io_header header;

  idl_s = (IDL_STRING *) argv[0];
  OutputDir = idl_s->s;
  Num = *(int *) argv[1];
  idl_s = (IDL_STRING *) argv[2];
  Snapbase = idl_s->s;

  Files = *(int *) argv[3];
  NumPart = *(int *) argv[4];
  NumGas = *(int *) argv[5];
  NumDm = *(int *) argv[6];
  NumMass = *(int *) argv[7];
  NumStars = *(int *) argv[8];
  Pos = (float *) argv[9];
  Vel = (float *) argv[10];
  ID = (unsigned long *) argv[11];
  Mass =(float *) argv[12];
  U = (float *) argv[13];
  Rho = (float *) argv[14];
  Ne = (float *) argv[15];
  NH = (float *) argv[16];
  Hsml = (float *) argv[17];
  SFR = (float *) argv[18];
  stellarage = (float *) argv[19];
  Metallicity = (float *) argv[20];
  Type = (int *) argv[21];
  RankID =  (int *) argv[22];
  RankIDgas =  (int *) argv[23];
  RankIDmass = (int *) argv[24];
  RankIDstars = (int *) argv[25];

  for(fnr = 0, count = 0, gascount=0, masscount=0, starscount=0; fnr < Files; fnr++)
    {
      if(Files > 1)
	sprintf(buf, "%s/%s_%03d.%d", OutputDir, Snapbase, Num, fnr);
      else
	sprintf(buf, "%s/%s_%03d", OutputDir, Snapbase, Num);
      
      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s'\n", buf);
	  return -1;
	}
      fread(&dummy, sizeof(int), 1, fd);
      fread(&header, sizeof(io_header), 1, fd);
      fread(&dummy, sizeof(int), 1, fd);

      BoxSize = header.BoxSize;

      for(i=0, len=0; i < 6; i++)
	len += header.npart[i];

      gaslen = 0;
      gaslen = header.npart[0];

      starslen = 0;
      starslen = header.npart[4];

      for(i=0, masslen=0; i < 6; i++)
	if(header.mass[i] == 0)
	    masslen += header.npart[i];
   
      fread(&dummy, sizeof(int), 1, fd);
      fread(&Pos[3 * count], sizeof(float), 3 * len, fd);
      fread(&dummy, sizeof(int), 1, fd);

      fread(&dummy, sizeof(int), 1, fd);
      fread(&Vel[3 * count], sizeof(float), 3 * len, fd);
      fread(&dummy, sizeof(int), 1, fd);
      
      fread(&dummy, sizeof(int), 1, fd);
      fread(&ID[count], sizeof(int), len, fd);
      fread(&dummy, sizeof(int), 1, fd);

      if(masslen > 0)
	{
	  fread(&dummy, sizeof(int), 1, fd);
	  fread(&Mass[masscount], sizeof(float), masslen, fd);
	  fread(&dummy, sizeof(int), 1, fd);
	}

      if(gaslen > 0)
	{
	  fread(&dummy, sizeof(int), 1, fd);
	  fread(&U[gascount], sizeof(float), gaslen, fd);
	  fread(&dummy, sizeof(int), 1, fd);

	  fread(&dummy, sizeof(int), 1, fd);
	  fread(&Rho[gascount], sizeof(float), gaslen, fd);
	  fread(&dummy, sizeof(int), 1, fd);

	  fread(&dummy, sizeof(int), 1, fd);
	  fread(&Ne[gascount], sizeof(float), gaslen, fd);
	  fread(&dummy, sizeof(int), 1, fd);

	  fread(&dummy, sizeof(int), 1, fd);
	  fread(&NH[gascount], sizeof(float), gaslen, fd);
	  fread(&dummy, sizeof(int), 1, fd);

	  fread(&dummy, sizeof(int), 1, fd);
	  fread(&Hsml[gascount], sizeof(float), gaslen, fd);
	  fread(&dummy, sizeof(int), 1, fd);

	  if(starslen > 0)
	    {
	      fread(&dummy, sizeof(int), 1, fd);
	      fread(&SFR[gascount], sizeof(float), gaslen, fd);
	      fread(&dummy, sizeof(int), 1, fd);

	      fread(&dummy, sizeof(int), 1, fd);
	      fread(&stellarage[starscount], sizeof(float), starslen, fd);
	      fread(&dummy, sizeof(int), 1, fd);

	      fread(&dummy, sizeof(int), 1, fd);
	      fread(&Metallicity[gascount+starscount], sizeof(float), gaslen+starslen, fd);
	      fread(&dummy, sizeof(int), 1, fd);

	    }

	}

      for(i=0; i<6; i++)
	for(j=0; j< header.npart[i]; j++)
	  {
	    Type[count] = i;
	    RankID[count]= count;
	 
	    if (i == 4)
	      {
		RankIDstars[starscount] = count;
		starscount++;
	      }
	 
	    count++;
	    
	    if( (header.mass[i] == 0) && (header.npart[i] > 0) )
	      {
		RankIDmass[masscount]= masscount;
		masscount++;
	      }	    

	  }

      for(j=0; j< header.npart[0]; j++)
	{
	  RankIDgas[gascount] = gascount;
	  gascount++;
	}

      fclose(fd);
    }


  iddat = malloc(sizeof(idsort_data)*NumPart);
  iddat_gas = malloc(sizeof(idsort_gasdata)*NumGas);
  iddat_stars = malloc(sizeof(idsort_starsdata)*NumStars);

  for(i=0; i<NumPart;i++)
    {
      iddat[i].ID = ID[i];
      iddat[i].rank = i;
    }

  for(i=0; i<NumGas;i++)
    {
      iddat_gas[i].ID = ID[i];
      iddat_gas[i].rank = i;
    }

  for(i=0; i<NumStars;i++)
    {
      iddat_stars[i].ID = ID[NumGas + NumDm + i];
      iddat_stars[i].rank = i + NumGas + NumDm;
    }

  qsort(iddat, NumPart, sizeof(idsort_data), iddat_sort_compare);
  qsort(iddat_gas, NumGas, sizeof(idsort_gasdata), iddat_gas_sort_compare);
  qsort(iddat_stars, NumStars, sizeof(idsort_starsdata), iddat_stars_sort_compare);

  for(i=0; i<NumPart; i++)
    RankID[i] = iddat[i].rank;

  for(i=0; i<NumGas; i++)
    RankIDgas[i] = iddat_gas[i].rank;

  for(i=0; i<NumStars; i++)
    RankIDstars[i] = iddat_stars[i].rank;


  free(iddat);
  free(iddat_gas);
  free(iddat_stars);

  return 0;
}


int iddat_sort_compare(const void *a, const void *b)
{
  if(((idsort_data *) a)->ID < ((idsort_data *) b)->ID)
    return -1;

  if(((idsort_data *) a)->ID > ((idsort_data *) b)->ID)
    return +1;

  return 0;
}

int iddat_gas_sort_compare(const void *a, const void *b)
{
  if(((idsort_gasdata *) a)->ID < ((idsort_gasdata *) b)->ID)
    return -1;

  if(((idsort_gasdata *) a)->ID > ((idsort_gasdata *) b)->ID)
    return +1;

  return 0;
}

int iddat_stars_sort_compare(const void *a, const void *b)
{
  if(((idsort_starsdata *) a)->ID < ((idsort_starsdata *) b)->ID)
    return -1;

  if(((idsort_starsdata *) a)->ID > ((idsort_starsdata *) b)->ID)
    return +1;

  return 0;
}


int get_group_indices(int argc, void *argv[])
{
  int Num, NTask;
  int i, j, n, NumPart, *RankID;
  unsigned long *ID, *RankIDstars, *Indicesstars;
  IDL_STRING *idl_s;
  char *OutputDir, *Snapbase;
  char buf[1000];
  FILE *fd;
  int GroupNr, FileNr, GroupLen, Ngroups, TotNgroups, Nids;
  int offset, *Indices;
  unsigned long *ids, *idsgas, *idsstars;
  int NumDm, NumGas, *RankIDgas, *Indicesgas, GroupLengas;
  int NumStars, GroupLenstars, GroupLendm;


  idl_s = (IDL_STRING *) argv[0];
  OutputDir = idl_s->s;
  Num = *(int *) argv[1];
  idl_s = (IDL_STRING *) argv[2];
  Snapbase = idl_s->s;

  GroupNr = *(int *) argv[3];
  FileNr = *(int *) argv[4];;
  GroupLen = *(int *) argv[5];
  GroupLengas = *(int *) argv[6];
  GroupLendm = *(int *) argv[7];
  GroupLenstars = *(int *) argv[8];
  NumPart  = *(int *) argv[9];
  NumGas = *(int *) argv[10];
  NumDm = *(int *) argv[11];
  NumStars = *(int *) argv[12];
  ID = (unsigned long *) argv[13];
  RankID =  (int *) argv[14];
  RankIDgas =  (int *) argv[15];
  RankIDstars =  (unsigned long *) argv[16];
  Indices = (int *) argv[17];
  Indicesgas = (int *) argv[18];
  Indicesstars = (unsigned long *) argv[19];

  sprintf(buf, "%s/groups_%03d/group_tab_%03d.%d", OutputDir, Num, Num, FileNr);
  if(!(fd = fopen(buf, "r")))
    {
      sprintf(buf, "%s/group_tab_%03d.%d", OutputDir, Num, FileNr);
      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s'\n", buf);
	  return -1;
	}
    }

  fread(&Ngroups, sizeof(int), 1, fd);
  fread(&Nids, sizeof(int), 1, fd);
  fread(&TotNgroups, sizeof(int), 1, fd);
  fread(&NTask, sizeof(int), 1, fd);
  /* skip group-len table */
  fseek(fd, sizeof(int) * Ngroups, SEEK_CUR);
  fseek(fd, sizeof(int) * GroupNr, SEEK_CUR);
  fread(&offset, sizeof(int), 1, fd);
  fclose(fd);

  ids = malloc(GroupLen * sizeof(unsigned long));
  idsgas = malloc(GroupLen * sizeof(unsigned long));
  idsstars = malloc(GroupLen * sizeof(unsigned long));

  sprintf(buf, "%s/groups_%03d/group_ids_%03d.%d", OutputDir, Num, Num, FileNr);
  if(!(fd = fopen(buf, "r")))
    {
      sprintf(buf, "%s/group_ids_%03d.%d", OutputDir, Num, FileNr);
      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s'\n", buf);
	  return -1;
	}
    }
  fread(&Ngroups, sizeof(int), 1, fd);
  fread(&Nids, sizeof(int), 1, fd);
  fread(&TotNgroups, sizeof(int), 1, fd);
  fread(&NTask, sizeof(int), 1, fd);

  fseek(fd, sizeof(int) * offset, SEEK_CUR);
  fread(ids, sizeof(unsigned long), GroupLen, fd);
  fclose(fd);

  for(i=0; i < GroupLen; i++)
    {
      idsgas[i]=ids[i];
      idsstars[i]=ids[i];
    }

   qsort(ids, GroupLen, sizeof(unsigned long), id_sort_compare_key);

  for(i = 0, n = 0; i < GroupLen; i++)
    {
      while(n < NumPart && ID[RankID[n]] != ids[i])
	n++;

      if(n >= NumPart)
	{
	  printf("We have a mismatch! (i=%d) something is wrong here.\n", i);
	  return -1;
	}
      
      Indices[i] = RankID[n];
    }

  for(i = 0, n = 0; i < GroupLengas; i++)
    {
      while(n < NumGas && ID[RankIDgas[n]] != idsgas[i])
	n++;

      if(n >= NumGas)
	{
	  printf("We have a mismatch again! (i=%d) something is wrong here.\n", i);
	  return -1;
	}
      
      Indicesgas[i] = RankIDgas[n];
    }

  for(i = 0, n = 0, j = NumGas+NumDm; i < GroupLenstars; i++)
    {
      while(n < NumStars && ID[RankIDstars[n]] != idsstars[i+GroupLengas+GroupLendm])
      	{
  	  n++;
	  j++;
  	}
      
      if(j >= NumStars+NumGas+NumDm+1)
     	{
            	  printf("We have a mismatch again!! (i=%d,n=%d) something is wrong here.\n",i, j);
     	  return -1;
     	}
      
      Indicesstars[i] = RankIDstars[n];
    }

  free(ids);
  free(idsgas);
  free(idsstars);
  return 0;
}



int id_sort_compare_key(const void *a, const void *b)
{
  if(*((int *) a) < *((int *) b))
    return -1;

  if(*((int *) a) > *((int *) b))
    return +1;

  return 0;
}
