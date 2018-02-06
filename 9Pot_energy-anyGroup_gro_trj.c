/* 
   This code:
             Analyze the potential energy of an ion at various
	     positions across the channel
   Input:
             1. .xtc or .trr file of the coordinates of SOL, NA/Cl
	     2. .dat file of the Wall atom position

   Note:                                           
           1. Energy in unit of kJ/mol                
	   2. Internal variables:                  

   requirements:
   1. order of atoms in trajectory file:  Na -> Cl -> SOL
   Rui Qiao
   
   v1.0 Oct, 23, 2002
   v2.0 Nov, 19, 2002  -> wall-ion electrostatic enabled
   v2.2 Nov, 25, 2002  -> ion-ion electrostatic enabled
                       -> energy unit now in kJ/mol
calculate the energy between any 2 groups   --------- 9.19.2010
from now on, each new code will have a name with 9 as the start.
*/

#include "string.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "copyrite.h"
#include "statutil.h"
#include "tpxio.h"
#include "Subfunc/create_matrix.h"
#include "Subfunc/gro_operate.h"
#include "Subfunc/itp_operate.h"
#include "Subfunc/useful_sf.h"

int main(int argc,char *argv[])
{
  static char *desc[] = {"Ion - Water/ion/wall interaction in a slit pore, cylindrical cutoff"};
  
  static int preframe  = 0;
  static int endframe  = 9999999;
  static int nStart1   = -1;        // starting index of reference atoms
  static int nStart2   = -1;    
  static int nm1       = 840;
  static int nm2       = 840;
  static int tN1       = 1;
  static int tN2       = 3;
  static int n_NCM     = 0;
  static int nbin      = 1;
  static int DIMN      = 2;
  static int combRule  = 1;
  static int nC6_12    = 1;

  static float cutoff=1.2;
  static float lowP = 0;
  static float upP  = 6;


  static char *ipgro ="start.gro",*frzNm = "a,b,c",*resSel = "a,b,c";
  static char *itpF ="FField_energy.itp";
  static char *itpS ="MSD_charge_para.dat";
  static float kbT   = 1.0;//6022.1367/(1.38*298);   /* use kJ/mol as unit */  
  t_pargs pa[] = {
    { "-gro", FALSE, etSTR,  {&ipgro},
      "The corresponding GRO file."  
    },
    { "-itpF", FALSE, etSTR,  {&itpF},
      "The itp file of Force Field."  
    },
    { "-itpS", FALSE, etSTR,  {&itpS},
      "The itp file of molecule structure."  
    },
    { "-nstart1",FALSE, etINT, {&nStart1},
      "Start number of atoms in trj file, if minus, in gro file"
    },
    { "-nstart2",FALSE, etINT, {&nStart2},
      "Start number of atoms in trj file, if minus, in gro file"
    },
    { "-nm1",FALSE, etINT, {&nm1},
      "Number of reference molecules"
    },
    { "-nm2",FALSE, etINT, {&nm2},
      "Number of selected molecules"
    },
    { "-tN1",FALSE, etINT, {&tN1},
      "the type of reference atom or molecule"
    },
    { "-tN2",FALSE, etINT, {&tN2},
      "the type of selected atom or molecule"
    },
    { "-NCM",FALSE, etINT, {&n_NCM},
      "center of molecule\n0: number; 1: charge; 2: mass"
    },
    { "-pre",FALSE, etINT, {&preframe},
      "Number of preframe"
    },
    { "-end",FALSE, etINT, {&endframe},
      "End number of trajectory set"
    },
    { "-dim",FALSE, etINT, {&DIMN},
      "the direction of bin sets made"
    },
    { "-c612",FALSE, etINT, {&nC6_12},
      "sigma/epsilon or C6/C12 in FF file\n0:sigma/epsilon\n1:C6/C12 "
    },
    { "-combR",FALSE, etINT, {&combRule},
      "combination rule\n1: mean by production; 2: mean by addation"
    },
    { "-nbin",FALSE, etINT, {&nbin},
      "number of bin sets along dim-direction"
    },
    { "-low",FALSE, etREAL, {&lowP},
      "lower position (nm)"
    },
    { "-up",FALSE, etREAL, {&upP},
      "upper position (nm)"
    },
    { "-Rcutoff",FALSE, etREAL, {&cutoff},
      "cutoff in LJ and Coulombic interactions(nm)"
    }
  };
   

  t_topology top;
  char       title[STRLEN];
  t_trxframe fr;
  rvec       *xtop;
  matrix     box;
  int        status;
  int        flags = TRX_READ_X | TRX_READ_V;   /* read velocity too! */
  

  t_filenm fnm[] = {
    { efTPS,  NULL,  NULL, ffREAD },   /* this is for the topology */
    { efTRX, "-f", NULL, ffREAD }      /* and this for the trajectory */
  };

#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

  //  ------------------------------------------------ my code -------------------------------------------//
  // read .gro file to know the molecule types
  FILE *inputgro,*inputFF, *inputStr, *file_energyPot;
  char name0[100];

  int  i,j,k,step,m,me,nType,max_atoms,totFrame = 1;
  double  ***pos1,***pos2,**groPos;
  double  *nCount,**energyPot,**Cij6,**Cij12;

  double dbin,dVol,dt,tmpt=0;
  double L[3];

  // for obtaining itp (structure and FF) information
  int  atoms1,atoms2,*numAt;
  char **itpMolResName,**itpMolscr, ***itpAtType, ***itpAtName;
  double ***atNCM, **mol_NCM,**sigma,**epsilon;
  // get information from gro file if necessary
  char **MolResName,**MolAtTy,**MolRName;
  int  ndxN,nTgro,ntotal,nSlt,*nStart,nMol[60],natoms[60];

  nType         = 10; 
  max_atoms     = 100; // max # of atoms in one molecule/ion
  itpMolResName = CreateMatrixChar(nType,10);               // residue name
  itpMolscr     = CreateMatrixChar(nType,10);               // describe of molecule/ion
  itpAtType     = CreateMatrixChar_3d(nType,max_atoms,10);  // atom type in FF
  itpAtName     = CreateMatrixChar_3d(nType,max_atoms,10);  // atom name
  atNCM  = CreateMatrix_3d(3,nType,max_atoms);              // number center, charge, mass of atom in each molecule/ion 
  mol_NCM= CreateMatrix(3,nType);                           // total number center, charge, mass of each molecule/ion 
  numAt  = CreateColumn(nType);                             // number of atoms in each molecule/ion 
  sigma  = CreateMatrix(nType,max_atoms);                   // sigma   in LJ-FF of each atom
  epsilon= CreateMatrix(nType,max_atoms);                   // epsilon in LJ-FF of each atom

  // get the parameter all types of molecule/ion in the structure file
  itp_str_FF_smart(itpF,itpS,itpMolResName,itpMolscr,itpAtType,itpAtName,atNCM,mol_NCM,sigma,epsilon,numAt,nType,nC6_12);
  if(tN1>nType || tN2>nType){
    printf("\nNOTE: there are NO types %d & %d both in MSD_charge_para.dat.\n",tN1,tN2);
    exit(0);
  }

  tN1 --; tN2 --; 
  atoms1 = numAt[tN1];  atoms2 = numAt[tN2];  // number of atoms in reference/selected molecule

  // get information from gro file if necessary
  if(nStart1 < 0 ||  nStart2 < 0){
    inputgro   = fopen(ipgro,"r"); 
    MolResName = CreateMatrixChar(nType,10);
    Read_gro(inputgro,MolResName,nMol,natoms,&nTgro); 
    rewind(inputgro); 
    nStart = CreateColumn(nTgro);
    for(i=0;i<nTgro-1;i++)   nStart[i+1] = nStart[i]+nMol[i]*natoms[i];
    ntotal =  nStart[nTgro-1]+nMol[nTgro-1]*natoms[nTgro-1];
    
    // get the position from gro file
    groPos = CreateMatrix(ntotal,3); 
    Read_gro_pos(inputgro,ntotal,groPos,L); 
   
    if(nStart1 < 0){      // reference group
      nSlt = 99999;
      for(i=0;i<nTgro;i++)
	if(strcmp(itpMolResName[tN1],str_rm_blank(MolResName[i]))==0) nSlt = i;
      if(nSlt == 99999){ printf("There is NO molecule named [%s]!!\n",itpMolResName[tN1]); exit(0);}
      nm1 = nMol[nSlt];
      if(atoms1 != natoms[nSlt]){ printf("The # of atoms in [%s] does NOT match!!\n",itpMolResName[tN1]); exit(0);}
      // get the coordinates of atoms in reference group
      pos1  = CreateMatrix_3d(nm1,atoms1,3);
      for(m=0;m<nm1;m++)
	for(j=0;j<atoms1;j++){
	  ndxN = nStart[nSlt] + atoms1*m+j;
	  for(k=0;k<3;k++)    pos1[m][j][k] = groPos[ndxN][k];
	}
    }
    
    if(nStart2 < 0){      // selected group
      nSlt = 99999;
      for(i=0;i<nTgro;i++)
	if(strcmp(itpMolResName[tN2],str_rm_blank(MolResName[i]))==0) nSlt = i;
      if(nSlt == 99999){ printf("There is NO molecule named [%s]!!\n",itpMolResName[tN2]); exit(0);}
      nm2 = nMol[nSlt];
      if(atoms2 != natoms[nSlt]){ printf("The # of atoms in [%s] does NOT match!!\n",itpMolResName[tN2]); exit(0);}
      // get the coordinates of atoms in reference group
      pos2  = CreateMatrix_3d(nm2,atoms2,3);
      for(m=0;m<nm2;m++)
	for(j=0;j<atoms2;j++){
	  ndxN = nStart[nSlt] + atoms2*m+j;
	  for(k=0;k<3;k++)    pos2[m][j][k] = groPos[ndxN][k];
	  //printf("coordinates:        %-8f    %-8f   %-8f\n",pos2[m][j][0],pos2[m][j][1],pos2[m][j][2]);
	}
    }
  }
  printf("\n            Molecule_type   number   atomsInMol   nStart\n");
  printf("reference:     %5s         %-8d    %-8d   %-8d\n",itpMolResName[tN1],nm1,atoms1,nStart1);
  printf("selected:      %5s         %-8d    %-8d   %-8d\n",itpMolResName[tN2],nm2,atoms2,nStart2);
     
  dbin  = (upP-lowP)/nbin; 

  Cij6        = CreateMatrix(atoms1,atoms2);  // note:  first is reference, second is selected
  Cij12       = CreateMatrix(atoms1,atoms2);
  energyPot   = CreateMatrix(2,nbin);
  nCount      = CreateVector(nbin);
  // initialize
  for(i=0;i<nbin;i++){
    nCount[i] = 0;
    for(j=0;j<2;j++)
      energyPot[j][i] = 0.0;
  }
  
  // get C6, C12 between atoms i and j
  computing_Cij_6_12(sigma[tN1],epsilon[tN1],sigma[tN2],epsilon[tN2],Cij6,Cij12,atoms1,atoms2,combRule);
  // compute the LJ/Coulombic potential between selected and reference groups
  if(nStart1 < 0 &&  nStart2 < 0)
    computing_pot(pos1,pos2,Cij6,Cij12,atNCM[n_NCM][tN1],atNCM[n_NCM][tN2],atNCM[1][tN1],atNCM[1][tN2],mol_NCM[n_NCM][tN1],mol_NCM[n_NCM][tN2],
		  energyPot,nCount,lowP,upP,dbin,cutoff,L,nm1,atoms1,nm2,atoms2,DIMN);
 

  if(nStart1 >= 0)   pos1  = CreateMatrix_3d(nm1,atoms1,3);
  if(nStart2 >= 0)   pos2  = CreateMatrix_3d(nm2,atoms2,3);
  
  if(nStart1 >= 0 || nStart2 >= 0){
    read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&xtop,NULL,box,TRUE);
    sfree(xtop);
    read_first_frame(&status,ftp2fn(efTRX,NFILE,fnm),&fr,flags);
    
    for(step=0;step<preframe;step++)     read_next_frame(status,&fr);
    tmpt= fr.time;
    do{ 
      L[0]   = fr.box[XX][XX]; 
      L[1]   = fr.box[YY][YY];
      L[2]   = fr.box[ZZ][ZZ];
      if(nStart1 >= 0) 
	load_position(fr,nm1,atoms1,nStart1,pos1,L); 
      
      if(tN1 == tN2){
	for(i=0;i<nm1;i++)
	  for(j=0;j<atoms1;j++)
	    for(k=0;k<3;k++)
	      pos2[i][j][k] = pos1[i][j][k];
      }else
	if(nStart2 >= 0) 
	  load_position(fr,nm2,atoms2,nStart2,pos2,L);
      // use the coordinates from trajectory file
      computing_pot(pos1,pos2,Cij6,Cij12,atNCM[n_NCM][tN1],atNCM[n_NCM][tN2],atNCM[1][tN1],atNCM[1][tN2],mol_NCM[n_NCM][tN1],mol_NCM[n_NCM][tN2],
		    energyPot,nCount,lowP,upP,dbin,cutoff,L,nm1,atoms1,nm2,atoms2,DIMN);

      step++;
      if(step==endframe) break;
    }while (read_next_frame(status,&fr));
    if(step<preframe+10){printf("\nInput trajactory file only has few frames!\n"); exit(0);}; 
    totFrame =    step-preframe;
    printf("\nTotal atom#:    %10d\n",fr.natoms);
    printf("bin size (nm): %10.3f\n",dbin);
    printf("\nTotal frame:    %10d\n",totFrame);
    printf("Time  step(ps): %10.3f\n",fr.tpf - fr.tppf);
    printf("Total time(ps): %10.3f\n",fr.time-tmpt);


  }
  dVol = L[0]*L[1]*L[2]*dbin/L[DIMN];
  //  printf("\nI'm HERE!\n");
  // report data
  //  sprintf(name0, "EnergyPot_bin%d_pos%.1f-%.1f_type%d-%d.dat",nbin,lowP,upP,tN1+1,tN2+1);
  sprintf(name0, "EnergyPot_bin%d_pos%.1f-%.1f_cutoff%.2f_%s_%s.dat",nbin,lowP,upP,cutoff,itpMolscr[tN1],itpMolscr[tN2]);
  printf("The output name is: %s! for m.s.d based on molecule center!!\n",name0);
  file_energyPot = fopen(name0,"w");
  for(i=0;i<nbin;i++)
    if(nCount[i] > 0)
      fprintf(file_energyPot,"%8.4e  %8.4e  %8.4e  %8.4e\n",dbin*(i+0.5)+lowP,nCount[i]/dVol/totFrame,energyPot[0][i]/nCount[i]/kbT,energyPot[1][i]/nCount[i]/kbT);
  fclose(file_energyPot);
  
  return 0;
}

