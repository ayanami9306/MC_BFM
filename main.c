#include <stdio.h>
#include "MC_Library.h"

#define ARM1_TYPE 1
#define ARM2_TYPE 2

char trajName[200], energyProfileName[200], datFileName[200], stateFileName[200];

unsigned int mcTimeMAX, mcInterval, mcTrial, mcTime, nMolecule;
unsigned int B_LENGTH, ARM_LENGTH[2];

__declspec(align(64)) SEGMENT molecule[20000];
__declspec(align(16)) int system_lattice_3D[BOUNDARY_MAX][BOUNDARY_MAX][BOUNDARY_Z];

double epsylon;
int nParticle;

FILE *out2;
FILE *fp_chiN;

void INPUT_PARAMETER(char *[]);
void SET_INITIAL_COORDINATE(char *[]);
inline void MetropolisMonteCarloEngine();
int Object_MCS = -1;

int main(int argc, char *argv[])
{
    char c;
    //load section
#if INPUT_MODE ==0
    printf("Is Continued? : Y or N");
    scanf("%c",&c);
#elif INPUT_MODE == 1
    c = *argv[1];
#endif
    if(c=='Y')
    {
#if INPUT_MODE == 0
        printf("Input file name : ");
        scanf("%s",stateFileName);
#elif INPUT_MODE == 1
        strcpy(stateFileName, argv[2]);
#endif
        INPUT_PARAMETER(argv);
        LOAD_CURRENT_STATE(stateFileName, molecule, &nParticle, &mcTime, &epsylon);
        sprintf(energyProfileName, "%s.csv",stateFileName);
        out2 = fopen(energyProfileName,"a");
        fclose(out2);
        for(int i=0;i<nParticle;i++) SET_SEGMENT_COORDINATE(molecule[i].coordinate[0], i, system_lattice_3D);
        INITIALIZE_PROGRAM();
        
        int sum_LJ = 0;
        for(int i=0;i<nParticle;i++)
        {
            sum_LJ += Calculate_eLJ_ON_LATTICE(molecule, i, 0, system_lattice_3D);
        }
        printf("%lf",(double)sum_LJ * (epsylon/2.0));
        
        SAVE_CURRENT_STATE(stateFileName, molecule, nParticle, mcTime, epsylon);
        PDB_File_Write(trajName, 0, molecule, nParticle, mcTime/mcInterval+1);
    }
    //new simulation section
    else if(c=='N')
    {
#if INPUT_MODE ==0
        printf("Input state save file name : ");
        scanf("%s", stateFileName);
#elif INPUT_MODE == 1
        strcpy(stateFileName, argv[2]);
#endif
        INITIALIZE_PROGRAM();
        INPUT_PARAMETER(argv);

        SET_INITIAL_COORDINATE(argv);
        mcTime = 0;
        sprintf(energyProfileName, "%s.csv",stateFileName);
        out2 = fopen(energyProfileName,"w");
        fclose(out2);
        
        SAVE_CURRENT_STATE(stateFileName, molecule, nParticle, mcTime, epsylon);
    }
    SET_EPSYLON_SET(fp_chiN, &Object_MCS, &epsylon, mcTime);
    while(mcTime < mcTimeMAX)
    {
        for(int submcTime = 0; submcTime < mcInterval; submcTime++) MetropolisMonteCarloEngine();
        mcTime += mcInterval;
        int sum_LJ = 0;
        for(int i=0;i<nParticle;i++)
        {
            sum_LJ += Calculate_eLJ_ON_LATTICE(molecule, i, 0, system_lattice_3D);
        }
        out2 = fopen(energyProfileName,"a");
        fprintf(out2, "%d,%lf,%lf\n", mcTime, epsylon, (double)sum_LJ * (epsylon/2.0));
        fclose(out2);
        //energyprofile
        SET_EPSYLON_SET(fp_chiN, &Object_MCS, &epsylon, mcTime);
        //PDB_File_Write(trajName, 0, molecule, nParticle, mcTime/mcInterval+1);
        SAVE_CURRENT_STATE(stateFileName, molecule, nParticle, mcTime, epsylon);
    }
    printf("End\n");
}

void INPUT_PARAMETER(char *argv[])
{
    //system lattice initialization
    for(int i=0;i<BOUNDARY_MAX;i++)
        for(int j=0;j<BOUNDARY_MAX;j++)
            for(int k=0;k<BOUNDARY_Z;k++)
                system_lattice_3D[i][j][k] = -1;
    //Input Variable
#if INPUT_MODE == 0
    printf("Input Total Simulation Time(MCS) : ");
    scanf("%d", &mcTimeMAX);
    printf("Input time intervals for recording property : ");
    scanf("%d", &mcInterval);
#elif INPUT_MODE == 1
    mcTimeMAX = atoi(argv[3]);
    mcInterval = atoi(argv[4]);
#endif
    strcpy(trajName, stateFileName);
    strcat(trajName, ".pdb");
    strcpy(datFileName, stateFileName);
    strcat(datFileName, ".dat");
    fp_chiN = fopen("chiN_profile.dat", "r");
}

//initial coordinate complex arm
void SET_INITIAL_COORDINATE(char *argv[])
{
    //set your molecule...
    PDB_File_Write(trajName, 1, molecule, nParticle, 0);
}

void MetropolisMonteCarloEngine()
{
    for(mcTrial = 0; mcTrial < nParticle; mcTrial++)
    {
        //select segment randomly
        int segment_num = (int)(WELLRNG512a()*nParticle);
        
        //random hopping
        Hopping_ON_LATTICE(molecule[segment_num].coordinate);
        
        //set empty in system box about original segment position
        SET_SEGMENT_COORDINATE(molecule[segment_num].coordinate[0], -1, system_lattice_3D);
        
        //is hopping position not occupied?
        if(!IS_OCCUPIED_LATTICE(molecule[segment_num].coordinate[1], system_lattice_3D))
        {
            //is bond lenth possible?
            if(IS_POSSIBLE_BOND_VECTOR(molecule, segment_num))
            {
                if(epsylon !=0)
                {
                    //calcualte difference eVDW
                    int sum_LJ = Calculate_eLJ_ON_LATTICE(molecule, segment_num, 1, system_lattice_3D) - Calculate_eLJ_ON_LATTICE(molecule, segment_num, 0, system_lattice_3D);
                    //calculate possiblity
                    if(sum_LJ <= 0)
                    {
                        molecule[segment_num].coordinate[0][0] = molecule[segment_num].coordinate[1][0];
                        molecule[segment_num].coordinate[0][1] = molecule[segment_num].coordinate[1][1];
                        molecule[segment_num].coordinate[0][2] = molecule[segment_num].coordinate[1][2];
                    }
                    else
                    {
                        double possibility = exp(-((double)sum_LJ)*epsylon);
                        //accept?
                        if(possibility>=WELLRNG512a())
                        {
                            molecule[segment_num].coordinate[0][0] = molecule[segment_num].coordinate[1][0];
                            molecule[segment_num].coordinate[0][1] = molecule[segment_num].coordinate[1][1];
                            molecule[segment_num].coordinate[0][2] = molecule[segment_num].coordinate[1][2];
                            
                        }
                    }
                }
                //if epsylon 0, non check energy
                else
                {
                    molecule[segment_num].coordinate[0][0] = molecule[segment_num].coordinate[1][0];
                    molecule[segment_num].coordinate[0][1] = molecule[segment_num].coordinate[1][1];
                    molecule[segment_num].coordinate[0][2] = molecule[segment_num].coordinate[1][2];
                }
            }
        }
        SET_SEGMENT_COORDINATE(molecule[segment_num].coordinate[0], segment_num, system_lattice_3D);
    }
}
