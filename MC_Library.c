#include "MC_Library.h"
#include <unistd.h>

int Half_Boundary_Max = BOUNDARY_MAX / 2;
int Half_Boundary_Z = BOUNDARY_Z / 2;

const int bond_possible_list[54][3] = {
    {-2, -1, -1},{-2, -1, 0},{-2, -1, 1},{-2, 0, -1},{-2, 0, 0},{-2, 0, 1},{-2, 1, -1},{-2, 1, 0},{-2, 1, 1},{-1, -2, -1},{-1, -2, 0},{-1, -2, 1},{-1, -1, -2},{-1, -1, 2},{-1, 0, -2},{-1, 0, 2},{-1, 1, -2},{-1, 1, 2},{-1, 2, -1},{-1, 2, 0},{-1, 2, 1},{0, -2, -1},{0, -2, 0},{0, -2, 1},{0, -1, -2},{0, -1, 2},{0, 0, -2},{0, 0, 2},{0, 1, -2},{0, 1, 2},{0, 2, -1},{0, 2, 0},{0, 2, 1},{1, -2, -1},{1, -2, 0},{1, -2, 1},{1, -1, -2},{1, -1, 2},{1, 0, -2},{1, 0, 2},{1, 1, -2},{1, 1, 2},{1, 2, -1},{1, 2, 0},{1, 2, 1},{2, -1, -1},{2, -1, 0},{2, -1, 1},{2, 0, -1},{2, 0, 0},{2, 0, 1},{2, 1, -1},{2, 1, 0},{2, 1, 1}
};

//set seed for random generator
void INITIALIZE_PROGRAM()
{
    //initialize WELL512(random generator) algorithm
    unsigned int seed[16];
    srand((unsigned)time(NULL));
    for(int i=0; i<16; i++) seed[i] = rand() + getpid();
    InitWELLRNG512a(seed);
}

//output segment info to pdb format
void PDB_File_Write(char file_name[50], bool is_new_file, SEGMENT *molecule, int nParticle, int time)
{
    FILE *fp_pdb;
    if(is_new_file)
    {
        fp_pdb = fopen(file_name, "w");
    }
    else fp_pdb = fopen(file_name, "a");
    fprintf(fp_pdb, "\n\nMODEL%9d\n", time);
    //print HETATM section
    for(int i=0; i<nParticle; i++)
    {
        fprintf(fp_pdb, "HETATM%5d %4s%c%3s %c%4d    %8d%8d%8d%25s\n",i+1, "C", 'A'+molecule[i].segment_type, "CAR", 'A'+molecule[i].segment_type, molecule[i].segment_type, molecule[i].coordinate[0][0], molecule[i].coordinate[0][1], molecule[i].coordinate[0][2]," ");
    }
    //print CONECT section
    for(int i=0; i<nParticle; i++)
    {
        int check = 0;
        for(int j=0; j<molecule[i].linked_segment_num; j++)
        {
            if(i<molecule[i].linked_segment[j])
            {
                if(!check)
                {
                    fprintf(fp_pdb, "CONECT%5d%5d", i+1, molecule[i].linked_segment[j]+1);
                    check=1;
                }
                else fprintf(fp_pdb, "%5d", molecule[i].linked_segment[j]+1);
            }
        }
        fprintf(fp_pdb, "\n");
    }
    //print ENDMDL
    fprintf(fp_pdb, "ENDMDL\n\n");
    fclose(fp_pdb);
}

//calculate LJ potentail in Onlattice(segment)
int Calculate_eLJ_ON_LATTICE(SEGMENT *molecule, unsigned int segment, int is_new_segment, int (*system_lattice_3D)[BOUNDARY_MAX][BOUNDARY_Z])
{
    int center_coordinate[3] = {molecule[segment].coordinate[is_new_segment][0], molecule[segment].coordinate[is_new_segment][1], molecule[segment].coordinate[is_new_segment][2]};
    int sum_LJ = 0;
    for(int i=0; i<54; i++)
    {
        //find near segment
        int
        x = center_coordinate[0] + bond_possible_list[i][0],
        y = center_coordinate[1] + bond_possible_list[i][1],
        z = center_coordinate[2] + bond_possible_list[i][2];
        
        //apply periodic boundary
        if(x>=BOUNDARY_MAX) x-=BOUNDARY_MAX;
        else if(x<0) x+=BOUNDARY_MAX;
        if(y>=BOUNDARY_MAX) y-=BOUNDARY_MAX;
        else if(y<0) y+=BOUNDARY_MAX;
        if(z>=BOUNDARY_Z) z-=BOUNDARY_Z;
        else if(z<0) z+=BOUNDARY_Z;
        
        //if near segment found?
        if(system_lattice_3D[x][y][z]>=0)
            if(molecule[segment].segment_type != molecule[system_lattice_3D[x][y][z]].segment_type) sum_LJ++;
    }
    
    //in short range, deduct connected segments
    for(int i=0; i < molecule[segment].linked_segment_num; i++)
    {
        //get bond vector
        SEGMENT * TEMP = &molecule[molecule[segment].linked_segment[i]];
        int bond_vector[3] =
        {
            center_coordinate[0] - TEMP->coordinate[0][0],
            center_coordinate[1] - TEMP->coordinate[0][1],
            center_coordinate[2] - TEMP->coordinate[0][2]
        };
        
        //periodic boundary revision
        if(Half_Boundary_Max<bond_vector[0]) bond_vector[0] -= BOUNDARY_MAX;
        else if((-Half_Boundary_Max)>bond_vector[0]) bond_vector[0] += BOUNDARY_MAX;
        if(Half_Boundary_Max<bond_vector[1]) bond_vector[1] -= BOUNDARY_MAX;
        else if((-Half_Boundary_Max)>bond_vector[1]) bond_vector[1] += BOUNDARY_MAX;
        if(Half_Boundary_Z<bond_vector[2]) bond_vector[2] -= BOUNDARY_Z;
        else if((-Half_Boundary_Z)>bond_vector[2]) bond_vector[2] += BOUNDARY_Z;
        
        int bond_length = bond_vector[0]*bond_vector[0] + bond_vector[1]*bond_vector[1] + bond_vector[2]*bond_vector[2];
        //energy count bond length^2 in SW : 4, 5, 6
        if(bond_length <= 6)
        {
            if(molecule[segment].segment_type != TEMP->segment_type) sum_LJ--;
        }
    }

    return sum_LJ;
}

//random hopping in on lattice
void Hopping_ON_LATTICE(int (*coordinate)[3])
{
    //move 1 way in 6way(+-1 in coordinate)
    unsigned int movable_decision_num = (int)(WELLRNG512a()*6);
    coordinate[1][0] = coordinate[0][0];
    coordinate[1][1] = coordinate[0][1];
    coordinate[1][2] = coordinate[0][2];
    switch(movable_decision_num)
    {
        case 0:
            coordinate[1][0]++;
            if(coordinate[1][0] >= BOUNDARY_MAX) coordinate[1][0] -= BOUNDARY_MAX;
            break;
        case 1:
            coordinate[1][0]--;
            if(coordinate[1][0] < 0) coordinate[1][0] += BOUNDARY_MAX;
            break;
        case 2:
            coordinate[1][1]++;
            if(coordinate[1][1] >= BOUNDARY_MAX) coordinate[1][1] -= BOUNDARY_MAX;
            break;
        case 3:
            coordinate[1][1]--;
            if(coordinate[1][1] < 0) coordinate[1][1] += BOUNDARY_MAX;
            break;
        case 4:
            coordinate[1][2]++;
            if(coordinate[1][2] >= BOUNDARY_Z) coordinate[1][2] -= BOUNDARY_Z;
            break;
        case 5:
            coordinate[1][2]--;
            if(coordinate[1][2] < 0) coordinate[1][2] += BOUNDARY_Z;
            break;
    }
}

//set segment in system box
void SET_SEGMENT_COORDINATE(int coordinate[3], int value, int (*system_lattice_3D)[BOUNDARY_MAX][BOUNDARY_Z])
{
    int other_value;
    if(value == -1) other_value = -1;
    else other_value = -5;
    int x = coordinate[0], y = coordinate[1], z = coordinate[2];
    int x_1 = coordinate[0]+1, y_1 = coordinate[1]+1, z_1 = coordinate[2]+1;
    if(x_1 >= BOUNDARY_MAX) x_1 -= BOUNDARY_MAX;
    if(y_1 >= BOUNDARY_MAX) y_1 -= BOUNDARY_MAX;
    if(z_1 >= BOUNDARY_Z) z_1 -= BOUNDARY_Z;
    
    //if empty, set -1
    //if occupy, center value is segment value, other value is -5
    system_lattice_3D[x][y][z] = value;
    system_lattice_3D[x][y][z_1] = other_value;
    system_lattice_3D[x][y_1][z] = other_value;
    system_lattice_3D[x][y_1][z_1] = other_value;
    system_lattice_3D[x_1][y][z] = other_value;
    system_lattice_3D[x_1][y][z_1] = other_value;
    system_lattice_3D[x_1][y_1][z] = other_value;
    system_lattice_3D[x_1][y_1][z_1] = other_value;

}

bool IS_OCCUPIED_LATTICE(int coordinate[3], int (*system_lattice_3D)[BOUNDARY_MAX][BOUNDARY_Z])
{
    int x = coordinate[0], y = coordinate[1], z = coordinate[2];
    int x_1 = coordinate[0]+1, y_1 = coordinate[1]+1, z_1 = coordinate[2]+1;
    if(x_1 >= BOUNDARY_MAX) x_1 -= BOUNDARY_MAX;
    if(y_1 >= BOUNDARY_MAX) y_1 -= BOUNDARY_MAX;
    if(z_1 >= BOUNDARY_Z) z_1 -= BOUNDARY_Z;
    
    if((system_lattice_3D[x][y][z]!=-1) |
       (system_lattice_3D[x][y][z_1]!=-1) |
       (system_lattice_3D[x][y_1][z]!=-1) |
       (system_lattice_3D[x][y_1][z_1]!=-1) |
       (system_lattice_3D[x_1][y][z]!=-1) |
       (system_lattice_3D[x_1][y][z_1]!=-1) |
       (system_lattice_3D[x_1][y_1][z]!=-1) |
       (system_lattice_3D[x_1][y_1][z_1]!=-1))
        return true;
    
    else return false;
}

//check possible bond vector if possible, return true
bool IS_POSSIBLE_BOND_VECTOR(SEGMENT *molecule, int segment_num)
{
    int criteria_coordinate[3] =
    {
        molecule[segment_num].coordinate[1][0],
        molecule[segment_num].coordinate[1][1],
        molecule[segment_num].coordinate[1][2]
    };
    
    int link_segment_num = molecule[segment_num].linked_segment_num;
    for(int i=0; i<link_segment_num; i++)
    {
        //get bond vector
        SEGMENT * TEMP = &molecule[molecule[segment_num].linked_segment[i]];
        int bond_vector[3] =
        {
            criteria_coordinate[0] - TEMP->coordinate[0][0],
            criteria_coordinate[1] - TEMP->coordinate[0][1],
            criteria_coordinate[2] - TEMP->coordinate[0][2]
        };
        
        //periodic boundary revision
        if(Half_Boundary_Max<bond_vector[0]) bond_vector[0] -= BOUNDARY_MAX;
        else if((-Half_Boundary_Max)>bond_vector[0]) bond_vector[0] += BOUNDARY_MAX;
        if(Half_Boundary_Max<bond_vector[1]) bond_vector[1] -= BOUNDARY_MAX;
        else if((-Half_Boundary_Max)>bond_vector[1]) bond_vector[1] += BOUNDARY_MAX;
        if(Half_Boundary_Z<bond_vector[2]) bond_vector[2] -= BOUNDARY_Z;
        else if((-Half_Boundary_Z)>bond_vector[2]) bond_vector[2] += BOUNDARY_Z;
        
        int bond_length = bond_vector[0]*bond_vector[0] + bond_vector[1]*bond_vector[1] + bond_vector[2]*bond_vector[2];
        //possible bond length = 4, 5, 6, 9, 10
        if(bond_length>10) return false;
        if(bond_length == 8) return false;
    }
    return true;
}

//save current segment state
void SAVE_CURRENT_STATE(char filename[50], SEGMENT *molecule, int nParticle, unsigned int mcTime, double epsylon)
{
    FILE *fp;
    char new_file_name[100];
    sprintf(new_file_name, "%s_MCS_%08d", filename, mcTime);
    fp = fopen(new_file_name, "w");
    //save particle number, qfactor, epsylon, mctime
    fprintf(fp,"%d %g %u ",nParticle, epsylon, mcTime);
    for(int i=0; i<nParticle; i++)
    {
        //save segment data
        SEGMENT *TEMP = &molecule[i];
        fprintf(fp, "%d %d %d %d %d ", TEMP->coordinate[0][0], TEMP->coordinate[0][1], TEMP->coordinate[0][2], TEMP->linked_segment_num, TEMP->segment_type);
        for(int j=0; j<TEMP->linked_segment_num; j++) fprintf(fp, "%d ", TEMP->linked_segment[j]);
        fprintf(fp,"\n");
    }
    fclose(fp);
}

//load saved data
void LOAD_CURRENT_STATE(char filename[50], SEGMENT *molecule, int *nParticle, unsigned int *mcTime, double *epsylon)
{
    FILE *fp;
    fp = fopen (filename, "r");
    fscanf(fp, "%d %lf %u ", nParticle, epsylon, mcTime);
    int local_nParticle = *nParticle;
    for(int i=0; i<local_nParticle; i++)
    {
        SEGMENT *TEMP = &molecule[i];
        fscanf(fp, "%d %d %d %d %d ", &(TEMP->coordinate[0][0]), &(TEMP->coordinate[0][1]), &(TEMP->coordinate[0][2]), &(TEMP->linked_segment_num), &(TEMP->segment_type));
        for(int j=0; j<TEMP->linked_segment_num; j++) fscanf(fp, "%d ", &(TEMP->linked_segment[j]));
    }
    fclose(fp);
}

void SET_EPSYLON_SET(FILE *fp_Energy_Profile, int *Object_MCS, double *epsylon, int mcTime)
{
    while(*Object_MCS < mcTime)
    {
        double temp_epsylon;
        fscanf(fp_Energy_Profile, "%d %lf", Object_MCS, &temp_epsylon);
        if(temp_epsylon >= 0) *epsylon = temp_epsylon;
    }
}

















