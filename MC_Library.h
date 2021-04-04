//header section
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "WELL512a.h"
#include <math.h>
#include <string.h>
#include <stdbool.h>

#define BOUNDARY_MAX 100 //system box size
#define BOUNDARY_Z 30
#define INPUT_MODE 1 //if 0 : console input mode, 1 : batch file input mode

//SEGMENT data structure
typedef struct SEGMENT_STRUCTURE
{
    int coordinate[2][3]; //[0 : old / 1 : new][x,y,z]
    int linked_segment[6]; //linked segment list
    int linked_segment_num; //for speed, save linked segment list number
    int segment_type; //segment charater type
}SEGMENT;

//function section

//math tool function
void INITIALIZE_PROGRAM(); //initialize program(set seed for WELL algorithm)

//file io function
void PDB_File_Write(char [], bool, SEGMENT *, int, int time); //output pdb file
void SAVE_CURRENT_STATE(char filename[50], SEGMENT *molecule, int nParticle, unsigned int mcTime, double epsylon); //output current state segment function
void LOAD_CURRENT_STATE(char filename[50], SEGMENT *molecule, int *nParticle, unsigned int *mcTime, double *epsylon); //input written state segment function

//segment move function
void Hopping_ON_LATTICE(int (*coordinate)[3]); //random hopping for on lattice

int Calculate_eLJ_ON_LATTICE(SEGMENT *, unsigned int segment, int is_new_segment, int (*system_lattice_3D)[BOUNDARY_MAX][BOUNDARY_Z]); //get LJ potential for on lattice

//for on lattice function
void SET_SEGMENT_COORDINATE(int coordinate[3], int value, int (*system_lattice_3D)[BOUNDARY_MAX][BOUNDARY_Z]); //set segment coordinate in system box
bool IS_OCCUPIED_LATTICE(int coordinate[3], int (*system_lattice_3D)[BOUNDARY_MAX][BOUNDARY_Z]); //is that lattice site is occupied? return bool value
bool IS_POSSIBLE_BOND_VECTOR(SEGMENT *molecule, int segment_num); //is bond vector is possible length? return bool value
void SET_EPSYLON_SET(FILE *fp_Energy_Profile, int *Object_MCS, double *epsylon, int mcTime);
























