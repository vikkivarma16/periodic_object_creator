#include <math.h>
#include <stdlib.h>
#include <stdint.h>

typedef struct {
    int count;
    int capacity;
    int *indices;
} Cell;

static inline int cell_hash(int ix, int iy, int iz, int nx, int ny, int nz) {
    return ix + nx * (iy + ny * iz);
}

/*
  Grid-based bond finder
  Arguments:
    obj: positions (x,y,z) array
    N: number of atoms
    tol: bond length tolerance
    body_ids: int array of body IDs
    neighbor_counts: int array of size N, counts of bonded neighbors
    neighbor_list: int** array of size N, each pointer points to array of bonded neighbor indices
*/
void bond_finder_grid(
    double *obj, int N,
    double tol,
    int *body_ids,
    int *neighbor_counts,
    int **neighbor_list
){
    double tol2 = tol*tol;

    // Compute bounds
    double xmin=1e300, xmax=-1e300;
    double ymin=1e300, ymax=-1e300;
    double zmin=1e300, zmax=-1e300;
    for(int i=0;i<N;i++){
        double x=obj[3*i], y=obj[3*i+1], z=obj[3*i+2];
        if(x<xmin) xmin=x; if(x>xmax) xmax=x;
        if(y<ymin) ymin=y; if(y>ymax) ymax=y;
        if(z<zmin) zmin=z; if(z>zmax) zmax=z;
    }

    int nx = (int)((xmax-xmin)/tol)+1;
    int ny = (int)((ymax-ymin)/tol)+1;
    int nz = (int)((zmax-zmin)/tol)+1;
    int total_cells = nx*ny*nz;

    Cell *grid = (Cell*)calloc(total_cells, sizeof(Cell));

    // Count occupancy
    for(int i=0;i<N;i++){
        int ix=(int)((obj[3*i]-xmin)/tol);
        int iy=(int)((obj[3*i+1]-ymin)/tol);
        int iz=(int)((obj[3*i+2]-zmin)/tol);
        int h = cell_hash(ix,iy,iz,nx,ny,nz);
        grid[h].capacity++;
    }

    // Allocate
    for(int i=0;i<total_cells;i++){
        if(grid[i].capacity>0){
            grid[i].indices = (int*)malloc(grid[i].capacity*sizeof(int));
            grid[i].count=0;
        }
    }

    // Fill grid
    for(int i=0;i<N;i++){
        int ix=(int)((obj[3*i]-xmin)/tol);
        int iy=(int)((obj[3*i+1]-ymin)/tol);
        int iz=(int)((obj[3*i+2]-zmin)/tol);
        int h = cell_hash(ix,iy,iz,nx,ny,nz);
        grid[h].indices[grid[h].count++] = i;
    }

    // Neighbor offsets (27 cells)
    int off[27][3];
    int c=0;
    for(int dx=-1;dx<=1;dx++)
        for(int dy=-1;dy<=1;dy++)
            for(int dz=-1;dz<=1;dz++){
                off[c][0]=dx; off[c][1]=dy; off[c][2]=dz; c++;
            }

    // Initialize neighbor counts
    for(int i=0;i<N;i++) neighbor_counts[i]=0;

    // Compute neighbors
    for(int i=0;i<N;i++){
        double xi=obj[3*i], yi=obj[3*i+1], zi=obj[3*i+2];
        int bi=body_ids[i];

        int ix=(int)((xi-xmin)/tol);
        int iy=(int)((yi-ymin)/tol);
        int iz=(int)((zi-zmin)/tol);

        for(int o=0;o<27;o++){
            int nx1=ix+off[o][0], ny1=iy+off[o][1], nz1=iz+off[o][2];
            if(nx1<0||ny1<0||nz1<0||nx1>=nx||ny1>=ny||nz1>=nz) continue;
            int h=cell_hash(nx1,ny1,nz1,nx,ny,nz);
            Cell *cell = &grid[h];

            for(int k=0;k<cell->count;k++){
                int j=cell->indices[k];
                if(j<=i) continue;
                if(body_ids[j]!=bi) continue;

                double dx=obj[3*j]-xi, dy=obj[3*j+1]-yi, dz=obj[3*j+2]-zi;
                double d2=dx*dx+dy*dy+dz*dz;
                if(d2 <= tol2){
                    neighbor_list[i][neighbor_counts[i]++] = j;
                    neighbor_list[j][neighbor_counts[j]++] = i; // symmetric
                }
            }
        }
    }

    // Free grid memory
    for(int i=0;i<total_cells;i++){
        if(grid[i].indices) free(grid[i].indices);
    }
    free(grid);
}

