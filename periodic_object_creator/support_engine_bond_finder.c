#include <math.h>
#include <stdlib.h>
#include <stdint.h>

/*
  Bond finder using spatial hashing grid.
  Arguments:
    obj: positions array (x,y,z)
    N: number of atoms
    tol: bond length tolerance
    neighbor_list: output array of size N*N (flattened) where 1 = bonded
*/

typedef struct {
    int count;
    int *indices;
} Cell;

static inline int cell_hash(int ix, int iy, int iz, int nx, int ny, int nz) {
    return ix + nx * (iy + ny * iz);
}

void bond_finder_grid(
    double *obj, int N,
    double tol,
    int *neighbor_list   // flattened N x N, 0/1
) {
    double tol2 = tol*tol;

    for(int i=0;i<N*N;i++) neighbor_list[i] = 0;

    // Compute bounding box
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
    for(int i=0;i<total_cells;i++){
        grid[i].count = 0;
        grid[i].indices = NULL;
    }

    // Count occupancy
    for(int i=0;i<N;i++){
        int ix = (int)((obj[3*i]   - xmin)/tol);
        int iy = (int)((obj[3*i+1] - ymin)/tol);
        int iz = (int)((obj[3*i+2] - zmin)/tol);

        int h = cell_hash(ix, iy, iz, nx, ny, nz);
        grid[h].count++;
    }

    // Allocate indices arrays
    for(int i=0;i<total_cells;i++){
        if(grid[i].count>0){
            grid[i].indices = (int*)malloc(grid[i].count*sizeof(int));
            grid[i].count = 0; // reuse as fill counter
        }
    }

    // Fill grid
    for(int i=0;i<N;i++){
        int ix = (int)((obj[3*i]   - xmin)/tol);
        int iy = (int)((obj[3*i+1] - ymin)/tol);
        int iz = (int)((obj[3*i+2] - zmin)/tol);

        int h = cell_hash(ix, iy, iz, nx, ny, nz);
        grid[h].indices[grid[h].count++] = i;
    }

    // Neighbor offsets (27 cells)
    int off[27][3];
    int c=0;
    for(int dx=-1; dx<=1; dx++)
        for(int dy=-1; dy<=1; dy++)
            for(int dz=-1; dz<=1; dz++){
                off[c][0]=dx; off[c][1]=dy; off[c][2]=dz;
                c++;
            }

    // Loop over atoms
    for(int i=0;i<N;i++){
        double xi=obj[3*i], yi=obj[3*i+1], zi=obj[3*i+2];
        int ix = (int)((xi - xmin)/tol);
        int iy = (int)((yi - ymin)/tol);
        int iz = (int)((zi - zmin)/tol);

        for(int o=0;o<27;o++){
            int nx1 = ix + off[o][0];
            int ny1 = iy + off[o][1];
            int nz1 = iz + off[o][2];
            if(nx1<0 || ny1<0 || nz1<0 || nx1>=nx || ny1>=ny || nz1>=nz)
                continue;
            int h = cell_hash(nx1, ny1, nz1, nx, ny, nz);
            Cell *cell = &grid[h];
            for(int k=0;k<cell->count;k++){
                int j = cell->indices[k];
                if(j <= i) continue; // avoid double count
                double dx = obj[3*j] - xi;
                double dy = obj[3*j+1] - yi;
                double dz = obj[3*j+2] - zi;
                if(dx*dx + dy*dy + dz*dz <= tol2){
                    neighbor_list[i*N+j] = 1;
                    neighbor_list[j*N+i] = 1;
                }
            }
        }
    }

    // Free memory
    for(int i=0;i<total_cells;i++)
        if(grid[i].indices) free(grid[i].indices);
    free(grid);
}

