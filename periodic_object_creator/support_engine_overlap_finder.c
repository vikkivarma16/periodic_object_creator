#include <math.h>
#include <stdlib.h>
#include <stdint.h>

/*
  Fast overlap finder using spatial hashing grid.
  Inspired by original C code given by user.

  Arguments:
    obj1, N1:   first object points (x,y,z)
    obj2, N2:   second object points
    delete_from: 0 = delete from obj1, 1 = delete from obj2
    tol:        tolerance distance

  Output:
    delete_mask array of size N1 or N2 filled with 0/1
*/

typedef struct {
    int count;
    int *indices;
} Cell;

static inline int cell_hash(int ix, int iy, int iz, int nx, int ny, int nz) {
    return ix + nx * (iy + ny * iz);
}

void overlap_eliminator(
    double *obj1, int N1,
    double *obj2, int N2,
    int delete_from,
    double tol,
    int *delete_mask
){
    int N = delete_from == 0 ? N1 : N2;
    double *A = delete_from == 0 ? obj1 : obj2;
    double *B = delete_from == 0 ? obj2 : obj1;
    int NA = delete_from == 0 ? N1 : N2;
    int NB = delete_from == 0 ? N2 : N1;

    double cell_size = tol;
    double tol2 = tol * tol;

    double xmin=1e300,xmax=-1e300, ymin=1e300,ymax=-1e300, zmin=1e300,zmax=-1e300;

    for (int i=0;i<NB;i++){
        double x=B[3*i], y=B[3*i+1], z=B[3*i+2];
        if (x<xmin) xmin=x;
        if (x>xmax) xmax=x;
        if (y<ymin) ymin=y;
        if (y>ymax) ymax=y;
        if (z<zmin) zmin=z;
        if (z>zmax) zmax=z;
    }

    int nx = (int)((xmax-xmin)/cell_size)+1;
    int ny = (int)((ymax-ymin)/cell_size)+1;
    int nz = (int)((zmax-zmin)/cell_size)+1;

    int total_cells = nx*ny*nz;

    Cell *grid = (Cell*)calloc(total_cells, sizeof(Cell));
    for(int i=0;i<total_cells;i++){
        grid[i].count = 0;
        grid[i].indices = NULL;
    }

    // Count occupancy
    for(int i=0;i<NB;i++){
        int ix = (int)((B[3*i]   - xmin)/cell_size);
        int iy = (int)((B[3*i+1] - ymin)/cell_size);
        int iz = (int)((B[3*i+2] - zmin)/cell_size);

        int h = cell_hash(ix,iy,iz,nx,ny,nz);
        grid[h].count++;
    }

    // Allocate memory
    for(int i=0;i<total_cells;i++){
        if(grid[i].count>0){
            grid[i].indices = (int*)malloc(grid[i].count*sizeof(int));
            grid[i].count = 0;
        }
    }

    // Fill the grid
    for(int i=0;i<NB;i++){
        int ix = (int)((B[3*i]   - xmin)/cell_size);
        int iy = (int)((B[3*i+1] - ymin)/cell_size);
        int iz = (int)((B[3*i+2] - zmin)/cell_size);

        int h = cell_hash(ix,iy,iz,nx,ny,nz);
        grid[h].indices[ grid[h].count++ ] = i;
    }

    // Clear delete mask
    for(int i=0;i<NA;i++) delete_mask[i] = 0;

    // Neighbor offsets
    int off[27][3];
    int c=0;
    for(int dx=-1;dx<=1;dx++)
        for(int dy=-1;dy<=1;dy++)
            for(int dz=-1;dz<=1;dz++){
                off[c][0]=dx;
                off[c][1]=dy;
                off[c][2]=dz;
                c++;
            }

    // Check overlaps
    for(int i=0;i<NA;i++){
        double x=A[3*i], y=A[3*i+1], z=A[3*i+2];

        int ix = (int)((x - xmin)/cell_size);
        int iy = (int)((y - ymin)/cell_size);
        int iz = (int)((z - zmin)/cell_size);

        for(int o=0;o<27;o++){
            int nx1=ix+off[o][0];
            int ny1=iy+off[o][1];
            int nz1=iz+off[o][2];

            if(nx1<0 || ny1<0 || nz1<0 || nx1>=nx || ny1>=ny || nz1>=nz)
                continue;

            int h = cell_hash(nx1,ny1,nz1,nx,ny,nz);
            Cell *cell = &grid[h];

            for(int k=0;k<cell->count;k++){
                int j = cell->indices[k];

                double dx = B[3*j]   - x;
                double dy = B[3*j+1] - y;
                double dz = B[3*j+2] - z;

                if(dx*dx + dy*dy + dz*dz < tol2){
                    delete_mask[i] = 1;
                    goto NEXT_POINT;
                }
            }
        }
        NEXT_POINT:;
    }

    // Free grid memory
    for(int i=0;i<total_cells;i++){
        if(grid[i].indices) free(grid[i].indices);
    }
    free(grid);
}

