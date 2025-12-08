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

void bond_finder_grid(
    double *coords, int N,
    double tol,
    int **neighbor_list,     // pre-allocated array of int* of length N
    int *neighbor_count      // pre-allocated array of int of length N
) {
    double tol2 = tol*tol;

    // Compute bounding box
    double xmin=1e300,xmax=-1e300, ymin=1e300,ymax=-1e300, zmin=1e300,zmax=-1e300;
    for(int i=0;i<N;i++){
        double x=coords[3*i], y=coords[3*i+1], z=coords[3*i+2];
        if(x<xmin) xmin=x; if(x>xmax) xmax=x;
        if(y<ymin) ymin=y; if(y>ymax) ymax=y;
        if(z<zmin) zmin=z; if(z>zmax) zmax=z;
    }

    double cell_size = tol;
    int nx = (int)((xmax-xmin)/cell_size)+1;
    int ny = (int)((ymax-ymin)/cell_size)+1;
    int nz = (int)((zmax-zmin)/cell_size)+1;
    int total_cells = nx*ny*nz;

    // Initialize grid
    Cell *grid = (Cell*)calloc(total_cells, sizeof(Cell));

    // First pass: count occupancy
    for(int i=0;i<N;i++){
        int ix = (int)((coords[3*i] - xmin)/cell_size);
        int iy = (int)((coords[3*i+1] - ymin)/cell_size);
        int iz = (int)((coords[3*i+2] - zmin)/cell_size);
        int h = cell_hash(ix,iy,iz,nx,ny,nz);
        grid[h].count++;
    }

    // Allocate cell indices
    for(int i=0;i<total_cells;i++){
        if(grid[i].count>0){
            grid[i].indices = (int*)malloc(grid[i].count*sizeof(int));
            grid[i].count = 0;
            grid[i].capacity = 0; // not used further
        }
    }

    // Fill grid
    for(int i=0;i<N;i++){
        int ix = (int)((coords[3*i] - xmin)/cell_size);
        int iy = (int)((coords[3*i+1] - ymin)/cell_size);
        int iz = (int)((coords[3*i+2] - zmin)/cell_size);
        int h = cell_hash(ix,iy,iz,nx,ny,nz);
        grid[h].indices[ grid[h].count++ ] = i;
    }

    // Neighbor offsets
    int off[27][3], c=0;
    for(int dx=-1;dx<=1;dx++)
        for(int dy=-1;dy<=1;dy++)
            for(int dz=-1;dz<=1;dz++)
                { off[c][0]=dx; off[c][1]=dy; off[c][2]=dz; c++; }

    // Initialize neighbor counts
    for(int i=0;i<N;i++) neighbor_count[i] = 0;

    // Check neighbors using cells
    for(int i=0;i<N;i++){
        double xi=coords[3*i], yi=coords[3*i+1], zi=coords[3*i+2];
        int ix = (int)((xi - xmin)/cell_size);
        int iy = (int)((yi - ymin)/cell_size);
        int iz = (int)((zi - zmin)/cell_size);

        for(int o=0;o<27;o++){
            int nx1=ix+off[o][0], ny1=iy+off[o][1], nz1=iz+off[o][2];
            if(nx1<0||ny1<0||nz1<0||nx1>=nx||ny1>=ny||nz1>=nz) continue;
            int h = cell_hash(nx1,ny1,nz1,nx,ny,nz);
            Cell *cell = &grid[h];

            for(int k=0;k<cell->count;k++){
                int j = cell->indices[k];
                if(j <= i) continue; // only upper triangle
                double dx=coords[3*j]-xi, dy=coords[3*j+1]-yi, dz=coords[3*j+2]-zi;
                double d2 = dx*dx + dy*dy + dz*dz;
                if(d2 <= tol2){
                    // Add neighbor for i
                    neighbor_list[i][ neighbor_count[i]++ ] = j;
                    // Add neighbor for j
                    neighbor_list[j][ neighbor_count[j]++ ] = i;
                }
            }
        }
    }

    // Free grid memory
    for(int i=0;i<total_cells;i++) if(grid[i].indices) free(grid[i].indices);
    free(grid);
}

