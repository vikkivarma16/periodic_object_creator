#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <stdio.h>

/* ============================================================
   Utilities
   ============================================================ */




static inline double randu(void){
    return (double)rand() / (double)RAND_MAX;
}


static inline double min_image(double dx, double box){
    if(dx >  0.5*box) dx -= box;
    if(dx < -0.5*box) dx += box;
    return dx;
}


static void random_unit(double v[3]){
    const double pi = M_PI;

    // generate random numbers
    double u = randu();           // uniform random in [0,1)
    double ra2 = randu();         // another uniform random in [0,1)

    // spherical coordinates
    double theta = acos(1.0 - 2.0*u);   // theta in [0, pi]
    double phi   = 2.0 * pi * ra2;      // phi in [0, 2*pi]

    // convert to Cartesian coordinates
    v[0] = sin(theta) * cos(phi);
    v[1] = sin(theta) * sin(phi);
    v[2] = cos(theta);
}



static void rotate_point(double p[3], const double a[3], double ang){
    double c = cos(ang), s = sin(ang);
    double dot = a[0]*p[0] + a[1]*p[1] + a[2]*p[2];
    double cross[3] = {
        a[1]*p[2] - a[2]*p[1],
        a[2]*p[0] - a[0]*p[2],
        a[0]*p[1] - a[1]*p[0]
    };
    for(int i=0;i<3;i++)
        p[i] = p[i]*c + cross[i]*s + a[i]*dot*(1.0-c);
}

/* ============================================================
   Grid
   ============================================================ */

typedef struct{
    int count, max_count;
    int *idx;
} Cell;

static inline int cell_hash(int x,int y,int z,int nx,int ny,int nz){
    return x + nx*(y + ny*z);
}

/* ============================================================
   Main routine
   ============================================================ */

void relax_spherical_particles(
    double *coords, int *mol_id, int *part_type, int N,
    int n_mol, int *mol_start, int *mol_count, int *mol_movable,
    double *sigma, int n_ptype, double box[3], double cell_size,
    int iter_max, double step_trans, double step_rot,
    int max_particles_per_cell, int grid_shifting_rate
){
    static int seeded = 0;
    if(!seeded){ srand(time(NULL)); seeded = 1; }

    /* ---------------- movable molecules ---------------- */
    int *movelist = malloc(n_mol*sizeof(int));
    int mcount = 0;
    for(int m=0;m<n_mol;m++)
        if(mol_movable[m]) movelist[mcount++] = m;
    if(mcount == 0){ free(movelist); return; }

    /* ---------------- overlap flags ---------------- */
    int *mol_overlap = calloc(n_mol,sizeof(int));
    for(int m=0;m<n_mol;m++) mol_overlap[m] = 1;

    /* ---------------- grid ---------------- */
    int nx = (int)(box[0]/cell_size);
    int ny = (int)(box[1]/cell_size);
    int nz = (int)(box[2]/cell_size);
    int nc = nx*ny*nz;

  
    Cell *grid = malloc(nc * sizeof(Cell));
    for(int i = 0; i < nc; i++){
        grid[i].count = 0;
        grid[i].max_count = max_particles_per_cell;
        grid[i].idx = malloc(max_particles_per_cell * sizeof(int));
        for(int j = 0; j < max_particles_per_cell; j++)
            grid[i].idx[j] = -1;  // mark all slots as "empty"
    }



    int *p_cell = malloc(N*sizeof(int));

    int off[27][3];
    int c = 0;
    for(int dx=-1;dx<=1;dx++)
        for(int dy=-1;dy<=1;dy++)
            for(int dz=-1;dz<=1;dz++){
                off[c][0]=dx; off[c][1]=dy; off[c][2]=dz; c++;
            }

    /* ---------------- buffers ---------------- */
    int max_cnt = 0;
    for(int m=0;m<n_mol;m++)
        if(mol_count[m] > max_cnt) max_cnt = mol_count[m];

    double *bak = malloc(3*max_cnt*sizeof(double));
    int    *new_cell = malloc(max_cnt*sizeof(int));

    /* ---------------- initial cell build ---------------- */
    for(int i=0;i<N;i++){
        int ix = (int)(coords[3*i]   / cell_size);
        int iy = (int)(coords[3*i+1] / cell_size);
        int iz = (int)(coords[3*i+2] / cell_size);
        int h = cell_hash(ix,iy,iz,nx,ny,nz);
        p_cell[i] = h;
        if(grid[h].count < grid[h].max_count)
            grid[h].idx[grid[h].count++] = i;
    }

    /* ---------------- temporaries (hoisted) ---------------- */
    double com[3], disp[3], axis[3], v[3], cs[3], co[3], cso[3], dx, dy, dz, cut2;
    int i,j,k,o,nov,accept;
    int any = 1;

    /* =======================================================
       Iterations
       ======================================================= */
    int bar_width = 50;
    
    int trials = 2*mcount;
    //printf("...number of trial moves per iteration %d, \n", trials);
    double sx=0,sy=0,sz=0;
    
    int ix, iy, iz, mol, start, cnt, cmid;
    
    double *sigma2 = malloc(n_ptype*n_ptype*sizeof(double));
    for(int a=0;a<n_ptype;a++)
        for(int b=0;b<n_ptype;b++)
            sigma2[a*n_ptype+b] = sigma[a*n_ptype+b] * sigma[a*n_ptype+b];

    
    printf ("\n\nOverlap removal on progress ...\n\n");
    
    for(int iter=0; iter<iter_max; iter++){

        
        for(int t=0;t<trials;t++){

            mol   = movelist[rand()%mcount];
            start = mol_start[mol];
            cnt   = mol_count[mol];
            
            
            //printf("molecular id is given as %d \n\n", mol);
            
         

            /* ---- COM ---- */
            
            /* ---- unwrap molecule ---- */
            
            memcpy(bak, &coords[3*start], 3*cnt*sizeof(double));
            
            double ref[3];
            ref[0] = coords[3*start];
            ref[1] = coords[3*start + 1];
            ref[2] = coords[3*start + 2];

            for(k = 1; k < cnt; k++){
                i = start + k;
                coords[3*i]     = ref[0] + min_image(coords[3*i]     - ref[0], box[0]);
                coords[3*i + 1] = ref[1] + min_image(coords[3*i + 1] - ref[1], box[1]);
                coords[3*i + 2] = ref[2] + min_image(coords[3*i + 2] - ref[2], box[2]);
            }

           
            com[0]=com[1]=com[2]=0.0;
            for(k=0;k<cnt;k++){
                i = start+k;
                com[0]+=coords[3*i];
                com[1]+=coords[3*i+1];
                com[2]+=coords[3*i+2];
            }
            com[0]/=cnt; com[1]/=cnt; com[2]/=cnt;

            

            /* ---- move ---- */
            if(drand48()<0.5){
                random_unit(disp);
                disp[0]*=step_trans; disp[1]*=step_trans; disp[2]*=step_trans;
                for(k=0;k<cnt;k++){
                    i=start+k;
                    coords[3*i]+=disp[0];
                    coords[3*i+1]+=disp[1];
                    coords[3*i+2]+=disp[2];
                }
            }else{
                random_unit(axis);
                double ang=(2.0*randu()-1.0)*step_rot;
                for(k=0;k<cnt;k++){
                    i=start+k;
                    double p[3]={coords[3*i]-com[0],coords[3*i+1]-com[1],coords[3*i+2]-com[2]};
                    rotate_point(p,axis,ang);
                    coords[3*i]=com[0]+p[0];
                    coords[3*i+1]=com[1]+p[1];
                    coords[3*i+2]=com[2]+p[2];
                }
            }

            cs[0]=cs[1]=cs[2]=0.0;
            cso[0]=cso[1]=cso[2]=0.0;
            co[0]=co[1]=co[2]=0.0;
            nov = 0;

            /* ---- overlap check ---- */
            for(k=0;k<cnt;k++){
                i=start+k;
                cmid  = mol_id[i];
                

                if(coords[3*i]<0.0) coords[3*i]+=box[0];
                else if (coords[3*i]>=box[0]) coords[3*i]-=box[0];
                if(coords[3*i+1]<0.0) coords[3*i+1]+=box[1];
                else if(coords[3*i+1]>=box[1]) coords[3*i+1]-=box[1];
                if(coords[3*i+2]<0.0) coords[3*i+2]+=box[2];
                else if(coords[3*i+2]>=box[2]) coords[3*i+2]-=box[2];

                 ix=(int)(coords[3*i]/cell_size);
                 iy=(int)(coords[3*i+1]/cell_size);
                 iz=(int)(coords[3*i+2]/cell_size);
                 
                if(ix >= nx) ix = nx - 1;
                if(iy >= ny) iy = ny - 1;
                if(iz >= nz) iz = nz - 1;

                

                new_cell[k]=cell_hash(ix,iy,iz,nx,ny,nz);

                for(o=0;o<27;o++){
                    int cx=ix+off[o][0], cy=iy+off[o][1], cz=iz+off[o][2];
                    
                    if(cx<0){cx=nx-1;sx=-1;} else if(cx>=nx){cx=0;sx=1;} else sx = 0;
                    if(cy<0){cy=ny-1;sy=-1;} else if(cy>=ny){cy=0;sy=1;} else sy = 0;
                    if(cz<0){cz=nz-1;sz=-1;} else if(cz>=nz){cz=0;sz=1;} else sz = 0; 

                    Cell *cell=&grid[cell_hash(cx,cy,cz,nx,ny,nz)];
                    for(j=0;j<cell->count;j++){
                        int p=cell->idx[j];
                        if(mol_id[p]==cmid) continue;

                        double dx=coords[3*i]-(coords[3*p]+sx*box[0]);
                        double dy=coords[3*i+1]-(coords[3*p+1]+sy*box[1]);
                        double dz=coords[3*i+2]-(coords[3*p+2]+sz*box[2]);
                        double cut2=sigma2[part_type[i]*n_ptype+part_type[p]];
                        
                        //printf ("%d,   %lf   %lf   %lf  %d \n\n", p, dx, dy, dz, mol_id[p]);
                        
                        if(dx*dx+dy*dy+dz*dz<cut2){
                            cs[0]+=coords[3*i]; cs[1]+=coords[3*i+1]; cs[2]+=coords[3*i+2];
                            cso[0]+=bak[3*k];  cso[1]+=bak[3*k+1];  cso[2]+=bak[3*k+2];
                            co[0]+=coords[3*p]+sx*box[0];
                            co[1]+=coords[3*p+1]+sy*box[1];
                            co[2]+=coords[3*p+2]+sz*box[2];
                            nov++;
                        }
                    }
                }
            }

            /* ---- accept / reject ---- */
            accept = 1;
            if(nov>0){
                if(!mol_overlap[mol]) accept=0;
                else{
                    cs[0]/=nov; cs[1]/=nov; cs[2]/=nov;
                    cso[0]/=nov; cso[1]/=nov; cso[2]/=nov;
                    co[0]/=nov; co[1]/=nov; co[2]/=nov;
                    disp[0]=cs[0]-cso[0];
                    disp[1]=cs[1]-cso[1];
                    disp[2]=cs[2]-cso[2];
                    v[0]=cs[0] - co[0];
                    v[1]=cs[1] - co[1];
                    v[2]=cs[2] - co[2];
                    if(disp[0]*v[0]+disp[1]*v[1]+disp[2]*v[2]<0.0)
                        accept=0;
                }
            }
           

            if(!accept){
                memcpy(&coords[3*start], bak, 3*cnt*sizeof(double));
            }else{
                mol_overlap[mol]=(nov>0);
                for(k = 0; k < cnt; k++){
                    i = start + k;
                    int old = p_cell[i], nw = new_cell[k];
                    if(old != nw){
                        // Remove i from old cell by shifting
                        int *idx = grid[old].idx;
                        int n = grid[old].count;
                        for(j = 0; j < n; j++){
                            if(idx[j] == i){
                                // shift all subsequent entries down
                                for(int tric = j; tric < n - 1; tric++)
                                    idx[tric] = idx[tric + 1];
                                idx[n - 1] = -1; // optional: mark last slot as empty
                                grid[old].count--;
                                break;
                            }
                        }

                        // Add i to new cell
                        if(grid[nw].count < grid[nw].max_count){
                            grid[nw].idx[grid[nw].count++] = i;
                        } else {
                            // handle overflow if needed
                            printf("Warning: cell overflow\n");
                        }

                        // Update particle's current cell
                        p_cell[i] = nw;
                    }
                }

            }
        }


        if(iter % grid_shifting_rate == 0){   // every 10 iterations
            int dx = (rand() % 3) - 1;
            int dy = (rand() % 3) - 1;
            int dz = (rand() % 3) - 1;

            if(dx || dy || dz){
                double shift[3] = {
                    dx * 0.5 * box[0],
                    dy * 0.5 * box[1],
                    dz * 0.5 * box[2]
                };

                for(i = 0; i < N; i++){
                    coords[3*i]     += shift[0];
                    coords[3*i + 1] += shift[1];
                    coords[3*i + 2] += shift[2];

                    // wrap back
                    if(coords[3*i] < 0) coords[3*i] += box[0];
                    if(coords[3*i] >= box[0]) coords[3*i] -= box[0];
                    if(coords[3*i+1] < 0) coords[3*i+1] += box[1];
                    if(coords[3*i+1] >= box[1]) coords[3*i+1] -= box[1];
                    if(coords[3*i+2] < 0) coords[3*i+2] += box[2];
                    if(coords[3*i+2] >= box[2]) coords[3*i+2] -= box[2];
                }

                /* rebuild grid safely */
                for(i = 0; i < nc; i++)
                    grid[i].count = 0;

                for(i = 0; i < N; i++){
                    int ix = (int)(coords[3*i]     / cell_size);
                    int iy = (int)(coords[3*i + 1] / cell_size);
                    int iz = (int)(coords[3*i + 2] / cell_size);
                    int h = cell_hash(ix, iy, iz, nx, ny, nz);
                    p_cell[i] = h;
                    grid[h].idx[grid[h].count++] = i;
                }
            }
        }

        
        
        




        if(iter % 10 == 0) { // update every trial step
        int overlap_count = 0;
        for(int m=0; m<n_mol; m++)
            if(mol_overlap[m]) overlap_count++;

        double progress = 1.0 - ((double)overlap_count / n_mol); // fraction of molecules resolved
        int pos = (int)(bar_width * progress);

        printf("\rOverlap progress: [");
        for(int i=0; i<bar_width; i++){
            if(i < pos) printf("=");
            else if(i == pos) printf(">");
            else printf(" ");
        }
        printf("] %3.0f%% Overlap-free: %d/%d", progress*100.0, n_mol - overlap_count, n_mol);
        fflush(stdout);
        }


        any=0;
        for(int m=0;m<n_mol;m++) if(mol_overlap[m]){ any=1; break; }
        if(!any) break;
    }

    if(any) printf("\nConvergence failed increase the iteration and wait for longer period... may the force be with you !!!!!\n");
    else 
      printf("\n Overlap removed and System converged thanks for the patient!!!!! \n\n");

    free(new_cell);
    free(bak);
    for(i=0;i<nc;i++) free(grid[i].idx);
    free(grid);
    free(p_cell);
    free(movelist);
    free(mol_overlap);
    free(sigma2);

}

