#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <signal.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>



/* ============================================================
   Utilities
   ============================================================ */


static volatile sig_atomic_t stop_requested = 0;
void handle_sigint(int sig){
    (void)sig;          // suppress unused warning
    stop_requested = 1;
}


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
    int count;
    int *idx;
} Cell;

static inline int cell_hash(int x,int y,int z,int nx,int ny,int nz){
    return x + nx*(y + ny*z);
}



static inline void cell_unhash(
    int h, int nx, int ny,
    int *x, int *y, int *z)
{
    *z = h / (nx * ny);
    int rem = h % (nx * ny);
    *y = rem / nx;
    *x = rem % nx;
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

    signal(SIGINT, handle_sigint);


    static int seeded = 0;
    if(!seeded){ srand(time(NULL)); seeded = 1; }

    /* ---------------- movable molecules ---------------- */
    int *movelist = malloc(n_mol*sizeof(int));
    int mcount = 0;
    for(int m=0;m<n_mol;m++)
        if(mol_movable[m]) movelist[mcount++] = m;
    if(mcount == 0){ free(movelist); return; }
    
    
    int *mol_id_c = malloc(N*sizeof(int));
    
    

    /* ---------------- overlap flags ---------------- */
    int *mol_overlap = calloc(n_mol,sizeof(int));
    for(int m=0;m<n_mol;m++) mol_overlap[m] = 1;

    /* ---------------- grid ---------------- */
    int nx = (int)floor(box[0]/cell_size);
    int ny = (int)floor(box[1]/cell_size);
    int nz = (int)floor(box[2]/cell_size);
    
    double cell_sizex =  box[0]/(float)nx;
    box[0] =  cell_sizex*(float)nx;
    
    double cell_sizey =  box[1]/(float)ny;
    box[1] =  cell_sizey*(float)ny;
    
    double cell_sizez =  box[2]/(float)nz;
    box[2] =  cell_sizez*(float)nz;
    
    
    int nc = nx*ny*nz;

  
    Cell *grid = malloc(nc * sizeof(Cell));
    for(int i = 0; i < nc; i++){
        grid[i].count = 0;
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
    int    *overlapping_mol = malloc(n_mol*sizeof(int));
    int n_over_mol ;
    
    for (int i =0; i<n_mol; i++)
    {
        overlapping_mol[i] = -1;
    }
    n_over_mol = 0;
    

    /* ---------------- initial cell build ---------------- */
    for(int i=0;i<N;i++){
        int ix = (int)floor(coords[3*i]   / cell_sizex);
        int iy = (int)floor(coords[3*i+1] / cell_sizey);
        int iz = (int)floor(coords[3*i+2] / cell_sizez);
        
        

        
        
        int h = cell_hash(ix,iy,iz,nx,ny,nz);
        p_cell[i] = h;
        if(grid[h].count < max_particles_per_cell){
            grid[h].idx[grid[h].count++] = i;
        } else {
            printf("Overflow during rebuild: cell %d\n", h);
        }

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
    
    int flag_move;
    
    
    
    
    
    /* ============================================================
   Remap part_type IDs to contiguous range [0 .. p_type-1]
   WITHOUT changing particle order
   ============================================================ */

  int *unique_types = malloc(N * sizeof(int));
  int p_type = 0;

  /* mark all as unused */
  for(int i = 0; i < N; i++)
      unique_types[i] = -1;

  for(int i = 0; i < N; i++){
      int t = part_type[i];
      int mapped = -1;

      /* check if this type already exists */
      for(int j = 0; j < p_type; j++){
          if(unique_types[j] == t){
              mapped = j;
              break;
          }
      }

      /* if new type, assign next ID */
      if(mapped == -1){
          unique_types[p_type] = t;
          part_type[i] = p_type;
          p_type++;
      } else {
          part_type[i] = mapped;
      }
  }

  /* p_type = number of unique particle types */
  /* unique_types[new_id] -> original ID */


    
    
    /// preparation part which will be executed before the simulation will start 
    /* ---- molecule COM and relative coordinates ---- */
    double (*mol_com)[3] = malloc(n_mol * sizeof(*mol_com));
    
    
    /* ============================================================
   PREPROCESSING: unwrap molecules, compute COM, store geometry
   ============================================================ */

    for(int m = 0; m < n_mol; m++){

        int start = mol_start[m];
        int cnt   = mol_count[m];

        /* ---- backup wrapped coordinates ---- */
        memcpy(bak, &coords[3 * start], 3 * cnt * sizeof(double));

        /* ---- unwrap molecule using minimum image ---- */
        double ref[3] = {
            coords[3 * start],
            coords[3 * start + 1],
            coords[3 * start + 2]
        };

        for(int k = 1; k < cnt; k++){
            int i = start + k;
            coords[3*i]     = ref[0] + min_image(coords[3*i]     - ref[0], box[0]);
            coords[3*i + 1] = ref[1] + min_image(coords[3*i + 1] - ref[1], box[1]);
            coords[3*i + 2] = ref[2] + min_image(coords[3*i + 2] - ref[2], box[2]);
        }
        
        for(int k = 0; k < cnt; k++){
            int i = start + k;
            mol_id_c[i] = m;
        }

        /* ---- compute COM in unwrapped space ---- */
        double com[3] = {0.0, 0.0, 0.0};
        for(int k = 0; k < cnt; k++){
            int i = start + k;
            com[0] += coords[3*i];
            com[1] += coords[3*i + 1];
            com[2] += coords[3*i + 2];
        }
        com[0] /= cnt;
        com[1] /= cnt;
        com[2] /= cnt;

        /* ---- wrap COM back into simulation box ---- */
        mol_com[m][0] = fmod(com[0] + box[0], box[0]);
        mol_com[m][1] = fmod(com[1] + box[1], box[1]);
        mol_com[m][2] = fmod(com[2] + box[2], box[2]);
        
        
        for(int k = 0; k < cnt; k++){
            int i = start + k;
            coords[3*i] = fmod(coords[3*i] + box[0], box[0]);
            coords[3*i+1] = fmod(coords[3*i+1] + box[1], box[1]);
            coords[3*i+2] = fmod(coords[3*i+2] + box[2], box[2]);
        }

      
        /* ---- restore original wrapped coordinates ---- */
        memcpy(&coords[3 * start], bak, 3 * cnt * sizeof(double));
        
    }

    /* ---- adaptive MC parameters ---- */
  double step_trans_max = step_trans;
  double step_rot_max   = step_rot;

  int trans_trials = 0, trans_accept = 0;
  int rot_trials   = 0, rot_accept   = 0;

  const int adapt_interval = 100;   // how often to adapt
  const double acc_low  = 0.30;
  const double acc_high = 0.50;
  const double step_min = 0.01;

    
    
    
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

            

           
            
            for (k = 0; k<n_over_mol; k++){
                overlapping_mol[k] = -1;
            }
            n_over_mol = 0;
            
            
            
            
            
            
            if(drand48()<0.5){
                /* ---------- TRANSLATION MOVE ---------- */
                flag_move = 0;
                trans_trials++;

                random_unit(disp);
                disp[0] *= step_trans;
                disp[1] *= step_trans;
                disp[2] *= step_trans;

                for(k = 0; k < cnt; k++){
                    i = start + k;
                    coords[3*i]     += disp[0];
                    coords[3*i + 1] += disp[1];
                    coords[3*i + 2] += disp[2];
                }
            }
            else{
                /* ---------- ROTATION MOVE ---------- */
                flag_move = 1;
                rot_trials++;

                random_unit(axis);
                double ang = (2.0 * randu() - 1.0) * step_rot;

                for(k = 0; k < cnt; k++){
                    i = start + k;
                    double p[3] = {
                        coords[3*i]     - com[0],
                        coords[3*i + 1] - com[1],
                        coords[3*i + 2] - com[2]
                    };
                    rotate_point(p, axis, ang);
                    coords[3*i]     = com[0] + p[0];
                    coords[3*i + 1] = com[1] + p[1];
                    coords[3*i + 2] = com[2] + p[2];
                }

                disp[0] = disp[1] = disp[2] = 0.0;
            }

           
            cs[0]=cs[1]=cs[2]=0.0;
            cso[0]=cso[1]=cso[2]=0.0;
            co[0]=co[1]=co[2]=0.0;
            nov = 0;

            /* ---- overlap check ---- */
            
            
            accept = 1;
            for(k=0;k<cnt;k++){
                i=start+k;
                cmid  = mol_id[i];
                

                coords[3*i] = fmod(coords[3*i] + box[0], box[0]);
                coords[3*i+1] = fmod(coords[3*i+1] + box[1], box[1]);
                coords[3*i+2] = fmod(coords[3*i+2] + box[2], box[2]);


              

                ix=(int)floor(coords[3*i]/cell_sizex);
                iy=(int)floor(coords[3*i+1]/cell_sizey);
                iz=(int)floor(coords[3*i+2]/cell_sizez);
                 
               
                new_cell[k]=cell_hash(ix,iy,iz,nx,ny,nz);

                for(o=0;o<27;o++){
                    int cx=ix+off[o][0], cy=iy+off[o][1], cz=iz+off[o][2];
                    
                    if(cx<0){cx=nx-1;sx=-1;} else if(cx>=nx){cx=0;sx=1;} else sx = 0;
                    if(cy<0){cy=ny-1;sy=-1;} else if(cy>=ny){cy=0;sy=1;} else sy = 0;
                    if(cz<0){cz=nz-1;sz=-1;} else if(cz>=nz){cz=0;sz=1;} else sz = 0; 

                    Cell *cell=&grid[cell_hash(cx,cy,cz,nx,ny,nz)];
                    for(j=0;j<cell->count;j++){
                        int p=cell->idx[j];
                        if(mol_id_c[p]==mol) continue;

                        double dx=coords[3*i]-(coords[3*p]+sx*box[0]);
                        double dy=coords[3*i+1]-(coords[3*p+1]+sy*box[1]);
                        double dz=coords[3*i+2]-(coords[3*p+2]+sz*box[2]);
                        double cut2=sigma2[part_type[i]*n_ptype+part_type[p]];
                        
                        //printf ("%d,   %lf   %lf   %lf  %d  %d \n\n", p, dx, dy, dz, mol_id_c[p], N);
                        
                        if(dx*dx+dy*dy+dz*dz<cut2){
                           if (mol_overlap[mol] ==0 || mol_overlap[mol_id_c[p]]==0 ){
                                  accept =0;
                                  break;
                           }
                           else {
                                  
                                  cs[0] += ref[0] + min_image(coords[3*i] - ref[0], box[0]);    
                                  cs[1] += ref[1] + min_image(coords[3*i + 1] - ref[1], box[1]); 
                                  cs[2] += ref[2] + min_image(coords[3*i + 2] - ref[2], box[2]);;
                                  
                                  cso[0]+= ref[0] + min_image(bak[3*k] - ref[0], box[0]);   
                                  cso[1]+= ref[1] + min_image(bak[3*k + 1] - ref[1], box[1]); 
                                  cso[2]+= ref[2] + min_image(bak[3*k + 2] - ref[2], box[2]);
                                  
                                  co[0] += ref[0] + min_image(coords[3*p] - ref[0], box[0]);    
                                  co[1] += ref[1] + min_image(coords[3*p + 1] - ref[1], box[1]); 
                                  co[2] += ref[2] + min_image(coords[3*p + 2] - ref[2], box[2]);;
                                  
                                  nov++;
                                  
                                  int found  = 0;
                                  for (int idx  = 0 ; idx <n_over_mol; idx++){
                                      if (overlapping_mol[idx] == mol_id_c[p]) { found = 1; break;}
                                  }
                                  if (found==0){  overlapping_mol[n_over_mol] =  mol_id_c[p]; n_over_mol++;}
                                  
                                  
                           }
                        }
                    }
                    
                    if (accept ==0) break;
                }
                if (accept ==0) break;
            }

            /* ---- accept / reject ---- */
            
            
            
            if(accept ==1 && nov>0){
            
                if (flag_move ==1){
            
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
                else {

                    /* -------------------------------------------
                       Compute COM of overlapping molecules
                       unwrapped w.r.t. moved molecule COM
                       ------------------------------------------- */
                    double ref_com[3] = {
                        mol_com[mol][0],
                        mol_com[mol][1],
                        mol_com[mol][2]
                    };

                    /* wrap reference COM just to be safe */
                    ref_com[0] = fmod(ref_com[0] + box[0], box[0]);
                    ref_com[1] = fmod(ref_com[1] + box[1], box[1]);
                    ref_com[2] = fmod(ref_com[2] + box[2], box[2]);

                    com[0] = com[1] = com[2] = 0.0;

                    for(j = 0; j < n_over_mol; j++){
                        int m2 = overlapping_mol[j];

                        /* unwrap overlapping molecule COM
                           relative to reference COM */
                        double ux = ref_com[0] +
                                    min_image(mol_com[m2][0] - ref_com[0], box[0]);
                        double uy = ref_com[1] +
                                    min_image(mol_com[m2][1] - ref_com[1], box[1]);
                        double uz = ref_com[2] +
                                    min_image(mol_com[m2][2] - ref_com[2], box[2]);

                        com[0] += ux;
                        com[1] += uy;
                        com[2] += uz;
                    }

                    com[0] /= n_over_mol;
                    com[1] /= n_over_mol;
                    com[2] /= n_over_mol;

                    /* -------------------------------------------
                       Directional acceptance test
                       ------------------------------------------- */
                    double dx = ref_com[0] - com[0];
                    double dy = ref_com[1] - com[1];
                    double dz = ref_com[2] - com[2];

                    if(disp[0]*dx + disp[1]*dy + disp[2]*dz < 0.0)
                        accept = 0;
                      
                }

                
            }
            
            
            if(!accept){
            
                memcpy(&coords[3*start], bak, 3*cnt*sizeof(double));
                 
            }
            else{
                mol_overlap[mol]=(nov>0);
            
                if(flag_move == 0) trans_accept++;
                else               rot_accept++;
                
                mol_com[mol][0] = fmod(mol_com[mol][0] + disp[0] + box[0], box[0]);
                mol_com[mol][1] = fmod(mol_com[mol][1] + disp[1] + box[1], box[1]);
                mol_com[mol][2] = fmod(mol_com[mol][2] + disp[2] + box[2], box[2]);
                
                for(k = 0; k < cnt; k++){
                    i = start + k;
                    int old = p_cell[i], nw = new_cell[k];
                    if(old != nw){
                        // Remove i from old cell by shifting
                        // int *idxx = grid[old].idx;
                        int n = grid[old].count;
                        
                        int found = 0;
                        for(j = 0; j < n; j++){
                            if(grid[old].idx[j] == i){
                                for(int tric = j; tric < n-1; tric++)
                                    grid[old].idx[tric] = grid[old].idx[tric + 1];
                                grid[old].idx[n - 1] = -1;
                                grid[old].count--;
                                found = 1;
                                break;
                            }
                        }

                        if(!found){
                            printf("ERROR: particle %d not found in old cell %d\n", i, old);
                            exit(0);
                            
                        }

                        
                        
                        
                        // Add i to new cell
                        if(grid[nw].count < max_particles_per_cell){
                            grid[nw].idx[grid[nw].count] = i;
                            grid[nw].count++;
                        } else {
                            // handle overflow if needed
                            printf("\n\nWarning: cell overflow: for max size cell ");
                            exit(0);
                        }
                        // Update particle's current cell
                        p_cell[i] = nw;
                    }
                }
            }
        }

        if(stop_requested){
            printf("\n\nCtrl+C detected — stopping relaxation cleanly.\n");
            break;
        }
     
        
        /* ---------- adaptive step control ---------- */
        if(iter % adapt_interval == 0 && (trans_trials + rot_trials) > 0){

            

           
                double acc = (double)trans_accept / trans_trials;
                if(acc > acc_high)
                    step_trans = fmin(step_trans * 1.1, step_trans_max);
                else if(acc < acc_low)
                    step_trans *= 0.9;

                /* enforce lower bound */
                step_trans = fmax(step_trans, step_min);
            

            
                acc = (double)rot_accept / rot_trials;
                if(acc > acc_high)
                    step_rot = fmin(step_rot * 1.1, step_rot_max);
                else if(acc < acc_low)
                    step_rot *= 0.9;

                /* enforce lower bound */
                step_rot = fmax(step_rot, step_min);
            

            /* reset counters */
            trans_trials = trans_accept = 0;
            rot_trials   = rot_accept   = 0;
        }




        // grid shifting in grid_shifting_rate steps     
        if(iter % grid_shifting_rate == 0){   // every grid_shifting_rate iterations
            // Choose one random direction: 0 = x, 1 = y, 2 = z
            int dir = rand() % 3;
            int dx = 0, dy = 0, dz = 0;

            if(dir == 0) dx = 1;   // shift in x-direction
            else if(dir == 1) dy = 1; // shift in y-direction
            else dz = 1;           // shift in z-direction


            if(dx || dy || dz){
            
               // printf("\nI am just checking |||||\n");
                double shift[3] = {
                    dx * 0.43 * box[0],
                    dy * 0.43 * box[1],
                    dz * 0.43 * box[2]
                };

                
                for(int i = 0; i < N; i++){
                    
                    coords[3*i]     += shift[0];
                    coords[3*i + 1] += shift[1];
                    coords[3*i + 2] += shift[2];

                    coords[3*i] = fmod(coords[3*i] + box[0], box[0]);
                    coords[3*i+1] = fmod(coords[3*i+1] + box[1], box[1]);
                    coords[3*i+2] = fmod(coords[3*i+2] + box[2], box[2]);

                    
                }

               
                for(int mol = 0; mol < n_mol; mol++){
                 
                      mol_com[mol][0] += shift[0];
                      mol_com[mol][1] += shift[1];
                      mol_com[mol][2] += shift[2];

                      mol_com[mol][0] = fmod(mol_com[mol][0] + box[0], box[0]);
                      mol_com[mol][1] = fmod(mol_com[mol][1] + box[1], box[1]);
                      mol_com[mol][2] = fmod(mol_com[mol][2] + box[2], box[2]);
                }


                // ---------- Clear all cells ----------
                
                
                for(int c = 0; c < nc; c++){
                    for(int j = 0; j < max_particles_per_cell; j++)
                        grid[c].idx[j] = -1;
                    grid[c].count = 0;
                }

               
                
                for(int it = 0; it < N; it++){
                    int ix = (int)floor(coords[3*it] / cell_sizex);
                    int iy = (int)floor(coords[3*it + 1] / cell_sizey);
                    int iz = (int)floor(coords[3*it + 2] / cell_sizez);

  

                    int h = cell_hash(ix, iy, iz, nx, ny, nz);

                   

                    p_cell[it] = h;

                    if(grid[h].count < max_particles_per_cell){
                        grid[h].idx[grid[h].count++] = it;
                    } else {
                        printf("Warning: Overflow during rebuild at cell %d\n", h);
                    }
                } 

               

                for(int m=0;m<n_mol;m++) mol_overlap[m] = 1;
            
            }
            
        }

        if(iter % 100 == 0) { // update every trial step
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
              printf("] %3.0f%% Overlap-free: %d/%d, iteration no: %d, with grid shifting interval: %d ", progress*100.0, n_mol - overlap_count, n_mol, iter, grid_shifting_rate);
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

