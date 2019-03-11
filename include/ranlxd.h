#ifndef RANLXD_H
#define RANLXD_H

extern void ranlxd(double r[],int n);
extern void rlxd_init(int level,int seed);
extern int rlxd_size(void);
extern void rlxd_get(int state[]);
extern void rlxd_reset(int state[]);

#endif
