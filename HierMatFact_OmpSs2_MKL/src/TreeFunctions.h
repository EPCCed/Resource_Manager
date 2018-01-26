#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <errno.h>
#include <math.h>
#include <sys/time.h>
#include <sys/times.h>
#include <malloc.h>

// Compute, recursively, the height of each node of the tree, where
// * chldtab[i] includes the head of the list of the children of node i-th
// * brthtab includes the rest of the lists which are headed in chldtab
// * hgthtab[i] includes the heigth of the node i-th
void ComputeHeigthNodes (int *chldtab, int *brthtab, int *hgthtab, int root);

// From the treetab, this routine computes other vectors related with the tree
// * chldtab[i] includes the head of the list in which the children of node i-th are
// * brthtab includes the rest of the lists which are headed in chldtab
// * nmchtab[i] includes the size of these lists
// * hgthtab[i] includes the height of each node
// The routine returns the number of nodes of the elimination tree
int ComputeEliminationTreeVectors (int *treetab, int *chldtab, int *rchldtab, int *nmchtab, int *brthtab, int *rbrthtab, int *hgthtab, int size);
