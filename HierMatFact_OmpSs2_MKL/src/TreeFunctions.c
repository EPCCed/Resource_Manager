#include "TreeFunctions.h"

// Compute, recursively, the height of each node of the tree, where
// * chldtab[i] includes the head of the list of the children of node i-th
// * brthtab includes the rest of the lists which are headed in chldtab
// * hgthtab[i] includes the heigth of the node i-th
void ComputeHeigthNodes (int *chldtab, int *brthtab, int *hgthtab, int root)
{
        // Definition of the local variables
        int val, node;
        // Parameters validation
        if ((chldtab == NULL) || (brthtab == NULL) || (hgthtab == NULL) || (root < 0))
        {
                printf ("Incorrect parameters in ComputeHeigthNodes (%d)\n", root);
                exit (-1);
        }
        // val --> the height of the children is equal to 1 plus the height of the root.
        // node --> begins with the first element of the list of the children.
        val = hgthtab[root] + 1; node = chldtab[root];
        while (node != -1)
        {
                hgthtab[node] = val;  // Fix the height of node
                ComputeHeigthNodes (chldtab, brthtab, hgthtab, node);
                node = brthtab[node]; // Move to the next element of the list
        }
}

// From the treetab, this routine computes other vectors related with the tree
// * chldtab[i] includes the head of the list in which the children of node i-th are
// * brthtab includes the rest of the lists which are headed in chldtab
// * nmchtab[i] includes the size of these lists
// * hgthtab[i] includes the height of each node
// The routine returns the number of nodes of the elimination tree
int ComputeEliminationTreeVectors (int *treetab, int *chldtab, int *rchldtab, int *nmchtab, int *brthtab, int *rbrthtab, int *hgthtab, int size)
{
        // Definition of the local vectors and variables
        int i, node = 0;
        // Parameters validation
        if ((treetab == NULL) || (chldtab == NULL) || (rchldtab == NULL) || (nmchtab == NULL) || (brthtab == NULL) || (rbrthtab == NULL) || (hgthtab == NULL) || (size < 1))
        {
                printf ("Incorrect parameters in ComputeEliminationTreeVectors (%d)\n", size);
                exit (-1);
        }
        // Initialize the vectors chldtab and brthtab
        for (i=0; i < size; i++)
        {
                chldtab[i] = -1; rchldtab[i] = -1; brthtab[i] = -1; rbrthtab[i] = -1; nmchtab[i] = 0;
        }
        // Visit of the nodes of the vector treetab, until the root is reached
        i=0;
        while ((i < size) && (node >= 0))
        {
                node = treetab[i];
                if (node >= 0)
                {  // the node is not the root
                        if (chldtab[node] < 0)
                        { // Create a new list of children related to node
                                chldtab[node] = i; nmchtab[node] = 1;
                                rchldtab[node] = i;
                        } else
                        { // Insert i as the header of the list related to node
                                brthtab[i] = chldtab[node];
                                rbrthtab[chldtab[node]] = i;
                                chldtab[node] = i; nmchtab[node]++;
                        }
                }
                i++;
        }
        // Compute the vector hgthtab
        hgthtab[i-1] = 1;
        ComputeHeigthNodes (chldtab, brthtab, hgthtab, i-1);
        // Return the number of nodes of the tree
        return i;
}


