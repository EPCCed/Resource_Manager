/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 **/
#include "plasma_workspace.h"
#include "plasma_internal.h"

#include <starpu.h>

/******************************************************************************/
int plasma_workspace_create(plasma_workspace_t *workspace, size_t lworkspace,
                            plasma_enum_t dtyp)
{
    // Allocate array of pointers.
    workspace->nthread = starpu_cpu_worker_get_count();
    workspace->lworkspace = lworkspace;
    workspace->dtyp  = dtyp;
    if ((workspace->spaces = (void**)calloc(workspace->nthread,
                                            sizeof(void*))) == NULL) {
        free(workspace->spaces);
        plasma_error("malloc() failed");
        return PlasmaErrorOutOfMemory;
    }

    // Each thread allocates its workspace.
    size_t size = (size_t)lworkspace * plasma_element_size(workspace->dtyp);
    int info = PlasmaSuccess;
    //int tid = starpu_worker_get_id();
    for (int tid = 0; tid < workspace->nthread; tid++) {
       if ((workspace->spaces[tid] = (void*)malloc(size)) == NULL) {
           info = PlasmaErrorOutOfMemory;
       }
    }
    if (info != PlasmaSuccess) {
        plasma_workspace_destroy(workspace);
    }
    return info;
}

/******************************************************************************/
int plasma_workspace_destroy(plasma_workspace_t *workspace)
{
    if (workspace->spaces != NULL) {
        for (int i = 0; i < workspace->nthread; ++i) {
            free(workspace->spaces[i]);
            workspace->spaces[i] = NULL;
        }
        free(workspace->spaces);
        workspace->spaces  = NULL;
        workspace->nthread = 0;
        workspace->lworkspace   = 0;
    }
    return PlasmaSuccess;
}
