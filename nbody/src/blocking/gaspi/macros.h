#ifndef GASPI_MACROS_H
#define GASPI_MACROS_H

#include <GASPI.h>

#include <stdlib.h>
#include <stdio.h>

#define   LOCAL_SEGMENT_ID 0
#define REMOTE1_SEGMENT_ID 1
#define REMOTE2_SEGMENT_ID 2

#define SUCCESS_OR_DIE(f...)                                              \
    {                                                                     \
        const gaspi_return_t __r = f;                                     \
        if (__r != GASPI_SUCCESS) {                                       \
            printf("Error: '%s' [%s:%i]: %i\n",#f,__FILE__,__LINE__,__r); \
            exit (EXIT_FAILURE);                                          \
        }                                                                 \
    }

#define WAIT_FOR_ENTRIES(queue, requested_entries)                   \
    {                                                                \
        gaspi_queue_id_t __queue = queue;                            \
        gaspi_number_t __requested_entries = requested_entries;      \
        gaspi_number_t __queue_size_max;                             \
        gaspi_number_t __queue_size;                                 \
                                                                     \
        SUCCESS_OR_DIE(gaspi_queue_size_max(&__queue_size_max));     \
        SUCCESS_OR_DIE(gaspi_queue_size(__queue, &__queue_size));    \
                                                                     \
        if (__queue_size + __requested_entries > __queue_size_max) { \
            SUCCESS_OR_DIE(gaspi_wait(__queue, GASPI_BLOCK));        \
        }                                                            \
    }

#define SEGMENT_ID(nbody, ptr)                   \
({                                               \
    const void *__ptr = (void *)(ptr);           \
    const nbody_t *__nbody = (nbody_t *)(nbody); \
    gaspi_segment_id_t __id;                     \
    if (__ptr == __nbody->local) {               \
        __id = LOCAL_SEGMENT_ID;                 \
    } else if (__ptr == __nbody->remote1) {      \
        __id = REMOTE1_SEGMENT_ID;               \
    } else {                                     \
        __id = REMOTE2_SEGMENT_ID;               \
    }                                            \
    __id;                                        \
})

#endif // GASPI_MACROS_H

