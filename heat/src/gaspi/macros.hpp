#ifndef MACROS_HPP
#define MACROS_HPP

#include <GASPI.h>
#include <cstdlib>
#include <cstdio>

#include "common/matrix.hpp"

// Compute the first compute (FCR) row offset
#define FCR_OFFSET(nbx, nby)                             \
    ({                                                   \
        gaspi_offset_t __offset = (nby)*sizeof(block_t); \
        __offset;                                        \
    })

// Compute the last compute row (LCR) offset
#define LCR_OFFSET(nbx, nby)                                    \
    ({                                                          \
        gaspi_offset_t __offset = ((BSX-1)*BSY)*sizeof(double); \
        __offset += ((nbx-2)*(nby))*sizeof(block_t);            \
        __offset;                                               \
    })

// Compute the upper border (UB) offset
#define UB_OFFSET(nbx, nby)                                     \
    ({                                                          \
        gaspi_offset_t __offset = ((BSX-1)*BSY)*sizeof(double); \
        __offset;                                               \
    })

// Compute the lower border (LB) offset
#define LB_OFFSET(nbx, nby)                                        \
    ({                                                             \
        gaspi_offset_t __offset = ((nbx-1)*(nby))*sizeof(block_t); \
        __offset;                                                  \
    })

// Compute the offset given a base and the current column block
#define OFFSET(base, by)                                       \
    ({                                                         \
        gaspi_offset_t __offset = (base)+(by)*sizeof(block_t); \
        __offset;                                              \
    })

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

#endif // MACROS_HPP

