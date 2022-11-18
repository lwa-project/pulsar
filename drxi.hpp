#ifndef __INCLUDE_DRXI_HPP_
#define __INCLUDE_DRXI_HPP_

#include <cstring>
#include "drx.hpp"

typedef struct __attribute__((packed)) {
    uint32_t sync_word;
    union {
        struct {
            uint32_t frame_count:24;
            uint8_t id:8;
        };
        /* Also: 
            struct {
                uint32_t frame_count:24;
                uint8_t  beam:3;
                uint8_t  tune:3;
                uint8_t  reserved:1;
                uint8_t  pol:1;
            };
        */
        uint32_t frame_count_word;
    };
    uint32_t second_count;
    uint16_t decimation;
    uint16_t time_offset;
} DRXIHeader;


typedef struct __attribute__((packed)) {
    uint64_t timetag;
    uint32_t tuning_word;
    uint32_t flags;
    uint8_t  bytes[8192];
} DRXIPayload;


typedef struct __attribute__((packed)) {
    DRXIHeader header;
    DRXIPayload payload;
} DRXIFrame;

void combine_drx_to_drxi(const DRXFrame *frameX, const DRXFrame *frameY, DRXIFrame *frameI);

#endif // __INCLUDE_DRXI_HPP_
