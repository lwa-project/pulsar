#include "drxi.hpp"

void combine_drx_to_drxi(const DRXFrame *frameX,
                         const DRXFrame *frameY,
                         DRXIFrame *frameI) {
  uint64_t timetag = __bswap_64(frameX->payload.timetag) - __bswap_16(frameX->header.time_offset);
  
  ::memcpy(frameI, frameX, sizeof(DRXIFrame)-8192);
  
  frameI->header.frame_count_word |= (1 << 6);
  frameI->header.time_offset = 0;
  frameI->payload.timetag = __bswap_64(timetag);
  
  for(int i=0; i<4096; i++) {\
    frameI->payload.bytes[2*i+0] = frameX->payload.bytes[i];
    frameI->payload.bytes[2*i+1] = frameY->payload.bytes[i];
  }
}
