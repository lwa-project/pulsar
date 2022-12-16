#include <iostream>
#include <cstring>

#include "drx.hpp"
#include "lwa.hpp"

DRXBuffer::DRXBuffer(std::string filename): _sample_rate(0), _timetag_skip(0), _last_timetag(0) {
  _filename = filename;
  _fh.open(filename, std::ios::in|std::ios::binary);
  if( !_fh.good() ) {
    throw(std::runtime_error("Failed to open file"));
  }
  
  // Find valid data
  DRXFrame frame;
  _fh.read(reinterpret_cast<char*>(&frame), sizeof(frame));
  while( _fh.good() && (   (frame.header.sync_word != 0x5CDEC0DE) \
                        || (frame.header.decimation == 0) ) ) {
    _fh.seekg(1-sizeof(frame), std::ios_base::cur);
    _fh.read(reinterpret_cast<char*>(&frame), sizeof(frame));
  }
  if( !_fh.good() ) {
    throw(std::runtime_error("Failed to find valid data"));
  }
  _start = timetag_to_lwatime(__bswap_64(frame.payload.timetag));
  _fh.seekg(-sizeof(frame), std::ios_base::cur);
  
  // Determine the sample rate and timetag skip
  _sample_rate = LWA_FS / __bswap_16(frame.header.decimation);
  _timetag_skip = 4096 * __bswap_16(frame.header.decimation);
  
  // Determine the how many tuning/polarization pairs are in the data
  long marker = _fh.tellg();
  for(int i=0; i<128; i++) {
    _fh.read(reinterpret_cast<char*>(&frame), sizeof(frame));
    if( !_fh.good() ) {
      break;
    }
    _frame_ids.insert(frame.header.frame_count_word);
  }
  _fh.seekg(marker, std::ios::beg);
  
  // Determine how many frames are in the file
  _fh.seekg (0, std::ios::end);
  _nframes = _fh.tellg() / sizeof(frame);
  _fh.seekg(marker, std::ios::beg);
}

double DRXBuffer::offset(double step) {
  LWATime t, t_start, t_offset;
  
  DRXFrame frame;
  _fh.read(reinterpret_cast<char*>(&frame), sizeof(frame));
  _fh.seekg(-sizeof(frame), std::ios_base::cur);
  
  t_start = t = timetag_to_lwatime(__bswap_64(frame.payload.timetag));
  t_offset = lwatime_offset(t_start, step);
  
  int nread = 0;
  while( _fh.good() && (lwatime_diff(t_offset, t) > 0.0) ) {
    _fh.read(reinterpret_cast<char*>(&frame), sizeof(frame));
    t = timetag_to_lwatime(__bswap_64(frame.payload.timetag));
    nread++;
  }
  _fh.seekg(-sizeof(frame), std::ios_base::cur);
  nread--;
  
  _start = t;
  _nframes = _nframes - nread;
  this->reset();
  
  return lwatime_diff(t, t_start);
}

std::list<DRXFrame> DRXBuffer::get() {
  DRXFrame frame;
  uint64_t timetag;
  while(   (_buffer.size() < 20) \
        && (_fh.read(reinterpret_cast<char*>(&frame), sizeof(frame)).good()) ) {
      // Validate the sync word
      if( frame.header.sync_word != 0x5CDEC0DE ) {
        throw SyncError("Invalid sync word");
      }
      
      timetag = __bswap_64(frame.payload.timetag);
      
      // See if we have already dumped this timetag
      if( timetag < std::begin(_buffer)->first ) {
        continue;
      }
      
      // Add it to the buffer
      // If we have a new timetag, create an entry for it
      if( _buffer.count(timetag) == 0 ) {
        _buffer.emplace(timetag, std::forward<std::map<uint32_t, DRXFrame> >({}));
      }
      
      // Save the frame
      _buffer[timetag][frame.header.frame_count_word] = frame;
  }
  
  std::list<DRXFrame> output;
  // Loop over ordered ID numbers
  if( !_buffer.empty() ) {
    _buff_it = std::begin(_buffer);
    
    // But first, check for any small gaps in the buffer by looking at the time
    // tag coming out vs. what we previously had
    if( _last_timetag > 0 ) {
      double missing = (_buff_it->first - _last_timetag - _timetag_skip) / _timetag_skip;
      
      // If it looks like we are missing something, fill the gap... up to a point
      if( missing > 0 ) {
        if( ((int) missing == missing) && (missing < 50) ) {
          DRXFrame dummy_frame;
          ::memcpy(&dummy_frame, &std::begin(_buff_it->second)->second, sizeof(dummy_frame));
          ::memset(&(dummy_frame.payload.bytes), 0, 4096);
          uint64_t dummy_timetag = _buff_it->first;
          
          for(int j=1; j<missing; j++) {
            _buffer.emplace(dummy_timetag + j*_timetag_skip, std::forward<std::map<uint32_t, DRXFrame> >({}));
            
            for(const uint32_t& frame_id: _frame_ids) {
              dummy_frame.header.frame_count_word = frame_id;
              dummy_frame.payload.timetag = __bswap_64(dummy_timetag + j*_timetag_skip);
              
              _buffer[dummy_timetag + j*_timetag_skip][frame_id] = dummy_frame;
            }
          }
          
          _buff_it = std::begin(_buffer);
        } else {
          throw std::runtime_error("Invalid timetag skip encountered");
        }
      }
    }
    
    _last_timetag = _buff_it->first;
    for(const uint32_t& frame_id: _frame_ids) {
      auto _set_it = _buff_it->second.find(frame_id);
      
      if( _set_it != std::end(_buff_it->second) ) {
        // Great, we have the frame
        frame = _set_it->second;
      } else {
        // Boo, we need to create a fake frame
        ::memcpy(&frame, &std::begin(_buff_it->second)->second, sizeof(frame));
        frame.header.frame_count_word = frame_id;
        ::memset(&(frame.payload.bytes), 0, 4096);
      }
      output.push_back(frame);
    }
    
    if( output.size() > 0 ) {
      _buffer.erase(_buff_it->first);
    }
  }
  
  return output;  
}

DRXReader::DRXReader(std::string filename): DRXBuffer(filename) {
  for(uint16_t i=0; i<256; i++) {
    for(uint16_t j=0; j<2; j++) {
      int16_t t = (i >> 4*(1-j)) & 15;
      _lut[i][j] = t;
      _lut[i][j] -= ((t&8)<<1);
    }
  }
}

int8_t* DRXReader::read(double t_read, LWATime *tStart, uint64_t *samples) {
  uint32_t nframes = round((t_read * this->sample_rate()) / 4096);
  nframes = std::max(nframes, 1u);
  
  int8_t *buffer = (int8_t*) calloc(sizeof(int8_t)*2, 4*nframes*4096);
  
  uint32_t j;
  for(j=0; j<nframes; j++) {
    std::list<DRXFrame> frames = this->get();
    if( frames.size() == 0 ) {
      break;
    }
    
    if( j == 0 ) {
      auto first_frame = std::begin(frames);
      *tStart = timetag_to_lwatime(__bswap_64(first_frame->payload.timetag) \
                                  - __bswap_16(first_frame->header.time_offset));
    }
    
    for(const DRXFrame& frame: frames) {
      uint8_t tune = DRX_GET_TUNE(frame);
      uint8_t pol = DRX_GET_POLN(frame);
      
      uint64_t i = 2*(tune - 1) + pol;
      const int8_t *fp;
      for(uint32_t k=0; k<4096; k++) {
        fp = _lut[ frame.payload.bytes[k] ];
        *(buffer + i*nframes*4096*2 + j*4096*2 + k*2 + 0) = fp[0];
        *(buffer + i*nframes*4096*2 + j*4096*2 + k*2 + 1) = fp[1];
      }
    }
  }
  
  if( (j < nframes) && (j > 0) ) {
    int8_t *new_buffer = (int8_t*) malloc(sizeof(int8_t)*2 * 4*j*4096);
    for(uint64_t i=0; i<4; i++) {
      ::memcpy(new_buffer + i*j*4096*2,
               buffer + i*nframes*4096*2, 
               sizeof(int8_t)*2 * j*4096);
    }
    ::free(buffer);
    buffer = new_buffer;
  }
  *samples = j*4096;
  
  return buffer;
}
