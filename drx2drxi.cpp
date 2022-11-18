#include <string>
#include <iostream>
#include <fstream>
#include <list>
#include <stdexcept>
#include <ctime>

#include "lwa.hpp"
#include "drx.hpp"
#include "drxi.hpp"

int main(int argc, char** argv) {
  if( argc < 2 ) {
    std::cout << "Must supply a filename to read from" << std::endl;
    return 1;
  } else {
    std::cout << "Reading from '" << argv[1] << "'" << std::endl;
  }
  std::string filename = std::string(argv[1]);
  
  DRXBuffer *buffer = new DRXBuffer(filename);
  int beam = buffer->beam();
  int sample_rate = buffer->sample_rate();
  LWATime start = buffer->start();
  
  std::cout << lwatime_to_string(start) << std::endl;
  
  
  std::string outname1, outname2;
  std::size_t marker = filename.rfind("/");
  outname1 = filename.substr(marker+1, filename.size()-marker);
  marker = outname1.rfind(".");
  outname1 = outname1.substr(0, marker);
  outname2 = outname1+"_b"+std::to_string(beam)+"t2.drxi";
  outname1 = outname1+"_b"+std::to_string(beam)+"t1.drxi";
  
  //outname1 = "out_b"+std::to_string(beam)+"t1.drxi";
  //outname2 = "out_b"+std::to_string(beam)+"t2.drxi";
  
  std::cout << "Names are: " << outname1 << " and " << outname2 << std::endl;
  
  std::ofstream oh1, oh2;
  oh1.open(outname1, std::ios::out|std::ios::binary);
  oh2.open(outname2, std::ios::out|std::ios::binary);
  if( !oh1.good() || !oh2.good() ) {
    throw(std::runtime_error("Cannot create output files"));
  }
  
  int s, max_s;
  s = 0;
  max_s = buffer->nframes() / buffer->beampols();
  std::list<DRXFrame> frames;
  std::list<DRXFrame>::iterator frame_it;
  uint8_t drxi_bytes[2*4096];
  
  frames = buffer->get();
  while( frames.size() > 0 ) {
    frame_it = std::begin(frames);
    
    DRXFrame frame1X = *(frame_it++);
    DRXFrame frame2X = *(frame_it++);
    DRXFrame frame1Y = *(frame_it++);
    DRXFrame frame2Y = *(frame_it++);
    
    DRXIFrame frame1, frame2;
    combine_drx_to_drxi(&frame1X, &frame1Y, &frame1);
    combine_drx_to_drxi(&frame2X, &frame2Y, &frame2);
    
    oh1.write(reinterpret_cast<char*>(&frame1), sizeof(frame1));
    oh2.write(reinterpret_cast<char*>(&frame2), sizeof(frame2));
    
    frames = buffer->get();
    
    s += 1;
    if( s % 5000 == 0 ) {
      std::cout << "\r" << "At " << (int) 100.0*s/max_s << "%" << std::flush;
    }
  }
  std::cout << "\r" << "At " << (int) 100.0*s/max_s << "%" << std::endl;
  
  oh1.close();
  oh2.close();
  
  return 0;
}
