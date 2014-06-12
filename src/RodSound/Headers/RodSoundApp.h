//
//  RodSoundApp.h
//  Visualizer
//
//  Created by eschweickart on 3/5/14.
//
//

#ifndef __Visualizer__SimulatorApp__
#define __Visualizer__SimulatorApp__

#include <iostream>

const int SampleRate = 44100;

template <typename T>
static uint16_t inline toSample(const T val, const T max) {
  return (uint16_t) ((val / max) * INT16_MAX);
}

// Tools for writing WAV files to disk.
// Code copied from http://joshparnell.com/blog/2013/03/21/how-to-write-a-wav-file-in-c/
#include <fstream>

template <typename T>
void write(std::ofstream& stream, const T& t) {
  stream.write((const char*)&t, sizeof(T));
}

template <typename T>
void writeFormat(std::ofstream& stream) {
  write<short>(stream, 1);
}

template <>
void writeFormat<float>(std::ofstream& stream) {
  write<short>(stream, 3);
}

template <typename SampleType>
void writeWAVData(
                  char const* outFile,
                  SampleType* buf,
                  size_t bufSize,
                  int sampleRate,
                  short channels)
{
  std::ofstream stream(outFile, std::ios::binary);
  stream.write("RIFF", 4);
  write<int>(stream, 36 + bufSize);
  stream.write("WAVE", 4);
  stream.write("fmt ", 4);
  write<int>(stream, 16);
  writeFormat<SampleType>(stream);                                // Format
  write<short>(stream, channels);                                 // Channels
  write<int>(stream, sampleRate);                                 // Sample Rate
  write<int>(stream, sampleRate * channels * sizeof(SampleType)); // Byterate
  write<short>(stream, channels * sizeof(SampleType));            // Frame size
  write<short>(stream, 8 * sizeof(SampleType));                   // Bits per sample
  stream.write("data", 4);
  stream.write((const char*)&bufSize, 4);
  stream.write((const char*)buf, bufSize);
  stream.close();
}


#endif /* defined(__Visualizer__SimulatorApp__) */
