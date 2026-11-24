#include "utils/commons.h"
#include "wavefront/wavefront_pcigar.h"

int main(void) {
  char pattern[42], text[43], output[48];
  for (int prefix = 0; prefix <= 32; ++prefix) {
    memset(pattern,'A',sizeof(pattern));
    memset(text,'A',sizeof(text));
    text[3+prefix] = 'C';
    wavefront_sequences_t sequences = {0};
    sequences.mode = wf_sequences_ascii;
    sequences.pattern = pattern+1;
    sequences.text = text+3;
    sequences.pattern_length = sequences.text_length = prefix+1;
    for (int offset = 0; offset < 8; ++offset) {
      memset(output,'?',sizeof(output));
      int v = 0, h = 0, length = 0;
      pcigar_unpack_linear(PCIGAR_MISMATCH,&sequences,&v,&h,output+offset,&length);
      if (v != prefix+1 || h != prefix+1 || length != prefix+1) return 1;
      for (int i = 0; i < (int)sizeof(output); ++i) {
        const char expected = (i < offset || i > offset+prefix) ? '?' :
                              i == offset+prefix ? 'X' : 'M';
        if (output[i] != expected) return 1;
      }
    }
  }
  return 0;
}
