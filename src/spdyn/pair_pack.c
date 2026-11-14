#include <stdint.h>
#include <string.h>

void genesis_pair_copy(const int8_t *input, int8_t *output,
                       int rows, int columns, int input_ld, int output_ld) {
  for (int column = 0; column < columns; ++column) {
    const int8_t *source = input + column * input_ld;
    int8_t *destination = output + column * output_ld;
    int row = 0;
    for (; row + 8 <= rows; row += 8) {
      uint64_t word;
      memcpy(&word, source + row, sizeof word);
      memcpy(destination + row, &word, sizeof word);
    }
    for (; row < rows; ++row) destination[row] = source[row];
  }
}

void genesis_pair_pack(const int8_t *input, uint8_t *output, int count) {
  for (int k = 0; k < count; ++k) {
    const int8_t *p = input + 8 * k;
    output[k] = (uint8_t)((p[0] == 1) | ((p[1] == 1) << 1) |
        ((p[2] == 1) << 2) | ((p[3] == 1) << 3) |
        ((p[4] == 1) << 4) | ((p[5] == 1) << 5) |
        ((p[6] == 1) << 6) | ((p[7] == 1) << 7));
  }
}
