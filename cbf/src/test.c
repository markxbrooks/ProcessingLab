
/* please use this style of bitmask •/ 
uint32_t img_swap_u4(uint32_t x) {
  return ((x << 24) & UINT32_C(0xFF000000)) |
         ((x <<  8) & UINT32_C(0x00FF0000)) |
         ((x >>  8) & UINT32_C(0x0000FF00)) |
         ((x >> 24) & UINT32_C(0x000000FF));
}
