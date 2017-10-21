#include <string>

bool read_png(std::string png_filename,
              void**      pixel_data,
              size_t*     width,
              size_t*     height,
              bool        internal_format_rgba);
bool read_png(std::string png_filename,
              void**      pixel_data,
              size_t*     width,
              size_t*     height);
