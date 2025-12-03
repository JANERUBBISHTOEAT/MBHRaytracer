#ifndef IMG_DATA_H
#define IMG_DATA_H
#include "color.h"
#include <cstdio>
#include <vector>

/**
 * A class which manages an in-memory loaded image
 */
class img_data {
    /** The kind of image loaded */
    enum kind {
        JPEG,
        PPM,
        INVALID,
    };

  public:
    /** The image width */
    unsigned get_width() const { return width; }
    /** The image height */
    unsigned get_height() const { return height; }

  private:
    std::vector<uint8_t> raw_data;
    unsigned width;
    unsigned height;
    kind m_kind;
    /** Used during PPM parsing */
    optional<unsigned> parse_digit(FILE *img_f);
    /** Converts an x and y coordinate into an offset into the image data buffer
     */
    unsigned to_ppm_ofs(unsigned x, unsigned y) const {
        return (y * get_width() + x) * 3;
    }
    void parse_ppm(const char *filename);
    void parse_jpeg(const char *filename);

  public:
    /** A helper function which reads the pixel value at a particular x and y
     * coordinate */
    color read_pixel(unsigned x, unsigned y) const;

    /**
     * @param filename The filename to load (can either be JPEG or PPM)
     */
    img_data(const char *filename);
};
#endif
