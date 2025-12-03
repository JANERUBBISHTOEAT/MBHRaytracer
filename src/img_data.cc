#include "img_data.h"
#include <jpeglib.h>

color img_data::read_pixel(unsigned x, unsigned y) const {
    // TODO: Assert px is in range.
    if (raw_data.empty()) {
        assert(false);
    }

    const uint8_t *col = raw_data.data() + to_ppm_ofs(x, y);
    return to_color(col[0], col[1], col[2]);
}
optional<unsigned> img_data::parse_digit(FILE *img_f) {
    unsigned dig = 0;
    char c;
    size_t ret = -1;
    while (ret != 0) {
        ret = fread(&c, sizeof(char), 1, img_f);

        if (ret != 1)
            return nullopt;

        unsigned single_dig = c - '0';

        if (single_dig > 10)
            break;

        dig *= 10;
        dig += single_dig;
    }
    return dig;
}

void img_data::parse_ppm(const char *filename) {
    FILE *img_f = fopen(filename, "r");

    char magic[3];
    size_t amt = fread(magic, sizeof(char), 3, img_f);
    assert(amt == 3);
    assert(magic[0] == 'P' and magic[1] == '6');
    optional<unsigned> l_width = parse_digit(img_f);
    optional<unsigned> l_height = parse_digit(img_f);
    optional<unsigned> max_value = parse_digit(img_f);
    // ignored for now
    (void)max_value;

    assert(l_width and l_height);
    this->width = *l_width;
    this->height = *l_height;

    unsigned size = width * height * 3;
    raw_data.resize(size);

    size_t res = fread(raw_data.data(), sizeof(char), size, img_f);
    assert(res == size);
}

img_data::img_data(const char *filename) {
    m_kind = kind::INVALID;
    if (ends_with(filename, "jpg")) {
        m_kind = kind::JPEG;
        parse_jpeg(filename);
    } else if (ends_with(filename, "ppm")) {
        m_kind = kind::PPM;
        parse_ppm(filename);
    }
    assert(m_kind != kind::INVALID);
}

void img_data::parse_jpeg(const char *filename) {
    FILE *infile = fopen(filename, "rb");
    assert(infile != nullptr);

    jpeg_decompress_struct cinfo;
    jpeg_error_mgr jerr;
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&cinfo);
    jpeg_stdio_src(&cinfo, infile);
    jpeg_read_header(&cinfo, TRUE);
    jpeg_start_decompress(&cinfo);

    width = cinfo.output_width;
    height = cinfo.output_height;
    int channels = cinfo.output_components;
    assert(channels >= 3);

    raw_data.resize(width * height * 3);
    JSAMPARRAY buffer = (*cinfo.mem->alloc_sarray)((j_common_ptr)&cinfo,
                                                   JPOOL_IMAGE,
                                                   width * channels, 1);
    while (cinfo.output_scanline < cinfo.output_height) {
        jpeg_read_scanlines(&cinfo, buffer, 1);
        unsigned int y = cinfo.output_scanline - 1;
        for (unsigned int x = 0; x < width; x++) {
            uint8_t *dest = raw_data.data() + to_ppm_ofs(x, y);
            dest[0] = buffer[0][x * channels + 0];
            dest[1] = buffer[0][x * channels + 1];
            dest[2] = buffer[0][x * channels + 2];
        }
    }

    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    fclose(infile);
}
