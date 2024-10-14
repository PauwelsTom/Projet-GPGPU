#include "filter_impl.h"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <stdlib.h>
#include <string.h>
#include <thread>
#include "logo.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <time.h>

struct rgb {
    uint8_t r, g, b;
};

struct lab {
    double l, a, b;
};

extern "C" {

    ////////////////////////////// first
    void rgb2xyz(uint8_t r, uint8_t g, uint8_t b,
                    double& x, double& y, double& z)
    {
        double _r = r / 255.0;
        double _g = g / 255.0;
        double _b = b / 255.0;

        _r = _r > 0.04045 ? pow((_r + 0.055) / 1.055, 2.4) : _r / 12.92;
        _g = _g > 0.04045 ? pow((_g + 0.055) / 1.055, 2.4) : _g / 12.92;
        _b = _b > 0.04045 ? pow((_b + 0.055) / 1.055, 2.4) : _b / 12.92;

        _r *= 100.0;
        _g *= 100.0;
        _b *= 100.0;

        x = _r * 0.4124 + _g * 0.3576 + _b * 0.1805;
        y = _r * 0.2126 + _g * 0.7152 + _b * 0.0722;
        z = _r * 0.0193 + _g * 0.1192 + _b * 0.9505;
    }

    void xyz2lab(double x, double y, double z,
                    double& l, double& a, double& b)
    {
        double _x = x / 95.047;
        double _y = y / 100.000;
        double _z = z / 108.883;

        _x = _x > 0.008856 ? pow(_x, 1.0/3.0) : (7.787 * _x) + (16.0 / 116.0);
        _y = _y > 0.008856 ? pow(_y, 1.0/3.0) : (7.787 * _y) + (16.0 / 116.0);
        _z = _z > 0.008856 ? pow(_z, 1.0/3.0) : (7.787 * _z) + (16.0 / 116.0);

        l = (116.0 * _y) - 16.0;
        a = 500.0 * (_x - _y);
        b = 200.0 * (_y - _z);
    }

    double* rgb2lab(uint8_t* buffer, int width, int height, int stride)
    {
        double* lab = (double*)malloc((stride / sizeof(rgb)) * height * sizeof(struct lab));

        for (int u = 0; u < height; ++u)
        {
            rgb* lineptr = (rgb*) (buffer + u * stride);
            struct lab* lablineptr = (struct lab*) (lab + u * (stride / sizeof(rgb)));
            for (int v = 0; v < width; ++v)
            {
                double x, y, z, l, a, b;
                rgb2xyz(lineptr[v].r, lineptr[v].g, lineptr[v].b, x, y, z);
                xyz2lab(x, y, z, l, a, b);
                lab[u * (stride / sizeof(rgb)) * 3 + v + 0] = l;
                lab[u * (stride / sizeof(rgb)) * 3 + v + 1] = a;
                lab[u * (stride / sizeof(rgb)) * 3 + v + 2] = b;
            }
        }
        return lab;
    }

    double deltaE_cie76(struct lab lab1, struct lab lab2)
    {
        double dl = lab1.l - lab2.l;
        double da = lab1.a - lab2.a;
        double db = lab1.b - lab2.b;
        return sqrt(dl * dl + da * da + db * db);
    }

    uint8_t* compute_deltaE(double* lab1, double* lab2,
                                        int width, int height, int stride)
    {
        stride = (stride / sizeof(rgb));
        double* tmpdeltaE = (double*)malloc(stride * height * sizeof(double));

        double minDeltaE = DBL_MAX;
        double maxDeltaE = 0;

        for (int y = 0; y < height; ++y)
        {
            struct lab* lab1lineptr = (struct lab*) (lab1 + y * stride);
            struct lab* lab2lineptr = (struct lab*) (lab2 + y * stride);
            for (int x = 0; x < width; ++x)
            {
                struct lab l1;
                l1.l = lab1[y * stride * 3 + x + 0];
                l1.a = lab1[y * stride * 3 + x + 1];
                l1.b = lab1[y * stride * 3 + x + 2];
                struct lab l2;
                l2.l = lab2[y * stride * 3 + x + 0];
                l2.a = lab2[y * stride * 3 + x + 1];
                l2.b = lab2[y * stride * 3 + x + 2];

                tmpdeltaE[y * stride + x] = deltaE_cie76(l1,
                                                        l2);
                if (tmpdeltaE[y * stride + x] < minDeltaE)
                    minDeltaE = tmpdeltaE[y * stride + x];
                if (tmpdeltaE[y * stride + x] > maxDeltaE)
                    maxDeltaE = tmpdeltaE[y * stride + x];
            }
        }

        // normalize
        uint8_t* deltaE = (uint8_t*)malloc(stride * height * sizeof(uint8_t));
        for (int y = 0; y < height; ++y)
        {
            for (int x = 0; x < width; ++x)
            {
                deltaE[y * stride + x] = (uint8_t)
                                (255.0 * (tmpdeltaE[y * stride + x] - minDeltaE) / (maxDeltaE - minDeltaE));
            }
        }
        free(tmpdeltaE);

        return deltaE;
    }

    ////////////////////////////// second
    void create_disk_element(int radius, int* diameter, int* element)
    {
        *diameter = radius * 2;
        int center = radius;
        for (int y = 0; y < *diameter; ++y) {
            for (int x = 0; x < *diameter; ++x) {
                int dx = x - center;
                int dy = y - center;
                if (sqrt(dx * dx + dy * dy) <= radius)
                {
                    element[y * (*diameter) + x] = 1;
                }
                else
                {
                    element[y * (*diameter) + x] = 0;
                }
            }
        }
    }

    uint8_t* erode(const uint8_t* input, int width, int height, int stride,
                    int* element, int diameter)
    {
        stride = (stride / sizeof(rgb));
        uint8_t* output = (uint8_t*)malloc(stride * height * sizeof(uint8_t));
        memset(output, 0, stride * height * sizeof(uint8_t));
        int radius = diameter / 2;
        for (int y = radius; y < height - radius; ++y)
        {
            for (int x = radius; x < width - radius; ++x)
            {
                uint8_t min_val = 255;
                for (int dy = -radius; dy <= radius; ++dy)
                {
                    for (int dx = -radius; dx <= radius; ++dx)
                    {
                        if (element[(dy + radius) * diameter + (dx + radius)])
                        {
                            uint8_t val = input[(y + dy) * stride + (x + dx)];
                            if (val < min_val)
                            {
                                min_val = val;
                            }
                        }
                    }
                }
                output[y * stride + x] = min_val;
            }
        }
        return output;
    }

    uint8_t* dilate(const uint8_t* input, int width, int height, int stride,
                        int* element, int diameter)
    {
        stride = (stride / sizeof(rgb));
        uint8_t* output = (uint8_t*)malloc(stride * height * sizeof(uint8_t));
        memset(output, 0, stride * height * sizeof(uint8_t));
        int radius = diameter / 2;
        
        for (int y = radius; y < height - radius; ++y)
        {
            for (int x = radius; x < width - radius; ++x)
            {
                uint8_t max_val = 0;
                for (int dy = -radius; dy <= radius; ++dy)
                {
                    for (int dx = -radius; dx <= radius; ++dx)
                    {
                        if (element[(dy + radius) * diameter + (dx + radius)])
                        {
                            uint8_t val = input[(y + dy) * stride + (x + dx)];
                            if (val > max_val)
                            {
                                max_val = val;
                            }
                        }
                    }
                }
                output[y * stride + x] = max_val;
            }
        }
        return output;
    }

    uint8_t* morphological_opening(const uint8_t* input, int width, int height, int stride,
                                         int* element, int diameter)
    {
        uint8_t* erroded = erode(input, width, height, stride, element, diameter);
        uint8_t* dilated = dilate(erroded, width, height, stride, element, diameter);
        free(erroded);
        return dilated;
    }


    ////////////////////////////// third


    void rec_propagate_high(uint8_t* input, uint8_t* res, int width, int height, int stride,
                                int x, int y, uint8_t low_threshold)
    {
        int index = y * stride + x;
        // Check 8-connected neighbors
        for (int dy = -1; dy <= 1; ++dy)
        {
            for (int dx = -1; dx <= 1; ++dx)
            {
                if (dx == 0 && dy == 0)
                    continue; // Skip the center pixel
                int nx = x + dx;
                int ny = y + dy;
                if (nx >= 0 && nx < width && ny >= 0 && ny < height)
                {
                    int neighbor_index = ny * stride + nx;
                    if (input[neighbor_index] >= low_threshold && res[neighbor_index] != 255)
                    {
                        res[neighbor_index] = 255;
                        rec_propagate_high(input, res, width, height, stride, nx, ny, low_threshold);
                    }
                }
            }
        }
    }

    uint8_t* apply_hysteresis_threshold(uint8_t* input, int width, int height, int stride,
                                            uint8_t low_threshold, uint8_t high_threshold)
    {
        stride = (stride / sizeof(rgb));
        uint8_t* res = (uint8_t*)malloc(stride * height * sizeof(uint8_t));
        // Initialize the output image
        memset(res, 0, stride * height * sizeof(uint8_t));

        // First pass: apply high threshold
        for (int y = 0; y < height; ++y)
        {
            for (int x = 0; x < width; ++x)
            {
                int index = y * stride + x;
                if (input[index] >= high_threshold)
                {
                    res[index] = 255;
                }
            }
        }

        // Second pass: apply low threshold and connect regions
        // doesn't work well because the background image is the first image:
        //      too many pixel with a small difference(bigger than low_threshold)
        for (int y = 0; y < height; ++y)
        {
            for (int x = 0; x < width; ++x)
            {
                int index = y * stride + x;
                if (res[index] == 255)
                {
                    // Propagate the high value to connected low-threshold regions
                    rec_propagate_high(input, res, width, height, stride, x, y, low_threshold);
                }
            }
        }
        return res;
    }

    ////////////////////////////// fourth

    // Function to create a new image by combining scaled and weighted images
    void combine_images(uint8_t* r3, uint8_t* output, int width, int height, int stride, double scale)
    {
        for (int y = 0; y < height; ++y)
        {
            struct rgb* lineptr = (struct rgb*) (output + y * stride);
            for (int x = 0; x < width; ++x)
            {
                double scaled_r3 = r3[y * (stride / sizeof(rgb)) + x] * 0.5;
                lineptr[x].r = (uint8_t) ((lineptr[x].r * scale) + scaled_r3);
                lineptr[x].g = (uint8_t) (lineptr[x].g * scale);  // Only the scaled 'a' value for the green channel
                lineptr[x].b = (uint8_t) (lineptr[x].b * scale);  // Only the scaled 'a' value for the blue channel
            }
        }
    }
    





    //////////////////// background
    // Function to compare two elements (for qsort)
    int compare(const void* a, const void* b)
    {
        return (*(uint8_t*)a - *(uint8_t*)b);
    }

    // Function to compute the median of an array of images
    uint8_t* median(uint8_t** input, int num_images, int width, int height, int stride) {
        uint8_t* output = (uint8_t*)malloc(stride * height * sizeof(uint8_t));

        // Temporary array to hold pixel values from all images for a specific position
        uint8_t* pixel_values_r = (uint8_t*)malloc(num_images * sizeof(uint8_t));
        uint8_t* pixel_values_g = (uint8_t*)malloc(num_images * sizeof(uint8_t));
        uint8_t* pixel_values_b = (uint8_t*)malloc(num_images * sizeof(uint8_t));

        for (int y = 0; y < height; ++y) {
            struct rgb* outputlineptr = (struct rgb*) (output + y * stride);
            for (int x = 0; x < width; ++x) {
                int index = y * stride + x;

                // Collect pixel values from all images at (x, y)
                for (int i = 0; i < num_images; ++i) {

                    struct rgb* lineptr = (struct rgb*) (input[i] + y * stride);
                    pixel_values_r[i] = lineptr[x].r;
                    pixel_values_g[i] = lineptr[x].g;
                    pixel_values_b[i] = lineptr[x].b;
                }

                // Sort the pixel values to find the median
                qsort(pixel_values_r, num_images, sizeof(uint8_t), compare);
                qsort(pixel_values_g, num_images, sizeof(uint8_t), compare);
                qsort(pixel_values_b, num_images, sizeof(uint8_t), compare);

                // Get the median value
                if (num_images % 2 == 0) {
                    outputlineptr[x].r = (pixel_values_r[num_images / 2 - 1] + pixel_values_r[num_images / 2]) / 2;
                    outputlineptr[x].g = (pixel_values_g[num_images / 2 - 1] + pixel_values_g[num_images / 2]) / 2;
                    outputlineptr[x].b = (pixel_values_b[num_images / 2 - 1] + pixel_values_b[num_images / 2]) / 2;
                } else {
                    outputlineptr[x].r = pixel_values_r[num_images / 2];
                    outputlineptr[x].g = pixel_values_g[num_images / 2];
                    outputlineptr[x].b = pixel_values_b[num_images / 2];
                }
            }
        }

        free(pixel_values_r);
        free(pixel_values_g);
        free(pixel_values_b);
        return output;
    }


    uint8_t* readPPM(const char* filename, int* width, int* height) {
        FILE* fp = fopen(filename, "rb");
        if (!fp) {
            fprintf(stderr, "Error: Unable to open file %s\n", filename);
            return NULL;
        }

        // Read magic number (first two characters)
        char magic[3];
        fgets(magic, sizeof(magic), fp);
        if (magic[0] != 'P' || magic[1] != '6') {
            fprintf(stderr, "Error: Not a valid PPM file\n");
            fclose(fp);
            return NULL;
        }

        // Read width, height, and max color value
        fscanf(fp, "%d %d\n", width, height);
        int max_val;
        fscanf(fp, "%d\n", &max_val);

        // Allocate memory for image data
        uint8_t* data = (uint8_t*)malloc((*width) * (*height) * 3);
        if (!data) {
            fprintf(stderr, "Error: Memory allocation failed\n");
            fclose(fp);
            return NULL;
        }

        // Read pixel data
        fread(data, 1, (*width) * (*height) * 3, fp);

        fclose(fp);
        return data;
    }

    uint8_t* loadppm(const char* filename, int *width, int *height)
    {
        uint8_t* imageData = readPPM(filename, width, height);
        if (!imageData) {
            return nullptr;
        }
        return imageData;
    }
    

    void filter_impl(uint8_t* buffer, int width, int height, int stride, int pixel_stride,
                        const char *bg, int opening_size, int th_low, int th_high,
                        int bg_sampling_rate, int bg_number_frame, struct timespec time)
    {        

        static double* backgroundLab = NULL;

        if (strcmp(bg, "") == 0)
        {
            // compute background image = median of few previous
            // will be final version (too slow for cpp)
            {
                static long long start_time = 0;

                static uint8_t** backgroundIms = (uint8_t**) malloc(bg_number_frame * sizeof(uint8_t*));
                static int nbbackgroundIm = 0;
                static int nbImageMet = 0;


                if (backgroundLab == NULL)
                {
                    start_time = (time.tv_nsec + time.tv_sec * 1000000000);
                }
                long long curr_time = (time.tv_nsec + time.tv_sec * 1000000000);
                long long val = (bg_sampling_rate * nbImageMet);
                val *= 1000000;
                long long diff = (curr_time - start_time);
                if (diff >= val)
                {
                    if (nbbackgroundIm == bg_number_frame)
                    {
                        free(backgroundIms[0]);
                        for (int i = 1; i < nbbackgroundIm; ++i)
                        {
                            backgroundIms[i - 1] = backgroundIms[i];
                        }

                        --nbbackgroundIm;
                    }
                    backgroundIms[nbbackgroundIm] = (uint8_t*) malloc(stride * height * sizeof(uint8_t));
                    memcpy(backgroundIms[nbbackgroundIm], buffer, (stride * height * sizeof(uint8_t)));
                    ++nbbackgroundIm;
                    ++nbImageMet;


                    // compute background
                    uint8_t* background = median(backgroundIms, nbbackgroundIm, width, height, stride);

                    if (backgroundLab != NULL)
                        free(backgroundLab);
                    backgroundLab = rgb2lab(background, width, height, stride);
                    free(background);
                }
            }

            // // set the background to the first image (work well in cpp)
            // {
            //     if (backgroundLab == NULL)
            //     {
            //         backgroundLab = rgb2lab(buffer, width, height, stride);
            //     }
            // }
        }
        else if (backgroundLab == NULL)
        {
            // uint8_t* background = loadJPG(bg, width, height);
            int w;
            int h;
            uint8_t* background = loadppm(bg, &w, &h);

            if (!background || w != width || h != height) {
                fprintf(stderr, "Error loading ppm file\n");
                return;
            }
            backgroundLab = rgb2lab(background, width, height, stride);
            free(background);
        }


        // compute difference of current im with background
        // modify the output
        {
            // first
            double* lab = rgb2lab(buffer, width, height, stride);
            uint8_t* deltaE = compute_deltaE(backgroundLab, lab, width, height, stride);
            // second
            int diameter;
            static int* element = (int*)malloc((opening_size * 2) * (opening_size * 2) * sizeof(int));
            create_disk_element(opening_size, &diameter, element);
            uint8_t* opened = morphological_opening(deltaE, width, height, stride, element, diameter);
            // third
            uint8_t* hThreshold = apply_hysteresis_threshold(opened, width, height, stride, th_low, th_high);
            // fourth: modify the output buffer
            combine_images(hThreshold, buffer, width, height, stride, 0.5);
            free(hThreshold);
            free(opened);
            free(deltaE);
            free(lab);
        }

        // You can fake a long-time process with sleep
        {
            using namespace std::chrono_literals;
            // std::this_thread::sleep_for(100ms);
        }
    }   
}
