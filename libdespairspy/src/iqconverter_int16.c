/*
Copyright (C) 2014, Youssef Touil <youssef@airspy.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#include "iqconverter_int16.h"

#include <stdlib.h>
#include <string.h>

#if defined(__MINGW32__) && !defined(__MINGW64_VERSION_MAJOR)
  #include <malloc.h>
  #define _aligned_malloc __mingw_aligned_malloc
  #define _aligned_free  __mingw_aligned_free
  #define _inline inline
#elif defined(__APPLE__)
  #include <malloc/malloc.h>
  #define _aligned_malloc(size, alignment) malloc(size)
  #define _aligned_free(mem) free(mem)
  #define _inline inline
#elif defined(__GNUC__) && !defined(__MINGW64_VERSION_MAJOR)
  #include <malloc.h>
  #define _aligned_malloc(size, alignment) memalign(alignment, size)
  #define _aligned_free(mem) free(mem)
  #define _inline inline
#endif

#define SIZE_FACTOR 16
#define DEFAULT_ALIGNMENT 16

#define SAMPLE_RESOLUTION 12
#define SAMPLE_ENCAPSULATION 15

#define SAMPLE_SHIFT (SAMPLE_ENCAPSULATION - SAMPLE_RESOLUTION)

int iqconverter_int16_init(iqconverter_int16_t *cnv, const int16_t *hb_kernel, int len)
{
    int ret = 0;
	int i;
	size_t buffer_size;

	cnv->len = len / 2 + 1;

	buffer_size = cnv->len * sizeof(int32_t);

	if (NULL == (cnv->fir_kernel = (int32_t *) _aligned_malloc(buffer_size, DEFAULT_ALIGNMENT))) {
        goto done;
    }

	if (NULL == (cnv->fir_queue = (int32_t *) _aligned_malloc(buffer_size * SIZE_FACTOR, DEFAULT_ALIGNMENT))) {
        goto done;
    }

	if (NULL == (cnv->delay_line = (int16_t *) _aligned_malloc(buffer_size / 4, DEFAULT_ALIGNMENT))) {
        goto done;
    }

	iqconverter_int16_reset(cnv);

	for (i = 0; i < cnv->len; i++)
	{
		cnv->fir_kernel[i] = hb_kernel[i * 2];
	}

done:
    if (0 != ret) {
        if (NULL != cnv->fir_kernel) {
            _aligned_free(cnv->fir_kernel);
        }
        if (NULL != cnv->fir_queue) {
            _aligned_free(cnv->fir_queue);
        }

        if (NULL != cnv->delay_line) {
            _aligned_free(cnv->delay_line);
        }
    }
	return ret;
}

void iqconverter_int16_free(iqconverter_int16_t *cnv)
{
	_aligned_free(cnv->fir_kernel);
	_aligned_free(cnv->fir_queue);
	_aligned_free(cnv->delay_line);
}

void iqconverter_int16_reset(iqconverter_int16_t *cnv)
{
	cnv->fir_index = 0;
	cnv->delay_index = 0;
	cnv->old_x = 0;
	cnv->old_y = 0;
	cnv->old_e = 0;
	memset(cnv->delay_line, 0, cnv->len * sizeof(int16_t) / 4);
	memset(cnv->fir_queue, 0, cnv->len * sizeof(int16_t) * SIZE_FACTOR);
}

static void fir_interleaved(iqconverter_int16_t *cnv, int16_t *samples, int len)
{
	int i;
	int j;
	int fir_index;
	int fir_len;
	int32_t *queue;
	int32_t acc;

	fir_len = cnv->len;
	fir_index = cnv->fir_index;

	for (i = 0; i < len; i += 2)
	{
		queue = cnv->fir_queue + fir_index;

		queue[0] = samples[i];

		acc = 0;

		// Auto vectorization works on VS2012, VS2013 and GCC
		for (j = 0; j < fir_len; j++)
		{
			acc += cnv->fir_kernel[j] * queue[j];
		}

		if (--fir_index < 0)
		{
			fir_index = cnv->len * (SIZE_FACTOR - 1);
			memcpy(cnv->fir_queue + fir_index + 1, cnv->fir_queue, (cnv->len - 1) * sizeof(int32_t));
		}

		samples[i] = acc >> 15;
	}

	cnv->fir_index = fir_index;
}

static void delay_interleaved(iqconverter_int16_t *cnv, int16_t *samples, int len)
{
	int i;
	int index;
	int half_len;
	int16_t res;

	half_len = cnv->len >> 1;
	index = cnv->delay_index;

	for (i = 0; i < len; i += 2)
	{
		res = cnv->delay_line[index];
		cnv->delay_line[index] = samples[i];
		samples[i] = res;

		if (++index >= half_len)
		{
			index = 0;
		}
	}

	cnv->delay_index = index;
}

/*
 * Given a sample, process an iteration of the DC removal IIR. 
 */
static inline
int16_t _remove_dc_sample(int32_t sample,
        int32_t old_e, int16_t old_x, int16_t old_y,
        int32_t *new_e, int16_t *new_x, int16_t *new_y)
{
	int32_t u;
	int16_t x, y, w, s;

    x = (sample - 2048) << SAMPLE_SHIFT;
    w = x - old_x;
    u = old_e + (int32_t) old_y * 32100;
    s = u >> 15;
    y = w + s;
    *new_e = u - (s << 15);
    *new_x = x;
    *new_y = y;
    return y;
}

/**
 * Apply the DC blocker filter to the input buffer
 */
static
void remove_dc(iqconverter_int16_t *cnv, uint16_t *samples_raw, int len)
{
	int i;
	int32_t old_e;
	int16_t old_x, old_y;
    int16_t *samples = (int16_t *)samples_raw;

	old_x = cnv->old_x;
	old_y = cnv->old_y;
	old_e = cnv->old_e;

    /* Apply the DC blocker filter, and update the output samples for I/Q calculation */
	for (i = 0; i < len; i += 4)
	{
        samples[i + 0] =
            -_remove_dc_sample(samples_raw[i + 0], old_e, old_x, old_y, &old_e, &old_x, &old_y);
        samples[i + 1] =
            (-_remove_dc_sample(samples_raw[i + 1], old_e, old_x, old_y, &old_e, &old_x, &old_y)) >> 1;
        samples[i + 2] =
            _remove_dc_sample(samples_raw[i + 2], old_e, old_x, old_y, &old_e, &old_x, &old_y);
        samples[i + 3] =
            _remove_dc_sample(samples_raw[i + 3], old_e, old_x, old_y, &old_e, &old_x, &old_y) >> 1;
	}

	cnv->old_x = old_x;
	cnv->old_y = old_y;
	cnv->old_e = old_e;
}

static void translate_fs_4(iqconverter_int16_t *cnv, int16_t *samples, int len)
{
	fir_interleaved(cnv, samples, len);
	delay_interleaved(cnv, samples + 1, len);
}

void iqconverter_int16_process(iqconverter_int16_t *cnv, uint16_t *samples, int len)
{
	remove_dc(cnv, samples, len);
	translate_fs_4(cnv, (int16_t *)samples, len);
}

