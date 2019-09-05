
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// core/film2d.cpp*
#include "film2d.h"
#include "paramset.h"
#include "imageio.h"
#include "stats.h"

namespace pbrt {

STAT_MEMORY_COUNTER("Memory/Film2d pixels", filmPixelMemory);

// Film2d Method Definitions

Film2d::Film2d(Film *film)
    : Film(film) {}

Bounds2i Film2d::GetSampleBounds() const {
    Bounds2f floatBounds(Floor(Point2f(croppedPixelBounds.pMin) +
                               Vector2f(0.5f, 0.0f) - Vector2f(filter->radius.x, 0)),
                         Ceil(Point2f(croppedPixelBounds.pMax) -
                              Vector2f(0.5f, 0.0f) + Vector2f(filter->radius.x, 0)));
    return (Bounds2i)floatBounds;
}

Bounds2f Film2d::GetPhysicalExtent() const {
    Float aspect = (Float)fullResolution.y / (Float)fullResolution.x;
    Float x = std::sqrt(diagonal * diagonal / (1 + aspect * aspect));
    Float y = 0;//aspect * x;
    return Bounds2f(Point2f(-x / 2, -y / 2), Point2f(x / 2, y / 2));
}

std::unique_ptr<FilmTile> Film2d::GetFilmTile(const Bounds2i &sampleBounds) {
    // Bound image pixels that samples in _sampleBounds_ contribute to
    Vector2f halfPixel = Vector2f(0.5f, 0.5f);
    Bounds2f floatBounds = (Bounds2f)sampleBounds;
    Point2i p0 = (Point2i)Ceil(floatBounds.pMin - halfPixel - filter->radius);
    Point2i p1 = (Point2i)Floor(floatBounds.pMax - halfPixel + filter->radius) +
                 Point2i(1, 1);
    Bounds2i tilePixelBounds = Intersect(Bounds2i(p0, p1), croppedPixelBounds);
    return std::unique_ptr<FilmTile>(new FilmTile(
        tilePixelBounds, filter->radius, filterTable, filterTableWidth,
        maxSampleLuminance));
}

void Film2d::Clear() {
    for (Point2i p : croppedPixelBounds) {
        Pixel &pixel = GetPixel(p);
        for (int c = 0; c < 3; ++c)
            pixel.splatXYZ[c] = pixel.xyz[c] = 0;
        pixel.filterWeightSum = 0;
    }
}

void Film2d::MergeFilmTile(std::unique_ptr<FilmTile> tile) {
    ProfilePhase p(Prof::MergeFilmTile);
    VLOG(1) << "Merging film tile " << tile->GetPixelBounds();
    std::lock_guard<std::mutex> lock(mutex);
    for (Point2i pixel : tile->GetPixelBounds()) {
        // Merge _pixel_ into _Film2d::pixels_
        const FilmTilePixel &tilePixel = tile->GetPixel(pixel);
        Pixel &mergePixel = GetPixel(pixel);
        Float xyz[3];
        tilePixel.contribSum.ToXYZ(xyz);
        for (int i = 0; i < 3; ++i) mergePixel.xyz[i] += xyz[i];
        mergePixel.L += tilePixel.contribSum;
        mergePixel.filterWeightSum += tilePixel.filterWeightSum;
    }
}

void Film2d::SetImage(const Spectrum *img) const {
    int nPixels = croppedPixelBounds.Area();
    for (int i = 0; i < nPixels; ++i) {
        Pixel &p = pixels[i];
        img[i].ToXYZ(p.xyz);
        p.filterWeightSum = 1;
        p.splatXYZ[0] = p.splatXYZ[1] = p.splatXYZ[2] = 0;
    }
}

void Film2d::AddSplat(const Point2f &p, Spectrum v) {
    ProfilePhase pp(Prof::SplatFilm);

    if (v.HasNaNs()) {
        LOG(ERROR) << StringPrintf("Ignoring splatted spectrum with NaN values "
                                   "at (%f, %f)", p.x, p.y);
        return;
    } else if (v.y() < 0.) {
        LOG(ERROR) << StringPrintf("Ignoring splatted spectrum with negative "
                                   "luminance %f at (%f, %f)", v.y(), p.x, p.y);
        return;
    } else if (std::isinf(v.y())) {
        LOG(ERROR) << StringPrintf("Ignoring splatted spectrum with infinite "
                                   "luminance at (%f, %f)", p.x, p.y);
        return;
    }

    if (!InsideExclusive((Point2i)p, croppedPixelBounds)) return;
    if (v.y() > maxSampleLuminance)
        v *= maxSampleLuminance / v.y();
    
    for (int k = 0; k < nTimeBins; k++) {
        Float xyz[3];
        v.ToXYZ(xyz, k);
        Pixel &pixel = GetPixel((Point2i)p);
        for (int i = 0; i < 3; ++i) pixel.splatXYZ[k*3+i].Add(xyz[i]);
    }

}

void Film2d::WriteImage(Float splatScale) {
    // Convert image to RGB and compute final pixel values
    LOG(INFO) <<
        "Converting image to RGB and computing final weighted pixel values";
    std::unique_ptr<Float[]> rgb(new Float[3 * nTimeBins * croppedPixelBounds.Area()]);
    int offset = 0;
    for (int k = 0; k < nTimeBins; k++) {
        for (Point2i p : croppedPixelBounds) {
            // Convert pixel XYZ color to RGB
            Pixel &pixel = GetPixel(p);
            Float xyz[3];
            pixel.L.ToXYZ(xyz, k);
            XYZToRGB(xyz, &rgb[3 * offset]);

            // Normalize pixel with weight sum
            Float filterWeightSum = pixel.filterWeightSum;
            if (filterWeightSum != 0) {
                Float invWt = (Float)1 / filterWeightSum;
                rgb[3 * offset] = std::max((Float)0, rgb[3 * offset] * invWt);
                rgb[3 * offset + 1] =
                    std::max((Float)0, rgb[3 * offset + 1] * invWt);
                rgb[3 * offset + 2] =
                    std::max((Float)0, rgb[3 * offset + 2] * invWt);
            }

            // Add splat value at pixel
            Float splatRGB[3];
            Float splatXYZ[3] = {pixel.splatXYZ[k*3+0], pixel.splatXYZ[k*3+1],
                                 pixel.splatXYZ[k*3+2]};
            XYZToRGB(splatXYZ, splatRGB);
            rgb[3 * offset] += splatScale * splatRGB[0];
            rgb[3 * offset + 1] += splatScale * splatRGB[1];
            rgb[3 * offset + 2] += splatScale * splatRGB[2];

            // Scale pixel value by _scale_
            rgb[3 * offset] *= scale;
            rgb[3 * offset + 1] *= scale;
            rgb[3 * offset + 2] *= scale;
            ++offset;
        }
    }
    // Write RGB image
    LOG(INFO) << "Writing image " << filename << " with bounds " <<
        croppedPixelBounds;
    char buffer[256];
    sprintf(buffer,"%s",filename.c_str());
    croppedPixelBounds.pMax.y *= nTimeBins;
    Point2i streakResolution = fullResolution;
    streakResolution.y *= nTimeBins;
    pbrt::WriteImage(buffer, &rgb[0], croppedPixelBounds, streakResolution);
}

}  // namespace pbrt
